sample346_path="IncludeSampleProband.txt"
DNV.1330TRIO_PATH="SHP2MC.txt"
ENST_OF_INTEREST_PATH="enst.coding.canonical.gencode.noalthaplo.txt"
chr_pos_depth.v4_path="DepthV41017M.txt"
chr_pos_depth.v5_path="DepthV51017M.txt"
chr_pos_depth.v6_path="DepthV61017M.txt"
sample_kit_path="Capture.txt"
ENST_NMD.NEG.BED_PATH="gencode.v19.comprehensive.ucsc.nonnmd.region.bed"
chr_start_end_enst_path="enst_start_end.txt"

library(tidyverse)
source("Variable.R")
read_tsv(sample346_path,col_names="sample") %>% mutate(trio346=TRUE) -> sample346
read_tsv(DNV.1330TRIO_PATH,col_types="ccdccccdddcd_") %>% mutate(chr=paste0("chr",chr)) -> dnv.1330trio
read_tsv(ENST_NMD.NEG.BED_PATH,col_names=c("chr.nmd.neg","start","end","enst")) -> enst_chr_start_end
read_tsv(chr_pos_depth.v4_path,col_types="cdd",col_names=c("chr","pos","depth.v4")) -> chr_pos_depth.v4
read_tsv(chr_pos_depth.v5_path,col_types="cdd",col_names=c("chr","pos","depth.v5")) -> chr_pos_depth.v5
read_tsv(chr_pos_depth.v6_path,col_types="cdd",col_names=c("chr","pos","depth.v6")) -> chr_pos_depth.v6
read_tsv(ENST_OF_INTEREST_PATH,col_names="enst") -> enst.of.interest
read_tsv(chr_start_end_enst_path,col_names=c("chr","start","end","enst")) %>% distinct(enst,.keep_all=TRUE) -> chr_start_end_enst
sample_kit = read_tsv(sample_kit_path,col_types="cc_____",skip=1,col_names=c("sample","kit")) %>% mutate(sample=paste0("Sample_",sample))

readRDS("exhaustive.vcf.rds") %>%
    left_join(chr_pos_depth.v4,by=c("chr","pos")) %>%
    left_join(chr_pos_depth.v5,by=c("chr","pos")) %>%
    left_join(chr_pos_depth.v6,by=c("chr","pos")) %>%
    mutate_at(vars(starts_with("depth.v")),funs(ifelse(is.na(.),0,.))) %>%
    modify_consequence("an") %>% 
    mutate(
        depth.v4=if_else(depth.v4>40,40,depth.v4),
        depth.v5=if_else(depth.v5>40,40,depth.v5),
        depth.v6=if_else(depth.v6>40,40,depth.v6)) %>% 
    mutate(rate.adj.for.346trio=16*2*mut.rate*depth.v4/40+196*2*mut.rate*depth.v5/40+134*2*mut.rate*depth.v6/40,
        rate.adj.for.1330trio=121*2*mut.rate*depth.v4/40+886*2*mut.rate*depth.v5/40+323*2*mut.rate*depth.v6/40) %>%
    group_by(enst,mod,nmd.neg) %>% 
    summarise(
        rate.adj.346trio.total=sum(rate.adj.for.346trio),
        rate.adj.1330trio.total=sum(rate.adj.for.1330trio)) %>%
    gather(key=cohort,value=rate.adj.cohort,rate.adj.346trio.total,rate.adj.1330trio.total) %>%
    spread(key=mod,value=rate.adj.cohort) %>%
    mutate(FS=1.25*NON,INF=FS/9) %>%
    gather(key=mod, value=rate.adj.cohort,MIS,NON,SYN,FS,INF,OTH) %>%
    mutate(cohort=if_else(cohort=="rate.adj.1330trio.total","trio1330","trio346"),
        rate.adj.cohort=if_else(is.na(rate.adj.cohort),0,rate.adj.cohort)
    ) %>% ungroup -> enst_mod_cohort_rate.adj.cohort

inner_join(dnv.1330trio,enst.of.interest) %>%
    left_join(sample346) %>%
    left_join(enst_chr_start_end) %>%
    mutate(nmd.neg=if_else(chr==chr.nmd.neg & pos >= start & pos <= end,"TRUE","FALSE")) %>%
    modify_consequence("ano") %>%
    mutate(filter=case_when(
        vqlsod>-7.18 & td>5.71 & dn>0.196 & dg>0.02 & df=="True" & mod %in% c("MIS","NON","SYN") ~ "PASS", #vqslod>-7.18
        vqlsod>-1.06 & td>5.5 & df=="True" & mod %in% c("FS","INF") ~ "PASS",
        TRUE ~ "FILTERED")) %>%
    group_by(enst,mod,nmd.neg) %>%
    filter(filter=="PASS") %>%
    summarise(trio1330=n(),trio346=sum(trio346=="TRUE",na.rm=T)) %>%
    gather(key=cohort,value=observed.count,trio1330,trio346) %>% ungroup -> enst_mod_cohort_observed.count

full_join(enst_mod_cohort_observed.count,enst_mod_cohort_rate.adj.cohort,by=c("enst","cohort","mod","nmd.neg")) %>%
    mutate(observed.count=if_else(is.na(observed.count),0,as.numeric(observed.count))) %>%
    nest(-c(enst,mod,cohort)) %>%
    mutate(data2=map(data,function(x)rbind(x,list("ALL",sum(x$observed.count),sum(x$rate.adj.cohort))))) %>% 
    unnest(data2) %>% 
    nest(-c(enst,cohort,nmd.neg)) %>%
    mutate(data2=map(data,function(x)rbind(x,list("LOF",
        filter(x,mod %in% c("NON","FS")) %>% select(observed.count) %>% sum,
        filter(x,mod %in% c("NON","FS")) %>% select(rate.adj.cohort) %>% sum)))) %>% 
    unnest(data2) %>% {
        group_by(.,cohort,mod,nmd.neg) %>% 
        summarise(observed.count.allenst=sum(observed.count),rate.adj.cohort.allenst=sum(rate.adj.cohort)) %>%
        mutate(low=map(rate.adj.cohort.allenst,function(x)poisson.test(round(x))$conf.int[1]),
            upp=map(rate.adj.cohort.allenst,function(x)poisson.test(round(x))$conf.int[2]),
            pvalue=map2(observed.count.allenst,rate.adj.cohort.allenst,function(x,y)1-ppois(x-1,y))
        ) %>% unnest(low,upp,pvalue) ->> mod_cohort_nmd.neg_ci_pvalue 
        
        mutate(.,low=map(rate.adj.cohort,function(x)poisson.test(round(x))$conf.int[1]),
            upp=map(rate.adj.cohort,function(x)poisson.test(round(x))$conf.int[2]),
            pvalue=map2(observed.count,rate.adj.cohort,function(x,y)1-ppois(x-1,y))
        ) %>% unnest(pvalue) ->> enst_mod_cohort_nmd.neg_ci_pvalue 
    }


### forrest plof for SYN,MIS,NON in whole coding or NMD(-) region in 1330 trios or NMD(-) region in 346 trios (Fig. 1B 1C) ###
filter(mod_cohort_nmd.neg_ci_pvalue,cohort=="trio1330" & nmd.neg=="ALL") %>% plot_forrest(c("SYN","MIS","NON"),"1017_1330_all.png")
filter(mod_cohort_nmd.neg_ci_pvalue,cohort=="trio1330" & nmd.neg=="TRUE") %>% plot_forrest(c("SYN","MIS","NON"),"1017_133_nmd.neg.png")
filter(mod_cohort_nmd.neg_ci_pvalue,cohort=="trio346" & nmd.neg=="TRUE") %>% plot_forrest(c("SYN","MIS","NON","FS","LOF"),"1017_346_nmd.neg.png")
1-ppois(5-1,2.75e-7*2.25*2*(1752+4293))

plot_forrest = function(DF,MOD,OUT){
    ggplot(DF,aes(x=mod,y=rate.adj.cohort.allenst)) +
        geom_point(size=4) +
        geom_point(aes(y=observed.count.allenst),size=4,color="red") +
        scale_x_discrete(limits=MOD) +
        #ylim(0,0) +
        geom_errorbar(aes(ymax=upp,ymin=low),width=0)
        ggsave(file=OUT)
}

### plot manhattan for gene-based enrichment analysis for LOF in nmd(-) regions (Fig. 1D) ###
library(qqman)
gene.n=20042
OUT="1017_346_nmd.neg_qqman.png"
par(ps=8)
png(OUT,width=3.5,height=3.5,units="in",res=100)
enst_mod_cohort_nmd.neg_ci_pvalue %>%
    unnest(pvalue) %>%
    inner_join(chr_start_end_enst,by="enst") %>%
    filter(cohort=="trio346" & nmd.neg=="TRUE" &  mod=="LOF" & !chr %in% c("chrY","chrM")) %>% 
    mutate(chr=if_else(chr=="chrX","chr23",chr),pos=(start+end)/2,
        chr=as.numeric(str_replace(chr,pattern=c("chr"),replacement=c("")))) %>%
    select(enst,chr,pos,pvalue) %>%
    rename(SNP=enst,CHR=chr,BP=pos,P=pvalue) %>%
    manhattan(suggestiveline=FALSE,genomewideline=-log10(0.05/gene.n),cex=1,ylim=c(0,10))
    dev.off()
