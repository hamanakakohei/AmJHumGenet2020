sample_path="proband.346ee.txt"
enst_path="enst_of_snpeff.canonical_gencode.protein.coding.txt" # snpeff canonical gencode HQ ENST
#obs_path="SHP2MC.bed" 
obs_path="SHP2MC_NonNMD.bed" 
#exp_path="chrallDMS.vcf-like" 
exp_path="chrallDMS_NonNMD.recode.vcf" 
v4_path="DepthV40417M.txt"
v5_path="DepthV50417M.txt"
v6_path="DepthV60417M.txt"

library(tidyverse)
source("Variable.R")
sample = read_tsv(sample_path,col_names="sample")
enst = read_tsv(enst_path,col_names="enst")
obs = read_tsv(obs_path,col_types="c_dccccdddcdcc")
exp = read_tsv(exp_path,col_types="cddcc__c_",comment="#",col_names=c("chr","pos","rate","ref","alt","an"))
v4 = read_tsv(v4_path,col_types="cdd",col_names=c("chr","pos","v4")) %>% filter(v4!=0)
v5 = read_tsv(v5_path,col_types="cdd",col_names=c("chr","pos","v5")) %>% filter(v5!=0)
v6 = read_tsv(v6_path,col_types="cdd",col_names=c("chr","pos","v6")) %>% filter(v6!=0)

obs2 = obs %>%
    #inner_join(sample,by="sample") %>% # for 346ee analysis
    inner_join(enst,by="enst") %>%
    modify_consequence("ano")

obs_snv = obs2 %>%
    filter(vqlsod>-7.18 & td>5.72 & dn>0.296 & df=="True" & dg>0.02) %>%
    filter(mod %in% c("MIS","NON","SYN")) %>%
    group_by(mod) %>%
    summarise(count=n()) 

obs_indel = obs2 %>%
    filter(vqlsod>-1.06 & td>5.5 & df=="True") %>%
    filter(mod %in% c("FS","INF")) %>%
    group_by(mod) %>%
    summarise(count=n())

exp2 = separate(exp,an,c("n1","an","n2","n3","n4","n5","enst"),"\\|") %>%
    select(-c(n1,n2,n3,n4,n5)) %>%
    inner_join(enst,by="enst") %>%
    left_join(v4,by=c("chr","pos")) %>%
    left_join(v5,by=c("chr","pos")) %>%
    left_join(v6,by=c("chr","pos")) %>%
    mutate_all(funs(ifelse(is.na(.),0,.))) %>%
    modify_consequence("an") 

exp3 = mutate(exp2,
        v4=if_else(v4>40,40,v4),
        v5=if_else(v5>40,40,v5),
        v6=if_else(v6>40,40,v6)) %>%
    #mutate(adj_rate=16*2*rate*v4/40+196*2*rate*v5/40+134*2*rate*v6/40) %>% # for 346ee analysis
    mutate(adj_rate=121*2*rate*v4/40+886*2*rate*v5/40+323*2*rate*v6/40) %>% # for 1342 trio analysis
    group_by(mod) %>%
    summarise(count=sum(adj_rate)) %>%
    filter(mod != "OTH") %>%
    spread(key=mod,value=count) %>%
    mutate(FS=1.25*NON,INF=FS/9) %>%
    gather(key=mod, value=count,MIS,NON,SYN,FS,INF) %>%
    mutate(type="exp")

# forrest plot
tmp = rbind(obs_snv,obs_indel) %>%
    mutate(type="obs") %>%
    rbind(exp3) %>%
    #transform(mod=factor(mod,levels=c("SYN","MIS","NON","INF","FS"))) %>%
    filter(!(mod %in% c("INF","FS"))) %>%
    #filter(!(mod %in% c("INF","MIS"))) %>%
    mutate(mod=as.character(mod)) %>%
    spread(key=type,value=count)

## lof option
#tmp=rbind(tmp,list("LOF",subset(tmp,mod=="NON")$exp+subset(tmp,mod=="FS")$exp,subset(tmp,mod=="NON")$obs+subset(tmp,mod=="FS")$obs))

low=c(); upp=c(); pva=c()
for( i in 1:nrow(tmp)){
    res = poisson.test(round(tmp$exp[i]))    
    low = c(low,res$conf.int[1])
    upp = c(upp,res$conf.int[2])
    pva = c(pva,1-ppois(tmp$obs[i]-1,tmp$exp[i]))
}

cbind(tmp,low,upp,pva) %>%
    #transform(mod=factor(mod,levels=c("SYN","NON","FS","LOF"))) %>%
    ggplot(aes(x=mod,y=exp)) +
    #theme(axis.text=element_text(size=6)) +
    geom_point(size=6) +
    geom_point(aes(y=obs),size=6,color="red") +
    scale_x_discrete(limits=c("SYN","MIS","NON")) +
    ylim(0,40) +
    geom_errorbar(aes(ymax=upp,ymin=low),width=0.01)
    ggsave(file="0626_forrest_nonnmd_1342.png")
