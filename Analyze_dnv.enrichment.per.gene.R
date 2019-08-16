sample_path="proband.346ee.txt"
enst_path="enst_of_snpeff.canonical_gencode.protein.coding.txt"
enst_pos_path="enst_gencode.v19.gff3.txt"
obs_path="SHP2MC_NonNMD.bed" 
exp_path="chrallDMS_NonNMD.recode.vcf" 
v4_path="DepthV40417M.txt"
v5_path="DepthV50417M.txt"
v6_path="DepthV60417M.txt"

library(tidyverse)
source("Variable.R")
sample = read_tsv(sample_path,col_names="sample")
enst = read_tsv(enst_path,col_names="enst")
enst_pos = read_tsv(enst_pos_path,col_types="cc__d_____",col_names=c("enst","chr","pos"))
obs = read_tsv(obs_path,col_types="c_dccccdddcdcc")
exp = read_tsv(exp_path,col_types="cddcc__c_",comment="#",col_names=c("chr","pos","rate","ref","alt","an"))
v4 = read_tsv(v4_path,col_types="cdd",col_names=c("chr","pos","v4")) %>% filter(v4!=0)
v5 = read_tsv(v5_path,col_types="cdd",col_names=c("chr","pos","v5")) %>% filter(v5!=0)
v6 = read_tsv(v6_path,col_types="cdd",col_names=c("chr","pos","v6")) %>% filter(v6!=0)

obs2 = inner_join(obs,sample,by="sample") %>%
inner_join(enst,by="enst") %>%
modify_consequence("ano")

obs_snv= obs2 %>%
filter(vqlsod>-7.18 & td>5.72 & dn>0.296 & df=="True" & dg>0.02) %>%
filter(mod %in% c("NON"))

obs_indel = obs2 %>%
filter(vqlsod>-1.06 & td>5.5 & df=="True") %>%
filter(mod %in% c("FS"))
                        
exp2 = separate(exp,an,c("n1","an","n2","n3","n4","n5","enst"),"\\|") %>%
    select(-c(n1,n2,n3,n4,n5)) %>%
    inner_join(enst,by="enst") %>%
    left_join(v4,by=c("chr","pos")) %>%
    left_join(v5,by=c("chr","pos")) %>%
    left_join(v6,by=c("chr","pos")) %>%
    mutate_all(funs(ifelse(is.na(.),0,.))) %>%
    modify_consequence("an") %>%
    filter(mod %in% c("NON"))

exp3 = mutate(exp2,
        v4=if_else(v4>40,40,v4),
        v5=if_else(v5>40,40,v5),
        v6=if_else(v6>40,40,v6)) %>%
    mutate(adj_rate=16*2*rate*v4/40+196*2*rate*v5/40+134*2*rate*v6/40) %>%
    group_by(enst) %>%
    summarise(count=sum(adj_rate*2.25)) %>%
    mutate(type="exp")

tmp = rbind(obs_snv,obs_indel) %>%
    group_by(enst) %>%
    summarise(count=n()) %>%
    mutate(type="obs") %>%
    rbind(exp3) %>%
    spread(key=type,value=count) %>%
    mutate_at(vars(obs),funs(ifelse(is.na(.),0,.))) %>%
    filter(!is.na(exp))

low=c(); upp=c(); pva=c()
for( i in 1:nrow(tmp)){
    res = poisson.test(round(tmp$exp[i]))    
    low = c(low,res$conf.int[1])
    upp = c(upp,res$conf.int[2])
    pva = c(pva,1-ppois(tmp$obs[i]-1,tmp$exp[i]))
}

distinct(enst_pos,enst,.keep_all=TRUE) %>%
    right_join(cbind(tmp,pva)) %>%
    write_tsv("0625_forQQman_346_nonNMD.txt",quote=F)
