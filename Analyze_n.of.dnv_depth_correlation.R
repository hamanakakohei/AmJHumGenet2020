# in R
library(tidyverse)
obs_path = "SHP2MC.txt" # obaserved denovo variant file
kit_path = "sample_capturekit.txt"
v4_path = "chr_pos_median.depth.v4.txt"
v5_path = "chr_pos_median.depth.v5.txt"
v6_path = "chr_pos_median.depth.v6.txt"
exp_path = "chrallDMS.vcf-like"

obs = read_tsv(obs_path,col_types="ccdcc_cdddcd_") %>% mutate(chr=paste0("chr",chr))
kit = read_tsv(kit_path,col_names=c("sample","kit"))
dep_v4 = read_tsv(v4_path,col_names=c("chr","pos","v4"),col_types="cdd")
dep_v5 = read_tsv(v5_path,col_names=c("chr","pos","v5"),col_types="cdd")
dep_v6 = read_tsv(v6_path,col_names=c("chr","pos","v6"),col_types="cdd")
exp = read_tsv(exp_path,comment="#",col_types="cddcc__c_",col_names=c("chr","pos","exp","ref","alt","ano"))

snv = c("splice_region_variant&synonymous_variant","synonymous_variant","splice_region_variant&initiator_codon_variant","initiator_codon_variant","stop_retained_variant")                                                                 

count_group = filter(obs, str_length(ref)==1 & str_length(alt)==1) %>%
    filter(vqlsod>-7.18 & td>5.71 & dn>0.196 & df=="True" & dg>0.02) %>%
    filter(ano %in% snv) %>%
    select(c(sample,chr,pos,ano)) %>%
    left_join(kit,by="sample") %>%
    left_join(dep_v4,by=c("chr","pos")) %>%
    left_join(dep_v5,by=c("chr","pos")) %>%
    left_join(dep_v6,by=c("chr","pos")) %>%
    mutate(dep=case_when(
        kit == "V4" ~ v4,
        kit == "V5" ~ v5,
        kit == "V6" ~ v6
    )) %>%
    mutate(group=((dep + 5)%/% 10)*10) %>%
    group_by(group) %>%
    summarize(count=n())

exp2 = separate(exp,ano,c("n","ano"),sep="\\|") %>%
    filter(ano %in% snv) %>%
    select(-c(n,ref,alt)) %>%
    left_join(dep_v4,by=c("chr","pos")) %>%
    left_join(dep_v5,by=c("chr","pos")) %>%
    left_join(dep_v6,by=c("chr","pos")) %>%
    mutate_all(funs(ifelse(is.na(.),0,.))) 

expv4_group = mutate(exp2,group=((v4 + 5)%/% 10)*10) %>% group_by(group) %>% summarize(expv4=121*2*sum(exp))
expv5_group = mutate(exp2,group=((v5 + 5)%/% 10)*10) %>% group_by(group) %>% summarize(expv5=886*2*sum(exp))
expv6_group = mutate(exp2,group=((v6 + 5)%/% 10)*10) %>% group_by(group) %>% summarize(expv6=323*2*sum(exp))

expt_group = full_join(expv4_group,expv5_group,by="group") %>%
    full_join(expv6_group,by="group") %>%
    mutate_all(funs(ifelse(is.na(.),0,.))) %>%
    mutate(expt=expv4+expv5+expv6)

left_join(expt_group,count_group,by="group") %>%
    mutate_at(vars(count),funs(ifelse(is.na(.),0,.))) %>%
    ggplot(aes(x=group,y=count/expt,size=count)) +
    geom_point() +
    scale_size_continuous(range=c(1,15)) +
    theme(legend.position="none") +
    xlim(0,100) + 
    ylim(0,1.5)

ggsave(file="depth_oe.ratio.syn.png")
