library(tidyverse)
DNV.1330TRIO_PATH="SHP2MC.txt"
kit_path="capture1342trio20190320.txt"
v4_path="DepthV41017M.txt"
v5_path="DepthV51017M.txt"
v6_path="DepthV61017M.txt"
SYN=c("splice_region_variant&synonymous_variant","synonymous_variant","stop_retained_variant","splice_region_variant&initiator_codon_variant","initiator_codon_variant","splice_region_variant&stop_retained_variant")                                                               
OUT="1017_depth.group.vs.count.1330trio.group.per.total.mut.rate.1330trio.group.png"
read_tsv(DNV.1330TRIO_PATH,col_types="ccdcc_cdddcd_") %>% mutate(chr=paste0("chr",chr)) -> dnv.1330trio
read_tsv(kit_path,col_names=c("sample","kit")) %>% mutate(kit=if_else(kit=="51Mb" | kit=="71Mb","V4",kit)) -> kit
read_tsv(v4_path,col_names=c("chr","pos","depth.v4"),col_types="cdd") -> chr_pos_depth.v4 
read_tsv(v5_path,col_names=c("chr","pos","depth.v5"),col_types="cdd") -> chr_pos_depth.v5
read_tsv(v6_path,col_names=c("chr","pos","depth.v6"),col_types="cdd") -> chr_pos_depth.v6

readRDS("tmp1.rds") %>% 
    filter(an %in% SYN) %>%
    left_join(chr_pos_depth.v4,by=c("chr","pos")) %>%
    left_join(chr_pos_depth.v5,by=c("chr","pos")) %>%
    left_join(chr_pos_depth.v6,by=c("chr","pos")) %>%
    gather(key=kit,value=depth,depth.v4,depth.v5,depth.v6) %>%
    mutate(total.mut.rate.each.kit=case_when(
        kit=="depth.v4" ~ mut.rate * 121 * 2,
        kit=="depth.v5" ~ mut.rate * 886 * 2,
        kit=="depth.v6" ~ mut.rate * 323 * 2)) %>%
    filter(!is.na(depth)) %>%
    group_by(depth) %>%
    summarise(total.mut.rate.1330trio=sum(total.mut.rate.each.kit)) -> depth_total.mut.rate.1330trio

filter(dnv.1330trio,vqlsod>-7.18 & td>5.71 & dn>0.196 & df=="True" & dg>0.02) %>%
    filter(ano %in% SYN) %>%
    select(c(sample,chr,pos,ano)) %>%
    left_join(kit,by="sample") %>%
    left_join(chr_pos_depth.v4,by=c("chr","pos")) %>%
    left_join(chr_pos_depth.v5,by=c("chr","pos")) %>%
    left_join(chr_pos_depth.v6,by=c("chr","pos")) %>%
    mutate(depth=case_when(
        kit == "V4" ~ depth.v4,
        kit == "V5" ~ depth.v5,
        kit == "V6" ~ depth.v6)) %>%
    filter(!is.na(depth)) %>%
    group_by(depth) %>%
    summarize(count.1330trio=n()) %>%
    mutate(count.1330trio=as.numeric(count.1330trio))-> depth_count.1330trio

full_join(depth_total.mut.rate.1330trio,depth_count.1330trio) %>%
    #mutate(count.1330trio=if_else(is.na(count.1330trio),0,count.1330trio),depth.group=((depth %/% 5) * 5 + 2.5)) %>%
    mutate(count.1330trio=if_else(is.na(count.1330trio),0,count.1330trio),depth.group=((depth %/% 10) * 10 + 5)) %>%
    group_by(depth.group) %>%
    summarise(count.1330trio.group=sum(count.1330trio),total.mut.rate.1330trio.group=sum(total.mut.rate.1330trio)) %>%
    ggplot(aes(x=depth.group,y=count.1330trio.group/total.mut.rate.1330trio.group)) + 
    geom_point(aes(size=count.1330trio.group)) +
    scale_size_continuous(range=c(0,10),name="Count") + 
    theme(legend.position=c(0.9,0.2)) +
    xlim(0,100) + ylim(0,1.5)
    ggsave(file=OUT)
                                                       
