library(tidyverse)
bed_path="gencode.v19.exon.bed"

# alt chr has >7 cols but no probelm here
bed = read_tsv(bed_path,col_types="cddc_c",col_names=c("chr","sta","end","enst","str")) %>%
    separate(enst,c("enst","n1","ex","n2","n3","n4","n5"),sep="_") %>%
    separate(enst,c("enst","n6"),sep="\\.") %>%
    select(-c(n1,n2,n3,n4,n5,n6))

# n of exons in each ENST (t_ex)
dup = table(bed$enst)
t_ex = tibble(t_ex=dup,enst=names(dup))

# remove alt haplotype
last_penu_1ex = left_join(bed,t_ex,key="enst") %>%
    mutate(
        an=case_when(
            .$t_ex == 1 ~ "1ex",
            .$t_ex != 1 & .$str == "+" & .$ex == (.$t_ex - 1) ~ "last",
            .$t_ex != 1 & .$str == "+" & .$ex == (.$t_ex - 2) ~ "penu",
            .$t_ex != 1 & .$str == "-" & .$ex == 0 ~ "last",
            .$t_ex != 1 & .$str == "-" & .$ex == 1 ~ "penu"
        )
    )%>%
    filter(an %in% c("1ex","last","penu")) %>%
    mutate(
        sta2=case_when(
            .$an %in% c("1ex","last") ~ .$sta,
            .$an == "penu" & .$str == "+" & (.$end - .$sta) >= 49 ~ .$end - 49,
            .$an == "penu" & .$str == "+" & (.$end - .$sta) < 49 ~ .$sta,
            .$an == "penu" & .$str == "-" ~ .$sta
        ),
        end2=case_when(
            .$an %in% c("1ex","last") ~ .$end,
            .$an == "penu" & .$str == "+" ~ .$end,
            .$an == "penu" & .$str == "-" & (.$end - .$sta) >= 49 ~ .$sta + 49,
            .$an == "penu" & .$str == "-" & (.$end - .$sta) < 49 ~ .$end
        )
    ) %>%
    filter(!grepl("_",chr))

last_penu_1ex %>%
    select(chr,sta2,end2,enst,str,an) %>% 
    write_tsv("GencodeV19_NonNMD.bed",quote=F,col_names=F)
