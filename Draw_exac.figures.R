library(tidyverse)

obs_path="EnstObs.txt"
exp_path="EnstExp.txt"
enst_path="enst_of_transcript_snpeff.canonical_gencode.protein.coding.txt"
exp_adj_path="EnstExp-DepAdj.txt"
obs_nmd_path="EnstObs_NonNMD.txt"
exp_nmd_path="EnstExp_NonNMDDepAdj.txt"

obs = read_tsv(obs_path,col_types="cd",col_names=c("enst","n"))
exp = read_tsv(exp_path,col_types="cd",col_names=c("enst","rate"))
ens = read_tsv(enst_path,col_types="c____",col_names=c("enst")) %>% distinct(enst,.keep_all=TRUE)
exp_adj = read_tsv(exp_adj_path,col_types="cd",col_names=c("enst","rate_adj"))
obs_nmd = read_tsv(obs_nmd_path,col_types="cd",col_names=c("enst","n_nmd"))
exp_nmd = read_tsv(exp_nmd_path,col_types="cd",col_names=c("enst","rate_nmdadj"))

enst_n_rate_len_adj = left_join(ens,obs,by="enst") %>% 
        left_join(exp,by="enst") %>%
        left_join(exp_adj,by="enst") %>%
        left_join(obs_nmd,by="enst") %>%
        left_join(exp_nmd,by="enst") %>%
        mutate_all(funs(ifelse(is.na(.),0,.)))
    
enst_n_rate_len_adj %>%
    ggplot(aes(x=rate_adj,y=n)) +
    geom_point(size=2)
    ggsave(file="DNV.rate_vs_depth.adj.n.png")

enst_n_rate_len_adj %>%
    ggplot(aes(x=rate_nmdadj,y=n_nmd)) +
    geom_point(size=2)
    ggsave(file="DNV.rate_depth.adj.n_at_nmd.negative.region.png")

