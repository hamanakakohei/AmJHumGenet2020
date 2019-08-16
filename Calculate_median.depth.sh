# region for depth check (20bp margin
awk '{print $1"\t"$2-21"\t"$3+20}' exon_of_transcript_of_snpeff.canonical_gencode.protein.coding.bed|sed -e 's/chr//g'|sort -k1,1 -k2,2n > interval.bed

# calc. depth
samtools depth -q 10 -Q 20 -b interval.bed -f bam.path.v4.txt > chr_pos_depths.v4.txt
samtools depth -q 10 -Q 20 -b interval.bed -f bam.path.v5.txt > chr_pos_depths.v5.txt
samtools depth -q 10 -Q 20 -b interval.bed -f bam.path.v6.txt > chr_pos_depths.v6.txt

### median calc. in R ##### 
library(tidyverse)
read_tsv("chr_pos_depths.v4.txt",col_names=c("chr","pos",1:104),col_types=paste(c("c",rep("d",105)),collapse="")) %>% 
    select(3:106) %>% apply(1,median) %>% cbind(select(v4,1:2)) %>% write.table("chr_pos_median.depth.v4.txt",quote=F,row.names=F,sep="\t")
read_tsv("chr_pos_depths.v5.txt",col_names=c("chr","pos",1:101),col_types=paste(c("c",rep("d",102)),collapse="")) %>% 
    select(3:103) %>% apply(1,median) %>% cbind(select(v5,1:2)) %>% write.table("chr_pos_median.depth.v5.txt",quote=F,row.names=F,sep="\t")
read_tsv("chr_pos_depths.v6.txt",col_names=c("chr","pos",1:103),,col_types=paste(c("c",rep("d",104)),collapse="")) %>% 
    select(3:105) %>% apply(1,median) %>% cbind(select(v6,1:2)) %>% write.table("chr_pos_median.depth.v6.txt",quote=F,row.names=F,sep="\t")
