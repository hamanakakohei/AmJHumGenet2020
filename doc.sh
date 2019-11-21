FUNCTION=function/
SAMTOOLS=samtools-1.9/bin/samtools
bed_path="gencode.v19.comprehensive.ucsc.bed"
CANONICAL_CODING_GENCODE_NOALTHAPLO="enst.coding.canonical.gencode.noalthaplo.txt"

# bed for depth calc.
join <(awk -F"[\t.]" 'NR>1{print $4,$1,$2,$3,$4}' $bed_path|sort) <(sort $CANONICAL_CODING_GENCODE_NOALTHAPLO)|awk '{print $2"\t"$3-21"\t"$4+20}'|sed -e 's/chr//g'|sort -k1,1 -k2,2n|tr ' ' '\t' > interval1016.bed

# bma path list
cut -f1 1342trioCapture.txt > 1342trioSample.txt
${FUNCTION}GetBamPath.sh 1342trioSample.txt 1017

# add path info to capture info 
#join -t$'\t' <(sort 1342trioCapture.txt) <(sort BamPathExistGetBamPath.sh04011400.txt) > 1342trioCaptureBam.txt
join -t$'\t' <(sort 1342trioCapture.txt) <(sort BamPathExistGetBamPath.sh1017.txt) > 1342trioCaptureBam.txt
grep -v -e titan -e phoenix 1342trioCaptureBam.txt > 1342trioCaptureBamNas.txt

# get 100 sample from each cap
#2019/10/16: V4: 353; V5: 2410; V6: 1238
awk '$2=="V4"' 1342trioCaptureBamNas.txt|awk 'NR%7==2' > DepthCheckV4.txt
awk '$2=="V5"' 1342trioCaptureBamNas.txt|awk 'NR%48==0' > DepthCheckV5.txt
awk '$2=="V6"' 1342trioCaptureBamNas.txt|awk 'NR%24==0' > DepthCheckV6.txt

# bam list
cut -f3 DepthCheckV4.txt > BamListV4.txt
cut -f3 DepthCheckV5.txt > BamListV5.txt
cut -f3 DepthCheckV6.txt > BamListV6.txt

# calc. depth (later, must add "chr"
$SAMTOOLS depth -q 10 -Q 20 -b interval1016.bed -f BamListV4.txt > DepthV41016.txt
$SAMTOOLS depth -q 10 -Q 20 -b interval1016.bed -f BamListV5.txt > DepthV51016.txt
$SAMTOOLS depth -q 10 -Q 20 -b interval1016.bed -f BamListV6.txt > DepthV61016.txt

# median calc. in R (and you need modify a bit, awk 'NR>1{print "chr"$2"\t"$3"\t"$1}' DepthV40417M.txt > tmpv4.txt
library(tidyverse)
v4 = read_tsv("DepthV41016.txt",col_names=c("chr","pos",1:51),col_types=paste(c("c",rep("d",52)),collapse=""))
v4 %>% select(3:53) %>% apply(1,median) %>% cbind(select(v4,1:2)) %>% mutate(chr=paste0("chr",chr)) %>% write.table("DepthV41017M.txt",quote=F,row.names=F,col.names=F,sep="\t")
v5 = read_tsv("DepthV51016.txt",col_names=c("chr","pos",1:50),col_types=paste(c("c",rep("d",51)),collapse=""))
v5 %>% select(3:52) %>% apply(1,median) %>% cbind(select(v5,1:2)) %>% mutate(chr=paste0("chr",chr)) %>% write.table("DepthV51017M.txt",quote=F,row.names=F,col.names=F,sep="\t")
v6 = read_tsv("DepthV61016.txt",col_names=c("chr","pos",1:51),col_types=paste(c("c",rep("d",52)),collapse=""))
v6 %>% select(3:53) %>% apply(1,median) %>% cbind(select(v6,1:2)) %>% mutate(chr=paste0("chr",chr)) %>% write.table("DepthV61017M.txt",quote=F,row.names=F,col.names=F,sep="\t")
