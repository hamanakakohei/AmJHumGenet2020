GATK="java -jar GenomeAnalysisTK.jar"
REF="hs37d5.fa"
SNPEFF="java -jar snpEff.jar"
source VcfFunctions.sh

# pass
$GATK -T SelectVariants -R $REF -ef --variant ExAC.r0.3.nonpsych.sites.vcf.gz --out P.vcf > log.txt 2>&1

# biallelic
$GATK -T SelectVariants -R $REF --restrictAllelesTo BIALLELIC --variant P.vcf --out PB.vcf > log2.txt 2>&1

# MAF filter
af_filter_for_exac PB.vcf 0.0001 > PBM.vcf

# snpeff
java -jar snpEff.jar -v -canon GRCh37.75 PBM.vcf 1> PBMS.vcf 2> log3.txt

# resolve multiple annos for each variant
for_multi_anno_norm PBMS.vcf MULTIPLE > PBMSN.vcf

# select syn
for_multi_anno_select PBMSN.vcf AnoSyn.txt > PBMSNS.vcf

# doc anno at 6th col. (PBMSNS.vcf has no "chr", so PosDep.txt dont need "chr" here
vcf_anno PBMSNS.vcf PosDep.txt 6 0.0 > PBMSNSDM.vcf

# only non-NMD region
vcftools --vcf PBMSNSDM.vcf --out PBMSNSDMN --bed GencodeV19_NonNMDNoChr.bed --recode --recode-INFO-all

# sum count according to ENST
for_multi_anno_sum_count PBMSNSDM.vcf > EnstObs.txt
for_multi_anno_sum_count PBMSNSDMN.recode.vcf > EnstObs_NonNMD.txt

# sum count according to depth
for_multi_anno_sum_count_depth PBMSNSDM.vcf DEPTH > DepthObs.txt

# sequence of exons of gencode transcript => all possible variant + mutation rate at 3rd column
sequence_to_all_variant2.R enst_chr_start_end_sequence.margin.40bp.txt tmp.vcf
awk 'BEGIN{OFS="\t"}!a[$1"-"$2"-"$4"-"$5]++' tmp.vcf > chrall.vcf-like

# doc at 6th col. chrall.vcf-like PosDepChr.txt 6 0.0 > chrallDM.vcf-like
vcf_anno chrall.vcf-like PosDepChr.txt 6 0.0 > chrallDM.vcf-like

# snpeff annotation
$SNPEFF -v -canon GRCh37.75 chrallDM.vcf-like > chrallDMS.vcf-like 2> log.txt

# resolve snpeff anno
for_multi_anno_norm chrallDMS.vcf-like MULTIPLE > chrallDMSN.vcf-like

### in R ####
library(tidyverse)
syn.anos = c("initiator_codon_variant","splice_region_variant&initiator_codon_variant","splice_region_variant&stop_retained_variant","splice_region_variant&synonymous_variant","stop_retained_variant","synonymous_variant")
read_tsv("chrallDMSN.vcf-like",col_types="cddccddccc",col_names=c("chr","pos","rate","ref","alt","dep","DNAm","an1","gt","an2"),comment="#") %>%
    mutate(an3=an2) %>%
    separate(an3,"an3","\\|") %>%
    filter(an3 %in% syn.anos)) %>%
    select(-an3) %>%
    write_tsv("chrallDMSNS2.vcf-like",col_names=FALSE)
##############

# only non-NMD region 
vcftools --vcf chrallDMSNS2.vcf-like --out chrallDMSNSN --bed GencodeV19_NonNMD.bed --recode --recode-INFO-all

# sum mutability according to ENST
for_multi_anno_sum_mutability_depth chrallDMSNS.vcf-like SIMPLE > EnstExp.txt
for_multi_anno_sum_mutability_depth chrallDMSNSN.recode.vcf SIMPLE > EnstExp_NonNMD.txt

# sum mutability according to depth
for_multi_anno_sum_mutability_depth chrallDMSNS.vcf-like DEPTH2 > DepthExp.txt

# sum mutability adjusted with depth accroding to ENST
for_multi_anno_sum_mutability_depth chrallDMSNS.vcf-like DEPTHADJ > EnstExp-DepAdj.txt
for_multi_anno_sum_mutability_depth chrallDMSNSN.recode.vcf DEPTHADJ > EnstExp_NonNMDDepAdj.txt
