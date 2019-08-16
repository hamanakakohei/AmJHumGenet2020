SORTVCF="java -jar picard-tools-2.10.0/picard.jar SortVcf"
BEDTOOLS="bedtools sort"
GATK="java -jar GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar"
REF="hs37d5.fa"
IDX="hs37d5.fa.fai"
VCF="joint_genotyped.vcf"
PED="all.trios.ped"
source VcfFunctions.sh

# sort
cut -f1 $IDX > names.txt
$BEDTOOLS sort -i $VCF -faidx names.txt > S.vcf

# add header
cat_vcf_header S.vcf $VCF > SH.vcf

# filter using parents' variants
vcf_parent_filter2 SH.vcf $PED > SHP2.vcf

# get varinats violating mendelian inheritance by GATK + extract vqslod
extract_variant_vcf SHP2.vcf > SHP2.txt

# triodenovo analysis
triodenovo SHP2.vcf $PED > SHPTd2.txt

# sort ped file and bam list file for DNMfilter

# DNMfilter analysis
dnmfilter SHPTd2.txt $PED $BAMLIST > SHPDn2.txt

# foramtting for denovofilter (this family file has no header, add by yourself
DenovofilterFormat.sh SHP2.vcf $PED SHPDfS2G.txt SHPDfI2G.txt SHPDfF.txt

# denovofilter anlaysis of SHPDfS2G.txt SHPDfI2G.txt by python
# follow github.com/jeremymcrae/denovoFilter
python setup.py install --user
python filter_de_novos.py --de-novos SHPDfS2G.txt --families SHPDfF.txt --annotate-only --output df.snv.txt
python filter_de_novos.py --de-novos-indels SHPDfI2G.txt --families SHPDfF.txt --annotate-only --output df.indel.txt

# formatting the results
cat df.snv.txt <(tail -n+2 df.indel.txt)|awk -F"\t" 'BEGIN{OFS="\t"}NR>1{print $6,$2,$7,$9,$1,$5}' > SHPDf2.txt

# denovogear
dng dnm auto --ped $PED --vcf SHP2.vcf > SHPDg.txt

############## in R #################
#formatting
source("Variable.R")
simplify_denovogear("SHPDg.txt","SHPDg.txt")

# merge denovo soft scores
merge_four_denovo_softs(ga="SHP2.txt",td="SHPTd2.txt",dn="SHPDn2.txt",df="SHPDf2.txt",dg="SHPDg2.txt",pc="Sanger.confirmed.variants.txt",out="SHP2M.txt")

# dist. of scores in R
plot_dnv_soft_score(dt="SHP2MC.txt")

