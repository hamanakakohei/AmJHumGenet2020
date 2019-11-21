GATK="java -jar GenomeAnalysisTK.jar"
REF=hs37d5.fa
VCFTOOLS=vcftools-0.1.15/bin/vcftools
BEDTOOLS=bedtools2-2.26.0/bin/bedtools
BCFTOOLS="bcftools-1.8/bin/bcftools"
FUNCTION=function/
SEQUENCE=gencodev19.comprehensive.40bp.seq.txt
CANONICAL_CODING_GENCODE_NOALTHAPLO="enst.coding.canonical.gencode.noalthaplo.txt"
EXAC="ExAC.r0.3.nonpsych.sites.vcf.gz"
SNPEFF="java -Xmx64g -jar snpEff.jar"
NMD_NEG_BED="gencode.v19.comprehensive.ucsc.nonnmd.region.bed"

source ${FUNCTION}Variable.sh
prep_seq_file $SEQUENCE|\
    awk '{split($1,col1,".");print col1[1],$2,$3,$4,$5}'|sort|\
    #join - <(sort $CANONICAL)|\
    join - <(sort $CANONICAL_CODING_GENCODE_NOALTHAPLO)|\
    tr ' ' '\t' > seq_canonical.protein.enst.txt

cat chr.txt|while read CHR; do awk '$2=="'"${CHR}"'"' seq_canonical.protein.enst.txt > seq_canonical.protein.enst.${CHR}.txt; done
cat chr.txt|while read CHR; do ${FUNCTION}SequenceToAllVariant2.R seq_canonical.protein.enst.${CHR}.txt tmp${CHR}.txt; done
cat tmpchr*|sort -u > v.vcf
$SNPEFF -v -onlyProtein -canon GRCh37.75 v.vcf > vs.vcf 2> log.txt
$VCFTOOLS --vcf vs.vcf --out vsn --bed $NMD_NEG_BED --recode --recode-INFO-all
#${FUNCTION}ForMultiAnnoNorm.sh vs.vcf MULTIPLE > vsn.vcf

# exac data
less $EXAC|awk -F"\t" 'BEGIN{OFS="\t"}$1 ~ /^#/{print $0}$1 !~ /^#/ && $7=="PASS"{$8=".";print $0}' > e1.vcf
$BCFTOOLS norm -m -any -f $REF e1.vcf > e2.vcf

#
library(tidyverse)
VCF_PATH="vs.vcf"
VCF.NMD.NEG_PATH="vsn.recode.vcf"
ENST_OF_INTEREST_PATH="enst.coding.canonical.gencode.noalthaplo.txt"
ENST_NMD.NEG.BED_PATH="gencode.v19.comprehensive.ucsc.nonnmd.region.bed"
chr_pos_mean_median_exac_PATH="Panel.chrall.coverage.txt"
EXAC_PATH="e2.vcf"
SYN=c("splice_region_variant&synonymous_variant","synonymous_variant","stop_retained_variant","splice_region_variant&initiator_codon_variant","initiator_codon_variant","splice_region_variant&stop_retained_variant")

read_tsv(ENST_OF_INTEREST_PATH,col_names="enst") -> enst.of.interest
read_tsv(ENST_NMD.NEG.BED_PATH,col_names=c("chr.nmd.neg","start","end","enst")) -> enst_chr_start_end
read_tsv(chr_pos_mean_median_exac_PATH,comment="#",col_types="cd_d",col_names=c("chr","pos","depth.exac")) %>% mutate(chr=paste0("chr",chr)) -> chr_pos_depth.exac
read_tsv(EXAC_PATH,comment="#",col_types="cd_cc",col_names=c("chr","pos","ref","alt")) %>% mutate(chr=paste0("chr",chr),exac="exac") -> exac
read_tsv(VCF_PATH,comment="#",col_types="cddcc__c_",col_names=c("chr","pos","mut.rate","ref","alt","snpeff")) %>%
    separate(snpeff,into=c("NA1","an","NA2","NA3","NA4","NA5","enst"),sep="\\|") %>%
    select(-c(NA1,NA2,NA3,NA4,NA5)) %>%
    inner_join(enst.of.interest,by="enst") %>% # tmp1.rds
    left_join(chr_pos_depth.exac,by=c("chr","pos")) %>%
    left_join(exac,by=c("chr","pos","ref","alt")) %>%
    left_join(enst_chr_start_end) %>%
    mutate(nmd.neg=if_else(chr==chr.nmd.neg & pos >= start & pos <= end,"TRUE","FALSE")) %>%
    saveRDS("exhaustive.vcf.rds")
    
readRDS("exhaustive.vcf.rds") %>% mutate(mut.rate.adj=case_when(
        depth.exac <= 12.5 ~ mut.rate * (0.04 * depth.exac + 0.125),
        depth.exac > 12.5 & depth.exac < 62.5 ~ mut.rate * (0.17 * log(depth.exac - 5.7) + 0.3),
        depth.exac >= 62.5 ~ mut.rate)) -> vcf # save as tmp2.rds

### table s1 s2, fig s2a, fig 1a ###
filter(vcf,nmd.neg=="TRUE") -> vcf.nmd.neg
make_mut.rate.table(vcf,enst.of.interest,"1008_enst_mut.rate.txt")
make_mut.rate.table(vcf.nmd.neg,enst.of.interest,"1008_enst_mut.rate.nmd.neg.txt")
plot_depth_vs_count.exac.per.total.mut.rate(vcf,SYN,"1008_depth_count.exac.per.total.mut.rate.png")
plot_total.mut.rate_vs_count.exac_at.each.enst(vcf,SYN,"1008_total.mut.rate_count.exac.png")
plot_total.mut.rate_vs_count.exac_at.each.enst(vcf.nmd.neg,SYN,"1008_total.mut.rate_count.exac_at.nmd.neg.png")

### functions ###
make_mut.rate.table = function(VCF,ENST.OF.INTEREST,OUT){
    group_by(VCF,enst,an) %>%
        summarise(total.mut.rate=sum(mut.rate)) %>%
        spread(key=an,value=total.mut.rate) %>%
        right_join(ENST.OF.INTEREST) %>%
        mutate_all(funs(ifelse(is.na(.),0,.))) %>%
        write_tsv(OUT,quote=F)
}

plot_depth_vs_count.exac.per.total.mut.rate = function(VCF,SYN,OUT){
    filter(VCF,an %in% SYN) %>%
        filter(!is.na(depth.exac)) %>%
        group_by(depth.exac) %>%
        summarise(total.mut.rate=sum(mut.rate),count.exac=sum(exac=="exac",na.rm=T)) %>%
        filter(count.exac>1000) %>%
        mutate(ratio=count.exac/total.mut.rate) %>% {
            filter(.,depth.exac==100) %>% select(ratio) %>% as.numeric ->> base
            ggplot(.,aes(x=depth.exac,y=ratio/base)) + geom_point(size=2) + ylim(0,1.0) +
                geom_abline(intercept=0.125,slope=0.04) +
                stat_function(fun=function(x)0.17*log(x-5.7)+0.3)
                ggsave(file=OUT)
    }
}

plot_total.mut.rate_vs_count.exac_at.each.enst = function(VCF,SYN,OUT){
    filter(VCF,an %in% SYN) %>%
        filter(!is.na(depth.exac)) %>%
        group_by(enst) %>%
        summarise(total.mut.rate.adj=sum(mut.rate.adj),count.exac=sum(exac=="exac",na.rm=T)) %>%
        ggplot(aes(x=total.mut.rate.adj,y=count.exac)) +
        geom_point(size=2)
        ggsave(file=OUT)
}

