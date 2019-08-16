NonNMD="GencodeV19_NonNMD.bed"
EXP="chrallDMS.vcf-like"
OBS="SHP2MC.txt"

# non-NMD two exon => non-NMD region
awk -F"\t" 'BEGIN{OFS="\t"}NR==1{ CHR[$4]=$1; STA[$4]=$2 ; end[$4]=$3 }NR>1{
    CHR[$4]=$1
    if($4 in STA){if($2 < STA[$4])STA[$4]=$2}else{STA[$4]=$2}
    if($4 in end){if($3 > end[$4])end[$4]=$3}else{end[$4]=$3}
}END{ for(i in CHR)print CHR[i],STA[i],end[i],i
}' $NonNMD > GencodeV19_NonNMD_region.bed

# non-NMD region filter 
vcftools --vcf $EXP --out chrallDMS_NonNMD --bed GencodeV19_NonNMD_region.bed --recode --recode-INFO-all

# observed in bed format + non-nmd region + manually add header
awk 'BEGIN{OFS="\t"}NR>1{print "chr"$2,$3-1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$1}' $OBS > SHP2MC.bed
bedtools intersect -a SHP2MC.bed -b GencodeV19_NonNMD_region.bed -wa|sort -u > SHP2MC_NonNMD.bed
