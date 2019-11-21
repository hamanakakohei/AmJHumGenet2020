NonNMD=gencode.v19.comprehensive.ucsc.nonnmd.bed

# non-NMD two exon => non-NMD region
awk -F"\t" 'BEGIN{OFS="\t"}NR==1{ CHR[$4]=$1; STA[$4]=$2 ; end[$4]=$3 }NR>1{
    CHR[$4]=$1
    if($4 in STA){if($2 < STA[$4])STA[$4]=$2}else{STA[$4]=$2}
    if($4 in end){if($3 > end[$4])end[$4]=$3}else{end[$4]=$3}
}END{ for(i in CHR)print CHR[i],STA[i],end[i],i
}' $NonNMD > gencode.v19.comprehensive.ucsc.nonnmd.region.bed

