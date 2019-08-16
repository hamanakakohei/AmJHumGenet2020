# 1st arg: vcf
# 2nd arg: ped
# 3rd arg: snv output name
# 4th arg: indel output name
# 5th arg: family file output name
# output: DenovoFilter format (sample,chr,pos,ref,alt,gene,dp_pt,dp_fa,dp_mo,anno,gt_pt,gt_fa,gt_mo

source VcfFunctions.sh
echo "WARNING: family id in ped file MUST be proband id!!"

# add gene at last
paste_last_for_snpeffvcf $1 GENE > tmp`basename $0`.txt

# add anno at last
paste_last_for_snpeffvcf tmp`basename $0`.txt ANNO > tmp`basename $0`2.txt

# sample chr;pos;ref;alt;gene;anno dp for trio memb.
# some sample has no AD even though the variant has AD col.
NROW=`head -n1000 tmp$(basename $0)2.txt|awk '$1=="#CHROM"{print NR}'`
awk -F"\t" 'NR=='"${NROW}"'{
    for(i=10;i<=NF;i++){SAMPLE[i]=$i}
} NR>'"${NROW}"' {
    for (i=10;i<=NF-2;i++){
        split($9,FORMAT,":")
        split($i,GT,":")
        for(j=1;j<=length(FORMAT);j++){
            if(FORMAT[j]=="AD"){
                if(GT[j]==""){
                    print SAMPLE[i]"\t"$1";"$2";"$4";"$5";"$(NF-1)";"$NF"\t"0;break
                }else{
                    print SAMPLE[i]"\t"$1";"$2";"$4";"$5";"$(NF-1)";"$NF"\t"GT[j];break
                }
            }
            if(j==length(FORMAT)){print SAMPLE[i]"\t"$1";"$2";"$4";"$5";"$(NF-1)";"$NF"\t""NoDP"}
        }
    }
}' tmp`basename $0`2.txt|sort > tmp`basename $0`2.5.txt

# dependent on vcf? AD: 0 => 0,0
awk 'BEGIN{OFS="\t"}$3!=0{print $0}$3==0{print $1,$2,"0,0"}' tmp`basename $0`2.5.txt > tmp`basename $0`3.txt

# make Sample_ID"\t"Family_ID list
awk '{print $2"\t"$1}' $2|sort > tmp`basename $0`4.txt

# make pt, fa, mo list
cut -f1 $2|sort -u > tmp`basename $0`pt.txt
cut -f3 $2|sort -u|grep -vw 0 > tmp`basename $0`fa.txt
cut -f4 $2|sort -u|grep -vw 0 > tmp`basename $0`mo.txt

# family;chr;pos;ref;alt;gene;anno dp for pt, fa, mo
join -t$'\t' tmp`basename $0`3.txt tmp`basename $0`pt.txt|awk '{print $1";"$2"\t"$3}'|sort > tmp`basename $0`pt2.txt
join -t$'\t' tmp`basename $0`3.txt tmp`basename $0`fa.txt|sort|join -t$'\t' tmp`basename $0`4.txt -|awk '{print $2";"$3"\t"$4}'|sort > tmp`basename $0`fa2.txt
join -t$'\t' tmp`basename $0`3.txt tmp`basename $0`mo.txt|sort|join -t$'\t' tmp`basename $0`4.txt -|awk '{print $2";"$3"\t"$4}'|sort > tmp`basename $0`mo2.txt

# family;chr;pos;ref;alt;gene;anno dp_pt dp_fa dp_mo
paste  tmp`basename $0`pt2.txt <(cut -f2 tmp`basename $0`fa2.txt) <(cut -f2 tmp`basename $0`mo2.txt) > tmp`basename $0`5.txt

# Denovofilter format (for SNV, max_af,pp_dnm,in_child_vcf,in_mother_vcf,in_father_vcf; for INDEL, only max_af
SNV="tmp$(basename $0)6.txt"
INDEL="tmp$(basename $0)7.txt"
awk -F";" -v SNV="${SNV}" -v INDEL="${INDEL}" 'BEGIN{OFS="\t"}{
    if(length($4)==length($5)){$1=$1;print $0,0,0,1,0,0 > SNV}else{$1=$1;print $0,0 > INDEL}    
}' tmp`basename $0`5.txt

# add header
echo -e "person_stable_id\tchrom\tpos\tref\talt\tsymbol\tconsequence\tdp4_child\tdp4_father\tdp4_mother\tmax_af" > tmp`basename $0`8.txt
cat <(paste tmp`basename $0`8.txt <(echo -e "pp_dnm\tin_child_vcf\tin_mother_vcf\tin_father_vcf")) tmp`basename $0`6.txt > $3
cat tmp`basename $0`8.txt tmp`basename $0`7.txt > $4

# make family file
awk '$6==2{print $1"\t"$1"\t"$5}' $2|sed -e 's/2$/F/' -e 's/1$/M/' > $5
if cut -f3 $5|grep -vw -e "F" -e "M">/dev/null; then echo "The family file has samples of unknown sex" 1>&2; fi
