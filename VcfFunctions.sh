PICARD="java -jar picard.jar"
GATK="java -jar GenomeAnalysisTK.jar"


cat_vcf_header(){
    # 1st arg: vcf with no header
    # 2nd arg: vcf with header
    # output: 1st vcf with header of 2nd vcf
    local NROW=`head -n1000 $2|awk '$1=="#CHROM"{print NR}'`
    cat <(head -n $NROW $2) $1
}

remove_vcf_header(){
    # 1st arg: vcf with no header
    # output: vcf with no header
    local NROW=`head -n1000 $1|awk '$1=="#CHROM"{print NR}'`
    tail -n+`expr $NROW + 1` $1
}

select_existing_variant(){
    local NROW_OF_HEADER=`awk '$1=="#CHROM"{print NR}' $1`
    tail -n +`expr 1 + $NROW_OF_HEADER` $1 | \
        awk '{
            for (i=10; i <= NF; i++) {
                if ($i ~ /^0\/1/ || $i ~ /^1\/0/ || $i ~ /^1\/1/) {
                    print $0
                    break
                }
            }
        }' | \
        cat <(head -n $NROW_OF_HEADER $1) -
}

subtract_vcf_from_another_vcf(){
    local NROW_OF_HEADER1=`awk '$1=="#CHROM"{print NR}' $1`
    local NROW_OF_HEADER2=`awk '$1=="#CHROM"{print NR}' $2`
    tail -n +`expr 1 + $NROW_OF_HEADER1` $1| awk -F"\t" 'BEGIN{OFS="\t"}{print $1"-"$2"-"$4"-"$5, $0}'| sort > Tmp`basename $0`1.txt
    tail -n +`expr 1 + $NROW_OF_HEADER2` $2| awk -F"\t" 'BEGIN{OFS="\t"}{print $1"-"$2"-"$4"-"$5}'| sort > Tmp`basename $0`2.txt
    join -t $'\t' -v 1 Tmp`basename $0`1.txt Tmp`basename $0`2.txt |\
        cut -f2- | \
        cat <(head -n $NROW_OF_HEADER1 $1) -
}

vcf_parent_filter2(){
    # 1st arg: vcf
    # 2nd arg: ped
    if [ $# != 2 ]; then echo "ERROR: two args required" 1>&2; fi
    echo "WARNING: not correct about multi-allelic sites" 1>&2
    
    # parent samples
    awk '$4==0{print $2}' $2 > tmp${FUNCNAME[0]}.txt 
    
    # parent vcf
    vcftools --vcf $1 --keep tmp${FUNCNAME[0]}.txt --out tmp${FUNCNAME[0]} --recode --recode-INFO-all   
    
    # parent variant vcf (0/1, 1/0, or 1/1
    select_exising_variant tmp${FUNCNAME[0]}.recode.vcf > tmp${FUNCNAME[0]}2.vcf
    
    # subtract from original vcf
    subtract_vcf_from_another_vcf $1 tmp${FUNCNAME[0]}2.vcf
}

paste_last_for_snpeffvcf(){
    # body
    local NROW=`head -n1000 $1|awk '$1=="#CHROM"{print NR}'`
    tail -n+`expr 1 + $NROW` $1 > tmp`basename $0`.txt

    # something
    if [ $2 == "GENE" ]; then
        awk '{print $8}' tmp`basename $0`.txt|tr ";" "\n"|grep ^ANN=|awk -F"|" '{print $4}' > tmp`basename $0`1.txt
        awk '{if($1==""){print "NoGeneForTFsite"}else{print $0}}' tmp`basename $0`1.txt > tmp`basename $0`2.txt
    elif [ $2 == "ANNO" ]; then
        awk '{print $8}' tmp`basename $0`.txt|tr ";" "\n"|grep ^ANN=|awk -F"|" '{print $2}' > tmp`basename $0`2.txt
    elif [ $2 == "ANNO_IMPACT" ]; then
        awk '{print $8}' tmp`basename $0`.txt|tr ";" "\n"|grep ^ANN=|awk -F"|" '{print $2"_"$3}' > tmp`basename $0`2.txt
    elif [ $2 == "AC" ]; then
        awk '{print $8}' tmp`basename $0`.txt|tr ";" "\n"|grep ^AC= |sed -e 's/AC=//g' > tmp`basename $0`2.txt
    elif [ $2 == "ENST" ]; then
        awk '{print $8}' tmp`basename $0`.txt|tr ";" "\n"|grep ^ANN=|awk -F"|" '{print $7}' > tmp`basename $0`1.txt
        awk '{if($1==""){print "NoENST"}else{print $0}}' tmp`basename $0`1.txt > tmp`basename $0`2.txt
    else
        echo "2ND ARG: GENE, ANNO, ANNO_IMPACT, or ENST" 1>&2
    fi
    
    # body + something + header
    paste tmp`basename $0`.txt tmp`basename $0`2.txt|cat <(head -n $NROW $1) -
}

extract_variant_vcf(){
    # 1st arg: vcf
    # output: sample,chr,pos,ref,alt,enst,anno
    paste_last_for_snpeffvcf $1 ANNO > tmp${FUNCNAME[0]}.vcf
    paste_last_for_snpeffvcf tmp${FUNCNAME[0]}.vcf ENST > tmp${FUNCNAME[0]}2.vcf
    extract_vqslod tmp${FUNCNAME[0]}2.vcf > tmp${FUNCNAME[0]}3.txt
    remove_vcf_header tmp${FUNCNAME[0]}2.vcf > tmp${FUNCNAME[0]}4.vcf
    paste tmp${FUNCNAME[0]}4.vcf tmp${FUNCNAME[0]}3.txt > tmp${FUNCNAME[0]}5.vcf
    cat_vcf_header tmp${FUNCNAME[0]}5.vcf $1 > tmp${FUNCNAME[0]}6.vcf
    local NROW=`head -n1000 tmp${FUNCNAME[0]}6.vcf|awk '$1=="#CHROM"{print NR}'`
    awk -F"\t" 'BEGIN{OFS="\t"}NR=='"${NROW}"'{for(i=10;i<=NF;i++)SAMPLE[i]=$i}NR>'"${NROW}"'{
        for(i=10;i<=(NF-2);i++){if($i ~ /^0\/1/ || $i ~ /^1\/0/ || $i ~ /^1\/1/){print SAMPLE[i],$1,$2,$4,$5,$(NF-1),$(NF-2),$NF}}
    }' tmp${FUNCNAME[0]}6.vcf
}

extract_vqslod(){ #1st arg: vcf; output: status"\t"vqslod score
    remove_vcf_header $1 > tmp${FUNCNAME[0]}.txt
    awk -F"\t" '{
        n=split($8,A,";")
        for(i=1;i<=n;i++){
            if(A[i] ~ /^VQSLOD=/){print A[i]; break}
            if(i == n){print "NoScore"}
        }
    }' tmp${FUNCNAME[0]}.txt > tmp${FUNCNAME[0]}2.txt   # PASS"\t"Score
    sed -e 's/VQSLOD=//g' < tmp${FUNCNAME[0]}2.txt
}

triodenovo(){
    TRIODENOVO="/path/to/triodenovo"
    $TRIODENOVO --ped $2  --minDQ 0.00 --in_vcf $1 --out_vcf tmp`basename $0`.vcf > log`basename $0`.txt 2>&1
    awk '$6==2{print $2}' $2 > tmp`basename $0`Pt.txt
    triodenovo_to_score tmp`basename $0`.vcf tmp`basename $0`Pt.txt
}

triodenovo_to_score(){
    local NROW=`awk '$1=="#CHROM"{print NR}' $1`
    awk 'NR=='"${NROW}"' {
        for (i=10; i<=NF; i++) {
            SAMPLE[i]=$i
        }
    }   NR>'"${NROW}"' {
        for (i=10; i<=NF; i++) {
            if ( $i !~ /^\./ ){
                split($i,ARRAY,":")
                print SAMPLE[i], $1, $2, $4, $5, ARRAY[2]
            }
        }
    }' $1 > tmp`basename $0`.txt
    
    # only for proband
    join <(sort tmp`basename $0`.txt) <(sort $2)|tr ' ' '\t'
}

dnmfilter(){
    awk '{print $1","$2","$3}' $1|sed -e 's/Sample_//g' -e 's/chr//g' -e 's/X/999X/g' -e 's/Y/999Y/g'|\
        sort -t"," -k1,1n -k2,2n -k3,3n|sed -e  's/999X/X/g' -e 's/999Y/Y/g'|awk '{print "Sample_"$0}' > tmp`basename $0`.csv
    
    java -jar DNMFilter.jar gbm \
    --reference hs37d5.fa \
    --pedigree $2 \
    --bam $3 \
    --training Training.264epi4ktrios.csv \
    --candidate tmp`basename $0`.csv \
    --configuration Features.conf \
    --cutoff 0 \
    --output tmp`basename $0`3.csv
    
    awk -F"," 'BEGIN{OFS="\t"}{print $1,$2,$3,"NA","NA",$4}' tmp`basename $0`3.csv
}

af_filter_for_exac(){
    local NROW=`head -n1000 $1|awk '$1=="#CHROM"{print NR}'`
    local MAF=`echo $2`
    
    cut -f8 $1|tr ';' '\n'|grep '^AF=' |sed -e 's/AF=//g' > tmp`basename $0`AF.txt
    head -n $NROW $1 > tmp`basename $0`He.txt
    paste tmp`basename $0`AF.txt <(tail -n+`expr 1 + $NROW` $1) |awk -F"\t" '$1<'"${MAF}"'' |cut -f2-| cat tmp`basename $0`He.txt -
}

for_multi_anno_norm(){
    NROW=$(expr $(head -n1000 $1 |awk '$1 !~ /^#/{print NR}'|head -n1) - 1)
    tail -n+`expr 1 + $NROW` $1 > tmp`basename $0`.txt
    
    if [ $2 == "MULTIPLE" ]; then
        awk -F"\t" '{
            split($8,ARRAY,"|")
            for(i=1;i<=length(ARRAY);i++){
                if(ARRAY[i]=="HIGH" || ARRAY[i]=="MODERATE" || ARRAY[i]=="LOW" || ARRAY[i]=="MODIFIER"){
                    if(1<=(i-1) && (i+4)<=length(ARRAY)){
                        print $0"\t"ARRAY[i-1]"|"ARRAY[i]"|"ARRAY[i+1]"|"ARRAY[i+4]
                    }
                }
            }
        }' tmp`basename $0`.txt > tmp`basename $0`2.txt
    elif [ $2 == "SIMPLE" ]; then
        awk -F"\t" '{
            split($8,ARRAY,"|")
            for(i=1;i<=length(ARRAY);i++){
                if(ARRAY[i]=="HIGH" || ARRAY[i]=="MODERATE" || ARRAY[i]=="LOW" || ARRAY[i]=="MODIFIER"){
                    if(1<=(i-1) && (i+4)<=length(ARRAY)){
                        print $0"\t"ARRAY[i-1]"|"ARRAY[i]"|"ARRAY[i+1]"|"ARRAY[i+4]
                    break
                    }
                }
            }
        }' tmp`basename $0`.txt > tmp`basename $0`2.txt
    fi
    
    cat <(head -n $NROW $1) tmp`basename $0`2.txt
}

for_multi_anno_select(){
    local NROW=`head -n1000 $1 | awk '$1=="#CHROM"{print NR}'`
    tail -n +`expr 1 + $NROW` $1 > tmp`basename $0`.txt
    awk '{print $NF}' tmp`basename $0`.txt|cut -d"|" -f1|paste - tmp`basename $0`.txt > tmp`basename $0`2.txt
    join -t$'\t' <(sort tmp`basename $0`2.txt) <(sort $2)|cut -f2- > tmp`basename $0`3.txt
    cat <(head -n $NROW $1) tmp`basename $0`3.txt
}

vcf_anno(){
    head -n1000 $1 > tmp`basename $0`0.txt
    if grep "#CHROM" tmp`basename $0`0.txt > /dev/null; then
        local NROW=`head -n1000 $1|awk '$1=="#CHROM"{print NR}'`
    else
        local NROW=0
    fi
    
    tail -n+`expr 1 + $NROW` $1 > tmp`basename $0`.txt
    awk '{print $1"_"$2"\t"$0}' tmp`basename $0`.txt > tmp`basename $0`2.txt
    awk '{print $1"_"$2"\t"$0}' $2 > tmp`basename $0`3.txt
    join -a 1 -e $4 -o 2.4 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 -t$'\t' <(sort tmp`basename $0`2.txt) <(sort tmp`basename $0`3.txt) > tmp`basename $0`4.txt
    paste <(cut -f2-$3 tmp`basename $0`4.txt) <(cut -f1 tmp`basename $0`4.txt) <(cut -f`expr 2 + $3`- tmp`basename $0`4.txt) > tmp`basename $0`5.txt
    cat <(head -n $NROW $1) tmp`basename $0`5.txt
}

for_multi_anno_sum_count(){
    local NROW=`head -n1000 $1 | awk '$1=="#CHROM"{print NR}'`
    tail -n +`expr 1 + $NROW` $1 > tmp`basename $0`.txt
    awk '{print $NF}' tmp`basename $0`.txt|cut -d"|" -f4|paste tmp`basename $0`.txt - > tmp`basename $0`2.txt
    awk '{SUM[$NF]++}END{for(i in SUM){print i"\t"SUM[i]}}' tmp`basename $0`2.txt
}

for_multi_anno_sum_count_depth(){
    local NROW=`head -n1000 $1 | awk '$1=="#CHROM"{print NR}'`
    tail -n +`expr 1 + $NROW` $1 > tmp`basename $0`.txt
    
    awk '{print $NF}' tmp`basename $0`.txt|cut -d"|" -f4|paste tmp`basename $0`.txt - > tmp`basename $0`2.txt
    
    if [ $2 == "DEPTH" ]; then
        awk '{SUM[$6]++}END{for(i in SUM){print i"\t"SUM[i]}}' tmp`basename $0`2.txt
    elif [ $2 == "DEPTH_ENST" ]; then
        awk '{SUM[$NF"_"$6]++}END{for(i in SUM){print i, SUM[i]}}' tmp`basename $0`2.txt|awk -F'[_ ]' '{print $1"\t"$2"\t"$3}'
    fi
}

for_multi_anno_sum_mutability_depth(){
    local NROW=`head -n1000 $1 | awk '$1=="#CHROM"{print NR}'`
    tail -n +`expr 1 + $NROW` $1 > tmp`basename $0`.txt
    
    awk '{print $NF}' tmp`basename $0`.txt|cut -d"|" -f4|cut -d":" -f1|paste tmp`basename $0`.txt - > tmp`basename $0`2.txt 
    
    if [ $2 == "SIMPLE" ]; then
        awk '{SUM[$NF]+=$3}END{for(i in SUM){print i"\t"SUM[i]}}' tmp`basename $0`2.txt
    elif [ $2 == "DEPTH2" ]; then
        awk '{SUM[$6]+=$3}END{for(i in SUM){print i"\t"SUM[i]}}' tmp`basename $0`2.txt
    elif [ $2 == "DEPTHADJ" ]; then
        # sum mutability adusted with depth(=$6) for each ENST
        awk '{
            if($6<12.5){
                SUM[$NF]+=$3*($6*0.04+0.125)
            }else if($6<55){
                SUM[$NF]+=$3*(0.177*log($6-6.29)+0.302)
            }else{
                SUM[$NF]+=$3
            }
        }END{for(i in SUM){print i"\t"SUM[i]}}' tmp`basename $0`2.txt
    fi
}
