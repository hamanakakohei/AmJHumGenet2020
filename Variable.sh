prep_seq_file(){
    # from ">" line, extract necessary info
    awk '$1 ~ /^>/{print $2,$1,$5}' $1|sed -e 's/range=//g'|\
        awk '{split($1,col1,":");split(col1[2],col1_pos,"-");split($2,col2,"_");print col1[1],col1_pos[1],col1_pos[2],col2[3],$3}' > tmphead.txt

    # connect pieces to one sequence
    awk 'NR>1{if($1 ~ /^>/){print ">"}else{print $0}}END{printf ">"}' $1|tr '\n' ' '|sed -e 's/ //g'|tr '>' '\n' > tmpseq.txt

    paste tmphead.txt tmpseq.txt > tmpheadseq.txt

    # complement seq
    cat tmpheadseq.txt|while read LINE; do
        if [ `echo $LINE|awk '{print $5}'` = "strand=+" ]; then
            echo $LINE|awk '{print $6}' >> tmpcmp.txt
        else
            echo $LINE|awk '{print $6}' |rev|tr acgt tgca >> tmpcmp.txt
        fi
    done
    
    awk 'BEGIN{OFS="\t"}{print $4,$1,$2,$3}' tmpheadseq.txt|paste - tmpcmp.txt
}
