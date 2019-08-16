#! /usr/bin/Rscript
library(tidyverse)

merge_four_denovo_softs= function(ga,td,dn,df,dg,pc,out){
    # ga: gatk; td: triodenovo; dn: DNMfilter; df: denovofilter
	ga <- read_tsv(ga,col_names=c("sample","chr","pos","ref","alt","gene","anno","vqlsod"),col_types="ccdccccd")
	td <- read_tsv(td,col_names=c("sample","chr","pos","ref","alt","td"),col_types="ccdccc")
	dn <- read_tsv(dn,col_names=c("sample","chr","pos","ref","alt","dn"),col_types="ccdccc")
	df <- read_tsv(df,col_names=c("sample","chr","pos","ref","alt","df"),col_types="ccdccc")
	dg <- read_tsv(dg,col_names=c("sample","chr","pos","ref","alt","dg"),col_types="ccdccc")
	pc <- read_tsv(pc,col_names=c("sample","chr","pos","ref","alt","pc"),col_types="ccdccc")
	 
    ## edit
	#ga$chr <- as.character(ga$chr)
	#ga <- select(ga,sample,chr,pos,ref,alt)
	dn <- select(dn,sample,chr,pos,dn)
    
    # merge
	ga_td 		        <- left_join(ga,            td,by=c("sample","chr","pos","ref","alt"))
	ga_td_dn 	        <- left_join(ga_td,         dn,by=c("sample","chr","pos"))
	ga_td_dn_df	        <- left_join(ga_td_dn,      df,by=c("sample","chr","pos","ref","alt"))
	ga_td_dn_df_dg	    <- left_join(ga_td_dn_df,   dg,by=c("sample","chr","pos","ref","alt"))
	ga_td_dn_df_dg_pc	<- full_join(ga_td_dn_df_dg,pc,by=c("sample","chr","pos","ref","alt"))

	write_tsv(ga_td_dn_df_dg_pc,out)	
}


simplify_denovogear = function(dg_path,output){
    dg <- read_delim(dg_path,delim=" ",col_names=as.character(c(1:45),44),col_types=cols(.default="c")) %>%
        select(seq(3,37,2),40,42,44)
    colnames(dg)=c("CHILD_ID","chr","pos","ref","alt","maxlike_null","pp_null","tgt_null(child/mom/dad)","snpcode","code","maxlike_dnm","pp_dnm","tgt_dnm(child/mom/dad)","lookup","flag","DEPTH_ch","DEPTH_da","DEPTH_mo","MQ_ch","MQ_da","MQ_mo")
    dg %>%
        select(c(CHILD_ID,chr,pos,ref,alt,pp_dnm)) %>%
        write.table(output,row.names=F,col.names=F,quote=F,sep="\t")
}
