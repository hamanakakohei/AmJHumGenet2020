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

sequence_to_all_variant2  = function(){
	options(scipen=100)
	print("WARNING: 'chr' col must contain 'chr'")
	SEQ = read_tsv(commandArgs(trailingOnly=TRUE)[1],col_names=c("enst","chr","start","end","seq")) %>% mutate(seq=tolower(seq)) %>% mutate(len=nchar(seq))
	OUTPUT <- commandArgs(trailingOnly=TRUE)[2]
	MU_SNP <- read.table("mut_rate.txt", header=TRUE)
	
	n_variant = 3*sum(SEQ$len-2)
	variant_list <- vector("list",n_variant)
	POS_LIST=0
	for(NROW in 1:nrow(SEQ)){
	    TMP = SEQ[NROW,]
	    CHR = as.character(TMP$chr)
	    for(NCHAR in 2:(TMP$len-1)){
	        POS = TMP$start + NCHAR - 1
	        CONTEXT = substr(TMP$seq,NCHAR-1,NCHAR+1)
	        ALT_CONTEXT = CONTEXT
	        POS_LIST=POS_LIST+3
	        if(substring(CONTEXT,2,2)=="a"){
	            substring(ALT_CONTEXT,2,2) <- "t"
	            variant_list[[POS_LIST-2]] = data.frame(chr=CHR,pos=POS,context=CONTEXT,alt_context=ALT_CONTEXT)
	            substring(ALT_CONTEXT,2,2) <- "g"
	            variant_list[[POS_LIST-1]] = data.frame(chr=CHR,pos=POS,context=CONTEXT,alt_context=ALT_CONTEXT)
	            substring(ALT_CONTEXT,2,2) <- "c"
	            variant_list[[POS_LIST]] = data.frame(chr=CHR,pos=POS,context=CONTEXT,alt_context=ALT_CONTEXT)
	            next
	        } else if (substring(CONTEXT,2,2)=="t"){
	            substring(ALT_CONTEXT,2,2) <- "a"
	            variant_list[[POS_LIST-2]] = data.frame(chr=CHR,pos=POS,context=CONTEXT,alt_context=ALT_CONTEXT)
	            substring(ALT_CONTEXT,2,2) <- "g"
	            variant_list[[POS_LIST-1]] = data.frame(chr=CHR,pos=POS,context=CONTEXT,alt_context=ALT_CONTEXT)
	            substring(ALT_CONTEXT,2,2) <- "c"
	            variant_list[[POS_LIST]] = data.frame(chr=CHR,pos=POS,context=CONTEXT,alt_context=ALT_CONTEXT)
	            next
	        } else if (substring(CONTEXT,2,2)=="g"){
	            substring(ALT_CONTEXT,2,2) <- "a"
	            variant_list[[POS_LIST-2]] = data.frame(chr=CHR,pos=POS,context=CONTEXT,alt_context=ALT_CONTEXT)
	            substring(ALT_CONTEXT,2,2) <- "t"
	            variant_list[[POS_LIST-1]] = data.frame(chr=CHR,pos=POS,context=CONTEXT,alt_context=ALT_CONTEXT)
	            substring(ALT_CONTEXT,2,2) <- "c"
	            variant_list[[POS_LIST]] = data.frame(chr=CHR,pos=POS,context=CONTEXT,alt_context=ALT_CONTEXT)
	            next
	        } else if (substring(CONTEXT,2,2)=="c"){
	            substring(ALT_CONTEXT,2,2) <- "a"
	            variant_list[[POS_LIST-2]] = data.frame(chr=CHR,pos=POS,context=CONTEXT,alt_context=ALT_CONTEXT)
	            substring(ALT_CONTEXT,2,2) <- "t"
	            variant_list[[POS_LIST-1]] = data.frame(chr=CHR,pos=POS,context=CONTEXT,alt_context=ALT_CONTEXT)
	            substring(ALT_CONTEXT,2,2) <- "g"
	            variant_list[[POS_LIST]] = data.frame(chr=CHR,pos=POS,context=CONTEXT,alt_context=ALT_CONTEXT)
	        } else {
	            print("ERROR:not any of a,t,g,c,,,,")
	            variant_list[[POS_LIST-2]] = data.frame(chr=NA,pos=NA,context=NA,alt_context=NA)
	            variant_list[[POS_LIST-1]] = data.frame(chr=NA,pos=NA,context=NA,alt_context=NA)
	            variant_list[[POS_LIST]] = data.frame(chr=NA,pos=NA,context=NA,alt_context=NA)
	        }
	    }
	
	    if(NROW==100){
	        print(paste(Sys.time(),"100finished"))
	    }else if(NROW==1000){
	        print(paste(Sys.time(),"1000finished"))
	    }else if(NROW==5000){
	        print(paste(Sys.time(),"5000finished"))
	    }
	}
	variant = bind_rows(variant_list)
	variant$context <- toupper(variant$context)
	variant$alt_context <- toupper(variant$alt_context)
	
	# mutability add, order, ref, alt, "." column add
	variant = merge(variant, MU_SNP, by.x = c("context", "alt_context"), by.y = c("from", "to"))
	variant <- variant[,c(3,4,5,1,2)]
	variant$context <- substring(variant$context,2,2)
	variant$alt_context <- substring(variant$alt_context,2,2)
	variant <- data.frame(variant,Ad1=".",Ad2=".",Ad3=".",Ad4=".")
	variant = variant[order(variant$chr,variant$pos),]
	write.table(variant,OUTPUT,sep="\t",quote=F,col.names=F,row.names=F,append=F)
}
