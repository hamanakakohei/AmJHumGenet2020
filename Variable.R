#! /usr/bin/Rscript
library(tidyverse)

modify_consequence = function(df,col_name){
    # 1st arg: dataframe
    # 2nd arg: name of col. for variant consequence
    # return dataframe with modified consequence as "mod"
	
    MIS=c("missense_variant"
        ,"missense_variant&splice_region_variant"
        ,"protein_protein_contact"
        ,"rare_amino_acid_variant"
    )
    NON=c("splice_acceptor_variant&intron_variant"
        ,"splice_acceptor_variant&splice_donor_variant&intron_variant"
        ,"splice_acceptor_variant&splice_region_variant&intron_variant"
        ,"splice_donor_variant&intron_variant"
        ,"splice_donor_variant&splice_region_variant&intron_variant"
        ,"start_lost"
        ,"start_lost&splice_region_variant"
        ,"stop_gained"
        ,"stop_gained&splice_region_variant"
        ,"stop_lost"
        ,"stop_lost&splice_region_variant"
        ,"splice_acceptor_variant&splice_donor_variant&splice_region_variant&intron_variant"
    )
    SYN=c("stop_retained_variant"
        ,"synonymous_variant"
        ,"initiator_codon_variant"
        ,"splice_region_variant&initiator_codon_variant"
        ,"splice_region_variant&stop_retained_variant"
        ,"splice_region_variant&synonymous_variant"
    )
    INF=c("disruptive_inframe_deletion"
        ,"disruptive_inframe_deletion&splice_region_variant"
        ,"disruptive_inframe_insertion"
        ,"inframe_deletion"
        ,"inframe_insertion"
    )
    FS=c(
        "frameshift_variant"
        ,"frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant"
        ,"frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant"
        ,"frameshift_variant&splice_region_variant"
        ,"frameshift_variant&start_lost"
        ,"frameshift_variant&stop_gained"
        ,"frameshift_variant&stop_gained&splice_region_variant"
        ,"splice_acceptor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant"
        ,"splice_donor_variant&inframe_deletion&splice_region_variant&intron_variant"
        ,"stop_gained&inframe_insertion"
    )

    n_col=which(colnames(df)==col_name)
    df$mod = apply(
        df[,n_col],1,
        function(x){
            if(x %in% MIS){return("MIS")
            } else if (x %in% NON){return("NON")
            } else if (x %in% SYN){return("SYN")
            } else if (x %in% INF){return("INF")
            } else if (x %in% FS){return("FS")
            } else {return("OTH")
            }
        }
    )
	return(df)	
}
