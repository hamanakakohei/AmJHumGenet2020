## depth file -> bedgraph file ##
library(tidyverse)
depth_path = "DepthV50417M.txt"
OUT = "v5.bedgraph"
tmp_v4 = read_tsv(depth_path,col_names=c("chr","pos","v4"),col_types="cdd")

## depth filt => bedgraph format ##
v4 = filter(tmp_v4,chr=="chr19" & pos>4500000 & pos<4600000)

sta_ary <- c(v4$pos[1])
end_ary <- c()
val_ary <- c(v4$v4[1])
for(i in 2:nrow(v4)){
    if(v4$pos[i]!=v4$pos[i-1]+1){
        end <- v4$pos[i-1]
        sta <- v4$pos[i] - 1
        val <- v4$v4[i]
        end_ary <- c(end_ary,end,sta)       
        sta_ary <- c(sta_ary,end,sta)       
        val_ary <- c(val_ary,0,v4$v4[i])  
    } else if(v4$pos[i]==v4$pos[i-1]+1 && v4$v4[i]!=v4$v4[i-1]){
        end <- v4$pos[i-1]
        sta <- v4$pos[i] - 1
        val <- v4$v4[i]
        end_ary <- c(end_ary,end)
        sta_ary <- c(sta_ary,sta)
        val_ary <- c(val_ary,v4$v4[i])
    }
}
end_ary <- c(end_ary, v4$pos[nrow(v4)])
df <- data.frame(chr="chr19",start=sta_ary,end=end_ary,value=val_ary)
write_tsv(df,OUT,quote=F)


## gviz plot ##
library(data.table)
library(Gviz)

grtrack = readRDS('grtrack.gencode.rds')
dt <- fread('v5.bedGraph',col.names = c('chromosome', 'start', 'end', 'value'))
st <- 4542000
en <- 4564000

dt=mutate(dt,chromosome='chr19')

dtrack <- DataTrack(
    range = dt,
    type = "histogram",
    genome = 'hg19',
    name = "Depth",
    chr = "chr19"
)

gtrack <- GenomeAxisTrack()

plotTracks(
    list(gtrack,dtrack,grtrack),
     chromosome='chr19', from = st, to = en
 )
