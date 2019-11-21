## depth file -> bedgraph file ##


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
