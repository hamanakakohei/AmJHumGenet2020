TypeOfInterest = c(
"TEGLU7",
"TEGLU8",
"TEGLU10",
"TEGLU3",
"TEGLU20",
"TEGLU2",
"TEINH4",
"TEINH5",
"TEINH6",
"TEINH10",
"TEINH11",
"TEINH12",
"TEINH14",
"TEINH15",
"TEINH16",
"TEINH17",
"TEINH18",
"TEINH19",
"TEINH21",
"CBINH1",
"CBPC",
"CBGRC",
"CBINH2",
"ACBG",
"CBNBL1",
"CBNBL2")

# access to loom file
library(loomR)
lfile <- connect(filename="l5_all.loom",mode="r+")

# check all cell types
unique(lfile$col.attrs$ClusterName[1:160796])
TypeNames = lfile$col.attrs$ClusterName[1:160796]

# get index of cell types of interest
Index = TypeNames %in% TypeOfInterest

# get exp matrix of interest cell types
data <- lfile$matrix[Index,]
row.names(data) <- lfile$col.attrs$CellID[Index]
colnames(data) <- lfile$row.attrs$Gene[1:27998]

# remove duplicate
data = data[unique(row.names(data)),unique(colnames(data))]

# get scaled data
library(Seurat)
library(dplyr)
pbmc <- CreateSeuratObject(raw.data = t(data),min.cells=3,min.genes=1,project="FemaleFetalGonad")
mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
pbmc <- NormalizeData(object=pbmc,normalization.method="LogNormalize",scale.factor=10000)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

# table of cell name and cluster name
Df1 = data.frame(Ord=1:length(pbmc@ident),ID=names(pbmc@ident)) # 
Df2 =data.frame(ID=lfile$col.attrs$CellID[1:160796],orig.ident=lfile$col.attrs$ClusterName[1:160796])
Df2 = Df2[!duplicated(Df2$ID),]
Df3 <- merge(Df1,Df2,by="ID",all.x=T)
Df3 <- Df3[order(Df3$Ord),]

# dot plot
CellIdent = as.factor(Df3$orig.ident)
names(CellIdent) = Df3$ID
CellIdent = gsub("TEGLU7",  "01TEGLU7",  CellIdent)
CellIdent = gsub("TEGLU8",  "02TEGLU8",  CellIdent)
CellIdent = gsub("TEGLU10", "03TEGLU10", CellIdent)
CellIdent = gsub("TEGLU3",  "04TEGLU3",  CellIdent)
CellIdent = gsub("TEGLU20", "05TEGLU20", CellIdent)
CellIdent = gsub("TEGLU2",  "06TEGLU2",  CellIdent)
CellIdent = gsub("TEINH4",  "07TEINH4",  CellIdent)
CellIdent = gsub("TEINH5",  "08TEINH5",  CellIdent)
CellIdent = gsub("TEINH6",  "09TEINH6",  CellIdent)
CellIdent = gsub("TEINH10", "10TEINH10", CellIdent)
CellIdent = gsub("TEINH11", "11TEINH11", CellIdent)
CellIdent = gsub("TEINH12", "12TEINH12", CellIdent)
CellIdent = gsub("TEINH14", "13TEINH14", CellIdent)
CellIdent = gsub("TEINH15", "14TEINH15", CellIdent)
CellIdent = gsub("TEINH16", "15TEINH16", CellIdent)
CellIdent = gsub("TEINH17", "16TEINH17", CellIdent)
CellIdent = gsub("TEINH18", "17TEINH18", CellIdent)
CellIdent = gsub("TEINH19", "18TEINH19", CellIdent)
CellIdent = gsub("TEINH21", "19TEINH21", CellIdent)
CellIdent = gsub("CBINH1",  "20CBINH1",  CellIdent)
CellIdent = gsub("CBPC",    "21CBPC",    CellIdent)
CellIdent = gsub("CBGRC",   "22CBGRC",   CellIdent)
CellIdent = gsub("CBINH2",  "23CBINH2",  CellIdent)
CellIdent = gsub("ACBG",    "24ACBG",    CellIdent)
CellIdent = gsub("CBNBL1",  "25CBNBL1",  CellIdent)
CellIdent = gsub("CBNBL2",  "26CBNBL2",  CellIdent)
pbmc@ident = as.factor(CellIdent)

# sema6b plot
for (x in c(500,600,700,800,900,1000,1100,1200,1300)) {
    png(paste("mousebrain.org.Sema6b190205-",x,".png"),width=500,height=x)
    DotPlot(pbmc,genes=c("Sema6b"),plot.legend=TRUE,dot.scale=20)
    dev.off()
}

# for plotting neurotransmitter
GENES_ARY=c("Slc17a7","Slc17a6","Slc17a8","Slc6a1","Gabbr1","Gabbr2")

# for plotting gene markers
GENES_ARY=c(
"Rasgrf2",
"Pvrl3",
"Cux2",
"Rorb",
"Plcxd2",
"Thsd7a",
"Kcnk2",
"Sulf2",
"Foxp2",
"Syt6",
"Rprm",
"Nr4a2",
"Synpr",
"Pcp4",
"Cplx3")

# for plotting gene markers
GENES_ARY=c(
"Rorb",
"Pvrl3",
"Nrgn",
"Grm2",
"Tfap2b",
"Eomes",
"Fabp7",
"Pvalb",
"Zeb2")

png(paste("xxxxxx.png",sep=""),width=1000,height=500)
DotPlot(object=pbmc, genes.plot=GENES_ARY, plot.legend=TRUE,dot.scale=45)
dev.off()
