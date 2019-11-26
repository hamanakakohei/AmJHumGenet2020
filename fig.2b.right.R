library("R.matlab")

data <- readMat("GSM3017261_150000_CNS_nuclei.mat")
# DGE: 156049 * 26894
# sample.type: 156049; p2_brain, p11_brain, p2_spine, p11_spine
# spinal.cluster.assignment: 156049
# genes: 26894
# cluster.assignment: 156049
# barcodes: 156049

# cluster.assignment
"53 Unresolved                  "
"44 Migrating Int Lhx6          "
"25 CB Granule Precursor        "
"2 OB Mitral/Tufted Ms4a15      "
"17 CTX PyrL6                   "
"28 CB Granule                  "
"69 Astro Prdm16                "
"58 Oligo NFOL1                 "
"7 CTX PyrL2/L3 Met             "
"47 Migrating Int Foxp2         "
"61 OPC                         "
"4 Medium Spiny Neurons         "
"9 CTX PyrL2/L3/L4 Mef2c        "
"46 Migrating Int Cpa6          "
"11 CTX PyrL4/L5                "
"20 THAL Glut                   "
"37 HIPP Pyr Precursor          "
"64 Endothelia                  "
"14 CTX PyrL6a                  "
"22 Purkinje Early              "
"35 HIPP Pyr Crym               "
"48 Migrating Int Pbx3          "
"23 Purkinje Late               "
"6 CTX PyrL2/L3/L4 Ntf3         "
"18 CLAU Pyr                    "
"68 Astro Slc7a10               "
"26 CB Int Stellate/Basket      "
"30 MD Glyc Int                 "
"57 Oligo MOL                   "
"16 CTX PyrL5/L6 Npr3           "
"38 HIPP Pyr Grik4              "
"13 CTX PyrL5 Fezf2             "
"40 HIPP Granule/PyrCA3         "
"60 Oligo COP2                  "
"21 THAL Int Six3               "
"45 Migrating Int Trdn          "
"51 SVZ Stem                    "
"49 Migrating Int Lgr6          "
"1 OB Mitral/Tufted Eomes       "
"31 MD Int Rxfp2                "
"8 CTX PyrL4 Wnt5b              "
"39 HIPP Granule Nrp2           "
"10 CTX PyrL4 Rorb              "
"15 CTX PyrL5/L6 Sulf1          "
"55 Oligo MFOL2                 "
"19 MTt Glut                    "
"27 CB Int Golgi/Stellate/Basket"
"12 CTX PyrL5 Itgb3             "
"54 Unresolved Kcng1            "
"65 SMC                         "
"71 Bergmann Glia               "
"43 SC Glut Gna14               "
"29 CB Int Precursor            "
"36 HIPP Granule Mki67          "
"56 Oligo MFOL1                 "
"66 VLMC Slc6a13                "
"24 CB Int Progenitor           "
"59 Oligo COP1                  "
"67 VLMC Slc47a1                "
"70 Astro Gfap                  "
"72 Ependyma                    "
"63 Microglia                   "
"73 OEC                         "
"34 SUB Pyr                     "
"62 Macrophage                  "
"5 CTX PyrL2/L3 Pappa2          "
"42 SC Glut Hmga2               "
"33 HIPP Pyr Cr2                "


medium spiny neuron: striatum
OB: olfactory bulb
CLAU: claustrum
SVZ: subventricular zone
SUB: subiculum
SMC: smooth muscle cell
VLMC: vascula and leptomeningeal cell
OEC: olbfactory ensheathling cell
MTt: Midbrain (Rostral?
SC: superior colliculus?
MD: Medulla

# cell type of interest
TypeOfInterest=c(
"5 CTX PyrL2/L3 Pappa2          ",
"7 CTX PyrL2/L3 Met             ",
"9 CTX PyrL2/L3/L4 Mef2c        ",
"6 CTX PyrL2/L3/L4 Ntf3         ",
"8 CTX PyrL4 Wnt5b              ",
"10 CTX PyrL4 Rorb              ",
"11 CTX PyrL4/L5                ",
"12 CTX PyrL5 Itgb3             ",
"13 CTX PyrL5 Fezf2             ",
"15 CTX PyrL5/L6 Sulf1          ",
"16 CTX PyrL5/L6 Npr3           ",
"17 CTX PyrL6                   ",
"14 CTX PyrL6a                  ",
"44 Migrating Int Lhx6          ",
"45 Migrating Int Trdn          ",
"46 Migrating Int Cpa6          ",
"47 Migrating Int Foxp2         ",
"48 Migrating Int Pbx3          ",
"49 Migrating Int Lgr6          ",
"22 Purkinje Early              ",
"23 Purkinje Late               ",
"25 CB Granule Precursor        ",
"28 CB Granule                  ",
"24 CB Int Progenitor           ",
"29 CB Int Precursor            ",
"26 CB Int Stellate/Basket      ",
"27 CB Int Golgi/Stellate/Basket",
"71 Bergmann Glia               "
)

# matrix of interest types
dge <- t(data$DGE)
row.names(dge) <- data$genes
colnames(dge)  <- data$barcodes
TypeNames = data$cluster.assignment
TFArray = TypeNames %in% TypeOfInterest
dge = dge[,TFArray]

# expression calcu.
library(Seurat)
library(dplyr)
pbmc <- CreateSeuratObject(raw.data=dge, min.cells=1, min.genes=1, project="10X_PBMC")
pbmc <- NormalizeData(object=pbmc, normalization.method="LogNormalize",scale.factor=10000)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))

# edit pbmc@ident
Ident = data$cluster.assignment[TFArray]
names(Ident) = names(pbmc@ident)
Ident = gsub("5 CTX PyrL2/L3 Pappa2          ", "01CTX PyrL2/L3 Pappa2          ", Ident)
Ident = gsub("7 CTX PyrL2/L3 Met             ", "02CTX PyrL2/L3 Met             ", Ident)
Ident = gsub("9 CTX PyrL2/L3/L4 Mef2c        ", "03CTX PyrL2/L3/L4 Mef2c        ", Ident)
Ident = gsub("6 CTX PyrL2/L3/L4 Ntf3         ", "04CTX PyrL2/L3/L4 Ntf3         ", Ident)
Ident = gsub("8 CTX PyrL4 Wnt5b              ", "05CTX PyrL4 Wnt5b              ", Ident)
Ident = gsub("10 CTX PyrL4 Rorb              ", "06CTX PyrL4 Rorb              ", Ident)
Ident = gsub("11 CTX PyrL4/L5                ", "07CTX PyrL4/L5                ", Ident)
Ident = gsub("12 CTX PyrL5 Itgb3             ", "08CTX PyrL5 Itgb3             ", Ident)
Ident = gsub("13 CTX PyrL5 Fezf2             ", "09CTX PyrL5 Fezf2             ", Ident)
Ident = gsub("15 CTX PyrL5/L6 Sulf1          ", "10CTX PyrL5/L6 Sulf1          ", Ident)
Ident = gsub("16 CTX PyrL5/L6 Npr3           ", "11CTX PyrL5/L6 Npr3           ", Ident)
Ident = gsub("17 CTX PyrL6                   ", "12CTX PyrL6                   ", Ident)
Ident = gsub("14 CTX PyrL6a                  ", "13CTX PyrL6a                  ", Ident)
Ident = gsub("44 Migrating Int Lhx6          ", "14Migrating Int Lhx6          ", Ident)
Ident = gsub("45 Migrating Int Trdn          ", "15Migrating Int Trdn          ", Ident)
Ident = gsub("46 Migrating Int Cpa6          ", "16Migrating Int Cpa6          ", Ident)
Ident = gsub("47 Migrating Int Foxp2         ", "17Migrating Int Foxp2         ", Ident)
Ident = gsub("48 Migrating Int Pbx3          ", "18Migrating Int Pbx3          ", Ident)
Ident = gsub("49 Migrating Int Lgr6          ", "19Migrating Int Lgr6          ", Ident)
Ident = gsub("22 Purkinje Early              ", "20Purkinje Early              ", Ident)
Ident = gsub("23 Purkinje Late               ", "21Purkinje Late               ", Ident)
Ident = gsub("25 CB Granule Precursor        ", "22CB Granule Precursor        ", Ident)
Ident = gsub("28 CB Granule                  ", "23CB Granule                  ", Ident)
Ident = gsub("24 CB Int Progenitor           ", "24CB Int Progenitor           ", Ident)
Ident = gsub("29 CB Int Precursor            ", "25CB Int Precursor            ", Ident)
Ident = gsub("26 CB Int Stellate/Basket      ", "26CB Int Stellate/Basket      ", Ident)
Ident = gsub("27 CB Int Golgi/Stellate/Basket", "27CB Int Golgi/Stellate/Basket", Ident)
Ident = gsub("71 Bergmann Glia               ", "28Bergmann Glia               ", Ident)
pbmc@ident = as.factor(Ident)

## dotplot
#for (x in c(500,600,700,800,900,1000,1100,1200,1300)) {
#png(paste("ScienceSema6b190205-2-",x,".png"),width=500,height=x)
#DotPlot(pbmc,genes=c("Sema6b        "),plot.legend=TRUE,dot.scale=20)
#dev.off()
#}

# check neurotransmitter 
# GENES_ARY=c("Slc17a7       ","Slc17a6       ","Slc17a8       ","Slc6a1        ","Gabbr1        ","Gabbr2        ")
GENES_ARY=c(
#"Pcp2          ",
"Car8          ",
"Cbln3         ",
"Gabra6        ",
"Pvalb         ",
"Tfap2b        ",
"Grm2          ",
"Nrgn          ",
"Fabp7         ")
#"Zeb2          ")
png(paste("Science-TableS2.png",sep=""),width=1000,height=500)
DotPlot(object=pbmc, genes.plot=GENES_ARY, plot.legend=TRUE,dot.scale=20)
dev.off()
