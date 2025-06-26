#Written and updated by Alyssa R. Hamm, June 2025

#Make graphs to visualize data

#open your libraries
library(SeuratObject)
library(Seurat)
library(EnhancedVolcano)
library(RColorBrewer)
library(ggthemes)

#set your working diectory
wd <- "Insert/Working/Directory/Here"
setwd(wd)
getwd()


Aa <-readRDS(file = "AaSgIntegrated_processed.rds")




########base and subpopulation figures

Aa <- SetIdent(Aa, value = "seurat_clusters")


#UMAP of all barcodes
DimPlot(Aa, reduction = "umap")

#UMAP of all barcodes, split between monoculture and coculture
DimPlot(Aa, reduction = "umap", split.by = "stim")

#UMAP of all barcodes, colored based on expression of a single operon
FeaturePlot(Aa, features = "ACT74-operon1559") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#Violin plot of a single operon, split between mono- and coculture, specifically in a single subpopulation
VlnPlot(Aa, assay= "RNA", features = "ACT74-operon1559", pt.size=0.3, split.by= "stim", y.max=8, slot="data", idents= "18")

#Violin plot of a single operon, split between mono- and coculture, shown in all subpopulations
VlnPlot(Aa, assay= "RNA", features = "ACT74-operon1559", pt.size=0, split.by= "stim", y.max=8, slot="data")


########NEAR-seq figures

Aa <- SetIdent(Aa, value = "aggregate")

AggVCo <- read.csv2("NEARseq_Markers.csv")
EnhancedVolcano(AggVCo,
                lab = rownames(AggVCo),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                col=c("grey30", "#00513B", "#0072B2", "#D55E00"),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                labSize = 0)
