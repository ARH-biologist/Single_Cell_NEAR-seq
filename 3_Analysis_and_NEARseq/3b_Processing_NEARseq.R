#Written and updated by Alyssa R. Hamm, June 2025

#Run data analysis and NEAR-seq

#open your libraries
library(SeuratObject)
library(Seurat)

#set your working diectory
wd <- "Insert/Working/Directory/Here"
setwd(wd)
getwd()



Aa <-readRDS(file = "AaSgIntegrated.rds")


##########################Process Data#######################################

Aa <- ScaleData(Aa)
Aa <- RunPCA(Aa, npcs = 20)
Aa <- RunUMAP(Aa, reduction = "pca", dims = 1:20)
Aa <- FindNeighbors(Aa, reduction = "pca", dims = 1:20, k.param = 41, prune.SNN = 0.06)
Aa <- FindClusters(Aa, resolution = 2.2)


AllMarkers <-FindAllMarkers(Aa, assay="integrated", min.pct=0, logfc.threshold=0.000000001)
write.table(AllMarkers, "AaSgIntegrated_AllMarkers.txt")


##############################Add all metadata############################
#Need to have a document with barcode names in the first column and your data 
#in another column. Here, I show the NEAR-seq aggregation information, but you
#should run this code for all attributes.
metadata <- read.csv("MetaData.csv", header=TRUE)


#Add this new column into the metadata
Aa <- AddMetaData(Aa, metadata, col.name = 'Aggregates')
#metadata is the name of the csv file
#aggregate is what the feature will be called
#col.name lines up with the name of the column in the dataset

Aa <- SetIdent(Aa, value = "aggregate")

AggVCo <- FindMarkers(Aa, assay= "integrated", ident.1= "Aggregate", ident.2= "Co", min.pct=-Inf, logfc.threshold=-Inf, min.cells.feature = -Inf,
                      min.cells.group = -Inf)

write.csv(AggVCo, "NEARseq_Markers.csv")


##Save the Object
saveRDS(Aa, file = "AaSgIntegrated_processed.rds")
