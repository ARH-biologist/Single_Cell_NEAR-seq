#Written and updated by Alyssa R. Hamm, May 2025

#Making the AaSg Integrated Dataset in Seurat

#open your libraries
library(SeuratObject)
library(Seurat)


#set your working diectory
wd <- "Insert/Working/Directory/Here"
setwd(wd)
getwd()

#############################M1############################################
#######################Prep the data###########################################
M1_Matrix <-read.csv2("M1_QC_Matrix.txt", sep=";",header=TRUE)

rownames(M1_Matrix) <- M1_Matrix[,1] #set the row names
M1_Data <-M1_Matrix[,-c(1)] #remove the repeated column


########################Make Seurat Object######################################

M1 <- CreateSeuratObject(counts = M1_Data, assay="RNA", min.cells = 5, min.features = 15)
M1

M1 <- NormalizeData(M1, normalization.method = "LogNormalize", scale.factor = 10000)
M1 <- FindVariableFeatures(M1, selection.method = "vst", nfeatures = 200)
M1$stim <- "ctrl"

saveRDS(M1, file = "M1_Raw.rds")
#############################M2############################################
#######################Prep the data###########################################
M2_Matrix <-read.csv2("M2_QC_Matrix.txt", sep=";",header=TRUE)

rownames(M2_Matrix) <- M2_Matrix[,1] #set the row names
M2_Data <-M2_Matrix[,-c(1)] #remove the repeated column


########################Make Seurat Object######################################

M2 <- CreateSeuratObject(counts = M2_Data, assay="RNA", min.cells = 5, min.features = 15)
M2

M2 <- NormalizeData(M2, normalization.method = "LogNormalize", scale.factor = 10000)
M2 <- FindVariableFeatures(M2, selection.method = "vst", nfeatures = 200)
M2$stim <- "ctrl"

saveRDS(M2, file = "M2_raw.rds")
#############################M3############################################
#######################Prep the data###########################################
M3_Matrix <-read.csv2("M3_QC_Matrix.txt", sep=";",header=TRUE)

rownames(M3_Matrix) <- M3_Matrix[,1] #set the row names
M3_Data <-M3_Matrix[,-c(1)] #remove the repeated column


########################Make Seurat Object######################################

M3 <- CreateSeuratObject(counts = M3_Data, assay="RNA", min.cells = 5, min.features = 15)
M3

M3 <- NormalizeData(M3, normalization.method = "LogNormalize", scale.factor = 10000)
M3 <- FindVariableFeatures(M3, selection.method = "vst", nfeatures = 200)
M3$stim <- "ctrl"

saveRDS(M3, file = "M3_raw.rds")
#############################C1############################################
#######################Prep the data###########################################
C1_Matrix <-read.csv2("C1_QC_Matrix.txt", sep=";",header=TRUE)

rownames(C1_Matrix) <- C1_Matrix[,1] #set the row names
C1_Data <-C1_Matrix[,-c(1)] #remove the repeated column


########################Make Seurat Object######################################

C1 <- CreateSeuratObject(counts = C1_Data, assay="RNA", min.cells = 5, min.features = 15)
C1

C1 <- NormalizeData(C1, normalization.method = "LogNormalize", scale.factor = 10000)
C1 <- FindVariableFeatures(C1, selection.method = "vst", nfeatures = 200)
C1$stim <- "stim"


saveRDS(C1, file = "C1_raw.rds")
#############################C2############################################
#######################Prep the data###########################################
C2_Matrix <-read.csv2("C2_QC_Matrix.txt", sep=";",header=TRUE)

rownames(C2_Matrix) <- C2_Matrix[,1] #set the row names
C2_Data <-C2_Matrix[,-c(1)] #remove the repeated column


########################Make Seurat Object######################################

C2 <- CreateSeuratObject(counts = C2_Data, assay="RNA", min.cells = 5, min.features = 15)
C2

C2 <- NormalizeData(C2, normalization.method = "LogNormalize", scale.factor = 10000)
C2 <- FindVariableFeatures(C2, selection.method = "vst", nfeatures = 200)
C2$stim <- "stim"

saveRDS(C2, file = "C2_raw.rds")
#############################C3############################################
#######################Prep the data###########################################
C3_Matrix <-read.csv2("C3_QC_Matrix.txt", sep=";",header=TRUE)

rownames(C3_Matrix) <- C3_Matrix[,1] #set the row names
C3_Data <-C3_Matrix[,-c(1)] #remove the repeated column


########################Make Seurat Object######################################

C3 <- CreateSeuratObject(counts = C3_Data, assay="RNA", min.cells = 5, min.features = 15)
C3

C3 <- NormalizeData(C3, normalization.method = "LogNormalize", scale.factor = 10000)
C3 <- FindVariableFeatures(C3, selection.method = "vst", nfeatures = 200)
C3$stim <- "stim"


saveRDS(C3, file = "C3_raw.rds")



############################Run the integration############################
#Create a list of Seurat objects and select features for integration
AaSg_datasets <- list(M1, M2, M3, C1, C2, C3)
AaSg_features <- SelectIntegrationFeatures(object.list = AaSg_datasets, nfeatures = 2000)


#Perform anchor-based integration
AaSg.anchors <- FindIntegrationAnchors(object.list = AaSg_datasets, anchor.features = AaSg_features, scale = TRUE, 
                                       normalization.method = "LogNormalize", reduction = "cca",dims = 1:20,)

###
all_genes <- lapply(AaSg_datasets, row.names) %>% Reduce(intersect, .) # get gene names present in ALL SCTransform'd datasets

AaSgIntegrated <- IntegrateData(anchorset = AaSg.anchors, dims = 1:20, normalization.method = "LogNormalize",
                              features.to.integrate= all_genes)


#Specify the "integrated" assay as default for further analysis
DefaultAssay(AaSgIntegrated) <- "integrated"

Aa <- AaSgIntegrated

##Save the Object
saveRDS(Aa, file = "AaSgIntegrated.rds")
