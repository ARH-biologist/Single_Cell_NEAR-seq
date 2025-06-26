#Written and updated by Alyssa R. Hamm, May 2025

#Now you can do some manual trimming of the data to prepare for Seurat
#Here you change some formatting and remove rRNA counts if wanted


#set your working diectory
wd <- "Insert/Working/Directory/Here"
setwd(wd)
getwd()

#import libraries
library(tidyverse)
library(dplyr)
library(naniar)
library(janitor)

#import and name operon matrix
M1 <-read.table("M1_Gene_Matrix_Alloperons.csv", sep ="," , header=TRUE)

#Make a data frame
Matrix <-data.frame(M1)

#####Remove cell barcodes for now
#Make sure you have your barcode names saved from the last script
Matrix <-Matrix[,-c(1)] 

#Ensure all number values are numeric and not characters
Matrix[Matrix == '1'] <- 1
Matrix[Matrix == '2'] <- 2
Matrix[Matrix == '3'] <- 3
Matrix[Matrix == '4'] <- 4
Matrix[Matrix == '5'] <- 5
Matrix[Matrix == '6'] <- 6
Matrix[Matrix == '7'] <- 7
Matrix[Matrix == '8'] <- 8
Matrix[Matrix == '9'] <- 9
Matrix[Matrix == '0'] <- 0
Matrix1 <-Matrix #save your work

#Make your data frame into a tibble so we can alter it
Matrix2 <- Matrix1 %>% 
  as_tibble()

#########Zero out rRNA from the table (unless you want to keep it)###########
#####This is a list of operons known to code for rRNA.
#We also made two extra columns at the end of the .csv file with the sum of
#Aa and Sg rRNA reads. They are also zeroed out here
rRNA <- Matrix2[,c("SGO_operon308","SGO_operon310","SGO_operon311","SGO_operon833",
                   "SGO_operon834","SGO_operon836","SGO_operon949","SGO_operon950","SGO_operon952",
                   "SGO_operon1024","SGO_operon1025","SGO_operon1027","ACT74_operon1204",
                   "ACT74_operon1459","ACT74_operon1499","ACT74_operon1820","ACT74_operon1976",
                   "ACT74_operon2047", "Sg_rRNA", "Aa_rRNA")]

mRNA <-Matrix2 %>%
  replace_with_na(Matrix2,replace = rRNA)

######Since you replaced rRNA reads with NAs,
#Use this line to replace any NAs with 0s
mRNA <- mRNA %>% replace(is.na(.), 0)

#####Remove SGO values
#In our analysis, we only want to look at Aa reads, so zero out Sg
SGO <- mRNA[,c(1:1017)]
Aa <-mRNA %>%
  replace_with_na(mRNA,replace = SGO)
Aa <- Aa %>% replace(is.na(.), 0)

#Add the barcode identities back into the dataframe
Barcodes <-read.table("M1_Barcodes.csv", sep ="," , header=TRUE)
Barcodes <- Barcodes %>% 
  as_tibble()

Matrix8 <- Aa %>%
  add_column(Barcodes, .before=("SGO_operon1"))

#####Transpose
#Subsequent steps need barcodes on top and operons on the side

Matrix9 <-t(Matrix8)

Matrix10 <-Matrix9

Matrix10 <-data.frame(Matrix10)

Matrix10[Matrix10 == '1'] <- 1
Matrix10[Matrix10 == '2'] <- 2
Matrix10[Matrix10 == '3'] <- 3
Matrix10[Matrix10 == '4'] <- 4
Matrix10[Matrix10 == '5'] <- 5
Matrix10[Matrix10 == '6'] <- 6
Matrix10[Matrix10 == '7'] <- 7
Matrix10[Matrix10 == '8'] <- 8
Matrix10[Matrix10 == '9'] <- 9
Matrix10[Matrix10 == '0'] <- 0
Matrix11 <-Matrix10 #save your work

Matrix11 <- Matrix11 %>% 
  as_tibble()
Matrix11 <- Matrix11 %>% row_to_names(row_number = 1)


Matrix11 <- mutate_all(Matrix11, function(x) as.numeric(as.character(x)))

#remove barcodes with over 5000 reads
#before doing this, ensure you have barcode names at the top

Matrix12 <- Matrix11[,colSums(Matrix11) < 5001]



#####Open a csv file that has all of the operons listed on the left side
#Use Excel, the first column needs operon names starting in cell A1
#save as operons_col.csv in your working directory
operons <-read.csv2("operons_col.csv", sep ="," , header=FALSE)
operons <- operons %>% 
  as_tibble()

Matrix13 <- Matrix12 %>%
  add_column(operons, .before=("M1_bc1_10_bc2_10_bc3_20")) 


#Transpose to get into Seurat format
Matrix14 <-t(Matrix13)

Matrix15 <- Matrix14 %>% row_to_names(row_number = 1)

Matrix16 <-t(Matrix15)

#Save file
write.csv2(Matrix16, "M1_QC_Matrix.txt")


###Move on to the next R script,


