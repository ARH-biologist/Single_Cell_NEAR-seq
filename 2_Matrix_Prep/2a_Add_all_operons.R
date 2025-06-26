#Written and updated by Alyssa R. Hamm, May 2025

#Add all operons to gene matrix
#Your dataset needs to include a column for every single operon possible



#set your working diectory
wd <- "Insert/Working/Directory/Here"
setwd(wd)
getwd()


#import libraries
library(dplyr)
library(tidyverse)


#####Open the csv file that has all of the operons listed on the top
#Use Excel, the top row needs operon names starting in cell A1
#save as operons_row.csv in your working directory
operons <-read.csv2("operons_row.csv", sep ="," , header=TRUE)
operons <- operons %>% 
  as_tibble()

#####Open your specific gene matrix 
#Here, our sample name is M1 and our species is abbreviated Aa
Aa <-read.csv2("M1_gene_matrix.csv", sep ="," , header=TRUE)
Aa <- Aa %>% 
  as_tibble()

#####Save your barcode names
#Open you gene matrix in Excel
#Copy and paste the barcode names into a new Excel sheet 
#Save as M1_barcodes.csv in your working directory


Aa <-Aa[,-c(1)] #remove the barcodes for now
Aa <- Aa %>% 
  as_tibble()

#Combine the two datasets
New <-full_join(operons,Aa)


#Since you made new columns, you may have random NAs
#Use this line to replace any NAs with 0s
New <- New %>% replace(is.na(.), 0)


#Add the barcode identities back into the dataframe
Barcodes <-read.table("M1_Barcodes.csv", sep ="," , header=TRUE)
Barcodes <- Barcodes %>% 
  as_tibble()


GeneMatrix <- New %>%
  add_column(Barcodes, .before=("SGO_operon1")) ###SGO-operon1 is the name of our first operon in the csv
View(GeneMatrix)

#Save the gene matrix
write.csv(GeneMatrix, "M1_Gene_Matrix_Alloperons.csv")

###Move on to the script 2b_Trimming.R

