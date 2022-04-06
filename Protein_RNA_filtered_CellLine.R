####### Filter datasets to get only data from cell lines with RNA, Protein and arm call data ####
### 210201
## Author: Klaske M. Schukken
## Title: "Filter" data: get cell lines with both RNA and Protein Data only. 371 cells

library('ggplot2')
library('tidyverse')
library('xlsx')
library(readxl)
library(reshape2)
library('BBmisc')

## Step 1: get Data
## Step 2: Get cell names from RNA and Protein Data. Merge. 
## Step 3: Merge cell data with RNA and Protein Dataframes-- 
# getting data only for cells with both RNA and Protein Data. Save data frames. 


####### Step 1: Get data.  ####
### Step 1: Get data. 
## Get  Protein expression and RNA expression Data. 
## Protein expression data
##  https://www.cell.com/cell/fulltext/S0092-8674(19)31385-6#secsectitle0190
## Nusinow et al. Cell, 2020, Quantitative proteomics of the cancer cell line encyclopedia
## and chromosome arm data : Cohen-Sharir et al. (2021) Nature paper 

## go to depmap.org and download the protein expression and RNA expression files. 
## and the cell line data files 


## Import Protein expression data and  chromosome arm data
###
setwd()#set working directory

Protein_Expression.1<-read_excel("mmc2.xlsx", 
                                 sheet= "Normalized Protein Expression")
Protein_Expression.1<-Protein_Expression.1[,c(1:426)] #delete empty collumns
#Name dataframe with protein_ID and corresponding Protein name, this way we have both gene ID and Protein name
Protein_ProID<-Protein_Expression.1[,c(1:2, 6)] 
Protein_ProID$Protein_Id<- str_replace_all(Protein_ProID$Protein_Id, "[|]", ".")
Protein_ProID$Protein_Id<- str_replace_all(Protein_ProID$Protein_Id, "[-]", ".")

write.csv(Protein_ProID, file="Protein_ID_info.csv")

#
#setwd()
Protein_Expression<-read.delim2("Protein_Effect_ID_numeric.csv", 
                                dec=".", header = TRUE, sep=",")
# Protein expression data (no aneuploidy data). 
## Proteins
Protein_Expression<-Protein_Expression[-c(346,90,29),] #non-unique values for 3 row names. repeat cell lines removed.
Protein_Expression_norm <- as.data.frame(Protein_Expression[,-c(1:3)]) #make new data frame
rownames(Protein_Expression_norm) <- Protein_Expression[,3] #add cell line names (rows)
Protein_Expression_norm2<-Protein_Expression_norm
colnames(Protein_Expression_norm2) <- colnames(Protein_Expression_norm) #add gene names
Protein_Expression_norm2$Cell_line <- rownames(Protein_Expression_norm2) #add cell names back in to data frame

# 
#setwd()
aneuploid<- read.csv("arm_calls_data_bendavid.csv", header=TRUE)


#Protein Info. names, chrm location, etc. 
#setwd()
Protein_Info<-read_excel("HGNC_Protein_Info_edited.xlsx", sheet ="HGNC names")
Protein_Info2<-Protein_Info[,c(2,3,11,12,22)]

## Import RNA expression data (no aneuploidy data). CCLE depmap data.
# 1 303 cell lines
# 19 145 RNAs
setwd("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/Quantile.Norm/210121_NoNorm/RNA/")
RNA_Expression_norm2<- read.csv("RNA.Expression.NoNorm.csv", header=TRUE)

# RNA Info. names, chrm location, etc. 
# setwd()
RNA_Info<-read_excel("HGNC_Protein_Info_edited.xlsx", sheet ="HGNC names")
RNA_Info2<-RNA_Info[,c(2,3,11,12,22,14)]




####### Step 2: Get cell names from RNA and Protein Data. Merge.  #####
## Step 2: Get cell names from RNA and Protein Data. Merge. 

RNA_Protein_Merged<-merge(y= Protein_Expression_norm2, x= RNA_Expression_norm2, 
                    by.y="Cell_line", by.x="Cell_line", 
                    sort = TRUE)# 371 cells, 31901  RNA/Proteins
RNA_Protein_CellLines<- RNA_Protein_Merged[,1:2]
RNA_Protein_CellLines<- subset(RNA_Protein_CellLines, Cell_line %in% aneuploid$DepMap_ID)
length(RNA_Protein_CellLines$Cell_line) #cell lines= 367

# setwd()
write.csv(RNA_Protein_CellLines, file= "RNA_Protein_Shared_cells.csv")
RNA_Protein_CellLines<- read_delim(file="RNA_Protein_Shared_cells.csv", 
                                   delim=";")

####### Step 3: Merge cell data with RNA and Protein Dataframes #####
## Step 3: Merge cell data with RNA and Protein Dataframes-- 
# getting data only for cells with both RNA and Protein Data. Save data frames. 

# Protein Data: 
Protein.Expression.filtered<-merge(y= Protein_Expression_norm2, x= RNA_Protein_CellLines, 
                             by.y="Cell_line", by.x="Cell_line", 
                             sort = TRUE)# 371 cells, 12756 Protein 
#Protein.Expression.filtered<- Protein.Expression.filtered[ , -which(names(Protein.Expression.filtered) %in% c("X"))]
length(Protein.Expression.filtered$Cell_line)
write.csv(Protein.Expression.filtered, file= "Protein_Expression_filtered.csv")

# RNA Data: 
RNA.Expression.filtered<-merge(y= RNA_Expression_norm2, x= RNA_Protein_CellLines, 
                                   by.y="Cell_line", by.x="Cell_line", 
                                   sort = TRUE)# 371 cells, 19145 RNA 
#RNA.Expression.filtered<- RNA.Expression.filtered[ , -which(names(RNA.Expression.filtered) %in% c("X.x","X.y"))]
length(RNA.Expression.filtered)
write.csv(RNA.Expression.filtered, file= "RNA_Expression_filtered.csv")

####### Step 4: Make a list of how many cells are in each chomosome arm category (gain/neutral/loss) #####
## Step 4: make a list of how many cells are in each category
# category: per chromosome arm: -1, 0, 1

RNA_Protein_CellLines$Cell_line

RNA_Protein_Cell_Info<-merge(y= RNA_Protein_CellLines, x= aneuploid, 
                         by.y="Cell_line", by.x="DepMap_ID", 
                         sort = TRUE, all=FALSE) 

AllChrmArms<- data.frame(Chrm=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,14,15,16,16,17,17,18,18,19,19,20,20,21,22), 
                         Arm=c("p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q",
                               "q","q","q","p","q","p","q","p","q","p","q","p","q","q","q"))
Cells.per.aneuploid.arm<-data.frame(Chrm.Num.Arm=as.character(), Num.cell.Mono=as.numeric(), Num.cell.Di=as.numeric(), Num.cell.Tri=as.numeric())

for (i in 1:length(AllChrmArms$Chrm)) {
  TestChrm<-AllChrmArms$Chrm[i]
  TestArm<-AllChrmArms$Arm[i]
  
  #subset 1: Find cells with gain or no gain of chromosome X
  testchrm.percell <- filter(RNA_Protein_Cell_Info, chrom==TestChrm) 
  testchrm.percell <- filter(testchrm.percell, arm==TestArm)
  
  Cells.Tri.forchrm<- filter(testchrm.percell, arm_call==1) #all cells triploid for testProt chrm
  Cells.Di.forchrm<- filter(testchrm.percell, arm_call==0) #all cells diploid for testProt chrm
  Cells.Mono.forchrm<- filter(testchrm.percell, arm_call==-1) #all cells diploid for testProt chrm
  
  Cells.per.aneuploid.arm <- rbind(Cells.per.aneuploid.arm, 
                                   data.frame(Chrm.Num.Arm= paste(TestChrm, TestArm),
                                              Num.cell.Mono=length(Cells.Mono.forchrm$DepMap_ID),
                                              Num.cell.Di=length(Cells.Di.forchrm$DepMap_ID),
                                              Num.cell.Tri= length(Cells.Tri.forchrm$DepMap_ID)))
  
}

write.csv(Cells.per.aneuploid.arm, file= "RNA_Protein_filtered_cells_perAneuploidArm.csv")

####### Get only ovarian/cancer_type cancer cells #####

RNA_Protein_CellLines <- read_delim(file="RNA_Protein_Shared_cells.csv", 
                                   delim=";")

Cell_Line_Info <- read_delim(file="sample_info.csv", 
                                   delim=",")
##Ovarian: 
Ovarian_cells <- subset(Cell_Line_Info, lineage=="ovary")

RNA_Protein_Ovary<-subset(RNA_Protein_CellLines, Cell_line %in% Ovarian_cells$DepMap_ID)

## lung=273, blood=133, skin=116, lymphocyte=109, CNS=108: 
Lung_cells <- subset(Cell_Line_Info, lineage_subtype=="NSCLC")
RNA_Protein_Lung<-subset(RNA_Protein_CellLines, Cell_line %in% Lung_cells$DepMap_ID)

blood_cells <- subset(Cell_Line_Info, lineage=="blood")
RNA_Protein_blood<-subset(RNA_Protein_CellLines, Cell_line %in% blood_cells$DepMap_ID)

skin_cells <- subset(Cell_Line_Info, lineage=="skin")
RNA_Protein_skin<-subset(RNA_Protein_CellLines, Cell_line %in% skin_cells$DepMap_ID)

lymphocyte_cells <- subset(Cell_Line_Info, lineage=="lymphocyte")
RNA_Protein_lymphocyte<-subset(RNA_Protein_CellLines, Cell_line %in% lymphocyte_cells$DepMap_ID)

CNS_cells <- subset(Cell_Line_Info, lineage=="central_nervous_system")
RNA_Protein_CNS<-subset(RNA_Protein_CellLines, Cell_line %in% CNS_cells$DepMap_ID)


