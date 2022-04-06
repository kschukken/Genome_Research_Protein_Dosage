##### Per specific cell line (ex MCF7): DNA RNA and Protein ratio relative to mean neutral ploidy ####
## 210215
## 210502-- update, only use genes with 10+ cells per category, in RNA and Protein
## Author: Klaske Schukken
## for Cell line X: Plot copy number per arm, average RNA & protein expression per arm
## Using only cells that have both RNA and Protein expression

## Protein expression data
##  https://www.cell.com/cell/fulltext/S0092-8674(19)31385-6#secsectitle0190
## Nusinow et al. Cell, 2020, Quantitative proteomics of the cancer cell line encyclopedia
## and chromosome arm data : Uri Ben-David paper 
## Data already normalized.

library('ggplot2')
library('tidyverse')
library(xlsx)
library(readxl)
library(reshape2)
library('BBmisc')

### Steps: 
# 1) Get data
# 2) get copy number per arm for specific cell line
# 3) get mean RNA expression of cell line, relative to diploid, per arm 
# 4) get mean Protein expression of cell line, relative to diploid, per arm 
# 5) Plot Copy number, RNA and Protein for 3-5 cells. seperate graphs. 

##### Step 1: Get data #####
### Step 1: Get data.  
## Get raw Protein and RNA Data. 

## Import Protein expression data and uri ben-david chromosome arm data
###
#setwd()#set working directory 
## Get depmap.org protein expression data and cell line information data. 
#              https://depmap.org/portal/download/
## Normalized Protein Expression data 
## Also get cell line info from depmap.org 
## https://www.cell.com/cell/fulltext/S0092-8674(19)31385-6#secsectitle0190
## Nusinow et al. Cell, 2020, Quantitative proteomics of the cancer cell line encyclopedia
## get expression data from depmap.org. 


Protein_Expression.1<-read_excel("mmc2.xlsx", 
                                 sheet= "Normalized Protein Expression")
Protein_Expression.1<-Protein_Expression.1[,c(1:426)] #delete empty collumns
#Name dataframe with protein_ID and corresponding Protein name, this way we have both gene ID and Protein name
Protein_ProID<-Protein_Expression.1[,c(1:2, 6)] 
Protein_ProID$Protein_Id<- str_replace_all(Protein_ProID$Protein_Id, "[|]", ".")
Protein_ProID$Protein_Id<- str_replace_all(Protein_ProID$Protein_Id, "[-]", ".")


#Protein Info. names, chrm location, etc. 
Protein_Info<-read_excel("HGNC_Protein_Info_edited.xlsx", sheet ="HGNC names")
Protein_Info2<-Protein_Info[,c(2,3,11,12,22)]

Protein_Info3<-Protein_Info2
Protein_Info3$arm <- gsub('[0-9]+', '', Protein_Info3$'Chromosome band') #Find test protein chromosome arm
Protein_Info3$arm <- gsub('[.]', '', Protein_Info3$arm) #also remove period, if needed
Protein_Info3$arm <- str_sub(Protein_Info3$arm, -1, -1) #start & end on last character. get only p or q
Protein_Info3$ChrmArm<- paste0(Protein_Info3$Chromosome, Protein_Info3$arm)
#filter cell data and get only aneuploidy scores for chrm arm of test protein location

write.csv(Protein_Info3, "Protein_location_info.csv")

###Get updated Protein info with only 371 cell lines (cell lines also have RNA Data)
# Get Protein.Expression.filtered and RNA.Expression.filtered from "Protein_RNA_filtered_CellLine.R" *
Protein.Expression.filtered<-read.delim2("Protein_Expression_filtered.csv", 
                                         dec=",", header = TRUE, sep=";")


####get cell line aneuploidy data:
#setwd()#set working directory
aneuploid<- read.csv("arm_calls_data_bendavid.csv", header=TRUE)


####
# Import RNA expression data (no aneuploidy data). CCLE depmap data.
# 1,303 cell lines
# 19,145 RNAs

#RNA Info. names, chrm location, etc. 
RNA_Info<-read_excel("HGNC_Protein_Info_edited.xlsx", sheet ="HGNC names")
RNA_Info2<-RNA_Info[,c(2,3,11,12,22,14)]

###Get updated Protein info with only 371 cell lines (cell lines also have RNA Data)
# Get Protein.Expression.filtered and RNA.Expression.filtered from "Protein_RNA_filtered_CellLine.R" *
RNA.Expression.filtered<-read.delim2("RNA_Expression_filtered.csv", 
                                     dec=",", header = TRUE, sep=";")




##### Filter Protein/RNA expression for genes that have a minimum of 10 cells/category
### Now filter Protein/RNA expression for genes that have a minimum of 10 cells per 
## category (gain, neutral, loss; in both RNA and protein): 
## 9414 genes in filtered data set 
## filter (10+ cells/category dataset from: Protein_RNA_Corr_min10.Filtered.R)

## This file is retrieved from / generated in Protein_RNA_expression.PerCell_v2.R
## can also be retrieved from supplementary data files
CN.Diff.xRNA.yProt.ThreeGroups<- read.csv("RNA.Protein_Mono.Di.Tri_Difference_Pvalue_min10points_3category.csv")

# Proteins: get collumns that have genes with 10+ cells/condition. 
#collumn names = "sp.Q9NQ94.A1CF_HUMAN" Protein ID format. same as CN.Diff.xRNA.yProt$Protein_ID
# select all collumns that have same ID as in filtered list. 
Protein.Expression.filtered_min10Cells<- Protein.Expression.filtered %>% select(one_of(CN.Diff.xRNA.yProt.ThreeGroups$Protein_ID))
# Now add Cell_lines back in. 
Protein.Expression.filtered_min10Cells$Cell_line<- Protein.Expression.filtered$Cell_line

# RNA: get collumns that have genes with 10+ cells/condition. 
# collumn names = "TSPAN6..7105." RNA ID format. same as CN.Diff.xRNA.yProt$RNA_ID
# select all collumns that have same ID as in filtered list. 
# length(unique(CN.Diff.xRNA.yProt$RNA_ID)), there are 9094 unique RNAs, 
# (9413 unique proteins, can have multiple protien isoforms per RNA)

RNA.Expression.filtered_min10Cells<- RNA.Expression.filtered %>% select(one_of(CN.Diff.xRNA.yProt.ThreeGroups$RNA_ID))
# Now add Cell_lines back in. 
RNA.Expression.filtered_min10Cells$Cell_line<- RNA.Expression.filtered$Cell_line




##### Step 2: get chrm arm copy number for specific cell line ####
## Step 2: get copy number per arm for specific cell line (~3 cell lines?)

## DLD1 = ACH-001061 --> do not have RNA & Protein data
# Hela = ACH-001086 --> no Protein
# MCF10A = ACH-001357 --> no Protein
## A2780 = ACH-000657 
## HCT116 = ACH-000971 
# MCF7 = ACH-000019
# LS 180 = ACH-000957
# SW48 = ACH-000958
# SW620= ACH-000651
# SKMEL2 = ACH-001190
# ACH-001307 = 8505C
# ACH-000755 = HCC2218
# ACH-000589 = NCIH1437
# NCIH838 =  ACH-000416 

# Make list of cell lines that 
List_Cell_ID<-c( "ACH-000019", "ACH-000416", "ACH-000657", "ACH-000957","ACH-000958", "ACH-000971", "ACH-001190")

#Check to see if cell lines are in the filtered cell line dataset (rna and prot info)
List_Cell_ID %in% Protein.Expression.filtered_min10Cells$Cell_line

An_Cell_data<-aneuploid[1,]# Make data frame with aneuploidy data
An_Cell_data[-1,]

#Now get all aneuploidy data for specific cells 
An_Cell_data<- aneuploid %>% filter(DepMap_ID %in% as.factor(List_Cell_ID))

#Now make a Chrm_Arm collumn
An_Cell_data$Chrm_Arm<- paste0(An_Cell_data$chrom, An_Cell_data$arm)

write.csv(An_Cell_data, 
          file= "Aneuploidy_Data_per_testcell.csv")
#An_Cell_data<- read_delim(delim=",",
#  file= "Aneuploidy_Data_per_testcell.csv")


An_Cell_data$DNA_ratio<- (An_Cell_data$wmed_CN/as.integer(An_Cell_data$ploidy))-1

#Sort by chrm arm level, not alphabetically
An_Cell_data$Chrm_Arm <- factor(An_Cell_data$Chrm_Arm, 
                                                 levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                            "4p", "4q", "5p", "5q", "6p", "6q",
                                                            "7p", "7q", "8p", "8q", "9p", "9q",
                                                            "10p", "10q", "11p", "11q", "12p", "12q",
                                                            "13p", "13q", "14p", "14q", "15p", "15q",
                                                            "16p", "16q", "17p", "17q", "18p", "18q",
                                                            "19p", "19q", "20p", "20q", "21p", "21q",
                                                            "22p", "22q", "Xp", "Xq"))
An_Cell_data$DepMap_ID <- factor(An_Cell_data$DepMap_ID, levels= unique(An_Cell_data$DepMap_ID))

ggplot(An_Cell_data, aes(x=DepMap_ID, y=Chrm_Arm))+ 
  geom_raster(aes(fill = DNA_ratio), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  xlab("Cell line")+
  ylab("Chromosome arms")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="DNA ratio") + #Chrm CN/ploidy -1
  coord_flip()+ #flip x and Y axis to arm on x.
  ggtitle("DNA arm call per cell line")

##### Step 3: get mean RNA expression of cell line, relative to diploid, per arm #####
## Step 3: get mean RNA expression of cell line, relative to diploid, per arm 
## Substep: 1) find cells with no gain/loss of chromosome X, arm p/q, and test cells
## substep: 2) Find difference in RNA expression per RNA
## substep: 3) find location of each RNA, chromosome and arm location
## substep: 4) find mean of difference in RNA exp per chromosome arm--> into datadrame
## substep: 5) repeat for all chromosome arms.
## substep: 6) plot mean change in RNA per chromosome, per chromosome gain. 

#Set up data
AllChrmArms<- data.frame(Chrm=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,14,15,16,16,17,17,18,18,19,19,20,20,21,22), 
                         Arm=c("p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q",
                               "q","q","q","p","q","p","q","p","q","p","q","p","q","q","q"))

#Before running, make RNA data for all RNA Expression data 
#with chrom arm info and RNA_ID equal to colnames on RNA expression data. 
#Make RNA expression names into data frame with ID and name
RNA.colnames<-data.frame(RNA_ID=colnames(RNA.Expression.filtered_min10Cells))
RNA.colnames$RNA_Name<-sub("[..].*", "", as.character(RNA.colnames$RNA_ID))
RNA.colnames$RNA_num<-sub(".*[.]{2}", "", as.character(RNA.colnames$RNA_ID))
RNA.colnames$RNA_num<-sub("[.]", "", as.character(RNA.colnames$RNA_num))

#add more RNA info per RNA. Add Uniprot_Acc and Gene_Symbol
RNA.onChrm1<-merge(x=RNA.colnames, y= RNA_Info2, 
                   by.x="RNA_Name", by.y="Protein_Name", 
                   sort = TRUE)# 8 collumns, 9090 genes



# make Chrm number & arm categories, and group by chrm num/arm
RNA.onChrm1$arm <- gsub('[0-9]+', '', RNA.onChrm1$'Chromosome band')
RNA.onChrm1$arm <- gsub('[.]', '', RNA.onChrm1$arm)
RNA.onChrm1$arm <- str_sub(RNA.onChrm1$arm, -1, -1) #start & end on character 1. get only p or q

RNA.onChrm1$Chrm.Num.Arm<- paste0(RNA.onChrm1$Chromosome, RNA.onChrm1$arm)
RNA.onChrm1$Chrm.Num.Arm<-as.factor(RNA.onChrm1$Chrm.Num.Arm) #as factor. 

#RNA.onChrm1 is used below in loop to subsection RNA_ID by chrm arm

#Set up empty datat frame to fill in with loop:
RNA.Diff.byChromosome<-data.frame(Chrm.Num.Arm=as.character(), 
                                  Difference=as.numeric(), 
                                  #PercDiff=as.numeric(), 
                                  Cell_ID=as.character())


###Run data below
##Repeat task for all chromosome arm gains:
for (i in 1:length(AllChrmArms$Chrm)) {
  TestChrm<-AllChrmArms$Chrm[i]
  TestArm<-AllChrmArms$Arm[i]
  TestChrmArm<-paste0(TestChrm, TestArm)
  
  #subset 1: Find cells with gain or no gain of chromosome X
  testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
  testchrm.percell <- filter(testchrm.percell, arm==TestArm)
  
  Cells.testCell.forchrm<- filter(testchrm.percell, (DepMap_ID %in% List_Cell_ID)) #all cells triploid for testRNA chrm
  Cells.Di.forchrm<- filter(testchrm.percell, arm_call==0) #all cells diploid for testRNA chrm
  
  #Substep 2: Get list of RNA on test chromosome
  ## take RNA names from RNA.Expression.filtered_min10Cells, and get only RNA on textChrmArm-->RNA.Expression.subset
  
  # per arm, only need RNA/Protein on that arm. 
  RNA.onChrm4<-RNA.onChrm1 #this is RNA info only, not RNA expression
  RNA.onChrm4<-subset(RNA.onChrm4, Chrm.Num.Arm == TestChrmArm) #get only RNA on chrm arm
  
  #Now subset RNA expression data with only chromosome arm of interest. 
  RNA.Expression.subset <- RNA.Expression.filtered_min10Cells %>% select(RNA.onChrm4$RNA_ID)
  RNA.Expression.subset$Cell_line<- RNA.Expression.filtered_min10Cells$Cell_line #add back cell info
  
  #Substep 3: get difference in RNA expression between diploid & test cells
  #Now get list of RNA expression in cells trisomic and disomic for each chrm arm
  RNA.testCell.Cells<-merge(y= RNA.Expression.subset, x= Cells.testCell.forchrm, 
                            by.y="Cell_line", by.x="DepMap_ID", 
                            sort = TRUE)
  RNA.Di.Cells<-merge(y= RNA.Expression.subset, x= Cells.Di.forchrm, 
                      by.y="Cell_line", by.x="DepMap_ID", 
                      sort = TRUE)
  
  # for each RNA, find mean Tri, mean Di, difference
  Diff.PerRNA<- data.frame(RNA_ID=as.character(), #Set up data.frame with 3 collumn
                           Difference=numeric(), 
                           #PercDiff=numeric(),
                           Cell_ID=as.character(),
                           stringsAsFactors = TRUE) #this is needed to prevent errors
  
  for (w in 1:length(List_Cell_ID)){ #for all cell lines in list:
    #Get difference in exp between diploid & test cells, per RNA
    for (t in 8:length(RNA.Di.Cells)) { #for all RNA in list
      Diff.PerRNA<- rbind(Diff.PerRNA, 
                          data.frame(RNA_ID=colnames(RNA.testCell.Cells[t]),
                                     Difference=RNA.testCell.Cells[w,t] - mean(RNA.Di.Cells[,t], na.rm=TRUE), 
                                     #PercDiff= (RNA.testCell.Cells[w,t] - mean(RNA.Di.Cells[,t], na.rm=TRUE))/mean(RNA.Di.Cells[,t], na.rm=TRUE),
                                     Cell_ID=List_Cell_ID[w]))
    }
  }
  ##
  # Get mean RNA expression difference by chrm arm, and per cell line
  Mean.Diff<-aggregate( Difference ~ Cell_ID, Diff.PerRNA, mean) 
  Mean.Diff$Chrm.Num.Arm<-as.factor(TestChrmArm)
  RNA.Diff.byChromosome<-rbind(RNA.Diff.byChromosome, Mean.Diff)
  
  print(paste0("Finished with chrm arm ", TestChrmArm))
}


### Now plot the heatmap for chromosome RNA chrm Gain changes. 
length(RNA.Diff.byChromosome$Chrm.Num.Arm) #273
write.csv(RNA.Diff.byChromosome, 
           file= "RNA_diff_PerCell_filtered.csv")
#RNA.Diff.byChromosome<- read_csv(
#          file= "/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/DNA_RNA_Protein_percell/RNA_diff_PerCell_filtered_210502.csv")


RNA.Diff.byChromosome$Chrm.Num.Arm <- factor(RNA.Diff.byChromosome$Chrm.Num.Arm, 
                                             levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                        "4p", "4q", "5p", "5q", "6p", "6q",
                                                        "7p", "7q", "8p", "8q", "9p", "9q",
                                                        "10p", "10q", "11p", "11q", "12p", "12q",
                                                        "13p", "13q", "14p", "14q", "15p", "15q",
                                                        "16p", "16q", "17p", "17q", "18p", "18q",
                                                        "19p", "19q", "20p", "20q", "21p", "21q",
                                                        "22p", "22q", "Xp", "Xq"))

RNA.Diff.byChromosome$Cell_ID <- factor(RNA.Diff.byChromosome$Cell_ID, levels= unique(RNA.Diff.byChromosome$Cell_ID))


ggplot(RNA.Diff.byChromosome, aes(x=Cell_ID, y=Chrm.Num.Arm))+ #3x10
  geom_raster(aes(fill = Difference), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  xlab("Cell lines")+
  ylab("RNA expression on chromosome arm")+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", 
                       midpoint=0, 
                       name="Difference") +
  coord_flip()+ #flip x and Y axis to arm on x
  ggtitle("Difference in RNA expresssion per chromosome arm between test cell \n and cells with no no gain/loss for chrm arm")



##### Step 4: get mean Protein expression of cell line, relative to diploid, per arm ####
### Step 4: get mean Protein expression of cell line, relative to diploid, per arm 
## List of cells: List_Cell_ID

### Using only 371 cell lines that have both RNA and Protein data
## Substep: 1) find cells with no gain/loss of chromosome X, arm p/q, and cell of interest
## substep: 2) find location of each protein, chromosome and arm location
## substep: 3) find difference in protein expression per protein
## substep: 4) find mean of difference in Prot exp per chromosome arm--> into dataframe
## substep: 5) repeat for all chromosome arms. then repeat per cell line
## substep: 6) plot mean change in prot per chromosome, per chromosome gain. 

#Set up data
AllChrmArms<- data.frame(Chrm=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,14,15,16,16,17,17,18,18,19,19,20,20,21,22), 
                         Arm=c("p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q",
                               "q","q","q","p","q","p","q","p","q","p","q","p","q","q","q"))

## make a list of protein info per protein_ID, so I can sort by prot location later.
# only need to do this once, then use data in for loop
Protein.colnames<-data.frame(Protein_ID=colnames(Protein.Expression.filtered_min10Cells)) #1 col, 12757 genes
# add more protein info per protein. Add Uniprot_Acc and Gene_Symbol
Protein.onChrm1<-merge(x= Protein.colnames, y= Protein_ProID, 
                       by.x="Protein_ID", by.y="Protein_Id", 
                       sort = TRUE)# 7 collumns, 12755 genes
#Find genes with same gene_Symbol
Protein.onChrm2<-merge(x= Protein.onChrm1, y= Protein_Info2, 
                       by.x="Gene_Symbol", by.y="Protein_Name", 
                       sort = TRUE)# 10 collumns, 12467 genes
#Of genes without matching gene symbol...
No_Gene_Symbol <- anti_join(Protein.onChrm1, Protein_Info2, #finding genes_Symbol without match
                            by = c("Gene_Symbol" = "Protein_Name"))
#...find the genes with matching uniprot_ID
Protein.onChrm3<-merge(x= No_Gene_Symbol, y= Protein_Info2, 
                       by.x="Uniprot_Acc", by.y="UniProt ID(supplied by UniProt)", 
                       sort = TRUE)# 10 collumns, 253 genes, only no gene-symbol genes with uniprot_ID
# Merge genes with gene symbol and genes with only uniprot ID. 
Protein.onChrm4<-merge(x= Protein.onChrm2, y= Protein.onChrm3, 
                       all=TRUE)# 12 collumns, 12720 genes


# make Chrm number & arm categories, and group by chrm num/arm
Protein.onChrm4$arm <- gsub('[0-9]+', '', Protein.onChrm4$'Chromosome band')
Protein.onChrm4$arm <- gsub('[.]', '', Protein.onChrm4$arm)
Protein.onChrm4$arm <- str_sub(Protein.onChrm4$arm, -1, -1) #start & end on character 1. get only p or q

Protein.onChrm4$Chrm.Num.Arm<- paste0(Protein.onChrm4$Chromosome, Protein.onChrm4$arm)
# 10 collumns, 12720 genes


# Set up empty datat frame to fill in with loop:
protein.Diff.byChromosome<-data.frame(Chrm.Num.Arm=as.character(), 
                                      Difference=as.numeric(), 
                                      #PercDiff=as.numeric(),
                                      Cell_ID=as.character())

###Run data below
##Repeat task for all chromosome arm gains:
for (i in 1:length(AllChrmArms$Chrm)) {
  TestChrm<-AllChrmArms$Chrm[i]
  TestArm<-AllChrmArms$Arm[i]
  TestChrmArm<-paste0(TestChrm, TestArm)
  
  #subset 1: Find cells with no gain/loss of ChrmArm, and test cells
  testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
  testchrm.percell <- filter(testchrm.percell, arm==TestArm)
  
  Cells.testcell.forchrm<- filter(testchrm.percell, (DepMap_ID %in% List_Cell_ID)) #all cells triploid for testProt chrm
  Cells.Di.forchrm<- filter(testchrm.percell, arm_call==0) #all cells diploid for testProt chrm
  
  # SubStep 2: get proteins located on TestChrmArm
  # First get protein info, loci, then subset protein expression data for testarm loci only
  # then link Protein values with location 
  #
  Protein.onChrm5<-Protein.onChrm4 #this is RNA info only, not RNA expression
  Protein.onChrm5<-subset(Protein.onChrm5, Chrm.Num.Arm == TestChrmArm) #get only RNA on chrm arm
  
  #Now subset RNA expression data with only chromosome arm of interest. 
  Protein.Expression.subset <- Protein.Expression.filtered_min10Cells %>% select(Protein.onChrm5$Protein_ID)
  Protein.Expression.subset$Cell_line<- Protein.Expression.filtered_min10Cells$Cell_line #add back cell info
  
  
  #Substep 3: get difference in protein expression between diploid & test cells
  #Now get list of protein expression in cells disomic & test cells for test ChrmArm
  Protein.TestCell.Cells<-merge(y= Protein.Expression.subset, x= Cells.testcell.forchrm, 
                                by.y="Cell_line", by.x="DepMap_ID", 
                                sort = TRUE)# 368 cells, 12757 proteins 
  Protein.Di.Cells<-merge(y= Protein.Expression.subset, x= Cells.Di.forchrm, 
                          by.y="Cell_line", by.x="DepMap_ID", 
                          sort = TRUE)# 368 cells, 12757 proteins 
  
  
  # for each protein, find mean Tri, mean Di, difference
  Diff.PerProtein<- data.frame(Protein_ID=as.character(), #Set up data.frame with 3 collumn
                               Difference=numeric(), 
                               #PercDiff=numeric(),
                               Cell_ID=as.character(),
                               stringsAsFactors = TRUE) #this is needed to prevent errors
  
  for (w in 1:length(List_Cell_ID)){ #for all cell lines in list:
    for (t in 8:length(Protein.Di.Cells)) { #Get difference in exp between cell line x and mean no gain/loss cells cells, per protein
      Diff.PerProtein<- rbind(Diff.PerProtein, 
                              data.frame(Protein_ID=colnames(Protein.TestCell.Cells[t]),
                                         Difference=Protein.TestCell.Cells[w,t] - mean(Protein.Di.Cells[,t], na.rm=TRUE), 
                                         #PercDiff= (Protein.TestCell.Cells[w,t] - mean(Protein.Di.Cells[,t], na.rm=TRUE))/(mean(Protein.Di.Cells[,t], na.rm=TRUE)),
                                         Cell_ID=List_Cell_ID[w]))
    }
  }
  ##
  # Get mean protein expression difference by chrm num/arm category
  Mean.Diff<-aggregate( Difference ~ Cell_ID, Diff.PerProtein, mean ) 
  Mean.Diff$Chrm.Num.Arm<-TestChrmArm
  
  protein.Diff.byChromosome<-rbind(protein.Diff.byChromosome, Mean.Diff)
  print(paste0("finished with chrm arm ", TestChrmArm))
}



### Now plot the heatmap for chromosome Protein chrm Gain changes. 
length(protein.Diff.byChromosome$Chrm.Num.Arm)
#setwd()
write.csv(protein.Diff.byChromosome, 
          file= "Protein_diff_PerCell_filtered.csv")
#protein.Diff.byChromosome<-read.csv(file="Protein_diff_PerCell_filtered_210502.csv", 
#                                     dec=".", header = TRUE, sep=",")


protein.Diff.byChromosome$Chrm.Num.Arm <- factor(protein.Diff.byChromosome$Chrm.Num.Arm, 
                                                 levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                            "4p", "4q", "5p", "5q", "6p", "6q",
                                                            "7p", "7q", "8p", "8q", "9p", "9q",
                                                            "10p", "10q", "11p", "11q", "12p", "12q",
                                                            "13p", "13q", "14p", "14q", "15p", "15q",
                                                            "16p", "16q", "17p", "17q", "18p", "18q",
                                                            "19p", "19q", "20p", "20q", "21p", "21q",
                                                            "22p", "22q", "Xp", "Xq"))

protein.Diff.byChromosome$Cell_ID <- factor(protein.Diff.byChromosome$Cell_ID, levels= unique(protein.Diff.byChromosome$Cell_ID))

ggplot(protein.Diff.byChromosome, aes(x=Cell_ID, y=Chrm.Num.Arm))+ #3x10 figure
  geom_raster(aes(fill = Difference), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  xlab("Cell line")+
  ylab("Protein expression per chrm arm")+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference") +
  coord_flip() + #flip x and Y axis to arm on x.
  ggtitle("Difference in protein expresssion per chromosome arm\n between test cell line and cells with no gain")

##### Step 5: Combining DNA, RNA and Protein Data, plotting per cell  #####
### Combining DNA, RNA and Protein Data, plotting per cell 
## 1) Reanme difference with "DNA", "RNA" or "Protein" to data frames, add cell&arm ID
## 2) Combine data frames
## 3) Make data frame per cell.
## 4) Plot data per cell: DNA, RNA and Protein

## MCF7
#1)
#An_Cell_data
RNA.Diff.byChromosome$RNADifference<- RNA.Diff.byChromosome$Difference
RNA.Diff.byChromosome$Cell.arm<- paste0(RNA.Diff.byChromosome$Cell_ID,RNA.Diff.byChromosome$Chrm.Num.Arm)
protein.Diff.byChromosome$ProteinDifference<- protein.Diff.byChromosome$Difference
protein.Diff.byChromosome$Cell.arm<- paste0(protein.Diff.byChromosome$Cell_ID,protein.Diff.byChromosome$Chrm.Num.Arm)


#2) DNA.ACH000019 -- MCF7
DNA.ACH000019<-subset(An_Cell_data, DepMap_ID == "ACH-000019")
RNA.ACH000019<-subset(RNA.Diff.byChromosome, Cell_ID == "ACH-000019")
protein.ACH000019<-subset(protein.Diff.byChromosome, Cell_ID == "ACH-000019")

Cell.DNA.RNA<-merge(x=DNA.ACH000019, y= RNA.ACH000019, 
                    by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                    sort = TRUE)# 
Cell.DNA.RNA.Protein<-merge(x=Cell.DNA.RNA, y= protein.ACH000019, 
                            by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                            sort = TRUE)# 

Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein %>% 
  rename( #new = old
    RNADifference = Difference.x,
    ProteinDifference = Difference.y
  )

Cell.DNA.RNA.Protein$arm_call<-Cell.DNA.RNA.Protein$arm_call/2
Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein[ , which(names(Cell.DNA.RNA.Protein) 
                                   %in% c("Chrm_Arm", "RNADifference", "ProteinDifference", "DNA_ratio"))]
Cell.DNA.RNA.Protein.m <- melt(Cell.DNA.RNA.Protein)

Cell.DNA.RNA.Protein.m$variable <- factor(Cell.DNA.RNA.Protein.m$variable, 
                                levels = c("ProteinDifference", "RNADifference", "DNA_ratio"))

cor.test(Cell.DNA.RNA.Protein$DNA_ratio, Cell.DNA.RNA.Protein$RNADifference, method="pearson")
# p= 2E-06, r=0.68
cor.test(Cell.DNA.RNA.Protein$DNA_ratio, Cell.DNA.RNA.Protein$ProteinDifference, method="pearson")
# p= 1E-05, r=0.63

#Heatmap DNA, RNA, Protein
cell_name= "MCF7"
pdf(file = paste0("plot.Heatmap.", cell_name, ".DNA.RNA.Protein_min10cell.pdf"),
    width = 8, 
    height = 3)
ggplot(data= Cell.DNA.RNA.Protein.m, 
       aes(x=variable, y=Chrm_Arm))+ #3x10
  geom_raster(aes(fill = value), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.7, 0.7)) + 
  xlab("MCF7")+
  ylab("Chromosome arm")+
  coord_flip()+ #flip x and Y axis to arm on x
  ggtitle("DNA ratio, difference in RNA and Protein expression \n per chromosome arm in MCF7, ACH-000019")
# plot.Heatmap.MCF7.DNA.RNA.Protein_min10cell
# 8x3
dev.off()

pdf(file = paste0("plot.Scatter.", cell_name, ".DNA.RNA.pdf"),
    width = 4, 
    height = 4)
# Scatterplot: DNA RNA
ggplot(data= Cell.DNA.RNA.Protein, 
       aes(x=DNA_ratio, y=RNADifference))+ #4x4
  geom_point(size=3)+
  coord_cartesian(xlim=c(-.7,.7), ylim=c(-.7,.7))+
  theme_classic()+
  geom_smooth(color="Red", method="lm") +
  xlab("DNA ratio")+
  ylab("Difference in RNA expression")
#ggtitle("DNA ratio vs difference in Protein expression \n per chromosome arm in SKMEL2, ACH-001190")
#4x4
# plot.Scatter.MCF7.DNA.RNA
dev.off()

pdf(file = paste0("plot.Scatter.", cell_name, ".DNA.Protein.pdf"),
    width = 4, 
    height = 4)
# Scatterplot: DNA Protein
ggplot(data= Cell.DNA.RNA.Protein, 
       aes(x=DNA_ratio, y=ProteinDifference))+ #4x4
  geom_point(size=3)+
  coord_cartesian(xlim=c(-.7,.7), ylim=c(-.7,.7))+
  theme_classic()+
  geom_smooth(color="Red", method="lm") +
  xlab("DNA ratio")+
  ylab("Difference in Protein expression")
#4x4
# plot.Scatter.MCF7.DNA.Protein
dev.off()


###

#2) HCT116 = ACH-000971
DNA.ACH000019<-subset(An_Cell_data, DepMap_ID == "ACH-000971")
RNA.ACH000019<-subset(RNA.Diff.byChromosome, Cell_ID == "ACH-000971")
protein.ACH000019<-subset(protein.Diff.byChromosome, Cell_ID == "ACH-000971")

Cell.DNA.RNA<-merge(x=DNA.ACH000019, y= RNA.ACH000019, 
                    by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                    sort = TRUE)# 
Cell.DNA.RNA.Protein<-merge(x=Cell.DNA.RNA, y= protein.ACH000019, 
                            by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                            sort = TRUE)# 

Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein %>% 
  rename( #new = old
    RNADifference = Difference.x,
    ProteinDifference = Difference.y
  )

Cell.DNA.RNA.Protein$arm_call<-Cell.DNA.RNA.Protein$arm_call/2
Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein[ , which(names(Cell.DNA.RNA.Protein) 
                                                     %in% c("Chrm_Arm", "RNADifference", "ProteinDifference", "DNA_ratio"))]
Cell.DNA.RNA.Protein.m <- melt(Cell.DNA.RNA.Protein)

Cell.DNA.RNA.Protein.m$variable <- factor(Cell.DNA.RNA.Protein.m$variable, 
                                          levels = c("ProteinDifference", "RNADifference", "DNA_ratio"))

cor.test(Cell.DNA.RNA.Protein$DNA_ratio, Cell.DNA.RNA.Protein$RNADifference, method="pearson")
# p= NS, r=0.25
cor.test(Cell.DNA.RNA.Protein$DNA_ratio, Cell.DNA.RNA.Protein$ProteinDifference, method="pearson")
# p= NS, r=0.38

cell_name= "HCT116"
pdf(file = paste0("plot.Heatmap.", cell_name, ".DNA.RNA.Protein_min10cell.pdf"),
    width = 8, 
    height = 3)
ggplot(data= Cell.DNA.RNA.Protein.m, 
       aes(x=variable, y=Chrm_Arm))+ #3x10
  geom_raster(aes(fill = value), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.7, 0.7)) + 
  xlab("HCT116")+
  ylab("Chromosome arm")+
  coord_flip()+ #flip x and Y axis to arm on x
  ggtitle("DNA ratio, difference in RNA and Protein expression \n per chromosome arm in HCT116, ACH-000971")
# 8x3 
# plot.Heatmap.HCT116.DNA.RNA.Protein_min10cell
dev.off()

pdf(file = paste0("plot.Scatter.", cell_name, ".DNA.RNA.pdf"),
    width = 4, 
    height = 4)
# Scatterplot: DNA RNA
ggplot(data= Cell.DNA.RNA.Protein, 
       aes(x=DNA_ratio, y=RNADifference))+ #4x4
  geom_point(size=3)+
  coord_cartesian(xlim=c(-.7,.7), ylim=c(-.7,.7))+
  theme_classic()+
  geom_smooth(color="Red", method="lm") +
  xlab("DNA ratio")+
  ylab("Difference in RNA expression")
#ggtitle("DNA ratio vs difference in Protein expression \n per chromosome arm in SKMEL2, ACH-001190")
#4x4
# plot.Scatter.HCT116.DNA.RNA
dev.off()

pdf(file = paste0("plot.Scatter.", cell_name, ".DNA.Protein.pdf"),
    width = 4, 
    height = 4)
# Scatterplot: DNA Protein
ggplot(data= Cell.DNA.RNA.Protein, 
       aes(x=DNA_ratio, y=ProteinDifference))+ #4x4
  geom_point(size=3)+
  coord_cartesian(xlim=c(-.7,.7), ylim=c(-.7,.7))+
  theme_classic()+
  geom_smooth(color="Red", method="lm") +
  xlab("DNA ratio")+
  ylab("Difference in Protein expression")
#4x4
# plot.Scatter.HCT116.DNA.Protein
dev.off()


###
#2) ACH000657- A2780
DNA.ACH000019<-subset(An_Cell_data, DepMap_ID == "ACH-000657")
RNA.ACH000019<-subset(RNA.Diff.byChromosome, Cell_ID == "ACH-000657")
protein.ACH000019<-subset(protein.Diff.byChromosome, Cell_ID == "ACH-000657")

Cell.DNA.RNA<-merge(x=DNA.ACH000019, y= RNA.ACH000019, 
                    by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                    sort = TRUE)# 
Cell.DNA.RNA.Protein<-merge(x=Cell.DNA.RNA, y= protein.ACH000019, 
                            by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                            sort = TRUE)# 
Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein %>% 
  rename( #new = old
    RNADifference = Difference.x,
    ProteinDifference = Difference.y
  )

Cell.DNA.RNA.Protein$arm_call<-Cell.DNA.RNA.Protein$arm_call/2
Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein[ , which(names(Cell.DNA.RNA.Protein) 
                                                     %in% c("Chrm_Arm", "RNADifference", "ProteinDifference", "DNA_ratio"))]
Cell.DNA.RNA.Protein.m <- melt(Cell.DNA.RNA.Protein)

Cell.DNA.RNA.Protein.m$variable <- factor(Cell.DNA.RNA.Protein.m$variable, 
                                          levels = c("ProteinDifference", "RNADifference", "DNA_ratio"))

cor.test(Cell.DNA.RNA.Protein$DNA_ratio, Cell.DNA.RNA.Protein$RNADifference, method="pearson")
# p= 7E-3, r=0.42
cor.test(Cell.DNA.RNA.Protein$DNA_ratio, Cell.DNA.RNA.Protein$ProteinDifference, method="pearson")
# p= 0.01, r=0.40

cell_name= "A2780"
pdf(file = paste0("plot.Heatmap.", cell_name, ".DNA.RNA.Protein_min10cell.pdf"),
    width = 8, 
    height = 3)
ggplot(data= Cell.DNA.RNA.Protein.m, 
       aes(x=variable, y=Chrm_Arm))+ #3x10
  geom_raster(aes(fill = value), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.7, 0.7)) + 
  xlab("A2780")+
  ylab("Chromosome arm")+
  coord_flip()+ #flip x and Y axis to arm on x
  ggtitle("Arm call, difference in RNA and Protein expression \n per chromosome arm in A2780 ACH-000657")
# 8x3
# plot.Heatmap.A2780.DNA.RNA.Protein_min10cell
dev.off()

pdf(file = paste0("plot.Scatter.", cell_name, ".DNA.RNA.pdf"),
    width = 4, 
    height = 4)
# Scatterplot: DNA RNA
ggplot(data= Cell.DNA.RNA.Protein, 
       aes(x=DNA_ratio, y=RNADifference))+ #4x4
  geom_point(size=3)+
  coord_cartesian(xlim=c(-.7,.7), ylim=c(-.7,.7))+
  theme_classic()+
  geom_smooth(color="Red", method="lm") +
  xlab("DNA ratio")+
  ylab("Difference in RNA expression")
#ggtitle("DNA ratio vs difference in Protein expression \n per chromosome arm in SKMEL2, ACH-001190")
#4x4
# plot.Scatter.A2780.DNA.RNA
dev.off()

pdf(file = paste0("plot.Scatter.", cell_name, ".DNA.Protein.pdf"),
    width = 4, 
    height = 4)
# Scatterplot: DNA Protein
ggplot(data= Cell.DNA.RNA.Protein, 
       aes(x=DNA_ratio, y=ProteinDifference))+ #4x4
  geom_point(size=3)+
  coord_cartesian(xlim=c(-.7,.7), ylim=c(-.7,.7))+
  theme_classic()+
  geom_smooth(color="Red", method="lm") +
  xlab("DNA ratio")+
  ylab("Difference in Protein expression")
#4x4
# plot.Scatter.A2780.DNA.Protein
dev.off()

#####          Plot cells ####

#2) LS 180 - ACH-000957
DNA.ACH000019<-subset(An_Cell_data, DepMap_ID == "ACH-000957")
RNA.ACH000019<-subset(RNA.Diff.byChromosome, Cell_ID == "ACH-000957")
protein.ACH000019<-subset(protein.Diff.byChromosome, Cell_ID == "ACH-000957")

Cell.DNA.RNA<-merge(x=DNA.ACH000019, y= RNA.ACH000019, 
                    by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                    sort = TRUE)# 
Cell.DNA.RNA.Protein<-merge(x=Cell.DNA.RNA, y= protein.ACH000019, 
                            by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                            sort = TRUE)# 
Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein %>% 
  rename( #new = old
    RNADifference = Difference.x,
    ProteinDifference = Difference.y
  )

Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein[ , which(names(Cell.DNA.RNA.Protein) 
                                                     %in% c("Chrm_Arm", "RNADifference", "ProteinDifference", "DNA_ratio"))]
Cell.DNA.RNA.Protein.m <- melt(Cell.DNA.RNA.Protein)

Cell.DNA.RNA.Protein.m$variable <- factor(Cell.DNA.RNA.Protein.m$variable, 
                                          levels = c("ProteinDifference", "RNADifference", "DNA_ratio"))

cor.test(Cell.DNA.RNA.Protein$DNA_ratio, Cell.DNA.RNA.Protein$RNADifference, method="pearson")
# p=NS , r=-0.09
cor.test(Cell.DNA.RNA.Protein$DNA_ratio, Cell.DNA.RNA.Protein$ProteinDifference, method="pearson")
# p=NS , r=-0.03

cell_name= "LS180"
pdf(file = paste0("plot.Heatmap.", cell_name, ".DNA.RNA.Protein_min10cell.pdf"),
    width = 8, 
    height = 3)
ggplot(data= Cell.DNA.RNA.Protein.m, 
       aes(x=variable, y=Chrm_Arm))+ #3x10
  geom_raster(aes(fill = value), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.7, 0.7)) + 
  xlab("LS 180")+
  ylab("Chromosome arm")+
  coord_flip()+ #flip x and Y axis to arm on x
  ggtitle("Arm call, difference in RNA and Protein expression \n per chromosome arm in LS 180 ACH-000957")
# 8x3
# plot.Heatmap.LS180.DNA.RNA.Protein_min10cell
dev.off()

pdf(file = paste0("plot.Scatter.", cell_name, ".DNA.RNA.pdf"),
    width = 4, 
    height = 4)
# Scatterplot: DNA RNA
ggplot(data= Cell.DNA.RNA.Protein, 
       aes(x=DNA_ratio, y=RNADifference))+ #4x4
  geom_point(size=3)+
  coord_cartesian(xlim=c(-.7,.7), ylim=c(-.7,.7))+
  theme_classic()+
  geom_smooth(color="Red", method="lm") +
  xlab("DNA ratio")+
  ylab("Difference in RNA expression")
#ggtitle("DNA ratio vs difference in Protein expression \n per chromosome arm in SKMEL2, ACH-001190")
#4x4
# plot.Scatter.LS180.DNA.RNA
dev.off()

pdf(file = paste0("plot.Scatter.", cell_name, ".DNA.Protein.pdf"),
    width = 4, 
    height = 4)
# Scatterplot: DNA Protein
ggplot(data= Cell.DNA.RNA.Protein, 
       aes(x=DNA_ratio, y=ProteinDifference))+ #4x4
  geom_point(size=3)+
  coord_cartesian(xlim=c(-.7,.7), ylim=c(-.7,.7))+
  theme_classic()+
  geom_smooth(color="Red", method="lm") +
  xlab("DNA ratio")+
  ylab("Difference in Protein expression")
#4x4
# plot.Scatter.LS180.DNA.Protein
dev.off()



###
#2) SKMEL2 = ACH-001190
DNA.ACH000019<-subset(An_Cell_data, DepMap_ID == "ACH-001190")
RNA.ACH000019<-subset(RNA.Diff.byChromosome, Cell_ID == "ACH-001190")
protein.ACH000019<-subset(protein.Diff.byChromosome, Cell_ID == "ACH-001190")

Cell.DNA.RNA<-merge(x=DNA.ACH000019, y= RNA.ACH000019, 
                    by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                    sort = TRUE)# 
Cell.DNA.RNA.Protein<-merge(x=Cell.DNA.RNA, y= protein.ACH000019, 
                            by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                            sort = TRUE)# 

Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein %>% 
  rename( #new = old
    RNADifference = Difference.x
    #ProteinDifference = Difference.y
  )

Cell.DNA.RNA.Protein$arm_call<-Cell.DNA.RNA.Protein$arm_call/2
Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein[ , which(names(Cell.DNA.RNA.Protein) 
                                                     %in% c("Chrm_Arm", "RNADifference", "ProteinDifference", "DNA_ratio"))]
Cell.DNA.RNA.Protein.m <- melt(Cell.DNA.RNA.Protein)

Cell.DNA.RNA.Protein.m$variable <- factor(Cell.DNA.RNA.Protein.m$variable, 
                                          levels = c("ProteinDifference", "RNADifference", "DNA_ratio"))

cor.test(Cell.DNA.RNA.Protein$DNA_ratio, Cell.DNA.RNA.Protein$RNADifference, method="pearson")
# p= 4E-08, r=0.75
cor.test(Cell.DNA.RNA.Protein$DNA_ratio, Cell.DNA.RNA.Protein$ProteinDifference, method="pearson")
# p= 3E-08, r=0.76

cell_name= "SKMEL2"
pdf(file = paste0("plot.Heatmap.", cell_name, ".DNA.RNA.Protein_min10cell.pdf"),
    width = 8, 
    height = 3)
#4) Plot data per cell: DNA, RNA and Protein
ggplot(data= Cell.DNA.RNA.Protein.m, 
       aes(x=variable, y=Chrm_Arm))+ #3x10
  geom_raster(aes(fill = value), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.71, 0.7)) + 
  xlab("SKMEL2")+
  ylab("Chromosome arm")+
  coord_flip()+ #flip x and Y axis to arm on x
  ggtitle("DNA ratio, difference in RNA and Protein expression \n per chromosome arm in SKMEL2, ACH-001190")
# 8x3
# plot.Heatmap.SKMEL2.DNA.RNA.Protein_min10cell
dev.off()

pdf(file = paste0("plot.Scatter.", cell_name, ".DNA.RNA.pdf"),
    width = 4, 
    height = 4)
# Scatterplot: DNA RNA
ggplot(data= Cell.DNA.RNA.Protein, 
       aes(x=DNA_ratio, y=RNADifference))+ #4x4
  geom_point(size=3)+
  coord_cartesian(xlim=c(-.7,.7), ylim=c(-.71,.7))+
  theme_classic()+
  geom_smooth(color="Red", method="lm") +
  xlab("DNA ratio")+
  ylab("Difference in RNA expression")
#ggtitle("DNA ratio vs difference in Protein expression \n per chromosome arm in SKMEL2, ACH-001190")
#4x4
# plot.Scatter.SKMEL2.DNA.RNA
dev.off()

pdf(file = paste0("plot.Scatter.", cell_name, ".DNA.Protein.pdf"),
    width = 4, 
    height = 4)
# Scatterplot: DNA Protein
ggplot(data= Cell.DNA.RNA.Protein, 
       aes(x=DNA_ratio, y=ProteinDifference))+ #4x4
  geom_point(size=3)+
  coord_cartesian(xlim=c(-.7,.7), ylim=c(-.7,.7))+
  theme_classic()+
  geom_smooth(color="Red", method="lm") +
  xlab("DNA ratio")+
  ylab("Difference in Protein expression")
#4x4
# plot.Scatter.SKMEL2.DNA.Protein
dev.off()

####
# ACH-000416 = NCIH838
DNA.ACH000019<-subset(An_Cell_data, DepMap_ID == "ACH-000416")
RNA.ACH000019<-subset(RNA.Diff.byChromosome, Cell_ID == "ACH-000416")
protein.ACH000019<-subset(protein.Diff.byChromosome, Cell_ID == "ACH-000416")

Cell.DNA.RNA<-merge(x=DNA.ACH000019, y= RNA.ACH000019, 
                    by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                    sort = TRUE)# 
Cell.DNA.RNA.Protein<-merge(x=Cell.DNA.RNA, y= protein.ACH000019, 
                            by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                            sort = TRUE)# 

Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein %>% 
  rename( #new = old
    RNADifference = Difference.x
    #ProteinDifference = Difference.y
  )

Cell.DNA.RNA.Protein$arm_call<-Cell.DNA.RNA.Protein$arm_call/2
Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein[ , which(names(Cell.DNA.RNA.Protein) 
                                                     %in% c("Chrm_Arm", "RNADifference", "ProteinDifference", "DNA_ratio"))]
Cell.DNA.RNA.Protein.m <- melt(Cell.DNA.RNA.Protein)

Cell.DNA.RNA.Protein.m$variable <- factor(Cell.DNA.RNA.Protein.m$variable, 
                                          levels = c("ProteinDifference", "RNADifference", "DNA_ratio"))

cor.test(Cell.DNA.RNA.Protein$DNA_ratio, Cell.DNA.RNA.Protein$RNADifference, method="pearson")
# p= 5E-07, r=0.70
cor.test(Cell.DNA.RNA.Protein$DNA_ratio, Cell.DNA.RNA.Protein$ProteinDifference, method="pearson")
# p= 4E-03, r=0.45


cell_name= "NCIH838"
pdf(file = paste0("plot.Heatmap.", cell_name, ".DNA.RNA.Protein_min10cell.pdf"),
    width = 8, 
    height = 3)
ggplot(data= Cell.DNA.RNA.Protein.m, 
       aes(x=variable, y=Chrm_Arm))+ #3x10
  geom_raster(aes(fill = value), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.7, 0.7)) + 
  xlab("NCIH838")+
  ylab("Chromosome arm")+
  coord_flip()+ #flip x and Y axis to arm on x
  ggtitle("DNA ratio, difference in RNA and Protein expression \n per chromosome arm in NCIH838, ACH-000416")
# 8x3
# plot.Heatmap.NCIH838.DNA.RNA.Protein_min10cell
dev.off()

pdf(file = paste0("plot.Scatter.", cell_name, ".DNA.RNA.pdf"),
    width = 4, 
    height = 4)
# Scatterplot: DNA RNA
ggplot(data= Cell.DNA.RNA.Protein, 
       aes(x=DNA_ratio, y=RNADifference))+ #4x4
  geom_point(size=3)+
  coord_cartesian(xlim=c(-.7,.7), ylim=c(-.7,.7))+
  theme_classic()+
  geom_smooth(color="Red", method="lm") +
  xlab("DNA ratio")+
  ylab("Difference in RNA expression")
#ggtitle("DNA ratio vs difference in Protein expression \n per chromosome arm in SKMEL2, ACH-001190")
#4x4
# plot.Scatter.NCIH838.DNA.RNA
dev.off()

pdf(file = paste0("plot.Scatter.", cell_name, ".DNA.Protein.pdf"),
    width = 4, 
    height = 4)
# Scatterplot: DNA Protein
ggplot(data= Cell.DNA.RNA.Protein, 
       aes(x=DNA_ratio, y=ProteinDifference))+ #4x4
  geom_point(size=3)+
  coord_cartesian(xlim=c(-.7,.7), ylim=c(-.7,.7))+
  theme_classic()+
  geom_smooth(color="Red", method="lm") +
  xlab("DNA ratio")+
  ylab("Difference in Protein expression")
#4x4
# plot.Scatter.NCIH838.DNA.Protein
dev.off()

#####          Optional: More cells #####
### Combining DNA, RNA and Protein Data, plotting per cell - continued
# ACH-000958 = SW48 
# ACH-000651 = SW620
# ACH-001307 = 8505C
# ACH-000755 = HCC2218
# ACH-000589 = NCIH1437


####
###
#2) SW48 = ACH-000958
DNA.ACH000019<-subset(An_Cell_data, DepMap_ID == "ACH-000958")
RNA.ACH000019<-subset(RNA.Diff.byChromosome, Cell_ID == "ACH-000651")
protein.ACH000019<-subset(protein.Diff.byChromosome, Cell_ID == "ACH-000651")

Cell.DNA.RNA<-merge(x=DNA.ACH000019, y= RNA.ACH000019, 
                    by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                    sort = TRUE)# 
Cell.DNA.RNA.Protein<-merge(x=Cell.DNA.RNA, y= protein.ACH000019, 
                            by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                            sort = TRUE)# 

Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein %>% 
  rename( #new = old
    RNADifference = Difference.x,
    ProteinDifference = Difference.y
  )

Cell.DNA.RNA.Protein$arm_call<-Cell.DNA.RNA.Protein$arm_call/2
Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein[ , which(names(Cell.DNA.RNA.Protein) 
                                                     %in% c("Chrm_Arm", "RNADifference", "ProteinDifference", "DNA_ratio"))]
Cell.DNA.RNA.Protein.m <- melt(Cell.DNA.RNA.Protein)

Cell.DNA.RNA.Protein.m$variable <- factor(Cell.DNA.RNA.Protein.m$variable, 
                                          levels = c("ProteinDifference", "RNADifference", "DNA_ratio"))

cell_name= "SW48"
pdf(file = paste0("plot.Heatmap.", cell_name, ".DNA.RNA.Protein_min10cell.pdf"),
    width = 8, 
    height = 3)
ggplot(data= Cell.DNA.RNA.Protein.m, 
       aes(x=variable, y=Chrm_Arm))+ #3x10
  geom_raster(aes(fill = value), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.7, 0.7)) + 
  xlab("SW48")+
  ylab("Chromosome arm")+
  coord_flip()+ #flip x and Y axis to arm on x
  ggtitle("DNA ratio, difference in RNA and Protein expression \n per chromosome arm in SW48, ACH-000958")
# 8x3
# plot.Heatmap.SW48.DNA.RNA.Protein_min10cell
dev.off()

pdf(file = paste0("plot.Scatter.", cell_name, ".DNA.RNA.pdf"),
    width = 4, 
    height = 4)
# Scatterplot: DNA RNA
ggplot(data= Cell.DNA.RNA.Protein, 
       aes(x=DNA_ratio, y=RNADifference))+ #4x4
  geom_point(size=3)+
  coord_cartesian(xlim=c(-.7,.7), ylim=c(-.7,.7))+
  theme_classic()+
  geom_smooth(color="Red", method="lm") +
  xlab("DNA ratio")+
  ylab("Difference in RNA expression")
#ggtitle("DNA ratio vs difference in Protein expression \n per chromosome arm in SKMEL2, ACH-001190")
#4x4
# plot.Scatter.SW48.DNA.RNA
dev.off()

pdf(file = paste0("plot.Scatter.", cell_name, ".DNA.Protein.pdf"),
    width = 4, 
    height = 4)
# Scatterplot: DNA Protein
ggplot(data= Cell.DNA.RNA.Protein, 
       aes(x=DNA_ratio, y=ProteinDifference))+ #4x4
  geom_point(size=3)+
  coord_cartesian(xlim=c(-.7,.7), ylim=c(-.7,.7))+
  theme_classic()+
  geom_smooth(color="Red", method="lm") +
  xlab("DNA ratio")+
  ylab("Difference in Protein expression")
#4x4
# plot.Scatter.SW48.DNA.Protein
dev.off()
#SW620= ACH-000651
DNA.ACH000019<-subset(An_Cell_data, DepMap_ID == "ACH-000651")
RNA.ACH000019<-subset(RNA.Diff.byChromosome, Cell_ID == "ACH-000651")
protein.ACH000019<-subset(protein.Diff.byChromosome, Cell_ID == "ACH-000651")

Cell.DNA.RNA<-merge(x=DNA.ACH000019, y= RNA.ACH000019, 
                    by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                    sort = TRUE)# 
Cell.DNA.RNA.Protein<-merge(x=Cell.DNA.RNA, y= protein.ACH000019, 
                            by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                            sort = TRUE)# 

Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein %>% 
  rename( #new = old
    RNADifference = Difference.x,
    ProteinDifference = Difference.y
  )

Cell.DNA.RNA.Protein$arm_call<-Cell.DNA.RNA.Protein$arm_call/2
Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein[ , which(names(Cell.DNA.RNA.Protein) 
                                                     %in% c("Chrm_Arm", "RNADifference", "ProteinDifference", "DNA_ratio"))]
Cell.DNA.RNA.Protein.m <- melt(Cell.DNA.RNA.Protein)

Cell.DNA.RNA.Protein.m$variable <- factor(Cell.DNA.RNA.Protein.m$variable, 
                                          levels = c("ProteinDifference", "RNADifference", "DNA_ratio"))


#4) Plot data per cell: DNA, RNA and Protein
ggplot(data= Cell.DNA.RNA.Protein.m, 
       aes(x=variable, y=Chrm_Arm))+ #3x10
  geom_raster(aes(fill = value), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.7, 0.7)) + 
  xlab("SW620")+
  ylab("Chromosome arm")+
  coord_flip()+ #flip x and Y axis to arm on x
  ggtitle("DNA ratio, difference in RNA and Protein expression \n per chromosome arm in SW620, ACH-000651")
# 8x3
# plot.Heatmap.SW620.DNA.RNA.Protein_min10cell

####
#ACH-001307 = 8505C
DNA.ACH000019<-subset(An_Cell_data, DepMap_ID == "ACH-001307")
RNA.ACH000019<-subset(RNA.Diff.byChromosome, Cell_ID == "ACH-001307")
protein.ACH000019<-subset(protein.Diff.byChromosome, Cell_ID == "ACH-001307")

Cell.DNA.RNA<-merge(x=DNA.ACH000019, y= RNA.ACH000019, 
                    by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                    sort = TRUE)# 
Cell.DNA.RNA.Protein<-merge(x=Cell.DNA.RNA, y= protein.ACH000019, 
                            by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                            sort = TRUE)# 

Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein %>% 
  rename( #new = old
    RNADifference = Difference.x,
    ProteinDifference = Difference.y
  )

Cell.DNA.RNA.Protein$arm_call<-Cell.DNA.RNA.Protein$arm_call/2
Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein[ , which(names(Cell.DNA.RNA.Protein) 
                                                     %in% c("Chrm_Arm", "RNADifference", "ProteinDifference", "DNA_ratio"))]
Cell.DNA.RNA.Protein.m <- melt(Cell.DNA.RNA.Protein)

Cell.DNA.RNA.Protein.m$variable <- factor(Cell.DNA.RNA.Protein.m$variable, 
                                          levels = c("ProteinDifference", "RNADifference", "DNA_ratio"))


ggplot(data= Cell.DNA.RNA.Protein.m, 
       aes(x=variable, y=Chrm_Arm))+ #3x10
  geom_raster(aes(fill = value), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference") + 
  xlab("8505C")+
  ylab("Chromosome arm")+
  coord_flip()+ #flip x and Y axis to arm on x
  ggtitle("DNA ratio, difference in RNA and Protein expression \n per chromosome arm in 8505C, ACH-001307")
# 8x3
# plot.Heatmap.8505C.DNA.RNA.Protein_min10cell

####
# ACH-000755 = HCC2218
DNA.ACH000019<-subset(An_Cell_data, DepMap_ID == "ACH-000755")
RNA.ACH000019<-subset(RNA.Diff.byChromosome, Cell_ID == "ACH-000755")
protein.ACH000019<-subset(protein.Diff.byChromosome, Cell_ID == "ACH-000755")

Cell.DNA.RNA<-merge(x=DNA.ACH000019, y= RNA.ACH000019, 
                    by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                    sort = TRUE)# 
Cell.DNA.RNA.Protein<-merge(x=Cell.DNA.RNA, y= protein.ACH000019, 
                            by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                            sort = TRUE)# 

Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein %>% 
  rename( #new = old
    RNADifference = Difference.x,
    ProteinDifference = Difference.y
  )

Cell.DNA.RNA.Protein$arm_call<-Cell.DNA.RNA.Protein$arm_call/2
Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein[ , which(names(Cell.DNA.RNA.Protein) 
                                                     %in% c("Chrm_Arm", "RNADifference", "ProteinDifference", "DNA_ratio"))]
Cell.DNA.RNA.Protein.m <- melt(Cell.DNA.RNA.Protein)

Cell.DNA.RNA.Protein.m$variable <- factor(Cell.DNA.RNA.Protein.m$variable, 
                                          levels = c("ProteinDifference", "RNADifference", "DNA_ratio"))


ggplot(data= Cell.DNA.RNA.Protein.m, 
       aes(x=variable, y=Chrm_Arm))+ #3x10
  geom_raster(aes(fill = value), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference") + 
  xlab("HCC2218")+
  ylab("Chromosome arm")+
  coord_flip()+ #flip x and Y axis to arm on x
  ggtitle("DNA ratio, difference in RNA and Protein expression \n per chromosome arm in HCC2218, ACH-000755")
# 8x3
# plot.Heatmap.HCC2218.DNA.RNA.Protein_min10cell

####
# ACH-000589 = NCIH1437
DNA.ACH000019<-subset(An_Cell_data, DepMap_ID == "ACH-000589")
RNA.ACH000019<-subset(RNA.Diff.byChromosome, Cell_ID == "ACH-000589")
protein.ACH000019<-subset(protein.Diff.byChromosome, Cell_ID == "ACH-000589")

Cell.DNA.RNA<-merge(x=DNA.ACH000019, y= RNA.ACH000019, 
                    by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                    sort = TRUE)# 
Cell.DNA.RNA.Protein<-merge(x=Cell.DNA.RNA, y= protein.ACH000019, 
                            by.x="Chrm_Arm", by.y="Chrm.Num.Arm", 
                            sort = TRUE)# 

Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein %>% 
  rename( #new = old
    RNADifference = Difference.x,
    ProteinDifference = Difference.y
  )

Cell.DNA.RNA.Protein$arm_call<-Cell.DNA.RNA.Protein$arm_call/2
Cell.DNA.RNA.Protein<- Cell.DNA.RNA.Protein[ , which(names(Cell.DNA.RNA.Protein) 
                                                     %in% c("Chrm_Arm", "RNADifference", "ProteinDifference", "DNA_ratio"))]
Cell.DNA.RNA.Protein.m <- melt(Cell.DNA.RNA.Protein)

Cell.DNA.RNA.Protein.m$variable <- factor(Cell.DNA.RNA.Protein.m$variable, 
                                          levels = c("ProteinDifference", "RNADifference", "DNA_ratio"))


ggplot(data= Cell.DNA.RNA.Protein.m, 
       aes(x=variable, y=Chrm_Arm))+ #3x10
  geom_raster(aes(fill = value), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.7, 0.7)) + 
  xlab("NCIH1437")+
  ylab("Chromosome arm")+
  coord_flip()+ #flip x and Y axis to arm on x
  ggtitle("DNA ratio, difference in RNA and Protein expression \n per chromosome arm in NCIH1437, ACH-000589")
# 8x3
# plot.Heatmap.NCIH1437.DNA.RNA.Protein_min10cell

