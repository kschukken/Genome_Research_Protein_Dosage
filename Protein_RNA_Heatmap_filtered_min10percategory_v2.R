###### Plot heatmaps of protein and RNA difference upon chrm arm gain and loss ####
### 210128
## Update- 210201- only use 371 cell lines with Protein and RNA data
## Update 210427- only use genes that have 10+ cells per category, for both RNA and protein
## Update 210502 - only use genes that have 10+ cells/category (~9k genes), 

### by Klaske Schukken
### Graphs for Depmap Protein & RNA data: 
### Graph: heatmap of average change in RNA/Protein experssion per chromosome, per chrm gain/loss

library('ggplot2')
library('tidyverse')
library('xlsx')
library(readxl)
library(reshape2)
library('BBmisc')

###### Step 1: Get data ####
### Step 1: Get data. 
## Get  Protein  expression data from depmap.org and cell info data from depmap

## Import Protein expression data and uri ben-david lab, Cohen-Sharir nature 2021 chromosome arm data
###
setwd()# ! set working directory to correct location

Protein_Expression.1<-read_excel("mmc2.xlsx", 
                                 sheet= "Normalized Protein Expression")
Protein_Expression.1<-Protein_Expression.1[,c(1:426)] #delete empty collumns
#Name dataframe with protein_ID and corresponding Protein name, this way we have both gene ID and Protein name
Protein_ProID<-Protein_Expression.1[,c(1:2, 6)] 
Protein_ProID$Protein_Id<- str_replace_all(Protein_ProID$Protein_Id, "[|]", ".")
Protein_ProID$Protein_Id<- str_replace_all(Protein_ProID$Protein_Id, "[-]", ".")


###Get updated Protein info with only 371 cell lines (cell lines also have RNA Data)
# from R file:  Protein_RNA_filtered_CellLine.R
Protein.Expression.filtered<-read.delim2("Protein_Expression_filtered.csv", 
                                dec=",", header = TRUE, sep=";")

###
aneuploid<- read.csv("arm_calls_data_bendavid.csv", header=TRUE)

#Protein Info. names, chrm location, etc. 
Protein_Info<-read_excel("HGNC_Protein_Info_edited.xlsx", sheet ="HGNC names")
Protein_Info2<-Protein_Info[,c(2,3,11,12,22)]

###
# Import RNA expression data (no aneuploidy data). CCLE depmap data.
# 1 303 cell lines
# 19 145 RNAs

#RNA Info. names, chrm location, etc. 
RNA_Info<-read_excel("HGNC_Protein_Info_edited.xlsx", sheet ="HGNC names")
RNA_Info2<-RNA_Info[,c(2,3,11,12,22,14)]

###Get updated Protein info with only 371 cell lines (cell lines also have RNA Data)
# from R file:  Protein_RNA_filtered_CellLine.R
RNA.Expression.filtered<-read.delim2("RNA_Expression_filtered.csv", 
                                         dec=",", header = TRUE, sep=";")


### Now filter Protein/RNA expression for genes that have a minimum of 10 cells per 
## category (gain, neutral, loss; in both RNA and protein): 
## 9414 genes in filtered data set 
## filter (10+ cells/category dataset from: Protein_RNA_Corr_min10.Filtered.R)
# CN.Diff.xRNA.yProt.ThreeGroups
# this is from R file: Protein_RNA_expression.PerCell_v2.R
CN.Diff.xRNA.yProt.ThreeGroups<- read.csv("RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points_3cat.csv")


# Proteins: get collumns that have genes with 10+ cells/condition. 
#collumn names = "sp.Q9NQ94.A1CF_HUMAN" Protein ID format. same as CN.Diff.xRNA.yProt$Protein_ID
# select all collumns that have same ID as in filtered list. 
Protein.Expression.filtered_min10Cells<- Protein.Expression.filtered %>% select(one_of(CN.Diff.xRNA.yProt$Protein_ID))
# Now add Cell_lines back in. 
Protein.Expression.filtered_min10Cells$Cell_line<- Protein.Expression.filtered$Cell_line

# RNA: get collumns that have genes with 10+ cells/condition. 
#collumn names = "TSPAN6..7105." RNA ID format. same as CN.Diff.xRNA.yProt$RNA_ID
# select all collumns that have same ID as in filtered list. 
# length(unique(CN.Diff.xRNA.yProt$RNA_ID)), there are 9094 unique RNAs, 
# (9413 unique proteins, can have multiple protien isoforms per RNA)
RNA.Expression.filtered_min10Cells<- RNA.Expression.filtered %>% select(one_of(CN.Diff.xRNA.yProt$RNA_ID))
# Now add Cell_lines back in. 
RNA.Expression.filtered_min10Cells$Cell_line<- RNA.Expression.filtered$Cell_line


###  Optional: Get data for filtered protein expression data
## Only data for Non-Small Cell Lung Cancer cells
## Only for genes that have top 90% RNA expression genes
## only for genes in cells where that gene is not mutated
## data generated in Protein_RNA_expression.PerCell_v2.R

Protein.Expression_Lung_highRNA_noMut<-read.delim2("Protein_Expression_Lung_highRNA_noMut.csv", 
                                     dec=".", header = TRUE, sep=",")
RNA.Expression_Lung_highRNA_noMut<-read.delim2("RNA_Expression_Lung_highRNA_noMut.csv", 
                                                   dec=".", header = TRUE, sep=",")


###### Step 2: Heatmap for PROTEIN Chromosome GAIN ####
### Step 2: make heatmap function for PROTEIN Chromosome GAIN
### Using only 371 cell lines that have both RNA and Protein data
## Substep: 1) find cells with gain & nogain of chromosome X, arm p/q
## substep: 2) Find difference in protein expression per protein
## substep: 3) find location of each protein, chromosome and arm location
## substep: 4) find mean of difference in Prot exp per chromosome arm--> into datadrame
## substep: 5) repeat for all chromosome arms.
## substep: 6) plot mean change in prot per chromosome, per chromosome gain. 

#Set up data
AllChrmArms<- data.frame(Chrm=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,14,15,16,16,17,17,18,18,19,19,20,20,21,22), 
                         Arm=c("p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q",
                               "q","q","q","p","q","p","q","p","q","p","q","p","q","q","q"))
Protein.Gain.Diff.byChromosome<-data.frame(Chrm.Num.Arm=as.character(), Difference=as.numeric(), Chrm.Gained=as.character())

Protein_Info3<-Protein_Info2
Protein_Info3$arm <- gsub('[0-9]+', '', Protein_Info3$'Chromosome band') #Find test protein chromosome arm
Protein_Info3$arm <- gsub('[.]', '', Protein_Info3$arm) #also remove period, if needed
Protein_Info3$arm <- str_sub(Protein_Info3$arm, -1, -1) #start & end on last character. get only p or q
#filter cell data and get only aneuploidy scores for chrm arm of test protein location


###Run data below
##Repeat task for all chromosome arm gains:
for (i in 1:length(AllChrmArms$Chrm)) {
  TestChrm<-AllChrmArms$Chrm[i]
  TestArm<-AllChrmArms$Arm[i]
  
  #subset 1: Find cells with gain or no gain of chromosome X
  testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
  testchrm.percell <- filter(testchrm.percell, arm==TestArm)
  
  Cells.Tri.forchrm<- filter(testchrm.percell, arm_call==1) #all cells triploid for testProt chrm
  Cells.Di.forchrm<- filter(testchrm.percell, arm_call==0) #all cells diploid for testProt chrm
  
  #Substep2: get difference in protein expression between diploid & triploid cells
  #Now get list of protein expression in cells trisomic and disomic for each chrm arm
  Protein.Tri.Cells<-merge(y= Protein.Expression.filtered_min10Cells, x= Cells.Tri.forchrm, 
                           by.y="Cell_line", by.x="DepMap_ID", 
                           sort = TRUE)# 368 cells, 12757 proteins 
  Protein.Di.Cells<-merge(y= Protein.Expression.filtered_min10Cells, x= Cells.Di.forchrm, 
                          by.y="Cell_line", by.x="DepMap_ID", 
                          sort = TRUE)# 368 cells, 12757 proteins 
  
  # for each protein, find mean Tri, mean Di, difference
  Diff.PerProtein<- data.frame(Protein_ID=as.character(), #Set up data.frame with 3 collumn
                               Triploid_Exp=numeric(), 
                               Diploid_Exp=numeric(), 
                               Difference=numeric(), 
                               stringsAsFactors = TRUE) #this is needed to prevent errors
  
  for (w in 8:9420) { #Get difference in exp between diploid & triploid cells, per protein
    Diff.PerProtein<- rbind(Diff.PerProtein, 
                            data.frame(Protein_ID=colnames(Protein.Tri.Cells[w]),
                                       Triploid_Exp=mean(Protein.Tri.Cells[,w], na.rm=TRUE),
                                       Diploid_Exp=mean(Protein.Di.Cells[,w], na.rm=TRUE),
                                       Difference=mean(Protein.Tri.Cells[,w], na.rm=TRUE) - mean(Protein.Di.Cells[,w], na.rm=TRUE)))
    }
  
  # then link Protein values with location 
  
  ## Substep 3: add chromosome location to each gene
  # add more protein info per protein. Add Uniprot_Acc and Gene_Symbol
  Diff.PerProtein2<-merge(x= Diff.PerProtein, y= Protein_ProID, 
                          by.x="Protein_ID", by.y="Protein_Id", 
                          sort = TRUE)# 5 collumns, 12755 genes
  #Find genes with same gene_Symbol
  Diff.PerProtein3<-merge(x= Diff.PerProtein2, y= Protein_Info2, 
                          by.x="Gene_Symbol", by.y="Protein_Name", 
                          sort = TRUE)# 10 collumns, 12467 genes
  #Of genes without matching gene symbol...
  No_Gene_Symbol <- anti_join(Diff.PerProtein2, Protein_Info2, #finding genes_Symbol without match
                              by = c("Gene_Symbol" = "Protein_Name"))
  #...find the genes with matching uniprot_ID
  Diff.PerProtein4<-merge(x= No_Gene_Symbol, y= Protein_Info2, 
                          by.x="Uniprot_Acc", by.y="UniProt ID(supplied by UniProt)", 
                          sort = TRUE)# 10 collumns, 253 genes, only no gene-symbol genes with uniprot_ID
  # Merge genes with gene symbol and genes with only uniprot ID. 
  Diff.PerProtein5<-merge(x= Diff.PerProtein3, y= Diff.PerProtein4, 
                          all=TRUE)# 11 collumns, 12720 genes
  
  
  # make Chrm number & arm categories, and group by chrm num/arm
  Diff.PerProtein5$arm <- gsub('[0-9]+', '', Diff.PerProtein5$'Chromosome band')
  Diff.PerProtein5$arm <- gsub('[.]', '', Diff.PerProtein5$arm)
  Diff.PerProtein5$arm <- str_sub(Diff.PerProtein5$arm, -1, -1) #start & end on character 1. get only p or q
  
  Diff.PerProtein5$Chrm.Num.Arm<- paste0(Diff.PerProtein5$Chromosome, Diff.PerProtein5$arm)
  
  # Remove bad chromosome locations. ex "mitochondria m"--> not for this analysis
  Diff.PerProtein6<-subset(Diff.PerProtein5, Chrm.Num.Arm != "mitochondriaa")
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "NANA")
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "reservedd")
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "2cen-q")
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "6s") #no idea what "s" arm is...
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "7s")
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "22p") #only 1 gene
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "22s") #?
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "21p") #only 1 gene
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "Yp") #?
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "Yq") #only 1 gene
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "NA") 
  Diff.PerProtein6$Chrm.Num.Arm<-as.factor(Diff.PerProtein6$Chrm.Num.Arm) #as factor. 
  
  ##
  # Get mean protein expression difference by chrm num/arm category
  Mean.Diff<-aggregate( Difference ~ Chrm.Num.Arm, Diff.PerProtein6, mean ) 
  Mean.Diff$Chrm.Gained<-paste0(TestChrm, TestArm)
  
  Protein.Gain.Diff.byChromosome<-rbind(Protein.Gain.Diff.byChromosome, Mean.Diff)
  print(paste0("Finished with chrom arm ", TestChrm, TestArm))
}

### Now plot the heatmap for chromosome Protein chrm Gain changes. 
length(Protein.Gain.Diff.byChromosome$Chrm.Num.Arm)#1599, after filtering for 10+ genes/category still 1599
write.csv(Protein.Gain.Diff.byChromosome, 
          file= "Protein.Gain.Diff.PerChromosome.Filtered.Min10cells2.csv")
#setwd()#set working directory 
#Protein.Gain.Diff.byChromosome<-read.delim2("Protein.Gain.Diff.PerChromosome.Filtered.csv", 
#                                     dec=".", header = TRUE, sep=",")


Protein.Gain.Diff.byChromosome$Chrm.Num.Arm <- factor(Protein.Gain.Diff.byChromosome$Chrm.Num.Arm, 
                                                      levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                                 "4p", "4q", "5p", "5q", "6p", "6q",
                                                                 "7p", "7q", "8p", "8q", "9p", "9q",
                                                                 "10p", "10q", "11p", "11q", "12p", "12q",
                                                                 "13p", "13q", "14p", "14q", "15p", "15q",
                                                                 "16p", "16q", "17p", "17q", "18p", "18q",
                                                                 "19p", "19q", "20p", "20q", "21p", "21q",
                                                                 "22p", "22q", "Xp", "Xq"))
Protein.Gain.Diff.byChromosome$Chrm.Gained <- factor(Protein.Gain.Diff.byChromosome$Chrm.Gained, 
                                                     levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                                "4p", "4q", "5p", "5q", "6p", "6q",
                                                                "7p", "7q", "8p", "8q", "9p", "9q",
                                                                "10p", "10q", "11p", "11q", "12p", "12q",
                                                                "13p", "13q", "14p", "14q", "15p", "15q",
                                                                "16p", "16q", "17p", "17q", "18p", "18q",
                                                                "19p", "19q", "20p", "20q", "21p", "21q",
                                                                "22p", "22q", "Xp", "Xq"))


#Protein.Gain.Diff.byChromosome<-read.delim2("Protein.Gain.Diff.PerChromosome.Filtered.csv", 
#                                            dec=".", header = TRUE, sep=",")

ggplot(Protein.Gain.Diff.byChromosome, aes(x=Chrm.Num.Arm, y=Chrm.Gained))+ 
  geom_raster(aes(fill = Difference), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  xlab("Protein expression per chromosome arm")+
  ylab("Chromosome arm gained")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.4, 0.4)) +
  coord_flip()+
  ggtitle("Difference in protein expresssion upon \n chromosome gain")
#7x6
# plot.heatmap.Protein.Gain.DiffperArm_10cellmin_dodgerblue
#sky blue1 or dodgerblue3

###### Step 3: Heatmap for PROTEIN Chromosome LOSS #### 
### Step 3: make heatmap function for PROTEIN Chromosome LOSS
## Substep: 1) Find cells with loss & noloss of chromosome X, arm p/q
## substep: 2) Find difference in protein expression per protein
## substep: 3) Find location of each protein, chromosome and arm location
## substep: 4) Find mean of difference in Prot exp per chromosome arm--> into datadrame
## substep: 5) Repeat for all chromosome arms.
## substep: 6) Plot mean change in prot per chromosome, per chromosome loss 

#Set up data
AllChrmArms<- data.frame(Chrm=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,14,15,16,16,17,17,18,18,19,19,20,20,21,22), 
                         Arm=c("p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q",
                               "q","q","q","p","q","p","q","p","q","p","q","p","q","q","q"))
Protein.Loss.Diff.byChromosome<-data.frame(Chrm.Num.Arm=as.character(), Difference=as.numeric(), Chrm.Lost=as.character())

Protein_Info3<-Protein_Info2
Protein_Info3$arm <- gsub('[0-9]+', '', Protein_Info3$'Chromosome band') #Find test protein chromosome arm
Protein_Info3$arm <- gsub('[.]', '', Protein_Info3$arm) #also remove period, if needed
Protein_Info3$arm <- str_sub(Protein_Info3$arm, -1, -1) #start & end on last character. get only p or q
#filter cell data and get only aneuploidy scores for chrm arm of test protein location


###Run data below PROTEIN LOSS
##Repeat task for all chromosome arm losses:
for (i in 1:length(AllChrmArms$Chrm)) {
  TestChrm<-AllChrmArms$Chrm[i]
  TestArm<-AllChrmArms$Arm[i]
  
  #subset 1: Find cells with loss or no loss of chromosome X
  testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
  testchrm.percell <- filter(testchrm.percell, arm==TestArm)
  
  Cells.Mono.forchrm<- filter(testchrm.percell, arm_call==-1) #all cells monoploid for testProt chrm
  Cells.Di.forchrm<- filter(testchrm.percell, arm_call==0) #all cells diploid for testProt chrm
  
  #Substep2: get difference in protein expression between diploid & monoploid cells
  #Now get list of protein expression in cells monosomic and disomic for each chrm arm
  Protein.Mono.Cells<-merge(y= Protein.Expression.filtered_min10Cells, x= Cells.Mono.forchrm, 
                            by.y="Cell_line", by.x="DepMap_ID", 
                            sort = TRUE)# 368 cells, 12757 proteins 
  Protein.Di.Cells<-merge(y= Protein.Expression.filtered_min10Cells, x= Cells.Di.forchrm, 
                          by.y="Cell_line", by.x="DepMap_ID", 
                          sort = TRUE)# 368 cells, 12757 proteins 
  
  # for each protein, find mean mono, mean Di, difference
  Diff.PerProtein<- data.frame(Protein_ID=as.character(), #Set up data.frame with 3 collumn
                               Monoploid_Exp=numeric(), 
                               Diploid_Exp=numeric(), 
                               Difference=numeric(), 
                               stringsAsFactors = TRUE) #this is needed to prevent errors
  
  for (w in 8:9420) {#Get difference in exp between diploid & monoploid cells, per protein
    Diff.PerProtein<- rbind(Diff.PerProtein, 
                            data.frame(Protein_ID=colnames(Protein.Mono.Cells[w]),
                                       Monoploid_Exp=mean(Protein.Mono.Cells[,w], na.rm=TRUE),
                                       Diploid_Exp=mean(Protein.Di.Cells[,w], na.rm=TRUE),
                                       Difference=mean(Protein.Mono.Cells[,w], na.rm=TRUE) - mean(Protein.Di.Cells[,w], na.rm=TRUE)))
    } 
  # then link Protein values with location 
  
  ## Substep 3: add chromosome location to each gene
  #add more protein info per protein. Add Uniprot_Acc and Gene_Symbol
  Diff.PerProtein2<-merge(x= Diff.PerProtein, y= Protein_ProID, 
                          by.x="Protein_ID", by.y="Protein_Id", 
                          sort = TRUE)# 5 collumns, 12755 genes
  #Find genes with same gene_Symbol
  Diff.PerProtein3<-merge(x= Diff.PerProtein2, y= Protein_Info2, 
                          by.x="Gene_Symbol", by.y="Protein_Name", 
                          sort = TRUE)# 10 collumns, 12467 genes
  #Of genes without matching gene symbol...
  No_Gene_Symbol <- anti_join(Diff.PerProtein2, Protein_Info2, #finding genes_Symbol without match
                              by = c("Gene_Symbol" = "Protein_Name"))
  #...find the genes with matching uniprot_ID
  Diff.PerProtein4<-merge(x= No_Gene_Symbol, y= Protein_Info2, 
                          by.x="Uniprot_Acc", by.y="UniProt ID(supplied by UniProt)", 
                          sort = TRUE)# 10 collumns, 253 genes, only no gene-symbol genes with uniprot_ID
  # Merge genes with gene symbol and genes with only uniprot ID. 
  Diff.PerProtein5<-merge(x= Diff.PerProtein3, y= Diff.PerProtein4, 
                          all=TRUE)# 11 collumns, 12720 genes
  
  
  # make Chrm number & arm categories, and group by chrm num/arm
  Diff.PerProtein5$arm <- gsub('[0-9]+', '', Diff.PerProtein5$'Chromosome band')
  Diff.PerProtein5$arm <- gsub('[.]', '', Diff.PerProtein5$arm)
  Diff.PerProtein5$arm <- str_sub(Diff.PerProtein5$arm, -1, -1) #start & end on character 1. get only p or q
  
  Diff.PerProtein5$Chrm.Num.Arm<- paste0(Diff.PerProtein5$Chromosome, Diff.PerProtein5$arm)
  
  # Remove bad chromosome locations. ex "mitochondria m"--> not for this analysis
  Diff.PerProtein6<-subset(Diff.PerProtein5, Chrm.Num.Arm != "mitochondriaa")
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "NANA")
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "reservedd")
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "2cen-q")
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "6s") #no idea what "s" arm is...
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "7s")
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "22p") #only 1 gene
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "22s") #?
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "21p") #only 1 gene
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "Yp") #?
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "Yq") #only 1 gene
  Diff.PerProtein6<-subset(Diff.PerProtein6, Chrm.Num.Arm != "NA") 
  Diff.PerProtein6$Chrm.Num.Arm<-as.factor(Diff.PerProtein6$Chrm.Num.Arm) #as factor. 
  
  ##
  # Get mean protein expression difference by chrm num/arm category
  Mean.Diff<-aggregate( Difference ~ Chrm.Num.Arm, Diff.PerProtein6, mean ) 
  Mean.Diff$Chrm.Lost<-paste0(TestChrm, TestArm)
  
  Protein.Loss.Diff.byChromosome<-rbind(Protein.Loss.Diff.byChromosome, Mean.Diff)
  print(paste0("Finished chromosome ", TestChrm, TestArm))
}

### Now plot the heatmap for chromosome Protein chrm loss changes. 
length(Protein.Loss.Diff.byChromosome$Chrm.Num.Arm)#1599, after sorting for 10+ cells/gene: 
#setwd()
write.csv(Protein.Loss.Diff.byChromosome, file= "Protein.Loss.Diff.PerChromosome_filtered_Min10cells2.csv")
#setwd()
#Protein.Loss.Diff.byChromosome<-read.delim2("Protein.Loss.Diff.PerChromosome.csv", 
#                                        dec=",", header = TRUE, sep=";")
#8:30 to

Protein.Loss.Diff.byChromosome$Chrm.Num.Arm <- factor(Protein.Loss.Diff.byChromosome$Chrm.Num.Arm, 
                                                      levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                                 "4p", "4q", "5p", "5q", "6p", "6q",
                                                                 "7p", "7q", "8p", "8q", "9p", "9q",
                                                                 "10p", "10q", "11p", "11q", "12p", "12q",
                                                                 "13p", "13q", "14p", "14q", "15p", "15q",
                                                                 "16p", "16q", "17p", "17q", "18p", "18q",
                                                                 "19p", "19q", "20p", "20q", "21p", "21q",
                                                                 "22p", "22q", "Xp", "Xq"))
Protein.Loss.Diff.byChromosome$Chrm.Lost <- factor(Protein.Loss.Diff.byChromosome$Chrm.Lost, 
                                                   levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                              "4p", "4q", "5p", "5q", "6p", "6q",
                                                              "7p", "7q", "8p", "8q", "9p", "9q",
                                                              "10p", "10q", "11p", "11q", "12p", "12q",
                                                              "13p", "13q", "14p", "14q", "15p", "15q",
                                                              "16p", "16q", "17p", "17q", "18p", "18q",
                                                              "19p", "19q", "20p", "20q", "21p", "21q",
                                                              "22p", "22q", "Xp", "Xq"))

#setwd()#set working directory to depmap
#Protein.Loss.Diff.byChromosome<-read.delim2("Protein.Loss.Diff.PerChromosome_filtered_Min10cells.csv", 
#                                            dec=".", header = TRUE, sep=",")


ggplot(Protein.Loss.Diff.byChromosome, aes(x=Chrm.Num.Arm, y=Chrm.Lost))+ 
  geom_raster(aes(fill = Difference), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  xlab("Protein expression per chromosome arm")+
  ylab("Chromosome arm Lost")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.4, 0.4)) +
  coord_flip()+
  ggtitle("Difference in protein expresssion upon \n chromosome loss")
# 7x6
# plot.heatmap.Protein.Loss.DiffperArm_10cellmin_dodgerblue


###### Step 4: Heatmap for RNA Chromosome GAIN ####
### Step 2: make heatmap function for RNA Chromosome GAIN
## Substep: 1) find cells with gain & nogain of chromosome X, arm p/q
## substep: 2) Find difference in RNA expression per RNA
## substep: 3) find location of each RNA, chromosome and arm location
## substep: 4) find mean of difference in RNA exp per chromosome arm--> into datadrame
## substep: 5) repeat for all chromosome arms.
## substep: 6) plot mean change in RNA per chromosome, per chromosome gain. 

#Set up data
AllChrmArms<- data.frame(Chrm=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,14,15,16,16,17,17,18,18,19,19,20,20,21,22), 
                         Arm=c("p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q",
                               "q","q","q","p","q","p","q","p","q","p","q","p","q","q","q"))
RNA.Gain.Diff.byChromosome<-data.frame(Chrm.Num.Arm=as.character(), 
                                       Difference=as.numeric(), 
                                       Chrm.Gained=as.character())

RNA_Info3<-RNA_Info2
RNA_Info3$arm <- gsub('[0-9]+', '', RNA_Info3$'Chromosome band') #Find test RNA chromosome arm
RNA_Info3$arm <- gsub('[.]', '', RNA_Info3$arm) #also remove period, if needed
RNA_Info3$arm <- str_sub(RNA_Info3$arm, -1, -1) #start & end on last character. get only p or q
#filter cell data and get only aneuploidy scores for chrm arm of test RNA location

###Run data below
##Repeat task for all chromosome arm gains:
for (i in 1:length(AllChrmArms$Chrm)) {
  TestChrm<-AllChrmArms$Chrm[i]
  TestArm<-AllChrmArms$Arm[i]
  
  #subset 1: Find cells with gain or no gain of chromosome X
  testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
  testchrm.percell <- filter(testchrm.percell, arm==TestArm)
  
  Cells.Tri.forchrm<- filter(testchrm.percell, arm_call==1) #all cells triploid for testRNA chrm
  Cells.Di.forchrm<- filter(testchrm.percell, arm_call==0) #all cells diploid for testRNA chrm
  
  #Substep2: get difference in RNA expression between diploid & triploid cells
  #Now get list of RNA expression in cells trisomic and disomic for each chrm arm
  RNA.Tri.Cells<-merge(y= RNA.Expression.filtered_min10Cells, x= Cells.Tri.forchrm, 
                       by.y="Cell_line", by.x="DepMap_ID", 
                       sort = TRUE)
  RNA.Di.Cells<-merge(y= RNA.Expression.filtered_min10Cells, x= Cells.Di.forchrm, 
                      by.y="Cell_line", by.x="DepMap_ID", 
                      sort = TRUE)
  
  # for each RNA, find mean Tri, mean Di, difference
  Diff.PerRNA<- data.frame(RNA_ID=as.character(), #Set up data.frame with 3 collumn
                           Triploid_Exp=numeric(), 
                           Diploid_Exp=numeric(), 
                           Difference=numeric(), 
                           stringsAsFactors = TRUE) #this is needed to prevent errors

  for (w in 8:9101) {
    #Get difference in exp between diploid & triploid cells, per RNA
    Diff.PerRNA<- rbind(Diff.PerRNA, 
                        data.frame(RNA_ID=colnames(RNA.Tri.Cells[w]),
                                   Triploid_Exp=mean(RNA.Tri.Cells[,w], na.rm=TRUE),
                                   Diploid_Exp=mean(RNA.Di.Cells[,w], na.rm=TRUE),
                                   Difference=mean(RNA.Tri.Cells[,w], na.rm=TRUE) - mean(RNA.Di.Cells[,w], na.rm=TRUE)))
    }
  # then link RNA values with location 
  
  Diff.PerRNA$RNA_Name<-sub("[..].*", "", as.character(Diff.PerRNA$RNA_ID))#12755 genes
  Diff.PerRNA$RNA_num<-sub(".*[.]{2}", "", as.character(Diff.PerRNA$RNA_ID))
  Diff.PerRNA$RNA_num<-sub("[.]", "", as.character(Diff.PerRNA$RNA_num))
  
  
  ## Substep 3: add chromosome location to each gene
  #add more RNA info per RNA. Add Uniprot_Acc and Gene_Symbol
  Diff.PerRNA2<-merge(x= Diff.PerRNA, y= RNA_Info2, 
                      by.x="RNA_Name", by.y="Protein_Name", 
                      sort = TRUE)# 5 collumns, 12755 genes
  #Of genes without matching gene symbol...
  No_Gene_Symbol <- anti_join(Diff.PerRNA, RNA_Info2, #finding genes_Symbol without match
                              by = c("RNA_Name" = "Protein_Name"))
  #...find the genes with matching uniprot_ID
  Diff.PerRNA4<-merge(x= No_Gene_Symbol, y= RNA_Info2, 
                      by.x="RNA_num", by.y="Entrez Gene ID", 
                      sort = TRUE)# 10 collumns, 253 genes, only no gene-symbol genes with uniprot_ID
  # Merge genes with gene symbol and genes with only uniprot ID. 
  Diff.PerRNA5<-merge(x= Diff.PerRNA2, y= Diff.PerRNA4, 
                      all=TRUE)# 11 collumns, 12720 genes
  
  
  # make Chrm number & arm categories, and group by chrm num/arm
  Diff.PerRNA5$arm <- gsub('[0-9]+', '', Diff.PerRNA5$'Chromosome band')
  Diff.PerRNA5$arm <- gsub('[.]', '', Diff.PerRNA5$arm)
  Diff.PerRNA5$arm <- str_sub(Diff.PerRNA5$arm, -1, -1) #start & end on character 1. get only p or q
  
  Diff.PerRNA5$Chrm.Num.Arm<- paste0(Diff.PerRNA5$Chromosome, Diff.PerRNA5$arm)
  
  # Remove bad chromosome locations. ex "mitochondria m"--> not for this analysis
  Diff.PerRNA6<-subset(Diff.PerRNA5, Chrm.Num.Arm != "mitochondriaa")
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "NANA")
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "reservedd")
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "2cen-q")
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "6s") #no idea what "s" arm is...
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "7s")
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "NA")
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "22p") #only 1 gene
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "14p") #only 1 gene
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "22s") #?
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "21p") #only 1 gene
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "Yp") #not enough Y
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "Yq") #not enough Y
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "NA") 
  Diff.PerRNA6$Chrm.Num.Arm<-as.factor(Diff.PerRNA6$Chrm.Num.Arm) #as factor. 
  
  ##
  # Get mean RNA expression difference by chrm num/arm category
  Mean.Diff<-aggregate( Difference ~ Chrm.Num.Arm, Diff.PerRNA6, mean ) 
  Mean.Diff$Chrm.Gained<-paste0(TestChrm, TestArm)
  
  RNA.Gain.Diff.byChromosome<-rbind(RNA.Gain.Diff.byChromosome, Mean.Diff)
  print(paste0("finished chromosome ", TestChrm, TestArm))
}

### Now plot the heatmap for chromosome RNA chrm Gain changes. 
length(RNA.Gain.Diff.byChromosome$Chrm.Gained) #1599 

setwd()
write.csv(RNA.Gain.Diff.byChromosome, 
           file= "RNA.Gain.Diff.PerChromosome_min10Cells2.csv")

#RNA.Gain.Diff.byChromosome<-read.delim2("RNA.Gain.Diff.PerChromosome_min10Cells.csv", 
#                                        dec=".", header = TRUE, sep=",")

RNA.Gain.Diff.byChromosome$Chrm.Num.Arm <- factor(RNA.Gain.Diff.byChromosome$Chrm.Num.Arm, 
                                                  levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                             "4p", "4q", "5p", "5q", "6p", "6q",
                                                             "7p", "7q", "8p", "8q", "9p", "9q",
                                                             "10p", "10q", "11p", "11q", "12p", "12q",
                                                             "13p", "13q", "14p", "14q", "15p", "15q",
                                                             "16p", "16q", "17p", "17q", "18p", "18q",
                                                             "19p", "19q", "20p", "20q", "21p", "21q",
                                                             "22p", "22q", "Xp", "Xq"))
RNA.Gain.Diff.byChromosome$Chrm.Gained <- factor(RNA.Gain.Diff.byChromosome$Chrm.Gained, 
                                                 levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                            "4p", "4q", "5p", "5q", "6p", "6q",
                                                            "7p", "7q", "8p", "8q", "9p", "9q",
                                                            "10p", "10q", "11p", "11q", "12p", "12q",
                                                            "13p", "13q", "14p", "14q", "15p", "15q",
                                                            "16p", "16q", "17p", "17q", "18p", "18q",
                                                            "19p", "19q", "20p", "20q", "21p", "21q",
                                                            "22p", "22q", "Xp", "Xq"))

ggplot(RNA.Gain.Diff.byChromosome, aes(x=Chrm.Num.Arm, y=Chrm.Gained))+ 
  geom_raster(aes(fill = Difference), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  xlab("RNA expression per chromosome arm")+
  ylab("Chromosome arm gained")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.4, 0.4)) +
  coord_flip()+
  ggtitle("Difference in RNA expresssion upon \n chromosome gain")
# 7X6
# plot.heatmap.RNA.Gain.DiffperArm_10cellmin_dodgerblue


###### Step 5: Heatmap for RNA Chromosome LOSS #### 
### Step 5: Make heatmap function for RNA Chromosome LOSS
## Substep: 1) Find cells with loss & noloss of chromosome X, arm p/q
## substep: 2) Find difference in RNA expression per RNA
## substep: 3) Find location of each RNA, chromosome and arm location
## substep: 4) Find mean of difference in RNA exp per chromosome arm--> into datadrame
## substep: 5) Repeat for all chromosome arms.
## substep: 6) Plot mean change in RNA per chromosome, per chromosome loss 

#Set up data
AllChrmArms<- data.frame(Chrm=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,14,15,16,16,17,17,18,18,19,19,20,20,21,22), 
                         Arm=c("p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q",
                               "q","q","q","p","q","p","q","p","q","p","q","p","q","q","q"))
RNA.Loss.Diff.byChromosome<-data.frame(Chrm.Num.Arm=as.character(), 
                                       Difference=as.numeric(), 
                                       Chrm.Lost=as.character())

RNA_Info3<-RNA_Info2
RNA_Info3$arm <- gsub('[0-9]+', '', RNA_Info3$'Chromosome band') #Find test RNA chromosome arm
RNA_Info3$arm <- gsub('[.]', '', RNA_Info3$arm) #also remove period, if needed
RNA_Info3$arm <- str_sub(RNA_Info3$arm, -1, -1) #start & end on last character. get only p or q
#filter cell data and get only aneuploidy scores for chrm arm of test RNA location


###Run data below
##Repeat task for all chromosome arm losses:
for (i in 1:length(AllChrmArms$Chrm)) {
  TestChrm<-AllChrmArms$Chrm[i]
  TestArm<-AllChrmArms$Arm[i]
  
  #subset 1: Find cells with loss or no loss of chromosome X
  testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
  testchrm.percell <- filter(testchrm.percell, arm==TestArm)
  
  Cells.Mono.forchrm<- filter(testchrm.percell, arm_call==-1) #all cells monoploid for testRNA chrm
  Cells.Di.forchrm<- filter(testchrm.percell, arm_call==0) #all cells diploid for testRNA chrm
  
  #Substep2: get difference in RNA expression between diploid & monoploid cells
  #Now get list of RNA expression in cells monosomy and disomic for each chrm arm
  RNA.Mono.Cells<-merge(y= RNA.Expression.filtered_min10Cells, x= Cells.Mono.forchrm, 
                        by.y="Cell_line", by.x="DepMap_ID", 
                        sort = TRUE)# 368 cells, 12757 RNAs 
  RNA.Di.Cells<-merge(y= RNA.Expression.filtered_min10Cells, x= Cells.Di.forchrm, 
                      by.y="Cell_line", by.x="DepMap_ID", 
                      sort = TRUE)# 368 cells, 12757 RNAs 
  
  # for each RNA, find mean mono, mean Di, difference
  Diff.PerRNA<- data.frame(RNA_ID=as.character(), #Set up data.frame with 3 collumn
                           Monoploid_Exp=numeric(), 
                           Diploid_Exp=numeric(), 
                           Difference=numeric(), 
                           stringsAsFactors = TRUE) #this is needed to prevent errors

  for (w in 8:9101) { #Get difference in exp between diploid & monoploid cells, per RNA
    Diff.PerRNA<- rbind(Diff.PerRNA, 
                        data.frame(RNA_ID=colnames(RNA.Mono.Cells[w]),
                                   Monoploid_Exp=mean(RNA.Mono.Cells[,w], na.rm=TRUE),
                                   Diploid_Exp=mean(RNA.Di.Cells[,w], na.rm=TRUE),
                                   Difference=mean(RNA.Mono.Cells[,w], na.rm=TRUE) - mean(RNA.Di.Cells[,w], na.rm=TRUE)))
    }
  
  # then link RNA values with location 
  
  Diff.PerRNA$RNA_Name<-sub("[..].*", "", as.character(Diff.PerRNA$RNA_ID))#12755 genes
  Diff.PerRNA$RNA_num<-sub(".*[.]{2}", "", as.character(Diff.PerRNA$RNA_ID))
  Diff.PerRNA$RNA_num<-sub("[.]", "", as.character(Diff.PerRNA$RNA_num))
  
  
  ## Substep 3: add chromosome location to each gene
  #add more RNA info per RNA. Add Uniprot_Acc and Gene_Symbol
  Diff.PerRNA2<-merge(x= Diff.PerRNA, y= RNA_Info2, 
                      by.x="RNA_Name", by.y="Protein_Name", 
                      sort = TRUE)# 5 collumns, 12755 genes
  #Of genes without matching gene symbol...
  No_Gene_Symbol <- anti_join(Diff.PerRNA, RNA_Info2, #finding genes_Symbol without match
                              by = c("RNA_Name" = "Protein_Name"))
  #...find the genes with matching uniprot_ID
  Diff.PerRNA4<-merge(x= No_Gene_Symbol, y= RNA_Info2, 
                      by.x="RNA_num", by.y="Entrez Gene ID", 
                      sort = TRUE)# 10 collumns, 253 genes, only no gene-symbol genes with uniprot_ID
  # Merge genes with gene symbol and genes with only uniprot ID. 
  Diff.PerRNA5<-merge(x= Diff.PerRNA2, y= Diff.PerRNA4, 
                      all=TRUE)# 11 collumns, 12720 genes
  
  
  # make Chrm number & arm categories, and group by chrm num/arm
  Diff.PerRNA5$arm <- gsub('[0-9]+', '', Diff.PerRNA5$'Chromosome band')
  Diff.PerRNA5$arm <- gsub('[.]', '', Diff.PerRNA5$arm)
  Diff.PerRNA5$arm <- str_sub(Diff.PerRNA5$arm, -1, -1) #start & end on character 1. get only p or q
  
  Diff.PerRNA5$Chrm.Num.Arm<- paste0(Diff.PerRNA5$Chromosome, Diff.PerRNA5$arm)
  
  # Remove bad chromosome locations. ex "mitochondria m"--> not for this analysis
  Diff.PerRNA6<-subset(Diff.PerRNA5, Chrm.Num.Arm != "mitochondriaa")
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "NANA")
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "reservedd")
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "2cen-q")
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "6s") #no idea what "s" arm is...
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "7s")
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "NA")
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "22p") #only 1 gene
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "14p") #only 1 gene
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "22s") #?
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "21p") #only 1 gene
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "Yp") #?
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "Yq") #only 1 gene
  Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "NA") 
  Diff.PerRNA6$Chrm.Num.Arm<-as.factor(Diff.PerRNA6$Chrm.Num.Arm) #as factor. 
  
  ##
  # Get mean RNA expression difference by chrm num/arm category
  Mean.Diff<-aggregate( Difference ~ Chrm.Num.Arm, Diff.PerRNA6, mean ) 
  Mean.Diff$Chrm.Lost<-paste0(TestChrm, TestArm)
  
  RNA.Loss.Diff.byChromosome<-rbind(RNA.Loss.Diff.byChromosome, Mean.Diff)
  print(paste0("Chromosome arm finished: ",TestChrm, TestArm))
}


### Now plot the heatmap for chromosome RNA chrm loss changes. 
length(RNA.Loss.Diff.byChromosome$Chrm.Num.Arm)#1599  

setwd()
write.csv(RNA.Loss.Diff.byChromosome, file= "RNA.Loss.Diff.PerChromosome_Min10cells2.csv")
#RNA.Loss.Diff.byChromosome<-read.csv("RNA.Loss.Diff.PerChromosome_Min10cells.csv")

RNA.Loss.Diff.byChromosome$Chrm.Num.Arm <- factor(RNA.Loss.Diff.byChromosome$Chrm.Num.Arm, 
                                                  levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                             "4p", "4q", "5p", "5q", "6p", "6q",
                                                             "7p", "7q", "8p", "8q", "9p", "9q",
                                                             "10p", "10q", "11p", "11q", "12p", "12q",
                                                             "13p", "13q", "14p", "14q", "15p", "15q",
                                                             "16p", "16q", "17p", "17q", "18p", "18q",
                                                             "19p", "19q", "20p", "20q", "21p", "21q",
                                                             "22p", "22q", "Xp", "Xq"))
RNA.Loss.Diff.byChromosome$Chrm.Lost <- factor(RNA.Loss.Diff.byChromosome$Chrm.Lost, 
                                               levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                          "4p", "4q", "5p", "5q", "6p", "6q",
                                                          "7p", "7q", "8p", "8q", "9p", "9q",
                                                          "10p", "10q", "11p", "11q", "12p", "12q",
                                                          "13p", "13q", "14p", "14q", "15p", "15q",
                                                          "16p", "16q", "17p", "17q", "18p", "18q",
                                                          "19p", "19q", "20p", "20q", "21p", "21q",
                                                          "22p", "22q", "Xp", "Xq"))

ggplot(RNA.Loss.Diff.byChromosome, aes(x=Chrm.Num.Arm, y=Chrm.Lost))+ 
  geom_raster(aes(fill = Difference), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  xlab("RNA expression per chromosome arm")+
  ylab("Chromosome arm lost")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.4, 0.4)) +
  coord_flip()+
  ggtitle("Difference in RNA expresssion upon \n chromosome loss")
# 7X6
# plot.heatmap.RNA.Loss.DiffperArm_10cellmin_dodgerblue

###### Step 6: Make heatmap for Protein & RNA with difference per chrm arm upon gain/loss of that arm ####
### Step 5: Make heatmap for Protein & RNA with diff per chrm arm upon gain/loss of arm
# Substep 1) Get data generated above: change in expression per arm upon each arm gain/loss
# Substep 2) Isolate difference in Arm X per gain/loss of arm X in RNA & Protein
# Substep 3) Combine RNA & Protein data. 
# Substep 4) Heatmap, plot data of RNA and Protein together

### heatmap 1: Gain of chromosome arm 
Protein.Gain.Diff.byChromosome
RNA.Gain.Diff.byChromosome

RNA.Protein.Gain.DiffInArm<- data.frame(RNAorPROTEIN=factor(), #Set up data.frame with 4 collumn
                         Chrm.Num.Arm=factor(), 
                         Difference=numeric(), 
                         Chrm.Gained=factor(), 
                         stringsAsFactors = TRUE) #this is needed to prevent errors

#Now add add Protein data: if gain chrm Xq, change in Chrm Xq protein expression
for (i in 1:length(Protein.Gain.Diff.byChromosome$Chrm.Num.Arm)) {
  if (Protein.Gain.Diff.byChromosome$Chrm.Num.Arm[i]==Protein.Gain.Diff.byChromosome$Chrm.Gained[i]) {
  RNA.Protein.Gain.DiffInArm<- rbind(RNA.Protein.Gain.DiffInArm, 
                                     data.frame(
                                       RNAorPROTEIN="Protein",
                                       Chrm.Num.Arm=Protein.Gain.Diff.byChromosome$Chrm.Num.Arm[i],
                                       Difference= Protein.Gain.Diff.byChromosome$Difference[i],
                                       Chrm.Gained= Protein.Gain.Diff.byChromosome$Chrm.Gained[i]
                                     ))
  }
}

#To same data frame, add RNA data. 
#Now add add RNA data: if gain chrm Xq, change in Chrm Xq RNA expression
for (i in 1:length(RNA.Gain.Diff.byChromosome$Chrm.Num.Arm)) {
  if (RNA.Gain.Diff.byChromosome$Chrm.Num.Arm[i]==RNA.Gain.Diff.byChromosome$Chrm.Gained[i]) {
    RNA.Protein.Gain.DiffInArm<- rbind(RNA.Protein.Gain.DiffInArm, 
                                       data.frame(
                                         RNAorPROTEIN="RNA",
                                         Chrm.Num.Arm=RNA.Gain.Diff.byChromosome$Chrm.Num.Arm[i],
                                         Difference= RNA.Gain.Diff.byChromosome$Difference[i],
                                         Chrm.Gained= RNA.Gain.Diff.byChromosome$Chrm.Gained[i]
                                       ))
  }
}

###Now plot Chromosome arm gain changes in RNA& Protein expression. heatmap. 
RNA.Protein.Gain.DiffInArm$RNAorPROTEIN <- factor(RNA.Protein.Gain.DiffInArm$RNAorPROTEIN, 
                                               levels = c("Protein", "RNA"))#plot RNA first

ggplot(RNA.Protein.Gain.DiffInArm, aes(x=RNAorPROTEIN, y=Chrm.Gained))+ # size 7x3?
  geom_raster(aes(fill = Difference), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  xlab("Difference in expression")+
  ylab("Chromosome arm gained")+
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  theme_classic()+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.4, 0.4)) +
  coord_flip()+
  ggtitle("Difference in expresssion upon \n chromosome gain")
# 8x3
# Plot.Protein.RNA.Gain.Difference_10cellMin_dodgerblue3

####
### heatmap 2: Loss of chromosome arm 
Protein.Loss.Diff.byChromosome
RNA.Loss.Diff.byChromosome

RNA.Protein.Loss.DiffInArm<- data.frame(RNAorPROTEIN=factor(), #Set up data.frame with 4 collumn
                                        Chrm.Num.Arm=factor(), 
                                        Difference=numeric(), 
                                        Chrm.Lost=factor(), 
                                        stringsAsFactors = TRUE) #this is needed to prevent errors

#Now add add Protein data: if gain chrm Xq, change in Chrm Xq protein expression
for (i in 1:length(Protein.Loss.Diff.byChromosome$Chrm.Num.Arm)) {
  if (Protein.Loss.Diff.byChromosome$Chrm.Num.Arm[i]==Protein.Loss.Diff.byChromosome$Chrm.Lost[i]) {
    RNA.Protein.Loss.DiffInArm<- rbind(RNA.Protein.Loss.DiffInArm, 
                                       data.frame(
                                         RNAorPROTEIN="Protein",
                                         Chrm.Num.Arm=Protein.Loss.Diff.byChromosome$Chrm.Num.Arm[i],
                                         Difference= Protein.Loss.Diff.byChromosome$Difference[i],
                                         Chrm.Lost= Protein.Loss.Diff.byChromosome$Chrm.Lost[i]
                                       ))
  }
}

#To same data frame, add RNA data. 
#Now add add RNA data: if Loss chrm Xq, change in Chrm Xq RNA expression
for (i in 1:length(RNA.Loss.Diff.byChromosome$Chrm.Num.Arm)) {
  if (RNA.Loss.Diff.byChromosome$Chrm.Num.Arm[i]==RNA.Loss.Diff.byChromosome$Chrm.Lost[i]) {
    RNA.Protein.Loss.DiffInArm<- rbind(RNA.Protein.Loss.DiffInArm, 
                                       data.frame(
                                         RNAorPROTEIN="RNA",
                                         Chrm.Num.Arm=RNA.Loss.Diff.byChromosome$Chrm.Num.Arm[i],
                                         Difference= RNA.Loss.Diff.byChromosome$Difference[i],
                                         Chrm.Lost= RNA.Loss.Diff.byChromosome$Chrm.Lost[i]
                                       ))
  }
}

length(RNA.Protein.Loss.DiffInArm$RNAorPROTEIN)
###Now plot Chromosome arm Loss changes in RNA& Protein expression. heatmap. 
RNA.Protein.Loss.DiffInArm$RNAorPROTEIN <- factor(RNA.Protein.Loss.DiffInArm$RNAorPROTEIN, 
                                                  levels = c("Protein", "RNA"))#plot RNA first

ggplot(RNA.Protein.Loss.DiffInArm, aes(x=RNAorPROTEIN, y=Chrm.Lost))+ 
  geom_raster(aes(fill = Difference), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  xlab("Difference in expression")+
  ylab("Chromosome arm lost")+
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="blakc"), 
        axis.text.y = element_text(color="black"))+
  theme_classic()+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.4, 0.4)) +
  coord_flip()+
  ggtitle("Difference in expresssion upon \n chromosome loss")
# 8x3
# Plot.Protein.RNA.Loss.Difference_10cellMin_dodgerblue
