##### Analysis of stable aneuploid cell proteomics ####
## Edited: 210407-- updated: Stingele et al. 2012 aneuploid proteomics data
## Edited: 210609-- cleaned up code for publication

## Author: Klaske Schukken
## Compare chromosome arm loss/gain with protein/RNA expression changes. 
## Data from Stingele et al. 2012

## Data is Supplementary table 1 from Stingele et al 2012
## Complete dataset combining quantitative genome (CGH), 
# transcriptome (microarrays) and prteome (SILAC) data of 
# trisomic and tetrasomic cell lines derived from HCT116 and RPE1. 
# The values are presented as log2 of the calculated aneuploid/diploid ratios. 

library('ggplot2')
library('tidyverse')
library(xlsx)
library(readxl)
library(reshape2)
library('BBmisc')

###
## plan: 
## Step 1: get Proteomics and RNA data (gain only, chrm 3,5,12 and 21)
##      it's log2(Fold Change) of expression in aneuploid cell vs control. 
##      limited number of cells. lots of proteins/RNAs
##      ony chromosome 5 has multiple cell line data-- only use chrm 5 data
## Step 2: look at RNA expression difference on ANEUPLOID CHROMOSOME only
##      Split into groups by aneuploid chromosome (ex chrm 5, chrm 12, chrm 21)
##      Chrm 5 has more cells than all the others, mean chrm 5 difference is more reliable than other data
## Step 3: look at Protein expression difference on ANEUPLOID CHROMOSOME only
##      Same as above, but now with proteins


# Notes: 
# I ended up only using chromosome 5 because multiple cell lines analysed had chrm 5 gain- 
# meaning I could take the mean of gene differences and remove some of the "one data point noise" 
# that other chromosomes displayed. 


##### Step 1: Get Data. Stingele RPE1, HCT116 aneuploidy data ####
setwd() #location of RNA and protein difference upon chrm gain or loss data
# get this dataset from supplementary data. or re-analyze using Protein_RNA_expression.PerCell_v2.R
# or from supplementary file 2
CN.Diff.xRNA.yProt.ThreeGroups<- read.csv(CN.Diff.xRNA.yProt.ThreeGroups, 
                                          file =paste("RNA.Protein_Mono.Di.Tri_Difference_Pvalue_min10points_3cat.csv", sep=','), 
                                          row.names = TRUE)


# Supplementary table 1. here I labeled it "Stingele.etal.2012.Genomics.xls" 
Stingele.RPE1<- read_xls("Stingele.etal.2012.Genomics.xls", sheet=3)
Stingele.HCT116<- read_xls("Stingele.etal.2012.Genomics.xls", sheet=2)

Stingele.info1<-Stingele.RPE1[,c(1,3,18)]
Stingele.info2<-Stingele.HCT116[,c(1,3,36)]
Stingele.info<-rbind(Stingele.info1, Stingele.info2)


## Divide data: get Protein data per cell line with Protein ID

Stingele.RPE1.Prot<- Stingele.RPE1[,c(1,18,8,14)]
Stingele.HCT116.Prot<- Stingele.HCT116[,c(1,20,26,30,34,36)]

Stingele.Prot <-merge(Stingele.RPE1.Prot, Stingele.HCT116.Prot, by.x="Protein IDs", by.y="Protein IDs", all=TRUE)
#10435 genes, ~2k of them shared between RPE1 and HCT116. ok
Stingele.Prot[,3]<-as.numeric(Stingele.Prot[,3]) #make numeric
Stingele.Prot[,4]<-as.numeric(Stingele.Prot[,4]) #make numeric
Stingele.Prot[,5]<-as.numeric(Stingele.Prot[,5]) #make numeric
Stingele.Prot[,6]<-as.numeric(Stingele.Prot[,6]) #make numeric
Stingele.Prot[,7]<-as.numeric(Stingele.Prot[,7]) #make numeric
Stingele.Prot[,8]<-as.numeric(Stingele.Prot[,8]) #make numeric


# Get RNA data. note: not all HCT116 cells that have Protein data also had RNA data
Stingele.RPE1.RNA<- Stingele.RPE1[,c(1,18,7,13)]
Stingele.HCT116.RNA<- Stingele.HCT116[,c(1,9,36)]

Stingele.RNA <-merge(Stingele.RPE1.RNA, Stingele.HCT116.RNA, by.x="Protein IDs", by.y="Protein IDs", all=TRUE)
#10435 genes, ~2k of them shared between RPE1 and HCT116. ok
Stingele.RNA$`mRNA RPE-1 21/3`<- as.numeric(Stingele.RNA$`mRNA RPE-1 21/3`) #make values numeric
Stingele.RNA$`mRNA RPE-1 5/3 12/3`<- as.numeric(Stingele.RNA$`mRNA RPE-1 5/3 12/3`) #make values numeric
Stingele.RNA$`mRNA HCT116 5/4`<- as.numeric(Stingele.RNA$`mRNA HCT116 5/4`) #make values numeric



##### Step 2: Now look at RNA on ANEUPLOID CHROMOSOME only #####
### Now look at RNA on ANEUPLOID CHROMOSOME only

## plot  chromosome 5 (mean diff per gene):
#  get genes on chrm 5, get mean RNA expression from cells with chrm 5 gain: 
#  only 2 cell lines with gain of chrm 5 in RNA dataset: 
Stingele.RNA.Chrm5<-subset(Stingele.RNA, Chromosome.x==5)
Stingele.RNA.Chrm5$meanDiff<-rowMeans(Stingele.RNA.Chrm5[,c(4:5)], na.rm=TRUE)


Stingele.RNA.Chrm5.info = merge(
  Stingele.RNA.Chrm5,
  Stingele.info[!duplicated(Stingele.info[, c("Protein IDs")]), ],#remove duplicated gene info
  by=c("Protein IDs","Protein IDs"),
  all.x=TRUE
) #add gene name, be sure to not add duplicates

### Combine data with depmap gain loss difference data
#   this is probably not the most efficent way to extract gene names from a list of names, but it works.

# Step 1: get gene names? (first name is list, then second)
Stingele.RNA.Chrm5.info$Gene_Name1<-sub(";.*", "", Stingele.RNA.Chrm5.info$"Gene Names")
# 2nd name: get first of first name in list, then get rid of all genes after "new first" gene to get second gene in name list
Stingele.RNA.Chrm5.info$Gene_Name2<-sub(".*?;", "", Stingele.RNA.Chrm5.info$"Gene Names")# get rid of first gene, "?" makes sub "lazy", finding first, not last occurance
Stingele.RNA.Chrm5.info$Gene_Name2<-sub(";.*", "", Stingele.RNA.Chrm5.info$Gene_Name2)
# 3rd and fourth names on list did not match with any missing data
Stingele.RNA.Chrm5.info$Gene_Name3<-sub(".*?;", "", Stingele.RNA.Chrm5.info$"Gene Names")# get rid of first gene
Stingele.RNA.Chrm5.info$Gene_Name3<-sub(".*?;", "", Stingele.RNA.Chrm5.info$Gene_Name3)# get rid of first gene
Stingele.RNA.Chrm5.info$Gene_Name3<-sub(";.*", "", Stingele.RNA.Chrm5.info$Gene_Name3)
# 4th and fourth names on list did not match with any missing data
Stingele.RNA.Chrm5.info$Gene_Name4<-sub(".*?;", "", Stingele.RNA.Chrm5.info$"Gene Names")# get rid of first gene
Stingele.RNA.Chrm5.info$Gene_Name4<-sub(".*?;", "", Stingele.RNA.Chrm5.info$Gene_Name4)# get rid of first gene
Stingele.RNA.Chrm5.info$Gene_Name4<-sub(".*?;", "", Stingele.RNA.Chrm5.info$Gene_Name4)# get rid of first gene
Stingele.RNA.Chrm5.info$Gene_Name4<-sub(";.*", "", Stingele.RNA.Chrm5.info$Gene_Name4)
# 5th and fourth names on list did not match with any missing data
Stingele.RNA.Chrm5.info$Gene_Name5<-sub(".*?;", "", Stingele.RNA.Chrm5.info$"Gene Names")# get rid of first gene
Stingele.RNA.Chrm5.info$Gene_Name5<-sub(".*?;", "", Stingele.RNA.Chrm5.info$Gene_Name5)# get rid of first gene
Stingele.RNA.Chrm5.info$Gene_Name5<-sub(".*?;", "", Stingele.RNA.Chrm5.info$Gene_Name5)# get rid of first gene
Stingele.RNA.Chrm5.info$Gene_Name5<-sub(".*?;", "", Stingele.RNA.Chrm5.info$Gene_Name5)# get rid of first gene
Stingele.RNA.Chrm5.info$Gene_Name5<-sub(";.*", "", Stingele.RNA.Chrm5.info$Gene_Name5)
# 6th and fourth names on list did not match with any missing data
Stingele.RNA.Chrm5.info$Gene_Name6<-sub(".*?;", "", Stingele.RNA.Chrm5.info$"Gene Names")# get rid of first gene
Stingele.RNA.Chrm5.info$Gene_Name6<-sub(".*?;", "", Stingele.RNA.Chrm5.info$Gene_Name6)# get rid of first gene
Stingele.RNA.Chrm5.info$Gene_Name6<-sub(".*?;", "", Stingele.RNA.Chrm5.info$Gene_Name6)# get rid of first gene
Stingele.RNA.Chrm5.info$Gene_Name6<-sub(".*?;", "", Stingele.RNA.Chrm5.info$Gene_Name6)# get rid of first gene
Stingele.RNA.Chrm5.info$Gene_Name6<-sub(".*?;", "", Stingele.RNA.Chrm5.info$Gene_Name6)# get rid of first gene
Stingele.RNA.Chrm5.info$Gene_Name6<-sub(";.*", "", Stingele.RNA.Chrm5.info$Gene_Name6)
# 7th and fourth names on list did not match with any missing data
Stingele.RNA.Chrm5.info$Gene_Name7<-sub(".*;", "", Stingele.RNA.Chrm5.info$"Gene Names")# get rid of first gene



## Step two: merge Dapmap data with Stingele data by gene name
## I have several collumns with difference names for the same gene, 
# if stingele geen name 1 = depmap gene symbol combine those two datapoints, else, test gene name 2, etc. 
# then combine all the data frames (game name 1 datset, gene name 2 dataset, etc. )

#merge with first name in list
Stingele.RNA.Chrm5.1<- merge(x=Stingele.RNA.Chrm5.info, #272 genes
                                    y=CN.Diff.xRNA.yProt.ThreeGroups, 
                                    by.x="Gene_Name1", 
                                    by.y="RNA_Name") #144 genes
# now get all non-merged cells, merge with second name in list
nomerge<- anti_join(x=Stingele.RNA.Chrm5.info, # 135 genes
                    y=CN.Diff.xRNA.yProt.ThreeGroups, 
                    by=c("Gene_Name1"="RNA_Name")) # genes
Stingele.RNA.Chrm5.2<- merge(x=nomerge, # genes
                                    y=CN.Diff.xRNA.yProt.ThreeGroups, 
                                    by.y="RNA_Name", 
                                    by.x="Gene_Name2") #65 genes
nomerge2<- anti_join(x=nomerge, # 45 genes
                    y=CN.Diff.xRNA.yProt.ThreeGroups, 
                    by=c("Gene_Name2"="RNA_Name")) 
Stingele.RNA.Chrm5.3<- merge(x=nomerge2, #1075 genes
                             y=CN.Diff.xRNA.yProt.ThreeGroups, 
                             by.y="RNA_Name", 
                             by.x="Gene_Name3") #65 genes
nomerge3<- anti_join(x=nomerge2, # 24 genes
                     y=CN.Diff.xRNA.yProt.ThreeGroups, 
                     by=c("Gene_Name3"="RNA_Name")) # genes
Stingele.RNA.Chrm5.4<- merge(x=nomerge3, #1075 genes
                             y=CN.Diff.xRNA.yProt.ThreeGroups, 
                             by.y="RNA_Name", 
                             by.x="Gene_Name4") #65 genes
nomerge4<- anti_join(x=nomerge3, # 23 genes
                     y=CN.Diff.xRNA.yProt.ThreeGroups, 
                     by=c("Gene_Name4"="RNA_Name")) # genes
Stingele.RNA.Chrm5.5<- merge(x=nomerge4, #1075 genes
                             y=CN.Diff.xRNA.yProt.ThreeGroups, 
                             by.y="RNA_Name", 
                             by.x="Gene_Name5") #65 genes
nomerge5<- anti_join(x=nomerge4, # 22 genes
                     y=CN.Diff.xRNA.yProt.ThreeGroups, 
                     by=c("Gene_Name5"="RNA_Name")) # genes
Stingele.RNA.Chrm5.6<- merge(x=nomerge5, # genes
                             y=CN.Diff.xRNA.yProt.ThreeGroups, 
                             by.y="RNA_Name", 
                             by.x="Gene_Name6") #22 genes
nomerge6<- anti_join(x=nomerge5, # 21 genes
                     y=CN.Diff.xRNA.yProt.ThreeGroups, 
                     by=c("Gene_Name6"="RNA_Name")) # 21 genes not found back in dataset


Stingele.RNA.Chrm5.merge<-rbind(Stingele.RNA.Chrm5.1, Stingele.RNA.Chrm5.2)
Stingele.RNA.Chrm5.merge<-rbind(Stingele.RNA.Chrm5.merge, Stingele.RNA.Chrm5.3)
Stingele.RNA.Chrm5.merge<-rbind(Stingele.RNA.Chrm5.merge, Stingele.RNA.Chrm5.4)
Stingele.RNA.Chrm5.merge<-rbind(Stingele.RNA.Chrm5.merge, Stingele.RNA.Chrm5.5)
Stingele.RNA.Chrm5.merge<-rbind(Stingele.RNA.Chrm5.merge, Stingele.RNA.Chrm5.6) #263 genes


##Plots: 
## density plot of  expression difference in stingele and depmap/CCLE data: 
ggplot(Stingele.RNA.Chrm5.merge, aes(x=meanDiff, y=RNA.Diff.Gain))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  theme_classic()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  coord_cartesian(ylim=c(-0.1,1), xlim=c(-0.1,1.75))+ 
  geom_smooth(method="lm", color="red")+
  ylab("Depmap: Chrm 5 RNA difference\nupon chrm gain")+
  xlab("Stingele: Chrm 5 RNA difference\n mean RPE1 & HCT116 Chrm 5 gain")
#4x4
# plot.Stingele.Depmap.Chrm5MeanDiff.RNA_density

#pearson correlation coefficient and p-value
cor.test(Stingele.RNA.Chrm5.merge$RNA.Diff.Gain, Stingele.RNA.Chrm5.merge$meanDiff, method="pearson")
# Pearson corr= 0.29
# p-value= 2E-05***
# n= 263
# cells = 2


##### Step 3: Now look at proteins on ANEUPLOID CHROMOSOME only, plot density#####
### Now look at proteins on ANEUPLOID CHROMOSOME only

### Chromosome 5 mean: 
## plot expression difference for chromosome 5:

# get genes on chrm 5, get mean expression from cells with chrm 5 gain
Stingele.Prot.Chrm5<-subset(Stingele.Prot, Chromosome.x==5)
Stingele.Prot.Chrm5$meanDiff<-rowMeans(Stingele.Prot.Chrm5[,c(4:7)], na.rm=TRUE)

Stingele.Prot.Chrm5.info = merge(
  Stingele.Prot.Chrm5,
  Stingele.info[!duplicated(Stingele.info[, c("Protein IDs")]), ],#remove duplicated gene info
  by=c("Protein IDs","Protein IDs"),
  all.x=TRUE
) #add gene name, be sure to not add duplicates


### Combine data with depmap gain loss difference data
# Stingele data has list of gene names: test depmap gene_symbol with name 1, if no match, test name 2. 

#step 1: get gene names? (first name is list, then second)
Stingele.Prot.Chrm5.info$Gene_Name1<-sub(";.*", "", Stingele.Prot.Chrm5.info$"Gene Names")
# 2nd name: get first of first name in list, then get rid of all genes after "new first" gene to get second gene in name list
Stingele.Prot.Chrm5.info$Gene_Name2<-sub(".*;", "", Stingele.Prot.Chrm5.info$"Gene Names")# get rid of first gene
Stingele.Prot.Chrm5.info$Gene_Name2<-sub(";.*", "", Stingele.Prot.Chrm5.info$Gene_Name2)
# 3rd and fourth names on list did not match with any missing data


#merge with first name in list
Stingele.DepmapProt.Chrm5.1<- merge(x=Stingele.Prot.Chrm5.info, #1388 genes
                                        y=CN.Diff.xRNA.yProt.ThreeGroups, 
                                        by.y="RNA_Name", 
                                        by.x="Gene_Name1") #945 genes
# now get all non-merged cells, merge with second name in list
nomerge<- anti_join(x=Stingele.Prot.Chrm5.info, #1388 genes
                    y=CN.Diff.xRNA.yProt.ThreeGroups, 
                    by=c("Gene_Name1"="RNA_Name")) #497 genes
Stingele.DepmapProt.Chrm5.2<- merge(x=nomerge, #1075 genes
                                        y=CN.Diff.xRNA.yProt.ThreeGroups, 
                                        by.y="RNA_Name", 
                                        by.x="Gene_Name2") #229 genes

Stingele.DepmapProt.Chrm5.3<-rbind(Stingele.DepmapProt.Chrm5.1, Stingele.DepmapProt.Chrm5.2)


#PLOTS: 
#scatterplot
## Plot all aneuploid gene difference: Stingele vs. depmap
ggplot(Stingele.DepmapProt.Chrm5.3, aes(x=Protein.Diff.Gain, y=meanDiff))+
  geom_point(color="black", size=2)+
  theme_classic()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  #coord_cartesian(ylim=c(-2,4), xlim=c(-2,4))+ 
  geom_smooth(method="lm", color="red")+ #linear trendline
  xlab("Depmap: Chrm 5 protein difference\nupon chrm gain")+
  ylab("Stingele: Chrm 5 protein difference\n mean RPE1 & HCT116 Chrm 5 gain")
#4x4
# plot.Stingele.Depmap.Chrm5MeanDiff

#Density plot
ggplot(Stingele.DepmapProt.Chrm5.3, aes(x=meanDiff, y=Protein.Diff.Gain))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  theme_classic()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  coord_cartesian(ylim=c(-0.25,0.75), xlim=c(-0.5,1.75))+ 
  geom_smooth(method="lm", color="red")+ #linear trendline
  ylab("Depmap: Chrm 5 protein difference\nupon chrm gain")+
  xlab("Stingele: Chrm 5 protein difference\n mean RPE1 & HCT116 Chrm 5 gain")
#4x4
# plot.Stingele.Depmap.Chrm5MeanDiff_density

#pearson correlation: 
cor.test(Stingele.DepmapProt.Chrm5.3$Protein.Diff.Gain, Stingele.DepmapProt.Chrm5.3$meanDiff, method="pearson")
# Pearson corr= 0.43
# p-value= 4E-10***
# n=372
