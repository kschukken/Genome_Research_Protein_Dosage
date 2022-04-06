##### Protein_RNA 1d density plots, per gene group #####
## author: Klaske M. Schukken

## Protein data analysis: Chromosome gain or loss vs diploid
## 3 category Grouping Analysis (scaling, buffered, anti-scaling)- Protein
## get scatterplot Protein difference per chrm arm 
## and 1d density plot per gene group (ribosome, transmembrane, etc) 


## 200914 start analysis
## Edited: 201217 -- Normalized data
## Edited: 210106 -- Normalize by Gene, not Cell line. mean=0 only, not SD. 
## Edited: 210121 -- no normalization beyond what was already present in the dataset. 
## Edited: 210203: use filtered cell lines- only cells with both RNA and Protein data
##                also measure changes in protein expression by gain/loss chromosome ARM (not whole chrm)
## Edited 210502: used filtered genes so I only use genes that have minimum of 10 cells per category
## Edited 210520: Use CN.DiffxRNA.yProtein data-- JUST DIFF UPON LOSS/GAIN of SPECIFC CHRM gene is located on!
##                was using old data with expression difference for all genes upon chrm loss 
##                (loose chrm 1, diff in gene express on chrm 3, etc.)
## Edited 210609: cleaned up code for publication. "_v2" 

library('ggplot2')
library('tidyverse')
library(xlsx)
library(readxl)
library(reshape2)


## Steps: 
# Step 1) Get data
#   get difference & Pvalue upon chrom gain/loss per gene (filtered cells, filtered genes)
#   for all chromosomes, not just diff upon gain/loss. also diff upon other chrm gain/loss
#   Get protein info: chromosome arm info
# Step 2) Seperate data into gene groups 
#    calculate if gene group is more or less buffered upon chrm gain/loss than other genes. 
# Step 3) plot desity (1d) plots of genes in specific gene group vs average. 
#   gene groups of interest: ribosome, transmembrane, metabolism, etc. 
#   y-axis scales change with each gene, be careful. 


##### Step 1) Get data and define functions #####

### density plots for gene groups relative to general gene expression
# Isolate Go term genes and view Protein / RNA difference in expression upon chrm gain/loss 
# ex: signalling moleculars, Ribosomal genes, G-protein coupled receptors, etc. 

# Step 1: get data 
#this data has 12k (all genes, not filter for 10+ cells/category) difference upon chrm gain/loss data
# so diff in express upon chrm gain/loss for genes on that chrm. 

## Filter for genes with minimum of 10+ cells per category to eliminate outliers. 
## Now filter Protein/RNA expression for genes that have a minimum of 10 cells per 
## category (gain, neutral, loss; in both RNA and protein): 
## 9414 genes in filtered data set 
## filter (10+ cells/category dataset from: Protein_RNA_Corr_min10.Filtered.R)
## CN.Diff.xRNA.yProt.ThreeGroups
# set working directory to location difference upon chromosome gain/loss dataset 
# (Supplementary data 2, or dataset from Protein_RNA_expression.PerCell_v2.R)

setwd()

CN.Diff.xRNA.yProt.ThreeGroups<- read.csv("RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points_3cat.csv")
CN.Diff.xRNA.yProt.ThreeGroups<- CN.Diff.xRNA.yProt.ThreeGroups[,-c(1)]

#note: all proteins are unique, but some RNA are repeated because 1 RNA can lead to multiple proteins
# we do not want to bias RNA analysis by analyzing one RNA multiple times, so I will subset data and isolate unique RNAs. 
CN.Diff.xRNA.unique<- CN.Diff.xRNA.yProt.ThreeGroups[,c(1:8,16,17)]
CN.Diff.xRNA.unique<- unique(CN.Diff.xRNA.unique) 

## get uniprot IDs, and add them to data for CORUM gene grouping section
Protein_Info3<-read.csv("Protein_location_info.csv")
Protein_ProID<-read.csv("Protein_ID_info.csv")
CN.Diff.xRNA.yProt.ThreeGroups2<- merge(x=CN.Diff.xRNA.yProt.ThreeGroups, y=Protein_ProID, by.x="Protein_ID", by.y="Protein_Id")


##Gene groups: 
CORUM.all<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=2) 
Ribosome<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=3) 
HSPA<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=4) 
Spliceosome<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=5) 
Autophagy<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=6) 
CellCycle<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=7) 
ER.Membrane<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=8) 
rRNAProcessing<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=9) 
Small_Mol_Metabolism<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=10) 
ExtracellMatrixStructural<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=11) 

##Tumor supressor gene (TSG) and oncogene (OG) list
# Bailey et al. Comprehensive Characterization of Cancer Driver Genes and Mutations, cell, 2018
TSG.OG<-read_xlsx("NIHMS948705_S8_Bailey2019.xlsx", sheet = 2) 



## FUNCTIONS
## define function to get significance in Protein/RNA for gain/loss relative to non-gene group expression: 
# Protein: 
P.Buff.or.Scale.Loss<- function(gene.list){
  testGenes<-CN.Diff.xRNA.yProt.ThreeGroups[CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% gene.list,]
  OtherGenes<-CN.Diff.xRNA.yProt.ThreeGroups[!CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% gene.list,]
  t.test(testGenes$Protein.Diff.Loss,OtherGenes$Protein.Diff.Loss)
}
P.Buff.or.Scale.Gain<- function(gene.list){
  testGenes<-CN.Diff.xRNA.yProt.ThreeGroups[CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% gene.list,]
  OtherGenes<-CN.Diff.xRNA.yProt.ThreeGroups[!CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% gene.list,]
  t.test(testGenes$Protein.Diff.Gain,OtherGenes$Protein.Diff.Gain)
}

#RNA
R.Buff.or.Scale.Loss<- function(gene.list){
  testGenes<-CN.Diff.xRNA.unique[CN.Diff.xRNA.unique$RNA_Name %in% gene.list,]
  OtherGenes<-CN.Diff.xRNA.unique[!CN.Diff.xRNA.unique$RNA_Name %in% gene.list,]
  t.test(testGenes$RNA.Diff.Loss,OtherGenes$RNA.Diff.Loss)
}
R.Buff.or.Scale.Gain<- function(gene.list){
  testGenes<-CN.Diff.xRNA.unique[CN.Diff.xRNA.unique$RNA_Name %in% gene.list,]
  OtherGenes<-CN.Diff.xRNA.unique[!CN.Diff.xRNA.unique$RNA_Name %in% gene.list,]
  t.test(testGenes$RNA.Diff.Gain,OtherGenes$RNA.Diff.Gain)
}



##### Step 2) Subset data by group and test (Protein and RNA) ####

## Some gene terms are based on gene sybol (like ribosomal subunits), while others are based on go terms. 
# download the corresponding go term from http://www.informatics.jax.org, save in a folder and upload into R as needed. 
# not all of these terms were used in the final paper. 

# additionally, get uniprot_corum_mapping.txt, from CORUM website http://mips.helmholtz-muenchen.de/corum/#download

## Get all Ribosomal proteins: 
Protein.Ribosome  <- subset(CN.Diff.xRNA.yProt.ThreeGroups, grepl( "RPS", CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name, fixed = TRUE))
Protein.Ribosome  <- Protein.Ribosome[order(Protein.Ribosome$RNA_Name),]
Protein.Ribosome  <- Protein.Ribosome[32:69,] #remove mitochindrial ribosome genes, PRPSAP genes, and TRPS genes
Protein.Ribosome2 <- subset(CN.Diff.xRNA.yProt.ThreeGroups, grepl( "RPL", CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name, fixed = TRUE))
Protein.Ribosome2 <- Protein.Ribosome2[order(Protein.Ribosome2$RNA_Name),]
Protein.Ribosome2 <- Protein.Ribosome2[47:92,] #remove mitochindrial ribosome genes

Protein.Gain.Ribosome <- rbind(Protein.Ribosome, Protein.Ribosome2) #this gives 84 proteins
P.Ribosome.list <- Protein.Gain.Ribosome$RNA_Name

x<- P.Buff.or.Scale.Gain(P.Ribosome.list) # <2E-16 anti scale
P.Buff.or.Scale.Loss(P.Ribosome.list) # 0.016 buffer
R.Buff.or.Scale.Gain(P.Ribosome.list) # NS
x<-R.Buff.or.Scale.Loss(P.Ribosome.list) # 4E-08 scale


## Get HSPA proteins
## buffered upon chrm gain/loss
Protein.Gain.HSP <- subset(CN.Diff.xRNA.yProt.ThreeGroups, grepl( "HSPA", CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name, fixed = TRUE))
Protein.Loss.HSP <- subset(CN.Diff.xRNA.yProt.ThreeGroups, grepl( "HSPA", CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name, fixed = TRUE))
HSPA.list<-Protein.Loss.HSP$RNA_Name

P.Buff.or.Scale.Gain(HSPA.list) # NS
P.Buff.or.Scale.Loss(HSPA.list) # NS
R.Buff.or.Scale.Gain(HSPA.list) # NS
R.Buff.or.Scale.Loss(HSPA.list) # NS


###Cell Cycle, GO:0007049
## Scale with chrm gain loss
CellCycle<- read.delim2("GO_term_CellCycle.txt",
                        dec=".", header = TRUE, sep="\t", row.names=NULL) 

CellCycle.list<-unique(CellCycle$Symbol)#get gene names
CellCycle.list<-toupper(CellCycle.list)#make string and uppercase

P.Buff.or.Scale.Gain(CellCycle.list) # NS
P.Buff.or.Scale.Loss(CellCycle.list) # NS
x<- R.Buff.or.Scale.Gain(CellCycle.list) # 6E-04
x<- R.Buff.or.Scale.Loss(CellCycle.list) # 3E-09



### rRNA processing GO: 0006364
## buffetted upon chrm gain/loss
#rRNA<- read.delim2("GO_term_rRNAProcessing.txt",
#                             dec=".", header = TRUE, sep="\t", row.names=NULL) 
rRNA.list<-unique(rRNAProcessing$Symbol)#get gene names
rRNA.list<-toupper(rRNA.list)#make string and uppercase

P.Buff.or.Scale.Gain(rRNA.list) # 3E-06 buffering
P.Buff.or.Scale.Loss(rRNA.list) # 0.018 buffering
R.Buff.or.Scale.Gain(rRNA.list) # 0.017
R.Buff.or.Scale.Loss(rRNA.list) # 4E-14



#### Extracellular matrix structural constituant GO:0005201
#EMSC<- read.delim2("GO_term_EMSC.txt",
#                          dec=".", header = TRUE, sep="\t", row.names=NULL) 
EMSC.list<-unique(ExtracellMatrixStructural$Symbol)#get gene names
EMSC.list<-toupper(EMSC.list)#make string and uppercase

P.Buff.or.Scale.Gain(gene.list=EMSC.list) # 0.045 buffer
P.Buff.or.Scale.Loss(gene.list=EMSC.list) # 3E-09 anti-scale
R.Buff.or.Scale.Gain(gene.list=EMSC.list) # NS
R.Buff.or.Scale.Loss(gene.list=EMSC.list) # 6E-11

#ExtraMatrix.list<-CN.Diff.xRNA.yProt.ThreeGroups[CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% EMSC.list,]



### Small molecule metabolic process GO:0044281 ##
#SM.Metabolism<- read.delim2("GO_term_SmallMoleculemetabolicPathway.txt",
#                         dec=".", header = TRUE, sep="\t", row.names=NULL) 
SM.Metabolism.list<-unique(Small_Mol_Metabolism$Symbol)#get gene names
SM.Metabolism.list<-toupper(SM.Metabolism.list)#make string and uppercase

P.Buff.or.Scale.Gain(gene.list=SM.Metabolism.list) # 1E-04 scale
x<-P.Buff.or.Scale.Loss(gene.list=SM.Metabolism.list) # 5E-14 scale
R.Buff.or.Scale.Gain(gene.list=SM.Metabolism.list) # NS
R.Buff.or.Scale.Loss(gene.list=SM.Metabolism.list) # 7E-03




### mRNA splicing, via spliceosome  GO:0000398
#spliceosome<- read.delim2("GO_term_Spliceosome.txt",
#                         dec=".", header = TRUE, sep="\t", row.names=NULL) 
spliceosome.list<-unique(Spliceosome$Symbol)#get gene names
spliceosome.list<-toupper(spliceosome.list)#make string and uppercase

x<- P.Buff.or.Scale.Gain(gene.list=spliceosome.list) # 8E-05 buffer
x<- P.Buff.or.Scale.Loss(gene.list=spliceosome.list) # 4E-05 buffer
x<- R.Buff.or.Scale.Gain(gene.list=spliceosome.list) # NS
x<- R.Buff.or.Scale.Loss(gene.list=spliceosome.list) # 8E-06 scale



# autophagy GO:0006914
#Authophagy<- read.delim2("GO_term_Authophagy.txt",
#                         dec=".", header = TRUE, sep="\t", row.names=NULL) 
Authophagy.list<-unique(Autophagy$Symbol)#get gene names
Authophagy.list<-toupper(Authophagy.list)#make string and uppercase

P.Buff.or.Scale.Gain(gene.list=Authophagy.list) # 0.08 NS scale
P.Buff.or.Scale.Loss(gene.list=Authophagy.list) # 0.001131 scale
R.Buff.or.Scale.Gain(gene.list=Authophagy.list) # 0.038 scale
R.Buff.or.Scale.Loss(gene.list=Authophagy.list) # 0.085 NS scale



#### CORUM proteins (protein complexes)
setwd()

#CORUM.Uniprot <- read.delim2("uniprot_corum_mapping.txt") #all CORUM IDs 

CORUM.all2<- subset(CN.Diff.xRNA.yProt.ThreeGroups2, Uniprot_Acc %in% CORUM.all$Uniprot_ID)
CORUM.all3<-CORUM.all2$RNA_Name

x<- P.Buff.or.Scale.Gain(gene.list=CORUM.all3) # <2E-16 buffer
x<- P.Buff.or.Scale.Loss(gene.list=CORUM.all3) #  2E-07 buffer
x<- R.Buff.or.Scale.Gain(gene.list=CORUM.all3) #  0.0011 scale
x<- R.Buff.or.Scale.Loss(gene.list=CORUM.all3) #  8E-10 scale


CORUM.gain2 <- CN.Diff.xRNA.yProt.ThreeGroups[CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% CORUM.all3,]
# gain =  0.1284378 --> 9.31094
# loss = -0.1101569 --> -7.35127%
# all gain = 0.1603291 --> 11.7542% 
# all loss = -0.1268382 --> -8.416361%

###Of the 2507 CORUM genes in the dataset: 
## upon gain: 7.9% anti-scale(199), 69% are buffered (1723), and only 23% scale(585)
    ## mean proteins upon gain: 60% buffered (5572 total)
## upon loss: 7% anti-scale (177), 76% are buffered (1915), and only 17% scale(415)
    ## mean protein upon loss: 66% buffered (6201 total)


## See what percent of CORUM genes are in other terms/categories
## if we remove some of these enriched terms, are CORUM genes still buffered? 
## Find what % of all genes fall into the other components. 
length(P.Ribosome.list)*100/length(unique(CN.Diff.xRNA.yProt.ThreeGroups$Protein_ID)) # 0.89%
length(HSPA.list)*100/length(unique(CN.Diff.xRNA.yProt.ThreeGroups$Protein_ID)) # 0.14%
length(CellCycle.list)*100/length(unique(CN.Diff.xRNA.yProt.ThreeGroups$Protein_ID)) # 18.1%
length(Authophagy.list)*100/length(unique(CN.Diff.xRNA.yProt.ThreeGroups$Protein_ID)) # 4.6%
length(SM.Metabolism.list)*100/length(unique(CN.Diff.xRNA.yProt.ThreeGroups$Protein_ID)) # 17.8%
length(rRNA.list)*100/length(unique(CN.Diff.xRNA.yProt.ThreeGroups$Protein_ID)) # 2.2%
length(spliceosome.list)*100/length(unique(CN.Diff.xRNA.yProt.ThreeGroups$Protein_ID)) # 3.0%
length(EMSC.list)*100/length(unique(CN.Diff.xRNA.yProt.ThreeGroups$Protein_ID)) # 1.6%

## Find what % of CORUM genes fall into the other components. 
sum(unique(P.Ribosome.list) %in% CORUM.all3)*100/length(CORUM.all3) # 2.7%
sum(unique(HSPA.list) %in% CORUM.all3)*100/length(CORUM.all3) # 0.20%
sum(unique(CellCycle.list) %in% CORUM.all3)*100/length(CORUM.all3) # 19.32%
sum(unique(Authophagy.list) %in% CORUM.all3)*100/length(CORUM.all3) # 4.4%
sum(unique(SM.Metabolism.list) %in% CORUM.all3)*100/length(CORUM.all3) # 5.9%
sum(unique(rRNA.list) %in% CORUM.all3)*100/length(CORUM.all3) # 2.9%
sum(unique(spliceosome.list) %in% CORUM.all3)*100/length(CORUM.all3) # 5.4%
sum(unique(EMSC.list) %in% CORUM.all3)*100/length(CORUM.all3) # 0.69%

## Find ratio between percent term in CORUM genes, over percent term in total genes
ratio1<- length(unique(CN.Diff.xRNA.yProt.ThreeGroups$Protein_ID))/length(CORUM.all3)
sum(unique(P.Ribosome.list) %in% CORUM.all3)*ratio1/length(P.Ribosome.list) #3.1
sum(unique(HSPA.list) %in% CORUM.all3)*ratio1/length(HSPA.list) #1.5
sum(unique(CellCycle.list) %in% CORUM.all3)*ratio1/length(CellCycle.list) #1.07
sum(unique(Authophagy.list) %in% CORUM.all3)*ratio1/length(Authophagy.list) #0.94
sum(unique(SM.Metabolism.list) %in% CORUM.all3)*ratio1/length(SM.Metabolism.list) #0.33
sum(unique(rRNA.list) %in% CORUM.all3)*ratio1/length(rRNA.list) #1.3
sum(unique(spliceosome.list) %in% CORUM.all3)*ratio1/length(spliceosome.list) # 1.8
sum(unique(EMSC.list) %in% CORUM.all3)*ratio1/length(EMSC.list) #0.44

# Combine all buffered terms enriched in CORUM genes
enriched_terms<-append(P.Ribosome.list, HSPA.list) 
enriched_terms<-append(enriched_terms, spliceosome.list) 
enriched_terms<-append(enriched_terms, rRNA.list) 
enriched_terms<-append(enriched_terms, CellCycle.list) 

# get all CORUM genes not in the above buffered terms
Subset_CORUM<- CORUM.all3[CORUM.all3 %in% enriched_terms == FALSE]

# CORUM genes are still buffered, even when you remove the previously identified buffered terms
P.Buff.or.Scale.Gain(gene.list=Subset_CORUM) #  7E-13 buffer
P.Buff.or.Scale.Loss(gene.list=Subset_CORUM) #  2E-04 buffer


#####         TSG and Oncogene analysis ####
##Tumor supressor gene (TSG) and oncogene (OG) list
# Bailey et al. Comprehensive Characterization of Cancer Driver Genes and Mutations, cell, 2018


Oncogene.Bailey<- subset(TSG.OG, TSG.OG$`Tumor suppressor or oncogene prediction (by 20/20+)`=="oncogene")
Oncogene.Bailey.list<-unique(Oncogene.Bailey$Gene) #84

TSG.Bailey<- subset(TSG.OG, TSG.OG$`Tumor suppressor or oncogene prediction (by 20/20+)`=="tsg")
TSG.Bailey.list<-unique(TSG.Bailey$Gene) #99

P.Buff.or.Scale.Gain(Oncogene.Bailey.list) # NS --> x8 correction= NS
P.Buff.or.Scale.Loss(Oncogene.Bailey.list) # NS--> x8 correction= NS
R.Buff.or.Scale.Gain(Oncogene.Bailey.list) # 0.0014--> x8 correction= 0.01
R.Buff.or.Scale.Loss(Oncogene.Bailey.list) # NS--> x8 correction= NS

P.Buff.or.Scale.Gain(TSG.Bailey.list) # NS--> x8 correction= NS
P.Buff.or.Scale.Loss(TSG.Bailey.list) # NS--> x8 correction= NS
R.Buff.or.Scale.Gain(TSG.Bailey.list) # NS--> x8 correction= NS
R.Buff.or.Scale.Loss(TSG.Bailey.list) # 0.033--> x8 correction= NS

#TSG OG t-test: 

P.Buff.or.Scale.Gain(gene.list=TSG.Bailey.list) # P = NS
P.Buff.or.Scale.Loss(gene.list=TSG.Bailey.list) # P = NS

P.Buff.or.Scale.Gain(gene.list=Oncogene.Bailey.list) # P = NS
P.Buff.or.Scale.Loss(gene.list=Oncogene.Bailey.list) # P = NS


## Ran all buffeted & scaling gene through gprofiler
## Nothing significant. see ***.R

#TSG.Bailey.Loss
#TSG.Bailey.Gain
#Oncogene.Bailey.Loss
#Oncogene.Bailey.Gain



##### Step 3) Plots: 1d density, barplot, scatter: all Protein diff vs gene group ####
# Plot 1d density of gene group vs all genes not in group
# for RNA & protein, upon chrm gain and loss


#SM.Metabolism.list # scale
#Authophagy.list #scale
#P.Ribosome.list # buffer
#rRNA.list # buffer
#CORUM.all3 # buffer
#spliceosome.list #buffer
#CellCycle.list # buffer, NS from mean gene
#PlasmaMembrane.list # buffer upon loss
#ExtraMatrix.list #anti scale  upon loss
#IntMembrane.list #scaling upon gain, buffered upon loss
#HSPA.list # NS
#ExtraMatrix.list #
#EMSC.list #
#transmembrane.list #
#Metabolism.list #
#Sphingolipid.list #
#Oncogene.Bailey.list
#TSG.Bailey.list

# Set working directory to where you want the graphs to end up: 
#setwd()

##INSTRUCTIONS: add list of genes you want to examine, and add group name and ID accordingly
## groupID will be added to .pdf title of the file. 

Test.List= P.Ribosome.list
GroupName= "Ribosomal genes"
GroupID= "Ribosome_v3"


## density plot PROTEIN GAIN 
pdf(file = paste0("plot.Protein.gain.density.", GroupID, "_v2.pdf"),
    width = 4, 
    height = 4)
  # non-list genes
ggplot(CN.Diff.xRNA.yProt.ThreeGroups, 
       aes(x=Protein.Diff.Gain))+
  geom_density(data=CN.Diff.xRNA.yProt.ThreeGroups[!CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% Test.List,],
               color="black", size=1)+
  geom_vline(xintercept=mean(CN.Diff.xRNA.yProt.ThreeGroups[!CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% Test.List,]$Protein.Diff.Gain), 
             color="black", linetype="dotted", size=1)+
  #list genes
  geom_density(data=CN.Diff.xRNA.yProt.ThreeGroups[CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% Test.List,], 
               color="red", size=1)+
  geom_vline(xintercept=mean(CN.Diff.xRNA.yProt.ThreeGroups[CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% Test.List,]$Protein.Diff.Gain), 
             color="red", linetype="dotted", size=1)+
  
  ylab("Density") +
  xlab("Difference in protein expression upon chrm gain")+
  ggtitle (GroupName)+
  coord_cartesian(xlim=c(-1,1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()
dev.off()

## density plot PROTEIN LOSS
pdf(file = paste0("plot.Protein.loss.density.", GroupID, "_v2.pdf"),
    width = 4, 
    height = 4)
ggplot(CN.Diff.xRNA.yProt.ThreeGroups, 
       aes(x=Protein.Diff.Loss))+
  geom_density(data=CN.Diff.xRNA.yProt.ThreeGroups[!CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% Test.List,], 
               color="black", size=1)+
  geom_vline(xintercept=mean(CN.Diff.xRNA.yProt.ThreeGroups[!CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% Test.List,]$Protein.Diff.Loss), 
             color="black", linetype="dotted", size=1)+
  
  geom_density(data=CN.Diff.xRNA.yProt.ThreeGroups[CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% Test.List,], 
               color="dodgerblue3", size=1)+
  geom_vline(xintercept=mean(CN.Diff.xRNA.yProt.ThreeGroups[CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% Test.List,]$Protein.Diff.Loss), 
             color="dodgerblue3", linetype="dotted", size=1)+
  
  ylab("Density") + 
  xlab("Difference in protein expression upon chrm loss")+
  ggtitle (GroupName)+
  coord_cartesian(xlim=c(-1,1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()
dev.off()


### RNA
pdf(file = paste0("plot.RNA.gain.density.", GroupID, "_v2.pdf"),
    width = 4, 
    height = 4)
## density plot RNA GAIN 
ggplot(CN.Diff.xRNA.unique, 
       aes(x=RNA.Diff.Gain))+
  geom_density(data=CN.Diff.xRNA.unique[!CN.Diff.xRNA.unique$RNA_Name %in% Test.List,], 
               color="black", size=1)+
  geom_vline(xintercept=mean(CN.Diff.xRNA.unique[!CN.Diff.xRNA.unique$RNA_Name %in% Test.List,]$RNA.Diff.Gain), 
             color="black", linetype="dotted", size=1)+
  
  geom_density(data=CN.Diff.xRNA.unique[CN.Diff.xRNA.unique$RNA_Name %in% Test.List,], 
               color="red", size=1)+
  geom_vline(xintercept=mean(CN.Diff.xRNA.unique[CN.Diff.xRNA.unique$RNA_Name %in% Test.List,]$RNA.Diff.Gain), 
             color="red", linetype="dotted", size=1)+
  ylab("Density") +
  xlab("Difference in RNA expression upon chrm gain")+
  ggtitle (GroupName)+
  coord_cartesian(xlim=c(-1,1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()
dev.off()

## density plot RNA LOSS
pdf(file = paste0("plot.RNA.loss.density.", GroupID, "_v2.pdf"),
    width = 4, 
    height = 4)
ggplot(CN.Diff.xRNA.unique, 
       aes(x=RNA.Diff.Loss))+
  geom_density(data=CN.Diff.xRNA.unique[!CN.Diff.xRNA.unique$RNA_Name %in% Test.List,], 
               color="black", size=1)+
  geom_vline(xintercept=mean(CN.Diff.xRNA.unique[!CN.Diff.xRNA.unique$RNA_Name %in% Test.List,]$RNA.Diff.Loss), 
             color="black", linetype="dotted", size=1)+
  
  geom_density(data=CN.Diff.xRNA.unique[CN.Diff.xRNA.unique$RNA_Name %in% Test.List,], 
               color="dodgerblue3", size=1)+
  geom_vline(xintercept=mean(CN.Diff.xRNA.unique[CN.Diff.xRNA.unique$RNA_Name %in% Test.List,]$RNA.Diff.Loss), 
             color="dodgerblue3", linetype="dotted", size=1)+
  
  ylab("Density") +
  xlab("Difference in RNA expression upon chrm loss")+
  ggtitle (GroupName)+
  coord_cartesian(xlim=c(-1,1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()
dev.off()

print("End of graphs, Protein and RNA")

##### Step 4) look at low  aneuploidy only cells: ribosomes still buffered? #####

CN.Diff.RNA.Prot_LowPloidy # from Protein_RNA_expression.PerCell_v2.R

CN.Diff.xRNA.LowPloidy.unique<- CN.Diff.RNA.Prot_LowPloidy[,c(1:8,16,17)]
CN.Diff.xRNA.LowPloidy.unique<- unique(CN.Diff.xRNA.LowPloidy.unique) ## Unique RNA genes

## define function to look at difference upon chrm gain/loss in low aneuploid cells
P.Buff.or.Scale.Loss.lowAneu<- function(gene.list){
  testGenes<-CN.Diff.RNA.Prot_LowPloidy[CN.Diff.RNA.Prot_LowPloidy$RNA_Name %in% gene.list,]
  OtherGenes<-CN.Diff.RNA.Prot_LowPloidy[!CN.Diff.RNA.Prot_LowPloidy$RNA_Name %in% gene.list,]
  t.test(testGenes$Protein.Diff.Loss,OtherGenes$Protein.Diff.Loss)
}
P.Buff.or.Scale.Gain.lowAneu<- function(gene.list){
  testGenes<-CN.Diff.RNA.Prot_LowPloidy[CN.Diff.RNA.Prot_LowPloidy$RNA_Name %in% gene.list,]
  OtherGenes<-CN.Diff.RNA.Prot_LowPloidy[!CN.Diff.RNA.Prot_LowPloidy$RNA_Name %in% gene.list,]
  t.test(testGenes$Protein.Diff.Gain,OtherGenes$Protein.Diff.Gain)
}
R.Buff.or.Scale.Loss.lowAneu<- function(gene.list){
  testGenes<-CN.Diff.xRNA.LowPloidy.unique[CN.Diff.xRNA.LowPloidy.unique$RNA_Name %in% gene.list,]
  OtherGenes<-CN.Diff.xRNA.LowPloidy.unique[!CN.Diff.xRNA.LowPloidy.unique$RNA_Name %in% gene.list,]
  t.test(testGenes$RNA.Diff.Loss,OtherGenes$RNA.Diff.Loss)
}
R.Buff.or.Scale.Gain.lowAneu<- function(gene.list){
  testGenes<-CN.Diff.xRNA.LowPloidy.unique[CN.Diff.xRNA.LowPloidy.unique$RNA_Name %in% gene.list,]
  OtherGenes<-CN.Diff.xRNA.LowPloidy.unique[!CN.Diff.xRNA.LowPloidy.unique$RNA_Name %in% gene.list,]
  t.test(testGenes$RNA.Diff.Gain,OtherGenes$RNA.Diff.Gain)
}



## Get all Ribosomal proteins: 
Protein.Ribosome <- subset(CN.Diff.RNA.Prot_LowPloidy, grepl( "RPS", CN.Diff.RNA.Prot_LowPloidy$RNA_Name, fixed = TRUE))
Protein.Ribosome <- Protein.Ribosome[order(as.character(Protein.Ribosome$RNA_Name)),]
Protein.Ribosome <- Protein.Ribosome[10:31,] #remove mitochindrial ribosome genes, PRPSAP genes, and TRPS genes
Protein.Ribosome2 <- subset(CN.Diff.RNA.Prot_LowPloidy, grepl( "RPL", CN.Diff.RNA.Prot_LowPloidy$RNA_Name, fixed = TRUE))
Protein.Ribosome2 <- Protein.Ribosome2[order(as.character(Protein.Ribosome2$RNA_Name)),]
Protein.Ribosome2 <- Protein.Ribosome2[20:33,] #remove mitochindrial ribosome genes

Protein.Gain.Ribosome <- rbind(Protein.Ribosome, Protein.Ribosome2) #this gives 84 proteins


P.Buff.or.Scale.Gain.lowAneu(Protein.Gain.Ribosome$RNA_Name) # buffer 0.024
P.Buff.or.Scale.Loss.lowAneu(Protein.Gain.Ribosome$RNA_Name) # anti-scale 5E-05
R.Buff.or.Scale.Gain.lowAneu(Protein.Gain.Ribosome$RNA_Name) # scale 0.025
R.Buff.or.Scale.Loss.lowAneu(Protein.Gain.Ribosome$RNA_Name) # NS

# rRNA processing

#rRNA<- read.delim2("GO_term_rRNAProcessing.txt",
#                   dec=".", header = TRUE, sep="\t", row.names=NULL) 
rRNA.list<-unique(rRNAProcessing$Symbol)#get gene names
rRNA.list<-toupper(rRNA.list)#make string and uppercase

P.Buff.or.Scale.Gain.lowAneu(rRNA.list) # NS
P.Buff.or.Scale.Loss.lowAneu(rRNA.list) # NS
R.Buff.or.Scale.Gain.lowAneu(rRNA.list) # NS
R.Buff.or.Scale.Loss.lowAneu(rRNA.list) # 0.0019 scaling

##INSTRUCTIONS: add list of genes you want to examine, and add group name and ID accordingly
## groupID will be added to .pdf title of the file. 
#setwd()

Test.List= rRNA.list
GroupName= "rRNA processing genes (Low aneuploid cells)"
GroupID= "rRNA_LowPloidy"


## density plot PROTEIN GAIN 
pdf(file = paste0("plot.Protein.gain.density.", GroupID, "_LowCellPloidy.pdf"),
    width = 4, 
    height = 4)
# non-list genes
ggplot(CN.Diff.RNA.Prot_LowPloidy, 
       aes(x=Protein.Diff.Gain))+
  geom_density(data=CN.Diff.RNA.Prot_LowPloidy[!CN.Diff.RNA.Prot_LowPloidy$RNA_Name %in% Test.List,],
               color="black", size=1)+
  geom_vline(xintercept=mean(CN.Diff.RNA.Prot_LowPloidy[!CN.Diff.RNA.Prot_LowPloidy$RNA_Name %in% Test.List,]$Protein.Diff.Gain), 
             color="black", linetype="dotted", size=1)+
  #list genes
  geom_density(data=CN.Diff.RNA.Prot_LowPloidy[CN.Diff.RNA.Prot_LowPloidy$RNA_Name %in% Test.List,], 
               color="red", size=1)+
  geom_vline(xintercept=mean(CN.Diff.RNA.Prot_LowPloidy[CN.Diff.RNA.Prot_LowPloidy$RNA_Name %in% Test.List,]$Protein.Diff.Gain), 
             color="red", linetype="dotted", size=1)+
  
  ylab("Density") +
  xlab("Difference in protein expression upon chrm gain (Low ploidy cells)")+
  ggtitle (GroupName)+
  coord_cartesian(xlim=c(-1,1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()
dev.off()

## density plot PROTEIN LOSS
pdf(file = paste0("plot.Protein.loss.density.", GroupID, "_LowCellPloidy.pdf"),
    width = 4, 
    height = 4)
ggplot(CN.Diff.xRNA.yProt.ThreeGroups, 
       aes(x=Protein.Diff.Loss))+
  geom_density(data=CN.Diff.RNA.Prot_LowPloidy[!CN.Diff.RNA.Prot_LowPloidy$RNA_Name %in% Test.List,], 
               color="black", size=1)+
  geom_vline(xintercept=mean(CN.Diff.RNA.Prot_LowPloidy[!CN.Diff.RNA.Prot_LowPloidy$RNA_Name %in% Test.List,]$Protein.Diff.Loss), 
             color="black", linetype="dotted", size=1)+
  
  geom_density(data=CN.Diff.RNA.Prot_LowPloidy[CN.Diff.RNA.Prot_LowPloidy$RNA_Name %in% Test.List,], 
               color="dodgerblue3", size=1)+
  geom_vline(xintercept=mean(CN.Diff.RNA.Prot_LowPloidy[CN.Diff.RNA.Prot_LowPloidy$RNA_Name %in% Test.List,]$Protein.Diff.Loss), 
             color="dodgerblue3", linetype="dotted", size=1)+
  
  ylab("Density") + 
  xlab("Difference in protein expression upon chrm loss (Low ploidy cells)")+
  ggtitle (GroupName)+
  coord_cartesian(xlim=c(-1,1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()
dev.off()

### RNA
pdf(file = paste0("plot.RNA.gain.density.", GroupID, "_LowPloidy.pdf"),
    width = 4, 
    height = 4)
## density plot RNA GAIN 
ggplot(CN.Diff.xRNA.LowPloidy.unique, 
       aes(x=RNA.Diff.Gain))+
  geom_density(data=CN.Diff.xRNA.LowPloidy.unique[!CN.Diff.xRNA.LowPloidy.unique$RNA_Name %in% Test.List,], 
               color="black", size=1)+
  geom_vline(xintercept=mean(CN.Diff.xRNA.LowPloidy.unique[!CN.Diff.xRNA.LowPloidy.unique$RNA_Name %in% Test.List,]$RNA.Diff.Gain), 
             color="black", linetype="dotted", size=1)+
  
  geom_density(data=CN.Diff.xRNA.LowPloidy.unique[CN.Diff.xRNA.LowPloidy.unique$RNA_Name %in% Test.List,], 
               color="red", size=1)+
  geom_vline(xintercept=mean(CN.Diff.xRNA.LowPloidy.unique[CN.Diff.xRNA.LowPloidy.unique$RNA_Name %in% Test.List,]$RNA.Diff.Gain), 
             color="red", linetype="dotted", size=1)+
  ylab("Density") +
  xlab("Difference in RNA expression upon chrm gain (Low ploidy cells)")+
  ggtitle (GroupName)+
  coord_cartesian(xlim=c(-1,1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()
dev.off()

## density plot RNA LOSS
pdf(file = paste0("plot.RNA.loss.density.", GroupID, "_LowPloidy.pdf"),
    width = 4, 
    height = 4)
ggplot(CN.Diff.xRNA.LowPloidy.unique, 
       aes(x=RNA.Diff.Loss))+
  geom_density(data=CN.Diff.xRNA.LowPloidy.unique[!CN.Diff.xRNA.LowPloidy.unique$RNA_Name %in% Test.List,], 
               color="black", size=1)+
  geom_vline(xintercept=mean(CN.Diff.xRNA.LowPloidy.unique[!CN.Diff.xRNA.LowPloidy.unique$RNA_Name %in% Test.List,]$RNA.Diff.Loss), 
             color="black", linetype="dotted", size=1)+
  
  geom_density(data=CN.Diff.xRNA.LowPloidy.unique[CN.Diff.xRNA.LowPloidy.unique$RNA_Name %in% Test.List,], 
               color="dodgerblue3", size=1)+
  geom_vline(xintercept=mean(CN.Diff.xRNA.LowPloidy.unique[CN.Diff.xRNA.LowPloidy.unique$RNA_Name %in% Test.List,]$RNA.Diff.Loss), 
             color="dodgerblue3", linetype="dotted", size=1)+
  
  ylab("Density") +
  xlab("Difference in RNA expression upon chrm loss (Low ploidy cells)")+
  ggtitle (GroupName)+
  coord_cartesian(xlim=c(-1,1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()
dev.off()

print("End of graphs, Protein and RNA")

##### Step 5) look at high aneuploidy only cells: ribosomes still buffered? #####

CN.Diff.RNA.Prot_HighPloidy # from Protein_RNA_expression.PerCell_v2.R

CN.Diff.xRNA.HighPloidy.unique<- CN.Diff.RNA.Prot_HighPloidy[,c(1:8,16,17)]
CN.Diff.xRNA.HighPloidy.unique<- unique(CN.Diff.xRNA.HighPloidy.unique) 

## define function to look at difference upon chrm gain/loss in low aneuploid cells
P.Buff.or.Scale.Loss.highAneu<- function(gene.list){
  testGenes<-CN.Diff.RNA.Prot_HighPloidy[CN.Diff.RNA.Prot_HighPloidy$RNA_Name %in% gene.list,]
  OtherGenes<-CN.Diff.RNA.Prot_HighPloidy[!CN.Diff.RNA.Prot_HighPloidy$RNA_Name %in% gene.list,]
  t.test(testGenes$Protein.Diff.Loss,OtherGenes$Protein.Diff.Loss)
}
P.Buff.or.Scale.Gain.highAneu<- function(gene.list){
  testGenes<-CN.Diff.RNA.Prot_HighPloidy[CN.Diff.RNA.Prot_HighPloidy$RNA_Name %in% gene.list,]
  OtherGenes<-CN.Diff.RNA.Prot_HighPloidy[!CN.Diff.RNA.Prot_HighPloidy$RNA_Name %in% gene.list,]
  t.test(testGenes$Protein.Diff.Gain,OtherGenes$Protein.Diff.Gain)
}
R.Buff.or.Scale.Loss.highAneu<- function(gene.list){
  testGenes<-CN.Diff.xRNA.HighPloidy.unique[CN.Diff.xRNA.HighPloidy.unique$RNA_Name %in% gene.list,]
  OtherGenes<-CN.Diff.xRNA.HighPloidy.unique[!CN.Diff.xRNA.HighPloidy.unique$RNA_Name %in% gene.list,]
  t.test(testGenes$RNA.Diff.Loss,OtherGenes$RNA.Diff.Loss)
}
R.Buff.or.Scale.Gain.highAneu<- function(gene.list){
  testGenes<-CN.Diff.xRNA.HighPloidy.unique[CN.Diff.xRNA.HighPloidy.unique$RNA_Name %in% gene.list,]
  OtherGenes<-CN.Diff.xRNA.HighPloidy.unique[!CN.Diff.xRNA.HighPloidy.unique$RNA_Name %in% gene.list,]
  t.test(testGenes$RNA.Diff.Gain,OtherGenes$RNA.Diff.Gain)
}



## Get all Ribosomal proteins: 
Protein.Ribosome <- subset(CN.Diff.RNA.Prot_HighPloidy, grepl( "RPS", CN.Diff.RNA.Prot_HighPloidy$RNA_Name, fixed = TRUE))
Protein.Ribosome <- Protein.Ribosome[order(as.character(Protein.Ribosome$RNA_Name)),]
Protein.Ribosome <- Protein.Ribosome[31:65,] #remove mitochindrial ribosome genes, PRPSAP genes, and TRPS genes
Protein.Ribosome2 <- subset(CN.Diff.RNA.Prot_HighPloidy, grepl( "RPL", CN.Diff.RNA.Prot_HighPloidy$RNA_Name, fixed = TRUE))
Protein.Ribosome2 <- Protein.Ribosome2[order(as.character(Protein.Ribosome2$RNA_Name)),]
Protein.Ribosome2 <- Protein.Ribosome2[45:89,] #remove mitochindrial ribosome genes

Protein.Gain.Ribosome <- rbind(Protein.Ribosome, Protein.Ribosome2) #this gives 84 proteins


P.Buff.or.Scale.Gain.highAneu(Protein.Gain.Ribosome$RNA_Name) # buffer 3E-06
P.Buff.or.Scale.Loss.highAneu(Protein.Gain.Ribosome$RNA_Name) # buffer 5E-07
R.Buff.or.Scale.Gain.highAneu(Protein.Gain.Ribosome$RNA_Name) # NS
R.Buff.or.Scale.Loss.highAneu(Protein.Gain.Ribosome$RNA_Name) # NS

#rRNA<- read.delim2("GO_term_rRNAProcessing.txt",
#                   dec=".", header = TRUE, sep="\t", row.names=NULL) 
rRNA.list<-unique(rRNAProcessing$Symbol)#get gene names
rRNA.list<-toupper(rRNA.list)#make string and uppercase

P.Buff.or.Scale.Gain.highAneu(rRNA.list) # NS
P.Buff.or.Scale.Loss.highAneu(rRNA.list) # 5E-05 buffer
R.Buff.or.Scale.Gain.highAneu(rRNA.list) # 0.0047 scaling
R.Buff.or.Scale.Loss.highAneu(rRNA.list) # NS 

##INSTRUCTIONS: add list of genes you want to examine, and add group name and ID accordingly
## groupID will be added to .pdf title of the file. 
#setwd()

Test.List= Protein.Gain.Ribosome$RNA_Name
GroupName= "Ribosome genes (high aneuploid cells)"
GroupID= "Ribosome_HighPloidy"


## density plot PROTEIN GAIN 
pdf(file = paste0("plot.Protein.gain.density.", GroupID, "_HighCellPloidy.pdf"),
    width = 4, 
    height = 4)
# non-list genes
ggplot(CN.Diff.RNA.Prot_HighPloidy, 
       aes(x=Protein.Diff.Gain))+
  geom_density(data=CN.Diff.RNA.Prot_HighPloidy[!CN.Diff.RNA.Prot_HighPloidy$RNA_Name %in% Test.List,],
               color="black", size=1)+
  geom_vline(xintercept=mean(CN.Diff.RNA.Prot_HighPloidy[!CN.Diff.RNA.Prot_HighPloidy$RNA_Name %in% Test.List,]$Protein.Diff.Gain), 
             color="black", linetype="dotted", size=1)+
  #list genes
  geom_density(data=CN.Diff.RNA.Prot_HighPloidy[CN.Diff.RNA.Prot_HighPloidy$RNA_Name %in% Test.List,], 
               color="red", size=1)+
  geom_vline(xintercept=mean(CN.Diff.RNA.Prot_HighPloidy[CN.Diff.RNA.Prot_HighPloidy$RNA_Name %in% Test.List,]$Protein.Diff.Gain), 
             color="red", linetype="dotted", size=1)+
  
  ylab("Density") +
  xlab("Difference in protein expression upon chrm gain (High ploidy cells)")+
  ggtitle (GroupName)+
  coord_cartesian(xlim=c(-1,1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()
dev.off()

## density plot PROTEIN LOSS
pdf(file = paste0("plot.Protein.loss.density.", GroupID, "_HighCellPloidy.pdf"),
    width = 4, 
    height = 4)
ggplot(CN.Diff.xRNA.yProt.ThreeGroups, 
       aes(x=Protein.Diff.Loss))+
  geom_density(data=CN.Diff.RNA.Prot_HighPloidy[!CN.Diff.RNA.Prot_HighPloidy$RNA_Name %in% Test.List,], 
               color="black", size=1)+
  geom_vline(xintercept=mean(CN.Diff.RNA.Prot_HighPloidy[!CN.Diff.RNA.Prot_HighPloidy$RNA_Name %in% Test.List,]$Protein.Diff.Loss), 
             color="black", linetype="dotted", size=1)+
  
  geom_density(data=CN.Diff.RNA.Prot_HighPloidy[CN.Diff.RNA.Prot_HighPloidy$RNA_Name %in% Test.List,], 
               color="dodgerblue3", size=1)+
  geom_vline(xintercept=mean(CN.Diff.RNA.Prot_HighPloidy[CN.Diff.RNA.Prot_HighPloidy$RNA_Name %in% Test.List,]$Protein.Diff.Loss), 
             color="dodgerblue3", linetype="dotted", size=1)+
  
  ylab("Density") + 
  xlab("Difference in protein expression upon chrm loss (High ploidy cells)")+
  ggtitle (GroupName)+
  coord_cartesian(xlim=c(-1,1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()
dev.off()

### RNA
pdf(file = paste0("plot.RNA.gain.density.", GroupID, "_HighPloidy.pdf"),
    width = 4, 
    height = 4)
## density plot RNA GAIN 
ggplot(CN.Diff.xRNA.HighPloidy.unique, 
       aes(x=RNA.Diff.Gain))+
  geom_density(data=CN.Diff.xRNA.HighPloidy.unique[!CN.Diff.xRNA.HighPloidy.unique$RNA_Name %in% Test.List,], 
               color="black", size=1)+
  geom_vline(xintercept=mean(CN.Diff.xRNA.HighPloidy.unique[!CN.Diff.xRNA.HighPloidy.unique$RNA_Name %in% Test.List,]$RNA.Diff.Gain), 
             color="black", linetype="dotted", size=1)+
  
  geom_density(data=CN.Diff.xRNA.HighPloidy.unique[CN.Diff.xRNA.HighPloidy.unique$RNA_Name %in% Test.List,], 
               color="red", size=1)+
  geom_vline(xintercept=mean(CN.Diff.xRNA.HighPloidy.unique[CN.Diff.xRNA.HighPloidy.unique$RNA_Name %in% Test.List,]$RNA.Diff.Gain), 
             color="red", linetype="dotted", size=1)+
  ylab("Density") +
  xlab("Difference in RNA expression upon chrm gain (High ploidy cells)")+
  ggtitle (GroupName)+
  coord_cartesian(xlim=c(-1,1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()
dev.off()

## density plot RNA LOSS
pdf(file = paste0("plot.RNA.loss.density.", GroupID, "_HighPloidy.pdf"),
    width = 4, 
    height = 4)
ggplot(CN.Diff.xRNA.HighPloidy.unique, 
       aes(x=RNA.Diff.Loss))+
  geom_density(data=CN.Diff.xRNA.HighPloidy.unique[!CN.Diff.xRNA.HighPloidy.unique$RNA_Name %in% Test.List,], 
               color="black", size=1)+
  geom_vline(xintercept=mean(CN.Diff.xRNA.HighPloidy.unique[!CN.Diff.xRNA.HighPloidy.unique$RNA_Name %in% Test.List,]$RNA.Diff.Loss), 
             color="black", linetype="dotted", size=1)+
  
  geom_density(data=CN.Diff.xRNA.HighPloidy.unique[CN.Diff.xRNA.HighPloidy.unique$RNA_Name %in% Test.List,], 
               color="dodgerblue3", size=1)+
  geom_vline(xintercept=mean(CN.Diff.xRNA.HighPloidy.unique[CN.Diff.xRNA.HighPloidy.unique$RNA_Name %in% Test.List,]$RNA.Diff.Loss), 
             color="dodgerblue3", linetype="dotted", size=1)+
  
  ylab("Density") +
  xlab("Difference in RNA expression upon chrm loss (High ploidy cells)")+
  ggtitle (GroupName)+
  coord_cartesian(xlim=c(-1,1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()
dev.off()

print("End of graphs, Protein and RNA")


####
#### DID not end up using this data in paper
# Now corrolate High aneuploid cell expression differences and low aneuploid cell expression differences
#Step 1: prep data
# x= Low aneuploidy score cells, y= high aneuploidy score cells

High.Low.Aneuploidy.Difference<- merge(x=CN.Diff.RNA.Prot_LowPloidy, y=CN.Diff.RNA.Prot_HighPloidy, 
                                       by.x="Protein_ID", by.y="Protein_ID")


#correlate difference upon chrm gain per gene, between cells with low and high aneuploidy
ggplot(High.Low.Aneuploidy.Difference, 
       aes(x=High.Low.Aneuploidy.Difference$Protein.Diff.Gain.x, 
           y=High.Low.Aneuploidy.Difference$Protein.Diff.Gain.y))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  xlab("Difference in protein expression upon gain\nLow aneuploid cells")+
  ylab("Difference in protein expression upon gain\nhigh aneuploid cells")+
  theme_classic()+
  geom_hline(yintercept=0.0)+
  geom_vline(xintercept=0.0)+
  geom_smooth(method="lm", color="Red")+
  coord_cartesian(xlim=c(-1, 1), ylim=c(-1,1))+
  ggtitle("Chromosome gain, Protein")
#4x4
# plot.Protein.Gain.Diff.Low.high.aneuploidCells
cor.test(High.Low.Aneuploidy.Difference$Protein.Diff.Gain.x, 
         High.Low.Aneuploidy.Difference$Protein.Diff.Gain.y, method="pearson")
#cor= -0.059
#P = 0.009412


###
#correlate difference upon chrm loss per gene, between cells with low and high aneuploidy
# x= low aneuploid quartile cells only, y= high aneuploid quartile cells only
ggplot(High.Low.Aneuploidy.Difference, 
       aes(x=High.Low.Aneuploidy.Difference$Protein.Diff.Loss.x, 
           y=High.Low.Aneuploidy.Difference$Protein.Diff.Loss.y))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  #geom_point()+
  xlab("Difference in protein expression upon loss\nLow aneuploid cells")+
  ylab("Difference in protein expression upon loss\nhigh aneuploid cells")+
  theme_classic()+
  geom_hline(yintercept=0.0)+
  geom_vline(xintercept=0.0)+
  geom_smooth(method="lm", color="Red")+
  coord_cartesian(xlim=c(-1, 1), ylim=c(-1,1))+
  ggtitle("Chromosome loss, Protein")
#4x4
# plot.Protein.Loss.Diff.Low.high.aneuploidCells
cor.test(High.Low.Aneuploidy.Difference$Protein.Diff.Loss.x, 
             High.Low.Aneuploidy.Difference$Protein.Diff.Loss.y, 
         method="pearson")
#cor= 0.06
#P = 0.0007629

## RNA: 

#correlate difference upon chrm gain per gene, between cells with low and high aneuploidy
ggplot(High.Low.Aneuploidy.Difference, 
       aes(x=High.Low.Aneuploidy.Difference$RNA.Diff.Gain.x, 
           y=High.Low.Aneuploidy.Difference$RNA.Diff.Gain.y))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  xlab("Difference in RNA expression upon gain\nLow aneuploid cells")+
  ylab("Difference in RNA expression upon gain\nhigh aneuploid cells")+
  theme_classic()+
  geom_hline(yintercept=0.0)+
  geom_vline(xintercept=0.0)+
  geom_smooth(method="lm", color="Red")+
  coord_cartesian(xlim=c(-1, 1), ylim=c(-1,1))+
  ggtitle("Chromosome gain, RNA")
#4x4
# plot.RNA.Gain.Diff.Low.high.aneuploidCells
cor.test(High.Low.Aneuploidy.Difference$RNA.Diff.Gain.x, 
         High.Low.Aneuploidy.Difference$RNA.Diff.Gain.y, method="pearson")
#cor= 0.0017
#P = NS


###
#correlate difference upon chrm loss per gene, between cells with low and high aneuploidy
# x= low aneuploid quartile cells only, y= high aneuploid quartile cells only
ggplot(High.Low.Aneuploidy.Difference, 
       aes(x=High.Low.Aneuploidy.Difference$RNA.Diff.Loss.x, 
           y=High.Low.Aneuploidy.Difference$RNA.Diff.Loss.y))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  #geom_point()+
  xlab("Difference in RNA expression upon loss\nLow aneuploid cells")+
  ylab("Difference in RNA expression upon loss\nhigh aneuploid cells")+
  theme_classic()+
  geom_hline(yintercept=0.0)+
  geom_vline(xintercept=0.0)+
  geom_smooth(method="lm", color="Red")+
  coord_cartesian(xlim=c(-1, 1), ylim=c(-1,1))+
  ggtitle("Chromosome loss, RNA")
#4x4
# plot.RNA.Loss.Diff.Low.high.aneuploidCells
cor.test(High.Low.Aneuploidy.Difference$RNA.Diff.Loss.x, 
         High.Low.Aneuploidy.Difference$RNA.Diff.Loss.y, 
         method="pearson")
#cor= 0.2
#P < 2E-16


##### Step 6) Reviewer Response: Find percent significant per difference ammount #####

CN.Diff.percentSig<-CN.Diff.xRNA.yProt.ThreeGroups
CN.Diff.percentSig$Diff.Group.Gain<-NA
CN.Diff.percentSig$Diff.Group.Gain[CN.Diff.percentSig$Protein.Diff.Gain < -0.5] = "< -0.5"
CN.Diff.percentSig$Diff.Group.Gain[CN.Diff.percentSig$Protein.Diff.Gain >= -0.5 & CN.Diff.percentSig$Protein.Diff.Gain < -0.25 ] = "-0.5 to -0.25"
CN.Diff.percentSig$Diff.Group.Gain[CN.Diff.percentSig$Protein.Diff.Gain >= -0.25 & CN.Diff.percentSig$Protein.Diff.Gain < -0.1 ] = "-0.25 to -0.1"
CN.Diff.percentSig$Diff.Group.Gain[CN.Diff.percentSig$Protein.Diff.Gain >= -0.1 & CN.Diff.percentSig$Protein.Diff.Gain < 0.1 ] = "-0.1 to 0.1"
CN.Diff.percentSig$Diff.Group.Gain[CN.Diff.percentSig$Protein.Diff.Gain >= 0.1 & CN.Diff.percentSig$Protein.Diff.Gain < 0.25 ] = "0.1 to 0.25"
CN.Diff.percentSig$Diff.Group.Gain[CN.Diff.percentSig$Protein.Diff.Gain >= 0.25 & CN.Diff.percentSig$Protein.Diff.Gain < 0.5 ] = "0.25 to 0.5"
CN.Diff.percentSig$Diff.Group.Gain[CN.Diff.percentSig$Protein.Diff.Gain > 0.5 ] = "> 0.5"

CN.Diff.percentSig$Diff.Group.Gain<- factor(CN.Diff.percentSig$Diff.Group.Gain, levels=c("< -0.5", "-0.5 to -0.25", "-0.25 to -0.1", "-0.1 to 0.1",
                                                                                         "0.1 to 0.25", "0.25 to 0.5", "> 0.5"))

CN.Diff.percentSig$Diff.Group.Loss<-NA
CN.Diff.percentSig$Diff.Group.Loss[CN.Diff.percentSig$Protein.Diff.Loss < -0.5] = "< -0.5"
CN.Diff.percentSig$Diff.Group.Loss[CN.Diff.percentSig$Protein.Diff.Loss >= -0.5 & CN.Diff.percentSig$Protein.Diff.Loss < -0.25 ] = "-0.5 to -0.25"
CN.Diff.percentSig$Diff.Group.Loss[CN.Diff.percentSig$Protein.Diff.Loss >= -0.25 & CN.Diff.percentSig$Protein.Diff.Loss < -0.1 ] = "-0.25 to -0.1"
CN.Diff.percentSig$Diff.Group.Loss[CN.Diff.percentSig$Protein.Diff.Loss >= -0.1 & CN.Diff.percentSig$Protein.Diff.Loss < 0.1 ] = "-0.1 to 0.1"
CN.Diff.percentSig$Diff.Group.Loss[CN.Diff.percentSig$Protein.Diff.Loss >= 0.1 & CN.Diff.percentSig$Protein.Diff.Loss < 0.25 ] = "0.1 to 0.25"
CN.Diff.percentSig$Diff.Group.Loss[CN.Diff.percentSig$Protein.Diff.Loss >= 0.25 & CN.Diff.percentSig$Protein.Diff.Loss < 0.5 ] = "0.25 to 0.5"
CN.Diff.percentSig$Diff.Group.Loss[CN.Diff.percentSig$Protein.Diff.Loss > 0.5 ] = "> 0.5"

CN.Diff.percentSig$Diff.Group.Loss<- factor(CN.Diff.percentSig$Diff.Group.Loss, levels=c("< -0.5", "-0.5 to -0.25", "-0.25 to -0.1", "-0.1 to 0.1",
                                                                                         "0.1 to 0.25", "0.25 to 0.5", "> 0.5"))

CN.Diff.percentSig$RankP.Gain<-rank(CN.Diff.percentSig$Protein.Pvalue.Gain)
CN.Diff.percentSig$B.H.CriticalValue<-(0.05/length(CN.Diff.percentSig$Protein.Pvalue.Gain))*CN.Diff.percentSig$RankP.Gain
CN.Diff.percentSig$significant.Gain<- CN.Diff.percentSig$Protein.Pvalue.Gain< CN.Diff.percentSig$B.H.CriticalValue
CN.Diff.percentSig$significant.Gain<- factor(CN.Diff.percentSig$significant.Gain, levels=c("TRUE", "FALSE"))

CN.Diff.percentSig$RankP.Loss<-rank(CN.Diff.percentSig$Protein.Pvalue.Loss)
CN.Diff.percentSig$B.H.CriticalValue2<-(0.05/length(CN.Diff.percentSig$Protein.Pvalue.Loss))*CN.Diff.percentSig$RankP.Loss
CN.Diff.percentSig$significant.Loss<- CN.Diff.percentSig$Protein.Pvalue.Loss< CN.Diff.percentSig$B.H.CriticalValue2
CN.Diff.percentSig$significant.Loss<- factor(CN.Diff.percentSig$significant.Loss, levels=c("TRUE", "FALSE"))


## Plot Number of significant and percent significant per category of how much difference:
# number significant (B.H. corrected p-value)
ggplot(CN.Diff.percentSig, aes(x=Diff.Group.Loss, fill= significant.Loss ))+
  geom_bar(position= "stack")+
  xlab("Difference upon chromosome loss")+
  ylab("Number of proteins")+
  geom_vline(xintercept=c(2.5, 4.5), linetype="dotted", color="black", size=1)+
  scale_fill_manual(values=c("Gold2", "Black"))
# 7x4
# plot.Protein.Difference.Significant.Cat.Loss

ggplot(CN.Diff.percentSig, aes(x=Diff.Group.Gain, fill= significant.Gain ))+
  geom_bar(position= "stack")+
  xlab("Difference upon chromosome loss")+
  ylab("Number of proteins")+
  geom_vline(xintercept=c(3.5, 5.5), linetype="dotted", color="black", size=1)+
  scale_fill_manual(values=c("Gold2", "Black"))
# 7x4
# plot.Protein.Difference.Significant.Cat.Gain

## Percent significant (B.H. corrected p-value): 
ggplot(CN.Diff.percentSig, aes(x=Diff.Group.Loss, fill= significant.Loss ))+
  geom_bar(position= "fill")+
  xlab("Difference upon chromosome loss")+
  ylab("Number of proteins")+
  scale_fill_manual(values=c("Gold2", "Black"))

ggplot(CN.Diff.percentSig, aes(x=Diff.Group.Gain, fill= significant.Gain ))+
  geom_bar(position= "fill")+
  xlab("Difference upon chromosome gain")+
  ylab("Number of proteins")+
  scale_fill_manual(values=c("Gold2", "Black"))


# scatterplot log2 p-value and difference with significance (corrected)
ggplot(CN.Diff.percentSig, aes(x=Protein.Diff.Gain, y= -log2(Protein.Pvalue.Gain), 
                               color= significant.Gain ))+
  geom_point()+
  xlab("Difference upon chromosome gain")+
  ylab("-log2(p-value)")+
  geom_vline(xintercept=c(-0.1, .25), linetype="dotted", color="black", size=1)+
  scale_color_manual(values=c("Gold2", "Black"))
# 5x4
# plot.Protein.Difference.Significant.Scatter.Gain

ggplot(CN.Diff.percentSig, aes(x=Protein.Diff.Loss, y= -log2(Protein.Pvalue.Loss), 
                               color= significant.Loss ))+
  geom_point()+
  xlab("Difference upon chromosome loss")+
  ylab("-log2(p-value)")+
  geom_vline(xintercept=c(0.1, -.25), linetype="dotted", color="black", size=1)+
  scale_color_manual(values=c("Gold2", "Black"))
# 5x4
# plot.Protein.Difference.Significant.Scatter.Loss