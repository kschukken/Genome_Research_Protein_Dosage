##### Protein expression data correlation with cellular aneuploidy score ####
## Protein expression correlate with aneuploidy score (bp per chrm corrected for ploidy)
## 200803
## update: 210218- added gene-count based aneuploidy data
## Klaske Schukken

library(ggplot2)
library(tidyverse) 
library(readxl)
library("corrr")

DataFileLocation<- "Documents" # ! Set to location you want files to be written to


##### Step 1: get Data ####
## Get depmap.org protein expression data and cell line information data. 
#              https://depmap.org/portal/download/
## Normalized Protein Expression data 
## https://www.cell.com/cell/fulltext/S0092-8674(19)31385-6#secsectitle0190
## Nusinow et al. Cell, 2020, Quantitative proteomics of the cancer cell line encyclopedia

# get data and put into proper formatting
# setwd(DataFileLocation)

# Normalized protein expression per cell line
Protein_Expression.1<-read_excel("mmc2.xlsx", 
                               sheet= "Normalized Protein Expression")
Protein_Expression.1<-Protein_Expression.1[,c(1:426)] #delete empty collumns

#This will be filtered for genes that have a minimum of 10 cells per category

# test<- Protein_Expression.1[-grep("tr", Protein_Expression.1$Protein_Id), ]
# Delete all the potential but unvalidated proteins. "tr|". 

n <- Protein_Expression.1$Gene_Symbol   
Protein_Expression2<- as.data.frame(t(as.matrix(Protein_Expression.1[,-1] ))) #switch rows and collumns
colnames(Protein_Expression2) <- n

Protein_Expression2$Cell_Lines<-factor(rownames(Protein_Expression2)) 

for (i in 1:(length(Protein_Expression2)-1)){
  Protein_Expression2[,i]<- as.numeric(as.character(Protein_Expression2[,i]))
}# NEED AS.CHARACTER otherwise data is bad! NAs become numbers!!!

Protein_Expression2$Cell_Lines<-gsub("_TenPx.*","", Protein_Expression2$Cell_Lines)

#CCLE.names <-read_excel("mmc1.xlsx", 
#             sheet="Sample_Information")
#CCLE.names<-CCLE.names[,1:2] #get CCLE names and depmap ID correlated table
#CCLE.names$`CCLE Code`<-factor(CCLE.names$`CCLE Code`)

#Protein_effect_3<-merge(x= CCLE.names, y= Protein_Expression2, 
#                         by.x="CCLE Code", by.y="Cell_Lines", 
#                         sort = TRUE)

#setwd(DataFileLocation)

### get cell line information from depmap.org download file: 
CCLE.Depmap.names<-read.delim2("DepMap-2018q3-celllines.csv", 
                               dec=".", header = TRUE, sep=",")
CCLE.Depmap.name<-CCLE.Depmap.names[,1:2] #get CCLE names and depmap ID correlated table
# now add Depmap ID to Protein data, so depmap_ID matched CCLE cell name. 

Protein_effect_4<-merge(x= CCLE.Depmap.name, y= Protein_Expression2, 
                         by.x="CCLE_Name", by.y="Cell_Lines", 
                         sort = TRUE)

#write.csv(Protein_effect_4,
#          "Protein_Effect_notr.csv", row.names = TRUE)
#setwd()#set working directory 
#Protein_effect_4<-read.delim2("Protein_Effect_notr.csv", 
#                               dec=".", header = TRUE, sep=",")

##Add Aneuploidy Sore
## Get Score from Aneuploid_Score.R, get gene number data

Protein_Scores<-merge(x= Protein_effect_4, y= Score, 
                              by.x="Broad_ID", by.y="DepMap_ID", 
                              sort = TRUE)

#write.csv(Protein_Scores,
#          "Protein_Score_GeneAneuploid.csv", row.names = TRUE)


rRNAMet<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=9) 
CORUM.all<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=2) 
ERMem<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=8) 
Ribosome<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=3) 


### get only the cells in the lowest quantile (1/4th of cells with lowest aneuploidy)
#setwd()
#Protein_Scores<-read.delim2("Protein_Score_GeneAneuploid.csv", 
#                              dec=".", header = TRUE, sep=",")

quantile(Protein_Scores$Gene_ploidy_Score) #0, 1773, 2668, 3309, 6985
low.aneuploid<- subset(Protein_Scores, Gene_ploidy_Score<= 1773) #only 94 cells


####
## Now get all the lists of gene groups needed for gene group analysis: 
Oncogene.Bailey.list#84
TSG.Bailey.list#99
Protein.HSP
Protein.Ribosome
rRNAMet.list
CORUM.all3
ERMem.list


##### Step 2: examine raw data ####
# plot data to see how it's distributed. 

ggplot(Protein_Scores,aes(arm_Score)) +
  geom_bar()+
  theme_classic()

ggplot(Protein_Scores,aes(ploidy)) +
  geom_bar()+
  theme_classic()

ggplot(Protein_Scores,aes(bp_Score)) +
  geom_histogram(binwidth = 100000000)+
  theme_classic()

ggplot(Protein_Scores,aes(Gene_ploidy_Score)) +
  geom_histogram()+
  theme_classic()



##### Step 3: corrolate gene expression with aneuploidy scores ####
###NEXT:  corrolate protein expression with aneuploidy scores
#Using corrr, pearson corrolation. 
Protein.cor_effect<- Protein_Scores[sapply(Protein_Scores, 
                                                function(x) !is.factor(x))] %>% 
  correlate() %>% 
  focus(arm_Score, bp_Score, bp_ploidy_Score, Gene_Score, Gene_ploidy_Score)


Protein.cor_effect <- Protein.cor_effect[order(-Protein.cor_effect$arm_Score,
                                         -Protein.cor_effect$Gene_ploidy_Score),]

#write.csv(Protein.cor_effect,
#          DataFileLocation, row.names = TRUE)

#setwd(DataFileLocation)#set working directory to wherever you wrote above file to
Protein.cor_effect<-read.delim2("Protein.cor_effect_geneAneuploid.csv", 
                             dec=".", header = TRUE, sep=",")

length(Protein.cor_effect$X)#12314
Pro.neg.cor<-Protein.cor_effect[12304:12314,] #highest negative corrolation
Pro.neg.cor<- Pro.neg.cor[order(Pro.neg.cor$Gene_ploidy_Score, Pro.neg.cor$arm_Score),]
Pro.pos.cor<-Protein.cor_effect[1:10,] #highest positive corrolation



#####         Plot: examine correlation data #### 
###PLOTS
#Scatterplot of all genes & corrolation to aneuploid scores
ggplot(Protein.cor_effect[2:17311,], 
                          aes(x=bp_ploidy_Score, y=arm_Score))+
  geom_point()+
  ylab("Correlation: Protein expression & arm score") +
  xlab("Correlation: Protein expression & base-pair ploidy score")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()



ggplot(Protein.cor_effect[2:12314,], 
                     aes(x=Gene_ploidy_Score, y=arm_Score))+
  geom_point(color="black")+
  ylab("Correlation: Protein expression & Arm score") +
  xlab("Correlation: Protein expression & Gene ploidy score")+
  geom_smooth(method="lm", color="Red")+
  coord_cartesian(xlim=c(-0.31,0.31), ylim=c(-0.31,0.31))+
  theme_classic()
cor.test(Protein.cor_effect$arm_Score, Protein.cor_effect$Gene_ploidy_Score)
# corr= 0.896
# P < 2E-16

ggplot(Protein.cor_effect[2:17311,], aes(x=Gene_ploidy_Score, 
                                      y=bp_ploidy_Score))+
  geom_point(color="black")+
  ylab("Correlation: Protein expression & basepair ploidy score") +
  xlab("Correlation: Protein expression & gene ploidy score")+
  theme_classic()

##### Step 4: Get p-values for  correlation ####
#Protein.cor_effect.p<- Protein.gene_effect_Score[sapply(Protein.gene_effect_Score, function(x) !is.factor(x))] %>% 
#  corr.test() %>% 
#  focus(Arm_Score, bp_Score, bp_ploidy_Score)

#### find corrolation score for top 5 & bottom 5 genes
# after multiplying by number of genes
# use formula to find p-values
Pro.pos.cor
Pro.neg.cor
##Since many Protein experiments use few samples (high cor, but low pvalue) 
## I wanted to get the p-value for all samples, and add to corr dataframe.

# Step 1: get corrolation data between aneuploid score and expression per gene. 
df<-data.frame(Name=NA, Pvalue=NA)

for (i in 1:(length(Protein_Scores)-9)) {
  w=i+4
  t.test<-cor.test(Protein_Scores$Gene_ploidy_Score, 
                   Protein_Scores[,w], 
                   method = "pearson")
  pvalue<-t.test$p.value 
  name<-as.factor(colnames(Protein_Scores[w]))
  data<-data.frame(Name=name, Pvalue=pvalue)
  df <- rbind(df,data)
}

df$Name<-as.factor(df$Name)
#step 2: combine with corr data
Protein.corr.pvalue<-merge(x= Protein.cor_effect, y= df, 
                        by.x="rowname", by.y="Name", 
                        sort = TRUE)

#####         Plot: corr v. pvalue ####

## Significance scatter plot:
# want to add "significant" TRUE-FALSE statement
# to do Benjamini Hoecht correction, step 1 rank by p-value. 
# critical value=  (rank/#tests)*max p-value
# ex: ( 200/12311 ) *0.05
# if critical value is less than p-value, it's significant
Protein.corr.pvalue$RankPvalue<-rank(Protein.corr.pvalue$Pvalue)
Protein.corr.pvalue$B.H.CriticalValue<-(0.05/length(Protein.corr.pvalue$Pvalue))*Protein.corr.pvalue$RankPvalue
Protein.corr.pvalue$significant<- Protein.corr.pvalue$Pvalue< Protein.corr.pvalue$B.H.CriticalValue

plot.cor.pvalue<- ggplot(Protein.corr.pvalue2, 
                         aes(x=Gene_ploidy_Score, y=-log2(Pvalue), color=significant))+
  geom_point()+
  ylab("-log2 (P-Value)") +
  xlab("Correlation: Protein expression & aneuploidy score")+
  ggtitle ("Protein expression correlated \n with aneuploidy score")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  scale_color_manual(values=c("black", "goldenrod"))
plot.cor.pvalue


# Density plot
subset.significant <- subset(Protein.corr.pvalue2,significant==TRUE)
ggplot(Protein.corr.pvalue2, 
                         aes(x=Gene_ploidy_Score))+
  geom_density(fill="black", binwidth=0.025)+
  ylab("-log2 (P-Value)") +
  xlab("Correlation: Protein expression & aneuploidy score")+
  ggtitle ("Protein expression correlated with aneuploid score")+
  theme(plot.title = element_text(hjust = 0.5))+
  coord_cartesian(xlim=c(-0.3,0.3))+
  theme_classic()

##### Step 5: Sort by p-value, plot top 2 positive or negative correlated proteins ####
Protein.corr.pvalue2 <- Protein.corr.pvalue[order(Protein.corr.pvalue$Pvalue,
                                            -Protein.corr.pvalue$Gene_ploidy_Score),]

#setwd(DataFileLocation)
write.csv(Protein.corr.pvalue2,
          "Protein_Corr_pvalue_GeneAneuploid.csv", row.names = TRUE)

Protein.corr.pvalue2<-read.delim2("Protein_Corr_pvalue_GeneAneuploid.csv", 
                              dec=".", header = TRUE, sep=",")


#Get significant p-value, positive correlation with aneuploidy
Protein.corr.pos.sig <- Protein.corr.pvalue[order(Protein.corr.pvalue$Pvalue),]
Protein.corr.pos.sig <- subset(Protein.corr.pos.sig, Protein.corr.pos.sig$significant==TRUE)
Protein.corr.pos.sig <- subset(Protein.corr.pos.sig, Protein.corr.pos.sig$Gene_ploidy_Score >0)

# top 2 most sig up: HSPA6, HSPA1A, both Heat Shock proteins. 4 of 31 HSP were upregulated, none downreg

#write.csv(Protein.corr.pos.sig,
#          DataFileLocation, row.names = TRUE)

#Get significant p-value, negative correlation with aneuploidy
Protein.corr.neg.sig <- Protein.corr.pvalue[order(Protein.corr.pvalue$Pvalue),]
Protein.corr.neg.sig <- subset(Protein.corr.neg.sig, Protein.corr.neg.sig$significant==TRUE)
Protein.corr.neg.sig <- subset(Protein.corr.neg.sig, Protein.corr.neg.sig$Gene_ploidy_Score <0)

#write.csv(Protein.corr.neg.sig.csv,
#          DataFileLocation, row.names = TRUE)

#top 2 most signif down: RPL22L1, RNF138


###Now to plot correlation between aneuploidy score and top genes

#top p-value significant correlations
#all_proteins<- colnames(Protein_Scores)
#all_proteins<- gsub("_HUMAN","",all_proteins)

##Plot expression by cellular aneuploidy score, for pos/neg corr genes: 
ggplot(Protein_Scores, 
                  aes(x=Gene_ploidy_Score, y=RNF138))+ #
  geom_point(color="dodgerblue3", size=2)+
  ylab("Protein expression: RNF138") + 
  xlab("Gene ploidy score") +
  geom_smooth(color="black", method="lm") +
  #geom_text(x=6000, y=-2, label="p-value= 1.3E-09")+
  #geom_text(x=6000, y=-2.5, label="cor= 0.22")+ 
  theme_classic()
# plot.Protein.aneuscore.Express.RNF138_blue
# 4x4
## COLORS: Red for pos, dodgerblue3 for neg


##### Step 6: get gene groups, test sig more/loss corr than mean proteins ####

##Get all HSPA HSP70 proteins:
Protein.HSP <- subset(Protein.corr.pvalue2, grepl( "HSPA", Protein.corr.pvalue2$rowname, fixed = TRUE))
Protein.HSP <- Protein.HSP[order(Protein.HSP$Pvalue),]
Protein.HSP.NOT <- Protein.corr.pvalue2[!Protein.corr.pvalue2$rowname %in% Protein.HSP$rowname,]# 

t.test(Protein.HSP$Gene_ploidy_Score, Protein.HSP.NOT$Gene_ploidy_Score)


## Get all Ribosomal proteins: 
#Protein.corr.pvalue3<-Protein.corr.pvalue2[order(Protein.corr.pvalue2$rowname),]
#Protein.Ribosome <- subset(Protein.corr.pvalue3, grepl( "RPS", Protein.corr.pvalue3$rowname, fixed = TRUE))
#Protein.Ribosome <- Protein.Ribosome[33:75,] #get rid of mitochondrial ribosome genes
#Protein.Ribosome2 <- subset(Protein.corr.pvalue3, grepl( "RPL", Protein.corr.pvalue3$rowname, fixed = TRUE))
#Protein.Ribosome2 <- Protein.Ribosome2[46:93,] #get rid of mitochondrial ribosome genes
#Protein.Ribosome<- rbind(Protein.Ribosome, Protein.Ribosome2) #this gives 91 genes
Ribosome.list<-Ribosome$Gene_Symbol
Protein.Ribosome<- Protein.corr.pvalue2[Protein.corr.pvalue2$rowname %in% Ribosome.list,]# 
Protein.Ribosome.NOT <- Protein.corr.pvalue2[!Protein.corr.pvalue2$rowname %in% Ribosome.list,]# 

t.test(Protein.Ribosome$Gene_ploidy_Score, Protein.Ribosome.NOT$Gene_ploidy_Score)


## Get all rRNA processing proteins GO:0006364
## 8 downregulated and 1 upregulated rRNA processing
#rRNAMet<- read.delim2("GO_term_rRNAprocessing.txt",
#                     dec=".", header = TRUE, sep="\t", row.names=NULL) 
rRNAMet.list<-unique(rRNAMet$Symbol)#get gene names
rRNAMet.list<-toupper(rRNAMet.list)#make string and uppercase. 209 genes
Protein.aneuScore.rRNAMet <- Protein.corr.pvalue2[Protein.corr.pvalue2$rowname %in% rRNAMet.list,]# 
Protein.aneuScore.rRNAMet.NOT <- Protein.corr.pvalue2[!Protein.corr.pvalue2$rowname %in% rRNAMet.list,]#

t.test(Protein.aneuScore.rRNAMet$Gene_ploidy_Score, Protein.aneuScore.rRNAMet.NOT$Gene_ploidy_Score)
# very significantly downregulated


#### CORUM Protein complexes
CORUM.all2<- subset(CN.Diff.xRNA.yProt.ThreeGroups2, Uniprot_Acc %in% CORUM.all$Uniprot_ID)
CORUM.all3<-CORUM.all2$RNA_Name
Protein.aneuscore.Corum <- Protein.corr.pvalue2[Protein.corr.pvalue2$rowname %in% CORUM.all3,]# 
Protein.aneuscore.Corum.NOT <- Protein.corr.pvalue2[!Protein.corr.pvalue2$rowname %in% CORUM.all3,]# 

t.test(Protein.aneuscore.Corum$Gene_ploidy_Score, Protein.aneuscore.Corum.NOT$Gene_ploidy_Score)
# CORUM genes are very mildly negatively correlated with aneuploidy. p=3.532e-06 


####intrinsic component of endoplasmic reticulum membrane, GO:0031227
#ERMem<- read.delim2("GO_term_ER_Membrane.txt",
#                    dec=".", header = TRUE, sep="\t", row.names=NULL) 
ERMem.list<-unique(ERMem$Symbol)#get gene names
ERMem.list<-toupper(ERMem.list)#make string and uppercase. 209 genes
Protein.aneuScore.ERMem <- Protein.corr.pvalue2[Protein.corr.pvalue2$rowname %in% ERMem.list,]# 
Protein.aneuScore.ERMem.NOT <- Protein.corr.pvalue2[!Protein.corr.pvalue2$rowname %in% ERMem.list,]# 

t.test(Protein.aneuScore.ERMem$Gene_ploidy_Score, Protein.aneuScore.ERMem.NOT$Gene_ploidy_Score)
# p= 5E-05
#good, upreg 


### Other gene groups examined, feel free to ignore this part: 

## Get all DNA double strand break repair proteins, go term GO:0006302
## 11 upregulated, 4 downregulated
#DNABreakRepair<- read.delim2("GO_term_DNADoubleStrandBreakRepair.txt",
#                             dec=".", header = TRUE, sep="\t", row.names=NULL) 
#DNABreakRepair.list<-unique(DNABreakRepair$MGI)#get gene names
#DNABreakRepair.list<-toupper(DNABreakRepair.list)#make string and uppercase. 243 genes
#Protein.aneuScore.DNABreakRepair <- Protein.corr.pvalue2[Protein.corr.pvalue2$rowname %in% DNABreakRepair.list,]#183 found back
#Protein.aneuScore.DNABreakRepair.NOT <- Protein.corr.pvalue2[!Protein.corr.pvalue2$rowname %in% DNABreakRepair.list,]# 190 genes found back

#t.test(Protein.aneuScore.DNABreakRepair$Gene_ploidy_Score, Protein.aneuScore.DNABreakRepair.NOT$Gene_ploidy_Score)
#0.0092, upregulated


## Get all Nonsense mediated decay proteins GO:0000184
## two significantly downregulated. most no change
#NMD<- read.delim2("GO_term_NonsenseMediatedDecay.txt",
#                             dec=".", header = TRUE, sep="\t", row.names=NULL) 
#NMD.list<-unique(NMD$MGI)#get gene names
#NMD.list<-toupper(NMD.list)#make string and uppercase. 36 genes
#Protein.aneuScore.NMD <- Protein.corr.pvalue[Protein.corr.pvalue$rowname %in% NMD.list,]# 31 genes found back


## Get all Metabolism of RNA proteins GO:0016070
## about equally up and down regulated
#RNAMet<- read.delim2("GO_term_RNAMetabolism.txt",
#                  dec=".", header = TRUE, sep="\t", row.names=NULL) 
#RNAMet.list<-unique(RNAMet$MGI)#get gene names
#RNAMet.list<-toupper(RNAMet.list)#make string and uppercase. 4258 genes
#Protein.aneuScore.RNAMet <- Protein.corr.pvalue[Protein.corr.pvalue$rowname %in% RNAMet.list,]# 2701 genes found back


## Plasma Membrane: slight upregulation in aneuploid, slight.
#PlasmaMembrane.list 
#Protein.aneuScore.PM <- Protein.corr.pvalue2[Protein.corr.pvalue2$rowname %in% PlasmaMembrane.list,]# 
#Protein.aneuScore.PM.NOT <- Protein.corr.pvalue2[!Protein.corr.pvalue2$rowname %in% PlasmaMembrane.list,]# 
#
#t.test(Protein.aneuScore.PM$Gene_ploidy_Score, Protein.aneuScore.PM.NOT$Gene_ploidy_Score)

### cell cell singaling , NS
#CCSignal.list
#Protein.aneuScore.CCSignal <- Protein.corr.pvalue[Protein.corr.pvalue$rowname %in% CCSignal.list,]# 

### Immune system , NS
#Immune.list
#Protein.aneuScore.Immune <- Protein.corr.pvalue[Protein.corr.pvalue$rowname %in% Immune.list,]# 

### TF, significantly upregulated *
#TF.list
#Protein.aneuScore.TF <- Protein.corr.pvalue2[Protein.corr.pvalue2$rowname %in% TF.list,]#
#Protein.aneuScore.TF.NOT <- Protein.corr.pvalue2[!Protein.corr.pvalue2$rowname %in% TF.list,]# 
#
#t.test(Protein.aneuScore.TF$Gene_ploidy_Score, Protein.aneuScore.TF.NOT$Gene_ploidy_Score)

### DNA repair, slight upregulation. NS
#DNARepair.list
#Protein.aneuScore.DNArepair <- Protein.corr.pvalue2[Protein.corr.pvalue2$rowname %in% DNARepair.list,]# 
#Protein.aneuScore.DNArepair.NOT <- Protein.corr.pvalue2[!Protein.corr.pvalue2$rowname %in% DNARepair.list,]# 
#
#t.test(Protein.aneuScore.DNArepair$Gene_ploidy_Score, Protein.aneuScore.DNArepair.NOT$Gene_ploidy_Score)

### Cell adhesion 
#CellAdhesion.list
#Protein.aneuScore.CellAdhesion <- Protein.corr.pvalue[Protein.corr.pvalue$rowname %in% CellAdhesion.list,]# 

### Tight Junction
#TightJunction.list
#Protein.aneuScore.TJ <- Protein.corr.pvalue[Protein.corr.pvalue$rowname %in% TightJunction.list,]# 

## Cellular response to heat, GO:0034605, NS 
#CellHeat.list
#Protein.aneuScore.cellheat <- Protein.corr.pvalue[Protein.corr.pvalue$rowname %in% CellHeat.list,]# 

# MRN complex (MRE11-RAD50-NBN complex), GO:0030870
#MRNcomplex.list
#Protein.aneuScore.MRN <- Protein.corr.pvalue2[Protein.corr.pvalue2$rowname %in% MRNcomplex.list,]# 
#Protein.aneuScore.MRN.NOT <- Protein.corr.pvalue2[!Protein.corr.pvalue2$rowname %in% MRNcomplex.list,]#

#t.test(Protein.aneuScore.MRN$Gene_ploidy_Score, Protein.aneuScore.MRN.NOT$Gene_ploidy_Score)

# Telomere lengthening. telomere maintenance via telomere lengthening, GO:0010833
#TeloLength.list
#Protein.aneuScore.TeloLength <- Protein.corr.pvalue[Protein.corr.pvalue$rowname %in% TeloLength.list,]# 



#### intrinsic component of plasma membrane, GO:0031226
#IntMem<- read.delim2("GO_term_Intrinsic.Plasma.Membrane.txt",
#                     dec=".", header = TRUE, sep="\t", row.names=NULL) 
#IntMem.list<-unique(IntMem$MGI)#get gene names
#IntMem.list<-toupper(IntMem.list)#make string and uppercase. 209 genes
#Protein.aneuScore.IntMem <- Protein.corr.pvalue2[Protein.corr.pvalue2$rowname %in% IntMem.list,]# 
#Protein.aneuScore.IntMem.NOT <- Protein.corr.pvalue2[!Protein.corr.pvalue2$rowname %in% IntMem.list,]# 
#
#t.test(Protein.aneuScore.IntMem$Gene_ploidy_Score, Protein.aneuScore.IntMem.NOT$Gene_ploidy_Score)
#good, upreg, but a lot of genes, hard to see trend. 



#####         TSG and OG subgroups ####
##Tumor supressor gene (TSG) and oncogene (OG) list
## Bailey Corr Oncogene: NS
## Bailey Corr TSG: NS

# Bailey et al. Comprehensive Characterization of Cancer Driver Genes and Mutations, cell, 2018
#Oncogene.Bailey.list#84
#TSG.Bailey.list#99

P.Corr.Oncogene.Bailey <- Protein.corr.pvalue2[Protein.corr.pvalue2$rowname %in% Oncogene.Bailey.list,]# 
P.Corr.TSG.Bailey <- Protein.corr.pvalue2[Protein.corr.pvalue2$rowname %in% TSG.Bailey.list,]# 

#Now get all genes EXCEPT for oncogene/TSG. for t-test purposes
P.Corr.Oncogene.Bailey.NOT <- Protein.corr.pvalue2[!Protein.corr.pvalue2$rowname %in% Oncogene.Bailey.list,]# 
P.Corr.TSG.Bailey.NOT <- Protein.corr.pvalue2[!Protein.corr.pvalue2$rowname %in% TSG.Bailey.list,]# 

## T-test between oncogenes/TSG vs all others: more likely to be pos or neg corr with aneuploidy? 
t.test(P.Corr.Oncogene.Bailey$Gene_ploidy_Score, P.Corr.Oncogene.Bailey.NOT$Gene_ploidy_Score)
#Oncogene NS. not correlated with aneuploidy score
t.test(P.Corr.TSG.Bailey$Gene_ploidy_Score, P.Corr.TSG.Bailey.NOT$Gene_ploidy_Score)
#TSG NS. not correlated with aneuploidy score



#####         Plot scatterplot correlation with aneuploidy score, per group ####
# CORUM  mildly negatively correlated with aneuploidy. p=3.532e-06 
# HSPA NS
# rRNA process neg corr p<2.2E-16
# Ribosome neg cor p<2.2E-16
# Intrinsic plasma membrane p=9E-07
# intrinsic component of endoplasmic reticulum membrane, 5E-05

##Plot scatterplot log2 correlation with aneuscore. 
## ! add correlation data subset with only gene group (made above) and name it: 
P.data.subset= Protein.aneuScore.ERMem
P.data.name = "Intrinsic component of endoplasmic reticulum membrane"

ggplot(Protein.corr.pvalue2,
                         aes(x=Gene_ploidy_Score, y=-log2(Pvalue)))+
  geom_point(size=1, color="black")+
  geom_point(size=2, color="red", data=P.data.subset,
             aes(x=Gene_ploidy_Score, y=-log2(Pvalue) ))+
  ylab("-log2 (P-Value)") +
  xlab("corr coefficient")+
  ggtitle (P.data.name)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  scale_color_manual(values=c("black", "red"))+
  scale_x_continuous(breaks=c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3), limits=c(-0.3, 0.3))+
  scale_y_continuous(limits=c(0, 35))

# 4x4 size
# plot.Protein.Corr.Pvalue.all.ERMembrane
# COLOR: "dodgerblue3" for negatively correlated, "red" for positively correlated. 



##### Step 7: Investigate low / high aneuploid cell quantile ####

### get only the cells in the lowest quantile (1/4th of cells with lowest aneuploidy)
#setwd(DataFileLocation)
#Protein_Scores<-read.delim2("Protein_Score_GeneAneuploid.csv", 
#                            dec=".", header = TRUE, sep=",")

#quantile(Protein_Scores$Gene_ploidy_Score) #0, 1773, 2668, 3309, 6985
#low.aneuploid<- subset(Protein_Scores, Gene_ploidy_Score<= 1773) #only 94 cells
#high.aneuploid<- subset(Protein_Scores, Gene_ploidy_Score> 3309) #only 94 cells
#low.aneuploid.cells<-low.aneuploid$Broad_ID
#high.aneuploid.cells<-high.aneuploid$Broad_ID
### This is continued in Protein_RNA_Expression.PerCell_V2.R, and Protein_RNA_1d.plot_v2.R
### where I calculate difference upon chrm arm gain or loss, and then make 1d plots per categories. 

##Make cell line info
Cell_info<- Protein_Scores[,2:3]
Cell_info$Cellular_Aneuploidy_Score<-Protein_Scores$Gene_ploidy_Score
Cell_info$Cellular_Aneuploidy_Score_NoPloidyCorrection<-Protein_Scores$Gene_Score
Cell_info$Lowest_Aneuploidy_Quartile<-Cell_info$Cellular_Aneuploidy_Score <= 1773
Cell_info$Highest_Aneuploidy_Quartile<-Cell_info$Cellular_Aneuploidy_Score > 3309
Cell_info<- unique(Cell_info)

write.csv(Cell_info, 
          file =paste("Cell.Info.AneuploidyScore.csv", sep=','), 
          row.names = TRUE)

### NEXT:  corrolate protein expression with aneuploidy scores

# Using corrr, pearson corrolation. 
Protein.cor_effect.low<- low.aneuploid[sapply(low.aneuploid, 
                                           function(x) !is.factor(x))] %>% 
  correlate() %>% 
  focus(Gene_ploidy_Score)


Protein.cor_effect.low <- Protein.cor_effect.low[order(-Protein.cor_effect.low$Gene_ploidy_Score,
                                               -Protein.cor_effect.low$Gene_ploidy_Score),]


length(Protein.cor_effect.low$Gene_ploidy_Score)#12762, of which 75=NA, 12687
Protein.cor_effect.low[12677:12687,] #highest negative corrolation
Protein.cor_effect.low[1:10,] #highest positive corrolation



#####         Get p-values for  correlation #####
#Protein.cor_effect.p<- Protein.gene_effect_Score[sapply(Protein.gene_effect_Score, function(x) !is.factor(x))] %>% 
#  corr.test() %>% 
#  focus(Arm_Score, bp_Score, bp_ploidy_Score)

#### find corrolation score for top 5 & bottom 5 genes
# after multiplying by number of genes
# use formula to find p-values
Pro.pos.cor
Pro.neg.cor
##Since many Protein experiments use few samples (high cor, but low pvalue) 
## I wanted to get the p-value for all samples, and add to corr dataframe.

# Step 1: get corrolation data between aneuploid score and expression per gene. ***
df<-data.frame(Name=NA, Pvalue=NA)

for (w in 4:(length(low.aneuploid)-75)) {
  if ( sum(is.na(low.aneuploid[,w])) < 81 ) {
  t.test<-cor.test(as.numeric(low.aneuploid$Gene_ploidy_Score), 
                   as.numeric(low.aneuploid[,w]), 
                   method = "pearson")
  pvalue<-t.test$p.value 
  name<-as.factor(colnames(low.aneuploid[w]))
  data<-data.frame(Name=name, Pvalue=pvalue)
  df <- rbind(df,data)
}
}

df$Name<-as.factor(df$Name)
#step 2: combine with corr data
Protein.cor_effect.Pvalue.low<-merge(x= Protein.cor_effect.low, y= df, 
                           by.x="rowname", by.y="Name", 
                           sort = TRUE) #11523 genes left

#####         Plot: corr v. pvalue ####

## Significance scatter plot:
# want to add "significant" TRUE-FALSE statement
# to do Benjamini Hoechberg correction, step 1 rank by p-value. 
# critical value=  (rank/#tests)*max p-value
# ex: ( 200/12311 ) *0.05
# if critical value is less than p-value, it's significant
Protein.cor_effect.Pvalue.low$RankPvalue<-rank(Protein.cor_effect.Pvalue.low$Pvalue)
Protein.cor_effect.Pvalue.low$B.H.CriticalValue<-(0.05/length(Protein.cor_effect.Pvalue.low$Pvalue))*Protein.cor_effect.Pvalue.low$RankPvalue
Protein.cor_effect.Pvalue.low$significant<- Protein.cor_effect.Pvalue.low$Pvalue< Protein.cor_effect.Pvalue.low$B.H.CriticalValue

ggplot(Protein.cor_effect.Pvalue.low, 
                         aes(x=Gene_ploidy_Score, y=-log2(Pvalue), color=significant))+
  geom_point()+
  ylab("-log2 (P-Value)") +
  xlab("Correlation: Protein expression & aneuploidy score")+
  ggtitle ("Protein expression correlated \n with aneuploidy score")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  scale_color_manual(values=c("black", "goldenrod"))



#####         Sort by p-value, plot top 2 positive or negative correlated proteins ####
Protein.cor_effect.Pvalue.low2 <- Protein.cor_effect.Pvalue.low[order(Protein.cor_effect.Pvalue.low$Pvalue,
                                                  -Protein.cor_effect.Pvalue.low$Gene_ploidy_Score),]



#Get significant p-value, positive correlation with aneuploidy
Protein.corr.pos.sig <- Protein.cor_effect.Pvalue.low[order(Protein.cor_effect.Pvalue.low$Pvalue),]
Protein.corr.pos.sig <- subset(Protein.corr.pos.sig, Protein.corr.pos.sig$significant==TRUE)
Protein.corr.pos.sig <- subset(Protein.corr.pos.sig, Protein.corr.pos.sig$Gene_ploidy_Score >0)

# top 2 most sig up: RICTOR, RB1CC1


#Get significant p-value, negative correlation with aneuploidy
Protein.corr.neg.sig <- Protein.cor_effect.Pvalue.low[order(Protein.cor_effect.Pvalue.low$Pvalue),]
Protein.corr.neg.sig <- subset(Protein.corr.neg.sig, Protein.corr.neg.sig$significant==TRUE)
Protein.corr.neg.sig <- subset(Protein.corr.neg.sig, Protein.corr.neg.sig$Gene_ploidy_Score <0)

# top 2 most sig neg correlation: CROCC and MGMT


###Now to plot correlation between aneuploidy score and top genes

Ribosome.list

##Plot expression by cellular aneuploidy score, for pos/neg corr genes: 
ggplot(low.aneuploid, 
       aes(x=Gene_ploidy_Score, y=RNF138))+ #
  geom_point(color="dodgerblue3", size=2)+
  ylab("Protein expression: RPL22L1") + 
  xlab("Gene ploidy score") +
  geom_smooth(color="black", method="lm") +
  #geom_text(x=6000, y=-2, label="p-value= 1.3E-09")+
  #geom_text(x=6000, y=-2.5, label="cor= 0.22")+ 
  theme_classic()
# plot.Protein.aneuscore.Express.RNF138_blue
# 4x4
## COLORS: Red for pos, dodgerblue3 for neg



#####         get gene groups, test sig more/loss corr than mean proteins  #####

## Get all Ribosomal proteins: 
Protein.cor_effect.Pvalue.low2<-Protein.cor_effect.Pvalue.low[order(Protein.cor_effect.Pvalue.low$rowname),]
Protein.Ribosome <- subset(Protein.cor_effect.Pvalue.low2, grepl( "RPS", Protein.cor_effect.Pvalue.low2$rowname, fixed = TRUE))
Protein.Ribosome <- Protein.Ribosome[34:76,] #get rid of mitochondrial ribosome genes
Protein.Ribosome2 <- subset(Protein.cor_effect.Pvalue.low2, grepl( "RPL", Protein.cor_effect.Pvalue.low2$rowname, fixed = TRUE))
Protein.Ribosome2 <- Protein.Ribosome2[46:93,] #get rid of mitochondrial ribosome genes
Protein.Ribosome<- rbind(Protein.Ribosome, Protein.Ribosome2) #this gives 91 genes
Ribosome.list<- Protein.Ribosome$rowname
Protein.Ribosome<- Protein.cor_effect.Pvalue.low2[Protein.cor_effect.Pvalue.low2$rowname %in% Ribosome.list,]# 
Protein.Ribosome.NOT <- Protein.cor_effect.Pvalue.low2[!Protein.cor_effect.Pvalue.low2$rowname %in% Ribosome.list,]# 

t.test(Protein.Ribosome$Gene_ploidy_Score, Protein.Ribosome.NOT$Gene_ploidy_Score)
# ribosomes upregulated in low aneuploid cells? 
# P-value = 3E-09
# 0.08 to 0.01


##  Plot scatterplot correlation with aneuploidy score, per group 

##Plot scatterplot log2 correlation with aneuscore. 
## ! add correlation data subset with only gene group (made above) and name it: 
P.data.subset= Protein.Ribosome
P.data.name = "Ribosomes_lowAneuploidCellsOnly"

ggplot(Protein.cor_effect.Pvalue.low2,
       aes(x=Gene_ploidy_Score, y=-log2(Pvalue)))+
  geom_point(size=1, color="black")+
  geom_point(size=2, color="red", data=P.data.subset,
             aes(x=Gene_ploidy_Score, y=-log2(Pvalue) ))+
  ylab("-log2 (P-Value)") +
  xlab("corr coefficient")+
  ggtitle (P.data.name)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  scale_color_manual(values=c("black", "red"))
  #scale_x_continuous(breaks=c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3), limits=c(-0.3, 0.3))+
  #scale_y_continuous(limits=c(0, 35))

# 4x4 size
# plot.Protein.Corr.Pvalue.all.ERMembrane
# COLOR: "dodgerblue3" for negatively correlated, "red" for positively correlated. 
