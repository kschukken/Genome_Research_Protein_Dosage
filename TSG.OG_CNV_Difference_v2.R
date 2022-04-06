#####  Tumor Supressor Ggene (TSG) and Oncogene (OG) gene copy number differences ####
### 210429
## Gene copy number correlation with expression
## for Oncogenes and Tumor Supressor Genes
## not Chromosome copy number, but GENE copy number. so we look at amplifications
## and other differences between whole chrm gain/loss and gene copy number changes
## By Klaske Schukken

## Note: you may need to change the name of the gene you want to examine occationally,
## because some oncogenes and tumor supressor genes are two protein varaints, so
## you'll have to edit the code slightly to make up for this, but not to difficult. 
## just find new name of $DNA.CopyNum

library('ggplot2')
library('tidyverse')
library('readxl')
library('reshape2')
library('BBmisc')
library('dplyr')
library("ROCR")
library("readODS")


## Protein expression data: 
## https://www.cell.com/cell/fulltext/S0092-8674(19)31385-6#secsectitle0190
## Nusinow et al. Cell, 2020, Quantitative proteomics of the cancer cell line encyclopedia
## Chromosome arm data: Cohen-Sharir 2021 Nature paper 
## RNA expression data:  CCLE depmap data.

## CCLE_gene_cn.csv
## Copy number variation on gene level
# CCLE gene copy number are log2(1+n) set relative to basal ploidy of the cell
# categories: Loss (<0.75), neutral(0.75-1.25), gain(1.25-2), multi-gain(2+)
# triploidy in diploid cell= Log2(3+1)/log2(2+1)= 1.26
# triploidy in tetraploid cell= Log2(3+1)/log2(4+1)= 0.86
# diploidy in diploid cell = Log2(2+1)/log2(2+1)= 1

### Steps to get CNV correlation with Oncogene & TSG expression
## Step 1: get Data
#   Bailey et al: Oncogene and TSG gene list
#   Protein expression
#   CNV per gene
#   Mean difference per gene upon gain/loss
#   Correlate Gene expression and CNV data
## Step 2: Plot OG and TSG difference upon aneuploidy category
#   boxplot average OG on gain, TSG on loss (also average genes)
#   barplot average expression difference per gene
## Step 3: Plot OG and TSG difference upon Copy number variation
#   scatterplot of CNV/ploidy vs expression per gene
#   categorical: lost, neutral, one-gain, multi-gain
## Step 4: Get differences in expression between gene copy number increase/decrease vs neutral
#   Categorical gene copy number: lost, neutral and gain
## Step 5: compare Gene CNV differences to chrm CNV diffference
#   correlate differnce in expression upon chrm and gene gain v neutral (& loss vs neutral)
#   find what percent of genes in chrm gain are not gained at gene level
## Step 6: Does gene amplification affect gene buffering in dataset? 


#####  Step 1: get data #####

### Set below directory to where you want your files to be saved: 
FigureResultsDirectory= "/Documents"
setwd(FigureResultsDirectory)

##Tumor supressor gene (TSG) and oncogene (OG) list
## Download Bailey et al. 2018 supplemental figure 8, upload table 1. example below: 
# Bailey et al. Comprehensive Characterization of Cancer Driver Genes and Mutations, cell, 2018

TSG.OG<-read_xlsx("NIHMS948705_S8_Bailey2019.xlsx", sheet = 2) 
Oncogene.Bailey<- subset(TSG.OG, TSG.OG$`Tumor suppressor or oncogene prediction (by 20/20+)`=="oncogene")
Oncogene.Bailey.list<-unique(Oncogene.Bailey$Gene) #84

TSG.Bailey<- subset(TSG.OG, TSG.OG$`Tumor suppressor or oncogene prediction (by 20/20+)`=="tsg")
TSG.Bailey.list<-unique(TSG.Bailey$Gene) #99


### Protein / RNA expression difference upon chrm arm gain/loss, per gene 
## Get aneuploidy difference data based on min 10 cells/gene,  
# From Protein_RNA_expression.PerCell_v2.R
setwd()
CN.Diff.xRNA.yProt.ThreeGroups <- read.csv("RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points_3cat.csv", header=TRUE) #9,414 genes
CN.Diff.xRNA.yProt.ThreeGroups<-CN.Diff.xRNA.yProt.ThreeGroups[,-c(1)] # Remove "X" collumn

# Difference upon chrm gain or loss (Diff & P-value)
Oncogene.Bailey.Diff<- CN.Diff.xRNA.yProt.ThreeGroups[CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% Oncogene.Bailey.list,]
TSG.Bailey.DIff <- CN.Diff.xRNA.yProt.ThreeGroups[CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% TSG.Bailey.list,]


###Protein expression per cell, filtered
## Protein info: 12755 genes
Protein_Info3<-read_csv("Protein_location_info.csv")
Protein_ProID<-read_csv2("Protein_ID_info.csv")
Protein_ProID<-Protein_ProID[,-c(1)]# Remove "X" collumn

#RNA Info. names, chrm location, etc. 
RNA_Info<-read_excel("HGNC_Protein_Info_edited.xlsx", sheet ="HGNC names")
RNA_Info2<-RNA_Info[,c(2,3,11,12,22,14)]

#Get updated Protein info with only 371 cell lines (cell lines also have RNA Data)
# protein data 12755 Proteins
Protein.Expression.filtered<-read.delim2("Protein_Expression_filtered.csv", 
                                         dec=",", header = TRUE, sep=";")
Protein.ID<-data.frame(Protein_ID=
                         colnames(Protein.Expression.filtered[3:length(Protein.Expression.filtered)]))

##Now merge depmap column name info with Protein info. 
## First add Uniprot ID and gene symbol info to depmap data. 
Depmap.Protein.info2<-merge(x= Protein.ID, y= Protein_ProID, 
                            by.x="Protein_ID", by.y="Protein_Id", 
                            sort = TRUE)# 3 collumns, 12755. all Proteins given Uniprot IDs.

###  Aneuploidy data:
## get  2021 cell line arm call data: 
aneuploid<- read.csv("arm_calls_data_bendavid.csv", header=TRUE)



### Copy Number Variance per gene
RNA_Protein_CellLines<- read_delim(file="RNA_Protein_Shared_cells.csv", 
                                   delim=";")


Gene.CNV<- read.csv("CCLE_gene_cn.csv", header=TRUE)
## Filter to get only cells with both RNA and Protein and aneuploidy data
# gene copy number variation, from only cell lines with RNA expression and protein expression  data.
Gene.CNV.filtered<-merge(y= Gene.CNV, x= RNA_Protein_CellLines, 
                                   by.y="X", by.x="Cell_line", 
                                   sort = TRUE)# 370 cells, 27562 genes 



#####  Step 2: Plot OG and TSG difference upon aneuploidy category #####
## Step 2: Plot OG and TSG difference upon aneuploidy category
#   boxplot average OG on gain, TSG on loss (also average genes)
#   barplot average expression difference per gene

## ONCOGENES: 
#Plot boxplot with average oncogene & other gene protein difference upon gain
CN.Diff.xRNA.yProt.ThreeGroups.onco<-CN.Diff.xRNA.yProt.ThreeGroups #make new df with collumn labeled oncogenes, categorical
CN.Diff.xRNA.yProt.ThreeGroups.onco$isOnco<-CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% Oncogene.Bailey.list
CN.Diff.xRNA.yProt.ThreeGroups.onco$isOnco[CN.Diff.xRNA.yProt.ThreeGroups.onco$isOnco==FALSE]<-"Other genes"
CN.Diff.xRNA.yProt.ThreeGroups.onco$isOnco[CN.Diff.xRNA.yProt.ThreeGroups.onco$isOnco==TRUE]<-"Oncogenes"
CN.Diff.xRNA.yProt.ThreeGroups.onco$isOnco<-factor(CN.Diff.xRNA.yProt.ThreeGroups.onco$isOnco, levels=c("Other genes", "Oncogenes"))

ggplot(CN.Diff.xRNA.yProt.ThreeGroups.onco, 
       aes(y=Protein.Diff.Gain, x=isOnco))+
  geom_boxplot()+
  geom_hline(yintercept=0)+
  ylab("Difference upon chrm gain")+
  xlab("Gene type")+
  theme_classic()
#Whiskers = 1.5*IQR
#4x4
#plot.Protein.Diff.Onco.all

# Now plot a few oncogene mean difference
# Difference upon chrm gain or loss (Diff & P-value)
Oncogene.Bailey.Diff<- CN.Diff.xRNA.yProt.ThreeGroups[CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% Oncogene.Bailey.list,]
Oncogene.Bailey.Diff<- Oncogene.Bailey.Diff[order(Oncogene.Bailey.Diff$Protein.Diff.Gain),]
write.csv(Oncogene.Bailey.Diff, file="Oncogene.ExpressionDifference.Gain.Loss.csv")
Oncogene.Bailey.Diff.subset<-Oncogene.Bailey.Diff[Oncogene.Bailey.Diff$RNA_Name %in% c("CTNNB1", "ERBB4", "ERBB2", "PPP2R1A", "MTOR", "MYC", 
                                                                                       "CDK4",  "RHOA", "MET", "PTPN11", "BRAF", "RHOB", 
                                                                                       "MAP2K1", "ERBB3", "KRAS", 
                                                                                        "EGFR", "MAPK1", "HRAS","SMAD4", "RXRA", "DACH1"),]
Oncogene.Bailey.Diff.subset<-Oncogene.Bailey.Diff.subset[-12,]
Oncogene.Bailey.Diff.subset$RNA_Name<-factor(Oncogene.Bailey.Diff.subset$RNA_Name, levels=Oncogene.Bailey.Diff.subset$RNA_Name)

ggplot(Oncogene.Bailey.Diff.subset, 
       aes(y=Protein.Diff.Gain, x=RNA_Name))+
  geom_bar(stat="identity")+
  theme_classic()+
  geom_hline(yintercept=0)+
  geom_hline(yintercept=mean(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain,na.rm=TRUE), color="black", size =1)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"))+
  coord_cartesian(ylim=c(-.75,0.75))+
  ylab("Difference upon chrm gain")+
  xlab("Gene name")
# plot.barplot.Onco.Difference.Mean_v2
# 6x4


### Tumor Suppressor genes (TSG): 
#Plot boxplot with average TSG & other gene protein difference upon Loss
CN.Diff.xRNA.yProt.ThreeGroups.TSG<-CN.Diff.xRNA.yProt.ThreeGroups #make new df with collumn labeled oncogenes, categorical
CN.Diff.xRNA.yProt.ThreeGroups.TSG$isTSG<-CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in%TSG.Bailey.list
CN.Diff.xRNA.yProt.ThreeGroups.TSG$isTSG[CN.Diff.xRNA.yProt.ThreeGroups.TSG$isTSG==FALSE]<-"Other genes"
CN.Diff.xRNA.yProt.ThreeGroups.TSG$isTSG[CN.Diff.xRNA.yProt.ThreeGroups.TSG$isTSG==TRUE]<-"Oncogenes"
CN.Diff.xRNA.yProt.ThreeGroups.TSG$isTSG<-factor(CN.Diff.xRNA.yProt.ThreeGroups.TSG$isTSG, levels=c("Other genes", "Oncogenes"))

ggplot(CN.Diff.xRNA.yProt.ThreeGroups.TSG, 
       aes(y=Protein.Diff.Loss, x=isTSG))+
  geom_boxplot()+
  geom_hline(yintercept=0)+
  ylab("Difference upon chrm loss")+
  xlab("Gene type")+
  theme_classic()
#Whiskers = 1.5*IQR
#4x4
# plot.Protein.Diff.TSG.all


# Now plot a few TSG mean difference
# Difference upon chrm gain or loss (Diff & P-value)
TSG.Bailey.DIff <- CN.Diff.xRNA.yProt.ThreeGroups[CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name %in% TSG.Bailey.list,]
TSG.Bailey.DIff<- TSG.Bailey.DIff[order(TSG.Bailey.DIff$Protein.Diff.Loss),]
#write.csv(TSG.Bailey.DIff, file="TSG.ExpressionDifference.Gain.Loss.csv")
TSG.Bailey.DIff.subset<-TSG.Bailey.DIff[TSG.Bailey.DIff$RNA_Name %in% c("SMAD4", "NF2", "STK11", "PTEN", "CASP8", "ATM", 
                                                                        "RASA1", "NOTCH1", "NF1", "APC", "EP300", "CREBBP", 
                                                                        "CDH1", "SOX9", "BRCA1", "KANSL1", 
                                                                        "CDKN1A", "CDK12","TP53", "TGFBR2", "CTNND1", "CDKN2A"),]
TSG.Bailey.DIff.subset<-TSG.Bailey.DIff.subset[-23,]
TSG.Bailey.DIff.subset$RNA_Name<-factor(TSG.Bailey.DIff.subset$RNA_Name, levels=TSG.Bailey.DIff.subset$RNA_Name)

ggplot(TSG.Bailey.DIff.subset, 
       aes(y=Protein.Diff.Loss, x=RNA_Name))+
  geom_bar(stat="identity")+
  theme_classic()+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(-.75,0.75))+
  geom_hline(yintercept=mean(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss,na.rm=TRUE), color="black", size =1)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"))+
  ylab("Difference upon chrm loss")+
  xlab("Gene name")
# plot.barplot.TSG.Difference.Mean_v2
# 6x4




#####  Step 3: Plot OG and TSG difference upon gene Copy number variation ####
# Set working directory to where you want to safe figures to. 
setwd(FigureResultsDirectory)


## Step 3: Plot OG and TSG difference upon Copy number variation
#   scatterplot of CNV/ploidy vs expression per gene
#   categorical: not amplified, one-gain, multi-gain

# INSTRUCTIONS: run each gene you are interested in through below set of analtsis/plotting. 
# note that some genes have multiple protein forms, and you'll have to make minor alterations 
# to the code to correctly select the protein version you want to analyse. 

TestGene<-"ERBB2" ### !!! Change the name of the test gene to the gene you want to test!!! 
# MAPK1, MTOR, MYC, SMAD4, KANSL1, TP53, CTNND1, CDKN2A, CDK12
# MYC->  buffered, CNV corr
# EGFR -> scaling , CNV corr
# KRAS --> buffered , CNV corr
# ERBB2 --> anti scaling , CNV corr
# CDK4 --> buffered, CNV corr
##Get  cell ID and CNV/cell for a specific gene 
CNV.col.testgene<-Gene.CNV.filtered[,sub("\\..*","",colnames(Gene.CNV.filtered))==TestGene] 
#find collumn for test gene
CNV.testGene<- data.frame(Cell_Line=Gene.CNV.filtered$Cell_line,
                        DNA.CopyNum= CNV.col.testgene )

##Get Cell ID and Protein expression/cell for specific gene
TestGene.info<-subset(Protein_ProID, Protein_ProID$Gene_Symbol==TestGene)#get protein name for gene

Protein.testGene<- data.frame( Cell_Line=Protein.Expression.filtered$Cell_line,
  Protein.expression= Protein.Expression.filtered[,colnames(Protein.Expression.filtered)==TestGene.info$Protein_Id[1]] )

## Combine CNV and expression data for specific gene
TestGene.Prot.CNV <- merge(x=Protein.testGene, y=CNV.testGene, 
                         by.x="Cell_Line", by.y="Cell_Line")

## Combine CNV, Protein expression data with aneuploidy data
TestGene.info2<-subset(Protein_Info3, Protein_Info3$`Approved Symbol`==TestGene)#get protein name for gene

aneuploid2<-subset(aneuploid, aneuploid$chrom==TestGene.info2$Chromosome) #Get aneu data for correct chrm
aneuploid2<-subset(aneuploid2, aneuploid2$arm==TestGene.info2$arm) #Get aneu data for correct chrm arm
aneuploid2<-subset(aneuploid2, aneuploid2$DepMap_ID %in% RNA_Protein_CellLines$Cell_line) #Get aneu data for correct chrm arm


TestGene.Prot.CNV2<-merge(x=TestGene.Prot.CNV, y=aneuploid2, 
                         by.x="Cell_Line", by.y="DepMap_ID")
TestGene.Prot.CNV2$arm_call[TestGene.Prot.CNV2$arm_call=="-1"]<-"Loss"
TestGene.Prot.CNV2$arm_call[TestGene.Prot.CNV2$arm_call=="0"]<-"Neutral"
TestGene.Prot.CNV2$arm_call[TestGene.Prot.CNV2$arm_call=="1"]<-"Gain"
TestGene.Prot.CNV2$arm_call<-factor(TestGene.Prot.CNV2$arm_call, levels = c("Loss", "Neutral", "Gain"))

#TestGene.Prot.CNV3<- subset(TestGene.Prot.CNV2, DNA.CopyNum<1.75)

##Plot x=CNV vs. y=expression scatterplot
pdf(file = paste0("plot.Protein.scatter.expression.GeneCN.", TestGene, ".pdf"),
    width = 4, 
    height = 4)
ggplot(TestGene.Prot.CNV3, aes(x=DNA.CopyNum, y=Protein.expression))+#, color=factor(arm_call)))+
  geom_point(size=2)+  
  #scale_color_manual(values=c("Cyan", "Grey50", "Red"), name="Chromosome arm CN")+
  xlab("Relative gene copy number")+
  ylab("Protein expression")+
  geom_smooth(method="lm", color="Red")+
  #coord_cartesian(xlim=c(0,8), ylim=c(-3.5,6))+
  theme_classic()+
  ggtitle(paste0(TestGene, " expression vs relative gene copy number"))
# 4x4
# plot.Protein.scatter.expression.GeneCN.SMAD4
dev.off()

# pearson correlation coefficient and p-value
cor.test(TestGene.Prot.CNV2$Protein.expression, TestGene.Prot.CNV2$DNA.CopyNum, method="pearson")



## Plot y=expresion vs. x=CNV (category) as boxplot
# CCLE gene copy number are log2(1+n) set relative to basal ploidy of the cell
# categories: Loss (<0.9), neutral(0.9-1.1), gain(1.1-1.75), amplified(1.75+)
TestGene.Prot.CNV$DNAcat<-cut(TestGene.Prot.CNV$DNA.CopyNum, 
                              breaks=c(0, 0.9, 1.1, 1.75, Inf), 
                              labels=c("Loss", "Neutral", "Gain", "Amplified"))
TestGene.Prot.CNV$DNAcat<-factor(TestGene.Prot.CNV$DNAcat, levels=c("Loss", "Neutral", "Gain", "Amplified"))
n<-subset(TestGene.Prot.CNV, DNAcat=="Neutral")
l<-subset(TestGene.Prot.CNV, DNAcat=="Loss")
g<-subset(TestGene.Prot.CNV, DNAcat=="Gain")
a<-subset(TestGene.Prot.CNV, DNAcat=="Amplified")
t.test(n$Protein.expression,l$Protein.expression) #Loss vs neutral
t.test(n$Protein.expression,g$Protein.expression) #gain vs neutral
t.test(n$Protein.expression,a$Protein.expression) #amp vs neutral


pdf(file = paste0("plot.Protein.box.expression.GeneCN.", TestGene, ".pdf"),
    width = 5, 
    height = 4)
ggplot(TestGene.Prot.CNV, 
       aes(y=Protein.expression, x=DNAcat, fill=DNAcat))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  scale_fill_manual(values=c("dodgerblue3", "grey80", "red", "dark red"), name="Gene CN")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(-6,6))+
  ylab("Protein expression")+
  xlab("Gene copy number category")+
  theme_classic()+
  ggtitle(paste0(TestGene," expression vs gene copy number"))
# 5x4
# plot.Protein.box.expression.GeneCN.MAPK1
dev.off()


## To check that arm call and expression level are as they were listed. 
pdf(file = paste0("plot.Protein.box.expression.ChrmCN.", TestGene, ".pdf"),
    width = 5, 
    height = 4)
ggplot(TestGene.Prot.CNV2, 
       aes(y=Protein.expression, x=arm_call, fill=arm_call))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  scale_fill_manual(values=c("dodgerblue3", "grey80", "red", "dark red"), name="Chromosome CN")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(-10,10))+
  ylab("Protein expression")+
  xlab("Chromosome arm relative copy number")+
  theme_classic()+
  ggtitle(paste0(TestGene," expression vs chromosome copy number"))
# 5x4
# plot.Protein.box.expression.ChrmCN.CTNND1
#skyblue1, dodgerblue3
dev.off()


#####  Step 4: Find overlap between gene copy number increase/decrease and chromosome copy number gain/loss ####
#   Categorical gene copy number: lost, neutral and gain
#   find what percent of genes in chrm gain are not gained at gene level
# DNA copy number is set relative to cell mean ploidy
# categorical cutoffs from above: 
# cut(TestGene.Prot.CNV$DNA.CopyNum, 
#                              breaks=c(0, 0.9, 1.1, 1.75, Inf), 
#                              labels=c("Gene Loss", "Gene Neutral", "Gene Gain"))

# Substep 1: Get data
Gene.CNV.filtered # see above, gene copy number per cell line data


# Substep 2: Cut data into categories based on gene copy number 
Gene.CNV.Category<- data.frame(Cell_Line= Gene.CNV.filtered$Cell_line, 
                               X=Gene.CNV.filtered$X)
Gene.CNV.Category[3:length(Gene.CNV.filtered)]<- NA
colnames(Gene.CNV.Category)<-colnames(Gene.CNV.filtered)

for (i in 3:length(Gene.CNV.filtered)){
  Gene.CNV.Category[i]<-cut(Gene.CNV.filtered[,i],
                       breaks=c(0, 0.9, 1.1, Inf),
                       labels=c("Gene_Loss", "Gene_Neutral", "Gene_Gain"))
}

# Substep 3: Find chrm copy number for each gene per cell
# Depmap.Protein.info5
# Protein.Expression.filtered
# format collumn headers as "Gene_Symbol..EnterezID." so that it matches Gene.CNV.Category headers
# Gene info has GENE info, not Protein info. 
Chrm.CNV.Category<- data.frame(DepMap_ID= unique(aneuploid$DepMap_ID))

for (i in 1:length(Depmap.Protein.info5$Gene_Symbol)){
  gene<-Depmap.Protein.info5[i,]
  #use If statement to check that Protein Info has Protein of interest. 
    TestArm <- gene$arm #Find test protein chromosome arm
    TestChrm<- as.numeric(gene$Chromosome) #Find test protein chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test protein location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm & arm==TestArm) 
    
    Chrm.CNV.Category<- merge(Chrm.CNV.Category, testchrm.percell[,c(2,7)], #only need cell line and arm call
                                 by.x="DepMap_ID", by.y= "DepMap_ID")
    colnames(Chrm.CNV.Category)[length(Chrm.CNV.Category)]<-paste0(gene$Gene_Symbol, "..", gene$`Entrez Gene ID`, ".") #label arm_call collumn with protein ID
}

Chrm.CNV.Category2<- Chrm.CNV.Category[,!duplicated(names(Chrm.CNV.Category))] #remove duplicated genes
# may have duplicated genes because multiple proteins can be made from one gene. 
# 9124 collumns, 9123 unique genes

## Substep 4: Find what percent of genes have different Gene category and Chromosome category
# make new data frame with info: for gene & cell: 
# chrm gain-    number gene loss, neutral and gain 
# chrm neutral- number gene loss, neutral and gain 
# chrm loss-    number gene loss, neutral and gain 
# look at sum(chrm gain/gene gain, chrm N/gene N, Chrn loss/gene loss) / all --> percent same

#Chrm.CNV.Category
#Gene.CNV.Category

# get list of cells shared by both Gene CNV data and Chromosome CNV data
# get list of genes shared by both Gene CNV data and Chromosome CNV data
ListofCells_Gene_Chrm<- Chrm.CNV.Category2$DepMap_ID[Gene.CNV.Category$Cell_line] #367 cell lines

#subset above Chrmosome and Gene copy number info for shared cells and genes
Chrm.CNV.Category_filtered<-Chrm.CNV.Category2[ListofCells_Gene_Chrm,]
Gene.CNV.Category_filtered<-Gene.CNV.Category[ListofCells_Gene_Chrm,]

colnames(Chrm.CNV.Category_filtered)<- sub("\\.\\..*", "", colnames(Chrm.CNV.Category_filtered))
colnames(Gene.CNV.Category_filtered)<- sub("\\.\\..*", "", colnames(Gene.CNV.Category_filtered))
ListofGenes<- colnames(Gene.CNV.Category_filtered)[colnames(Gene.CNV.Category_filtered) %in% colnames(Chrm.CNV.Category_filtered)]
ListofGenes<-append(ListofGenes, c("Cell_line", "DepMap_ID"))

# 9078 genes shared between databases

Chrm.CNV.Category_filtered<-subset(Chrm.CNV.Category_filtered, select = ListofGenes)
Gene.CNV.Category_filtered<-Gene.CNV.Category_filtered[,as.character(ListofGenes)]


#make empty vectors
ChrmGain.GeneGain<- 0
ChrmGain.GeneNeutral<- 0
ChrmGain.GeneLoss<- 0
ChrmNeutral.GeneGain<- 0
ChrmNeutral.GeneNeutral<- 0
ChrmNeutral.GeneLoss<- 0
ChrmLoss.GeneGain<- 0
ChrmLoss.GeneNeutral<- 0
ChrmLoss.GeneLoss<- 0


for (i in 1:length(ListofGenes)){
  Gene.CNV1<-subset(Gene.CNV.Category_filtered, select=c(ListofGenes[i], "Cell_line")) #get Gene CNV info per gene
  Chrm.CNV1<-subset(Chrm.CNV.Category_filtered, select=c(ListofGenes[i], "DepMap_ID")) #get Chrm CNV info per gene
  Chrm.Gene.CNV1<- merge(Chrm.CNV1, Gene.CNV1, by.x="DepMap_ID", by.y="Cell_line")
  
  ChrmGain.GeneGain<- ChrmGain.GeneGain + sum(Chrm.Gene.CNV1[,2]==1 & Chrm.Gene.CNV1[,3]=="Gene_Gain")
  ChrmGain.GeneNeutral<- ChrmGain.GeneNeutral + sum(Chrm.Gene.CNV1[,2]==1 & Chrm.Gene.CNV1[,3]=="Gene_Neutral")
  ChrmGain.GeneLoss<- ChrmGain.GeneLoss + sum(Chrm.Gene.CNV1[,2]==1 & Chrm.Gene.CNV1[,3]=="Gene_Loss")
  ChrmNeutral.GeneGain<- ChrmNeutral.GeneGain + sum(Chrm.Gene.CNV1[,2]==0 & Chrm.Gene.CNV1[,3]=="Gene_Gain")
  ChrmNeutral.GeneNeutral<- ChrmNeutral.GeneNeutral + sum(Chrm.Gene.CNV1[,2]==0 & Chrm.Gene.CNV1[,3]=="Gene_Neutral")
  ChrmNeutral.GeneLoss<- ChrmNeutral.GeneLoss + sum(Chrm.Gene.CNV1[,2]==0 & Chrm.Gene.CNV1[,3]=="Gene_Loss")
  ChrmLoss.GeneGain<- ChrmLoss.GeneGain + sum(Chrm.Gene.CNV1[,2]== -1 & Chrm.Gene.CNV1[,3]=="Gene_Gain")
  ChrmLoss.GeneNeutral<- ChrmLoss.GeneNeutral + sum(Chrm.Gene.CNV1[,2]== -1 & Chrm.Gene.CNV1[,3]=="Gene_Neutral")
  ChrmLoss.GeneLoss<- ChrmLoss.GeneLoss + sum(Chrm.Gene.CNV1[,2]== -1 & Chrm.Gene.CNV1[,3]=="Gene_Loss")
}

Chrm.Gene.CopyNumber.Overlap<-data.frame(Chrm_Gain<-c(ChrmGain.GeneGain, ChrmGain.GeneNeutral, ChrmGain.GeneLoss),
                                         Chrm_Neutral<-c(ChrmNeutral.GeneGain, ChrmNeutral.GeneNeutral, ChrmNeutral.GeneLoss), 
                                         Chrm_Loss<-c(ChrmLoss.GeneGain, ChrmLoss.GeneNeutral, ChrmLoss.GeneLoss))
rownames(Chrm.Gene.CopyNumber.Overlap)<-c("Gene_Gain", "Gene_Neutral", "Gene_Loss")
colnames(Chrm.Gene.CopyNumber.Overlap)<-c("Chrm_Gain", "Chrm_Neutral", "Chrm_Loss")
Chrm.Gene.CopyNumber.Overlap

Percent_overlap<-100*sum(ChrmGain.GeneGain, ChrmNeutral.GeneNeutral, ChrmLoss.GeneLoss)/ sum(Chrm.Gene.CopyNumber.Overlap)
Percent_overlap #70.34% of genes are in same category as chromosome gain/loss
Percent_Anti<-100*sum(ChrmGain.GeneLoss, ChrmLoss.GeneGain)/ sum(Chrm.Gene.CopyNumber.Overlap)
Percent_Anti # 2.96% of genes are in opposite category as chromosome gain/loss, mostly gained upon chrm loss

## cutoffs 0.9 to 1.1-We examined whether or not the extensive protein buffering we find in our data is caused by gene mutations, gene “flipping”, low expression artifacts, unreliable protein measurements or mass spec specific measurement effects. -> 70% gene/chrm correct, and 2.96% opposite of expected
Chrm.Gene.CopyNumber.Overlap_0.9_1.1<-Chrm.Gene.CopyNumber.Overlap
Percent_ChrmGain.GeneLoss<-100*sum(ChrmGain.GeneLoss)/ sum(Chrm.Gene.CopyNumber.Overlap$Chrm_Gain )
## opposite of expected: 2.40% of genes on gained chromosomes are lost at gene level
Percent_ChrmLoss.GeneGain<-100*sum(ChrmLoss.GeneGain)/ sum(Chrm.Gene.CopyNumber.Overlap$Chrm_Loss )
## opposite of expected: 10.38% of genes on lost chromosomes are gained at gene level



## Other cutoff attempts: 
## cutoffs 0.8 to 1.2--> 71.6% gene/chrm correct, and 1.66% opposite of expected
#Chrm.Gene.CopyNumber.Overlap_0.8_1.2<-Chrm.Gene.CopyNumber.Overlap
Percent_ChrmGain.GeneLoss<-100*sum(ChrmGain.GeneLoss)/ sum(Chrm.Gene.CopyNumber.Overlap$Chrm_Gain )
## opposite of expected: 1.59% of genes on gained chromosomes are lost at gene level
Percent_ChrmLoss.GeneGain<-100*sum(ChrmLoss.GeneGain)/ sum(Chrm.Gene.CopyNumber.Overlap$Chrm_Loss )
## opposite of expected: 5.68% of genes on lost chromosomes are gained at gene level

## cutoffs 0.9 to 1.2--> 75.4% gene/chrm correct, and 1.77% opposite of expected 
#  Chrm.Gene.CopyNumber.Overlap_0.9_1.2<-Chrm.Gene.CopyNumber.Overlap
Percent_ChrmGain.GeneLoss<-100*sum(ChrmGain.GeneLoss)/ sum(Chrm.Gene.CopyNumber.Overlap$Chrm_Gain )
## opposite of expected: 2.40% of genes on gained chromosomes are lost at gene level
Percent_ChrmLoss.GeneGain<-100*sum(ChrmLoss.GeneGain)/ sum(Chrm.Gene.CopyNumber.Overlap$Chrm_Loss )
## opposite of expected: 5.68% of genes on lost chromosomes are gained at gene level


#####  Step 5: Find list of genes/cell that has flipped GCN and CCN #####
#   correlate differnce in expression upon chrm and gene gain v neutral (& loss vs neutral)

Chrm.CNV.Category_filtered<-subset(Chrm.CNV.Category_filtered, select = ListofGenes)
Gene.CNV.Category_filtered<-Gene.CNV.Category_filtered[,as.character(ListofGenes)]

List.ChrmGain.GeneLoss<-data.frame(Cell_line=as.character,
                                   Gene_ID= as.character)
List.ChrmLoss.GeneGain<-data.frame(Cell_line=as.character,
                                   Gene_ID= as.character)

for (i in 1:(length(ListofGenes)-2)){
  Gene.CNV1<-subset(Gene.CNV.Category_filtered, select=c(ListofGenes[i], "Cell_line")) #get Gene CNV info per gene
  Chrm.CNV1<-subset(Chrm.CNV.Category_filtered, select=c(ListofGenes[i], "DepMap_ID")) #get Chrm CNV info per gene
  Chrm.Gene.CNV1<- merge(Chrm.CNV1, Gene.CNV1, by.x="DepMap_ID", by.y="Cell_line")
  
  if (sum(Chrm.Gene.CNV1[,2]== 1 & Chrm.Gene.CNV1[,3]=="Gene_Loss")!=0){
  List.ChrmGain.GeneLoss<- rbind(List.ChrmGain.GeneLoss, data.frame(Cell_line= Chrm.Gene.CNV1[(Chrm.Gene.CNV1[,2]== 1 & Chrm.Gene.CNV1[,3]=="Gene_Loss"),1], 
                                                                    Gene_ID= ListofGenes[i]) )
  }
  if (sum(Chrm.Gene.CNV1[,2]== -1 & Chrm.Gene.CNV1[,3]=="Gene_Gain")!=0){
  List.ChrmLoss.GeneGain<- rbind(List.ChrmLoss.GeneGain, data.frame(Cell_line= Chrm.Gene.CNV1[(Chrm.Gene.CNV1[,2]== -1 & Chrm.Gene.CNV1[,3]=="Gene_Gain"),1], 
                                                                    Gene_ID= ListofGenes[i]) )
  }
}

List.GeneCN.flipped.ChrmCN<-rbind(List.ChrmLoss.GeneGain, List.ChrmGain.GeneLoss)
# 33085 genes gained upon chrm loss or lost upon chrm gain
write.csv(List.GeneCN.flipped.ChrmCN, file="List.Flipped.GeneCN.ChrmCN.csv")



#####  Step 6: Does gene amplification percent (of dataset) affect buffering? DID NOT END UP USING THIS ####

# did not end up using this in the final paper. 
# But percent of cells in dataset with gene amplification does predict buffering upon chrm loss. 
# make dataset about percent of cells with gene amplification that I used in buffering factor analysis. 

## Find percent of cell with amplified (>1.75 relative CN) gene Copy number, per gene. 
## See if this significantly affects Buffering ROC AUC. 
GeneAmplification<-data.frame(Gene_ID=character(), 
                              Gene_Name=character(),
                              AmplificationRatio=numeric(), 
                              stringsAsFactors = TRUE)

for (i in 4:length(Gene.CNV.filtered)){#first three columns are cell line data
  # select collumn, each collumn = a gene, each row a cell
  # CCLE gene copy number are log2(1+n) set relative to basal ploidy of the cell
  # set cutoff for amplified at 1.75
  # would need 6+ copies in diploid, or 11+ copies in triploid, 16+ in tetra
  percentAmp<-sum(Gene.CNV.filtered[,i]>1.75)/370 #370 cell lines in this data. 
  GeneAmplification<-rbind(GeneAmplification, 
                           data.frame(Gene_ID=colnames(Gene.CNV.filtered[i]),
                                      Gene_Name=sub("\\..*","",colnames(Gene.CNV.filtered[i])),
                                      AmplificationRatio= percentAmp))
  
}
write.csv(GeneAmplification, 
          file="AmplificationRatio_perGene1.75.csv")


CN.Diff.RNA.Prot.3Groups.Amp<- merge(x=CN.Diff.xRNA.yProt.ThreeGroups, y=GeneAmplification, 
                                     by.x="RNA_ID", by.y="Gene_ID")
Amp.LossDiff.Corr<-cor.test(CN.Diff.RNA.Prot.3Groups.Amp$Protein.Diff.Loss, CN.Diff.RNA.Prot.3Groups.Amp$AmplificationRatio, method="pearson")
# p= 8E-11
# corr= 0.067
Amp.GainDiff.Corr<-cor.test(CN.Diff.RNA.Prot.3Groups.Amp$Protein.Diff.Gain, CN.Diff.RNA.Prot.3Groups.Amp$AmplificationRatio, method="pearson")
#p= 1E-2
#corr= 0.026


# Corr of protein diff upon LOSS with amplification: 
ggplot(CN.Diff.RNA.Prot.3Groups.Amp, 
       aes(y=Protein.Diff.Loss, x=100*AmplificationRatio))+
  geom_point(size=2)+
  theme_classic()+
  coord_cartesian(ylim=c(-2,2), xlim=c(0,20))+
  geom_smooth(method="lm", color="Red")+
  ylab("Protein difference upon chrm loss")+
  xlab("Percent of cell lines with amplifications")
#4x4
# p= 1E-10
# corr= 0.066
# plot.Protein.LossDiff.Amp.sig
CN.Diff.RNA.Prot.3Groups.Amp$Amplification.Group<-cut(CN.Diff.RNA.Prot.3Groups.Amp$AmplificationRatio, breaks=c(-Inf, 0.01,0.05, Inf), 
                                                      labels=c("0-1%", "1-5%", ">5%"))
amp0.loss<-subset(CN.Diff.RNA.Prot.3Groups.Amp, Amplification.Group=="0-1%")
amp1.loss<-subset(CN.Diff.RNA.Prot.3Groups.Amp, Amplification.Group=="1-5%")
amp5.loss<-subset(CN.Diff.RNA.Prot.3Groups.Amp, Amplification.Group==">5%")
t.test(amp0.loss$Protein.Diff.Gain, amp1.loss$Protein.Diff.Gain)

# Groups of protein diff upon LOSS with amplification: 
ggplot(CN.Diff.RNA.Prot.3Groups.Amp, 
       aes(y=Protein.Diff.Loss, x=Amplification.Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  coord_cartesian(ylim=c(-1,1))+
  #geom_smooth(method="lm", color="Red")+
  ylab("Protein difference upon chrm loss")+
  xlab("Percent of cell lines with amplifications")
#4x4
# plot.boxplot.Protein.loss.Diff.percentAmplifications
# 0 to 1: 7E-05;  0 to 5: 4E-03

## Corr of protein difference upon gain with amplifications
ggplot(CN.Diff.RNA.Prot.3Groups.Amp, 
       aes(y=Protein.Diff.Gain, x=100*AmplificationRatio))+
  geom_point(size=2)+
  theme_classic()+
  coord_cartesian(ylim=c(-2,2), xlim=c(0,20))+
  geom_smooth(method="lm", color="Red")+
  ylab("Protein difference upon chrm gain")+
  xlab("Percent of cell lines with amplifications")
#4x4
#p= 1E-2
#corr= 0.027
# plot.Protein.GainDiff.Amp.sig
# top gene= PRSS2
# MYC is at 11, near 0 difference

# Groups of protein diff upon Gain with amplification: 
ggplot(CN.Diff.RNA.Prot.3Groups.Amp, 
       aes(y=Protein.Diff.Gain, x=Amplification.Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  coord_cartesian(ylim=c(-1,1))+
  #geom_smooth(method="lm", color="Red")+
  ylab("Protein difference upon chrm gain")+
  xlab("Percent of cell lines with amplifications")
#4x4
# plot.boxplot.Protein.gain.Diff.percentAmplifications
# 0 to 1: 0.021;  0 to 5: NS
##

### Look at TSG and Oncogenes:
TSG.Bailey.DIff.Amp <- CN.Diff.RNA.Prot.3Groups.Amp[CN.Diff.RNA.Prot.3Groups.Amp$RNA_Name %in% TSG.Bailey.list,]
TSG.Bailey.DIff.Amp<- TSG.Bailey.DIff.Amp[order(TSG.Bailey.DIff.Amp$Protein.Diff.Loss),]

ggplot(TSG.Bailey.DIff.Amp, 
       aes(y=Protein.Diff.Loss, x=100*AmplificationRatio))+
  geom_point()+
  theme_classic()+
  coord_cartesian(ylim=c(-.75,0.75), xlim=c(0,11))+
  geom_smooth(method="lm", color="Red")+
  ylab("Protein difference upon chrm loss")+
  xlab("Percent of cell lines with amplifications")
# plot.TSG.Difference.GeneAmp.NS
# 4x4


Onco.Bailey.DIff.Amp <- CN.Diff.RNA.Prot.3Groups.Amp[CN.Diff.RNA.Prot.3Groups.Amp$RNA_Name %in% Oncogene.Bailey.list,]
Onco.Bailey.DIff.Amp<- Onco.Bailey.DIff.Amp[order(Onco.Bailey.DIff.Amp$Protein.Diff.Gain),]

ggplot(Onco.Bailey.DIff.Amp, 
       aes(y=Protein.Diff.Gain, x=(100*AmplificationRatio)))+
  geom_point()+
  theme_classic()+
  coord_cartesian(ylim=c(-.75,0.75), xlim=c(0,11))+
  geom_smooth(method="lm", color="Red")+
  ylab("Protein difference upon chrm gain")+
  xlab("Percent of cell lines with amplifications")
# plot.Onco.Difference.GeneAmp.NS
# 4x4

