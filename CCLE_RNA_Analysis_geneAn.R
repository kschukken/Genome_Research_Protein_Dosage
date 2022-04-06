##### RNA expression data correlation with cellular aneuploidy score ####
### RNA expression correlate with aneuploidy score (bp per chrm corrected for ploidy)
### 200518
### Updated: 210217- now use gene based aneuploidy scores
### Klaske Schukken
### RNA dropout analysis
### RNA correlation to aneuploidy scores
### Using AVANA 20Q2 CCLE data. 
### downloaded data on 200512

library(ggplot2)
library(tidyverse)
library(readxl)
library("corrr")

setwd()#set working directory 

DataFileLocation<- "Documents" # ! Set to location you want files to be written to


###INFO: 
# ploidy is mean chromosome count
# arm_Score is sum of aneuploid chromosome arms. 
# bp_Score is aneuploid abs(CN of chrm -  round(ploidy)) * length of arm (in bp)
#         this way it gives you 2 if cell has gain of 2 chrm. so #aneuploid bp double
# bp_ploidy_Score is aneuploid bp_Score/ ploidy 
#         so that it accounts for gene dosage imbalance: 
#         ex. tetraploid needs 2x aneuploid chrm for same imbalance. 
# Gene_Score is aneuploid chrm arm * genes on arm
#         only genes with known protein coding genes
#         downloaded data 210217
#         from HGNC. https://www.genenames.org/
# Gene_ploidy_Score is aneuploid chrm arm * genes on arm/ ploidy 
#         so that it accounts for gene dosage imbalance: 
#         ex. tetraploid needs 2x aneuploid chrm for same imbalance.
#         !!! This is the score used as the "Cellular aneuploidy score" in manuscript

#CCLE is RNA seq file expression data. only the protein coding genes
#RSEM log2 transformed 
#19144 genes
#1304 cell lines 


##### Step 1: Get CCLE RNA expression data and merge ####
### from depmap.org
### check data frame correct, and rename cell lines as "Cell_line"
CCLE_RNA<-read.delim2("CCLE_expression.csv", 
                             dec=".", header = TRUE, sep="\t")

names(CCLE_RNA)[names(CCLE_RNA) == "X"] <- "Cell_line"

length(CCLE_RNA[,1])
CCLE_RNA[5,19]

# Generated in Aneuploid_Score.R
Score<- read.csv("Score2.aneuploid.cell_line.csv", header=TRUE)

##Gene groups: 
rRNAMet<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=9) 
CORUM.all<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=2) 
ERMem<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=8) 


##### Now add aneuploidy scores per cell lines
### get Score data from Aneuploidy_Score.R
CCLE_RNA_Score<-merge(x= CCLE_RNA, y= Score, 
                             by.x="Cell_line", by.y="DepMap_ID", 
                             sort = TRUE)
length(CCLE_RNA_Score[1,])

write.csv(CCLE_RNA_Score,
          "CCLE_RNA_Score2_GeneAneuploid.csv", row.names = TRUE)

#setwd()
#CCLE_RNA_Score<-read.delim2("CCLE_RNA_Score2_GeneAneuploid.csv", 
#                            dec=".", header = TRUE, sep=",")


#step 1.5- make some plots/overview of data 
plot.ploidy<- ggplot(CCLE_RNA_Score,aes(ploidy)) +
  geom_bar()

plot.bp.Score<- ggplot(CCLE_RNA_Score,aes(bp_Score)) +
  geom_histogram(binwidth = 100000000)

plot.bp.ploidy.Score<- ggplot(CCLE_RNA_Score,aes(bp_ploidy_Score)) +
  geom_histogram(binwidth = 25000000)

plot.ploidy
plot.bp.Score
plot.bp.ploidy.Score



##### Step 2: Corrolate RNA expression with aneuploidy score ####
### Using corrr, pearson corrolation 

CCLE_RNA_corr<- CCLE_RNA_Score[sapply(CCLE_RNA_Score, 
                                              function(x) !is.factor(x))] %>% 
  correlate() %>% 
  focus(arm_Score, bp_Score, bp_ploidy_Score, Gene_Score, Gene_ploidy_Score)


CCLE_RNA_corr <- CCLE_RNA_corr[order(-CCLE_RNA_corr$Gene_ploidy_Score),]

#write.csv(CCLE_RNA_corr,
#          "CCLE_RNA_correlation_GeneAneuploidy.csv", row.names = TRUE)

x<-CCLE_RNA_corr$rowname
CCLE_RNA_corr<- CCLE_RNA_corr %>% separate(rowname, c("Name", "ID"))
CCLE_RNA_corr$rowname<-x

write.csv(CCLE_RNA_corr,
          "CCLE_RNA_Gene_Names.csv", row.names = TRUE)


#setwd()#set working directory 
#CCLE_RNA_corr<-read.delim2("CCLE_RNA_correlation_GeneAneuploidy.csv", 
#                            dec=".", header = TRUE, sep=",")
CCLE_RNA_corr<-subset(CCLE_RNA_corr, select= -c(X))#re-run data if you want arm score


##### Step 3: Order genes by p-value. find top pos and negative correlated #### 
### then save data in that excel sheet. 
### then plot overall data. 

length(CCLE_RNA_corr)# 6 collumns
CCLE.pos.cor<-CCLE_RNA_corr[1:10,] #highest  corrolation
CCLE.neg.cor<-CCLE_RNA_corr[19105:19115,] #highest negative corrolation
CCLE.neg.cor<- CCLE.neg.cor[order(CCLE.neg.cor$bp_ploidy_Score),]

CCLE.neg.cor
CCLE.pos.cor


#ordered by absolute bp ploidy score. so pos & neg mixed. 
#because "lowest" scores, after normal sorting, are NA. 

CCLE_GenePloidy_corr<- CCLE_RNA_corr[order(-abs(CCLE_RNA_corr$Gene_ploidy_Score)),]
head(CCLE_GenePloidy_corr)

#####         Plot correlation data ####
#plot scatter plot and ploidy. 
#Bp ploidy vs arm score
plot.cor.scatter<- ggplot(RNA.corr.pvalue2[2:19145,], #remove ploidy
                          aes(x=bp_ploidy_Score, y=bp_Score))+
  geom_point()+
  ylab("Correlation: RNA expression & arm score") +
  xlab("Correlation: RNA expression & bp ploidy score")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()

plot.cor.scatter# Removed 30 rows with missing values

#Bp ploidy vs gene ploidy score
plot.gene<- ggplot(CCLE_bp_corr[2:19145,], 
                   aes(x=Gene_ploidy_Score, y=Gene_Score))+
  geom_point(color="black")+
  ylab("Correlation: RNA expression & gene score") +
  xlab("Correlation: RNA expression & gene ploidy score")+
  theme_classic()

plot.gene

#Bp ploidy vs gene ploidy score
plot.gene<- ggplot(RNA.corr.pvalue2[2:19145,], 
                   aes(x=Gene_ploidy_Score, y=bp_ploidy_Score))+
  geom_point(color="black")+
  ylab("Correlation: RNA expression & base-pair ploidy score") +
  xlab("Correlation: RNA expression & gene ploidy score")+
  theme_classic()

plot.gene


##### Step 4: Plot top 2 pos/neg correlated genes ####
### Step 4: calculate top p-values and plot genes of interest
#Top 10 total correlations
#Note that the top 20 are definately significant. with 20th gene having
#a p-value of 1.1e-29. 
head(CCLE_GenePloidy_corr)

RNA.corr.test<-cor.test(CCLE_RNA_Score$Gene_ploidy_Score, 
                       CCLE_RNA_Score$LAMA5..3911., 
                       method = "pearson")
RNA.corr.test
RNA.corr.pvalue<- RNA.corr.test$p.value 

# RPL3..6122. 
# PPP3CC..5533. 
# RPL34..6164.
# KIAA1522 p= 8.7e-37. r=0.417

ggplot(CCLE_RNA_Score, 
                    aes(x=Gene_ploidy_Score, y=RPL34..6164.))+
  geom_point(color="dodgerblue3", size=2)+
  ylab("RNA expresion: RPL34") +
  xlab("Aneuploidy Score: gene ploidy score")+
  geom_smooth(color="black", method="lm")+
  #geom_text(x=7000, y=2, label=paste0("p= 7.7E-36"))+
  #geom_text(x=7000, y=3, label= paste0("cor= 0.39"))+
  theme_classic()
# plot.RNA.aneuscore.Express.RPL34.Blue
#4x4
# red for pos, dodgerblue3 for neg corr

##### Side Step 1: Look at whole Genome doubling vs. no doubling ####
## Did not end up using this. 
#head(Gene_ploidy_Score)
#Look at whole Genome doubling vs. no doubling t-test:
#CCLE_RNA_Score$WGD<-CCLE_RNA_Score$ploidy>2.25

#ggplot(CCLE_RNA_Score, aes(WGD, TSPAN15..23555.)) +
#  geom_boxplot()+
#  ggtitle("TSPAN15 RNA: t-test p-value= 6E-15, Ploidy>2.25 is WGD+")+
#  theme_classic()
#t.test(TSPAN15..23555. ~ WGD, data = CCLE_RNA_Score)

#ggplot(CCLE_RNA_Score, aes(WGD, KIAA1522..57648.)) +
#  geom_boxplot()+
#  ggtitle("KIAA1522 RNA: t-test p-value= 1E-15, Ploidy>2.25 is WGD+")+
#  theme_classic()
#t.test(KIAA1522..57648. ~ WGD, data = CCLE_RNA_Score)

#ggplot(CCLE_RNA_Score, aes(WGD, CTSV..1515.)) +
#  geom_boxplot()+
#  ggtitle("CTSV RNA: t-test p-value= 2E-16, Ploidy>2.25 is WGD+")+
#  theme_classic()
#t.test(CTSV..1515. ~ WGD, data = CCLE_RNA_Score)

#ggplot(CCLE_RNA_Score, 
#       aes(x=Gene_ploidy_Score, y=ST6GALNAC6..30815.))+
#  geom_point(color="black")+
#  ylab("CCLE RNA expresion: ST6GALNAC6..30815.") +
#  xlab("Aneuploidy Score: bp ploidy")+
#  theme_classic()


#potentially good: TMEM184A..202915., LAMA5..3911., TSPAN15..23555., 
#SHTN1..57698., RIPK4..54101., SOX13..9580., 
#CCLE_RNA_corr[19115,]
#R.top500.up<-CCLE_RNA_corr[1:501,] #highest  corrolation
#R.top500.down<-CCLE_RNA_corr[18615:19115,] 
#head(R.top500.up)

#a<- R.top500.up$rowname
#R.top500.up<- R.top500.up %>% separate(rowname, c("Name", "ID"))
#R.top500.up$rowname<-a

#b<-R.top500.down$rowname
#R.top500.down<- R.top500.down %>% separate(rowname, c("Name", "ID"))
#R.top500.down$rowname<-b

#write.csv(R.top500.up,
#          "CCLE_RNA_top500.csv", row.names = TRUE)
#write.csv(R.top500.down,
#          "CCLE_RNA_bottom500.csv", row.names = TRUE)

##### Step 5: Calculate p-values for RNA correlation with aneuploidy score #####
## because I calculated p-value for protein, I will also do so for RNA. 
## I wanted to get the p-value for all samples, and add to corr dataframe. 
# Step 1: get corrolation data between aneuploid score and expression per gene. 
df2<-data.frame(Name=NA, Pvalue=NA)

for (i in 3:(length(CCLE_RNA_Score)-6)) {
  t.test<-cor.test(CCLE_RNA_Score$Gene_ploidy_Score, 
                   CCLE_RNA_Score[,i], 
                   method = "pearson")
  pvalue<-t.test$p.value 
  name<-as.factor(colnames(CCLE_RNA_Score[i]))
  data<-data.frame(Name=name, Pvalue=pvalue)
  df2 <- rbind(df2,data)
}

df2$Name<-as.factor(df2$Name)
#step 2: combine with corr data
RNA.corr.pvalue<-merge(x= CCLE_RNA_corr, y= df2, 
                           by.x="rowname", by.y="Name", 
                           sort = TRUE)

# step 3: plot cor v. pvalue
RNA.corr.pvalue$RankPvalue<-rank(RNA.corr.pvalue$Pvalue)
RNA.corr.pvalue$B.H.CriticalValue<-(0.05/length(RNA.corr.pvalue$Pvalue))*RNA.corr.pvalue$RankPvalue
RNA.corr.pvalue$significant<- RNA.corr.pvalue$Pvalue< RNA.corr.pvalue$B.H.CriticalValue

RNA.corr.pvalue$RNA_Name<-sapply(strsplit(as.character(RNA.corr.pvalue$rowname),"[..]"), `[`, 1) #split at .., only take before.. character to get RNA name

RNA.corr.pvalue2 <- RNA.corr.pvalue[order(RNA.corr.pvalue$Pvalue,
                                          -RNA.corr.pvalue$Gene_ploidy_Score),]

write.csv(RNA.corr.pvalue2,
          "RNA_Corr_pvalue_GeneAneuploid.csv", row.names = TRUE)


ggplot(RNA.corr.pvalue2, 
                         aes(x=Gene_ploidy_Score, y=-log2(Pvalue), color=significant))+
  geom_point()+
  ylab("-log2 (P-Value)") +
  xlab("Correlation: RNA expression & aneuploidy score")+
  ggtitle ("Correlating aneuploidy scores and RNA expression difference, per gene")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  scale_color_manual(values=c("black", "goldenrod"))
#5x4

ggplot(RNA.corr.pvalue2, 
                      aes(x=Gene_ploidy_Score, y=Gene_Score, color=significant))+
  geom_density2d(size=0.5)+
  ylab("Correlation: RNA expression & aneuploidy score\nnot corrected for ploidy") +
  xlab("Correlation: RNA expression & aneuploidy score")+
  ggtitle ("Correlating aneuploidy scores and RNA expression difference, per gene")+
  theme(plot.title = element_text(hjust = 0.5))+
  coord_cartesian(xlim=c(-0.4, 0.4), ylim=c(-0.4,0.4))+
  theme_classic()+
  scale_color_manual(values=c("black", "goldenrod"))


ggplot(RNA.corr.pvalue2, 
       aes(x=Gene_ploidy_Score))+
  geom_density(size=1)+
  ylab("Density") +
  xlab("Correlation: RNA expression & aneuploidy score")+
  ggtitle ("Density plot of RNA aneuploidy score correlation")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()
#plot.density.RNA.AneuScore.Corr.all


# step 4: sort by p-value
RNA.corr.pvalue$RNA_Name<-sapply(strsplit(as.character(RNA.corr.pvalue$rowname),"[..]"), `[`, 1) #split at .., only take before.. character to get RNA name

RNA.corr.pvalue2 <- RNA.corr.pvalue[order(RNA.corr.pvalue$Pvalue,
                                                  -RNA.corr.pvalue$Gene_ploidy_Score),]

write.csv(RNA.corr.pvalue2,
          "RNA_Corr_pvalue_GeneAneuploid.csv", row.names = TRUE)

#Get significant p-value, positive correlation with aneuploidy
RNA.corr.pos.sig <- RNA.corr.pvalue[order(RNA.corr.pvalue$Pvalue),]
RNA.corr.pos.sig <- subset(RNA.corr.pos.sig, RNA.corr.pos.sig$significant==TRUE)
RNA.corr.pos.sig <- subset(RNA.corr.pos.sig, RNA.corr.pos.sig$Gene_ploidy_Score >0)

write.csv(RNA.corr.pos.sig,
          "RNA_Corr_pvalue_GeneAneuploid_PosSig.csv", row.names = TRUE)

# top2 most signif pos cor: TSPAN15, KIAA1522

#Get significant p-value, negative correlation with aneuploidy
RNA.corr.neg.sig <- RNA.corr.pvalue[order(RNA.corr.pvalue$Pvalue),]
RNA.corr.neg.sig <- subset(RNA.corr.neg.sig, RNA.corr.neg.sig$significant==TRUE)
RNA.corr.neg.sig <- subset(RNA.corr.neg.sig, RNA.corr.neg.sig$Gene_ploidy_Score <0)

write.csv(RNA.corr.neg.sig,
          "RNA_Corr_pvalue_GeneAneuploid_NegSig.csv", row.names = TRUE)
# RPL3..6122. 
# PPP3CC..5533. 
# RPL34..6164

# top 2 most sig neg corr: RPL3, PPP3CC

#setwd()#set working directory to depmap
#RNA.corr.pvalue2<-read.delim2("RNA_Corr_pvalue_GeneAneuploid.csv", 
#                              dec=".", header = TRUE, sep=",")


### Get aneuploidy correlation score data by chromosome
## sort my chromosome, give 1q data to Asad
RNA.corr.Chrm<-merge(x=Protein_Info2, y=RNA.corr.pvalue2, by.x="Approved Symbol", by.y="RNA_Name")
# Lost 513 genes when adding location data. somehow no location data? or names don't match up
RNA.corr.Chrm1<-subset(RNA.corr.Chrm, Chromosome==1) #1200 genes
#write.csv(RNA.corr.Chrm1,
#          "RNA_AneuploidCancer.Corr_1q.csv", row.names = TRUE)


#### plot RNA expression vs aneuploidy score
ggplot(CCLE_RNA_Score, 
       aes(x=Gene_ploidy_Score, y=RPL3..6122.))+ #
  geom_point(color="red", size=2)+
  ylab("RNA expression: RPL3") + 
  xlab("Gene ploidy score") +
  geom_smooth(color="black", method="lm") +
  #geom_text(x=6000, y=2, label="p-value= 2.4E-35")+
  #geom_text(x=6000, y=1.5, label="cor= 0.39")+ 
  theme_classic()


##### Step 6: Get gene groups, test sig more/loss corr than mean RNA ####
## to get gene groups, I downloaded GO terms from http://www.informatics.jax.org
## you can also get them from the supplemental data (excel format)

## Looking at HSPA and Ribosome RNAs
## top postively correlated RNA with aneuploidy score terms
##Get all HSPA / HSP70 RNA:
RNA.HSP <- subset(RNA.corr.pvalue2, grepl( "HSPA", RNA.corr.pvalue2$rowname, fixed = TRUE))
RNA.HSP <- RNA.HSP[order(RNA.HSP$Pvalue),]
RNA.HSP.NOT <- RNA.corr.pvalue2[!RNA.corr.pvalue2$RNA_Name %in% RNA.HSP$rowname,]# 

t.test(RNA.HSP$Gene_ploidy_Score, RNA.HSP.NOT$Gene_ploidy_Score)
#p=0.02907, upreg

## Get all Ribosomal RNA: 
RNA.corr.pvalue3<-RNA.corr.pvalue2[order(RNA.corr.pvalue2$rowname),]

RNA.Ribosome <- subset(RNA.corr.pvalue3, grepl( "RPS", RNA.corr.pvalue3$rowname, fixed = TRUE))
RNA.Ribosome <- RNA.Ribosome[35:75,]
RNA.Ribosome2 <- subset(RNA.corr.pvalue3, grepl( "RPL", RNA.corr.pvalue3$rowname, fixed = TRUE))
RNA.Ribosome2 <- RNA.Ribosome2[50:99,]
RNA.Ribosome<- rbind(RNA.Ribosome, RNA.Ribosome2) #this gives 91 genes
RNA.Ribosome.NOT <- RNA.corr.pvalue2[!RNA.corr.pvalue2$RNA_Name %in% RNA.Ribosome$rowname,]# 

t.test(RNA.Ribosome$Gene_ploidy_Score, RNA.Ribosome.NOT$Gene_ploidy_Score)
#P<2.2E-16



## Get all rRNA processing proteins GO:0006364
## 8 downregulated and 1 upregulated rRNA processing
#rRNAMet<- read.delim2("GO_term_rRNAprocessing.txt",
#                      dec=".", header = TRUE, sep="\t", row.names=NULL) 
#rRNAMet<- read_xlsx("Supplemental Data 6 gene lists.xlsx", sheet=9) 
rRNAMet.list<-unique(rRNAMet$Symbol)#get gene names
rRNAMet.list<-toupper(rRNAMet.list)#make string and uppercase. 209 genes
RNA.aneuScore.rRNAMet <- RNA.corr.pvalue2[RNA.corr.pvalue2$RNA_Name %in% rRNAMet.list,]# 
RNA.aneuScore.rRNAMet.NOT <- RNA.corr.pvalue2[!RNA.corr.pvalue2$RNA_Name %in% rRNAMet.list,]# 

t.test(RNA.aneuScore.rRNAMet$Gene_ploidy_Score, RNA.aneuScore.rRNAMet.NOT$Gene_ploidy_Score)
#NS



## CORUM gene list (protein complex)

CORUM.all2<- subset(CN.Diff.xRNA.yProt.ThreeGroups2, Uniprot_Acc %in% CORUM.all$Uniprot_ID)
CORUM.all3<-CORUM.all2$RNA_Name
RNA.Corum.Corr <- RNA.corr.pvalue2[RNA.corr.pvalue2$RNA_Name %in% CORUM.all3,]# 
RNA.Corum.Corr.NOT <- RNA.corr.pvalue2[!RNA.corr.pvalue2$RNA_Name %in% CORUM.all3,]# 

t.test(RNA.Corum.Corr$Gene_ploidy_Score, RNA.Corum.Corr.NOT$Gene_ploidy_Score)
# CORUM RNAs are negatively correlated with aneuploidy score
#(CORUM mean -0.18, all gene mean = 0.007, p<2.2E-16)


####intrinsic component of endoplasmic reticulum membrane, GO:0031227
#ERMem<- read.delim2("GO_term_ER_Membrane.txt",
#                    dec=".", header = TRUE, sep="\t", row.names=NULL) 
ERMem.list<-unique(ERMem$Symbol)#get gene names
ERMem.list<-toupper(ERMem.list)#make string and uppercase. 209 genes
R.aneuScore.ERMem <- RNA.corr.pvalue2[RNA.corr.pvalue2$RNA_Name %in% ERMem.list,]# 
R.aneuScore.ERMem.NOT <- RNA.corr.pvalue2[!RNA.corr.pvalue2$RNA_Name %in% ERMem.list,]# 

t.test(R.aneuScore.ERMem$Gene_ploidy_Score, R.aneuScore.ERMem.NOT$Gene_ploidy_Score)
# p= 5E-05
# good, upreg, but a lot. 
# 0.0022


### Othere gene groups we did not end up using: 
## Get all DNA double strand break repair proteins, go term GO:0006302
## 11 upregulated, 4 downregulated
#DNABreakRepair<- read.delim2("GO_term_DNADoubleStrandBreakRepair.txt",
#                             dec=".", header = TRUE, sep="\t", row.names=NULL) 
#DNABreakRepair.list<-unique(DNABreakRepair$MGI)#get gene names
#DNABreakRepair.list<-toupper(DNABreakRepair.list)#make string and uppercase. 243 genes
#RNA.aneuScore.DNABreakRepair <- RNA.corr.pvalue2[RNA.corr.pvalue2$RNA_Name %in% DNABreakRepair.list,]#183 found back
#RNA.aneuScore.DNABreakRepair.NOT <- RNA.corr.pvalue2[!RNA.corr.pvalue2$RNA_Name %in% DNABreakRepair.list,]# 

#t.test(RNA.aneuScore.DNABreakRepair$Gene_ploidy_Score, RNA.aneuScore.DNABreakRepair.NOT$Gene_ploidy_Score)
#P=1.6E-7, downregulated
## Get all Nonsense mediated decay proteins GO:0000184
## two significantly downregulated. most no change
#NMD<- read.delim2("GO_term_NonsenseMediatedDecay.txt",
#                  dec=".", header = TRUE, sep="\t", row.names=NULL) 
#NMD.list<-unique(NMD$MGI)#get gene names
#NMD.list<-toupper(NMD.list)#make string and uppercase. 36 genes
#RNA.aneuScore.NMD <- RNA.corr.pvalue[RNA.corr.pvalue$RNA_Name %in% NMD.list,]# 31 genes found back

## Plasma Membrane# slight upregulation in aneuploid, slight.
#PlasmaMembrane.list 
#RNA.aneuScore.PM <- RNA.corr.pvalue[RNA.corr.pvalue$RNA_Name %in% PlasmaMembrane.list,]# 

### cell cell singaling , NS
#CCSignal.list
#RNA.aneuScore.CCSignal <- RNA.corr.pvalue[RNA.corr.pvalue$RNA_Name %in% CCSignal.list,]# 

### Immune system , NS
#Immune.list
#RNA.aneuScore.Immune <- RNA.corr.pvalue[RNA.corr.pvalue$RNA_Name %in% Immune.list,]# 

### TF, NS 
#TF.list
#RNA.aneuScore.TF <- RNA.corr.pvalue[RNA.corr.pvalue$RNA_Name %in% TF.list,]# 

## Get all DNA double strand break repair proteins, go term GO:0006302
## 11 upregulated, 4 downregulated
#DNABreakRepair<- read.delim2("GO_term_DNADoubleStrandBreakRepair.txt",
#                             dec=".", header = TRUE, sep="\t", row.names=NULL) 
#DNABreakRepair.list<-unique(DNABreakRepair$MGI)#get gene names
#DNABreakRepair.list<-toupper(DNABreakRepair.list)#make string and uppercase. 243 genes
#RNA.aneuScore.DNArepair <- RNA.corr.pvalue2[RNA.corr.pvalue2$RNA_Name %in% DNABreakRepair.list,]# 
#RNA.aneuScore.DNArepair.NOT <- RNA.corr.pvalue2[!RNA.corr.pvalue2$RNA_Name %in% DNABreakRepair.list,]# 

#t.test(RNA.aneuScore.DNArepair$Gene_ploidy_Score, RNA.aneuScore.DNArepair.NOT$Gene_ploidy_Score)


### Cell adhesion 
#CellAdhesion.list
#RNA.aneuScore.CellAdhesion <- RNA.corr.pvalue[RNA.corr.pvalue$RNA_Name %in% CellAdhesion.list,]# 

### Tight Junction
#TightJunction.list
#RNA.aneuScore.TJ <- RNA.corr.pvalue[RNA.corr.pvalue$RNA_Name %in% TightJunction.list,]# 

## Cellular response to heat, GO:0034605, NS
#CellHeat.list
#RNA.aneuScore.Heat <- RNA.corr.pvalue[RNA.corr.pvalue$RNA_Name %in% CellHeat.list,]# 

# MRN complex (MRE11-RAD50-NBN complex), GO:0030870, NS
#MRNcomplex.list
#RNA.aneuScore.MRN <- RNA.corr.pvalue2[RNA.corr.pvalue2$RNA_Name %in% MRNcomplex.list,]# 
#RNA.aneuScore.MRN.NOT <- RNA.corr.pvalue2[!RNA.corr.pvalue2$RNA_Name %in% MRNcomplex.list,]# 

#t.test(RNA.aneuScore.MRN$Gene_ploidy_Score, RNA.aneuScore.MRN.NOT$Gene_ploidy_Score)
#NS


#### intrinsic component of plasma membrane, GO:0031226
#IntMem<- read.delim2("GO_term_Intrinsic.Plasma.Membrane.txt",
#                     dec=".", header = TRUE, sep="\t", row.names=NULL) 
#IntMem.list<-unique(IntMem$MGI)#get gene names
#IntMem.list<-toupper(IntMem.list)#make string and uppercase. 209 genes
#R.aneuScore.IntMem <- RNA.corr.pvalue2[RNA.corr.pvalue2$RNA_Name %in% IntMem.list,]# 
#R.aneuScore.IntMem.NOT <- RNA.corr.pvalue2[!RNA.corr.pvalue2$RNA_Name %in% IntMem.list,]# 

#t.test(R.aneuScore.IntMem$Gene_ploidy_Score, R.aneuScore.IntMem.NOT$Gene_ploidy_Score)
#good, upreg, but a lot. 
# 1E-10




#####         Tumor supressor gene (TSG) and oncogene (OG) groups ####
## Bailey Corr Oncogene: NS
## Vigano Corr Oncogene: NS
## Devoli Corr Oncogene: NS

## Bailey Corr TSG: neg corr P= 1.467E-5
## Vigano Corr TSG: neg corr P= 6.266E-7
## Devoli Corr TSG: neg corr P= 1.004E-7

# Bailey et al. Comprehensive Characterization of Cancer Driver Genes and Mutations, cell, 2018
Oncogene.Bailey.list#84
TSG.Bailey.list#99

RNA.Corr.Oncogene.Bailey <- RNA.corr.pvalue2[RNA.corr.pvalue2$RNA_Name %in% Oncogene.Bailey.list,]# 
RNA.Corr.TSG.Bailey <- RNA.corr.pvalue2[RNA.corr.pvalue2$RNA_Name %in% TSG.Bailey.list,]# 

#Now get all genes EXCEPT for oncogene/TSG. for t-test purposes
RNA.Corr.Oncogene.Bailey.NOT <- RNA.corr.pvalue2[!RNA.corr.pvalue2$RNA_Name %in% Oncogene.Bailey.list,]# 
RNA.Corr.TSG.Bailey.NOT <- RNA.corr.pvalue2[!RNA.corr.pvalue2$RNA_Name %in% TSG.Bailey.list,]# 

## T-test between oncogenes/TSG vs all others: more likely to be pos or neg corr with aneuploidy? 
t.test(RNA.Corr.Oncogene.Bailey$Gene_ploidy_Score, RNA.Corr.Oncogene.Bailey.NOT$Gene_ploidy_Score)
#Oncogene NS. not correlated with aneuploidy score
t.test(RNA.Corr.TSG.Bailey$Gene_ploidy_Score, RNA.Corr.TSG.Bailey.NOT$Gene_ploidy_Score)
#TSG more downregulated in more aneuploid cells. (higher aneuploidy, lower TSG expression)
# P= 1.467E-5


#####         Plot scatterplot correlation with aneuploidy score, per group ####
#RNA.aneuScore.DNABreakRepair   DNA Double strand break repair neg cor with aneuploidy, p=1.169e-07 (471 genes)
#RNA.Corum.Corr                 CORUM negatively correlated with aneuploidy. p<2.2E-16
#RNA.HSP                        HSPA 0.02907, upreg*
#RNA.aneuScore.MRN              MRN NS
#RNA.aneuScore.rRNAMet          rRNA process neg corr p 1.187e-12
#RNA.Ribosome                   Ribosome neg cor p<2.2E-16
#R.aneuScore.ERMem              intrinsic component of endoplasmic reticulum membrane, 0.0022

data.subset= R.aneuScore.ERMem
data.name = "intrinsic component of endoplasmic reticulum membrane"

ggplot(RNA.corr.pvalue2, 
                         aes(x=Gene_ploidy_Score, y=-log2(Pvalue), color=significant))+
  geom_point(size=1, color="black")+
  geom_point(size=2, color="dodgerblue3", data=data.subset,
             aes(x=Gene_ploidy_Score, y=-log2(Pvalue) ))+
  #geom_point(size=3, color="Cyan", data=data.subset, 
  #           aes(x=mean(Gene_ploidy_Score), y=mean(-log2(Pvalue)) ))+ #mean 
  ylab("-log2 (P-Value)") +
  xlab("correlation coefficient")+
  ggtitle (data.name)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  scale_color_manual(values=c("black", "red"))+
  scale_x_continuous(breaks=c(-0.4, -0.2, 0.0, 0.2, 0.4), limits=c(-0.4, 0.4))+
  scale_y_continuous(limits=c(0, 115))
#4x4 size  
# plot.RNA.aneuScoreCorr.pvalue.mean.all.Oncogene.Vigano
# plot.RNA.aneuScoreCorr.logPvalue.mean.CORUM


## Plot for gene group, but now with corr x corr
ggplot(RNA.corr.pvalue2, 
       aes(x=Gene_ploidy_Score, y=Gene_Score))+
  geom_point(size=1, color="black")+
  geom_point(size=2, color="red", data=data.subset,
             aes(x=Gene_ploidy_Score, y=Gene_Score ))+
  #geom_point(size=3, color="Cyan", data=data.subset, 
  #           aes(x=mean(Gene_ploidy_Score), y=mean(Gene_Score) ))+
  ylab("Correlation coefficient \nnot corrected for ploidy") +
  xlab("correlation coefficient")+
  ggtitle (data.name)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  scale_color_manual(values=c("black", "red"))+
  scale_x_continuous(breaks=c(-0.4, -0.2, 0.0, 0.2, 0.4), limits=c(-0.4, 0.4))+
  scale_y_continuous(breaks=c(-0.4, -0.2, 0.0, 0.2, 0.4), limits=c(-0.4, 0.4))

# 4x4 size  
# plot.RNA.aneuScoreCorr.CorrNoPloidy.all.ERMembrane.blue

# colors: negative = "dodgerblue3", positive = "Red"
