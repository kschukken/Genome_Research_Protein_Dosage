###### Comparing Mass Spec to RPPA data ####
### 211122
### by Klaske Schukken
### Compare RPPA (antibody data) and mass spec data
### 

library('ggplot2')
library('tidyverse')
library('xlsx')
library(readxl)
library(reshape2)
library('BBmisc')
library('tidyr')
library(tidyverse) 

###### Step 1: Get data ####
### Step 1: Get data. 

### Get RPPA data
### RPPA 
### download the RPPA protein expression data and format it according to the methods section.
setwd("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Controls_Lung_noMut_highRNA/RPPA control")
RPPA_data<- read.xlsx("RPPA aneuploidy.xlsx", sheetIndex = 1)

## Get CCLE raw Protein and RNA Data. 

## Import Protein expression data and uri ben-david lab, Cohen-Sharir nature 2021 chromosome arm data
##

DataFileLocation= "/Documents" ## ! Set file path to place you downloaded supplementary files and other suggested files

##New Protein data. Filtered for only cells with both RNA and Protein data: 
## generated in Protein_RNA_filtered_CellLine.R
setwd(DataFileLocation)
Protein.Expression.filtered<-read.delim2("Protein_Expression_filtered.csv", 
                                         dec=",", header = TRUE, sep=";")

#Name dataframe with protein_ID and corresponding Protein name, this way we have both gene ID and Protein name
Protein_Expression.1<-read_excel("mmc2.xlsx", 
                                 sheet= "Normalized Protein Expression")
Protein_Expression.1<-Protein_Expression.1[,c(1:426)] #delete empty collumns

Protein_ProID<-Protein_Expression.1[,c(1:2, 6)] 
Protein_ProID$Protein_Id<- str_replace_all(Protein_ProID$Protein_Id, "[|]", ".")
Protein_ProID$Protein_Id<- str_replace_all(Protein_ProID$Protein_Id, "[-]", ".")

#Protein Info. names, chrm location, etc. 
#setwd()
Protein_Info<-read_excel("HGNC_Protein_Info_edited.xlsx", sheet ="HGNC names")

#Add chrm arm data
Protein_Info4<-Protein_Info
Protein_Info4$arm <- gsub('[0-9]+', '', Protein_Info4$'Chromosome band') #Find test protein chromosome arm
Protein_Info4$arm <- gsub('[.]', '', Protein_Info4$arm) #also remove period, if needed
Protein_Info4$arm <- str_sub(Protein_Info4$arm, -1, -1) #start & end on last character. get only p or q
#filter cell data and get only aneuploidy scores for chrm arm of test protein location

###
#setwd()
aneuploid<- read.csv("arm_calls_data_bendavid.csv", header=TRUE)

###
# Import RNA expression data (no aneuploidy data). CCLE depmap data.

## RNA data. Filtered for only cells with both RNA and Protein data: 
#  RNA data for only cells that have both RNA and Protein Data available. (371 cell lines)
## a few extra cells are filtered out when we add aneuploidy data--> 367 cells total
## RNA expression differences generated in Protein_RNA_filtered_CellLine.R
#setwd()
setwd("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison")
RNA.Expression.filtered<-read.delim2("RNA_Expression_filtered.csv", 
                                     dec=",", header = TRUE, sep=";")

#RNA Info. names, chrm location, etc. 
#setwd()
RNA_Info<-read_excel("HGNC_Protein_Info_edited.xlsx", sheet ="HGNC names")
RNA_Info2<-RNA_Info[,c(2,3,11,12,22,14)]

RNA_Info3<-RNA_Info2
RNA_Info3$arm <- gsub('[0-9]+', '', RNA_Info3$'Chromosome band') #Find test RNA chromosome arm
RNA_Info3$arm <- gsub('[.]', '', RNA_Info3$arm) #also remove period, if needed
RNA_Info3$arm <- str_sub(RNA_Info3$arm, -1, -1) #start & end on last character. get only p or q

CN.Diff.xRNA.yProt.ThreeGroups<- read.csv("RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points_3cat.csv") 

setwd(DataFileLocation) # set wd to where you want the graphs to go to.


##### Correlate RPPA and mass spec difference data ####
colnames(RPPA_data) <- c("Protein_ID", "Gene_Symbol", "Protein.Diff.Gain.RPPA", "Protein.Diff.Loss.RPPA")
RPPA_data_MassSpec<- merge(RPPA_data, CN.Diff.xRNA.yProt.ThreeGroups, 
                           by.x="Gene_Symbol", by.y="RNA_Name")

## Protein difference upon chrm gain RPPA vs. mass spec
ggplot(RPPA_data_MassSpec, aes(x=Protein.Diff.Gain.RPPA, y=Protein.Diff.Gain))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  xlab("Difference in protein upon gain\nRPPA data")+
  ylab("Difference in protein upon gain\nMassSpec data")+
  theme_classic()+
  geom_hline(yintercept=0.0)+
  geom_vline(xintercept=0.0)+
  geom_smooth(method="lm", color="Red")+
  coord_cartesian(xlim=c(-0.7, 0.7), ylim=c(-0.7,0.7))+
  ggtitle("RPPA vs MassSpec: Gain")
# 4x4
# plot.2D.Diff.Gain.RPPA.MassSpec

## Protein difference upon chrm loss RPPA vs. mass spec
ggplot(RPPA_data_MassSpec, aes(x=Protein.Diff.Loss.RPPA, y=Protein.Diff.Loss))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  xlab("Difference in protein upon loss\nRPPA data")+
  ylab("Difference in protein upon loss\nMassSpec data")+
  theme_classic()+
  geom_hline(yintercept=0.0)+
  geom_vline(xintercept=0.0)+
  geom_smooth(method="lm", color="Red")+
  coord_cartesian(xlim=c(-0.7, 0.7), ylim=c(-0.7,0.7))+
  ggtitle("RPPA vs MassSpec: Loss")
# 4x4
# plot.2D.Diff.Loss.RPPA.MassSpec

# Corr gain: p= 1E-08, coef= 0.464
cor.test(RPPA_data_MassSpec$Protein.Diff.Gain.RPPA, RPPA_data_MassSpec$Protein.Diff.Gain)
# Corr loss: p=8E-11, coef= 0.517
cor.test(RPPA_data_MassSpec$Protein.Diff.Loss.RPPA, RPPA_data_MassSpec$Protein.Diff.Loss)

write.csv(RPPA_data_MassSpec, "RPPA.MassSpec.difference.csv")


# significant difference between data types? (NOPE, no significant difference)
t.test(RPPA_data_MassSpec$Protein.Diff.Gain.RPPA, RPPA_data_MassSpec$Protein.Diff.Gain)
# NS 
t.test(RPPA_data_MassSpec$Protein.Diff.Loss.RPPA, RPPA_data_MassSpec$Protein.Diff.Loss)
# NS 


##### Double buffered: RPPA and Mass Spec ####

RPPA.MassSpec.doubleBuffer.gain<-subset(RPPA_data_MassSpec, Protein.Diff.Gain.RPPA<0.25 & Protein.Diff.Gain<0.25)
# 74 genes double buffered
RPPA.MassSpec.doubleBuffer.loss<-subset(RPPA_data_MassSpec, Protein.Diff.Loss.RPPA> -0.25 & Protein.Diff.Loss> -0.25)
# 98 genes double buffered
## Ran this data through g:profiler GSEA (see gprofiler_bargraph_Var.R), 
# but there wasn't enough data to create significant enrichment

##### Bargraph of protein AS buffering Scaling profile 
RPPA_data.ThreeGroups<- RPPA_data

RPPA_data.ThreeGroups$Three.Protein.Gain<- cut(RPPA_data.ThreeGroups$Protein.Diff.Gain,
                                                        breaks=c(-Inf,-0.1,0.25,Inf),
                                                        include.lowest=TRUE,
                                                        labels=c("Anti-Scaling","Buffering","Scaling"))
RPPA_data.ThreeGroups$Three.Protein.Loss<- cut(RPPA_data.ThreeGroups$Protein.Diff.Loss,
                                                        breaks=c(-Inf,-0.25,0.1,Inf),
                                                        include.lowest=TRUE,
                                                        labels=c("Scaling", "Buffering","Anti-Scaling"))
RPPA_data.ThreeGroups$Three.Protein.Loss<-factor(RPPA_data.ThreeGroups$Three.Protein.Loss, 
                                                          levels=c("Anti-Scaling","Buffering","Scaling"))
#format data so I can plot both groups: Prot gain/loss
dat.m <- melt(RPPA_data.ThreeGroups, id.vars='Gene_Symbol', measure.vars=c('Three.Protein.Gain', 'Three.Protein.Loss'))

## RPPA Protein Gain
ggplot(data= dat.m, aes(x=variable, fill=value)) + 
  geom_bar(position="fill")+ #dodge (next to eachother)
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  #coord_cartesian(ylim=c(0,6000))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Percent of genes buffered upon gain and loss")
# 4x4
# plot.Protein.Bargraph.ThreeCat.Gain_RPPA

## Percent of genes AS buffer scaling upon gain: 
table(subset(dat.m, variable=="Three.Protein.Gain")$value)/length(subset(dat.m, variable=="Three.Protein.Gain")$value)
# Percent: 5.5, 74.4, 20
## Percent of genes AS buffer scaling upon loss: 
table(subset(dat.m, variable=="Three.Protein.Loss")$value)/length(subset(dat.m, variable=="Three.Protein.Gain")$value)
# Percent: 11.7, 75.9, 12.4


##### Boxplot of RPPA and mass spec protein expression difference #####

## Plot 2: barplot of difference: 
# RPPA_data_MassSpec
Diff.Data.RPPA<- data.frame(Dataset= as.character(), 
                                  DifferenceType= as.character(), 
                                  Difference = as.numeric())
# get data from RPPA and mass spec 
Diff.Data.RPPA<- rbind(Diff.Data.RPPA, 
                             data.frame(Dataset="Mass Spec", 
                                        DifferenceType= "Gain", 
                                        Difference= RPPA_data_MassSpec$Protein.Diff.Gain ))
Diff.Data.RPPA<- rbind(Diff.Data.RPPA, 
                             data.frame(Dataset="RPPA", 
                                        DifferenceType= "Gain", 
                                        Difference= RPPA_data_MassSpec$Protein.Diff.Gain.RPPA ))
Diff.Data.RPPA<- rbind(Diff.Data.RPPA, 
                             data.frame(Dataset="Mass Spec", 
                                        DifferenceType= "Loss", 
                                        Difference= RPPA_data_MassSpec$Protein.Diff.Loss ))
Diff.Data.RPPA<- rbind(Diff.Data.RPPA, 
                             data.frame(Dataset="RPPA", 
                                        DifferenceType= "Loss", 
                                        Difference= RPPA_data_MassSpec$Protein.Diff.Loss.RPPA ))

pdf(height=4, width=5, file="plot.bar.MeanDiff.MS.RPPA.pdf")
ggplot(Diff.Data.RPPA, aes(x=factor(DifferenceType), y=Difference))+ 
  geom_boxplot(aes(fill = Dataset), outlier.shape=NA)+ # no outliers
  xlab("")+
  ylab("Difference in expression")+
  theme_classic()+
  coord_cartesian(ylim=c(-0.6,.75))+
  geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  scale_fill_manual(name="Dataset", 
                    values=c("grey25", "grey75"))+
  ggtitle("Difference in protein expression by data type")
dev.off()



