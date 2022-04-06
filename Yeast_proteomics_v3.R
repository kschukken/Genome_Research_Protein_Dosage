##### Yeast Proteomics and transcriptomics #####
## 210507 
## Yeast Proteomics and transcriptomics- relative to CCLE/depmap dosage compensation
## Author: Klaske Schukken
## Compare chromosome loss/gain difference with protein expression changes. 
## Data from Dephoure et al. 2014, eLife; 3:e03023
## Quantitative proteomic analysis reveals posttranslational responses to aneuploidy in yeast


library('ggplot2')
library('tidyverse')
library(xlsx)
library(readxl)
library(reshape2)
library('BBmisc')

##### Step by step plan: ####

## Step 1: get Proteomics  data 
##      Yeast: several different aneuploid yeast clones relative to euploid
##      Data type: proteomics (log2 transformed difference between aneuploid & euploid)
##      Protein in aneuploid chromosomes (Depmap)
##      also get RNA data
## Step 2: scatterplot Protein difference in yeast vs cancer. --> density plot clearer
## Step 3: scatterplot RNA difference in yeast vs cancer. --> density plot clearer
## Step 4: Plot individual genes: yeast & CCLE protein & RNA expression diff upon gain and loss 


##### Step 1: get data #####


## This is the data Jason M. Sheltzer prepared. 
## data available in supplementary file 3, sheet named "Yeast_Data" 
# Yeast TMT means: 
setwd()
Yeast.Prot<- read_xlsx("yeast-human aneuploidy.xlsx") 


##### Step 2: scatterplot Protein difference in yeast vs cancer. #####
###Plot difference upon gain in yeast vs human
## Plot all aneuploid gene difference: Dephoure vs. depmap
# scatterplot
ggplot(Yeast.Prot, aes(y=`Human Tri vs. Di - protein`, x=`Yeast TMT - protein`))+
  geom_point(color="black")+
  theme_classic()+
  geom_hline(yintercept = 0)+ #add lines at y=0
  geom_vline(xintercept = 0)+ #add lines at x=0
  #coord_cartesian(ylim=c(-2,3), xlim=c(-2,3))+
  geom_smooth(method="lm", color="red")+ #linear trendline
  ylab("Depmap: Protein difference\nupon chrm gain")+
  xlab("Dephoure: Protein difference\nupon chrm gain in yeast")
# size 4x4 
# plot.Dephoure.Depmap.YeastMeanDiff_v2


## Plot all aneuploid gene difference: Dephoure vs. depmap
# same as above, but density plot
ggplot(Yeast.Prot, aes(y=`Human Tri vs. Di - protein`, x=`Yeast TMT - protein`))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  theme_classic()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  coord_cartesian(ylim=c(-.2,.5), xlim=c(-0.2,1.5))+
  geom_smooth(method="lm", color="red")+ #linear trendline
  ylab("Depmap: Protein difference\nupon chrm gain")+
  xlab("Dephoure: Protein difference\nupon chrm gain in yeast")
# 4x4
# plot.Dephoure.Depmap.YeastMeanDiff_density

# Get pearson correlation coefficient and p-value: 
YeasCorr<-cor.test(Yeast.Prot$`Yeast TMT - protein`, Yeast.Prot$`Human Tri vs. Di - protein`, method="pearson")
# Pearson corr= 0.214
# p-value = 4E-09
# n=738


##### Step 3: scatterplot RNA difference in yeast vs cancer. #####
## Plot all aneuploid gene difference: Dephoure vs. depmap
ggplot(Yeast.Prot, aes(y=`Human RNA`, x=`Yeast RNA`))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  theme_classic()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  coord_cartesian(ylim=c(0,.6), xlim=c(0,2))+
  geom_smooth(method="lm", color="red")+
  ylab("Depmap: Protein difference\nupon chrm gain")+
  xlab("Dephoure: Protein difference\nupon chrm gain in yeast")
# 4x4
# plot.Dephoure.Depmap.YeastMeanDiff.RNA_density

# Get pearson correlation coefficient and p-value: 
YeasCorr<-cor.test(Yeast.Prot$`Yeast RNA`, Yeast.Prot$`Human RNA`, method="pearson")
# Pearson corr= 0.073
# p-value = 0.049
# n=725


##### Step 4: Plot individual genes: yeast upon gain and loss ####

#prep data
Yeast.Prot$Protein.CCLE<-Yeast.Prot$`Human Tri vs. Di - protein`/mean(Yeast.Prot$`Human Tri vs. Di - protein`)
Yeast.Prot$Protein.Yeast<-Yeast.Prot$`Yeast TMT - protein`/mean(Yeast.Prot$`Yeast TMT - protein`)
Yeast.Prot$RNA.CCLE<-Yeast.Prot$`Human RNA`/mean(Yeast.Prot$`Human RNA`, na.rm=TRUE)
Yeast.Prot$RNA.Yeast<-Yeast.Prot$`Yeast RNA`/mean(Yeast.Prot$`Yeast RNA`, na.rm=TRUE)


gene=Yeast.Prot$`Human Gene`[58] # get a gene based on index (note: genes sorted by difference upon gain)
gene= "RPL3" #OR add name of gene you want to study
Yeast.Prot1<-subset(Yeast.Prot, `Human Gene` == gene) #only get data for gene of interest
Yeast.Prot2<-melt(Yeast.Prot1[,c(1,7:10)]) #prep data frame for plotting
Yeast.Prot2$variable<- factor(Yeast.Prot2$variable, level=c("Human RNA", "Yeast RNA", "Human Tri vs. Di - protein", "Yeast TMT - protein"))

#plot specific gene differences (RNA & protein, human & yeast)
ggplot(Yeast.Prot2, aes(x=variable, y=as.numeric(value)))+
  geom_col(position="dodge")+
  ggtitle(gene)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  ylab("Difference in gene expression")+
  xlab("")
# size: 4x4
# plot.Diff.yeast.human.perGene.RPL38

# some interesting genes: 
# EMC3
# PFDN4
# RPL3
# RPS15
# STX16
# SF3A1
# MRPS7
# YME1L1
# POLR2K, RPL38, or UQCRQ

# plot.YME1L1.yeast.human.protein.RNA
# 4x4




