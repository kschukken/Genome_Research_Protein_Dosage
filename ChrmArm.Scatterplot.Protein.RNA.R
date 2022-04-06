##### Plot scatterplots: gene expression difference for a chrm arm #####
## author: Klaske M. Schukken

## Protein data analysis: Chromosome gain or loss vs neutral ploidy expression difference
## get scatterplot Protein difference per chrm arm 

## 200914
## Edited: 201217 -- Normalized data
## Edited: 210106 -- Normalize by Gene, not Cell line. mean=0 only, not SD. 
## Edited: 210121 -- no normalization beyond what was already present in the dataset. 
## Edited: 210203: use filtered cell lines- only cells with both RNA and Protein data
##                also measure changes in protein expression by gain/loss chromosome ARM (not whole chrm)
## Edited 210502: used filtered genes so I only use genes that have minimum of 10 cells per category


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
# Step 2) scatter plot protein diff difference (black & red): all vs chrm 5q
#   plot RNA & Protein gain & loss: diff per gene, red=chrm arm 5q


##### Step 1: Get Protein difference for all genes per chromosome arm ####
###Step 1: Get data
## Protein changes upon chrm gain
## Get difference data per chrm arm: difference for all proteins (chrm 1-22) upon chrm X gain/loss: 
## Compare difference between proteins on chrm, vs not on chrm. 

setwd()

### Below files were generated in Protein_filtered_analysis_3cat.R*
## Below data sets are difference in all proteins (12k genes) (not filtered for 10+ cell per condition)
## upon chrm gain and loss. per chrm arm 

## sub-Step 1: get all the datasets for chrm Loss
P.Loss.Diff.Chrm1p<-read_csv(file= "ChrmArmLoss.1p.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm1q<-read_csv(file= "ChrmArmLoss.1q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm2p<-read_csv(file= "ChrmArmLoss.2p.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm2q<-read_csv(file= "ChrmArmLoss.2q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm3p<-read_csv(file= "ChrmArmLoss.3p.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm3q<-read_csv(file= "ChrmArmLoss.3q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm4p<-read_csv(file= "ChrmArmLoss.4p.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm4q<-read_csv(file= "ChrmArmLoss.4q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm5p<-read_csv(file= "ChrmArmLoss.5p.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm5q<-read_csv(file= "ChrmArmLoss.5q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm6p<-read_csv(file= "ChrmArmLoss.6p.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm6q<-read_csv(file= "ChrmArmLoss.6q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm7p<-read_csv(file= "ChrmArmLoss.7p.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm7q<-read_csv(file= "ChrmArmLoss.7q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm8p<-read_csv(file= "ChrmArmLoss.8p.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm8q<-read_csv(file= "ChrmArmLoss.8q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm9p<-read_csv(file= "ChrmArmLoss.9p.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm9q<-read_csv(file= "ChrmArmLoss.9q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm10p<-read_csv(file= "ChrmArmLoss.10p.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm10q<-read_csv(file= "ChrmArmLoss.10q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm11p<-read_csv(file= "ChrmArmLoss.11p.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm11q<-read_csv(file= "ChrmArmLoss.11q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm12p<-read_csv(file= "ChrmArmLoss.12p.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm12q<-read_csv(file= "ChrmArmLoss.12q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm13q<-read_csv(file= "ChrmArmLoss.13q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm14q<-read_csv(file= "ChrmArmLoss.14q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm15q<-read_csv(file= "ChrmArmLoss.15q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm16p<-read_csv(file= "ChrmArmLoss.16p.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm16q<-read_csv(file= "ChrmArmLoss.16q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm17p<-read_csv(file= "ChrmArmLoss.17p.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm17q<-read_csv(file= "ChrmArmLoss.17q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm18p<-read_csv(file= "ChrmArmLoss.18p.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm18q<-read_csv(file= "ChrmArmLoss.18q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm19p<-read_csv(file= "ChrmArmLoss.19p.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm19q<-read_csv(file= "ChrmArmLoss.19q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm20p<-read_csv(file= "ChrmArmLoss.20p.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm20q<-read_csv(file= "ChrmArmLoss.20q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm21q<-read_csv(file= "ChrmArmLoss.21q.Protein.Diff_min10cells.csv")
P.Loss.Diff.Chrm22q<-read_csv(file= "ChrmArmLoss.22q.Protein.Diff_min10cells.csv")

##sub-Step 2: repeat for chromosome gain changes
P.Gain.Diff.Chrm1p<-read_csv(file= "ChrmArmGain.1p.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm1q<-read_csv(file= "ChrmArmGain.1q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm2p<-read_csv(file= "ChrmArmGain.2p.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm2q<-read_csv(file= "ChrmArmGain.2q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm3p<-read_csv(file= "ChrmArmGain.3p.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm3q<-read_csv(file= "ChrmArmGain.3q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm4p<-read_csv(file= "ChrmArmGain.4p.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm4q<-read_csv(file= "ChrmArmGain.4q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm5p<-read_csv(file= "ChrmArmGain.5p.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm5q<-read_csv(file= "ChrmArmGain.5q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm6p<-read_csv(file= "ChrmArmGain.6p.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm6q<-read_csv(file= "ChrmArmGain.6q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm7p<-read_csv(file= "ChrmArmGain.7p.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm7q<-read_csv(file= "ChrmArmGain.7q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm8p<-read_csv(file= "ChrmArmGain.8p.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm8q<-read_csv(file= "ChrmArmGain.8q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm9p<-read_csv(file= "ChrmArmGain.9p.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm9q<-read_csv(file= "ChrmArmGain.9q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm10p<-read_csv(file= "ChrmArmGain.10p.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm10q<-read_csv(file= "ChrmArmGain.10q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm11p<-read_csv(file= "ChrmArmGain.11p.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm11q<-read_csv(file= "ChrmArmGain.11q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm12p<-read_csv(file= "ChrmArmGain.12p.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm12q<-read_csv(file= "ChrmArmGain.12q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm13q<-read_csv(file= "ChrmArmGain.13q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm14q<-read_csv(file= "ChrmArmGain.14q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm15q<-read_csv(file= "ChrmArmGain.15q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm16p<-read_csv(file= "ChrmArmGain.16p.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm16q<-read_csv(file= "ChrmArmGain.16q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm17p<-read_csv(file= "ChrmArmGain.17p.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm17q<-read_csv(file= "ChrmArmGain.17q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm18p<-read_csv(file= "ChrmArmGain.18p.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm18q<-read_csv(file= "ChrmArmGain.18q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm19p<-read_csv(file= "ChrmArmGain.19p.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm19q<-read_csv(file= "ChrmArmGain.19q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm20p<-read_csv(file= "ChrmArmGain.20p.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm20q<-read_csv(file= "ChrmArmGain.20q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm21q<-read_csv(file= "ChrmArmGain.21q.Protein.Diff.all_min10cells.csv")
P.Gain.Diff.Chrm22q<-read_csv(file= "ChrmArmGain.22q.Protein.Diff.all_min10cells.csv")


##### Step 2: Plot Protein expression change per chrm arm gain/loss #####
## Step 2: Plot Protein expression change per chrm arm gain/loss

# Plot mean Difference in expression of Protein between with gain chm vs no gain. 
# dotplot. one dot = one Protein Difference (gain - no gain) 
# color based on gain or no gain


data=P.Gain.Diff.Chrm5q
testChrmArm="5q" 

ggplot(data) + 
  geom_point(data = subset(data, Chrm.Num.Arm!=testChrmArm), 
             aes(x = Difference, y = -log2(p.value), 
                 col= Chrm.Num.Arm==testChrmArm))+
  geom_point(data = subset(data, Chrm.Num.Arm==testChrmArm), 
             aes(x = Difference, y = -log2(p.value), 
                 col= Chrm.Num.Arm==testChrmArm))+
  xlab(paste0("Protein Expression: \n(Chrm ", testChrmArm, " gain) - (no Chrm gain)"))+
  ylab("-log2(p-value)") +
  labs(col=paste0("Protein on Chrm Arm ", testChrmArm))+
  coord_cartesian(xlim=c(-4,4), ylim=c(0,25))+
  scale_color_manual(values=c("Black", "Red"))+
  theme_classic()+
  theme(axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))+
  ggtitle(paste0("Chromosome Arm ", testChrmArm, " gain vs no gain: \nProtein Difference vs. p-value"))
## dimentions: 5x4
## plot.Protein.Chrm5q.Gaun.scatter_min10cells_v3_dodgerblue3
## make loss dodgerblue3, gain =red


data=P.Loss.Diff.Chrm5q
testChrmArm="5q" 

ggplot(data) + 
  geom_point(data = subset(data, Chrm.Num.Arm!=testChrmArm), 
             aes(x = Difference, y = -log2(p.value), 
                 col= Chrm.Num.Arm==testChrmArm))+
  geom_point(data = subset(data, Chrm.Num.Arm==testChrmArm), 
             aes(x = Difference, y = -log2(p.value), 
                 col= Chrm.Num.Arm==testChrmArm))+
  xlab(paste0("Protein Expression: \n(Chrm ", testChrmArm, " loss) - (no Chrm loss)"))+
  ylab("-log2(p-value)") +
  labs(col=paste0("Protein on Chrm Arm ", testChrmArm))+
  coord_cartesian(xlim=c(-4,4), ylim=c(0,25))+
  scale_color_manual(values=c("Black", "dodgerblue3"))+
  theme_classic()+
  theme(axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))+
  ggtitle(paste0("Chromosome Arm ", testChrmArm, " loss vs no loss: \nProtein Difference vs. p-value"))
## dimentions: 5x4
## plot.Protein.Chrm5q.Loss.scatter_min10cells_v3_dodgerblue3
## make loss dodgerblue3, gain =red


##### Step 3: Get RNA difference for all genes per chromosome arm ####
### Step 3: Get data
### RNA changes upon chrm gain
### Get difference data per chrm arm: difference for all proteins (chrm 1-22) upon chrm X gain/loss: 
### Compare difference between proteins on chrm, vs not on chrm. 

# below files were make in RNA_filtered_analysis_3cat.R
# Below data sets are difference in all proteins (12k genes) (not filtered for 10+ cell per condition)
# upon chrm gain and loss. per chrm arm 

##Step 3: get all the datasets for chrm Loss
R.Loss.Diff.Chrm1p<-read_csv(file= "ChrmArmLoss.1p.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm1q<-read_csv(file= "ChrmArmLoss.1q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm2p<-read_csv(file= "ChrmArmLoss.2p.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm2q<-read_csv(file= "ChrmArmLoss.2q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm3p<-read_csv(file= "ChrmArmLoss.3p.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm3q<-read_csv(file= "ChrmArmLoss.3q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm4p<-read_csv(file= "ChrmArmLoss.4p.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm4q<-read_csv(file= "ChrmArmLoss.4q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm5p<-read_csv(file= "ChrmArmLoss.5p.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm5q<-read_csv(file= "ChrmArmLoss.5q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm6p<-read_csv(file= "ChrmArmLoss.6p.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm6q<-read_csv(file= "ChrmArmLoss.6q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm7p<-read_csv(file= "ChrmArmLoss.7p.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm7q<-read_csv(file= "ChrmArmLoss.7q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm8p<-read_csv(file= "ChrmArmLoss.8p.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm8q<-read_csv(file= "ChrmArmLoss.8q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm9p<-read_csv(file= "ChrmArmLoss.9p.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm9q<-read_csv(file= "ChrmArmLoss.9q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm10p<-read_csv(file= "ChrmArmLoss.10p.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm10q<-read_csv(file= "ChrmArmLoss.10q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm11p<-read_csv(file= "ChrmArmLoss.11p.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm11q<-read_csv(file= "ChrmArmLoss.11q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm12p<-read_csv(file= "ChrmArmLoss.12p.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm12q<-read_csv(file= "ChrmArmLoss.12q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm13q<-read_csv(file= "ChrmArmLoss.13q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm14q<-read_csv(file= "ChrmArmLoss.14q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm15q<-read_csv(file= "ChrmArmLoss.15q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm16p<-read_csv(file= "ChrmArmLoss.16p.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm16q<-read_csv(file= "ChrmArmLoss.16q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm17p<-read_csv(file= "ChrmArmLoss.17p.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm17q<-read_csv(file= "ChrmArmLoss.17q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm18p<-read_csv(file= "ChrmArmLoss.18p.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm18q<-read_csv(file= "ChrmArmLoss.18q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm19p<-read_csv(file= "ChrmArmLoss.19p.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm19q<-read_csv(file= "ChrmArmLoss.19q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm20p<-read_csv(file= "ChrmArmLoss.20p.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm20q<-read_csv(file= "ChrmArmLoss.20q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm21q<-read_csv(file= "ChrmArmLoss.21q.RNA.Diff.all.filtered_min10cells.csv")
R.Loss.Diff.Chrm22q<-read_csv(file= "ChrmArmLoss.22q.RNA.Diff.all.filtered_min10cells.csv")

##sub-Step 2: repeat for chromosome gain changes
R.Gain.Diff.Chrm1p<-read_csv(file= "ChrmArmGain.1p.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm1q<-read_csv(file= "ChrmArmGain.1q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm2p<-read_csv(file= "ChrmArmGain.2p.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm2q<-read_csv(file= "ChrmArmGain.2q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm3p<-read_csv(file= "ChrmArmGain.3p.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm3q<-read_csv(file= "ChrmArmGain.3q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm4p<-read_csv(file= "ChrmArmGain.4p.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm4q<-read_csv(file= "ChrmArmGain.4q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm5p<-read_csv(file= "ChrmArmGain.5p.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm5q<-read_csv(file= "ChrmArmGain.5q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm6p<-read_csv(file= "ChrmArmGain.6p.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm6q<-read_csv(file= "ChrmArmGain.6q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm7p<-read_csv(file= "ChrmArmGain.7p.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm7q<-read_csv(file= "ChrmArmGain.7q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm8p<-read_csv(file= "ChrmArmGain.8p.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm8q<-read_csv(file= "ChrmArmGain.8q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm9p<-read_csv(file= "ChrmArmGain.9p.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm9q<-read_csv(file= "ChrmArmGain.9q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm10p<-read_csv(file= "ChrmArmGain.10p.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm10q<-read_csv(file= "ChrmArmGain.10q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm11p<-read_csv(file= "ChrmArmGain.11p.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm11q<-read_csv(file= "ChrmArmGain.11q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm12p<-read_csv(file= "ChrmArmGain.12p.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm12q<-read_csv(file= "ChrmArmGain.12q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm13q<-read_csv(file= "ChrmArmGain.13q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm14q<-read_csv(file= "ChrmArmGain.14q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm15q<-read_csv(file= "ChrmArmGain.15q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm16p<-read_csv(file= "ChrmArmGain.16p.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm16q<-read_csv(file= "ChrmArmGain.16q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm17p<-read_csv(file= "ChrmArmGain.17p.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm17q<-read_csv(file= "ChrmArmGain.17q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm18p<-read_csv(file= "ChrmArmGain.18p.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm18q<-read_csv(file= "ChrmArmGain.18q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm19p<-read_csv(file= "ChrmArmGain.19p.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm19q<-read_csv(file= "ChrmArmGain.19q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm20p<-read_csv(file= "ChrmArmGain.20p.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm20q<-read_csv(file= "ChrmArmGain.20q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm21q<-read_csv(file= "ChrmArmGain.21q.RNA.Diff.all.filtered_min10cells.csv")
R.Gain.Diff.Chrm22q<-read_csv(file= "ChrmArmGain.22q.RNA.Diff.all.filtered_min10cells.csv")



##### Step 4: Scatterplot RNA difference upon chm arm gain: on chrm X vs not on chrm X ####
### Step 2: plot RNA difference upon chm arm gain: on chrm X vs not on chrm X
## Plot RNA expression change per chrm gain/loss

# Plot mean Difference in expression of RNA between with gain chm vs no gain. 
# dotplot. one dot = one RNA Difference (gain - no gain) 
# color based on gain or no gain

data=R.Gain.Diff.Chrm5q
testChrmArm="5q"

ggplot(data) + 
  geom_point(aes(x = Difference, y = -log2(p.value), col= Chrm.Num.Arm==testChrmArm))+
  geom_point(data = subset(data, Chrm.Num.Arm==testChrmArm), 
             aes(x = Difference, y = -log2(p.value), col= Chrm.Num.Arm==testChrmArm))+
  xlab(paste0("Difference in RNA Expression: \n(Chrm ", testChrmArm, " gain) - (no Chrm gain)"))+
  ylab("-log2(p-value)") +
  labs(col=paste0("RNA on Chrm Arm ", testChrmArm))+
  coord_cartesian(xlim=c(-2,2), ylim=c(0,50)) +
  scale_color_manual(values=c("Black", "Red"))+
  theme_classic()+
  theme(axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))+
  ggtitle(paste0("Chromosome Arm ", testChrmArm, " gain vs no gain: \nRNA Difference vs. p-value"))
## dimentions: 5x4
## plot.RNA.Chrm5q.Gain.scatter_min10cells_v3_dodgerblue
## loss= skyblue1, gain = red


data=R.Loss.Diff.Chrm5q
testChrmArm="5q"

ggplot(data) + 
  geom_point(aes(x = Difference, y = -log2(p.value), col= Chrm.Num.Arm==testChrmArm))+
  geom_point(data = subset(data, Chrm.Num.Arm==testChrmArm), 
             aes(x = Difference, y = -log2(p.value), col= Chrm.Num.Arm==testChrmArm))+
  xlab(paste0("Difference in RNA Expression: \n(Chrm ", testChrmArm, " loss) - (no Chrm loss)"))+
  ylab("-log2(p-value)") +
  labs(col=paste0("RNA on Chrm Arm ", testChrmArm))+
  coord_cartesian(xlim=c(-2,2), ylim=c(0,50)) +
  scale_color_manual(values=c("Black", "dodgerblue3"))+
  theme_classic()+
  theme(axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"))+
  ggtitle(paste0("Chromosome Arm ", testChrmArm, " loss vs no loss: \nRNA Difference vs. p-value"))
## dimentions: 5x4
## plot.RNA.Chrm5q.Loss.scatter_min10cells_v3_dodgerblue
## loss= skyblue1, gain = red


