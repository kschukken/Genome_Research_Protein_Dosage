###### TCGA DATA OVARIAN SAMPLES: Difference upon gain/loss dataset(s) ####
### and protein expression vs DNA copy number boxplots 
### 210127
### by Klaske Schukken
### Graphs for Depmap Protein & RNA data: 
### Graph: protein/RNA expression per cell in mono/disomy/trisomy categories
### Edited: 210204-- only in cell lines with both RNA and Protein data. 
##                -- and look at chromosome arm level gains/losses.
## Edited: 210630-- Look at TSGA data for Ovarian, ovarian and colorectal cancer. 

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
## Get raw Protein Data and TCGA aneuploidy scores.

DataFileLocation= "" 
setwd(DataFileLocation)

# TCGA aneuploidy data/ arm scores per sample from Taylor et al. 2018, supplementary table 2
TCGA_armScores<-read_xlsx("NIHMS958047-supplement-1_ArmCall.xlsx")


### Get Proteins: 
# CPTAC proteomics data that corresponds to TCGA samples. 
# from NIH NCI CPTAC
# data for Breat, ovarian, and colorectal data. 
# 111 Breast, 181 Ovarian, 84 Colorectal samples
# for Breast and Ovarian: use iTRAQ data- log transformed and normalized
# for colorectal: use Spectral Counts, normalize to total protein expression and log2 transform data. 

CPTAC_Ovarian<-read_tsv("TCGA_Ovarian_JHU_Proteome.itraq.tsv") # 132 samples, 8600 genes
CPTAC_Ovarian2<-read_tsv("TCGA_Ovarian_PNNL_Proteome.itraq.tsv") # 84 samples, 7483 genes

CPTAC_Protein_Expression<- merge(CPTAC_Ovarian[,1:265], CPTAC_Ovarian2, 
                                          by.x="Gene", by.y="Gene", all=TRUE)

#protein reported as "log ratio" and "log ratio unshared", remove the unshared collumns. repeats. 
CPTAC_Protein_Expression_rename <- CPTAC_Protein_Expression %>% select(-contains("Unshared"))
#now remove "Log Ratio" and/or "Spectral Counts", to just get sample name. 
# also, there are some repeat samples, use only first instance of that sample:
#      so remove .x, so that sample name with match armScore sample name list. 
colnames(CPTAC_Protein_Expression_rename) <- gsub(' Log Ratio', '', colnames(CPTAC_Protein_Expression_rename))
colnames(CPTAC_Protein_Expression_rename) <- gsub('.x', '', colnames(CPTAC_Protein_Expression_rename))
colnames(CPTAC_Protein_Expression_rename) <- gsub('-01A', '-01', colnames(CPTAC_Protein_Expression_rename))

# Filter CPTAC_Protein_Expression for samples with aneuploidy data
SampleColNames<-gsub('TCGA-', '', TCGA_armScores$Sample)
SampleColNames<- append(SampleColNames, c("Gene", "NCBIGeneID", "Authority", "Description", "Organism", "Chromosome", "Locus"))

CPTAC_Protein_Expression_ovarian<- CPTAC_Protein_Expression_rename[,(names(CPTAC_Protein_Expression_rename) %in% SampleColNames)]
# 168 samples have protein data and aneuploidy arm length data
# 8911 genes

TCGA_armScores_ovarian<- TCGA_armScores[gsub('TCGA-', '', TCGA_armScores$Sample) %in% colnames(CPTAC_Protein_Expression_ovarian),]
TCGA_armScores_ovarian$Sample<- gsub('TCGA-', '', TCGA_armScores_ovarian$Sample)
# 168 ovarian cancer samples with chromosome arm data. 

#Add chrm number + arm collumn
CPTAC_Protein_Expression_ovarian$arm<- gsub("[^a-zA-Z]", "", CPTAC_Protein_Expression_ovarian$Locus)
CPTAC_Protein_Expression_ovarian$arm<- gsub("X", "", CPTAC_Protein_Expression_ovarian$arm)
CPTAC_Protein_Expression_ovarian$arm<- gsub("qqq", "q", CPTAC_Protein_Expression_ovarian$arm)
CPTAC_Protein_Expression_ovarian$arm<- gsub("qq", "q", CPTAC_Protein_Expression_ovarian$arm)
CPTAC_Protein_Expression_ovarian$arm<- gsub("ppp", "p", CPTAC_Protein_Expression_ovarian$arm)
CPTAC_Protein_Expression_ovarian$arm<- gsub("pp", "p", CPTAC_Protein_Expression_ovarian$arm)
CPTAC_Protein_Expression_ovarian$arm<- gsub("pandYp", "p", CPTAC_Protein_Expression_ovarian$arm)

CPTAC_Protein_Expression_ovarian$ChrmNumArm<- paste0(CPTAC_Protein_Expression_ovarian$Chromosome, CPTAC_Protein_Expression_ovarian$arm) 
 #8911 genes

### Add BRCA1 chromosome location
#CPTAC_Protein_Expression_ovarian[,CPTAC_Protein_Expression_rename$Gene=="BRCA1"]
CPTAC_Protein_Expression_ovarian$Locus[837]<-"17q21.31"
CPTAC_Protein_Expression_ovarian$Chromosome[837]<-"17"
CPTAC_Protein_Expression_ovarian$ChrmNumArm[697]<-"17q"

subset(CPTAC_Protein_Expression_rename, Gene=="BRCA2")
CPTAC_Protein_Expression_ovarian$Locus[838]<-"13q13.1"
CPTAC_Protein_Expression_ovarian$Chromosome[838]<-"13"
CPTAC_Protein_Expression_ovarian$ChrmNumArm[838]<-"13q"


### Get RNA expression: 
# TCGA transcriptomic data that corresponds to TCGA samples. 
# from 
# RNA expressed as log2 lowess normalized (cy5/cy3) collapsed by gene symbol
setwd(DataFileLocation)
RNA_Ovarian<-read_tsv("OV.transcriptome__agilentg4502a_07_2_unc_edu_Level_3_data.txt") # 132 samples, 8600 genes
RNA_Ovarian2<-read_tsv("OV.transcriptome__agilentg4502a_07_3_unc_edu_Level_3_data.txt") # 84 samples, 7483 genes
#RNA_Ovarian<- RNA_Ovarian[-1,]
#RNA_Ovarian2<- RNA_Ovarian2[-1,]
colnames(RNA_Ovarian)[1]<-"Gene"
colnames(RNA_Ovarian2)[1]<-"Gene"
RNA_Expression<- merge(RNA_Ovarian, RNA_Ovarian2, 
                                 by.x="Gene", by.y="Gene", all=TRUE) 
#length= 598, 17814 genes


#protein reported as "log ratio" and "log ratio unshared", remove the unshared collumns. repeats. 
RNA_Expression_rename <- RNA_Expression %>% select(-contains("Unshared"))
#now remove "Log Ratio" and/or "Spectral Counts", to just get sample name. 
# also, there are some repeat samples, use only first instance of that sample:
#      so remove .x, so that sample name with match armScore sample name list. 
colnames(RNA_Expression_rename) <- gsub('.x', '', colnames(RNA_Expression_rename))
colnames(RNA_Expression_rename) <- gsub('-01A.*', '-01', colnames(RNA_Expression_rename))
colnames(RNA_Expression_rename) <- gsub('-01B.*', '-01', colnames(RNA_Expression_rename))
colnames(RNA_Expression_rename) <- gsub('-01C.*', '-01', colnames(RNA_Expression_rename))
colnames(RNA_Expression_rename) <- gsub('-01D.*', '-01', colnames(RNA_Expression_rename))
colnames(RNA_Expression_rename) <- gsub('-01R.*', '-01', colnames(RNA_Expression_rename))
colnames(RNA_Expression_rename) <- gsub('-02A.*', '-02', colnames(RNA_Expression_rename))
colnames(RNA_Expression_rename) <- gsub('-11A.*', '-11', colnames(RNA_Expression_rename))


# Filter CPTAC_Protein_Expression for samples with aneuploidy data
SampleColNames<-TCGA_armScores$Sample
SampleColNames<- append(SampleColNames, c("Gene", "NCBIGeneID", "Authority", "Description", "Organism", "Chromosome", "Locus"))

RNA_Expression_ovarian<- RNA_Expression_rename[,(names(RNA_Expression_rename) %in% SampleColNames)]
# 168 samples have protein data and aneuploidy arm length data
# 8911 genes

TCGA_armScores_ovarian<- TCGA_armScores[TCGA_armScores$Sample %in% colnames(RNA_Expression_ovarian),]
TCGA_armScores_ovarian$Sample<- TCGA_armScores_ovarian$Sample
# 168 ovarian cancer samples with chromosome arm data. 

GeneInfo <- CPTAC_Protein_Expression_ovarian[,c(1,176)] #gene names and ChrmNumArm

#Add chrm number + arm collumn
RNA_Expression_ovarian<- merge(RNA_Expression_ovarian, GeneInfo, 
                               by.x="Gene", by.y="Gene")
# 7427 genes shared with protein data. lost 1484 genes
# 542 cell lines (544 with gene and chrmNumArm)

### add BRCA1 chromosome location
#subset(RNA_Expression_ovarian, Gene=="BRCA1")
RNA_Expression_ovarian$ChrmNumArm[697]<-"17q"
#subset(RNA_Expression_ovarian, Gene=="BRCA2")
RNA_Expression_ovarian$ChrmNumArm[698]<-"13q"


###### Find Protein expression difference & t-test list, all ####
## make list of pvalue and difference between gain-neutral-loss cells per gene 
## for each gene in Protein data. using only cells that have Protein and aneuploidy data


t.test_TCGA_ovarian<- data.frame(Gene=character(),
                              Pvlaue.Gain= numeric(),
                              Diff.Gain= numeric(), 
                              CV.Gain=numeric(),
                              Pvlaue.Loss= numeric(),
                              Diff.Loss= numeric(),
                              CV.Loss= numeric(),
                              Pvalue.GainvLoss= numeric(),
                              Diff.GainvLoss= numeric())



for (i in 1:length(CPTAC_Protein_Expression_ovarian$Gene)){ #exclude mean, median st dev rows
  TestProt = CPTAC_Protein_Expression_ovarian$Gene[i]
  
  testProtChrm <- CPTAC_Protein_Expression_ovarian$ChrmNumArm[i] #get test protein data
  
  Protein_Expression_TestProt <- CPTAC_Protein_Expression_ovarian[i,c(2:168)] #sample collumns only
  
  ## Get samples ID with gain, loss or neutral for text chrm arm
  if (testProtChrm %in% colnames(TCGA_armScores_ovarian)) {
  AneuploidData<- TCGA_armScores_ovarian %>%
    select("Sample", all_of(testProtChrm))
  
  GainChrm <- subset(AneuploidData, AneuploidData[,2]== 1)
  NeutralChrm <- subset(AneuploidData, AneuploidData[,2]== 0)
  LossChrm <- subset(AneuploidData, AneuploidData[,2]== -1)
  
  ## put data into categories based on chrm arm number
  Chrm.tri<- Protein_Expression_TestProt %>%
    select(GainChrm$Sample)
  Chrm.di<- Protein_Expression_TestProt %>%
    select(NeutralChrm$Sample)
  Chrm.mono<- Protein_Expression_TestProt %>%
    select(LossChrm$Sample)
  
  ##Make data frame with t-test info about protein expression per arm number category.
  ## Return this data at the end of the function. 
  if (rowSums(!is.na(Chrm.tri))>=10 & #[8] is column with protein expression data. check it has values. 
      #Added if statement so only genes with 2+ values per condition are analyzed
      #if I don't do this, the t-test crashes and I get no values. 
      rowSums(!is.na(Chrm.di))>=10 &
      rowSums(!is.na(Chrm.mono))>=10) {
    # Trisomy vs disomy
    di.tri<-t.test(Chrm.tri, 
                   Chrm.di, 
                   mu = 0, 
                   alt = "two.sided",
                   conf.level = 0.99) #get p-value
    Diff.Tri.Di<-rowMeans(Chrm.tri, na.rm=TRUE) - rowMeans(Chrm.di, na.rm=TRUE) #get difference
    Variance.Gain<- var(as.numeric(Chrm.tri), na.rm=TRUE)/rowMeans(Chrm.tri, na.rm=TRUE)
    #Disomy vs monosome
    Di.Mono<-t.test(Chrm.mono, 
                    Chrm.di, 
                    mu = 0, 
                    alt = "two.sided",
                    conf.level = 0.99) #get p-value
    Diff.Di.Mono<- rowMeans(Chrm.mono, na.rm=TRUE) - rowMeans(Chrm.di, na.rm=TRUE)  #get difference
    Variance.Loss<- var(as.numeric(Chrm.mono), na.rm=TRUE)/rowMeans(Chrm.mono, na.rm=TRUE)
    # Tri vs monosomy
    Tri.Mono<-t.test(Chrm.mono, 
                     Chrm.tri, 
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
    Diff.Tri.Mono<-rowMeans(Chrm.tri, na.rm=TRUE) - rowMeans(Chrm.mono, na.rm=TRUE) #get difference
    
    t.test_TCGA_ovarian<- rbind(t.test_TCGA_ovarian, 
                             data.frame(Gene=TestProt,
                                        Pvalue.Gain= di.tri$p.value,
                                        Diff.Gain= Diff.Tri.Di, 
                                        CV.Gain= Variance.Gain,
                                        Pvalue.Loss= Di.Mono$p.value,
                                        Diff.Loss= Diff.Di.Mono, 
                                        CV.Loss= Variance.Loss,
                                        Pvalue.GainvLoss= Tri.Mono$p.value,
                                        Diff.GainvLoss= Diff.Tri.Mono))
  } else {
    t.test_TCGA_ovarian<-t.test_TCGA_ovarian
  }
  }
}


t.test_TCGA_ovarian<-t.test_TCGA_ovarian[order(t.test_TCGA_ovarian$Pvalue.GainvLoss),]
length(t.test_TCGA_ovarian$Gene) #1509 genes pass 10+ cells/category test

# 5216 genes pass the 3+ sample/category test
# CONTROL normalized: 1196 genes pass the 10+ sample/category test

setwd(DataFileLocation)
#write.csv(t.test_TCGA_ovarian, 
#          file =paste("Protein_TCGA_OvarianCancer_Diff_Pvalue_min10points.csv", sep=','), 
#          row.names = TRUE)

setwd()
t.test_TCGA_ovarian<-read.delim2("Protein_TCGA_OvarianCancer_Diff_Pvalue_min10points.csv", 
                                  dec=".", header = TRUE, sep=",")


###### Find RNA expression difference & t-test list, all ####
## make list of pvalue and difference between gain-neutral-loss cells per gene 
## for each gene in Protein data. using only cells that have Protein and aneuploidy data


t.test_TCGA_ovarian.R<- data.frame(Gene=character(),
                                 Pvlaue.Gain.R= numeric(),
                                 Diff.Gain.R= numeric(), 
                                 CV.Gain.R=numeric(),
                                 Pvlaue.Loss.R= numeric(),
                                 Diff.Loss.R= numeric(),
                                 CV.Loss.R= numeric(),
                                 Pvalue.GainvLoss.R= numeric(),
                                 Diff.GainvLoss.R= numeric())



for (i in 1:length(RNA_Expression_ovarian$Gene)){ #exclude mean, median st dev rows
  TestProt = RNA_Expression_ovarian$Gene[i]
  
  testProtChrm <- RNA_Expression_ovarian$ChrmNumArm[i] #get test protein data
  
  Protein_Expression_TestProt <- RNA_Expression_ovarian[i,c(2:543)] #sample collumns only
  
  ## Get samples ID with gain, loss or neutral for text chrm arm
  if (testProtChrm %in% colnames(TCGA_armScores_ovarian)) {
    AneuploidData<- TCGA_armScores_ovarian %>%
      select("Sample", all_of(testProtChrm))
    
    GainChrm <- subset(AneuploidData, AneuploidData[,2]== 1)
    NeutralChrm <- subset(AneuploidData, AneuploidData[,2]== 0)
    LossChrm <- subset(AneuploidData, AneuploidData[,2]== -1)
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- as.numeric(Protein_Expression_TestProt %>%
      select(GainChrm$Sample) )
    Chrm.di<- as.numeric( Protein_Expression_TestProt %>%
      select(NeutralChrm$Sample) )
    Chrm.mono<- as.numeric( Protein_Expression_TestProt %>%
      select(LossChrm$Sample) )
    
    ##Make data frame with t-test info about protein expression per arm number category.
    ## Return this data at the end of the function. 
    if (length(!is.na(Chrm.tri))>=10 & #[8] is column with protein expression data. check it has values. 
        #Added if statement so only genes with 10+ values per condition are analyzed 
        length(!is.na(Chrm.di))>=10 &
        length(!is.na(Chrm.mono))>=10) {
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri, 
                     Chrm.di, 
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri, na.rm=TRUE) - mean(Chrm.di, na.rm=TRUE) #get difference
      Variance.Gain<- var(as.numeric(Chrm.tri), na.rm=TRUE)/mean(Chrm.tri, na.rm=TRUE)
      #Disomy vs monosome
      Di.Mono<-t.test(Chrm.mono, 
                      Chrm.di, 
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Di.Mono<- mean(Chrm.mono, na.rm=TRUE) - mean(Chrm.di, na.rm=TRUE)  #get difference
      Variance.Loss<- var(as.numeric(Chrm.mono), na.rm=TRUE)/mean(Chrm.mono, na.rm=TRUE)
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono, 
                       Chrm.tri, 
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri, na.rm=TRUE) - mean(Chrm.mono, na.rm=TRUE) #get difference
      
      t.test_TCGA_ovarian.R<- rbind(t.test_TCGA_ovarian.R, 
                                  data.frame(Gene=TestProt,
                                             Pvalue.Gain.R= di.tri$p.value,
                                             Diff.Gain.R= Diff.Tri.Di, 
                                             CV.Gain.R= Variance.Gain,
                                             Pvalue.Loss.R= Di.Mono$p.value,
                                             Diff.Loss.R= Diff.Di.Mono, 
                                             CV.Loss.R= Variance.Loss,
                                             Pvalue.GainvLoss.R= Tri.Mono$p.value,
                                             Diff.GainvLoss.R= Diff.Tri.Mono))
    } else {
      t.test_TCGA_ovarian.R<-t.test_TCGA_ovarian.R
    }
  }
}


t.test_TCGA_ovarian.R<-t.test_TCGA_ovarian.R[order(t.test_TCGA_ovarian.R$Pvalue.GainvLoss),]
length(t.test_TCGA_ovarian.R$Gene) #1509 genes pass 10+ cells/category test

# 4994 genes pass the 10+ sample/category test

setwd(DataFileLocation)
write.csv(t.test_TCGA_ovarian.R, 
          file =paste("RNA_TCGA_OvarianCancer_Diff_Pvalue_min10points.csv", sep=','), 
          row.names = TRUE)

#t.test_TCGA_ovarian.R<-read.delim2("RNA_TCGA_OvarianCancer_Diff_Pvalue_min10points.csv", 
#                                  dec=".", header = TRUE, sep=",")

write.csv(t.test_TCGA_ovarian.ThreeGroups.R, 
          file =paste("RNA_TCGA_Difference_Pvalue_min10points_3cat.csv", sep=','), 
          row.names = TRUE)

t.test_TCGA_ovarian.ThreeGroups.R<- read.delim2("RNA_TCGA_Difference_Pvalue_min10points_3cat.csv", 
                                             dec=".", header = TRUE, sep=",")

###### Combine Protein and RNA data #####

TCGA_ovarian_Diff.Pvalue<- merge(t.test_TCGA_ovarian.R, t.test_TCGA_ovarian, 
                                 by.x="Gene", by.y="Gene")
# 1304 genes have both RNA and Protein data

#setwd(DataFileLocation)
#write.csv(t.test_TCGA_ovarian.ThreeGroups, 
#          file="TCGA_ovarian.diff.pvalue.RNA.Protein.ThreeGroups.csv")

t.test_TCGA_ovarian.ThreeGroups<- read.csv("TCGA_ovarian.diff.pvalue.RNA.Protein.ThreeGroups.csv")

###### Categorize genes by significance and difference (not used in paper) #####
## Categorize genes by if they scale or buffer upon chrm gain/loss. plot. 
## Did not end up using these graphs in the paper. but can use this to find highly 
##    significant scaling/buffering proteins/RNAs. 

## Find genes per category: Per chrm gain: 
## Categories based on Significance and +/- Difference,
## can be used to find genes that are significantly Scaling/anti-scaling (or NS) for gain & loss

CN.gain.ProtNS<- subset(t.test_TCGA_ovarian, 
                                    (Pvalue.Gain>0.05) )
CN.gain.ProtNS<- CN.gain.ProtNS[order(CN.gain.ProtNS$Pvalue.Gain),]
CN.gain.Protscale<- subset(t.test_TCGA_ovarian, 
                                      (Pvalue.Gain<0.05 & Diff.Gain>0) )
CN.gain.Protscale<- CN.gain.Protscale[order(CN.gain.Protscale$Pvalue.Gain),]

CN.gain.ProtAS<- subset(t.test_TCGA_ovarian, 
                                     (Pvalue.Gain<0.05 & Diff.Gain<0) )
CN.gain.ProtAS<- CN.gain.ProtAS[order(CN.gain.ProtAS$Pvalue.Gain),]

length(CN.gain.ProtNS$Gene)  # 2885
CN.gain.ProtNS$Gene[1:5]
length(CN.gain.Protscale$Gene) # 336
CN.gain.Protscale$Gene[1:5]
length(CN.gain.ProtAS$Gene) # 32
CN.gain.ProtAS$Gene[1:5]

## Upon chromosome gain: 
##Minimum of 10 data points per condition (-1, 0, +1):
# NS 89% of genes
# gain significant: 10%
#   GSPT1 USP7  PARN  MKL2  PGAM5
# Loss significant: 1%
# SCIN   CLC    SCGN   SLC1A4 CD300A

####
##Find categories for chromosome loss: 
CN.loss.ProtNS<- subset(t.test_TCGA_ovarian, 
                        (Pvalue.Loss>0.05) )
CN.loss.ProtNS<- CN.loss.ProtNS[order(CN.loss.ProtNS$Pvalue.Loss),]
CN.loss.Protscale<- subset(t.test_TCGA_ovarian, 
                           (Pvalue.Loss<0.05 & Diff.Loss<0) )
CN.loss.Protscale<- CN.loss.Protscale[order(CN.loss.Protscale$Pvalue.Loss),]

CN.loss.ProtAS<- subset(t.test_TCGA_ovarian, 
                        (Pvalue.Loss<0.05 & Diff.Loss>0) )
CN.loss.ProtAS<- CN.loss.ProtAS[order(CN.loss.ProtAS$Pvalue.Loss),]

length(CN.gain.ProtNS$Gene)  # 2657
CN.gain.ProtNS$Gene[1:5]
length(CN.gain.Protscale$Gene) # 474
CN.gain.Protscale$Gene[1:5]
length(CN.gain.ProtAS$Gene) # 122
CN.gain.ProtAS$Gene[1:5]

## Upon chromosome gain: 
##Minimum of 10 data points per condition (-1, 0, +1):
# NS 82% of genes
# Loss significant: 15%
#  PCBD2    COL4A3BP KRT24    PLXDC1   MTAP 
# Gain significant: 4%
# KPNA2  RPL23  EFTUD2 UTP15  PSMC5 


###Now to look at genes that are buffeted/scaling in BOTH gain and loss 
CN.aneu.gain.loss.Scale<- subset(t.test_TCGA_ovarian, 
                                    ((Pvalue.Loss<0.05 & Diff.Loss<0) & (Pvalue.Gain<0.05 & Diff.Gain>0) ))  # sig and scale
CN.aneu.gain.loss.Scale<- CN.aneu.gain.loss.Scale[order(CN.aneu.gain.loss.Scale$Pvalue.GainvLoss),]

CN.aneu.gain.loss.AS<- subset(t.test_TCGA_ovarian, 
                                 ((Pvalue.Loss<0.05 & Diff.Loss>0) & (Pvalue.Gain<0.05 & Diff.Gain<0) ))  # sig and scale
CN.aneu.gain.loss.AS<- CN.aneu.gain.loss.AS[order(CN.aneu.gain.loss.AS$Pvalue.GainvLoss),]


length(CN.aneu.gain.loss.Scale$Gene)  #42
CN.aneu.RNAbuff.Protbuff$Gene[1:5]
length(CN.aneu.gain.loss.AS$Gene)  #0
CN.aneu.gain.loss.AS$Gene[1:5]



###### Categorize genes: Scaling, buffering. anti-scaling ####
## Make THREE Categories: scaling, anti-scale and buffering, based on difference. 
## Three groups: Scaling, buffering and anti-scaling
## Three categories, upon gain: sclaing (>0.25), buffering (-0.1 to 0.25), anti-scaling (< -0.1)
## 0.25 because thats halfway between .5 fold change, expected from 1.5-fold increase. 
## -0.1 to account for non-significant downregulation upon gain
## vice versa for chrm loss. 

### Protein
t.test_TCGA_ovarian.ThreeGroups<- TCGA_ovarian_Diff.Pvalue
t.test_TCGA_ovarian.ThreeGroups$Three.Protein.Gain<- cut(t.test_TCGA_ovarian.ThreeGroups$Diff.Gain,
                                                    breaks=c(-Inf,-0.1,0.25,Inf),
                                                    include.lowest=TRUE,
                                                    labels=c("Anti-Scaling","Buffering","Scaling"))
t.test_TCGA_ovarian.ThreeGroups$Three.Protein.Loss<- cut(t.test_TCGA_ovarian.ThreeGroups$Diff.Loss,
                                                    breaks=c(-Inf,-0.25,0.1,Inf),
                                                    include.lowest=TRUE,
                                                    labels=c("Scaling", "Buffering","Anti-Scaling"))
t.test_TCGA_ovarian.ThreeGroups$Three.RNA.Gain<- cut(t.test_TCGA_ovarian.ThreeGroups$Diff.Gain.R,
                                                         breaks=c(-Inf,-0.1,0.25,Inf),
                                                         include.lowest=TRUE,
                                                         labels=c("Anti-Scaling","Buffering","Scaling"))
t.test_TCGA_ovarian.ThreeGroups$Three.RNA.Loss<- cut(t.test_TCGA_ovarian.ThreeGroups$Diff.Loss.R,
                                                         breaks=c(-Inf,-0.25,0.1,Inf),
                                                         include.lowest=TRUE,
                                                         labels=c("Scaling", "Buffering","Anti-Scaling"))

t.test_TCGA_ovarian.ThreeGroups$Three.Protein.Loss<-factor(t.test_TCGA_ovarian.ThreeGroups$Three.Protein.Loss, 
                                                           levels=c("Anti-Scaling","Buffering","Scaling"))
t.test_TCGA_ovarian.ThreeGroups$Three.Protein.Gain<-factor(t.test_TCGA_ovarian.ThreeGroups$Three.Protein.Gain, 
                                                           levels=c("Anti-Scaling","Buffering","Scaling"))
t.test_TCGA_ovarian.ThreeGroups$Three.RNA.Loss<-factor(t.test_TCGA_ovarian.ThreeGroups$Three.RNA.Loss, 
                                                          levels=c("Anti-Scaling","Buffering","Scaling"))
t.test_TCGA_ovarian.ThreeGroups$Three.RNA.Gain<-factor(t.test_TCGA_ovarian.ThreeGroups$Three.RNA.Gain, 
                                                      levels=c("Anti-Scaling","Buffering","Scaling"))



### Gain count percentages: 6% AS, 86% Buffer, 8% Scale
sum(t.test_TCGA_ovarian.ThreeGroups$Three.Protein.Gain=="Anti-Scaling")/length(t.test_TCGA_ovarian.ThreeGroups$Gene) #4%
sum(t.test_TCGA_ovarian.ThreeGroups$Three.Protein.Gain=="Buffering")/length(t.test_TCGA_ovarian.ThreeGroups$Gene) #80%
sum(t.test_TCGA_ovarian.ThreeGroups$Three.Protein.Gain=="Scaling")/length(t.test_TCGA_ovarian.ThreeGroups$Gene) # 16%

### Loss count percentages: 4% AS, 80% Buffer, 16% Scale
sum(t.test_TCGA_ovarian.ThreeGroups$Three.Protein.Loss=="Anti-Scaling")/length(t.test_TCGA_ovarian.ThreeGroups$Gene) #4%
sum(t.test_TCGA_ovarian.ThreeGroups$Three.Protein.Loss=="Buffering")/length(t.test_TCGA_ovarian.ThreeGroups$Gene) #80%
sum(t.test_TCGA_ovarian.ThreeGroups$Three.Protein.Loss=="Scaling")/length(t.test_TCGA_ovarian.ThreeGroups$Gene) # 16%

### RNA Gain count percentages: 4% AS, 52% Buffer, 44% Scale
sum(t.test_TCGA_ovarian.ThreeGroups$Three.RNA.Gain=="Anti-Scaling")/length(t.test_TCGA_ovarian.ThreeGroups$Gene) #4%
sum(t.test_TCGA_ovarian.ThreeGroups$Three.RNA.Gain=="Buffering")/length(t.test_TCGA_ovarian.ThreeGroups$Gene) #80%
sum(t.test_TCGA_ovarian.ThreeGroups$Three.RNA.Gain=="Scaling")/length(t.test_TCGA_ovarian.ThreeGroups$Gene) # 16%

### Loss count percentages: 3% AS, 27% Buffer, 71% Scale
sum(t.test_TCGA_ovarian.ThreeGroups$Three.RNA.Loss=="Anti-Scaling")/length(t.test_TCGA_ovarian.ThreeGroups$Gene) #4%
sum(t.test_TCGA_ovarian.ThreeGroups$Three.RNA.Loss=="Buffering")/length(t.test_TCGA_ovarian.ThreeGroups$Gene) #80%
sum(t.test_TCGA_ovarian.ThreeGroups$Three.RNA.Loss=="Scaling")/length(t.test_TCGA_ovarian.ThreeGroups$Gene) # 16%


##Scatterplots with categories colored in: 
# Chrm gain graphs
ggplot(t.test_TCGA_ovarian.ThreeGroups, aes(x=Diff.Gain, y=-log2(Pvalue.Gain), color=Three.Protein.Gain))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in Protein expression")+
  ylab("log2(p-value)")+
  labs(color = "Category:")+ #legend title
  theme_classic()+
  coord_cartesian(xlim=c(-1, 1), ylim=c(0,45))+
  ggtitle("Ovarian cancer: Protein upon chrm gain")
# plot.Protein.TCGA.ovarian.Diff.pvalue.category.Gain

# Loss graphs
ggplot(t.test_TCGA_ovarian.ThreeGroups, aes(x=Diff.Loss, y=-log2(Pvalue.Loss), color=Three.Protein.Loss))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in Protein expression")+
  ylab("log2(p-value)")+
  labs(color = "Category:")+ #legend title
  theme_classic()+
  coord_cartesian(xlim=c(-1, 1), ylim=c(0,45))+
  ggtitle("Ovarian cancer: Protein upon chrm loss")
#plot.Protein.TCGA.ovarian.Diff.pvalue.category.Loss
# 5x4

##Scatterplots with categories colored in: RNA
# Chrm gain graphs
ggplot(t.test_TCGA_ovarian.ThreeGroups, aes(x=Diff.Gain.R, y=-log2(Pvalue.Gain.R), color=Three.RNA.Gain))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in RNA expression")+
  ylab("log2(p-value)")+
  labs(color = "Category:")+ #legend title
  theme_classic()+
  coord_cartesian(xlim=c(-1, 1), ylim=c(0,160))+
  ggtitle("Ovarian cancer: RNA upon chrm gain")
# plot.TCGA.ovarian.RNA.Diff.pvalue.category.Gain

# RNA Loss graphs
ggplot(t.test_TCGA_ovarian.ThreeGroups, aes(x=Diff.Loss.R, y=-log2(Pvalue.Loss.R), color=Three.RNA.Loss))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in RNA expression")+
  ylab("log2(p-value)")+
  labs(color = "Category:")+ #legend title
  theme_classic()+
  coord_cartesian(xlim=c(-1, 1), ylim=c(0,160))+
  ggtitle("Ovarian cancer: RNA upon chrm loss")
#plot.TCGA.ovarian.RNA.Diff.pvalue.category.Loss
# 5x4

### variance
## did not end up using variance graphs
# gain graphs
ggplot(t.test_TCGA_ovarian.ThreeGroups, aes(x=Diff.Gain, y=-log(abs(CV.Gain)), color= -log(abs(CV.Gain))< 0 ))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  scale_color_manual(values=c("Black", "Gold2"))+
  xlab("Difference in Protein expression")+
  ylab("-log2 coefficient of variance")+
  labs(color = "Category:")+ #legend title
  theme_classic()+
  geom_vline(xintercept = -0.1)+
  geom_vline(xintercept = 0.25)+
  #coord_cartesian(xlim=c(-1, 1), ylim=c(0,45))+
  ggtitle("Ovarian cancer: Protein upon chrm gain")
#plot.Protein.TCGA.ovarian.Diff.CVar.gain
# 5x4
cor.test(t.test_TCGA_ovarian.ThreeGroups$Diff.Gain, -log2(abs(t.test_TCGA_ovarian.ThreeGroups$CV.Gain))) 
# corr = 0.236, p<2E-16

# Loss graphs
ggplot(t.test_TCGA_ovarian.ThreeGroups, aes(x=Diff.Loss, y=-log(abs(CV.Loss)), color=-log(abs(CV.Loss))< 0 ))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  scale_color_manual(values=c("black", "Gold2"))+
  xlab("Difference in Protein expression")+
  ylab("coefficient of variance")+
  labs(color = "Category:")+ #legend title
  theme_classic()+
  geom_vline(xintercept = 0.1)+
  geom_vline(xintercept = -0.25)+
  #coord_cartesian(xlim=c(-1, 1), ylim=c(0,45))+
  ggtitle("Ovarian cancer: Protein upon chrm loss")
#plot.Protein.TCGA.ovarian.Diff.pvalue.category.Loss
# 5x4
cor.test(t.test_TCGA_ovarian.ThreeGroups$Diff.Loss, -log2(abs(t.test_TCGA_ovarian.ThreeGroups$CV.Loss))) 
# corr = -0.45, p<2E-16

# conclusion, the larger the difference in expression upon gain/loss, 
# the larger the variance within gain/loss cell population. 

#var var graph
ggplot(t.test_TCGA_ovarian.ThreeGroups, aes(x=-log(abs(CV.Gain)), y=-log(abs(CV.Loss))))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  scale_color_manual(values=c("black", "Gold2"))+
  xlab("Coefficient of variance: Gain")+
  ylab("Coefficient of variance: Loss")+
  labs(color = "Category:")+ #legend title
  theme_classic()+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  #coord_cartesian(xlim=c(-1, 1), ylim=c(0,45))+
  ggtitle("Ovarian cancer: Protein upon chrm gain")
#plot.Protein.TCGA.ovarian.Diff.pvalue.category.Loss
# 5x4

###Bar graph of Protein Scaling/Buffered

#format data so I can plot all 2 groups: Prot gain/loss
dat.m <- melt(t.test_TCGA_ovarian.ThreeGroups, id.vars='Gene', 
              measure.vars=c('Three.RNA.Gain', 'Three.Protein.Gain', 
                             'Three.RNA.Loss', 'Three.Protein.Loss'
                             ))

## Percent of genes per category
ggplot(data= dat.m, aes(x=variable, fill=value)) + 
  geom_bar(position="fill")+ #dodge (next to eachother)
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Protein Difference upon\nchromosome arm gain: per Category")
# 4x4
# plot.TCGA.ovarian.Bargraph.R.P.Gain.Loss


####

## Protein Gain & Loss
ggplot(data= t.test_TCGA_ovarian.ThreeGroups, aes(x=Three.Protein.Gain, fill=Three.Protein.Loss)) + 
  geom_bar(position="fill")+
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Protein difference upon chrm arm gain: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference upon loss:\nCategory")+ #legend title
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Relationship between protein genes \nupon chromosome arm gain or loss")
# 5x4
# plot.TCGA.Protein.Ovarian.Bargraph.ThreeCat.Gain.Loss

setwd(DataFileLocation)
write.csv(t.test_TCGA_ovarian.ThreeGroups, 
          file =paste("Protein_TCGA_Difference_Pvalue_min10points_3cat.csv", sep=','), 
          row.names = TRUE)

t.test_TCGA_ovarian.ThreeGroups<- read.csv("Protein_TCGA_Difference_Pvalue_min10points_3cat.csv")


###### Define function: plot Protein/RNA expression by chrm CN ####
### Step 2: make function to plot protein expression in neutral, gian and loss cells 
## ProtExp.ChrmCN.Ovarian 
## for Protein Data:
ProtExp.ChrmCN.Ovarian<- function(Protein){
  
  TestProt=Protein
  geneData<- subset(CPTAC_Protein_Expression_ovarian, Gene==TestProt)
  TestChrm<- geneData$ChrmNumArm
  colnames(geneData)<-gsub('^', 'TCGA-', colnames(geneData))
  
  #filter cell data and get only aneuploidy scores for chrm arm of test protein location
  if (TestChrm %in% colnames(TCGA_armScores_ovarian)) {
    AneuploidData<- TCGA_armScores_ovarian %>%
      select("Sample", all_of(TestChrm))
    
    GainChrm <- subset(AneuploidData, AneuploidData[,2]== 1)
    NeutralChrm <- subset(AneuploidData, AneuploidData[,2]== 0)
    LossChrm <- subset(AneuploidData, AneuploidData[,2]== -1)
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- (geneData[, (colnames(geneData) %in% GainChrm$Sample)])
    Chrm.di<- (geneData[, (colnames(geneData) %in% NeutralChrm$Sample)])
    Chrm.mono<- (geneData[, (colnames(geneData) %in% LossChrm$Sample)])
    
  ##Make data frame with t-test info about protein expression per arm number category.
  ## Return this data at the end of the function. 
  if (length(Chrm.tri[!is.na(Chrm.tri)])>=3 & 
      #Added if statement so only genes with 2+ values per condition are analyzed
      #if I don't do this, the t-test crashes and I get no values. 
      length(Chrm.di[!is.na(Chrm.di)])>=3 &
      length(Chrm.di[!is.na(Chrm.mono)])>=3) {
    # Trisomy vs disomy
    di.tri<-t.test(Chrm.tri, 
                   Chrm.di, 
                   mu = 0, 
                   alt = "two.sided",
                   conf.level = 0.99) #get p-value
    Diff.Tri.Di<-rowMeans(Chrm.tri, na.rm=TRUE) - rowMeans(Chrm.di, na.rm=TRUE) #get difference
    #Disomy vs monosome
    Mono.Di<-t.test(Chrm.mono, 
                    Chrm.di, 
                    mu = 0, 
                    alt = "two.sided",
                    conf.level = 0.99) #get p-value
    Diff.Mono.Di<-rowMeans(Chrm.mono, na.rm=TRUE) - rowMeans(Chrm.di, na.rm=TRUE)  #get difference
    # Tri vs monosomy
    Tri.Mono<-t.test(Chrm.mono, 
                     Chrm.tri, 
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
    Diff.Tri.Mono<-rowMeans(Chrm.tri, na.rm=TRUE) - rowMeans(Chrm.mono, na.rm=TRUE) #get difference
    
    t.test_TCGA_ovarian<- data.frame(Gene=TestProt,
                                  Pvalue.Gain= di.tri$p.value,
                                  Diff.Gain= Diff.Tri.Di, 
                                  Pvalue.Loss= Mono.Di$p.value,
                                  Diff.Loss= Diff.Mono.Di, 
                                  Pvalue.GainvLoss= Tri.Mono$p.value,
                                  Diff.GainvLoss= Diff.Tri.Mono)
    
    
  } else {
    print("Not enough cells to perform t-test on all conditions") 
  }
    
    ##Make Plot of protein expression per aneuploid category
    Data <- data.frame(expression= t(Chrm.mono), 
                       arm_call="Loss")
    Data<- rbind(Data, 
                 data.frame(expression= t(Chrm.di), 
                            arm_call="Neutral"))
    Data<- rbind(Data, 
                 data.frame(expression= t(Chrm.tri), 
                            arm_call="Gain"))
    
    setwd(DataFileLocation) 
    pdf(paste("Plot.",TestProt,".Protein.Expression.per.Gain.Neutral.Loss.pdf", sep=''), width=4, height=4)## Save as PDF
    Plot.testProt<- ggplot(Data,
                           aes(x = as_factor(arm_call), y=Data[,1])) + 
      geom_boxplot(fill= c("Loss"="dodgerblue3", "Neutral"="grey80", "Gain"="red"), outlier.shape = NA) +
      geom_jitter() +
      xlab(paste("Chromosome arm call of ", TestProt," location: Chrm", TestChrm))+
      ylab("Protein Expression per cell line") +
      theme_classic()+
      ggtitle(paste(TestProt," protein expression per cell\n per corresponding chromosome copy number"))
    print(Plot.testProt)
    dev.off()##stop saving PDF
    return(t.test_TCGA_ovarian)
    }
}
RNAExp.ChrmCN.Ovarian<- function(RNA){
  
  TestProt=RNA
  geneData<- subset(RNA_Expression_ovarian, Gene==TestProt)
  TestChrm<- geneData$ChrmNumArm
  
  #filter cell data and get only aneuploidy scores for chrm arm of test protein location
  if (TestChrm %in% colnames(TCGA_armScores_ovarian)) {
    AneuploidData<- TCGA_armScores_ovarian %>%
      select("Sample", all_of(TestChrm))
    
    GainChrm <- subset(AneuploidData, AneuploidData[,2]== 1)
    NeutralChrm <- subset(AneuploidData, AneuploidData[,2]== 0)
    LossChrm <- subset(AneuploidData, AneuploidData[,2]== -1)
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- as.numeric(geneData[, (colnames(geneData) %in% GainChrm$Sample)])
    Chrm.di<- as.numeric(geneData[, (colnames(geneData) %in% NeutralChrm$Sample)])
    Chrm.mono<- as.numeric(geneData[, (colnames(geneData) %in% LossChrm$Sample)])
    
    ##Make data frame with t-test info about protein expression per arm number category.
    ## Return this data at the end of the function. 
    if (length(Chrm.tri[!is.na(Chrm.tri)])>=3 & 
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        length(Chrm.di[!is.na(Chrm.di)])>=3 &
        length(Chrm.di[!is.na(Chrm.mono)])>=3) {
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri, 
                     Chrm.di, 
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri, na.rm=TRUE) - mean(Chrm.di, na.rm=TRUE) #get difference
      #Disomy vs monosome
      Mono.Di<-t.test(Chrm.mono, 
                      Chrm.di, 
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Mono.Di<-mean(Chrm.mono, na.rm=TRUE) - mean(Chrm.di, na.rm=TRUE)  #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono, 
                       Chrm.tri, 
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri, na.rm=TRUE) - mean(Chrm.mono, na.rm=TRUE) #get difference
      
      t.test_TCGA_ovarian<- data.frame(Gene=TestProt,
                                       Pvalue.Gain= di.tri$p.value,
                                       Diff.Gain= Diff.Tri.Di, 
                                       Pvalue.Loss= Mono.Di$p.value,
                                       Diff.Loss= Diff.Mono.Di, 
                                       Pvalue.GainvLoss= Tri.Mono$p.value,
                                       Diff.GainvLoss= Diff.Tri.Mono)
      
      
    } else {
      print("Not enough cells to perform t-test on all conditions") 
    }
    
    ##Make Plot of protein expression per aneuploid category
    Data <- data.frame(expression= (Chrm.mono), 
                       arm_call="Loss")
    Data<- rbind(Data, 
                 data.frame(expression= (Chrm.di), 
                            arm_call="Neutral"))
    Data<- rbind(Data, 
                 data.frame(expression= (Chrm.tri), 
                            arm_call="Gain"))
    
    setwd(DataFileLocation) 
    pdf(paste("Plot.",TestProt,".RNA.Expression.per.Gain.Neutral.Loss.pdf", sep=''), width=4, height=4)## Save as PDF
    Plot.testRNA<- ggplot(Data,
                          aes(x = as_factor(arm_call), y=Data[,1])) + 
      geom_boxplot(fill= c("Loss"="dodgerblue3", "Neutral"="grey80", "Gain"="red"), outlier.shape = NA) +
      geom_jitter() +
      xlab(paste("Chromosome arm call of ", TestProt," location: Chrm", TestChrm))+
      ylab("RNA Expression per cell line") +
      theme_classic()+
      ggtitle(paste(TestProt," RNA expression per cell\n per corresponding chromosome copy number"))
    print(Plot.testRNA)
    dev.off()##stop saving PDF
    return(t.test_TCGA_ovarian)
  }
}

###### Boxplots of specific gene expression per gain/neutral/loss ####
DataFileLocation<-"/Volumes/Schukken_SSD/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison"
## plot specific protein expression per category
## plot and get p-values for Protein expression changes by chromosome CN
ProtExp.ChrmCN.Ovarian("KRAS") #NS gain & loss
ProtExp.ChrmCN.Ovarian("ERBB2") #NS gain, Sig loss
ProtExp.ChrmCN.Ovarian("EGFR") #Sig gain
ProtExp.ChrmCN.Ovarian("TP53") #not enough datapoints, but also: pvalue loss= 0.5712, Diff.loss= 0.1080914
ProtExp.ChrmCN.Ovarian("NF1")  #NS gain & Loss
ProtExp.ChrmCN.Ovarian("BRCA1") #no location data, even then, not enough data
ProtExp.ChrmCN.Ovarian("BRCA2") #no location data, even then, not enough data
ProtExp.ChrmCN.Ovarian("NRAS")  #NS gain & loss
ProtExp.ChrmCN.Ovarian("RAD21")  #NS gain & loss
ProtExp.ChrmCN.Ovarian("PTEN")  # NS gain and loss
ProtExp.ChrmCN.Ovarian("GAB2")  #NS gain and loss
ProtExp.ChrmCN.Ovarian("RHOA")  #NS gain and loss
ProtExp.ChrmCN.Ovarian("CDK4")  #NS gain and loss
ProtExp.ChrmCN.Ovarian("HLA-B")  #NS
ProtExp.ChrmCN.Ovarian("ATM")  #gain 0.00156, loss=0.03 NS
ProtExp.ChrmCN.Ovarian("PPP2R1A")  #NS
ProtExp.ChrmCN.Ovarian("CASP8")  # NS / NS
ProtExp.ChrmCN.Ovarian("MEN1")  # NS / NS
ProtExp.ChrmCN.Ovarian("KIF1A")  # gain sig (0.011), loss NS
ProtExp.ChrmCN.Ovarian("MTOR")  # NS/NS

## plot and get p-values for RNA expression changes by chromosome CN
RNAExp.ChrmCN.Ovarian("KRAS") #Very Significant
X<- RNAExp.ChrmCN.Ovarian("ERBB2") #Sig loss
RNAExp.ChrmCN.Ovarian("TP53") #NS
RNAExp.ChrmCN.Ovarian("EGFR") #Sig gain
RNAExp.ChrmCN.Ovarian("NF1")  #Significant upon loss
RNAExp.ChrmCN.Ovarian("BRCA1") #Significant
RNAExp.ChrmCN.Ovarian("BRCA2") #Significant upon loss
RNAExp.ChrmCN.Ovarian("NRAS")  #Significant
RNAExp.ChrmCN.Ovarian("CDK4")  #Significant
RNAExp.ChrmCN.Ovarian("RAD21")  #Significant
RNAExp.ChrmCN.Ovarian("HLA-B")  #NS
RNAExp.ChrmCN.Ovarian("PPP2R1A")  #Significant
RNAExp.ChrmCN.Ovarian("ATM")  #Significant
RNAExp.ChrmCN.Ovarian("CASP8")  # Sig
RNAExp.ChrmCN.Ovarian("KIF1A")  # sig gain p= 5E-04
RNAExp.ChrmCN.Ovarian("CTNND1")  # sig loss
RNAExp.ChrmCN.Ovarian("MEN1")  # Sig
RNAExp.ChrmCN.Ovarian("NIPBL")  #sig
RNAExp.ChrmCN.Ovarian("MTOR")  #??


## Bargraph of difference upon chromosome gain or loss
## Since so few cells, many onc are excluded due to not having 10+ cells/category
Onco<-c("RAD21", "RHOB", "PPP2R1A", "CTNNB1", "RHOA", "MTOR", 
        "KRAS", "CDK4", "ERBB2", "PTPN11", "MAP2K1", "GAB2", "SMAD4", "EGFR", "NRAS"
           )
Subset.Ovarian.Onco<- data.frame(Gene=character(), 
                                 Diff.Gain= as.numeric())
for (i in 1:length(Onco)){
  Oncogene<-Onco[i]
  X<-ProtExp.ChrmCN.Ovarian(Oncogene)
  Subset.Ovarian.Onco<- rbind(Subset.Ovarian.Onco, 
                              data.frame(Gene=X$Gene,
                                         Diff.Gain= X$Diff.Gain))
}
Subset.Ovarian.Onco<- Subset.Ovarian.Onco[order(Subset.Ovarian.Onco$Diff.Gain),]

ggplot(Subset.Ovarian.Onco, 
       aes(y=Diff.Gain, x=Gene))+
  geom_bar(stat="identity")+
  theme_classic()+
  geom_hline(yintercept=0)+
  geom_hline(yintercept=mean(TCGA_ovarian_Diff.Pvalue$Diff.Gain,na.rm=TRUE), color="black", size =1)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"))+
  coord_cartesian(ylim=c(-.75,0.75))+
  ylab("Difference upon chrm gain")+
  xlab("Gene name")
# plot.barplot.Onco.Diff.Ovarian_v3
# 6x4



## I tried making bargraph of TSG TP53, BRCA1 and BRCA2, but there were not enough datapoints to graph. 
## Since so few cells, many tsg are excluded due to not having 10+ cells/category
## TSG.Bailey.list 
TSG_test<- subset(t.test_TCGA_ovarian.ThreeGroups, Gene %in% TSG.Bailey.list)

TSG<- c("PTEN" , "NF1", "SMAD4", "ARID2", "ATM", "CASP8", "CDKN2A", "CTNND1", "CUL3", 
         "KIF1A", "MEN1", "NIPBL", "PSIP1", "RUNX1", "TP53"
) 

Subset.Ovarian.TSG<- data.frame(Gene=character(), 
                                 Diff.Loss= as.numeric())
for (i in 1:length(TSG)){
  TSGene<-TSG[i]
  X<-ProtExp.ChrmCN.Ovarian(TSGene)
  Subset.Ovarian.TSG<- rbind(Subset.Ovarian.TSG, 
                              data.frame(Gene=X$Gene,
                                         Diff.Loss= X$Diff.Loss))
}
Subset.Ovarian.TSG<- Subset.Ovarian.TSG[order(Subset.Ovarian.TSG$Diff.Loss),]
Subset.Ovarian.TSG$Gene<- factor(Subset.Ovarian.TSG$Gene, levels= Subset.Ovarian.TSG$Gene)


ggplot(Subset.Ovarian.TSG, 
       aes(y=Diff.Loss, x=Gene))+
  geom_bar(stat="identity")+
  theme_classic()+
  geom_hline(yintercept=0)+
  geom_hline(yintercept=mean(t.test_TCGA_ovarian.ThreeGroups$Diff.Loss,na.rm=TRUE), color="black", size =1)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"))+
  coord_cartesian(ylim=c(-.75,0.75))+
  ylab("Difference upon chrm loss")+
  xlab("Gene name")
# plot.barplot.TSG.Diff.Ovarian_v3
# 6x4



##Tumor supressor gene (TSG) and oncogene (OG) list
## go to TSG.OG_CNV_Difference.R
## Download Bailey et al. 2018 supplemental figure 8, upload table 1. example below: 
# Bailey et al. Comprehensive Characterization of Cancer Driver Genes and Mutations, cell, 2018

Oncogene.Bailey.list #get all oncogenes in Bailey et al 
TSG.Bailey.list # get all TSG in Bailey et al 

## Oncogenes
## Bailey
Oncog.Protein.ttest<-subset(t.test_TCGA_ovarian, Protein_Name %in% Oncogene.Bailey.list)
Oncog.Protein.ttest<- Oncog.Protein.ttest[order(Oncog.Protein.ttest$Diff.Di_Mono),]

# 25 out of 39 (64%) genes not differently expressed between chrm gain or loss (protein level)

p.BRAF<-ProtExp.ChrmCN.Ovarian("BRAF") #Sig (Most sig onco protein)
r.BRAF<-RNAExp.ChrmCN.filtered("BRAF") #Sig (tri)
p.MAPK1<-ProtExp.ChrmCN.Ovarian("MAPK1") #SIG
r.MAPK1<-RNAExp.ChrmCN.filtered("MAPK1") #SIG 
p.RRAS2<-ProtExp.ChrmCN.Ovarian("RRAS2") #SIG (Tri)
r.RRAS2<-RNAExp.ChrmCN.filtered("RRAS2") #SIG (Tri)
p.MTOR<-ProtExp.ChrmCN.Ovarian("MTOR") #NS
r.MTOR<-RNAExp.ChrmCN.filtered("MTOR") #SIG!
p.MYC<-ProtExp.ChrmCN.Ovarian("MYC") #NS
r.MYC<-RNAExp.ChrmCN.filtered("MYC") #NS

## TSG
# Bailey


###### Compare protein/RNA dosage compensation: TCGA vs CCLE #####
### Combine data TCGA Difference and Depmap Difference
### plot (& calculate corr) for difference upon gain/loss.
t.test_TCGA_ovarian.ThreeGroups   # ovarian cancer TCGA Protein & RNA

setwd()
CN.Diff.xRNA.yProt.ThreeGroups<- read.csv(
  file ="RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points_3cat.csv", 
  header = TRUE)


### Protein
Diff.pvalue.TCGAovarianx.CCLEy<-merge(x=t.test_TCGA_ovarian.ThreeGroups, 
                                     y=CN.Diff.xRNA.yProt.ThreeGroups, 
                                     by.x="Gene", 
                                     by.y= "RNA_Name") #2991 genes

# Gain graphs
ggplot(Diff.pvalue.TCGAovarianx.CCLEy, aes(x=Diff.Gain, 
                                          y=Protein.Diff.Gain))+
  #geom_point(size=2)+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  xlab("Difference upon gain: TGCA Ovarian")+
  ylab("Difference upon gain: CCLE")+
  theme_classic()+
  geom_hline(yintercept =0)+
  geom_vline(xintercept=0)+
  geom_smooth(method="lm", color="red")+ #linear trendline
  coord_cartesian(xlim=c(-0.2, 0.5), ylim=c(-0.2, 0.5))+
  ggtitle("ovarian cancer: Protein upon chrm gain")
# plot.Protein.TCGA.Ovarian.CCLE.Gain2
# 4x4

cor.test(Diff.pvalue.TCGAovarianx.CCLEy$Diff.Gain, 
         Diff.pvalue.TCGAovarianx.CCLEy$Protein.Diff.Gain, 
         method="pearson")
# p = 2E-15
# corr= 0.216
# n= 1313

# Min 15 samples: P<2E-16, cor=0.286, n=979
# Min 3 samples : P=4E-13, cor=0.103, n= 4936
# CONTROL normalized, p=2.6E-13, cor=0.210, n= 1191
# min 10 samples, all proteins: # p < 2E-16, corr= 0.230,  n= 1481

# Loss graphs Protein
ggplot(Diff.pvalue.TCGAovarianx.CCLEy, aes(x=Diff.Loss, 
                                           y=Protein.Diff.Loss))+
  #geom_point(size=2)+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  xlab("Difference upon loss: TGCA Ovarian")+
  ylab("Difference upon loss: CCLE")+
  theme_classic()+
  geom_hline(yintercept =0)+
  geom_vline(xintercept=0)+
  geom_smooth(method="lm", color="red")+ #linear trendline
  coord_cartesian(xlim=c(-0.5, 0.2), ylim=c(-0.5, 0.2))+
  ggtitle("ovarian cancer: Protein upon chrm loss")
# plot.Protein.TCGA.Ovarian.CCLE.Loss2
# 4x4

cor.test(Diff.pvalue.TCGAovarianx.CCLEy$Diff.Loss, 
         Diff.pvalue.TCGAovarianx.CCLEy$Protein.Diff.Loss, 
         method="pearson")
# p < 2E-16
# corr= 0.232
# n= 1313

# Min 15 samples: P<2E-16, cor=0.269, n=979
# Min 3 samples : P<2E-16, cor=0.186, n= 4936
# CONTROL normalized, p<2E-16, cor=0.258, n= 1191
# min 10 samples, all proteins: # p < 2E-16, corr= 0.234, n= 1481


### RNA 
# Gain graphs RNA
ggplot(Diff.pvalue.TCGAovarianx.CCLEy, aes(x=Diff.Gain.R, 
                                           y=RNA.Diff.Gain))+
  #geom_point(size=2)+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  xlab("Difference upon gain: TGCA Ovarian")+
  ylab("Difference upon gain: CCLE")+
  theme_classic()+
  geom_hline(yintercept =0)+
  geom_vline(xintercept=0)+
  geom_smooth(method="lm", color="red")+ #linear trendline
  coord_cartesian(xlim=c(-0.1, 0.7), ylim=c(-0.1, 0.7))+
  ggtitle("ovarian cancer: RNA upon chrm gain")
# plot.RNA.TCGA.Ovarian.CCLE.Gain2
# 4x4

cor.test(Diff.pvalue.TCGAovarianx.CCLEy$Diff.Gain.R, 
         Diff.pvalue.TCGAovarianx.CCLEy$RNA.Diff.Gain, 
         method="pearson")
# p = 4E-3
# corr= 0.198
# n= 1313

# Loss graphs RNA
ggplot(Diff.pvalue.TCGAovarianx.CCLEy, aes(x=Diff.Loss.R, 
                                           y=RNA.Diff.Loss))+
  #geom_point(size=2)+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  xlab("Difference upon loss: TGCA Ovarian")+
  ylab("Difference upon loss: CCLE")+
  theme_classic()+
  geom_hline(yintercept =0)+
  geom_vline(xintercept =0)+
  geom_smooth(method="lm", color="red")+ #linear trendline
  coord_cartesian(xlim=c(-0.7, 0.2), ylim=c(-0.7, 0.2))+
  ggtitle("ovarian cancer: RNA upon chrm loss")
# plot.RNA.TCGA.Ovarian.CCLE.Loss2
# 4x4

cor.test(Diff.pvalue.TCGAovarianx.CCLEy$Diff.Loss.R, 
         Diff.pvalue.TCGAovarianx.CCLEy$RNA.Diff.Loss, 
         method="pearson")
# p < 2E-16
# corr= 0.282
# n= 1313

###### Ovarian TCGA vs CCLE (CCLE ovary only: failed, no genes)  #####
### Combine data TCGA Difference and Depmap Difference
### plot (& calculate corr) for difference upon gain/loss.
t.test_TCGA_ovarian.ThreeGroups   # ovarian cancer TCGA Protein & RNA
CN.Diff.RNA.Prot_Ovary # CCLE ovary data (3 datapoints each)


# I make a CCLE difference dataset looking at only ovarian cancer cell lines
# but there were no genes present in that dataset also present in the TCGA dataset
# No genes overlap. :(  
# no genes with data (min 3 points) difference in TCGA and CCLE data. 


###### Buffering Factor Boxplots for Ovarian tumors ######
### Find values for ubiquitination, PPI and protein complex membership, 
### and dependency score per protein
### plot factor score per AS Buffering and Scaling group. 
### Get all factors from Protein_buffering_factors_v2.R, or from supplementary file
setwd()
All.Factors.Diff<- read.csv(file= "Protein.AllFactors.csv")
All.Factors.Diff.Ovarian<-All.Factors.Diff[-c(1,3,4, 6:21)]

All.Factors.Diff.Ovarian2<-merge(x=t.test_TCGA_ovarian.ThreeGroups, y=All.Factors.Diff.Ovarian, 
                                 by.x= "Gene", by.y= "Gene_Symbol")

### boxplots of mean factor score/value in Buffered, scaling and antiscaling genes upon either gain and loss

y.lab="NED" #Name of factor you are looking at
testcollumn<- All.Factors.Diff.Ovarian2$Non.exponential.decay.delta #get the collumn of factor you want to look at


## Chrm gain categories
ggplot(All.Factors.Diff.Ovarian2, aes(x=Three.Protein.Gain.x, y=testcollumn, fill=Three.Protein.Gain.x))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("darkorange2", "Gold2", "palegreen4"))+
  theme_classic()+
  ggtitle("Gain")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(-0.2,0.6))+
  xlab("")+
  ylab(y.lab)
# 5x4
# boxplot.mean.factor.Ovarian.Gain.NED

#Chrm loss categories
ggplot(All.Factors.Diff.Ovarian2, aes(x=Three.Protein.Loss.x, y=testcollumn, fill=Three.Protein.Loss.x))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("darkorange2", "Gold2", "palegreen4"))+
  theme_classic()+
  ggtitle("Loss")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(-0.2,0.6))+
  xlab("")+
  ylab(y.lab)
# 5x4
# boxplot.mean.factor.Ovarian.Loss.NED

G.AS<- subset(All.Factors.Diff.Ovarian2, Three.Protein.Gain.x=="Anti-Scaling") #chrm gain, category anti scaling
G.B<- subset(All.Factors.Diff.Ovarian2, Three.Protein.Gain.x=="Buffering") #chrm gain, category buffered
G.S<- subset(All.Factors.Diff.Ovarian2, Three.Protein.Gain.x=="Scaling") #chrm gain, category scaling
L.AS<- subset(All.Factors.Diff.Ovarian2, Three.Protein.Loss.x=="Anti-Scaling")
L.B<- subset(All.Factors.Diff.Ovarian2, Three.Protein.Loss.x=="Buffering")
L.S<- subset(All.Factors.Diff.Ovarian2, Three.Protein.Loss.x=="Scaling")

## Test significance between groups
# change the collumn to the factor you want to test. 
t.test(G.AS$Dependency.Score, G.B$Dependency.Score)
t.test(G.S$Dependency.Score,  G.B$Dependency.Score)
t.test(L.AS$Dependency.Score, L.B$Dependency.Score)
t.test(L.S$Dependency.Score,  L.B$Dependency.Score)

#Aggregation, scale 0-300
# Gain Buffer v AS, NS (almost significant)
# Gain Buffer v Scale, NS
# loss Buffer v AS, NS
# loss Buffer v Scale, NS

#PPI, scale 0-200
# Gain Buffer v AS, NS
# Gain Buffer v Scale, NS
# loss Buffer v AS, NS
# loss Buffer v Scale, 0.01

#Ubiquitination, scale 0-50
# Gain Buffer v AS, p= 0.00024
# Gain Buffer v Scale, NS
# loss Buffer v AS, NS
# loss Buffer v Scale, NS

# CORUM, scale 0-8
# Gain Buffer v AS, NS
# Gain Buffer v Scale, 3E-05
# loss Buffer v AS, NS
# loss Buffer v Scale, 5E-08

# Dependency, -2 to 1
# Gain Buffer v AS, 9E-05
# Gain Buffer v Scale, <2E-16
# loss Buffer v AS, <2E-16
# loss Buffer v Scale, 6E-11

# NED, -0.2 to 0.6
# Gain Buffer v AS, 0.0007761
# Gain Buffer v Scale, NS
# loss Buffer v AS, NS
# loss Buffer v Scale, NS

# RNA variance, -6 to 3
# Gain Buffer v AS, 2E-12
# Gain Buffer v Scale, 2E-14
# loss Buffer v AS, 2E-16
# loss Buffer v Scale, 1E-04


###### Intro: Plot heatmaps of protein difference upon chrm arm gain and loss in TCGA Ovarian samples####
### 210128
## Update- 210201- only use 371 cell lines with Protein and RNA data
## Update 210427- only use genes that have 10+ cells per category, for both RNA and protein
## Update 210502 - only use genes that have 10+ cells/category (~9k genes), 
## Update 210701 - Use TCGA Ovarian cancer models, 10+ cells/category ()

### by Klaske Schukken
### Graphs for Depmap Protein & RNA data: 
### Graph: heatmap of average change in RNA/Protein experssion per chromosome, per chrm gain/loss

library('ggplot2')
library('tidyverse')
library('xlsx')
library(readxl)
library(reshape2)
library('BBmisc')

###### Step 1: Get data for difference per chrm arm heatmaps ####
### Step 1: Get data. 
## Get  Protein  expression data from depmap.org and cell info data from depmap

## Import Protein expression data and uri ben-david lab, Cohen-Sharir nature 2021 chromosome arm data
###
setwd()# ! set working directory to correct location


DataFileLocation= "" ## set file location to where you need it to be!
setwd(DataFileLocation)

# TCGA aneuploidy data/ arm scores per sample from Taylor et al. 2018, supplementary table 2
TCGA_armScores<-read_xlsx("NIHMS958047-supplement-1_ArmCall.xlsx")

# CPTAC proteomics data that corresponds to TCGA samples. 
# from NIH NCI CPTAC
# data for Breat, ovarian, and colorectal data. 
# 111 Ovarian, 181 Ovarian, 84 Colorectal samples
# for Ovarian and ovarian: use iTRAQ data- log transformed and normalized
# for colorectal: use Spectral Counts, normalize to total protein expression and log2 transform data. 

CPTAC_Ovarian<-read_tsv("TCGA_Ovarian_JHU_Proteome.itraq.tsv") # 132 samples, 8600 genes
CPTAC_Ovarian2<-read_tsv("TCGA_Ovarian_PNNL_Proteome.itraq.tsv") # 84 samples, 7483 genes

CPTAC_Protein_Expression<- merge(CPTAC_Ovarian[,1:265], CPTAC_Ovarian2, 
                                 by.x="Gene", by.y="Gene", all=TRUE)

#protein reported as "log ratio" and "log ratio unshared", remove the unshared collumns. repeats. 
CPTAC_Protein_Expression_rename <- CPTAC_Protein_Expression %>% select(-contains("Unshared"))
#now remove "Log Ratio" and/or "Spectral Counts", to just get sample name. 
# also, there are some repeat samples, use only first instance of that sample:
#      so remove .x, so that sample name with match armScore sample name list. 
colnames(CPTAC_Protein_Expression_rename) <- gsub(' Log Ratio', '', colnames(CPTAC_Protein_Expression_rename))
colnames(CPTAC_Protein_Expression_rename) <- gsub('.x', '', colnames(CPTAC_Protein_Expression_rename))
colnames(CPTAC_Protein_Expression_rename) <- gsub('-01A', '-01', colnames(CPTAC_Protein_Expression_rename))


# Filter CPTAC_Protein_Expression for samples with aneuploidy data
SampleColNames<-gsub('TCGA-', '', TCGA_armScores$Sample)
SampleColNames<- append(SampleColNames, c("Gene", "NCBIGeneID", "Authority", "Description", "Organism", "Chromosome", "Locus"))

CPTAC_Protein_Expression_Ovarian<- CPTAC_Protein_Expression_rename[,(names(CPTAC_Protein_Expression_rename) %in% SampleColNames)]
# 99 samples have protein data and aneuploidy arm length data
# 10628 genes

TCGA_armScores_Ovarian<- TCGA_armScores[gsub('TCGA-', '', TCGA_armScores$Sample) %in% colnames(CPTAC_Protein_Expression_Ovarian),]
TCGA_armScores_Ovarian$Sample<- gsub('TCGA-', '', TCGA_armScores_Ovarian$Sample)
# 99 Ovarian cancer samples with chromosome arm data. 


#Add chrm number + arm collumn
CPTAC_Protein_Expression_Ovarian$arm<- gsub("[^a-zA-Z]", "", CPTAC_Protein_Expression_Ovarian$Locus)
CPTAC_Protein_Expression_Ovarian$arm<- gsub("X", "", CPTAC_Protein_Expression_Ovarian$arm)
CPTAC_Protein_Expression_Ovarian$arm<- gsub("qqq", "q", CPTAC_Protein_Expression_Ovarian$arm)
CPTAC_Protein_Expression_Ovarian$arm<- gsub("qq", "q", CPTAC_Protein_Expression_Ovarian$arm)
CPTAC_Protein_Expression_Ovarian$arm<- gsub("ppp", "p", CPTAC_Protein_Expression_Ovarian$arm)
CPTAC_Protein_Expression_Ovarian$arm<- gsub("pp", "p", CPTAC_Protein_Expression_Ovarian$arm)
CPTAC_Protein_Expression_Ovarian$arm<- gsub("pandYp", "p", CPTAC_Protein_Expression_Ovarian$arm)

CPTAC_Protein_Expression_Ovarian$ChrmNumArm<- paste0(CPTAC_Protein_Expression_Ovarian$Chromosome, CPTAC_Protein_Expression_Ovarian$arm) 



### Get RNA expression: 
# TCGA transcriptomic data that corresponds to TCGA samples. 
# from 
# RNA expressed as log2 lowess normalized (cy5/cy3) collapsed by gene symbol
setwd(DataFileLocation)
RNA_Ovarian<-read_tsv("OV.transcriptome__agilentg4502a_07_2_unc_edu_Level_3_data.txt") # 132 samples, 8600 genes
RNA_Ovarian2<-read_tsv("OV.transcriptome__agilentg4502a_07_3_unc_edu_Level_3_data.txt") # 84 samples, 7483 genes
#RNA_Ovarian<- RNA_Ovarian[-1,]
#RNA_Ovarian2<- RNA_Ovarian2[-1,]
colnames(RNA_Ovarian)[1]<-"Gene"
colnames(RNA_Ovarian2)[1]<-"Gene"
RNA_Expression<- merge(RNA_Ovarian, RNA_Ovarian2, 
                       by.x="Gene", by.y="Gene", all=TRUE) 
#length= 598, 17814 genes


#protein reported as "log ratio" and "log ratio unshared", remove the unshared collumns. repeats. 
RNA_Expression_rename <- RNA_Expression %>% select(-contains("Unshared"))
#now remove "Log Ratio" and/or "Spectral Counts", to just get sample name. 
# also, there are some repeat samples, use only first instance of that sample:
#      so remove .x, so that sample name with match armScore sample name list. 
colnames(RNA_Expression_rename) <- gsub('.x', '', colnames(RNA_Expression_rename))
colnames(RNA_Expression_rename) <- gsub('-01A.*', '-01', colnames(RNA_Expression_rename))
colnames(RNA_Expression_rename) <- gsub('-01B.*', '-01', colnames(RNA_Expression_rename))
colnames(RNA_Expression_rename) <- gsub('-01C.*', '-01', colnames(RNA_Expression_rename))
colnames(RNA_Expression_rename) <- gsub('-01D.*', '-01', colnames(RNA_Expression_rename))
colnames(RNA_Expression_rename) <- gsub('-01R.*', '-01', colnames(RNA_Expression_rename))
colnames(RNA_Expression_rename) <- gsub('-02A.*', '-02', colnames(RNA_Expression_rename))
colnames(RNA_Expression_rename) <- gsub('-11A.*', '-11', colnames(RNA_Expression_rename))


# Filter CPTAC_Protein_Expression for samples with aneuploidy data
SampleColNames<-TCGA_armScores$Sample
SampleColNames<- append(SampleColNames, c("Gene", "NCBIGeneID", "Authority", "Description", "Organism", "Chromosome", "Locus"))

RNA_Expression_ovarian<- RNA_Expression_rename[,(names(RNA_Expression_rename) %in% SampleColNames)]
# 168 samples have protein data and aneuploidy arm length data
# 8911 genes

TCGA_armScores_ovarian<- TCGA_armScores[TCGA_armScores$Sample %in% colnames(RNA_Expression_ovarian),]
TCGA_armScores_ovarian$Sample<- TCGA_armScores_ovarian$Sample
# 168 ovarian cancer samples with chromosome arm data. 

GeneInfo <- CPTAC_Protein_Expression_ovarian[,c(1,176)] #gene names and ChrmNumArm

#Add chrm number + arm collumn
RNA_Expression_ovarian<- merge(RNA_Expression_ovarian, GeneInfo, 
                               by.x="Gene", by.y="Gene")
# 7427 genes shared with protein data. lost 1484 genes
# 542 cell lines (544 with gene and chrmNumArm)


### Now filter Protein & RNA expression for genes that have a minimum of 10 cells per
## category (gain, neutral, loss; in protein and RNA): 
## 3253 genes in filtered data set 
## filter (10+ cells/category dataset from: TCGA_ProteinDiff_Ovarian.R)
# t.test_TCGA_ovarian.ThreeGroups
t.test_TCGA_ovarian.ThreeGroups<-read.delim2("TCGA_ovarian.diff.pvalue.RNA.Protein.ThreeGroups.csv", 
                                             dec=".", header = TRUE, sep=",")

# Proteins: get collumns that have genes with 10+ cells/condition. 
# select all collumns that have same ID as in filtered list. 
CPTAC_Protein_Expression_Ovarian_min10Cells<- CPTAC_Protein_Expression_Ovarian[CPTAC_Protein_Expression_Ovarian$Gene %in% t.test_TCGA_ovarian.ThreeGroups$Gene,]

# RNA: get collumns that have genes with 10+ cells/condition. 
# select all collumns that have same ID as in filtered list. 
RNA_Expression_Ovarian_min10Cells<- RNA_Expression_ovarian[RNA_Expression_ovarian$Gene %in% t.test_TCGA_ovarian.ThreeGroups$Gene,]



###### Step 2: Heatmap for PROTEIN Chromosome GAIN and LOSS####
### Step 2: make heatmap function for PROTEIN Chromosome GAIN or LOSS
### Using only 99 Ovarian cancer samples with Protein data
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
Protein.Gain.Diff.Ovarian<-data.frame(ChrmArm=as.character(), Difference=as.numeric(), Chrm.Gained=as.character())
Protein.Loss.Diff.Ovarian<-data.frame(ChrmArm=as.character(), Difference=as.numeric(), Chrm.Lost=as.character())


###Run data below
##Repeat task for all chromosome arm gains:
for (i in 1:length(AllChrmArms$Chrm)) {
  TestChrm<-AllChrmArms$Chrm[i]
  TestArm<-AllChrmArms$Arm[i]
  TestChrmArm<-paste0(TestChrm, TestArm)
  
  #subset 1: Find cells with gain or no gain of chromosome X
  
  ## Get samples ID with gain, loss or neutral for test chrm arm
  AneuploidData<- TCGA_armScores_Ovarian %>%
    select("Sample", all_of(TestChrmArm))
  
  GainChrm <- subset(AneuploidData, AneuploidData[,2]== 1)
  NeutralChrm <- subset(AneuploidData, AneuploidData[,2]== 0)
  LossChrm <- subset(AneuploidData, AneuploidData[,2]== -1)
  
  ## put data into categories based on chrm arm number
  Chrm.tri<- CPTAC_Protein_Expression_Ovarian_min10Cells %>%
    select(ChrmNumArm, Gene, GainChrm$Sample)
  Chrm.di<- CPTAC_Protein_Expression_Ovarian_min10Cells %>%
    select(ChrmNumArm, Gene, NeutralChrm$Sample)
  Chrm.mono<- CPTAC_Protein_Expression_Ovarian_min10Cells %>%
    select(ChrmNumArm, Gene, LossChrm$Sample)
  
  #Substep2: get difference in protein expression between gain & neutral cells
  if (length(Chrm.tri)>12 &
      length(Chrm.di)>12 &   ## Minimum of 10 samples per 
      length(Chrm.mono)>12){
    #Now get list of protein expression in cells trisomic and disomic for each chrm arm
    # for each protein, find mean Tri, mean Di, difference
    Diff.PerProtein<- data.frame(Protein_ID=as.character(), #Set up data.frame with 3 collumn
                                 ChrmArm= as.character(),
                                 Gain_Exp=numeric(), 
                                 Neutral_Exp=numeric(), 
                                 Loss_Exp=numeric(),
                                 Difference.Gain=numeric(), 
                                 Difference.Loss=numeric(),
                                 stringsAsFactors = TRUE) #this is needed to prevent errors
    
    for (w in 1:1509) { #Get difference in exp between diploid & triploid cells, per protein
      Diff.PerProtein<- rbind(Diff.PerProtein, 
                              data.frame(Protein_ID=Chrm.tri$Gene[w],
                                         ChrmArm= Chrm.tri$ChrmNumArm[w], 
                                         Gain_Exp=rowMeans(Chrm.tri[w,3:length(Chrm.tri)], na.rm=TRUE),
                                         Neutral_Exp=rowMeans(Chrm.di[w,3:length(Chrm.di)], na.rm=TRUE),
                                         Loss_Exp=rowMeans(Chrm.mono[w,3:length(Chrm.mono)], na.rm=TRUE),
                                         Difference.Gain=rowMeans(Chrm.tri[w,3:length(Chrm.tri)], na.rm=TRUE) - rowMeans(Chrm.di[w,3:length(Chrm.di)], na.rm=TRUE), 
                                         Difference.Loss=rowMeans(Chrm.mono[w,3:length(Chrm.mono)], na.rm=TRUE) - rowMeans(Chrm.di[w,3:length(Chrm.di)], na.rm=TRUE)  ))
    }
    
    # then link Protein values with location 
    
    ## Substep 3: add chromosome location to each gene
    # make Chrm number & arm categories, and group by chrm num/arm
    
    # Remove bad chromosome locations. ex "mitochondria m"--> not for this analysis
    Diff.PerProtein6<-subset(Diff.PerProtein, ChrmArm != "mitochondriaa")
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "NANA")
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "reservedd")
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "2cen-q")
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "6s") #no idea what "s" arm is...
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "7s")
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "22p") #only 1 gene
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "22s") #?
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "21p") #only 1 gene
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "Yp") #?
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "Yq") #only 1 gene
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "NA") 
    Diff.PerProtein6$ChrmArm<-as.factor(Diff.PerProtein6$ChrmArm) #as factor. 
    
    ##
    # Get mean protein expression difference by chrm num/arm category
    Mean.Gain.Diff<-aggregate( Difference.Gain ~ ChrmArm, Diff.PerProtein6, mean ) 
    Mean.Gain.Diff$Chrm.Gained<- TestChrmArm
    
    Mean.Loss.Diff<-aggregate( Difference.Loss ~ ChrmArm, Diff.PerProtein6, mean ) 
    Mean.Loss.Diff$Chrm.Lost<- TestChrmArm
    
    
    Protein.Gain.Diff.Ovarian<-rbind(Protein.Gain.Diff.Ovarian, Mean.Gain.Diff)
    Protein.Loss.Diff.Ovarian<-rbind(Protein.Loss.Diff.Ovarian, Mean.Loss.Diff)
    print(paste0("Finished with chrom arm ", TestChrmArm))
  }
}

### Now plot the heatmap for chromosome Protein chrm Gain changes. 
length(Protein.Gain.Diff.Ovarian$ChrmArm) 
write.csv(Protein.Gain.Diff.Ovarian, 
          file= "Protein.Gain.Diff.PerChromosome.TCGA.Ovarian.csv")

length(Protein.Loss.Diff.Ovarian$ChrmArm) 
write.csv(Protein.Loss.Diff.Ovarian, 
          file= "Protein.Loss.Diff.PerChromosome.TCGA.Ovarian.csv")
#setwd()#set working directory 
Protein.Gain.Diff.Ovarian<-read.delim2("Protein.Gain.Diff.PerChromosome.TCGA.Ovarian.csv", 
                                     dec=".", header = TRUE, sep=",")
Protein.Loss.Diff.Ovarian<-read.delim2("Protein.Loss.Diff.PerChromosome.TCGA.Ovarian.csv", 
                                     dec=".", header = TRUE, sep=",")


Protein.Gain.Diff.Ovarian$ChrmArm <- factor(Protein.Gain.Diff.Ovarian$ChrmArm, 
                                            levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                       "4p", "4q", "5p", "5q", "6p", "6q",
                                                       "7p", "7q", "8p", "8q", "9p", "9q",
                                                       "10p", "10q", "11p", "11q", "12p", "12q",
                                                       "13p", "13q", "14p", "14q", "15p", "15q",
                                                       "16p", "16q", "17p", "17q", "18p", "18q",
                                                       "19p", "19q", "20p", "20q", "21p", "21q",
                                                       "22p", "22q", "Xp", "Xq"))
Protein.Gain.Diff.Ovarian$Chrm.Gained <- factor(Protein.Gain.Diff.Ovarian$Chrm.Gained, 
                                                levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                           "4p", "4q", "5p", "5q", "6p", "6q",
                                                           "7p", "7q", "8p", "8q", "9p", "9q",
                                                           "10p", "10q", "11p", "11q", "12p", "12q",
                                                           "13p", "13q", "14p", "14q", "15p", "15q",
                                                           "16p", "16q", "17p", "17q", "18p", "18q",
                                                           "19p", "19q", "20p", "20q", "21p", "21q",
                                                           "22p", "22q", "Xp", "Xq"))


#Protein.Gain.Diff.Ovarian<-read.delim2("Protein.Gain.Diff.PerChromosome.Filtered.csv", 
#                                            dec=".", header = TRUE, sep=",")

ggplot(Protein.Gain.Diff.Ovarian, aes(x=ChrmArm, y=Chrm.Gained))+ 
  geom_raster(aes(fill = Difference.Gain), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  xlab("Protein expression per chromosome arm")+
  ylab("Chromosome arm gained")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.6, 0.6)) +
  coord_flip()+
  ggtitle("Difference in protein expresssion upon \n chromosome gain, TCGA Ovarian samples")
# 4x5
# plot.heatmap.Protein.Gain.TCGA.Ovarian
# sky blue1 or dodgerblue3


## Plot Lost: 
Protein.Loss.Diff.Ovarian$ChrmArm <- factor(Protein.Loss.Diff.Ovarian$ChrmArm, 
                                            levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                       "4p", "4q", "5p", "5q", "6p", "6q",
                                                       "7p", "7q", "8p", "8q", "9p", "9q",
                                                       "10p", "10q", "11p", "11q", "12p", "12q",
                                                       "13p", "13q", "14p", "14q", "15p", "15q",
                                                       "16p", "16q", "17p", "17q", "18p", "18q",
                                                       "19p", "19q", "20p", "20q", "21p", "21q",
                                                       "22p", "22q", "Xp", "Xq"))
Protein.Loss.Diff.Ovarian$Chrm.Lost <- factor(Protein.Loss.Diff.Ovarian$Chrm.Lost, 
                                              levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                         "4p", "4q", "5p", "5q", "6p", "6q",
                                                         "7p", "7q", "8p", "8q", "9p", "9q",
                                                         "10p", "10q", "11p", "11q", "12p", "12q",
                                                         "13p", "13q", "14p", "14q", "15p", "15q",
                                                         "16p", "16q", "17p", "17q", "18p", "18q",
                                                         "19p", "19q", "20p", "20q", "21p", "21q",
                                                         "22p", "22q", "Xp", "Xq"))


#Protein.Gain.Diff.Ovarian<-read.delim2("Protein.Gain.Diff.PerChromosome.Filtered.csv", 
#                                            dec=".", header = TRUE, sep=",")

ggplot(Protein.Loss.Diff.Ovarian, aes(x=ChrmArm, y=Chrm.Lost))+ 
  geom_raster(aes(fill = Difference.Loss), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  xlab("Protein expression per chromosome arm")+
  ylab("Chromosome arm lost")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.6, 0.6)) +
  coord_flip()+
  ggtitle("Difference in protein expresssion upon \n chromosome loss, TCGA Ovarian samples")
# 5x4
# plot.heatmap.Protein.Lost.TCGA.Ovarian
# sky blue1 or dodgerblue3

###### Step 2: Heatmap for RNA Chromosome GAIN and LOSS####
### Step 2: make heatmap function for RNA Chromosome GAIN or LOSS
### Using only ** Ovarian cancer samples with RNA data
## Substep: 1) find cells with gain & nogain of chromosome X, arm p/q
## substep: 2) Find difference in protein expression per protein
## substep: 3) find location of each protein, chromosome and arm location
## substep: 4) find mean of difference in Prot exp per chromosome arm--> into datadrame
## substep: 5) repeat for all chromosome arms.
## substep: 6) plot mean change in prot per chromosome, per chromosome gain. 

#Set up data
AllChrmArms <- data.frame(Chrm=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,14,15,16,16,17,17,18,18,19,19,20,20,21,22), 
                          Arm=c("p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q",
                                "q","q","q","p","q","p","q","p","q","p","q","p","q","q","q"))
RNA.Gain.Diff.Ovarian<-data.frame(ChrmArm=as.character(), Difference=as.numeric(), Chrm.Gained=as.character())
RNA.Loss.Diff.Ovarian<-data.frame(ChrmArm=as.character(), Difference=as.numeric(), Chrm.Lost=as.character())


###Run data below
##Repeat task for all chromosome arm gains:
for (i in 1:length(AllChrmArms$Chrm)) {
  TestChrm<-AllChrmArms$Chrm[i]
  TestArm<-AllChrmArms$Arm[i]
  TestChrmArm<-paste0(TestChrm, TestArm)
  
  #subset 1: Find cells with gain or no gain of chromosome X
  
  ## Get samples ID with gain, loss or neutral for test chrm arm
  AneuploidData<- TCGA_armScores_Ovarian %>%
    select("Sample", all_of(TestChrmArm))
  AneuploidData$Sample<- gsub('^', 'TCGA-', AneuploidData$Sample)
  
  GainChrm <- subset(AneuploidData, AneuploidData[,2]== 1)
  NeutralChrm <- subset(AneuploidData, AneuploidData[,2]== 0)
  LossChrm <- subset(AneuploidData, AneuploidData[,2]== -1)
  
  ## put data into categories based on chrm arm number
  Chrm.tri<- RNA_Expression_Ovarian_min10Cells %>%
    select(ChrmNumArm, Gene, GainChrm$Sample)
  Chrm.di<- RNA_Expression_Ovarian_min10Cells %>%
    select(ChrmNumArm, Gene, NeutralChrm$Sample)
  Chrm.mono<- RNA_Expression_Ovarian_min10Cells %>%
    select(ChrmNumArm, Gene, LossChrm$Sample)
  
  #Substep2: get difference in protein expression between gain & neutral cells
  if (length(Chrm.tri)>12 &
      length(Chrm.di)>12 &   ##min 10 samples, plus gene names and chrm location columns
      length(Chrm.mono)>12){
    # make numeric
    Chrm.tri[,3:length(Chrm.tri)]<- data.frame(lapply(Chrm.tri[,3:length(Chrm.tri)],as.numeric))
    Chrm.di[,3:length(Chrm.di)]<- data.frame(lapply(Chrm.di[,3:length(Chrm.di)],as.numeric))
    Chrm.mono[,3:length(Chrm.mono)]<- data.frame(lapply(Chrm.mono[,3:length(Chrm.mono)],as.numeric))
    
    #Now get list of protein expression in cells trisomic and disomic for each chrm arm
    # for each protein, find mean Tri, mean Di, difference
    Diff.PerProtein<- data.frame(Protein_ID=as.character(), #Set up data.frame with 3 collumn
                                 ChrmArm= as.character(),
                                 Gain_Exp=numeric(), 
                                 Neutral_Exp=numeric(), 
                                 Loss_Exp=numeric(),
                                 Difference.Gain=numeric(), 
                                 Difference.Loss=numeric(),
                                 stringsAsFactors = TRUE) #this is needed to prevent errors
    
    for (w in 1:1509) { #Get difference in exp between diploid & triploid cells, per protein
      Diff.PerProtein<- rbind(Diff.PerProtein, 
                              data.frame(Protein_ID=Chrm.tri$Gene[w],
                                         ChrmArm= Chrm.tri$ChrmNumArm[w], 
                                         Gain_Exp=rowMeans(Chrm.tri[w,3:length(Chrm.tri)], na.rm=TRUE),
                                         Neutral_Exp=rowMeans(Chrm.di[w,3:length(Chrm.di)], na.rm=TRUE),
                                         Loss_Exp=rowMeans(Chrm.mono[w,3:length(Chrm.mono)], na.rm=TRUE),
                                         Difference.Gain=rowMeans(Chrm.tri[w,3:length(Chrm.tri)], na.rm=TRUE) - rowMeans(Chrm.di[w,3:length(Chrm.di)], na.rm=TRUE), 
                                         Difference.Loss=rowMeans(Chrm.mono[w,3:length(Chrm.mono)], na.rm=TRUE) - rowMeans(Chrm.di[w,3:length(Chrm.di)], na.rm=TRUE)  ))
    }
    
    # then link Protein values with location 
    
    ## Substep 3: add chromosome location to each gene
    # make Chrm number & arm categories, and group by chrm num/arm
    
    # Remove bad chromosome locations. ex "mitochondria m"--> not for this analysis
    Diff.PerProtein6<-subset(Diff.PerProtein, ChrmArm != "mitochondriaa")
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "NANA")
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "reservedd")
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "2cen-q")
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "6s") #no idea what "s" arm is...
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "7s")
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "22p") #only 1 gene
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "22s") #?
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "21p") #only 1 gene
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "Yp") #?
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "Yq") #only 1 gene
    Diff.PerProtein6<-subset(Diff.PerProtein6, ChrmArm != "NA") 
    Diff.PerProtein6$ChrmArm<-as.factor(Diff.PerProtein6$ChrmArm) #as factor. 
    
    ##
    # Get mean protein expression difference by chrm num/arm category
    Mean.Gain.Diff<-aggregate( Difference.Gain ~ ChrmArm, Diff.PerProtein6, mean ) 
    Mean.Gain.Diff$Chrm.Gained<- TestChrmArm
    
    Mean.Loss.Diff<-aggregate( Difference.Loss ~ ChrmArm, Diff.PerProtein6, mean ) 
    Mean.Loss.Diff$Chrm.Lost<- TestChrmArm
    
    
    RNA.Gain.Diff.Ovarian<-rbind(RNA.Gain.Diff.Ovarian, Mean.Gain.Diff)
    RNA.Loss.Diff.Ovarian<-rbind(RNA.Loss.Diff.Ovarian, Mean.Loss.Diff)
    print(paste0("Finished with chrom arm ", TestChrmArm))
  }
}

### Now plot the heatmap for chromosome RNA chrm Gain changes. 
length(RNA.Gain.Diff.Ovarian$ChrmArm) 
write.csv(RNA.Gain.Diff.Ovarian, 
          file= "RNA.Gain.Diff.PerChromosome.TCGA.Ovarian.csv")

length(RNA.Loss.Diff.Ovarian$ChrmArm) 
write.csv(RNA.Loss.Diff.Ovarian, 
          file= "RNA.Loss.Diff.PerChromosome.TCGA.Ovarian.csv")
#setwd() #set working directory 
RNA.Gain.Diff.Ovarian<-read.delim2("RNA.Gain.Diff.PerChromosome.TCGA.Ovarian.csv", 
                                     dec=".", header = TRUE, sep=",")
RNA.Loss.Diff.Ovarian<-read.delim2("RNA.Loss.Diff.PerChromosome.TCGA.Ovarian.csv", 
                                     dec=".", header = TRUE, sep=",")


RNA.Gain.Diff.Ovarian$ChrmArm <- factor(RNA.Gain.Diff.Ovarian$ChrmArm, 
                                        levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                   "4p", "4q", "5p", "5q", "6p", "6q",
                                                   "7p", "7q", "8p", "8q", "9p", "9q",
                                                   "10p", "10q", "11p", "11q", "12p", "12q",
                                                   "13p", "13q", "14p", "14q", "15p", "15q",
                                                   "16p", "16q", "17p", "17q", "18p", "18q",
                                                   "19p", "19q", "20p", "20q", "21p", "21q",
                                                   "22p", "22q", "Xp", "Xq"))
RNA.Gain.Diff.Ovarian$Chrm.Gained <- factor(RNA.Gain.Diff.Ovarian$Chrm.Gained, 
                                            levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                       "4p", "4q", "5p", "5q", "6p", "6q",
                                                       "7p", "7q", "8p", "8q", "9p", "9q",
                                                       "10p", "10q", "11p", "11q", "12p", "12q",
                                                       "13p", "13q", "14p", "14q", "15p", "15q",
                                                       "16p", "16q", "17p", "17q", "18p", "18q",
                                                       "19p", "19q", "20p", "20q", "21p", "21q",
                                                       "22p", "22q", "Xp", "Xq"))


#RNA.Gain.Diff.Ovarian<-read.delim2("RNA.Gain.Diff.PerChromosome.Filtered.csv", 
#                                            dec=".", header = TRUE, sep=",")

ggplot(RNA.Gain.Diff.Ovarian, aes(x=ChrmArm, y=Chrm.Gained))+ 
  geom_raster(aes(fill = Difference.Gain), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  xlab("RNA expression per chromosome arm")+
  ylab("Chromosome arm gained")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.6, 0.6)) +
  coord_flip()+
  ggtitle("Difference in RNA expresssion upon \n chromosome gain, TCGA Ovarian samples")
# 5x4
# plot.heatmap.RNA.Gain.TCGA.Ovarian
# sky blue1 or dodgerblue3


## Plot Lost: 
RNA.Loss.Diff.Ovarian$ChrmArm <- factor(RNA.Loss.Diff.Ovarian$ChrmArm, 
                                        levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                   "4p", "4q", "5p", "5q", "6p", "6q",
                                                   "7p", "7q", "8p", "8q", "9p", "9q",
                                                   "10p", "10q", "11p", "11q", "12p", "12q",
                                                   "13p", "13q", "14p", "14q", "15p", "15q",
                                                   "16p", "16q", "17p", "17q", "18p", "18q",
                                                   "19p", "19q", "20p", "20q", "21p", "21q",
                                                   "22p", "22q", "Xp", "Xq"))
RNA.Loss.Diff.Ovarian$Chrm.Lost <- factor(RNA.Loss.Diff.Ovarian$Chrm.Lost, 
                                          levels = c("1p", "1q", "2p", "2q", "3p", "3q",
                                                     "4p", "4q", "5p", "5q", "6p", "6q",
                                                     "7p", "7q", "8p", "8q", "9p", "9q",
                                                     "10p", "10q", "11p", "11q", "12p", "12q",
                                                     "13p", "13q", "14p", "14q", "15p", "15q",
                                                     "16p", "16q", "17p", "17q", "18p", "18q",
                                                     "19p", "19q", "20p", "20q", "21p", "21q",
                                                     "22p", "22q", "Xp", "Xq"))


#RNA.Gain.Diff.Ovarian<-read.delim2("RNA.Gain.Diff.PerChromosome.Filtered.csv", 
#                                            dec=".", header = TRUE, sep=",")

ggplot(RNA.Loss.Diff.Ovarian, aes(x=ChrmArm, y=Chrm.Lost))+ 
  geom_raster(aes(fill = Difference.Loss), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  xlab("RNA expression per chromosome arm")+
  ylab("Chromosome arm lost")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.6, 0.6)) +
  coord_flip()+
  ggtitle("Difference in RNA expresssion upon \n chromosome loss, TCGA Ovarian samples")
# 5x4
# plot.heatmap.RNA.Lost.TCGA.Ovarian
# sky blue1 or dodgerblue3

###### Step 3: Make heatmap for Protein with difference per chrm arm upon gain/loss of that arm ####
### Step 3: Make heatmap for Protein with diff per chrm arm upon gain/loss of arm
# Substep 1) Get data generated above: change in expression per arm upon each arm gain/loss
# Substep 2) Isolate difference in Arm X per gain/loss of arm X in Ovarian cancer samples
# Substep 3) Combine RNA & Protein data. 
# Substep 4) Heatmap, plot data of RNA and Protein together

### heatmap 1: Gain of chromosome arm 
Protein.Gain.Diff.Ovarian
Protein.Loss.Diff.Ovarian
RNA.Gain.Diff.Ovarian
RNA.Loss.Diff.Ovarian

Gene.Gain.Loss.DiffInArm.ovarian<- data.frame(GAINorLOSS=factor(), #Set up data.frame with 4 collumn
                                              ChrmArm=factor(), 
                                              Difference=numeric(), 
                                              Aneuploid.Chrm=factor(), 
                                              stringsAsFactors = TRUE) #this is needed to prevent errors

#Now add add Protein data: if gain chrm Xq, change in Chrm Xq protein expression
for (i in 1:length(Protein.Gain.Diff.Ovarian$ChrmArm)) {
  if (Protein.Gain.Diff.Ovarian$ChrmArm[i]==Protein.Gain.Diff.Ovarian$Chrm.Gained[i]) {
    Gene.Gain.Loss.DiffInArm.ovarian<- rbind(Gene.Gain.Loss.DiffInArm.ovarian, 
                                             data.frame(
                                               GAINorLOSS="Protein Gain",
                                               ChrmArm=Protein.Gain.Diff.Ovarian$ChrmArm[i],
                                               Difference= Protein.Gain.Diff.Ovarian$Difference.Gain[i],
                                               Aneuploid.Chrm= Protein.Gain.Diff.Ovarian$Chrm.Gained[i]
                                             ))
  }
}

for (i in 1:length(Protein.Loss.Diff.Ovarian$ChrmArm)) {
  if (Protein.Loss.Diff.Ovarian$ChrmArm[i]==Protein.Loss.Diff.Ovarian$Chrm.Lost[i]) {
    Gene.Gain.Loss.DiffInArm.ovarian<- rbind(Gene.Gain.Loss.DiffInArm.ovarian, 
                                             data.frame(
                                               GAINorLOSS="Protein Loss",
                                               ChrmArm=Protein.Loss.Diff.Ovarian$ChrmArm[i],
                                               Difference= Protein.Loss.Diff.Ovarian$Difference.Loss[i],
                                               Aneuploid.Chrm= Protein.Loss.Diff.Ovarian$Chrm.Lost[i]
                                             ))
  }
}

for (i in 1:length(RNA.Gain.Diff.Ovarian$ChrmArm)) {
  if (RNA.Gain.Diff.Ovarian$ChrmArm[i]==RNA.Gain.Diff.Ovarian$Chrm.Gained[i]) {
    Gene.Gain.Loss.DiffInArm.ovarian<- rbind(Gene.Gain.Loss.DiffInArm.ovarian, 
                                             data.frame(
                                               GAINorLOSS="RNA Gain",
                                               ChrmArm=RNA.Gain.Diff.Ovarian$ChrmArm[i],
                                               Difference= RNA.Gain.Diff.Ovarian$Difference.Gain[i],
                                               Aneuploid.Chrm= RNA.Gain.Diff.Ovarian$Chrm.Gained[i]
                                             ))
  }
}

for (i in 1:length(RNA.Loss.Diff.Ovarian$ChrmArm)) {
  if (RNA.Loss.Diff.Ovarian$ChrmArm[i]==RNA.Loss.Diff.Ovarian$Chrm.Lost[i]) {
    Gene.Gain.Loss.DiffInArm.ovarian<- rbind(Gene.Gain.Loss.DiffInArm.ovarian, 
                                             data.frame(
                                               GAINorLOSS="RNA Loss",
                                               ChrmArm=RNA.Loss.Diff.Ovarian$ChrmArm[i],
                                               Difference= RNA.Loss.Diff.Ovarian$Difference.Loss[i],
                                               Aneuploid.Chrm= RNA.Loss.Diff.Ovarian$Chrm.Lost[i]
                                             ))
  }
}

###Now plot Chromosome arm gain changes in RNA& Protein expression. heatmap. 
Gene.Gain.Loss.DiffInArm.ovarian$GAINorLOSS <- factor(Gene.Gain.Loss.DiffInArm.ovarian$GAINorLOSS, 
                                                      levels = c("Protein Gain", "RNA Gain", "Protein Loss", "RNA Loss"))#plot RNA first

ggplot(Gene.Gain.Loss.DiffInArm.ovarian, aes(x=GAINorLOSS, y=Aneuploid.Chrm))+ # size 7x3?
  geom_raster(aes(fill = Difference), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  xlab("Difference in expression")+
  ylab("Chromosome arm gained/lost")+
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  theme_classic()+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.6, 0.6)) +
  coord_flip()+
  ggtitle("Difference in ovarian sample protein expresssion \n upon chromosome gain or loss")
# 8x3
# Plot.Protein.RNA.Gain.Loss.Ovarian.TCGA.Difference
