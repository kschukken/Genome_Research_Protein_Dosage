### Calculate aneuploidy score(s) per cell line
### 200508
### Klaske Schukken
### Goal: take aneuploidy  data per cell line from Uri Ben-David
### and give each cell line an aneuploidy score
### Make a .csv or .txt file of score. 
### than later (in other script) compare to sgRNA dropouts. and 
### find corrolation between aneuploidy & genes? 
### update 210216 - added Gene # Score



library('plot.matrix')# for making colored plots of aneuploid cells 
library('ggplot2')
library('tidyverse')
library('xlsx')

### Step 1: get data ####
#Import aneuploid arm call data from Cohen-Sharir et al. 2021 Nature paper
setwd()
aneuploid<- read.csv("arm_calls_data_bendavid.csv", header=TRUE)

# make dataframe of chrm+arm, bp_arm, and gene_arm
# make this in excel, save as .csv, then upload to R
# bp data estimated from NIH genome viewer
# gene data gotten from HGNC, genes with protein products only. https://www.genenames.org/
bp_arm<- read.csv("bp_per_arm.csv", header=TRUE)


# gene data gotten from HGNC, genes with protein products only. https://www.genenames.org/
#get HGNC gene data. info per gene, including chrm loci and type of gene
HGCN_gene_data<- read.csv("HGNC.csv", header=TRUE) 



### Step 2: calculate aneuploidy score(s) per chrm arm ####
## add score for aneuploid chrm arms. column 8. 
## count number of aneuploid arms per cell. note: only last row has number of aneuploid arms
#There are 39 chromosome arms per cell line, 997 cell lines

# to do: make table of arm score: sum of aneuploid arms per cell line
# then merge by cell line, adding proper names
arm_Score<-data.frame(tapply(abs(aneuploid$arm_call), aneuploid$DepMap_ID, FUN=sum))
arm_Score <- cbind(rownames(arm_Score), data.frame(arm_Score, row.names=NULL))
arm_Score<- arm_Score %>% rename(DepMap_ID='rownames(arm_Score)', 
                                 arm_Score='tapply.abs.aneuploid.arm_call...aneuploid.DepMap_ID..FUN...sum.')

aneuploid<-merge(x= aneuploid, y= arm_Score, 
                 by.x="DepMap_ID", by.y="DepMap_ID", sort = FALSE)
aneuploid <- aneuploid[order(aneuploid$DepMap_ID, aneuploid$chrom),]


head(aneuploid)
nrow(aneuploid)/39


### Now I want to make a score based on basepairs of chrmosome. 
## so an aneuploid chrm 1 is bigger than aneuploid chrm 23
## this will be bp_Score
## and I want to take into account ploidy
## so that diploid + 1p is bigger than octoploid+1p. 
## this will be bp_ploidy_Score

# Step 1: add collumn to aneuploid dataframe:  chrom + arm 
aneuploid$Chrm_Arm<- paste0(aneuploid$chrom, aneuploid$arm) #new collumn with chrm_arm: "1_p" 

# Step 2: make dataframe of chrm+arm, bp_arm, and gene_arm
# make this in excel, save as .csv, then upload to R... probably easier. 
# bp data estimated from NIH genome viewer
# gene data gotten from HGNC, genes with protein products only. https://www.genenames.org/
bp_arm$Chrm_Arm<- gsub("_", "", bp_arm$Chrm_Arm)

# Step 3: combine dataframes by chrm_arm. so aneuploid has row of arm bp
aneuploid<-merge(x= aneuploid, y= bp_arm, 
                             by.x="Chrm_Arm", by.y="Chrm_Arm", sort = FALSE)
aneuploid <- aneuploid[order(aneuploid$DepMap_ID, aneuploid$Chrm_Arm),]

###Step 4: bp_aneuploid: arm_call*bp_arm per row
# sum bp_aneuploid per cell line for bp_Score 
aneuploid <- aneuploid %>% mutate(bp_aneuploid = abs(wmed_CN - round(ploidy)) * bp_arm)
head(aneuploid)

# to do: make table of bp score: sum of aneuploid bp per cell line
# then merge by cell line, adding proper names
bp_Score<-data.frame(tapply(aneuploid$bp_aneuploid, aneuploid$DepMap_ID, FUN=sum))
bp_Score <- cbind(rownames(bp_Score), data.frame(bp_Score, row.names=NULL))
bp_Score<- bp_Score %>% rename(DepMap_ID='rownames(bp_Score)', 
                               bp_Score='tapply.aneuploid.bp_aneuploid..aneuploid.DepMap_ID..FUN...sum.')

aneuploid<-merge(x= aneuploid, y= bp_Score, 
                 by.x="DepMap_ID", by.y="DepMap_ID", sort = FALSE)
aneuploid <- aneuploid[order(aneuploid$DepMap_ID, aneuploid$Chrm_Arm),]

head(aneuploid)

### Step 5: bp_ploidy_aneuploid: 
#( (wmed_CN - round(ploidy)) * bp_arm)/round(ploidy)
# sum bp_ploidy_aneuploid per cell line for bp_ploidy_Score 
aneuploid <- aneuploid %>% mutate(bp_ploidy_aneuploid = abs( (wmed_CN - round(ploidy)) * bp_arm)/round(ploidy))

# to do: make table of bp anueploid score: sum of aneuploid bp (divided by ploidy) per cell line
# then merge by cell line, adding proper names
bp_ploidy_Score<-data.frame(tapply(aneuploid$bp_ploidy_aneuploid, aneuploid$DepMap_ID, FUN=sum))
bp_ploidy_Score <- cbind(rownames(bp_ploidy_Score), data.frame(bp_ploidy_Score, row.names=NULL))
bp_ploidy_Score<- bp_ploidy_Score %>% rename(DepMap_ID='rownames(bp_ploidy_Score)', 
                                                   bp_ploidy_Score='tapply.aneuploid.bp_ploidy_aneuploid..aneuploid.DepMap_ID..FUN...sum.')

aneuploid<-merge(x= aneuploid, y= bp_ploidy_Score, 
                 by.x="DepMap_ID", by.y="DepMap_ID", sort = FALSE)
aneuploid <- aneuploid[order(aneuploid$DepMap_ID, aneuploid$Chrm_Arm),]

head(aneuploid)

###
### Step 6: gene_aneuploid: arm_call*gene_arm per row

#sort data: only genes with protein products, and make chrm arm collumn
ProteinCodingGenes<- HGCN_gene_data %>% filter(Locus.Type== "gene with protein product")

ProteinCodingGenes$arm <- gsub('[0-9]+', '', ProteinCodingGenes$'Chromosome.band')
ProteinCodingGenes$arm <- gsub('[.]', '', ProteinCodingGenes$arm)
ProteinCodingGenes$arm <- str_sub(ProteinCodingGenes$arm, -1, -1) #start & end on character 1. get only p or q

ProteinCodingGenes$Chrm.Num.Arm<- paste0(ProteinCodingGenes$chromosome, ProteinCodingGenes$arm)

chrmList<-c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
        "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", 
        "21", "22") #list all chrm locations I want (get rid of "NA" and "random")
armList<- c("p", "q") # list all chrm arms I want (to get rid of "NA", and "s", and "d")
GenesOnChrm<- subset(ProteinCodingGenes, arm %in% armList) #19012 prot coding genes. filter good arm
GenesOnChrm<- subset(GenesOnChrm, chromosome %in% chrmList) #18048 prot coding genes, filter good chrm
GenesOnChrm<- GenesOnChrm[GenesOnChrm$Chrm.Num.Arm != "14p", ]
GenesOnChrm<- GenesOnChrm[GenesOnChrm$Chrm.Num.Arm != "21p", ] #18043 prot coding genes, filter out 14p and 21p


gene_arm<- data.frame(table(GenesOnChrm$Chrm.Num.Arm)) #make data frame with number of genes per chrm arm
gene_arm<- gene_arm %>% rename(Chrm.arm=Var1, Gene.Num=Freq)

### Good, now I have a list of protein coding genes per chromosome arm
## can now make aneuploidy score based on number of genes. 

#first: add number of genes per chrm arm to data
aneuploid<-merge(x= aneuploid, y= gene_arm, 
                 by.x="Chrm_Arm", by.y="Chrm.arm", sort = FALSE)
aneuploid <- aneuploid[order(aneuploid$DepMap_ID, aneuploid$Chrm_Arm),]

#make new category: aneuploid genes per aneuploid chrm arm. (0 if not aneuploid)
# I am basing this on arm call, so 2 chrm gain is same as 1. but 2 chrm gain usually only happens in 
# polyploid background, so this takes ploidy into account more? 
aneuploid <- aneuploid %>% mutate(gene_aneuploid = abs(arm_call) * Gene.Num)


# make list of aneuploid score per cell, sum of aneuploid genes per Cell ID
Gene_Score_perCell<-data.frame(tapply(aneuploid$gene_aneuploid, aneuploid$DepMap_ID, FUN=sum))
Gene_Score_perCell <- cbind(rownames(Gene_Score_perCell), data.frame(Gene_Score_perCell, row.names=NULL))
Gene_Score_perCell<- Gene_Score_perCell %>% rename(DepMap_ID='rownames(Gene_Score_perCell)', 
                                                   Gene_Score='tapply.aneuploid.gene_aneuploid..aneuploid.DepMap_ID..FUN...sum.')

aneuploid<-merge(x= aneuploid, y= Gene_Score_perCell, 
                 by.x="DepMap_ID", by.y="DepMap_ID", sort = FALSE)
aneuploid <- aneuploid[order(aneuploid$DepMap_ID, aneuploid$Chrm_Arm),]

aneuploid$Gene_ploidy_Score<- aneuploid$Gene_Score / aneuploid$ploidy

### Step 3: Make new dataframe with aneuploidy scores/cell #####
###Make new dataframe from the last row of each cell line,
## with the total aneuploidy score 
## then cut out un-needed collumns

list<- c()
for (i in 1:997){
  list<- c(list, i*39)
}# Make a list of values for row numbers with true aneuploid score
#the last row for each cell line (39th row), (for chrmosome 9q)

Score<- aneuploid[c(list),] #make a list of aneuploid scores/cell (only one row per cell)
head(Score)

#ggplot(Score,aes(ploidy)) +
#  geom_bar() # checked ploidy distribution. mostly diploid & triploid. 
#make new data frame with them
#then, delete unneccessary rows: chrm, arm, wmed_CN, arm_call

#If you want to try categorical data, you can. but then it's harder
#to do statistical analysis

#Score$Category<- "Aneuploid"
#Score$Category[Score$bp_ploidy_Score<250000000]= "near euploid"
#Score$Category#235 near euploid, 762 aneuploid. ~24%
#Score$Strict<-  "Aneuploid"
#Score$Strict[Score$bp_ploidy_Score==0]= "euploid"
#Strict : 34 euploid, 963 aneuploid, ~3% euploid

Score<- Score[ , -which(names(Score) %in% c("wmed_CN", "Chrm_Arm",
                                            "arm_call", "chrom", "arm", 
                                            "bp_arm", "bp_ploidy_aneuploid",
                                            "bp_aneuploid", "X", "Chrm_Arm", 
                                            "Gene.Num", "gene_aneuploid"))]
head(Score)

write.csv(Score,
          "Aneuploidy.Score.perCell.csv", row.names = TRUE)

#or upload xlsx file Aneuploidy.Score.perCell.xlsx 
Score<- read_xlsx("Aneuploidy.Score.perCell.xlsx")
# 6.6% of cell lines tested are euploid. 93% are aneuploid to varying degrees. 
#       of which, 32 (3.2% total) have partial chrm arm gain
# all euploid cells are diploid.
# Score[Score$Gene_ploidy_Score==0,]

### Score data info: 
# ploidy is mean chromosome count
# arm_Score is sum of aneuploid chromosome arms. 
# bp_Score is aneuploid abs(CN of chrm -  round(ploidy)) * length of arm (in bp)
#         This will give you a score if you have only part of chromosome gained
# bp_ploidy_Score is aneuploid bp_Score/ ploidy 
#         so that it accounts for gene dosage imbalance: 
#         ex. tetraploid needs 2x aneuploid chrm for same imbalance. 
# Gene_Score is aneuploid chrm arm * genes on arm
#         only genes with known protein coding genes
#         downloaded data 210217
#         from HGNC. https://www.genenames.org/
# Gene_ploidy_Score is aneuploid chrm arm * genes on arm/ ploidy **!!! This is the one used as cellular aneuploidy score
#         so that it accounts for gene dosage imbalance: 
#         ex. tetraploid needs 2x aneuploid chrm for same imbalance.
