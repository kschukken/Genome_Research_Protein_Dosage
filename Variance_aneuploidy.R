###### Calculating gene neutral ploidy coefficient variance (RNA and Protein) #####
## Variance aneuploidy.R
## Protien/RNA Variance: corrolate with aneuploidy-induce expression difference? 
## By Klaske Schukken
## 210302
## Update 210504: only do genes with minimum 10 cells per category. 

###
# I will calculate Protein (and/or RNA) neutral-ploidy coefficient of variance
# coefficient of variance= standard deviation of protein (and/or RNA), divided by mean expression, 
#    in only cells with "neutral ploidy" for the chrm the gene is located on
# (cells not-aneuploid for arm of chromosome gene is located on)
# then corrolate variance score with expression difference upon chrm gain and/or loss.

# Final result: correlation coefficient and p-value for difference vs. variance scores
# also plot difference vs variance scores. 
# Expect: low variance == low difference upon aneuploidy (buffered)
#         high varaince == high difference upon aneuploidy (Scaling and anti scaling)

## Protein expression data
## https://www.cell.com/cell/fulltext/S0092-8674(19)31385-6#secsectitle0190
## Nusinow et al. Cell, 2020, Quantitative proteomics of the cancer cell line encyclopedia

## chromosome arm data : Cohen-Sharir (2021) Nature paper 
## RNA expression data:  CCLE depmap data.

library('ggplot2')
library('tidyverse')
library('xlsx')
library('readxl')
library('reshape2')
library('BBmisc')

### Steps: 
# 1) Get filtered data (RNA & Protein & aneuploidy)
#       only cells that have aneuploudy, RNA and Protein data
#       get data about location of each gene
#       get data about difference upon gain or loss, per gene
# 2) Per gene: calculate varaince: 
#       i) find gene location (chrm & arm)
#       ii) get cells neutral-ploidy for that arm
#       iii) calculate coefficient of variance in neutral-ploidy cells, per gene
# 3) Correlate variance with difference per gene
#       for loss and gain
#       get p-value and correlation coefficient 
# 4) Plot difference (loss or gain) vs. variance (within diploid)


###### Step 1: get data #####
### Get filtered data (RNA & Protein & aneuploidy)
#       only cells that have aneuploudy, RNA and Protein data
#       get data about location of each gene? 
#       get data about difference upon chrm gain or diploid, per gene? 

#### Protein expression data 
# filtered to only get data from cells that have RNA and Protein expression
## Protein info: 12755 genes
setwd("/Documents")# ! Set working directory
Protein_Info3<-read_csv("Protein_location_info.csv")
Protein_ProID<-read_csv("Protein_ID_info.csv")

#Get updated Protein info with only 371 cell lines (cell lines also have RNA Data)
# protein data 12755 Proteins
Protein.Expression.filtered<-read.delim2("Protein_Expression_filtered.csv", 
                                         dec=",", header = TRUE, sep=";")


### RNA expression data 
# filtered to only get data from cells that have RNA and Protein expression
# Import RNA expression data (no aneuploidy data). CCLE depmap data.
# 371 cell lines
# 19 144 RNAs

# RNA Info. names, chrm location, etc. 
# 40 674 RNAs info 
RNA_Info<-read_excel("HGNC_Protein_Info_edited.xlsx", sheet ="HGNC names")
RNA_Info2<-RNA_Info[,c(2,3,11,12,22,14)]

# Get updated RNA expression with only 371 cell lines (cell lines also have RNA Data)
# 19 144 RNAs expression
RNA.Expression.filtered<-read.delim2("RNA_Expression_filtered.csv", 
                                     dec=",", header = TRUE, sep=";")


####get cell line aneuploidy data:
aneuploid<- read.csv("arm_calls_data_bendavid.csv", header=TRUE)


# Protein / RNA expression difference upon chrm arm gain/loss, per gene 
# from Protein_RNA_expression.PerCell_v2.R 
# or from supplementary data, difference data
CN.Diff.xRNA.yProt.ThreeGroups <- read.csv("RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points_3cat.csv", 
                                            header=TRUE) #9,414 genes

# Note: loss of about 2k genes (protein ID)-- 
# due to needing both RNA and Protein expression/gene
# and due to needing 3+ cell line with chrm gain/loss for chrm arm of gene location.

## Notes! : 
# VarMean is Variance/Mean, also known as coefficient of variance
# Var.Mean.noAn = is coefficient of variance taken from only cells with neutral ploidy for the gene's chrm location



###### Step 2) Per Protein: calculate coefficient of varaince (neutral ploidy cells only): ####
# 2) Per Protein: calculate varaince: 
#       i) find gene location (chrm & arm)
#       ii) get cells diploid for that arm
#       iii) calculate variance (&SD?) of expression in diploid cells, per gene
# var(x) # x= data to get variance of
# sd (x) # x= data to get standard deviation of. 
# var(x)/mean(x) = coeffiecient of variance
# sapply(data, var) # do var() to each collumn? of data

### Protein variance: 

### Using only 371 cell lines that have both RNA and Protein data (only 367 also with aneuploidy data)
## Substep: 1) find cells with no gain/loss of chromosome X, arm p/q 
## substep: 2) find location of each protein, chromosome and arm location
## substep: 3) find variance, coefficient of varance, & st. dev. in protein expression per protein, in diploid cells
## substep: 4) add variance, coefficient of varance & st. dev. per protein into a dataframe
## substep: 5) repeat for all chromosome arms.
## substep: 6) plot coefficient of varance vs st. dev. of all proteins

#Set up data
AllChrmArms<- data.frame(Chrm=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,14,15,16,16,17,17,18,18,19,19,20,20,21,22), 
                         Arm=c("p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q",
                               "q","q","q","p","q","p","q","p","q","p","q","p","q","q","q"))

## make a list of protein info per protein_ID, so I can sort by prot location later.
# only need to do this once, then use data in for loop
Protein.colnames<-data.frame(Protein_ID=colnames(Protein.Expression.filtered)) #1 col, 12757 genes
# add more protein info per protein. Add Uniprot_Acc and Gene_Symbol
Protein.onChrm1<-merge(x= Protein.colnames, y= Protein_ProID, 
                       by.x="Protein_ID", by.y="Protein_Id", 
                       sort = TRUE)# 7 collumns, 12755 genes

#Find genes with same gene_Symbol
Protein.onChrm2<-merge(x= Protein.onChrm1, y= Protein_Info3, 
                       by.x="Gene_Symbol", by.y="Approved Symbol", 
                       sort = TRUE)# 10 collumns, 12467 genes
#Of genes without matching gene symbol...
No_Gene_Symbol <- anti_join(Protein.onChrm1, Protein_Info3, #finding genes_Symbol without match
                            by = c("Gene_Symbol" = "Approved Symbol"))
#...find the genes with matching uniprot_ID
Protein.onChrm3<-merge(x= No_Gene_Symbol, y= Protein_Info3, 
                       by.x="Uniprot_Acc", by.y="UniProt ID(supplied by UniProt)", 
                       sort = TRUE)# 10 collumns, 253 genes, only no gene-symbol genes with uniprot_ID
# Merge genes with gene symbol and genes with only uniprot ID. 
Protein.onChrm4<-merge(x= Protein.onChrm2, y= Protein.onChrm3, 
                       all=TRUE)# 12 collumns, 12720 genes


# make Chrm number & arm categories, and group by chrm num/arm
Protein.onChrm4$arm <- gsub('[0-9]+', '', Protein.onChrm4$'Chromosome band')
Protein.onChrm4$arm <- gsub('[.]', '', Protein.onChrm4$arm)
Protein.onChrm4$arm <- str_sub(Protein.onChrm4$arm, -1, -1) #start & end on character 1. get only p or q

Protein.onChrm4$Chrm.Num.Arm<- paste0(Protein.onChrm4$Chromosome, Protein.onChrm4$arm)
# 10 collumns, 12720 genes


# Set up empty datat frame to fill in with loop:
Protein.var.sd<-data.frame(Protein_ID=as.character(),
                           Protein_ID.All=as.character(),
                           Variance.noAn=numeric(),
                           Variance.Mean.noAn=numeric(),
                           Variance.All=numeric(),
                           St.Dev.noAn= numeric(),
                           St.Dev.All= numeric(),
                           Chrm.Num.Arm= as.character(),
                           stringsAsFactors = TRUE)

###Run data below
##Repeat task for all chromosome arm gains:
for (i in 1:length(AllChrmArms$Chrm)) {
  TestChrm<-AllChrmArms$Chrm[i]
  TestArm<-AllChrmArms$Arm[i]
  TestChrmArm<-paste0(TestChrm, TestArm)
  
  #subset 1: Find cells with no gain/loss of ChrmArm, and test cells
  testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
  testchrm.percell <- filter(testchrm.percell, arm==TestArm)
  
  Cells.Di.forchrm<- filter(testchrm.percell, arm_call==0) #all cells diploid for testProt chrm
  
  # SubStep 2: get proteins located on TestChrmArm
  # First get protein info, loci, then subset protein expression data for testarm loci only
  # then link Protein values with location 
  
  Protein.onChrm5<-Protein.onChrm4 #this is Protein info only, not Prot expression
  Protein.onChrm5<-subset(Protein.onChrm5, Chrm.Num.Arm == TestChrmArm) #get only Prot on chrm arm
  
  #Now subset RNA expression data with only chromosome arm of interest. 
  Protein.Expression.subset <- Protein.Expression.filtered %>% select(Protein.onChrm5$Protein_ID)
  Protein.Expression.subset$Cell_line<- Protein.Expression.filtered$Cell_line #add back cell info
  
  
  #Substep 3: get difference in protein expression between diploid & test cells
  #Now get list of protein expression in cells disomic & test cells for test ChrmArm
  Protein.Di.Cells<-merge(y= Protein.Expression.subset, x= Cells.Di.forchrm, 
                          by.y="Cell_line", by.x="DepMap_ID", 
                          sort = TRUE)# 368 cells, 12757 proteins 
  
  # for each protein, find mean gain, mean Di, difference
  Var.PerProtein <- data.frame(Protein_ID=as.character(),
                               Protein_ID.All=as.character(),
                               Variance.noAn=numeric(),
                               Variance.Mean.noAn=numeric(),
                               Variance.All=numeric(),
                               St.Dev.noAn= numeric(),
                               St.Dev.All= numeric(),
                               Chrm.Num.Arm= as.character(),
                               stringsAsFactors = TRUE) #this is needed to prevent errors
  
    for (t in 8:length(Protein.Di.Cells)) { 
      #Get difference in exp between cell line x and mean no gain/loss cells cells, per protein
      #protein express subset does not have 7 extra info columns, so add -7 to get same Protein_ID
      Var.PerProtein<- rbind(Var.PerProtein, 
                              data.frame(Protein_ID=colnames(Protein.Di.Cells[t]),
                                         Protein_ID.All=colnames(Protein.Expression.subset[t-7]),
                                         Variance.noAn=var(Protein.Di.Cells[,t], na.rm=TRUE), 
                                         Variance.Mean.noAn=var(Protein.Di.Cells[,t], na.rm=TRUE)/abs(mean(Protein.Di.Cells[,t], na.rm=TRUE)),
                                         Variance.All=var(Protein.Expression.subset[,t-7], na.rm=TRUE),
                                         St.Dev.noAn= sd(Protein.Di.Cells[,t], na.rm=TRUE),
                                         St.Dev.All= sd(Protein.Expression.subset[,t-7], na.rm=TRUE),
                                         Chrm.Num.Arm= TestChrmArm))
  
  }
  ##
  # Get mean protein expression difference by chrm num/arm category
  #Variance.Mean.noAn= variance divided by mean expression on non-aneuploid chromosomes
  Protein.var.sd<-rbind(Protein.var.sd, Var.PerProtein)
  print(paste0("finished with chrm arm ", TestChrmArm))
}


### Now plot the heatmap for chromosome Protein chrm Gain changes. 
length(Protein.var.sd$Protein_ID) #12 149 genes, 502 X, 42 Y, 15 other, 12 missing?
write.csv(Protein.var.sd, 
          file= "Protein.Variance.SD.VarMean.csv")
#Protein.var.sd<-read.delim2("Protein.Variance.SD.VarMean.csv", 
#                                     dec=".", header = TRUE, sep=",")


ggplot(Protein.var.sd, aes(x=Variance.noAn, y=Variance.Mean.noAn))+ #5x5 figure
  geom_point()+
  xlab("Protein Varance \n(cells without gene loci gain/loss)")+
  ylab("Protein Coefficient of Variance \n(cells without gene loci gain/loss)")+
  theme_classic()+
  ggtitle("Protein variance and standard deviation \nwithin cells without gene loci aneuploidy")


###### Step 2.5) Per RNA: calculate coefficient of varaince (neutral ploidy cells only):  ####
# step 2.5) Per RNA: calculate varaince: 
#       i) find gene location (chrm & arm)
#       ii) get cells diploid for that arm
#       iii) calculate variance (&SD?) of expression in diploid cells, per gene

### RNA variance: 
## Step 2.5: get Variance and SD for RNA, in diploid-only and in all cells 
## Substep: 1) find cells with no gain/loss of chromosome X, arm p/q 
## substep: 2) find location of each protein, chromosome and arm location
## substep: 3) find variance & st. dev. in protein expression per protein, in diploid cells
## substep: 4) add variance & st. dev. per protein into a dataframe
## substep: 5) repeat for all chromosome arms.
## substep: 6) plot variance vs st. dev. of all proteins

#Set up data
AllChrmArms<- data.frame(Chrm=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,14,15,16,16,17,17,18,18,19,19,20,20,21,22), 
                         Arm=c("p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q",
                               "q","q","q","p","q","p","q","p","q","p","q","p","q","q","q"))

#Before running, make RNA data for all RNA Expression data 
#with chrom arm info and RNA_ID equal to colnames on RNA expression data. 
#Make RNA expression names into data frame with ID and name
RNA.colnames<-data.frame(RNA_ID=colnames(RNA.Expression.filtered))
RNA.colnames$RNA_Name<-sub("[..].*", "", as.character(RNA.colnames$RNA_ID))
RNA.colnames$RNA_num<-sub(".*[.]{2}", "", as.character(RNA.colnames$RNA_ID))
RNA.colnames$RNA_num<-sub("[.]", "", as.character(RNA.colnames$RNA_num))

#add more RNA info per RNA. Add Uniprot_Acc and Gene_Symbol
RNA.onChrm1<-merge(x=RNA.colnames, y= RNA_Info2, 
                   by.x="RNA_Name", by.y="Approved Symbol", 
                   sort = TRUE)# 8 collumns, 18591 genes
#Of genes without matching gene symbol...
No_Gene_Name <- anti_join(RNA.colnames, RNA_Info2, #finding genes_Symbol without match
                          by = c("RNA_Name" = "Approved Symbol"))
#...find the genes with matching uniprot_ID
RNA.onChrm2<-merge(x= No_Gene_Name, y= RNA_Info2, 
                   by.x="RNA_num", by.y="Entrez Gene ID", 
                   sort = TRUE)# 462 genes
# Merge genes with gene symbol and genes with only uniprot ID. 
RNA.onChrm3<-merge(x= RNA.onChrm1, y= RNA.onChrm2, 
                   all=TRUE)# 9 collumns, 19053 genes


# make Chrm number & arm categories, and group by chrm num/arm
RNA.onChrm3$arm <- gsub('[0-9]+', '', RNA.onChrm3$'Chromosome band')
RNA.onChrm3$arm <- gsub('[.]', '', RNA.onChrm3$arm)
RNA.onChrm3$arm <- str_sub(RNA.onChrm3$arm, -1, -1) #start & end on character 1. get only p or q

RNA.onChrm3$Chrm.Num.Arm<- paste0(RNA.onChrm3$Chromosome, RNA.onChrm3$arm)
RNA.onChrm3$Chrm.Num.Arm<-as.factor(RNA.onChrm3$Chrm.Num.Arm) #19 053 genes 

#RNA.onChrm3 is used below in loop to subsection RNA_ID by chrm arm


# Set up empty data frame to fill in with loop:
RNA.var.sd<-data.frame(RNA_ID=as.character(),
                           RNA_ID.All=as.character(),
                           Variance.noAn=numeric(),
                           Variance.Mean.noAn=numeric(),
                           Variance.All=numeric(),
                           St.Dev.noAn= numeric(),
                           St.Dev.All= numeric(),
                           Chrm.Num.Arm= as.character(),
                           stringsAsFactors = TRUE)

###Run data below
##Repeat task for all chromosome arm gains:
for (i in 1:length(AllChrmArms$Chrm)) {
  TestChrm<-AllChrmArms$Chrm[i]
  TestArm<-AllChrmArms$Arm[i]
  TestChrmArm<-paste0(TestChrm, TestArm)
  
  #subset 1: Find cells with no gain/loss of ChrmArm, and test cells
  testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
  testchrm.percell <- filter(testchrm.percell, arm==TestArm)
  
  Cells.Di.forchrm<- filter(testchrm.percell, arm_call==0) #all cells diploid for testProt chrm
  
  # SubStep 2: get RNAs located on TestChrmArm
  # First get RNA info, loci, then subset RNA expression data for testarm loci only
  # then link RNA values with location 
  
  RNA.onChrm4<-RNA.onChrm3 #this is RNA info only, not RNA expression
  RNA.onChrm4<-subset(RNA.onChrm4, Chrm.Num.Arm == TestChrmArm) #get only RNA on chrm arm
  
  #Now subset RNA expression data with only chromosome arm of interest. 
  RNA.Expression.subset <- RNA.Expression.filtered %>% select(RNA.onChrm4$RNA_ID)
  RNA.Expression.subset$Cell_line<- RNA.Expression.filtered$Cell_line #add back cell info
  
  
  #Substep 3: get difference in RNA expression between diploid & test cells
  #Now get list of RNA expression in cells disomic & test cells for test ChrmArm
  RNA.Di.Cells<-merge(y= RNA.Expression.subset, x= Cells.Di.forchrm, 
                          by.y="Cell_line", by.x="DepMap_ID", 
                          sort = TRUE)# 368 cells, 12757 RNAs 
  
  # for each RNA, find mean gain, mean of neutral variance cells only, difference
  Var.PerRNA <- data.frame(RNA_ID=as.character(),
                               RNA_ID.All=as.character(),
                               Variance.noAn=numeric(),
                               Variance.Mean.noAn=numeric(),
                               Variance.All=numeric(),
                               St.Dev.noAn= numeric(),
                               St.Dev.All= numeric(),
                               Chrm.Num.Arm= as.character(),
                               stringsAsFactors = TRUE) #this is needed to prevent errors
  
  for (t in 8:length(RNA.Di.Cells)) { 
    #Get difference in exp between cell line x and mean no gain/loss cells cells, per RNA
    #RNA express subset does not have 7 extra info columns, so add -7 to get same RNA_ID
    Var.PerRNA<- rbind(Var.PerRNA, 
                           data.frame(RNA_ID=colnames(RNA.Di.Cells[t]),
                                      RNA_ID.All=colnames(RNA.Expression.subset[t-7]),
                                      Variance.noAn=var(RNA.Di.Cells[,t], na.rm=TRUE), 
                                      Variance.Mean.noAn=var(RNA.Di.Cells[,t], na.rm=TRUE)/abs(mean(RNA.Di.Cells[,t], na.rm=TRUE)),
                                      Variance.All=var(RNA.Expression.subset[,t-7], na.rm=TRUE),
                                      St.Dev.noAn= sd(RNA.Di.Cells[,t], na.rm=TRUE),
                                      St.Dev.All= sd(RNA.Expression.subset[,t-7], na.rm=TRUE),
                                      Chrm.Num.Arm= TestChrmArm))
    
  }
  ##
  # Get mean RNA expression difference by chrm num/arm category
  RNA.var.sd<-rbind(RNA.var.sd, Var.PerRNA)
  print(paste0("finished with chrm arm ", TestChrmArm))
}



### Now plot the heatmap for chromosome RNA chrm Gain changes. 
#setwd()#set working directory 

length(RNA.var.sd$RNA_ID) #18 166 genes, 829 genes on X chromosome, 124 Y, 20 other/non-chrm 

RNA.var.sd$Variance.Mean.noAn<-abs(RNA.var.sd$Variance.Mean.noAn)
write.csv(RNA.var.sd,
           file= "RNA.Variance.SD.VarMean.csv")
#RNA.var.sd<-read.delim2("RNA.Variance.SD.VarMean.csv", 
#                                     dec=".", header = TRUE, sep=",")

ggplot(RNA.var.sd, aes(x=log2(St.Dev.noAn), y=log2(Variance.Mean.noAn)))+ #5x5 figure
  geom_point()+
  xlab("log2 RNA St. Dev. \n(in cells without chrm gain or loss)")+
  ylab("Log2 RNA Coefficient of Variance \n(in cells without chrm gain or loss)")+
  theme_classic()+
  ggtitle("RNA variance and standard deviation per gene \nwithin cells without chromosome gain or loss")



###### Step 3) Correlate variance with difference per gene ####
# 3) Correlate variance with difference per gene
#       for loss and gain
#       get p-value and correlation coefficient 

Protein.var.sd$Protein_Name<-as.factor(sub(".*[.]", "", as.character(Protein.var.sd$Protein_ID)))
Protein.var.sd$Protein_Name<-as.factor(sub("_HUMAN", "", as.character(Protein.var.sd$Protein_Name)))

### For Protein Coeff of variance vs Protein difference:
Prot.Gain.Diff.Var<-merge(x= CN.Diff.xRNA.yProt.ThreeGroups, y= Protein.var.sd, 
                       by.x="Protein_ID", by.y="Protein_ID", 
                       sort = TRUE) # 6 collumns, 11235 genes (all with difference score)

#Protein gain pearson correlation: low varance/SD vs difference in expression upon chrm gain/loss
Corr.Prot.Gain.Diff.CV<-cor.test(Prot.Gain.Diff.Var$Protein.Diff.Gain, log2(Prot.Gain.Diff.Var$Variance.Mean.noAn), method="pearson")
# Corr (log2CV to difference)= -0.03, p=0.00198

Corr.Prot.Loss.Diff.CV<-cor.test(Prot.Gain.Diff.Var$Protein.Diff.Loss, log2(Prot.Gain.Diff.Var$Variance.Mean.noAn), method="pearson")
# Corr(log2CV to difference)= 0.10, p-value<2E-16


###
### For RNA chromosome gain:
# CN.Diff.xRNA.yProt.ThreeGroups has 9414 unique RNA_IDs
RNA.var.sd$RNA_Name<-as.factor(sub("[..].*", "", as.character(RNA.var.sd$RNA_ID)))
RNA.var.sd$RNA_num<-as.factor(sub(".*[.]{2}", "", as.character(RNA.var.sd$RNA_ID)))
RNA.var.sd$RNA_num<-as.factor(sub("[.]", "", as.character(RNA.var.sd$RNA_num)))

RNA.Gain.Diff.Var<-merge(x= CN.Diff.xRNA.yProt.ThreeGroups, y=RNA.var.sd, #difference data for 11235 genes
                          by.x="RNA_Name", by.y="RNA_Name", 
                          sort = TRUE) # 6 collumns,  9414 unique RNA names

# Some proteins have same RNA, leading to "Unique" rows with protein data, but same RNA
# then we have repeat RNA rows that mess with difference vs. variable info. 
RNA.Gain.Diff.Var<-distinct(RNA.Gain.Diff.Var) #get only unique RNA data, get rid of duplicates

### RNA log2CV to RNA difference upon aneuploidy
#RNA Gain pearson correlation: low varance/SD vs difference in expression upon chrm gain/loss
Corr.RNA.CV.RNADiff.gain<-cor.test(RNA.Gain.Diff.Var$RNA.Diff.Gain, log2(RNA.Gain.Diff.Var$Variance.Mean.noAn), method="pearson")
# Corr= -0.133, p-value< 2E-16
Corr.RNA.CV.RNADiff.loss<-cor.test(RNA.Gain.Diff.Var$RNA.Diff.Loss, log2(RNA.Gain.Diff.Var$Variance.Mean.noAn), method="pearson")
# Corr= 0.412, p-value< 2E-16


### RNA log2CV to Protein difference upon aneuploidy
#RNA Gain pearson correlation: low varance/SD vs difference in expression upon chrm gain/loss
Corr.RNA.CV.ProteinDiff.gain<-cor.test(RNA.Gain.Diff.Var$Protein.Diff.Gain, log2(RNA.Gain.Diff.Var$Variance.Mean.noAn), method="pearson")
# Corr= -0.000311, p-value= NS
Corr.RNA.CV.ProteinDiff.loss<-cor.test(RNA.Gain.Diff.Var$Protein.Diff.Loss, log2(RNA.Gain.Diff.Var$Variance.Mean.noAn), method="pearson")
# Corr= 0.146, p-value< 2E-16


####
## RNA has distinct low variance group of genes. ~3.5k genes, from 11k total
RNA.LowCV<-subset(RNA.Gain.Diff.Var, log2(Variance.Mean.noAn)>=-1)
write.csv(RNA.LowCV, file= "RNA.lowCV.group.csv")


###### Step 4) Plot difference (gain & loss) vs. variance (within diploid) ####
# 4) Plot difference (gain & loss) vs. variance (within diploid) 
## RNA varaince vs. RNA difference (gain) 
ggplot(RNA.Gain.Diff.Var, aes(x=log2(Variance.Mean.noAn), y=(RNA.Diff.Gain)))+ #5x5 figure
  #geom_density2d(color="black")+
  #geom_point(data = subset(RNA.Gain.Diff.Var, RNA_Name %in% Protein.Loss.HSP$Gene_Symbol), 
  #           aes(x = log2(Variance.Mean.noAn), y=(Difference.y), color="Red"))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  scale_color_manual(values=c("a"="black", "b"=" Red")) +
  scale_fill_manual(values=c("a"="black", "b"="Red")) +
  xlab("log2(RNA coefficient of variance) \n (in cells without chrm gain/loss)")+
  ylab("RNA difference upon chrm gain")+
  theme_classic()+
  #scale_color_manual(values=c("black", "red"))+
  geom_smooth(method="lm", color="Red")+
  coord_cartesian(xlim=c(-6,3), ylim=c(-1,1))+
  ggtitle("RNA Relative Variance within cells without gene loci aneuploidy 
          vs RNA difference in expression upon chromosome gain:
          all genes")
#3x3 for density2d plot
# plot.RNA.Gain.Log2CV.RNADiff

## RNA varaince vs. RNA difference (gain) 
ggplot(RNA.Gain.Diff.Var, aes(x=log2(Variance.Mean.noAn), y=(RNA.Diff.Loss)))+ #5x5 figure
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  scale_color_manual(values=c("a"="black", "b"=" Red")) +
  scale_fill_manual(values=c("a"="black", "b"="Red")) +
  xlab("log2(RNA coefficient of variance) \n (in cells without chrm gain/loss)")+
  ylab("RNA difference upon chrm loss")+
  theme_classic()+
  geom_smooth(method="lm", color="Red")+
  coord_cartesian(xlim=c(-6,3), ylim=c(-1,1))+
  ggtitle("RNA Relative Variance within cells without gene loci aneuploidy 
          vs RNA difference in expression upon chromosome loss:
          all Proteins")
#3x3 for density2d plot
# plot.RNA.Loss.Log2CV.RNADiff

## RNA varaince vs. Protein difference (gain) 
ggplot(RNA.Gain.Diff.Var, aes(x=log2(Variance.Mean.noAn), y=(Protein.Diff.Loss)))+ #5x5 figure
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  scale_color_manual(values=c("a"="black", "b"=" Red")) +
  scale_fill_manual(values=c("a"="black", "b"="Red")) +
  xlab("log2(RNA coefficient of variance) \n (in cells without chrm gain/loss)")+
  ylab("Protein difference upon chrm loss")+
  theme_classic()+
  geom_smooth(method="lm", color="Red")+
  coord_cartesian(xlim=c(-6,3), ylim=c(-1,1))+
  ggtitle("RNA Relative Variance within cells without gene loci aneuploidy 
          vs Protein difference in expression upon chromosome loss:
          all Proteins")
#3x3 
# plot.RNA.Loss.Log2CV.ProtDiff

RNA.Gain.Diff.Var$LowRNA.CV<-log2(RNA.Gain.Diff.Var$Variance.Mean.noAn)< -1
RNA.Gain.Diff.Var$LowRNA.CV<- gsub("TRUE","Low", RNA.Gain.Diff.Var$LowRNA.CV)
RNA.Gain.Diff.Var$LowRNA.CV<- gsub("FALSE","High", RNA.Gain.Diff.Var$LowRNA.CV)
RNA.Gain.Diff.Var$LowRNA.CV<-factor(RNA.Gain.Diff.Var$LowRNA.CV, levels=c("Low", "High"))

# test to see if categorical (High RNA var vs low RNA var): affects RNA/Protein differnce upon aneuploidy
lowRNAvar<-subset(RNA.Gain.Diff.Var, log2(RNA.Gain.Diff.Var$Variance.Mean.noAn)< -1)
highRNAvar<-subset(RNA.Gain.Diff.Var, log2(RNA.Gain.Diff.Var$Variance.Mean.noAn)> -1)
t.test(lowRNAvar$Protein.Diff.Loss, highRNAvar$Protein.Diff.Loss) #p<2E-16 
t.test(lowRNAvar$Protein.Diff.Gain, highRNAvar$Protein.Diff.Gain) #p=0.028
t.test(lowRNAvar$RNA.Diff.Loss, highRNAvar$RNA.Diff.Loss) #p<2E-16 
t.test(lowRNAvar$RNA.Diff.Gain, highRNAvar$RNA.Diff.Gain) #p<2E-16

## RNA varaince vs. Protein difference in expression
ggplot(RNA.Gain.Diff.Var, aes(x=LowRNA.CV, y=Protein.Diff.Loss))+ 
  geom_boxplot(outlier.shape = NA)+
  geom_hline(yintercept = 0)+
  xlab("RNA neutral ploidy variance")+
  ylab("RNA difference upon chrm loss")+
  theme_classic()+
  coord_cartesian( ylim=c(-2,2))+
  ggtitle("RNA neutral ploidy variance  
          vs RNA difference upon chromosome loss")
#4x4
# plot.boxplot.RNALog2CV.RNADiff.Loss

###### Step 7) Save files  #####
##save as csv files. to run through g: profiler. 
write.csv(RNA.lowVar.withDiffData, file= "RNA.lowVar.Diff.chrmgainandloss.csv")

RNA.lowVar<-subset(RNA.var.sd, RNA.var.sd$Variance.Mean.noAn <= 0.1) 

RNA.var.sd #18 052 genes, 829 genes on X chromosome, 124 Y, 20 other, 28 missing?
#this is just RNA variance info, not difference upon aneuploidy info. 
Protein.var.sd #12 149 genes, 502 X, 42 Y, 15 other, 12 missing?

Prot.lowVar.withDiffData<- merge(Prot.Gain.lowVar, Prot.Loss.lowVar, 
          by.x="Protein_ID", by.y= "Protein_ID")#512 Proteins
RNA.lowVar.withDiffData<- merge(RNA.Gain.lowVar, RNA.Loss.lowVar, 
                                by.x="RNA_Name", by.y="RNA_Name") #80 RNAs


# Note: var-diff data is based on gene with both RNA and Protein data, 
# and only genes on autosomes (aneuploidy data) (so no X, Y, mitochondria, or centromere genes)
# so it is a smaller list than all RNA or all Protein data.  
# also 146 RNA_ID genes did not have corresponding RNA IDs in difference data. 