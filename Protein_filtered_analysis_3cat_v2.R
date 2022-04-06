##### Calculate protein difference for all proteins upon chrm arm gain/loss #####
## 200904
## Edited: 201217-- Normalized data
## Edited: 210106-- Normalize by Gene, not Cell line. mean=0 only, not SD. 
##         also Fold change is not just FC, not log2(FC)!
## Edited: 210121-- No More normalization. Data already normalized
##         Figured out I was messing up data during "as.numeric" step. Fixed.
## Edited: 210202-- Use filtered data: only cell lines that have both RNA and Protein data
##                -- also measure changes in protein expression by gain/loss chromosome ARM (not whole chrm)
## Edited: 210302: only use genes with minimum 10 cells per category (Protein & RNA, gain neutral, loss)
#                   also made 3 category split (gain neutral loss) instead of quartiles. 

## Author: Klaske Schukken
## Compare chromosome arm loss/gain with protein expression changes. 
## Note: no protein on chrm doubles or halves as a result of 
##       corresponding chrm gain or loss. 
## Protein expression data
##  https://www.cell.com/cell/fulltext/S0092-8674(19)31385-6#secsectitle0190
## Nusinow et al. Cell, 2020, Quantitative proteomics of the cancer cell line encyclopedia
## and chromosome arm data : Uri Ben-David paper 


## Data already normalized. quote from paper (Nusinow et al. Cell, 2020): 
# For each protein, the filtered peptide TMT values were summed to 
# create non-normalized protein quantifications. To control for 
# differential protein loading within a ten-plex, the summed protein 
# quantities were adjusted to be equal within a ten-plex. Following this, 
# values were log2-transformed, and within each ten-plex the bridge channel 
# protein quantity was subtracted from each sample quantity to create a ratio 
# to the bridge. Bridge samples, now 0, were removed. For each protein, there 
# is some measurement error in the measurement of the bridge sample. To account 
# for this, within each ten-plex, the mean protein expression was centered at 0. 
# Finally, ten-plexes were joined by protein identification to create the complete 
# dataset.

library('ggplot2')
library('tidyverse')
library(xlsx)
library(readxl)
library(reshape2)
library('BBmisc')


##### Step 1: Get data & prep data ####
# Prep data 
# Get protein expression data into proper format
# Only need to do this once, then go to next step and upload data.
dataLocation= c("/Documents") # !! set working directory as datalocation to write files into this location
setwd(dataLocation)#set working directory

Protein_Expression.1<-read_excel("mmc2.xlsx", 
                                 sheet= "Normalized Protein Expression")
Protein_Expression.1<-Protein_Expression.1[,c(1:426)] #delete empty collumns

#Name dataframe with protein_ID and corresponding Protein name, this way we have both gene ID and Protein name
Protein_ProID<-Protein_Expression.1[,c(1:2, 6)] 
Protein_ProID$Protein_Id<- str_replace_all(Protein_ProID$Protein_Id, "[|]", ".")
Protein_ProID$Protein_Id<- str_replace_all(Protein_ProID$Protein_Id, "[-]", ".")


##Step 1: Import Protein expression data and 
# and import Cohen-Sharir et al. 2021 Nature (lab of uri ben-david) chromosome arm data
# Protein expression data (no aneuploidy data). 
## Proteins
## Cell lines 

###Upload Protein expression data. only cell lines with Both RNA and Protein data
## Filtered Protein Data
### This file from Protein_RNA_filtered_CellLine.R
## 210202

Protein.Expression.filtered<-read.delim2("Protein_Expression_filtered.csv", 
                                         dec=".", header = TRUE, sep=",")

# arm call data from  Cohen-Sharir et al. 2021 Nature
aneuploid<- read.csv("arm_calls_data_bendavid.csv", header=TRUE)

#Protein Info. names, chrm location, etc. 
Protein_Info<-read_excel("HGNC_Protein_Info_edited.xlsx", sheet ="HGNC names")
Protein_Info2<-Protein_Info[,c(2,3,11,12,22)]

#Add chrm arm data
Protein_Info3<-Protein_Info2
Protein_Info3$arm <- gsub('[0-9]+', '', Protein_Info3$'Chromosome band') #Find test protein chromosome arm
Protein_Info3$arm <- gsub('[.]', '', Protein_Info3$arm) #also remove period, if needed
Protein_Info3$arm <- str_sub(Protein_Info3$arm, -1, -1) #start & end on last character. get only p or q
#filter cell data and get only aneuploidy scores for chrm arm of test protein location



# Proteins: get collumns that have genes with 10+ cells/condition. 
# collumn names = "sp.Q9NQ94.A1CF_HUMAN" Protein ID format. same as CN.Diff.xRNA.yProt$Protein_ID
# select all collumns that have same ID as in filtered list. 
Protein.Expression.filtered_min10Cells<- Protein.Expression.filtered %>% select(one_of(CN.Diff.xRNA.yProt.ThreeGroups$Protein_ID))
# Now add Cell_lines back in. 
Protein.Expression.filtered_min10Cells$Cell_line<- Protein.Expression.filtered$Cell_line

# RNA: get collumns that have genes with 10+ cells/condition. 
# collumn names = "TSPAN6..7105." RNA ID format. same as CN.Diff.xRNA.yProt.ThreeGroups$RNA_ID
# select all collumns that have same ID as in filtered list. 
# length(unique(CN.Diff.xRNA.yProt.ThreeGroups$RNA_ID)), there are 9094 unique RNAs, 
# (9413 unique proteins, can have multiple protien isoforms per RNA)
RNA.Expression.filtered_min10Cells<- RNA.Expression.filtered %>% select(one_of(CN.Diff.xRNA.yProt.ThreeGroups$RNA_ID))
# Now add Cell_lines back in. 
RNA.Expression.filtered_min10Cells$Cell_line<- RNA.Expression.filtered$Cell_line





##### Step 2: Define function: Add chromosome specific data to Protein Data (ex: chrm 12 GAIN or no) ####
## Step 2: Add chromosome specific data to Protein Data (ex: chrm 12 gain or no)
# Define function to write change in gene expression upon chromosome gain
# TestChrm is chromosome number you want to test. 
# nChrmArm is a number (1 or 2) corresponding to the number of arms gained/lost.
# P-value is corrected for number of test. P-value*12k
## Note: 13, 14, 15, 21, 22 chrm only one sided. 
#dataLocation is setwd () computer location for .csv and .pdf files. 


ChrmGainEffect<-function(nChrmArm=1, 
                         dataLocation=DataLocation)
{
  for (i in 1:length(AllChrmArms$Chrm)) {
    TestChrm<-AllChrmArms$Chrm[i]
    TestArm<-AllChrmArms$Arm[i]
    
    #subset 1: Find cells with gain or no gain of chromosome X
    TestChrm.percell <- filter(aneuploid, chrom==TestChrm) 
    TestChrm.percell <- filter(TestChrm.percell, arm==TestArm)
    
    Cells.Tri.forchrm<- filter(TestChrm.percell, arm_call==nChrmArm) #all cells triploid for testProt chrm
    Cells.Di.forchrm<- filter(TestChrm.percell, arm_call==0) #all cells diploid for testProt chrm
    
    #Substep2: get difference in protein expression between diploid & triploid cells
    #Now get list of protein expression in cells trisomic and disomic for each chrm arm
    Chrm.gain<-merge(y= Protein.Expression.filtered_min10Cells, x= Cells.Tri.forchrm, 
                             by.y="Cell_line", by.x="DepMap_ID", 
                             sort = TRUE)# 368 cells, 12757 proteins 
    Chrm.nogain<-merge(y= Protein.Expression.filtered_min10Cells, x= Cells.Di.forchrm, 
                            by.y="Cell_line", by.x="DepMap_ID", 
                            sort = TRUE)# 368 cells, 12757 proteins 
    
    
  ####
  ## Step 4: get p-values and difference between protein expression in cells with trisomy7 and without
  # length(Protein_Expression) #12757. so [4:12757] are Protein collumns
  #Note: Fold change is not just FOLD Change, NOT log2 FC!!
  
  t.test.Chrm<- data.frame(Protein_ID=as.character(), #Set up data.frame with 3 collumn
                           p.value=numeric(), 
                           Difference=numeric(), 
                           stringsAsFactors = TRUE) #this is needed to prevent errors
  
  for (w in 8:9420)  {
      x<-t.test(Chrm.gain[,w], 
                Chrm.nogain[,w], 
                mu = 0, 
                alt = "two.sided",
                conf.level = 0.99) #get p-value
      Diff<-mean(Chrm.gain[,w], na.rm=TRUE) - mean(Chrm.nogain[,w], na.rm=TRUE) #get difference
      t.test.Chrm<- rbind(t.test.Chrm, 
                          data.frame(Protein_ID=colnames(Chrm.gain[w]), 
                                     p.value=x$p.value, 
                                     Difference=Diff))
    } 
  } #Note: rbind works when you bind 2 data.frames. so need to structure data as dataframe to work
  
  
  ###
  ## Step 5: add chromosome location to each gene
  #add more protein info per protein. Add Uniprot_Acc and Gene_Symbol
  t.test.Chrm2<-merge(x= t.test.Chrm, y= Protein_ProID, 
                      by.x="Protein_ID", by.y="Protein_Id", 
                      sort = TRUE)# 5 collumns, 9410 genes
  
  #Find genes with same gene_Symbol
  t.test.Chrm.loci3<-merge(x= t.test.Chrm2, y= Protein_Info2, 
                           by.x="Gene_Symbol", by.y="Approved Symbol", 
                           sort = TRUE)# 10 collumns, 253 genes, only no gene-symbol genes with uniprot_ID

  #Of genes without matching gene symbol....
  No_Gene_Symbol <- anti_join(t.test.Chrm2, Protein_Info2, #finding genes_Symbol without match
                              by = c("Gene_Symbol" = "Approved Symbol"))
  #...find the genes with matching uniprot_ID
  t.test.Chrm.loci4<-merge(x= No_Gene_Symbol, y= Protein_Info2, 
                           by.x="Uniprot_Acc", by.y="UniProt ID(supplied by UniProt)", 
                           sort = TRUE)# 10 collumns, 12467 genes
  
  # Merge genes with gene symbol and genes with only uniprot ID. 
  t.test.Chrm.loci5<-merge(x= t.test.Chrm.loci3, y= t.test.Chrm.loci4, 
                           all=TRUE)# 11 collumns, 9413 genes
  # 35 genes are lost: they have no gene symbol or uniprot ID identicle to list. 
  # or they do not have a minimum of 10 data points per condition. 
  # Optional: data about each protein, in Excel: 
  #write.xlsx(x= t.test.Chrm.loci3, file= "/Users/user/Documents/Depmap Aneuploidy/Protein/Chrm.gain.deregulated/Protein_Info_Chrm7.xlsx")
  
  
  # make Chrm number & arm categories, and group by chrm num/arm
  t.test.Chrm.loci5$arm <- gsub('[0-9]+', '', t.test.Chrm.loci5$'Chromosome band')
  t.test.Chrm.loci5$arm <- gsub('[.]', '', t.test.Chrm.loci5$arm)
  t.test.Chrm.loci5$arm <- str_sub(t.test.Chrm.loci5$arm, -1, -1) #start & end on character 1. get only p or q
  
  t.test.Chrm.loci5$Chrm.Num.Arm<- paste0(t.test.Chrm.loci5$Chromosome, t.test.Chrm.loci5$arm)
  
  # Remove bad chromosome locations. ex "mitochondria m"--> not for this analysis
  t.test.Chrm.loci6<-subset(t.test.Chrm.loci5, Chrm.Num.Arm != "mitochondriaa")
  t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "NANA")
  t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "reservedd")
  t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "2cen-q")
  t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "6s") #no idea what "s" arm is...
  t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "7s")
  t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "22p") #only 1 gene
  t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "22s") #?
  t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "21p") #only 1 gene
  t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "Yp") #?
  t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "Yq") #only 1 gene
  t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "NA") 
  t.test.Chrm.loci6$Chrm.Num.Arm<-as.factor(t.test.Chrm.loci6$Chrm.Num.Arm) #as factor. 
  
  ##above gets rid of genes with no loci. 
  
  ###
  ## Step 5.2: order by chromosome location
  #Now seperate proteins on chrm of choice vs not on chrm. 
  TestChrmArm=paste0(TestChrm, TestArm)
  t.tst.Chrm.gain<- filter(t.test.Chrm.loci6, Chrm.Num.Arm==TestChrmArm)#Proteins on chrm 7
  
  t.tst.Chrm.nogain<- filter(t.test.Chrm.loci6, Chrm.Num.Arm!=TestChrmArm)#Proteins not on Chrm 7 
  
  #See if difference between protein expression on gained chrm vs not gained chrm is significant
  #mean(t.tst.Chrm.gain$FoldChange)
  #mean(t.tst.Chrm.nogain$FoldChange)
  
  Diff.test<-t.test(t.tst.Chrm.gain$Difference, t.tst.Chrm.nogain$Difference, mu = 0, 
                    alt = "two.sided",
                    conf.level = 0.99) #See if difference between protein expression on gained chrm vs not gained chrm is significant
  
  ###
  ## Step 6: Seperate into 3 categories.  
  ## Save all gene deregulated by chromsome
  ## save test chrm genes deregulated. with quartile number
  ## also plot quartiles if you want. 
  
  ### find top deregulated genes
  t.test.Chrm.loci_ordered<- t.test.Chrm.loci6[order(t.test.Chrm.loci6$Difference),]
  
  t.test.Chrm.loci_ordered$Protein.Gain.3cat<- cut(t.test.Chrm.loci_ordered$Difference,
                                            breaks=c(-Inf,-0.1,0.25,Inf),
                                            include.lowest=TRUE,
                                            labels=c("Anti-Scaling", "Buffering","Scaling"))
  
  setwd(dataLocation)
  
  write.csv(t.test.Chrm.loci_ordered, 
            file =paste0("ChrmArmGain.",TestChrm, TestArm, ".Protein.Diff.all_min10cells.csv"), 
            row.names = TRUE)
  
  ##Now plot quartiles of test chrm proteins, and fold change
  testProteinChange<-t.test.Chrm.loci6[which(t.test.Chrm.loci6$Chrm.Num.Arm == TestChrmArm),]
  
  testProteinChange<-testProteinChange[order(testProteinChange$Difference),]
  
  write.csv(testProteinChange, 
             file =paste0("ChmArmGain.",TestChrm, TestArm,".Protein.Diff.testChrm_min10cells.csv", sep=''), 
             row.names = TRUE)
  
  
  print(paste0("Finished calculating Protein changes for gain of ", TestChrm, TestArm))
}

##### Step 2.5: Define function: Add chromosome specific data to Protein Data (ex: chrm 12 LOSS or no) ####
#similar to ChrmGainEffect, but change text to say "Loss" 
ChrmLossEffect<-function(nChrmArm, 
                         dataLocation="/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Protein_filtered/210503_plot_chrmArm_gain.loss")
{
  for (i in 1:length(AllChrmArms$Chrm)) {
    TestChrm<-AllChrmArms$Chrm[i]
    TestArm<-AllChrmArms$Arm[i]
    
    #subset 1: Find cells with gain or no gain of chromosome X
    TestChrm.percell <- filter(aneuploid, chrom==TestChrm) 
    TestChrm.percell <- filter(TestChrm.percell, arm==TestArm)
    
    Cells.Tri.forchrm<- filter(TestChrm.percell, arm_call==nChrmArm) #all cells triploid for testProt chrm
    Cells.Di.forchrm<- filter(TestChrm.percell, arm_call==0) #all cells diploid for testProt chrm
    
    #Substep2: get difference in protein expression between diploid & triploid cells
    #Now get list of protein expression in cells trisomic and disomic for each chrm arm
    Chrm.gain<-merge(y= Protein.Expression.filtered_min10Cells, x= Cells.Tri.forchrm, 
                     by.y="Cell_line", by.x="DepMap_ID", 
                     sort = TRUE)# 368 cells, 12757 proteins 
    Chrm.nogain<-merge(y= Protein.Expression.filtered_min10Cells, x= Cells.Di.forchrm, 
                       by.y="Cell_line", by.x="DepMap_ID", 
                       sort = TRUE)# 368 cells, 12757 proteins 
    
    
    ####
    ## Step 4: get p-values and difference between protein expression in cells with trisomy7 and without
    #length(Protein_Expression) #12757. so [4:12757] are Protein collumns
    #Note: Fold change is not just FOLD Change, NOT log2 FC!!
    
    t.test.Chrm<- data.frame(Protein_ID=as.character(), #Set up data.frame with 3 collumn
                             p.value=numeric(), 
                             Difference=numeric(), 
                             FoldChange=numeric(), 
                             stringsAsFactors = TRUE) #this is needed to prevent errors
    
    for (i in 8:9420) { #Set up t.test function. also adds Difference and FC. 
      if (length(Chrm.gain[,i][!is.na(Chrm.gain[,i])])>=10 & 
          #Added if statement so only genes with 2+ values per condition are analyzed
          #if I don't do this, the t-test crashes and I get no values. 
          length(Chrm.nogain[,i][!is.na(Chrm.nogain[,i])])>=10) {
        x<-t.test(Chrm.gain[,i], 
                  Chrm.nogain[,i], 
                  mu = 0, 
                  alt = "two.sided",
                  conf.level = 0.99) #get p-value
        Diff<-mean(Chrm.gain[,i], na.rm=TRUE) - mean(Chrm.nogain[,i], na.rm=TRUE) #get difference
        Fold.Change<- mean(Chrm.gain[,i], na.rm=TRUE) / mean(Chrm.nogain[,i], na.rm=TRUE)
        t.test.Chrm<- rbind(t.test.Chrm, 
                            data.frame(Protein_ID=colnames(Chrm.gain[i]), 
                                       p.value=x$p.value, 
                                       Difference=Diff, 
                                       FoldChange=Fold.Change))
      } else {
        t.test.Chrm<-t.test.Chrm 
      }
    } #Note: rbind works when you bind 2 data.frames. so need to structure data as dataframe to work
    
    
    ###
    ## Step 5: add chromosome location to each gene
    #add more protein info per protein. Add Uniprot_Acc and Gene_Symbol

    t.test.Chrm2<-merge(x= t.test.Chrm, y= Protein_ProID, 
                        by.x="Protein_ID", by.y="Protein_Id", 
                        sort = TRUE)# 5 collumns, 9410 genes
    
    #Find genes with same gene_Symbol
    t.test.Chrm.loci3<-merge(x= t.test.Chrm2, y= Protein_Info2, 
                             by.x="Gene_Symbol", by.y="Approved Symbol", 
                             sort = TRUE)# 10 collumns, 7901 genes
    
    #Of genes without matching gene symbol....
    No_Gene_Symbol <- anti_join(x=t.test.Chrm2, y= Protein_Info2,  #finding genes_Symbol without match
                                by = c("Gene_Symbol" = "Approved Symbol"))
    
    #...find the genes with matching uniprot_ID
    t.test.Chrm.loci4<-merge(x= No_Gene_Symbol, y= Protein_Info2, 
                             by.x="Uniprot_Acc", by.y="UniProt ID(supplied by UniProt)", 
                             sort = TRUE)# 10 collumns, 253 genes, only no gene-symbol genes with uniprot_ID
    
    # Merge genes with gene symbol and genes with only uniprot ID. 
    t.test.Chrm.loci5<-merge(x= t.test.Chrm.loci3, y= t.test.Chrm.loci4, 
                             all=TRUE)# 11 collumns, 9410 genes
    #Optional: data about each protein, in Excel: 
    #write.xlsx(x= t.test.Chrm.loci3, file= "/Users/user/Documents/Depmap Aneuploidy/Protein/Chrm.gain.deregulated/Protein_Info_Chrm7.xlsx")
    
    
    # make Chrm number & arm categories, and group by chrm num/arm
    t.test.Chrm.loci5$arm <- gsub('[0-9]+', '', t.test.Chrm.loci5$'Chromosome band')
    t.test.Chrm.loci5$arm <- gsub('[.]', '', t.test.Chrm.loci5$arm)
    t.test.Chrm.loci5$arm <- str_sub(t.test.Chrm.loci5$arm, -1, -1) #start & end on character 1. get only p or q
    
    t.test.Chrm.loci5$Chrm.Num.Arm<- paste0(t.test.Chrm.loci5$Chromosome, t.test.Chrm.loci5$arm)
    
    # Remove bad chromosome locations. ex "mitochondria m"--> not for this analysis
    t.test.Chrm.loci6<-subset(t.test.Chrm.loci5, Chrm.Num.Arm != "mitochondriaa")
    t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "NANA")
    t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "reservedd")
    t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "2cen-q")
    t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "6s") #no idea what "s" arm is...
    t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "7s")
    t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "22p") #only 1 gene
    t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "22s") #?
    t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "21p") #only 1 gene
    t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "Yp") #?
    t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "Yq") #only 1 gene
    t.test.Chrm.loci6<-subset(t.test.Chrm.loci6, Chrm.Num.Arm != "NA") 
    t.test.Chrm.loci6$Chrm.Num.Arm<-as.factor(t.test.Chrm.loci6$Chrm.Num.Arm) #as factor. 
    
    ##above gets rid of genes with no loci. 
    
    ###
    ## Step 5.2: order by chromosome location
    #Now seperate proteins on chrm of choice vs not on chrm. 
    TestChrmArm=paste0(TestChrm, TestArm)
    t.tst.Chrm.gain<- filter(t.test.Chrm.loci6, Chrm.Num.Arm==TestChrmArm)#Proteins on chrm 7
    
    t.tst.Chrm.nogain<- filter(t.test.Chrm.loci6, Chrm.Num.Arm!=TestChrmArm)#Proteins not on Chrm 7 
    
    #See if difference between protein expression on gained chrm vs not gained chrm is significant
    #mean(t.tst.Chrm.gain$FoldChange)
    #mean(t.tst.Chrm.nogain$FoldChange)
    
    Diff.test<-t.test(t.tst.Chrm.gain$Difference, t.tst.Chrm.nogain$Difference, mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #See if difference between protein expression on gained chrm vs not gained chrm is significant
    
    ###
    ## Step 6: Seperate into quartiles. 
    ## Save all gene deregulated by chromsome
    ## save test chrm genes deregulated. with quartile number
    ## also plot quartiles if you want. 
    
    ### find top deregulated genes
    t.test.Chrm.loci_ordered<- t.test.Chrm.loci6[order(t.test.Chrm.loci6$Difference),]
    
    setwd(dataLocation)
    write.csv(t.test.Chrm.loci_ordered, 
              file =paste0("ChrmArmLoss.",TestChrm, TestArm, ".Protein.Diff_min10cells.csv"), 
              row.names = TRUE)
    
    
    ##Now plot quartiles of test chrm proteins, and fold change
    testProteinChange<-t.test.Chrm.loci6[which(t.test.Chrm.loci6$Chrm.Num.Arm == TestChrmArm),]
    
    testProteinChange<-testProteinChange[order(testProteinChange$Difference),]
    
    testProteinChange$Protein.Loss.3cat<- cut(testProteinChange$Difference,
                                              breaks=c(-Inf,-0.25,0.1,Inf),
                                              include.lowest=TRUE,
                                              labels=c("Scaling", "Buffering","Anti-Scaling"))
    
    
    write.csv(testProteinChange, 
               file =paste0("ChmArmLoss.",TestChrm, TestArm,".Protein.Diff.3cat_min10cells.csv"), 
               row.names = TRUE)
    
    
    print(paste0("Finished calculating Protein changes for loss of ", TestChrm, TestArm))
  }
}
##### Step 3: run functions ####
## run functions
## Note: 13, 14, 15, 21, 22 chrm only one sided. 

AllChrmArms<- data.frame(Chrm=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,14,15,16,16,17,17,18,18,19,19,20,20,21,22), 
                         Arm=c("p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q",
                               "q","q","q","p","q","p","q","p","q","p","q","p","q","q","q"))

##Chromosome gain
ChrmGainEffect(nChrmArm = 1)

##Chromosome loss
## Note: 13, 14, 15, 21, 22 chrm only one sided. 
ChrmLossEffect(nChrmArm = -1)
