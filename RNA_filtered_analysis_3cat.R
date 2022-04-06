##### Calculate RNA difference for all proteins upon chrm arm gain/loss ####
## Filter RNA expression for cells with protein & aneuploidy data
## filter genes for genes with minimum 10+ cells/category
# Get RNA difference upon chrm gain or loss for all genes (9k)
#     subset: gene on chrm x upon chrm x gain

## 201014 
## edited 201217 -- Normalize RNA expression data per cell line
## 210108 Update: Normalize RNA expression by gene, not cell line
## 200121 Update: No Normalization. RNA expression already normalized. 
## 210202 Update: use filtered cell lines- only cells with both RNA and Protein data
##                also measure changes in protein expression by gain/loss chromosome ARM (not whole chrm)
## 210503 Update: filter genes for 10+ cells/category


## By: Klaske M. Schukken
## Compare  chromosome arm loss/gain with RNA expression changes. 
## and group into Anti Scale, Buffer, Scaling categories
## RNA expression data
## Depmap RNA expression database
## and chromosome arm data : Uri Ben-David lab, Cohen-Sharir 2021 nature paper paper 

## Edited: No Extra Normailization. RNA Data has already been normalized and log2 transformed
## "The RNA data has been log normalized: "Log2 transformed, using a pseudo-count of 1."

# To do: find way to correlate RNA name with chrm location. 
# then find cells with chrm gain/loss and find RNA from that chrm. 
# then sort into 3 categories: scaling, buffering, anti-scaling 
# then corrolate with protein expression categories 


library('ggplot2')
library('tidyverse')
library('xlsx')
library(readxl)
library(reshape2)
library('BBmisc')

##### Step 1: Get data & prep data ####
# Prep data 
#Get RNA expression data into proper format

## Get aneuploidy info
dataLocation= c("/Documents") # !! set working directory as datalocation to write files into this location
setwd(dataLocation)#set working directory
aneuploid<- read.csv("arm_calls_data_bendavid.csv", header=TRUE)

#RNA Info. names, chrm location, etc. 
#setwd()#set working directory to depmap
RNA_Info<-read_excel("HGNC_Protein_Info_edited.xlsx", sheet ="HGNC names")
RNA_Info2<-RNA_Info[,c(2,3,11,12,22,14)]

RNA_Info3<-RNA_Info2
RNA_Info3$arm <- gsub('[0-9]+', '', RNA_Info3$'Chromosome band') #Find test RNA chromosome arm
RNA_Info3$arm <- gsub('[.]', '', RNA_Info3$arm) #also remove period, if needed
RNA_Info3$arm <- str_sub(RNA_Info3$arm, -1, -1) #start & end on last character. get only p or q

#RNA data for only cells that have both RNA and Protein Data available. (371 cell lines)
# from Protein_RNA_filtered_CellLine.R 
RNA.Expression.filtered<-read.delim2("RNA_Expression_filtered.csv", 
                                         dec=".", header = TRUE, sep=",")


### Now filter Protein/RNA expression for genes that have a minimum of 10 cells per 
## category (gain, neutral, loss; in both RNA and protein): 
## 9414 genes in filtered data set 
## filter (10+ cells/category dataset from: Protein_RNA_Corr_min10.Filtered.R)
# CN.Diff.xRNA.yProt.ThreeGroups
# from Protein_RNA_expression.PerCell_v2.R 
CN.Diff.xRNA.yProt.ThreeGroups<- read.csv("RNA.Protein_Mono.Di.Tri_Difference_Pvalue_min10points.csv")

# RNA: get collumns that have genes with 10+ cells/condition. 
#collumn names = "TSPAN6..7105." RNA ID format. same as CN.Diff.xRNA.yProt$RNA_ID
# select all collumns that have same ID as in filtered list. 
# length(unique(CN.Diff.xRNA.yProt$RNA_ID)), there are 9094 unique RNAs, 
# (9413 unique proteins, can have multiple protien isoforms per RNA)
RNA.Expression.filtered_min10Cells<- RNA.Expression.filtered %>% select(one_of(CN.Diff.xRNA.yProt.ThreeGroups$RNA_ID))
# Now add Cell_lines back in. 
RNA.Expression.filtered_min10Cells$Cell_line<- RNA.Expression.filtered$Cell_line




##### Step 2: Define function: Add chromosome specific data per cell to RNA Data (ex: chrm 12 GAIN or no) #####
## Step 2: Add chromosome specific data per cell to RNA Data (ex: chrm 12 gain or no)
## ex. does this cell have chrm 12 gain Y/N. (not RNA location data)
## Note: 13, 14, 15, 21, 22 chrm only one sided. 
# Define function to write change in gene expression upon chromosome gain (chrm loss is below)
# TestChrm is chromosome number you want to test. 
# nChrmArm is a number (1 or 2) corresponding to the number of arms gained.
# P-value is corrected for number of test. P-value*12k
## Note: 13, 14, 15, 21, 22 chrm only one sided. 
#dataLocation is setwd () computer location for .csv and .pdf files. 

ChrmGainEffect.RNA<-function(nChrmArm=1, 
                             dataLocation="/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/RNA_filtered/210503_plot_chrmArm_gain.loss")
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
    Chrm.gain<-merge(y= RNA.Expression.filtered_min10Cells, x= Cells.Tri.forchrm, 
                     by.y="Cell_line", by.x="DepMap_ID", 
                     sort = TRUE)# 368 cells, 12757 proteins 
    Chrm.nogain<-merge(y= RNA.Expression.filtered_min10Cells, x= Cells.Di.forchrm, 
                       by.y="Cell_line", by.x="DepMap_ID", 
                       sort = TRUE)# 368 cells, 12757 proteins 
    
  #####
  ## Step 4: get p-values and difference between RNA expression in cells with trisomy7 and without
  #length(RNA.Expression.filtered_min10Cells) #12758. so [4:12758] are RNA collumns
  
  t.test.Chrm<- data.frame(RNA_ID=as.character(), #Set up data.frame with 3 collumn
                           p.value=numeric(), 
                           Difference=numeric(), 
                           stringsAsFactors = TRUE) #this is needed to prevent errors
  
  for (w in 8:9101) {
    x<-t.test(Chrm.gain[,w], 
              Chrm.nogain[,w], 
              mu = 0, 
              alt = "two.sided",
              conf.level = 0.99) #get p-value
    Diff<-mean(Chrm.gain[,w], na.rm=TRUE) - mean(Chrm.nogain[,w], na.rm=TRUE) #get difference
    #note difference is FC. as Log2(A)-Log2(B)=Log2(A/B)
    #Difference: aneuploid - euploid
    t.test.Chrm<- rbind(t.test.Chrm, 
                        data.frame(RNA_ID=colnames(Chrm.gain[w]), 
                                   p.value=x$p.value, 
                                   Difference=Diff))
  }#Note: rbind works when you bind 2 data.frames. so need to structure data as dataframe to work
  
  
  ###
  ## Step 5: add chromosome location to each gene
  #add more RNA info per RNA. Add Uniprot_Acc and Gene_Symbol
    
    t.test.Chrm$RNA_Name<-sub("[..].*", "", as.character(t.test.Chrm$RNA_ID))#9094 genes
    t.test.Chrm$RNA_num<-sub(".*[.]{2}", "", as.character(t.test.Chrm$RNA_ID))
    t.test.Chrm$RNA_num<-sub("[.]", "", as.character(t.test.Chrm$RNA_num))
    
    
    ## Substep 3: add chromosome location to each gene
    #9094 genes to start
    #add more RNA info per RNA. Add Uniprot_Acc and Gene_Symbol
    Diff.PerRNA2<-merge(x= t.test.Chrm, y= RNA_Info3, 
                        by.x="RNA_Name", by.y="Protein_Name", 
                        sort = TRUE)# 5 collumns,  genes
    #Of genes without matching gene symbol...
    No_Gene_Symbol <- anti_join(t.test.Chrm, RNA_Info3, #finding genes_Symbol without match
                                by = c("RNA_Name" = "Protein_Name")) #4 genes
    #...find the genes with matching uniprot_ID
    Diff.PerRNA4<-merge(x= No_Gene_Symbol, y= RNA_Info3, 
                        by.x="RNA_num", by.y="Entrez Gene ID", 
                        sort = TRUE)# 10 collumns,4 genes
    # Merge genes with gene symbol and genes with only uniprot ID. 
    Diff.PerRNA5<-merge(x= Diff.PerRNA2, y= Diff.PerRNA4, 
                        all=TRUE)# 11 collumns, 9094 genes
    
    
    # make Chrm number & arm categories, and group by chrm num/arm

    Diff.PerRNA5$Chrm.Num.Arm<- paste0(Diff.PerRNA5$Chromosome, Diff.PerRNA5$arm)
    
    # Remove bad chromosome locations. ex "mitochondria m"--> not for this analysis
    Diff.PerRNA6<-subset(Diff.PerRNA5, Chrm.Num.Arm != "mitochondriaa")
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "NANA")
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "reservedd")
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "2cen-q")
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "6s") #no idea what "s" arm is...
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "7s")
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "NA")
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "22p") #only 1 gene
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "14p") #only 1 gene
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "22s") #?
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "21p") #only 1 gene
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "Yp") #?
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "Yq") #only 1 gene
    Diff.PerRNA6$Chrm.Num.Arm<-as.factor(Diff.PerRNA6$Chrm.Num.Arm) #as factor. 
  
  ##Remove values of negative infinity. this messes up values downstream
    Diff.PerRNA6 <- Diff.PerRNA6[!is.infinite(Diff.PerRNA6$Difference),]
  #12745 genes left (negative 1 infinity)
  
  ###
  ## Step 5.2: order by chromosome location
  #Now seperate RNAs on chrm of choice vs not on chrm. 
  TestChrmArm=paste0(TestChrm, TestArm)
  t.tst.Chrm.gain<- filter(Diff.PerRNA6, Chrm.Num.Arm==TestChrmArm)#RNAs on chrm 7
  
  t.tst.Chrm.nogain<- filter(Diff.PerRNA6, Chrm.Num.Arm!=TestChrmArm)#RNAs not on Chrm 7 
  
  #See if difference between RNA expression on gained chrm vs not gained chrm is significant
  #mean(t.tst.Chrm.gain$Difference)
  #mean(t.tst.Chrm.nogain$Difference)
  
  Differencetest<-t.test(t.tst.Chrm.gain$Difference, t.tst.Chrm.nogain$Difference, mu = 0, 
                         alt = "two.sided",
                         conf.level = 0.99) 
  #See if difference between RNA expression on gained chrm vs not gained chrm is significant
  
  
  ###
  ## Step 6: Seperate into quartiles. 
  ## Save all gene deregulated by chromsome
  ## save test chrm genes deregulated. with quartile number
  ## also plot quartiles if you want. 
  
  ### find top deregulated genes
  t.test.Chrm.loci_ordered<- Diff.PerRNA6[order(Diff.PerRNA6$Difference),]
  
  t.test.Chrm.loci_ordered$RNA.Gain.3cat<- cut(t.test.Chrm.loci_ordered$Difference,
                                                   breaks=c(-Inf,-0.1,0.25,Inf),
                                                   include.lowest=TRUE,
                                                   labels=c("Anti-Scaling", "Buffering","Scaling"))
  
  setwd(dataLocation)
  write.csv(t.test.Chrm.loci_ordered, 
            file =paste0("ChrmArmGain.",TestChrm, TestArm,".RNA.Diff.all.filtered_min10cells.csv"), 
            row.names = TRUE)
  
  
  ##Now plot quartiles of test chrm RNAs, and log2-fold change
  testRNAChange<-Diff.PerRNA6[which(Diff.PerRNA6$Chrm.Num.Arm==TestChrmArm),]
  
  testRNAChange<-testRNAChange[order(testRNAChange$Difference),]
  
  write.csv2(testRNAChange, 
             file =paste0("ChmArmGain.",TestChrm, TestArm,".RNAOnChrm.3cat.Filtered_min10cells.csv"), 
             row.names = TRUE)
  
  print(paste0("Finished calculating RNA changes for gain of ", TestChrm, TestArm))
  }
}

##### Step 2.5: Define function: Add chromosome specific data per cell to RNA Data (ex: chrm 12 LOSS or no) #####
#Chromosome Loss effects on RNA expression
#similar to ChrmGainEffect, but change text to say "Loss" 
#nChrmArms should be "-1" or "-2". 
ChrmLossEffect.RNA<-function(nChrmArm=-1, 
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
    Chrm.gain<-merge(y= RNA.Expression.filtered_min10Cells, x= Cells.Tri.forchrm, 
                     by.y="Cell_line", by.x="DepMap_ID", 
                     sort = TRUE)# 368 cells, 12757 proteins 
    Chrm.nogain<-merge(y= RNA.Expression.filtered_min10Cells, x= Cells.Di.forchrm, 
                       by.y="Cell_line", by.x="DepMap_ID", 
                       sort = TRUE)# 368 cells, 12757 proteins 
    
    #####
    ## Step 4: get p-values and difference between RNA expression in cells with trisomy7 and without
    #length(RNA.Expression.filtered_min10Cells) #12758. so [4:12758] are RNA collumns
    
    t.test.Chrm<- data.frame(RNA_ID=as.character(), #Set up data.frame with 3 collumn
                             p.value=numeric(), 
                             Difference=numeric(), 
                             stringsAsFactors = TRUE) #this is needed to prevent errors
    
    for (w in 8:9101) {
      x<-t.test(Chrm.gain[,w], 
                Chrm.nogain[,w], 
                mu = 0, 
                alt = "two.sided",
                conf.level = 0.99) #get p-value
      Diff<-mean(Chrm.gain[,w], na.rm=TRUE) - mean(Chrm.nogain[,w], na.rm=TRUE) #get difference
      #Difference: aneuploid - euploid
      log2Fold<- log2(mean(Chrm.gain[,w], na.rm=TRUE) / mean(Chrm.nogain[,w], na.rm=TRUE))
      t.test.Chrm<- rbind(t.test.Chrm, 
                          data.frame(RNA_ID=colnames(Chrm.gain[w]), 
                                     p.value=x$p.value, 
                                     Difference=Diff))
    }#Note: rbind works when you bind 2 data.frames. so need to structure data as dataframe to work
    
    
    ###
    ## Step 5: add chromosome location to each gene
    #add more RNA info per RNA. Add Uniprot_Acc and Gene_Symbol
    
    t.test.Chrm$RNA_Name<-sub("[..].*", "", as.character(t.test.Chrm$RNA_ID))#12755 genes
    t.test.Chrm$RNA_num<-sub(".*[.]{2}", "", as.character(t.test.Chrm$RNA_ID))
    t.test.Chrm$RNA_num<-sub("[.]", "", as.character(t.test.Chrm$RNA_num))
    
    
    ## Substep 3: add chromosome location to each gene
    # 9094 genes to start with
    # add more RNA info per RNA. Add Uniprot_Acc and Gene_Symbol
    Diff.PerRNA2<-merge(x= t.test.Chrm, y= RNA_Info3, 
                        by.x="RNA_Name", by.y="Protein_Name",
                        sort = TRUE)# 5 collumns, 9090 genes
    # Of genes without matching gene symbol... 
    No_Gene_Symbol <- anti_join(x=t.test.Chrm, y=RNA_Info3, #finding genes_Symbol without match
                                by = c("RNA_Name" = "Protein_Name")) #4 genes
    #...find the genes with matching uniprot_ID
    Diff.PerRNA4<-merge(x= No_Gene_Symbol, y= RNA_Info3, 
                        by.x="RNA_num", by.y="Entrez Gene ID", 
                        sort = TRUE)# 10 collumns, 4 genes
    # Merge genes with gene symbol and genes with only uniprot ID. 
    Diff.PerRNA5<-merge(x= Diff.PerRNA2, y= Diff.PerRNA4, 
                        all=TRUE)# 11 collumns, 9094 genes
    
    
    # make Chrm number & arm categories, and group by chrm num/arm
    
    Diff.PerRNA5$Chrm.Num.Arm<- paste0(Diff.PerRNA5$Chromosome, Diff.PerRNA5$arm)
    
    # Remove bad chromosome locations. ex "mitochondria m"--> not for this analysis
    Diff.PerRNA6<-subset(Diff.PerRNA5, Chrm.Num.Arm != "mitochondriaa")
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "NANA")
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "reservedd")
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "2cen-q")
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "6s") #no idea what "s" arm is...
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "7s")
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "NA")
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "22p") #only 1 gene
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "14p") #only 1 gene
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "22s") #?
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "21p") #only 1 gene
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "Yp") #?
    Diff.PerRNA6<-subset(Diff.PerRNA6, Chrm.Num.Arm != "Yq") #only 1 gene
    Diff.PerRNA6$Chrm.Num.Arm<-as.factor(Diff.PerRNA6$Chrm.Num.Arm) #as factor. 
    
    
    ##Remove values of negative infinity. this messes up values downstream
    Diff.PerRNA6 <- Diff.PerRNA6[!is.infinite(Diff.PerRNA6$Difference),]
    
    ###
    ## Step 5.2: order by chromosome location
    #Now seperate RNAs on chrm of choice vs not on chrm. 
    TestChrmArm=paste0(TestChrm, TestArm)
    t.tst.Chrm.gain<- filter(Diff.PerRNA6, Chrm.Num.Arm==TestChrmArm)#RNAs on chrm 7
    
    t.tst.Chrm.nogain<- filter(Diff.PerRNA6, Chrm.Num.Arm!=TestChrmArm)#RNAs not on Chrm 7 
    
    #See if difference between RNA expression on gained chrm vs not gained chrm is significant
    #mean(t.tst.Chrm.gain$Difference)
    #mean(t.tst.Chrm.nogain$Difference)
    
    Differencetest<-t.test(t.tst.Chrm.gain$Difference, t.tst.Chrm.nogain$Difference, mu = 0, 
                           alt = "two.sided",
                           conf.level = 0.99) 
    #See if difference between RNA expression on gained chrm vs not gained chrm is significant
    
    
    ###
    ## Step 6: Seperate into categories based on expression 
    ## Save all gene deregulated by chromsome
    ## save test chrm genes deregulated. with quartile number
    ## also plot quartiles if you want. 
    
    ### find top deregulated genes
    t.test.Chrm.loci_ordered<- Diff.PerRNA6[order(Diff.PerRNA6$Difference),]
    
    t.test.Chrm.loci_ordered$RNA.Loss.3cat<- cut(t.test.Chrm.loci_ordered$Difference,
                                                 breaks=c(-Inf,-0.25,0.1,Inf),
                                                 include.lowest=TRUE,
                                                 labels=c("Scaling", "Buffering","Anti-Scaling"))
    
    setwd(dataLocation)
    write.csv(t.test.Chrm.loci_ordered, 
              file =paste0("ChrmArmLoss.",TestChrm, TestArm,".RNA.Diff.all.filtered_min10cells.csv"), 
              row.names = TRUE)
    
    
    ##Now plot quartiles of test chrm RNAs, and log2-fold change
    testRNAChange<-Diff.PerRNA6[which(Diff.PerRNA6$Chrm.Num.Arm==TestChrmArm),]
    
    testRNAChange<-testRNAChange[order(testRNAChange$Difference),]

    write.csv2(testRNAChange, 
               file =paste0("ChmArmLoss.",TestChrm, TestArm,".RNAOnChrm.3cat.Filtered_min10cells.csv"), 
               row.names = TRUE)
    
    print(paste0("Finished calculating RNA changes for loss of ", TestChrm, TestArm))
  }
}
##### Step 3: Run functions #####

AllChrmArms<- data.frame(Chrm=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,14,15,16,16,17,17,18,18,19,19,20,20,21,22), 
                         Arm=c("p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q","p","q",
                               "q","q","q","p","q","p","q","p","q","p","q","p","q","q","q"))

## run functions for RNA chrm gain effect analysis
## Note: 13, 14, 15, 21, 22 chrm only one sided. 
##Chromosome gain effect on RNA: 
ChrmGainEffect.RNA(nChrmArm = 1)


##Chromosome loss
## Note: 13, 14, 15, 21, 22 chrm only one sided. 
ChrmLossEffect.RNA(nChrmArm = -1)
