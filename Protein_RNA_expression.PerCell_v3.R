###### Difference upon gain/loss dataset(s), and RNA and Protein expression vs DNA copy number boxplots ####
### 210127
### by Klaske Schukken

### Generate difference in expression upon chromosome gain and loss for CCLE dataset
###   also for diploid-only/ triploid-only cells
###   also for high-aneuploidy and low-aneuploidy cells only
###   also for difference control datasets: no mutations, no low reproducibility, etc. 

### also: plot protein/RNA expression per chromosome category (gain neutral loss)
### Graphs for Protein & RNA data: 
### Graph: protein/RNA expression per cell in gain/neutral/loss categories


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
## Get raw Protein and RNA Data. 

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

####
# Import RNA expression data (no aneuploidy data). CCLE depmap data.

## RNA data. Filtered for only cells with both RNA and Protein data: 
#  RNA data for only cells that have both RNA and Protein Data available. (371 cell lines)
## a few extra cells are filtered out when we add aneuploidy data--> 367 cells total
## RNA expression differences generated in Protein_RNA_filtered_CellLine.R
setwd(DataFileLocation)
RNA.Expression.filtered<-read.delim2("RNA_Expression_filtered.csv", 
                                     dec=",", header = TRUE, sep=";")

#RNA Info. names, chrm location, etc. 
#setwd()
RNA_Info<-Protein_Info
RNA_Info2<-RNA_Info[,c(2,3,11,12,22,14)]

RNA_Info3<-RNA_Info2
RNA_Info3$arm <- gsub('[0-9]+', '', RNA_Info3$'Chromosome band') #Find test RNA chromosome arm
RNA_Info3$arm <- gsub('[.]', '', RNA_Info3$arm) #also remove period, if needed
RNA_Info3$arm <- str_sub(RNA_Info3$arm, -1, -1) #start & end on last character. get only p or q


setwd(DataFileLocation) # set wd to where you want the graphs to go to.

# One of the datasets generated here is the main difference in RNA and protein expression upon chromosome gain/loss dataset
#CN.Diff.xRNA.yProt.ThreeGroups<- read.csv("RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points_3cat.csv")

## get list of cell lines used & info about cell lines
Cell_Lines_Used<- subset(Cell_Line_Info, DepMap_ID %in% Protein.Expression.filtered_min10Cells$Cell_line)
Cell_Lines_Used<- subset(Cell_Lines_Used, DepMap_ID %in% aneuploid$DepMap_ID)

write.xlsx(Cell_Lines_Used, "Cell_Line_Used_Info.xlsx")

###### Find RNA expression difference & t-test list, all ####
## Step 5: make list of all RNAs 
# MINIMUM OF 10 CELLS per cateogry (gain/neutral/loss), per gene. 
#colnames(RNA.Expression.filtered)<-sub("[..].*", "", as.character(colnames(RNA.Expression.filtered)))#12755 genes
#RNA.Expression.filtered<-RNA.Expression.filtered[,-c(1)]

t.test_RNA_category<- data.frame(RNA_ID=character(),
                                 RNA_Name=character(),
                                 Pvalue.Gain= numeric(),
                                 Diff.Gain= numeric(), 
                                 Pvalue.Loss= numeric(),
                                 Diff.Loss= numeric(),
                                 Pvalue.GainvLoss= numeric(),
                                 Diff.GainvLoss= numeric())

for (i in 2:length(RNA.Expression.filtered)){
  
  TestRNA_name=sub("[..].*", "", as.character(colnames(RNA.Expression.filtered[i])))
  TestRNA=colnames(RNA.Expression.filtered[i])
  
  testRNAChrm <- filter(RNA_Info3, RNA_Info3$"Approved Symbol"==TestRNA_name) #get test RNA data
  if (length(testRNAChrm$Chromosome)!=0){ #Make sure RNA data is in RNA_info3, else it crashes
    TestArm <- testRNAChrm$arm #Find test RNA chromosome arm
    TestChrm<- as.numeric(testRNAChrm$Chromosome) #Find test RNA chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test RNA location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    RNA_Expression_TestRNA <- RNA.Expression.filtered %>% select(Cell_line, all_of(TestRNA))
    
    #Combine RNA data with aneuploidy data for test RNA location
    RNA.effect.Chrm<-merge(y= RNA_Expression_TestRNA, x= testchrm.percell, 
                           by.y="Cell_line", by.x="DepMap_ID", 
                           sort = TRUE)# 368 cells, 12757 RNAs 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(RNA.effect.Chrm, arm_call==1) #all cells with chrm gain for  testRNA chrm
    Chrm.di<- filter(RNA.effect.Chrm, arm_call==0) #all cells with chrm neutral ploidy for  testRNA chrm
    Chrm.mono<- filter(RNA.effect.Chrm, arm_call==-1) #all cells with chrm loss for  testRNA chrm
    
    ##Make data frame with t-test info about RNA expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 &  #minimum 10 cells in "tri" category
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 &
        colSums(!is.na(Chrm.mono[8]))>=10) {
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with RNA data
                     Chrm.di[,8], # [,8] because that is column with RNA data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Mono.Di<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                      Chrm.di[,8], # [,8] because that is column with RNA data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Mono.Di<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                       Chrm.tri[,8], # [,8] because that is column with RNA data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_RNA_category<- rbind(t.test_RNA_category, data.frame(
        RNA_ID=colnames(RNA.effect.Chrm[8]),
        RNA_Name=TestRNA_name,
        Pvalue.Gain= di.tri$p.value,
        Diff.Gain= Diff.Tri.Di, 
        Pvalue.Loss= Mono.Di$p.value,
        Diff.Loss= Diff.Mono.Di, 
        Pvalue.GainvLoss= Tri.Mono$p.value,
        Diff.GainvLoss= Diff.Tri.Mono))
    }
  }
}
t.test_RNA_category<-distinct(t.test_RNA_category)#17774 genes
t.test_RNA_category<-t.test_RNA_category[order(t.test_RNA_category$Pvlaue.Tri.Mono),]

#17621 RNAs measured. 
setwd(DataFileLocation)
write.csv(t.test_RNA_category, 
          file =paste("RNA_Loss.Neutral.Gain_Difference_Pvalue_min10points.csv", sep=','), 
          row.names = TRUE)

t.test_RNA_category<-read.delim2("RNA_Loss.Neutral.Gain_Difference_Pvalue_min10points.csv", 
                                 dec=".", header = TRUE, sep=",")

###### Find Protein expression difference & t-test list, all ####
## make list of pvalue and difference between mono- di-tri cells per gene 
## for each gene in Protein data. using only cells that have both RNA and Protein data

##First step is to combine Protein Uniprot ID info with Depmap Protein_ID
## Using information given by initial mmcr protein database. see above. 
## Get Protein_ID as dataframe, combine with uniprot id and Protein symbol data. 
Protein.ID<-data.frame(Protein_ID=
                         colnames(Protein.Expression.filtered[3:length(Protein.Expression.filtered)]))

##Now merge depmap column name info with Protein info. 
## First add Uniprot ID and gene symbol info to depmap data. 
Depmap.Protein.info2<-merge(x= Protein.ID, y= Protein_ProID, 
                            by.x="Protein_ID", by.y="Protein_Id", 
                            sort = TRUE)# 3 collumns, 12755. all Proteins given Uniprot IDs.
# Now combine all genes with same uniprot ID: 
Depmap.Protein.info3<-merge(x= Depmap.Protein.info2, y= Protein_Info4, 
                            by.x="Uniprot_Acc", by.y="UniProt ID(supplied by UniProt)", 
                            sort = TRUE)# 10594 Proteins
# Find those genes without uniprot id: 
No_UniprotID <- anti_join(Depmap.Protein.info2, Protein_Info4, #finding genes_Symbol without match
                          by = c("Uniprot_Acc" = "UniProt ID(supplied by UniProt)"))
#...find the genes with matching gene symbols
Depmap.Protein.info4<-merge(x= No_UniprotID, y= Protein_Info4, 
                            by.x="Gene_Symbol", by.y="Approved Symbol", 
                            sort = TRUE)# 10 collumns, 253 genes, only no gene-symbol genes with uniprot_ID

# Merge genes with gene symbol and genes with only uniprot ID. 
Depmap.Protein.info5<-merge(x= Depmap.Protein.info3, y= Depmap.Protein.info4, 
                            all=TRUE)# 11 collumns, 12100 genes


t.test_prot_category<- data.frame(Protein_ID=character(),
                                  Protein_Name=character(),
                                  Pvlaue.Tri.Di= numeric(),
                                  Diff.Gain= numeric(), 
                                  Pvlaue.Di.Mono= numeric(),
                                  Diff.Di_Mono= numeric(),
                                  Pvlaue.Tri.Mono= numeric(),
                                  Diff.GainvLoss= numeric())



for (i in 3:length(Protein.Expression.filtered)){
  TestProt = colnames(Protein.Expression.filtered[i])
  
  testProtChrm <- filter(Depmap.Protein.info5, Protein_ID==TestProt) #get test protein data
  
  
  #use If statement to check that Protein Info has Protein of interest. 
  if (length(testProtChrm[,1])!=0){
    TestArm <- testProtChrm$arm #Find test protein chromosome arm
    TestChrm<- as.numeric(testProtChrm$Chromosome) #Find test protein chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test protein location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    Protein_Expression_TestProt <- Protein.Expression.filtered %>% select(Cell_line, all_of(TestProt))
    
    #Combine protein data with aneuploidy data for test protein location
    Protein.effect.Chrm<-merge(y= Protein_Expression_TestProt, x= testchrm.percell, 
                               by.y="Cell_line", by.x="DepMap_ID", 
                               sort = TRUE)# 368 cells, 12757 proteins 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(Protein.effect.Chrm, arm_call==1) #all cells triploid for testProt chrm
    Chrm.di<- filter(Protein.effect.Chrm, arm_call==0) #all cells diploid for testProt chrm
    Chrm.mono<- filter(Protein.effect.Chrm, arm_call==-1) #all cells mono for testProt chrm
    
    ##Make data frame with t-test info about protein expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 & #[8] is column with protein expression data. check it has values. 
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 &
        colSums(!is.na(Chrm.mono[8]))>=10) {
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with Protein data
                     Chrm.di[,8], # [,8] because that is column with Protein data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Di.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                      Chrm.di[,8], # [,8] because that is column with Protein data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Di.Mono<-mean(Chrm.di[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      Di.Mono.ploidy<- mean(Chrm.di[,6], na.rm=TRUE)
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                       Chrm.tri[,8], # [,8] because that is column with Protein data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      Tri.Mono.ploidy<- mean(Chrm.tri[,6], na.rm=TRUE)
      
      t.test_prot_category<- rbind(t.test_prot_category, 
                                   data.frame(Protein_ID=TestProt,
                                              Protein_Name=testProtChrm$Gene_Symbol,
                                              Pvlaue.Tri.Di= di.tri$p.value,
                                              Diff.Gain= Diff.Tri.Di, 
                                              Pvlaue.Di.Mono= Di.Mono$p.value,
                                              Diff.Di_Mono= Diff.Di.Mono, 
                                              Pvlaue.Tri.Mono= Tri.Mono$p.value,
                                              Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_prot_category<-t.test_prot_category
    }
  }
}


t.test_prot_category<-t.test_prot_category[order(t.test_prot_category$Pvlaue.Tri.Mono),]#11458 genes. only the proteins with 3+ cells per mono/di.tri category. lost ~1k proteins. 
t.test_prot_category<-distinct(t.test_prot_category)
t.test_prot_category$Diff.Di_Mono<-t.test_prot_category$Diff.Di_Mono*-1 #make this into Aneu-Diploid.

setwd(DataFileLocation)
write.csv(t.test_prot_category, 
          file =paste("Protein_Loss.Neutral.Gain_Difference_Pvalue_min10points.csv", sep=','), 
          row.names = TRUE)

#t.test_prot_category<-read.delim2("Protein_Loss.Neutral.Gain_Difference_Pvalue_min10points.csv", 
#                                  dec=".", header = TRUE, sep=",")


### Combine RNA and Protein Data: 
## This is with filtered cell lines (only cells with RNA and Protein and aneuploid data:)
## This is with filtered genes: only genes with 10+ cells in all categories: RNA & protein gain, neutral and loss
CN.Diff.xRNA.yProt<-merge(x=t.test_RNA_category, y=t.test_prot_category, by.x="RNA_Name", by.y="Protein_Name")
#CN.Diff.xRNA.yProt<- subset(CN.Diff.xRNA.yProt, select=-c(X,X.x, X.y, X.1,X.2))

colnames(CN.Diff.xRNA.yProt)<-c("RNA_Name", "RNA_ID", 
                                "RNA.Pvalue.Gain", "RNA.Diff.Gain", 
                                "RNA.Pvalue.Loss", "RNA.Diff.Loss", 
                                "RNA.Pvalue.GainvLoss", "RNA.Diff.GainvLoss", 
                                "Protein_ID", "Protein.Pvalue.Gain", 
                                "Protein.Diff.Gain", "Protein.Pvalue.Loss", 
                                "Protein.Diff.Loss", "Protein.Pvalue.GainvLoss", 
                                "Protein.Diff.GainvLoss")
#CN.Diff.xRNA.yProt<-CN.Diff.xRNA.yProt[,-c(2,10:11)]

write.csv(CN.Diff.xRNA.yProt, 
          file =paste("RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points.csv", sep=','), 
          row.names = TRUE)

#read.csv(CN.Diff.xRNA.yProt, 
#          file =paste("RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points.csv", sep=','), 
#          row.names = TRUE)

write.csv(CN.Diff.xRNA.yProt.ThreeGroups, 
          file =paste("RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points_3cat.csv", sep=','), 
          row.names = TRUE)

CN.Diff.xRNA.yProt.ThreeGroups<- read.csv("RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points_3cat.csv") 


###### Categorize genes by significance and difference #####
## Categorize genes by if they scale or buffer upon chrm gain/loss. plot. 
## Did not end up using these graphs in the paper. but can use this to find highly 
##    significant scaling/buffering proteins/RNAs. 

## Find genes per category: Per chrm gain: 
## Categories based on Significance and +/- 0 Difference,
## can be used to find genes that are significantly Scaling/Buffering (or NS) for gain & loss

CN.Diff.xRNA.yProt<-merge(x=t.test_RNA_category, y=t.test_prot_category, by.x="RNA_ID", by.y="Protein_Name")

CN.gain.RNAbuff.Protbuff<- subset(CN.Diff.xRNA.yProt, 
                                  (RNA.Pvalue.Gain>0.05 | RNA.Diff.Gain<0) & 
                                    (Protein.Pvalue.Gain>0.05 | Protein.Diff.Gain<0) )
CN.gain.RNAbuff.Protbuff<- CN.gain.RNAbuff.Protbuff[order(CN.gain.RNAbuff.Protbuff$Protein.Pvalue.Gain),]
CN.gain.RNAscale.Protscale<- subset(CN.Diff.xRNA.yProt, 
                                    (RNA.Pvalue.Gain<0.05 & RNA.Diff.Gain>0) & 
                                      (Protein.Pvalue.Gain<0.05 & Protein.Diff.Gain>0) )
CN.gain.RNAscale.Protscale<- CN.gain.RNAscale.Protscale[order(CN.gain.RNAscale.Protscale$Protein.Pvalue.Gain),]
CN.gain.RNAscale.Protbuff<- subset(CN.Diff.xRNA.yProt, 
                                   (RNA.Pvalue.Gain<0.05 & RNA.Diff.Gain>0) & 
                                     (Protein.Pvalue.Gain>0.05 | Protein.Diff.Gain<0) )
CN.gain.RNAscale.Protbuff<- CN.gain.RNAscale.Protbuff[order(CN.gain.RNAscale.Protbuff$Protein.Pvalue.Gain),]
CN.gain.RNAbuff.Protscale<- subset(CN.Diff.xRNA.yProt, 
                                   (RNA.Pvalue.Gain>0.05 | RNA.Diff.Gain<0) &
                                     (Protein.Pvalue.Gain<0.05 & Protein.Diff.Gain>0) )
CN.gain.RNAbuff.Protscale<- CN.gain.RNAbuff.Protscale[order(CN.gain.RNAbuff.Protscale$Protein.Pvalue.Gain),]

length(CN.gain.RNAbuff.Protbuff$RNA_ID) 
CN.gain.RNAbuff.Protbuff$RNA_Name[1:5]
length(CN.gain.RNAscale.Protscale$RNA_ID) 
CN.gain.RNAscale.Protscale$RNA_Name[1:5]
length(CN.gain.RNAscale.Protbuff$RNA_ID) 
CN.gain.RNAscale.Protbuff$RNA_Name[1:5]
length(CN.gain.RNAbuff.Protscale$RNA_ID) 
CN.gain.RNAbuff.Protscale$RNA_Name[1:5]

## Upon chromosome gain: 
##Minimum of 3 data points per condition (-1, 0, +1):
#RNA buff, Protein Buffered: 4737 genes 
#     top 5: IRF1    DTNA    CYBRD1  ZMYND19 FCGRT  
#RNA scale, Protein scales: 2309 genes 
#     top 5 most significant: PURB   NUDCD3 STAU1  RPRD1B CDK5
#RNA scale, Protein buff: 3450 genes 
#     top 5 most significant: ACTR5   MLKL    C9orf85 TMEM14C MTCH1  
#RNA buff, Protein scales: 445 genes 
#     top 5 most significant: IP6K2  MUM1   ABHD10 DOK3   RNF223

## Upon chromosome gain: 
##Minimum of 10 data points per condition (-1, 0, +1):
#RNA buff, Protein Buffered: 3849 genes 
#     top 5: LAMC2 AP1M2 DOK3  PCDH9 LAMB3  
#RNA scale, Protein scales: 2205 genes 
#     top 5 most significant: PURB   NUDCD3 STAU1  RPRD1B CDK5
#RNA scale, Protein buff: 3013 genes 
#     top 5 most significant: TNNC2    TRAPPC12 RFC2     HAUS4    PSMB2     
#RNA buff, Protein scales: 346 genes 
#     top 5 most significant: ABHD10   TMEM126A RABEP2   ALDOA    PGAM5


####
##Find categories for chromosome loss: 
CN.loss.RNAbuff.Protbuff<- subset(CN.Diff.xRNA.yProt, 
                                  (RNA.Pvalue.Loss>0.05 | RNA.Diff.Loss>0) & 
                                    (Protein.Pvalue.Loss>0.05 | Protein.Diff.Loss<0) )
CN.loss.RNAbuff.Protbuff<- CN.loss.RNAbuff.Protbuff[order(CN.loss.RNAbuff.Protbuff$Protein.Pvalue.Loss),]

CN.loss.RNAscale.Protscale<- subset(CN.Diff.xRNA.yProt, 
                                    (RNA.Pvalue.Loss<0.05 & RNA.Diff.Loss<0) & 
                                      (Protein.Pvalue.Loss<0.05 & Protein.Diff.Loss<0) )
CN.loss.RNAscale.Protscale<- CN.loss.RNAscale.Protscale[order(CN.loss.RNAscale.Protscale$Protein.Pvalue.Loss),]

CN.loss.RNAscale.Protbuff<- subset(CN.Diff.xRNA.yProt, 
                                   (RNA.Pvalue.Loss<0.05 & RNA.Diff.Loss<0) & 
                                     (Protein.Pvalue.Loss>0.05 | Protein.Diff.Loss>0) )
CN.loss.RNAscale.Protbuff<- CN.loss.RNAscale.Protbuff[order(CN.loss.RNAscale.Protbuff$Protein.Pvalue.Loss),]

CN.loss.RNAbuff.Protscale<- subset(CN.Diff.xRNA.yProt, 
                                   (RNA.Pvalue.Loss>0.05 | RNA.Diff.Loss>0) & 
                                     (Protein.Pvalue.Loss<0.05 & Protein.Diff.Loss<0) )
CN.loss.RNAbuff.Protscale<- CN.loss.RNAbuff.Protscale[order(CN.loss.RNAbuff.Protscale$Protein.Pvalue.Loss),]

length(CN.loss.RNAbuff.Protbuff$RNA_ID) 
CN.loss.RNAbuff.Protbuff$RNA_Name[1:5]
length(CN.loss.RNAscale.Protscale$RNA_ID) 
CN.loss.RNAscale.Protscale$RNA_Name[1:5]
length(CN.loss.RNAscale.Protbuff$RNA_ID) 
CN.loss.RNAscale.Protbuff$RNA_Name[1:5]
length(CN.loss.RNAbuff.Protscale$RNA_ID) 
CN.loss.RNAbuff.Protscale$RNA_Name[1:5]

## Upon chromosome loss: 
##Minimum of 3 data points per condition (-1, 0, +1):
#RNA buff, Protein Buffered: 4059 genes 
#     top 5: GLI2   EDIL3  GMPS   TXNRD3 CEBPB
#RNA scale, Protein scales: 3021 genes 
#     top 5 most significant: PDE12  TXNL1  HDHD2  NARS   PPP2CB
#RNA scale, Protein buff: 3503 genes 
#     top 5 most significant: RASSF2 ZNF622 ADD3   ATMIN  NHSL1 
#RNA buff, Protein scales: 358 genes 
#     top 5 most significant: DNAJB1   MEX3D    UBE2E2   CHMP2B   ARHGEF18

## Upon chromosome loss: 
##Minimum of 3 data points per condition (-1, 0, +1):
#RNA buff, Protein Buffered: 3386 genes 
#     top 5: DSG2     DSG3     CYBRD1   SERPINB5 S100A10 
#RNA scale, Protein scales: 2855 genes 
#     top 5 most significant: PDE12  TXNL1  HDHD2  NARS   PPP2CB
#RNA scale, Protein buff: 2911 genes 
#     top 5 most significant: CTNNB1 TAF13  TRIM11 STARD9 SPCS1  
#RNA buff, Protein scales: 261 genes 
#     top 5 most significant: DNAJB1   MEX3D    UBE2E2   CHMP2B   ARHGEF18

###Now to look at genes that are buffeted/scaling in BOTH gain and loss 
CN.aneu.RNAbuff.Protbuff<- subset(CN.Diff.xRNA.yProt, 
                                  ((RNA.Pvalue.Loss>0.05 | RNA.Diff.Loss>0) & (RNA.Pvalue.Gain>0.05 |RNA.Diff.Gain<0) & (RNA.Pvalue.GainvLoss>0.05| RNA.Diff.GainvLoss<0)) & # RNA not sig OR wrong direction
                                    ((Protein.Pvalue.Loss>0.05 | Protein.Diff.Loss>0) & (Protein.Pvalue.Gain>0.05 |Protein.Diff.Gain<0) & (Protein.Pvalue.GainvLoss>0.05| Protein.Diff.GainvLoss<0)))  # Protein not sig OR wrong direction
CN.aneu.RNAbuff.Protbuff<- CN.aneu.RNAbuff.Protbuff[order(CN.aneu.RNAbuff.Protbuff$Protein.Pvalue.GainvLoss),]

CN.aneu.RNAscale.Protscale<- subset(CN.Diff.xRNA.yProt, 
                                    ((RNA.Pvalue.Loss<0.05 & RNA.Diff.Loss<0) & (RNA.Pvalue.Gain<0.05 & RNA.Diff.Gain>0) & (RNA.Pvalue.GainvLoss<0.05 & RNA.Diff.GainvLoss>0)) & # RNA sig AND scales directionally
                                      ((Protein.Pvalue.Loss<0.05 & Protein.Diff.Loss<0) & (Protein.Pvalue.Gain<0.05 & Protein.Diff.Gain>0) & (Protein.Pvalue.GainvLoss<0.05 & Protein.Diff.GainvLoss>0)))# Protein sig AND scales directionally
CN.aneu.RNAscale.Protscale<- CN.aneu.RNAscale.Protscale[order(CN.aneu.RNAscale.Protscale$Protein.Pvalue.GainvLoss),]

CN.aneu.RNAscale.Protbuff<- subset(CN.Diff.xRNA.yProt, 
                                   ((RNA.Pvalue.Loss<0.05 & RNA.Diff.Loss<0) & (RNA.Pvalue.Gain<0.05 & RNA.Diff.Gain>0) & (RNA.Pvalue.GainvLoss<0.05 & RNA.Diff.GainvLoss>0)) & # RNA sig AND scales directionally
                                     ((Protein.Pvalue.Loss>0.05 | Protein.Diff.Loss>0) & (Protein.Pvalue.Gain>0.05 |Protein.Diff.Gain<0) & (Protein.Pvalue.GainvLoss>0.05| Protein.Diff.GainvLoss<0)))  # Protein not sig OR wrong direction
CN.aneu.RNAscale.Protbuff<- CN.aneu.RNAscale.Protbuff[order(CN.aneu.RNAscale.Protbuff$Protein.Pvalue.GainvLoss),]

CN.aneu.RNAbuff.Protscale<- subset(CN.Diff.xRNA.yProt, 
                                   ((RNA.Pvalue.Loss>0.05 | RNA.Diff.Loss>0) & (RNA.Pvalue.Gain>0.05 |RNA.Diff.Gain<0) & (RNA.Pvalue.GainvLoss>0.05| RNA.Diff.GainvLoss<0)) & # RNA not sig OR wrong direction
                                     ((Protein.Pvalue.Loss<0.05 & Protein.Diff.Loss<0) & (Protein.Pvalue.Gain<0.05 & Protein.Diff.Gain>0) & (Protein.Pvalue.GainvLoss<0.05 & Protein.Diff.GainvLoss>0)))# Protein sig AND scales directionally
CN.aneu.RNAbuff.Protscale<- CN.aneu.RNAbuff.Protscale[order(CN.aneu.RNAbuff.Protscale$Protein.Pvalue.GainvLoss),]

length(CN.aneu.RNAbuff.Protbuff$RNA_ID) 
CN.aneu.RNAbuff.Protbuff$RNA_Name[1:5]
length(CN.aneu.RNAscale.Protscale$RNA_ID) 
CN.aneu.RNAscale.Protscale$RNA_Name[1:5]
length(CN.aneu.RNAscale.Protbuff$RNA_ID) 
CN.aneu.RNAscale.Protbuff$RNA_Name[1:5]
length(CN.aneu.RNAbuff.Protscale$RNA_ID) 
CN.aneu.RNAbuff.Protscale$RNA_Name[1:5]

## Buffered upon BOTH gain and loss, scale with either one or other: 
##Minimum of 3 data points per condition (-1, 0, +1):
#RNA buff, Protein Buffered: 1489 genes 
#     top 5: IRF1   IAH1   NOS1AP RAB21  HPCAL4
#RNA scale, Protein scales: 884 genes 
#     top 5 most significant: PTPN2  SMCHD1 RNMT   USP14  PPP4R1
#RNA scale, Protein buff: 1425 genes 
#     top 5 most significant: DICER1  CARD8   RAPGEF1 DCLRE1C MDM4 
#RNA buff, Protein scales: 16 genes 
#     top 5 most significant: SIX5     STMN4    GOLGA3   GLB1L3   ARHGEF28

## Buffered upon BOTH gain and loss, scale with either one or other: 
##Minimum of 10 data points per condition (-1, 0, +1):
#RNA buff, Protein Buffered: 1705 genes 
#     top 5: DSG3     AP1M2    DSC2     SERPINB5 S100A10
#RNA scale, Protein scales: 860 genes 
#     top 5 most significant: PTPN2  SMCHD1 RNMT   USP14  PPP4R1
#RNA scale, Protein buff: 1014 genes 
#     top 5 most significant: RAB11B COPG1  HDAC3  TRIM11 HAUS4 
#RNA buff, Protein scales: 3 genes 
#     top 5 most significant: GLB1L3  ST3GAL2 SLC27A2 

ProtExp.ChrmCN.filtered("GLB1L3")
RNAExp.ChrmCN.filtered("GLB1L3")
ProtExp.ChrmCN.filtered("ST3GAL2")
RNAExp.ChrmCN.filtered("ST3GAL2")
ProtExp.ChrmCN.filtered("SLC27A2")
RNAExp.ChrmCN.filtered("SLC27A2")


###Some other genes that I had for some reason
ProtExp.ChrmCN.filtered("DSG3")
RNAExp.ChrmCN.filtered("DSG3")
ProtExp.ChrmCN.filtered("ABL2")
RNAExp.ChrmCN.filtered("ABL2")
ProtExp.ChrmCN.filtered("STMN4")
RNAExp.ChrmCN.filtered("STMN4")
ProtExp.ChrmCN.filtered("TROAP")
RNAExp.ChrmCN.filtered("TROAP")
ProtExp.ChrmCN.filtered("ERBB4")
RNAExp.ChrmCN.filtered("ERBB4")


#Make category vectors for buffer/scale upon gain, loss or both: 

CN.Diff.xRNA.yProt$Gain<-NA 
for (i in 1:length(CN.Diff.xRNA.yProt$RNA_ID)){
  if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.gain.RNAbuff.Protbuff$RNA_ID){
    CN.Diff.xRNA.yProt$Gain[i]<-"RNA buffered, Protein buffered"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.gain.RNAscale.Protscale$RNA_ID){
    CN.Diff.xRNA.yProt$Gain[i]<-"RNA scales, Protein scales"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.gain.RNAscale.Protbuff$RNA_ID){
    CN.Diff.xRNA.yProt$Gain[i]<-"RNA scales, Protein buffered"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.gain.RNAbuff.Protscale$RNA_ID){
    CN.Diff.xRNA.yProt$Gain[i]<-"RNA buffered, Protein scales"
  }
}
CN.Diff.xRNA.yProt$Gain<-factor(CN.Diff.xRNA.yProt$Gain, levels=c("RNA buffered, Protein buffered", "RNA scales, Protein buffered", "RNA scales, Protein scales", "RNA buffered, Protein scales"))

CN.Diff.xRNA.yProt$Loss<-NA 
for (i in 1:length(CN.Diff.xRNA.yProt$RNA_ID)){
  if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.loss.RNAbuff.Protbuff$RNA_ID){
    CN.Diff.xRNA.yProt$Loss[i]<-"RNA buffered, Protein buffered"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.loss.RNAscale.Protscale$RNA_ID){
    CN.Diff.xRNA.yProt$Loss[i]<-"RNA scales, Protein scales"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.loss.RNAscale.Protbuff$RNA_ID){
    CN.Diff.xRNA.yProt$Loss[i]<-"RNA scales, Protein buffered"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.loss.RNAbuff.Protscale$RNA_ID){
    CN.Diff.xRNA.yProt$Loss[i]<-"RNA buffered, Protein scales"
  }
}
CN.Diff.xRNA.yProt$Loss<-factor(CN.Diff.xRNA.yProt$Loss, levels=c("RNA buffered, Protein buffered", "RNA scales, Protein buffered", "RNA scales, Protein scales", "RNA buffered, Protein scales"))

CN.Diff.xRNA.yProt$Gain.Loss<-NA 
for (i in 1:length(CN.Diff.xRNA.yProt$RNA_ID)){
  if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.aneu.RNAbuff.Protbuff$RNA_ID){
    CN.Diff.xRNA.yProt$Gain.Loss[i]<-"RNA buffered, Protein buffered"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.aneu.RNAscale.Protscale$RNA_ID){
    CN.Diff.xRNA.yProt$Gain.Loss[i]<-"RNA scales, Protein scales"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.aneu.RNAscale.Protbuff$RNA_ID){
    CN.Diff.xRNA.yProt$Gain.Loss[i]<-"RNA scales, Protein buffered"
  } else if (CN.Diff.xRNA.yProt$RNA_ID[i] %in% CN.aneu.RNAbuff.Protscale$RNA_ID){
    CN.Diff.xRNA.yProt$Gain.Loss[i]<-"RNA buffered, Protein scales"
  }
}

CN.Diff.xRNA.yProt$Gain.Loss<-factor(CN.Diff.xRNA.yProt$Gain.Loss, levels=c("RNA buffered, Protein buffered", "RNA scales, Protein buffered", "RNA scales, Protein scales", "RNA buffered, Protein scales"))

table.diff.scale<-table(CN.Diff.xRNA.yProt$Gain)
table.diff.scale<-rbind(table.diff.scale,table(CN.Diff.xRNA.yProt$Loss))
table.diff.scale<-rbind(table.diff.scale,table(CN.Diff.xRNA.yProt$Gain.Loss))
rownames(table.diff.scale)<-c("Gain", "Loss", "Gain and Loss")

melt.table.diff.scale<- melt(table.diff.scale)

ggplot(data=melt.table.diff.scale, aes(x=Var1, y=value))+
  geom_bar(stat="identity", position="dodge", aes(fill=Var2))+
  xlab("Gene scaling or buffering upon chromosome aneuploidy")+
  ylab("Gene count")+ 
  scale_fill_manual(values=c("grey20", "grey40", "grey60", "grey80"))+
  labs(fill=c("Gene regulation type (significance):"))+
  theme_classic()




###### Categorize genes: Scaling, buffering. anti-scaling ####
## Make THREE Categories: scaling, anti-scale and buffering, based on difference. 
## Three groups: Scaling, buffering and anti-scaling
## Three categories, upon gain: sclaing (>0.25), buffering (-0.1 to 0.25), anti-scaling (< -0.1)
## 0.25 because thats halfway between .5 fold change, expected from 1.5-fold increase. 
## -0.1 to account for non-significant downregulation upon gain
## vice versa for chrm loss. 
CN.Diff.xRNA.yProt.ThreeGroups<- CN.Diff.xRNA.yProt
CN.Diff.xRNA.yProt.ThreeGroups$Three.RNA.Gain<- cut(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain,
                                                    breaks=c(-Inf,-0.1,0.25,Inf),
                                                    include.lowest=TRUE,
                                                    labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.xRNA.yProt.ThreeGroups$Three.RNA.Loss<- cut(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss,
                                                    breaks=c(-Inf,-0.25,0.1,Inf),
                                                    include.lowest=TRUE,
                                                    labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.xRNA.yProt.ThreeGroups$Three.Protein.Gain<- cut(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain,
                                                        breaks=c(-Inf,-0.1,0.25,Inf),
                                                        include.lowest=TRUE,
                                                        labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.xRNA.yProt.ThreeGroups$Three.Protein.Loss<- cut(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss,
                                                        breaks=c(-Inf,-0.25,0.1,Inf),
                                                        include.lowest=TRUE,
                                                        labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.xRNA.yProt.ThreeGroups$Three.Protein.Loss<-factor(CN.Diff.xRNA.yProt.ThreeGroups$Three.Protein.Loss, 
                                                          levels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.xRNA.yProt.ThreeGroups$Three.RNA.Loss<-factor(CN.Diff.xRNA.yProt.ThreeGroups$Three.RNA.Loss, 
                                                          levels=c("Anti-Scaling","Buffering","Scaling"))

### count percentages
sum(CN.Diff.xRNA.yProt.ThreeGroups$Three.Protein.Loss=="Anti-Scaling")/length(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name)
x<-subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.RNA.Gain=="Buffering")
sum(x$Three.Protein.Gain=="Buffering")/length(x$RNA_Name)
x<-subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.Protein.Loss=="Buffering")
sum(x$Three.Protein.Gain=="Buffering")/length(x$RNA_Name)


##Scatterplots with categories colored in: 
#RNA Chrm gain graphs
ggplot(CN.Diff.xRNA.yProt.ThreeGroups, aes(x=RNA.Diff.Gain, y=-log2(RNA.Pvalue.Gain), color=Three.RNA.Gain))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in RNA expression")+
  ylab("log2(p-value)")+
  labs(color = "Quartiles:")+ #legend title
  theme_classic()+
  coord_cartesian(xlim=c(-2.7, 2.7), ylim=c(0,100))+
  ggtitle("RNA Gain scatterplot quartiles")

#RNA Loss graphs
ggplot(CN.Diff.xRNA.yProt.ThreeGroups, aes(x=RNA.Diff.Loss, y=-log2(RNA.Pvalue.Loss), color=Three.RNA.Loss))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in RNA expression")+
  ylab("log2(p-value)")+
  labs(color = "Quartiles:")+ #legend title
  theme_classic()+
  coord_cartesian(xlim=c(-2.7, 2.7), ylim=c(0,100))+
  ggtitle("RNA Loss scatterplot quartiles")

#Protein Chrm gain graphs
ggplot(CN.Diff.xRNA.yProt.ThreeGroups, aes(x=Protein.Diff.Gain, y=-log2(Protein.Pvalue.Gain), color=Three.Protein.Gain))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in protein expression")+
  ylab("log2(p-value)")+
  labs(color = "Quartiles:")+ #legend title
  theme_classic()+
  coord_cartesian(xlim=c(-2.7, 2.7), ylim=c(0,100))+
  ggtitle("Protein Gain scatterplot quartiles")

#Protein Loss graphs
ggplot(CN.Diff.xRNA.yProt.ThreeGroups, aes(x=Protein.Diff.Loss, y=-log2(Protein.Pvalue.Loss), color=Three.Protein.Loss))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in protein expression")+
  ylab("log2(p-value)")+
  labs(color = "Quartiles:")+ #legend title
  theme_classic()+
  coord_cartesian(xlim=c(-2.7, 2.7), ylim=c(0,100))+
  ggtitle("protein Loss scatterplot quartiles")


###Bar graph of RNA and Protein quantiles Scaling/Buffered

#format data so I can plot all 4 groups: RNA/Prot gain/loss
dat.m <- melt(CN.Diff.xRNA.yProt.ThreeGroups,id.vars='RNA_ID', measure.vars=c('Three.RNA.Gain','Three.Protein.Gain','Three.RNA.Loss', 'Three.Protein.Loss'))

## RNA & Protein Gain
ggplot(data= dat.m, aes(x=variable, fill=value)) + 
  geom_bar(position="fill")+ #dodge (next to eachother)
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("RNA Difference upon chrm arm gain: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  #coord_cartesian(ylim=c(0,6000))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Relationship between RNA and Protein Difference\nupon chromosome arm gain: per Category")
#4x4
# plot.Protein.RNA.Bargraph.ThreeCat.Gain


####
###Bar graph of RNA and Protein quantiles Scaling/Buffered
## RNA & Protein Gain
ggplot(data= CN.Diff.xRNA.yProt.ThreeGroups, aes(x=Three.RNA.Gain, fill=Three.Protein.Gain)) + 
  geom_bar(position="fill")+
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("RNA Difference upon chrm arm gain: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Relationship between RNA and Protein Difference\nupon chromosome arm gain: per Category")
#4x4
# plot.Protein.RNA.Bargraph.ThreeCat.Gain

#RNA & Protein Loss
ggplot(data= CN.Diff.xRNA.yProt.ThreeGroups, aes(x=Three.RNA.Loss, fill=Three.Protein.Loss)) + 
  geom_bar(position="fill")+
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("RNA Difference upon chrm arm loss: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Relationship between RNA and Protein Difference\nupon chromosome arm loss: per Category")
#4x4
# plot.Protein.RNA.Bargraph.ThreeCat.Loss

## RNA Gain & Loss
ggplot(data= CN.Diff.xRNA.yProt.ThreeGroups, aes(x=Three.RNA.Gain, fill=Three.RNA.Loss)) + 
  geom_bar(position="fill")+
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("RNA Difference upon chrm arm gain: category")+
  ylab("Percent of genes")+
  labs(fill = "RNA Difference upon loss:\nCategory")+ #legend title
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Relationship between RNA genes \nupon chromosome arm gain or loss")
#4x4
# plot.Protein.RNA.Bargraph.ThreeCat.Gain

## Protein Gain & Loss
ggplot(data= CN.Diff.xRNA.yProt.ThreeGroups, aes(x=Three.Protein.Gain, fill=Three.Protein.Loss)) + 
  geom_bar(position="fill")+
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Protein difference upon chrm arm gain: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference upon loss:\nCategory")+ #legend title
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Relationship between protein genes \nupon chromosome arm gain or loss")
#4x4
# plot.Protein.RNA.Bargraph.ThreeCat.Gain


###### Protein vs RNA gain/loss correlation graphs ######
## density plot of change in RNA & Protein upon chrm gain/Loss 
## using only data with 10+ datapoints per condition. 

CN.Diff.xRNA.yProt.ThreeGroups
#9413 genes analyzed
# CHROMOSOME GAIN: correlate RNA and Protein expression difference
# Mimimum of 10 datapoints per condition (gain, no aneu, loss)
ggplot(CN.Diff.xRNA.yProt.ThreeGroups, aes(x=CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain, y=CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  xlab("Difference in RNA expression")+
  ylab("Difference in protein expression")+
  theme_classic()+
  geom_hline(yintercept=0.0)+
  geom_vline(xintercept=0.0)+
  geom_smooth(method="lm", color="Red")+
  coord_cartesian(xlim=c(-0.7, 0.7), ylim=c(-0.7,0.7))+
  ggtitle("Chromosome gain")
#4x4
# plot.Protein.RNA.expression.density.Gain
Chrm_Gain_min10<-cor.test(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain, CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain, method="pearson")
#cor=0.546
#P < 2E-16


# CHOMOSOME LOSS: correlate RNA and Protein expression difference
ggplot(CN.Diff.xRNA.yProt.ThreeGroups, aes(x=RNA.Diff.Loss, y=Protein.Diff.Loss))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  xlab("Difference in RNA expression")+
  ylab("Difference in protein expression")+
  theme_classic()+
  geom_hline(yintercept=0.0)+
  geom_vline(xintercept=0.0)+
  geom_smooth(method="lm", color="Red")+
  coord_cartesian(xlim=c(-0.7, 0.7), ylim=c(-0.7,0.7))+
  ggtitle("Chromosome loss")
#4x4
#plot.Protein.RNA.expression.density.Loss
Chrm_Loss_min10<-cor.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss, CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss, method="pearson")
Chrm_Loss_min10
#cor=0.554
#P<2E-16




###### Define function: plot Protein expression by chrm CN ####
### Step 2: make function to plot protein expression in mono, di and triploid cells 
## ProtExp.ChrmCN.filtered
## for Protein Data:
ProtExp.ChrmCN.filtered<- function(Protein){
  
  TestProt=Protein
  
  testProtChrm <- filter(Protein_Info4, Protein_Info4$"Approved Symbol"==TestProt) #get test protein data
  TestArm <- gsub('[0-9]+', '', testProtChrm$'Chromosome band') #Find test protein chromosome arm
  TestArm <- gsub('[.]', '', TestArm) #also remove period, if needed
  TestChrm<- as.numeric(testProtChrm$Chromosome) #Find test protein chromosome
  #filter cell data and get only aneuploidy scores for chrm arm of test protein location
  testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
  testchrm.percell <- filter(testchrm.percell, arm==TestArm)
  
  
  TestProteinInfo <- filter(Protein_ProID, Gene_Symbol==TestProt)
  TestProteinID<- TestProteinInfo$Protein_Id
  
  Protein_Expression_TestProt <- Protein.Expression.filtered %>% select(Cell_line, TestProteinID)
  
  #Combine protein data with aneuploidy data for test protein location
  Protein.effect.Chrm<-merge(y= Protein_Expression_TestProt, x= testchrm.percell, 
                             by.y="Cell_line", by.x="DepMap_ID", 
                             sort = TRUE)# 368 cells, 12757 proteins 
  
  ## put data into categories based on chrm arm number
  Chrm.tri<- filter(Protein.effect.Chrm, arm_call==1) #all cells triploid for testProt chrm
  Chrm.di<- filter(Protein.effect.Chrm, arm_call==0) #all cells diploid for testProt chrm
  Chrm.mono<- filter(Protein.effect.Chrm, arm_call==-1) #all cells mono for testProt chrm
  
  ##Make data frame with t-test info about protein expression per arm number category.
  ## Return this data at the end of the function. 
  if (length(Chrm.tri[!is.na(Chrm.tri)])>=3 & 
      #Added if statement so only genes with 2+ values per condition are analyzed
      #if I don't do this, the t-test crashes and I get no values. 
      length(Chrm.di[!is.na(Chrm.di)])>=3 &
      length(Chrm.di[!is.na(Chrm.mono)])>=3) {
    # Trisomy vs disomy
    di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with Protein data
                   Chrm.di[,8], # [,8] because that is column with Protein data
                   mu = 0, 
                   alt = "two.sided",
                   conf.level = 0.99) #get p-value
    Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
    #Disomy vs monosome
    Mono.Di<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                    Chrm.di[,8], # [,8] because that is column with Protein data
                    mu = 0, 
                    alt = "two.sided",
                    conf.level = 0.99) #get p-value
    Diff.Mono.Di<-mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE)  #get difference
    # Tri vs monosomy
    Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                     Chrm.tri[,8], # [,8] because that is column with Protein data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
    Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
    
    t.test_prot_category<- data.frame(Protein_ID=colnames(Protein.effect.Chrm[8]), 
                                      Protein_Name=TestProt,
                                      Pvlaue.Tri.Di= di.tri$p.value,
                                      Diff.Gain= Diff.Tri.Di, 
                                      Pvlaue.Loss= Mono.Di$p.value,
                                      Diff.Loss= Diff.Mono.Di, 
                                      Pvlaue.Tri.Mono= Tri.Mono$p.value,
                                      Diff.GainvLoss= Diff.Tri.Mono)
  } else {
    t.test_prot_category<- c("Not enough cells to perform t-test on all conditions") 
  }
  
  
  ##Make Plot of protein expression per aneuploid category
  setwd(DataFileLocation) 
  pdf(paste("Plot.",TestProt,".Protein.Expression.per.Tri.Di.Mono.pdf", sep=''), width=4, height=4)## Save as PDF
  Plot.testProt<- ggplot(Protein.effect.Chrm,
                         aes(x = as_factor(arm_call), y=Protein.effect.Chrm[,8])) + 
    geom_boxplot(fill= c("-1"="dodgerblue3", "0"="grey80", "1"="red"), outlier.shape = NA) +
    geom_jitter() +
    xlab(paste("Chromosome arm call of ", TestProt," location: Chrm",TestChrm, TestArm))+
    ylab("Protein Expression per cell line") +
    theme_classic()+
    ggtitle(paste(TestProt," protein expression per cell\n per corresponding chromosome copy number"))
  print(Plot.testProt)
  dev.off()##stop saving PDF
  
  return(t.test_prot_category)
}

###### Define function: plot RNA expression by chrm CN #####
### Step 3: make function to plot RNA expression in mono, di and triploid cells 
## for RNA Data:
RNAExp.ChrmCN.filtered<- function(RNA){
  
  TestRNA=RNA
  
  testRNAChrm <- filter(RNA_Info3, RNA_Info3$"Approved Symbol"==TestRNA) #get test RNA data
  TestArm <- gsub('[0-9]+', '', testRNAChrm$'Chromosome band') #Find test RNA chromosome arm
  TestArm <- gsub('[.]', '', TestArm) #also remove period, if needed
  TestChrm<- as.numeric(testRNAChrm$Chromosome) #Find test RNA chromosome
  #filter cell data and get only aneuploidy scores for chrm arm of test RNA location
  testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
  testchrm.percell <- filter(testchrm.percell, arm==TestArm)
  
  colnames(RNA.Expression.filtered)<-sub("[..].*", "", as.character(colnames(RNA.Expression.filtered)))#12755 genes
  
  RNA_Expression_TestRNA <- RNA.Expression.filtered %>% select(Cell_line, all_of(TestRNA))
  
  #Combine RNA data with aneuploidy data for test RNA location
  RNA.effect.Chrm<-merge(y= RNA_Expression_TestRNA, x= testchrm.percell, 
                         by.y="Cell_line", by.x="DepMap_ID", 
                         sort = TRUE)# 368 cells, 12757 RNAs 
  
  ## put data into categories based on chrm arm number
  Chrm.tri<- filter(RNA.effect.Chrm, arm_call==1) #all cells triploid for testRNA chrm
  Chrm.di<- filter(RNA.effect.Chrm, arm_call==0) #all cells diploid for testRNA chrm
  Chrm.mono<- filter(RNA.effect.Chrm, arm_call==-1) #all cells mono for testRNA chrm
  
  ##Make data frame with t-test info about RNA expression per arm number category.
  ## Return this data at the end of the function. 
  if (length(Chrm.tri[!is.na(Chrm.tri)])>=3 & 
      #Added if statement so only genes with 2+ values per condition are analyzed
      #if I don't do this, the t-test crashes and I get no values. 
      length(Chrm.di[!is.na(Chrm.di)])>=3 &
      length(Chrm.di[!is.na(Chrm.mono)])>=3) {
    # Trisomy vs disomy
    di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with RNA data
                   Chrm.di[,8], # [,8] because that is column with RNA data
                   mu = 0, 
                   alt = "two.sided",
                   conf.level = 0.99) #get p-value
    Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
    #Disomy vs monosome
    Mono.Di<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                    Chrm.di[,8], # [,8] because that is column with RNA data
                    mu = 0, 
                    alt = "two.sided",
                    conf.level = 0.99) #get p-value
    Diff.Mono.Di<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE)  #get difference
    # Tri vs monosomy
    Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                     Chrm.tri[,8], # [,8] because that is column with RNA data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
    Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
    
    t.test_RNA_category<- data.frame(RNA_ID=colnames(RNA.effect.Chrm[8]),
                                     RNA_Name=TestRNA,
                                     Pvlaue.Tri.Di= di.tri$p.value,
                                     Diff.Gain= Diff.Tri.Di, 
                                     Pvlaue.Loss= Mono.Di$p.value,
                                     Diff.Loss= Diff.Mono.Di, 
                                     Pvlaue.Tri.Mono= Tri.Mono$p.value,
                                     Diff.GainvLoss= Diff.Tri.Mono)
  } else {
    t.test_RNA_category<- c("Not enough cells to perform t-test on all conditions") 
  }
  
  
  ##Make Plot of RNA expression per aneuploid category
  setwd(DataFileLocation)
  pdf(paste("Plot.",TestRNA,".RNA.Expression.per.Tri.Di.Mono.pdf", sep=''), width=4, height=4)## Save as PDF
  Plot.testRNA<- ggplot(RNA.effect.Chrm,
                        aes(x = as_factor(arm_call), y=RNA.effect.Chrm[,8])) + 
    geom_boxplot(fill= c("-1"="dodgerblue3", "0"="grey80", "1"="red"), outlier.shape = NA) +
    geom_jitter() +
    xlab(paste("Chromosome arm call of ", TestRNA," location: Chrm",TestChrm, TestArm))+
    ylab("RNA Expression per cell line") +
    theme_classic()+
    ggtitle(paste(TestRNA," RNA expression per cell\n per corresponding chromosome copy number"))
  print(Plot.testRNA)
  dev.off()##stop saving PDF
  
  return(t.test_RNA_category)
}

###### Boxplots of specific gene expression per gain/neutral/loss ####

## plot specific protein expression per category
## plot and get p-values for Protein/RNA expression changes by chromosome CN
p.P53<-ProtExp.ChrmCN.filtered("TP53")
p.MDM2<-ProtExp.ChrmCN.filtered("MDM2")
p.MDM4<-ProtExp.ChrmCN.filtered("MDM4") #Significant Difference! 
p.CDK1<-ProtExp.ChrmCN.filtered("CDK1")
p.CDK13<-ProtExp.ChrmCN.filtered("CDK13") 
p.CDK12<-ProtExp.ChrmCN.filtered("CDK12")

## RNA expression
r.P53<-RNAExp.ChrmCN.filtered("TP53") #Significant Difference!
r.MDM2<-RNAExp.ChrmCN.filtered("MDM2") #Significant Difference!
r.MDM4<-RNAExp.ChrmCN.filtered("MDM4") #Significant Difference!
r.CDK1<-RNAExp.ChrmCN.filtered("CDK1") #Significant Difference!
r.CDK13<-RNAExp.ChrmCN.filtered("CDK13") #NS
r.CDK12<-RNAExp.ChrmCN.filtered("CDK12") #NS

## Found top 5 deregulated Proteins and RNA, by p-value, betwaan Mono and Tri. 
##plot Top 5 Protein
t.test_prot_category[1:5,2]
p.PTPN2<-ProtExp.ChrmCN.filtered("PTPN2")
r.PTPN2<-RNAExp.ChrmCN.filtered("PTPN2")
p.SMCHD1<-ProtExp.ChrmCN.filtered("SMCHD1")
r.SMCHD1<-RNAExp.ChrmCN.filtered("SMCHD1")
p.RNMT<-ProtExp.ChrmCN.filtered("RNMT")
r.RNMT<-RNAExp.ChrmCN.filtered("RNMT")
p.USP14<-ProtExp.ChrmCN.filtered("USP14")
r.USP14<-RNAExp.ChrmCN.filtered("USP14")
p.PPP4R1<-ProtExp.ChrmCN.filtered("PPP4R1")
r.PPP4R1<-RNAExp.ChrmCN.filtered("PPP4R1")

##plot Top 5 RNA
t.test_RNA_category[1:5,1]
p.NDUFAF5<-ProtExp.ChrmCN.filtered("NDUFAF5")
r.NDUFAF5<-RNAExp.ChrmCN.filtered("NDUFAF5")
p.NDUFV2<-ProtExp.ChrmCN.filtered("NDUFV2")
r.NDUFV2<-RNAExp.ChrmCN.filtered("NDUFV2")
p.MGME1<-ProtExp.ChrmCN.filtered("MGME1")
r.MGME1<-RNAExp.ChrmCN.filtered("MGME1")
p.ESF1<-ProtExp.ChrmCN.filtered("ESF1")
r.ESF1<-RNAExp.ChrmCN.filtered("ESF1")
p.MAVS<-ProtExp.ChrmCN.filtered("MAVS")
r.MAVS<-RNAExp.ChrmCN.filtered("MAVS")

## Plot TSG scaling w/gain loss at RNA
## Onco and TSG genes should be in Bailey lists

##Tumor supressor gene (TSG) and oncogene (OG) list
## Download Bailey et al. 2018 supplemental figure 8, upload table 1. example below: 
# Bailey et al. Comprehensive Characterization of Cancer Driver Genes and Mutations, cell, 2018

Oncogene.Bailey.list #get all oncogenes in Bailey et al 
TSG.Bailey.list # get all TSG in Bailey et al 

## Oncogenes
## Bailey
Oncog.Protein.ttest<-subset(t.test_prot_category, Protein_Name %in% Oncogene.Bailey.list)
Oncog.Protein.ttest<- Oncog.Protein.ttest[order(Oncog.Protein.ttest$Diff.Di_Mono),]
# most buffeted to most scaling (Tri and mono)--No sig, 
# Most scaling to most buffeted (Tri and mono)-- No Sig
# setwd(DataFileLocation)
# write_csv(Oncog.Protein.ttest, "Devoli.Onco.csv")
# Q1, Q4- No Sig
oncog.RNA.ttest<- subset(t.test_RNA_category, RNA_ID %in% Oncogene.Bailey.list)
oncog.RNA.ttest<- oncog.RNA.ttest[order(oncog.RNA.ttest$Diff.Gain),]
# 25 out of 39 (64%) genes not differently expressed between chrm gain or loss (protein level)
# 9 out of 37 (24%) genes not differenetly expressed between chrm gain or loss (RNA level)

p.BRAF<-ProtExp.ChrmCN.filtered("BRAF") #Sig (Most sig onco protein)
r.BRAF<-RNAExp.ChrmCN.filtered("BRAF") #Sig (tri)
p.MAPK1<-ProtExp.ChrmCN.filtered("MAPK1") #SIG
r.MAPK1<-RNAExp.ChrmCN.filtered("MAPK1") #SIG 
p.RRAS2<-ProtExp.ChrmCN.filtered("RRAS2") #SIG (Tri)
r.RRAS2<-RNAExp.ChrmCN.filtered("RRAS2") #SIG (Tri)
p.MTOR<-ProtExp.ChrmCN.filtered("MTOR") #NS
r.MTOR<-RNAExp.ChrmCN.filtered("MTOR") #SIG!
p.MYC<-ProtExp.ChrmCN.filtered("MYC") #NS
r.MYC<-RNAExp.ChrmCN.filtered("MYC") #NS

## TSG
# Bailey
#TSG.list2<-TSG.list[-c(34,35)] #Remove HLA-A and HLA-B

TSG.Protein.ttest<-subset(t.test_prot_category, Protein_Name %in% TSG.Devoli.list)
TSG.Protein.ttest<- TSG.Protein.ttest[order(TSG.Protein.ttest$Diff.Di_Mono),]
# most buffeted to most scaling (Tri and mono)-- NS, NS
# Most scaling to most buffeted (Tri and mono)-- NS, NS

TSG.RNA.ttest<- subset(t.test_RNA_category, RNA_ID %in% TSG.Bailey.list)
TSG.RNA.ttest<- TSG.RNA.ttest[order(TSG.RNA.ttest$Diff.Gain),]


# 35 of 66 genes (53%) are Not differently different expressed between gain or loss (protein level)
# 24 of 69 genes (35%) are Not differently different expressed between gain or loss (RNA level)
p.SMAD4<-ProtExp.ChrmCN.filtered("SMAD4") #Sig
r.SMAD4<-RNAExp.ChrmCN.filtered("SMAD4") #Sig
p.ATM<-ProtExp.ChrmCN.filtered("ATM") #Sig
r.ATM<-RNAExp.ChrmCN.filtered("ATM") #Sig
p.TP53<-ProtExp.ChrmCN.filtered("TP53") #NS
r.TP53<-RNAExp.ChrmCN.filtered("TP53") # Sig (mono)
p.MAP3K1<-ProtExp.ChrmCN.filtered("MAP3K1") #NS
r.MAP3K1<-RNAExp.ChrmCN.filtered("MAP3K1") #slight sig tri
p.CTNND1<-ProtExp.ChrmCN.filtered("CTNND1") #NS
r.CTNND1<-RNAExp.ChrmCN.filtered("CTNND1") #NS
p.KANSL1<-ProtExp.ChrmCN.filtered("KANSL1") #NS
r.KANSL1<-RNAExp.ChrmCN.filtered("KANSL1") #Sig
p.CDKN1A<-ProtExp.ChrmCN.filtered("CDKN1A") #NS
r.CDKN1A<-RNAExp.ChrmCN.filtered("CDKN1A") #NS
p.BRCA1<-ProtExp.ChrmCN.filtered("BRCA1") #NS
r.BRCA1<-RNAExp.ChrmCN.filtered("BRCA1") #Sig on mono
p.SMAD4<-ProtExp.ChrmCN.filtered("SMAD4") #NS
r.SMAD4<-RNAExp.ChrmCN.filtered("SMAD4") #Sig on mono

GainNS.LossNS<- subset(TSG.Protein.ttest, 
                       (Pvlaue.Tri.Di>0.05 &
                          Diff.Gain<0.2 &
                          Pvlaue.Di.Mono>0.05 &
                          Diff.Di_Mono> -0.2 ) )




### find gene increase upon gain, ns upon loss
GainSig.LossNS<- subset(CN.Diff.xRNA.yProt.ThreeGroups, 
                        (Protein.Pvalue.Gain<0.05 & #0r 0.05/9414
                           Protein.Diff.Gain>0.25 &
                           Protein.Pvalue.Loss>0.05 &
                           Protein.Diff.Loss>0 ) )
GainSig.LossNS<-GainSig.LossNS[order(GainSig.LossNS$Protein.Diff.Gain),]
p.KANSL1<-ProtExp.ChrmCN.filtered("TUBA3C") # gain sig, loss ns
r.KANSL1<-RNAExp.ChrmCN.filtered("TUBA3C") # gain sig, loss ns
p.KANSL1<-ProtExp.ChrmCN.filtered("GOLGA2") # gain sig, loss ns
r.KANSL1<-RNAExp.ChrmCN.filtered("GOLGA2") # gain sig, loss ns

FindPercent<- subset(CN.Diff.xRNA.yProt.ThreeGroups, 
                     (Protein.Pvalue.Gain<0.05 & #Sig 
                        Protein.Diff.Gain>0 & #increase
                        Protein.Pvalue.Loss<0.05 &
                        Protein.Diff.Loss<0 &
                        RNA.Pvalue.Gain<0.05 & 
                        RNA.Diff.Gain>0 &
                        RNA.Pvalue.Loss<0.05 &
                        RNA.Diff.Loss<0 ) )
(length(FindPercent$RNA_Name)*100)/length(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name)
# 19% no sig RNA nor Protein, gain and loss
# 9 % sig RNA and protein, gain and loss
# 13% sig RNA, not protein, gain and loss
# 6 % sig gain only, RNA and Protein
# 9 % sig loss only, RNA and Protein
# 7 % sig RNA, Protein gain only
# 11% sig RNA, Protein loss only
# 7 % sig RNA gain only
# 9 % sig RNA loss only
# 1.2 % sig Protein gain only
# 1.3 % sig Protein loss only
# 3.4 % Sig anti-scaling at RNA or protein, gain or loss

# 16% Protein gain sig, loss not (regardless of RNA)
# 22% Protein loss sig, gain not (regardless of RNA)


###### Make RNA/Protein difference & t-test list per Diploid, & triploid ####

### subset 367 cells into ploidy
# change this for diploid or triploid subsets: 
FilteredCellPloidyData <- subset(aneuploid, DepMap_ID %in% Protein.Expression.filtered$Cell_line)
FilteredCellPloidyData <- FilteredCellPloidyData[,c(2,6)]
FilteredCellPloidyData <- unique(FilteredCellPloidyData)

#MonoploidCells <- subset(FilteredCellPloidyData, ploidy<1.5) # 7, 2%
DiploidCells <- subset(FilteredCellPloidyData, ploidy>=1.5 & ploidy<2.5) # 185, 50%
TriploidCells <- subset(FilteredCellPloidyData, ploidy>=2.5 & ploidy<3.5) # 137, 37%
# TetraploidCells <- subset(FilteredCellPloidyData, ploidy>=3.5 & ploidy<4.5) # 24, 7%
#HighPloidyCells <- subset(FilteredCellPloidyData, ploidy>4.5) # 14, mean ploidy =5.7, 4%



### RNA difference & p-value in near-diploid, near-triploid and near-tetraploid aneuploidy
colnames(RNA.Expression.filtered)<-sub("[..].*", "", as.character(colnames(RNA.Expression.filtered)))#12755 genes
RNA.Expression.ByPloidy<-RNA.Expression.filtered[,-c(1)]

#set cell subgroup by ploidy: 
### !!! edit the below code to do Tiploid or diploid cells as needed:
RNA.Expression.ByPloidy<- subset(RNA.Expression.ByPloidy, Cell_line %in% TriploidCells$DepMap_ID)

t.test_RNA_category_byPloidy<- data.frame(RNA_ID=character(),
                                 RNA_Name=character(),
                                 Pvalue.Gain= numeric(),
                                 Diff.Gain= numeric(), 
                                 Pvalue.Loss= numeric(),
                                 Diff.Loss= numeric(),
                                 Pvalue.GainvLoss= numeric(),
                                 Diff.GainvLoss= numeric())

for (i in 2:length(RNA.Expression.ByPloidy)){
  
  TestRNA_name=sub("[..].*", "", as.character(colnames(RNA.Expression.ByPloidy[i])))
  TestRNA=colnames(RNA.Expression.ByPloidy[i])
  
  testRNAChrm <- filter(RNA_Info3, RNA_Info3$"Approved Symbol"==TestRNA_name) #get test RNA data
  if (length(testRNAChrm$Chromosome)!=0){ #Make sure RNA data is in RNA_info3, else it crashes
    TestArm <- testRNAChrm$arm #Find test RNA chromosome arm
    TestChrm<- as.numeric(testRNAChrm$Chromosome) #Find test RNA chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test RNA location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    RNA_Expression_TestRNA <- RNA.Expression.ByPloidy %>% select(Cell_line, all_of(TestRNA))
    
    #Combine RNA data with aneuploidy data for test RNA location
    RNA.effect.Chrm<-merge(y= RNA_Expression_TestRNA, x= testchrm.percell, 
                           by.y="Cell_line", by.x="DepMap_ID", 
                           sort = TRUE)# 368 cells, 12757 RNAs 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(RNA.effect.Chrm, arm_call==1) #all cells triploid for testRNA chrm
    Chrm.di<- filter(RNA.effect.Chrm, arm_call==0) #all cells diploid for testRNA chrm
    Chrm.mono<- filter(RNA.effect.Chrm, arm_call==-1) #all cells mono for testRNA chrm
    
    ##Make data frame with t-test info about RNA expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 &  
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 & 
        colSums(!is.na(Chrm.mono[8]))>=10) { 
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with RNA data
                     Chrm.di[,8], # [,8] because that is column with RNA data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Mono.Di<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                      Chrm.di[,8], # [,8] because that is column with RNA data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Mono.Di<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                       Chrm.tri[,8], # [,8] because that is column with RNA data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_RNA_category_byPloidy<- rbind(t.test_RNA_category_byPloidy, data.frame(
        RNA_ID=colnames(RNA.effect.Chrm[8]),
        RNA_Name=TestRNA_name,
        Pvalue.Gain= di.tri$p.value,
        Diff.Gain= Diff.Tri.Di, 
        Pvalue.Loss= Mono.Di$p.value,
        Diff.Loss= Diff.Mono.Di, 
        Pvalue.GainvLoss= Tri.Mono$p.value,
        Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_RNA_category_byPloidy<-t.test_RNA_category_byPloidy
    }
  }
}
t.test_RNA_category_byPloidy<-distinct(t.test_RNA_category_byPloidy) 
length(t.test_RNA_category_byPloidy$RNA_ID)

### minimum of 10 cells/category:  
# diploid: 10841 RNA ( genes in all data)
# triploid: 16434 RNA 
t.test_RNA_category_byPloidy<-t.test_RNA_category_byPloidy[order(t.test_RNA_category_byPloidy$Pvalue.GainvLoss),]

### if you run the above code with a minimum of 3 cells/condition, you get this data: 
# diploid: 10841 RNA (17774 genes in all data)
# triploid: 16434 RNA 
# tetraploid: 9579 RNA (min 3)


###
### Protein difference & p-value in near-diploid, near-triploid and near-tetraploid aneuploidy 
## make list of pvalue and difference between mono- di-tri cells per gene 
## for each gene in Protein data. using only cells that have both RNA and Protein data

## !! edit the below code as needed to get diploid or triploid data:
Protein.Expression.byPloidy<- subset(Protein.Expression.filtered, 
                                     Cell_line %in% TriploidCells$DepMap_ID)

Protein.ID<-data.frame(Protein_ID=
                         colnames(Protein.Expression.filtered[3:length(Protein.Expression.filtered)]))

##Now merge depmap column name info with Protein info. 
## First add Uniprot ID and gene symbol info to depmap data. 
Depmap.Protein.info2<-merge(x= Protein.ID, y= Protein_ProID, 
                            by.x="Protein_ID", by.y="Protein_Id", 
                            sort = TRUE)# 3 collumns, 12755. all Proteins given Uniprot IDs.
# Now combine all genes with same uniprot ID: 
Depmap.Protein.info3<-merge(x= Depmap.Protein.info2, y= Protein_Info4, 
                            by.x="Uniprot_Acc", by.y="UniProt ID(supplied by UniProt)", 
                            sort = TRUE)# 10594 Proteins
# Find those genes without uniprot id: 
No_UniprotID <- anti_join(Depmap.Protein.info2, Protein_Info4, #finding genes_Symbol without match
                          by = c("Uniprot_Acc" = "UniProt ID(supplied by UniProt)"))
#...find the genes with matching gene symbols
Depmap.Protein.info4<-merge(x= No_UniprotID, y= Protein_Info4, 
                            by.x="Gene_Symbol", by.y="Approved Symbol", 
                            sort = TRUE)# 10 collumns, 253 genes, only no gene-symbol genes with uniprot_ID

# Merge genes with gene symbol and genes with only uniprot ID. 
Depmap.Protein.info5<-merge(x= Depmap.Protein.info3, y= Depmap.Protein.info4, 
                            all=TRUE)# 11 collumns, 12100 genes


t.test_prot_category_byPloidy<- data.frame(Protein_ID=character(),
                                  Protein_Name=character(),
                                  Pvlaue.Tri.Di= numeric(),
                                  Diff.Gain= numeric(), 
                                  Pvlaue.Di.Mono= numeric(),
                                  Diff.Di_Mono= numeric(),
                                  Pvlaue.Tri.Mono= numeric(),
                                  Diff.GainvLoss= numeric())


for (i in 3:length(Protein.Expression.byPloidy)){
  TestProt = colnames(Protein.Expression.byPloidy[i])
  
  testProtChrm <- filter(Depmap.Protein.info5, Protein_ID==TestProt) #get test protein data
  
  
  #use If statement to check that Protein Info has Protein of interest. 
  if (length(testProtChrm[,1])!=0){
    TestArm <- testProtChrm$arm #Find test protein chromosome arm
    TestChrm<- as.numeric(testProtChrm$Chromosome) #Find test protein chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test protein location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    Protein_Expression_TestProt <- Protein.Expression.byPloidy %>% select(Cell_line, all_of(TestProt))
    
    #Combine protein data with aneuploidy data for test protein location
    Protein.effect.Chrm<-merge(y= Protein_Expression_TestProt, x= testchrm.percell, 
                               by.y="Cell_line", by.x="DepMap_ID", 
                               sort = TRUE)# 368 cells, 12757 proteins 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(Protein.effect.Chrm, arm_call==1) #all cells triploid for testProt chrm
    Chrm.di<- filter(Protein.effect.Chrm, arm_call==0) #all cells diploid for testProt chrm
    Chrm.mono<- filter(Protein.effect.Chrm, arm_call==-1) #all cells mono for testProt chrm
    
    ##Make data frame with t-test info about protein expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 & #[8] is column with protein expression data. check it has values. 
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 & #usually 10, now 3**
        colSums(!is.na(Chrm.mono[8]))>=10) { #usually 10, now 3**
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with Protein data
                     Chrm.di[,8], # [,8] because that is column with Protein data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Di.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                      Chrm.di[,8], # [,8] because that is column with Protein data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Di.Mono<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE)  #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                       Chrm.tri[,8], # [,8] because that is column with Protein data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_prot_category_byPloidy<- rbind(t.test_prot_category_byPloidy, 
                                   data.frame(Protein_ID=TestProt,
                                              Protein_Name=testProtChrm$Gene_Symbol,
                                              Pvlaue.Tri.Di= di.tri$p.value,
                                              Diff.Gain= Diff.Tri.Di, 
                                              Pvlaue.Di.Mono= Di.Mono$p.value,
                                              Diff.Di_Mono= Diff.Di.Mono, 
                                              Pvlaue.Tri.Mono= Tri.Mono$p.value,
                                              Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_prot_category_byPloidy<-t.test_prot_category_byPloidy
    }
  }
}

t.test_prot_category_byPloidy<-distinct(t.test_prot_category_byPloidy) 
## minimum of 10 cells/category: 
# diploid=  4685 proteins
# triploid = 7738 proteins
t.test_prot_category_byPloidy<-t.test_prot_category_byPloidy[order(t.test_prot_category_byPloidy$Pvlaue.Tri.Mono),]# *** genes

## if you run this code with a miminum of 3 cells/condition instead of 10, you get this data:
# diploid= 4685 proteins
# triploid = 7738 proteins
# tetraploid = 4147 proteins



######          Di-Triploid: Combine RNA and Protein Data, save and calculate diff per ploidy group ####
## This is with filtered cell lines (only cells with RNA and Protein and aneuploid data:)
## This is with filtered genes: only genes with 10+ cells in all categories: RNA & protein gain, neutral and loss
setwd(DataFileLocation)

CN.Diff.RNA.Prot_byPloidy<-merge(x=t.test_RNA_category_byPloidy, 
                                 y=t.test_prot_category_byPloidy, 
                                 by.x="RNA_Name", by.y="Protein_Name") #4147

colnames(CN.Diff.RNA.Prot_byPloidy)<-c("RNA_Name", "RNA_ID", 
                                "RNA.Pvalue.Gain", "RNA.Diff.Gain", 
                                "RNA.Pvalue.Loss", "RNA.Diff.Loss", 
                                "RNA.Pvalue.GainvLoss", "RNA.Diff.GainvLoss", 
                                "Protein_ID", "Protein.Pvalue.Gain", 
                                "Protein.Diff.Gain", "Protein.Pvalue.Loss", 
                                "Protein.Diff.Loss", "Protein.Pvalue.GainvLoss", 
                                "Protein.Diff.GainvLoss")

write.csv(CN.Diff.RNA.Prot_byPloidy, #change name as needed ! for triploid diploid
          file =paste("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_Triploid_min10cells.csv", sep=','), 
          row.names = TRUE)


### Save diploid and triploid datasets: 
### Diploid 
## run above code with diploid cells dataset, then save as "diploid"
CN.Diff.RNA.Prot_Diploid<-CN.Diff.RNA.Prot_byPloidy

#CN.Diff.RNA.Prot_Diploid<- read.csv("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_Diploid_min10cells.csv")


### Diploid: 
mean(CN.Diff.RNA.Prot_Diploid$RNA.Diff.Gain)
mean(CN.Diff.RNA.Prot_Diploid$Protein.Diff.Gain)
mean(CN.Diff.RNA.Prot_Diploid$RNA.Diff.Loss)
mean(CN.Diff.RNA.Prot_Diploid$Protein.Diff.Loss)

# now get percent change, not log2 fold change! 
100*2^(mean(CN.Diff.RNA.Prot_Diploid$RNA.Diff.Gain))-100
100*2^(mean(CN.Diff.RNA.Prot_Diploid$Protein.Diff.Gain))-100
100*2^(mean(CN.Diff.RNA.Prot_Diploid$RNA.Diff.Loss))-100
100*2^(mean(CN.Diff.RNA.Prot_Diploid$Protein.Diff.Loss))-100

# cells= 185
# genes (>10 points/category shared in diploid and tridploid)= 4489
# Mean Prot gain: 0.188957, 14%
# Mean RNA  gain: 0.328812, 26%
# Mean Prot loss: -0.1717, -11%
# Mean RNA  loss: -0.3156, -19%

length(CN.Diff.RNA.Prot_Diploid$Protein.Diff.Gain)# 4489
sum(CN.Diff.RNA.Prot_Diploid$Protein.Diff.Gain >= log2(3/2)) #318, 7.1% >= DNA CN change upon gain
sum(CN.Diff.RNA.Prot_Diploid$Protein.Diff.Gain <= log2(1/2)) #5, 0.1% <= DNA CN change upon loss


### Triploid 
## run above code with diploid cells dataset, then save as "diploid"
CN.Diff.RNA.Prot_Triploid<-CN.Diff.RNA.Prot_byPloidy

CN.Diff.RNA.Prot_Triploid<- read.csv("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_Triploid_min10cells.csv")

### Triploid : 
mean(CN.Diff.RNA.Prot_Triploid$RNA.Diff.Gain)
mean(CN.Diff.RNA.Prot_Triploid$Protein.Diff.Gain)
mean(CN.Diff.RNA.Prot_Triploid$RNA.Diff.Loss)
mean(CN.Diff.RNA.Prot_Triploid$Protein.Diff.Loss)

# cells= 137
# genes (>10 points/category)= 7411
# Mean Prot gain:  0.134, 9.7%
# Mean RNA  gain:  0.245, 19%
# Mean Prot loss: -0.151, -10%
# Mean RNA  loss: -0.288, -18%

# number of genes >= log2(4/3) expression (>= DNA CN change) **
length(CN.Diff.RNA.Prot_Triploid$Protein.Diff.Gain) #7411 proteins
sum(CN.Diff.RNA.Prot_Triploid$Protein.Diff.Gain >= log2(4/3)) #801, 11% >= DNA CN change upon gain
sum(CN.Diff.RNA.Prot_Triploid$Protein.Diff.Gain <= log2(2/3)) #82, 1.1% <= DNA CN change upon loss



### All 
 CN.Diff.xRNA.yProt.ThreeGroups
 100*2^(mean(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain))-100
 100*2^(mean(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain))-100
 100*2^(mean(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss))-100
 100*2^(mean(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss))-100

# All 
 #### NOTE!!! log2 fold change difference is not the same as percent! calculate back to percent
 ## percent change = 2^(log2 fold change)-1
# cells= 367
# genes (>10 points/category)= 9414
# Mean Prot gain: 12%
# Mean RNA  gain: 22%
# Mean Prot loss: -8.4%
# Mean RNA  loss: -15%


# Diploid (min 10/category)
# cells= 185
# genes (>10 points/category shared in diploid and tridploid)= 4489
# Mean Prot gain: 0.188957, 14%
# Mean RNA  gain: 0.328812, 26%
# Mean Prot loss: -0.1717, -11%
# Mean RNA  loss: -0.3156, -19%

# Triploid (min 10/category)
# cells= 137
# genes (>10 points/category)= 7411
# Mean Prot gain:  0.134, 9.7%
# Mean RNA  gain:  0.245, 19%
# Mean Prot loss: -0.151, -10%
# Mean RNA  loss: -0.288, -18%


####Heatmap for mean expression upon gain and loss, protein/RNA
#calculate mean gain loss in DNA copy number upon chrm gain and loss, account for ploidy
log2(3/2)*(185/367) + log2(4/3)*(137/367) + log2(5/4)*(24/367) + log2(2/1)*(7/367) + log2(6.7/5.7)*(14/367)
log2(1/2)*(185/367) + log2(2/3)*(137/367) + log2(3/4)*(24/367) + log2(4.7/5.7)*(14/367)


#make dataframe for heatmap
Di.Tri.all.meanDiff<-data.frame(
  Cells=as.factor(c("Diploid", "Diploid", "Diploid", "Diploid", "Diploid", "Diploid", 
                       "Triploid", "Triploid", "Triploid", "Triploid", "Triploid", "Triploid", 
                       #"Tetraploid", "Tetraploid", "Tetraploid", "Tetraploid", "Tetraploid", "Tetraploid", 
                       "All", "All", "All", "All", "All", "All")),
  Condition=as.factor(c("Gain.DNA","Gain.RNA", "Gain.Protein", "Loss.DNA", "Loss.RNA", "Loss.Protein", 
                           "Gain.DNA","Gain.RNA", "Gain.Protein", "Loss.DNA", "Loss.RNA", "Loss.Protein", 
                           #"Gain.DNA","Gain.RNA", "Gain.Protein", "Loss.DNA", "Loss.RNA", "Loss.Protein", 
                           "Gain.DNA","Gain.RNA", "Gain.Protein", "Loss.DNA", "Loss.RNA", "Loss.Protein")), 
  values=as.numeric(c(log2(3/2), mean(CN.Diff.RNA.Prot_Diploid$RNA.Diff.Gain), mean(CN.Diff.RNA.Prot_Diploid$Protein.Diff.Gain), 
                      log2(1/2), mean(CN.Diff.RNA.Prot_Diploid$RNA.Diff.Loss), mean(CN.Diff.RNA.Prot_Diploid$Protein.Diff.Loss), 
                      log2(4/3), mean(CN.Diff.RNA.Prot_Triploid$RNA.Diff.Gain), mean(CN.Diff.RNA.Prot_Triploid$Protein.Diff.Gain), 
                      log2(2/3), mean(CN.Diff.RNA.Prot_Triploid$RNA.Diff.Loss), mean(CN.Diff.RNA.Prot_Triploid$Protein.Diff.Loss), 
                      #log2(5/4), mean(CN.Diff.RNA.Prot_Tetraploid$RNA.Diff.Gain), mean(CN.Diff.RNA.Prot_Tetraploid$Protein.Diff.Gain), 
                      #log2(3/4), mean(CN.Diff.RNA.Prot_Tetraploid$RNA.Diff.Loss), mean(CN.Diff.RNA.Prot_Tetraploid$Protein.Diff.Loss), 
                      0.4988, mean(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain), mean(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain), 
                      -0.76, mean(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss), mean(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss) 
  )
  )
)

Di.Tri.all.meanDiff$Condition<- factor(Di.Tri.all.meanDiff$Condition, levels=c("Loss.Protein", "Loss.RNA", "Loss.DNA", "Gain.Protein", "Gain.RNA", "Gain.DNA"))
Di.Tri.all.meanDiff$Cells<- factor(Di.Tri.all.meanDiff$Cells, levels=c("Diploid", "Triploid", "Tetraploid", "All"))

#plot heatmap of mean difference by ploidy: 
ggplot(Di.Tri.all.meanDiff, aes(x=Cells, y=Condition))+
  geom_raster(aes(fill = values), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  xlab("Cancer cell line ploidy")+
  ylab("Gene expression difference")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-1, 1)) +
  ggtitle("Mean difference per cell ploidy")
# 5x4
# plot.heatmap.meanDiff.Ploidy.Min10
# plot.heatmap.meanDiff.Ploidy_withTetra

####
### now merge the diploid and triploid datasets by protein name
## and do a paired t-test between RNA and Protein expression
## see if gaining a ploidy significantly buffers gene overexpression upon gain and loss
## in both RNA and Protein
colnames(CN.Diff.RNA.Prot_Triploid)
colnames(CN.Diff.RNA.Prot_Diploid)

Merged_CN.Diff.Di_Triploid<-merge(CN.Diff.RNA.Prot_Diploid, CN.Diff.RNA.Prot_Triploid, 
                                  by.x="Protein_ID", by.y="Protein_ID")
## length 4175
## we do indeed only have samples that are present in both. Good. 

#now do paired t-tests
## percent change = 2^(log2 fold change)-1 * 100%
# x= diploid cells only
# y= triploid cells only

## RNA Diff Gain
t.test(paired=TRUE, Merged_CN.Diff.Di_Triploid$RNA.Diff.Gain.x, Merged_CN.Diff.Di_Triploid$RNA.Diff.Gain.y)
# mean diff =0.07826104, p-value <2.2E-16
# mean percent diff = 5.574473 %

## RNA Diff Loss
t.test(paired=TRUE, Merged_CN.Diff.Di_Triploid$RNA.Diff.Loss.x, Merged_CN.Diff.Di_Triploid$RNA.Diff.Loss.y)
# mean diff = -0.01046477, p-value 0.06568 NS
# mean percent diff = -0.7279997 % 


## Protein Diff Gain
t.test(paired=TRUE, Merged_CN.Diff.Di_Triploid$Protein.Diff.Gain.x, Merged_CN.Diff.Di_Triploid$Protein.Diff.Gain.y)
# mean diff =0.04656913, p-value <2.2E-16
# mean percent diff = 3.280589 %

## Protein Diff Loss
t.test(paired=TRUE, Merged_CN.Diff.Di_Triploid$Protein.Diff.Loss.x, Merged_CN.Diff.Di_Triploid$Protein.Diff.Loss.y)
# mean diff = -0.01456149, p-value 0.003147
# mean percent diff = -1.014436 %

### now get percent of diploid-only and triploid-only protein differences that are buffered: 
colnames(CN.Diff.RNA.Prot_Triploid) #7411
colnames(CN.Diff.RNA.Prot_Diploid) #4489

colnames(Merged_CN.Diff.Di_Triploid) # x= diploid, y=triploid

sum(Merged_CN.Diff.Di_Triploid$Protein.Diff.Gain.x> Merged_CN.Diff.Di_Triploid$Protein.Diff.Gain.y)
#2384 genes have larger log2FC in diploid than triploid upon chrm gain
# 2384/4175 = 57.1% of genes have greater difference in diploid

sum(Merged_CN.Diff.Di_Triploid$Protein.Diff.Loss.x < Merged_CN.Diff.Di_Triploid$Protein.Diff.Loss.y)
#2161 genes have lower (more negative) log2FC in diploid than triploid upon chrm gain
# 2161/4175 = 51.8% of genes have more negative diff in diploid


## chrm GAIN in diploid and triploid only cells: (PROTEIN)
#diploid
length(subset(Merged_CN.Diff.Di_Triploid, Protein.Diff.Gain.x >-0.1 & Protein.Diff.Gain.x <0.25)$Protein.Diff.Gain.x)
# 2192 buffered
# 2192/4175 = 52.5%

#Triploid
length(subset(Merged_CN.Diff.Di_Triploid, Protein.Diff.Gain.y >-0.1 & Protein.Diff.Gain.y <0.25)$Protein.Diff.Gain.y)
# 2440 buffered
# 2440/4175 = 58.4%


## chrm LOSS in diploid and triploid only cells: (PROTEIN)
#Diploid:
length(subset(Merged_CN.Diff.Di_Triploid, Protein.Diff.Loss.x >-0.25 & Protein.Diff.Loss.x <0.1)$Protein.Diff.Loss.x)
# 2147 buffered
# 2147/4175 = 51.4%

#Triploid: 
length(subset(Merged_CN.Diff.Di_Triploid, Protein.Diff.Loss.y >-0.25 & Protein.Diff.Loss.y <0.1)$Protein.Diff.Loss.y)
# 2715 buffered
# 2715/4175 = 65.0%


###### Low/high Aneuploid cells only (still see buffering of ribosomes?) ####

### Protein difference & p-value in near-diploid, near-triploid and near-tetraploid aneuploidy 
## make list of pvalue and difference between mono- di-tri cells per gene 
## for each gene in Protein data. using only cells that have both RNA and Protein data

### sub-step 1: Get protein and RNA expression data set for only cell lines with low aneuploidy: 

#List of low aneuploid cells from protein_expression_data_GeneScore.R
#quantile(Protein_Scores$Gene_ploidy_Score) #0, 1773, 2668, 3309, 6985
#low.aneuploid<- subset(Protein_Scores, Gene_ploidy_Score<= 1773) #only 94 cells
#low.aneuploid.cells<-low.aneuploid$Broad_ID
#high.aneuploid<- subset(Protein_Scores, Gene_ploidy_Score> 3309) #only 94 cells
#high.aneuploid.cells<-high.aneuploid$Broad_ID


Protein.Expression_lowaneuploidCells<- subset(Protein.Expression.filtered, Cell_line %in% low.aneuploid.cells)
#Protein.Expression_lowaneuploidCells<- subset(Protein.Expression.filtered, Cell_line %in% high.aneuploid.cells)

RNA.Expression_lowaneuploidCells<- subset(Protein.Expression.filtered, Cell_line %in% low.aneuploid.cells)
#RNA.Expression_lowaneuploidCells<- subset(Protein.Expression.filtered, Cell_line %in% high.aneuploid.cells)



### sub-step 2: get protein difference in cells with low aneuploidy
Protein.ID<-data.frame(Protein_ID=
                         colnames(Protein.Expression_lowaneuploidCells[3:length(Protein.Expression_lowaneuploidCells)]))

# Now merge depmap column name info with Protein info. 
# First add Uniprot ID and gene symbol info to depmap data. 
Depmap.Protein.info2<-merge(x= Protein.ID, y= Protein_ProID, 
                            by.x="Protein_ID", by.y="Protein_Id", 
                            sort = TRUE)# 3 collumns, 12755. all Proteins given Uniprot IDs.
# Now combine all genes with same uniprot ID: 
Depmap.Protein.info3<-merge(x= Depmap.Protein.info2, y= Protein_Info4, 
                            by.x="Uniprot_Acc", by.y="UniProt ID(supplied by UniProt)", 
                            sort = TRUE)# 10594 Proteins
# Find those genes without uniprot id: 
No_UniprotID <- anti_join(Depmap.Protein.info2, Protein_Info4, #finding genes_Symbol without match
                          by = c("Uniprot_Acc" = "UniProt ID(supplied by UniProt)"))
#...find the genes with matching gene symbols
Depmap.Protein.info4<-merge(x= No_UniprotID, y= Protein_Info4, 
                            by.x="Gene_Symbol", by.y="Approved Symbol", 
                            sort = TRUE)# 10 collumns, 253 genes, only no gene-symbol genes with uniprot_ID

# Merge genes with gene symbol and genes with only uniprot ID. 
Depmap.Protein.info5<-merge(x= Depmap.Protein.info3, y= Depmap.Protein.info4, 
                            all=TRUE)# 11 collumns, 12100 genes


t.test_prot_category_LowPloidy<- data.frame(Protein_ID=character(),
                                           Protein_Name=character(),
                                           Pvlaue.Tri.Di= numeric(),
                                           Diff.Gain= numeric(), 
                                           Pvlaue.Di.Mono= numeric(),
                                           Diff.Di_Mono= numeric(),
                                           Pvlaue.Tri.Mono= numeric(),
                                           Diff.GainvLoss= numeric())


for (i in 3:length(Protein.Expression_lowaneuploidCells)){
  TestProt = colnames(Protein.Expression_lowaneuploidCells[i])
  
  testProtChrm <- filter(Depmap.Protein.info5, Protein_ID==TestProt) #get test protein data
  
  
  #use If statement to check that Protein Info has Protein of interest. 
  if (length(testProtChrm[,1])!=0){
    TestArm <- testProtChrm$arm #Find test protein chromosome arm
    TestChrm<- as.numeric(testProtChrm$Chromosome) #Find test protein chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test protein location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    Protein_Expression_TestProt <- Protein.Expression_lowaneuploidCells %>% select(Cell_line, all_of(TestProt))
    
    #Combine protein data with aneuploidy data for test protein location
    Protein.effect.Chrm<-merge(y= Protein_Expression_TestProt, x= testchrm.percell, 
                               by.y="Cell_line", by.x="DepMap_ID", 
                               sort = TRUE)# 368 cells, 12757 proteins 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(Protein.effect.Chrm, arm_call==1) #all cells triploid for testProt chrm
    Chrm.di<- filter(Protein.effect.Chrm, arm_call==0) #all cells diploid for testProt chrm
    Chrm.mono<- filter(Protein.effect.Chrm, arm_call==-1) #all cells mono for testProt chrm
    
    ##Make data frame with t-test info about protein expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=3 & #[8] is column with protein expression data. check it has values. 
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=3 & #usually 10, now 3**
        colSums(!is.na(Chrm.mono[8]))>=3) { #usually 10, now 3**
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with Protein data
                     Chrm.di[,8], # [,8] because that is column with Protein data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Di.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                      Chrm.di[,8], # [,8] because that is column with Protein data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Di.Mono<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE)  #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                       Chrm.tri[,8], # [,8] because that is column with Protein data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_prot_category_LowPloidy<- rbind(t.test_prot_category_LowPloidy, 
                                            data.frame(Protein_ID=TestProt,
                                                       Protein_Name=testProtChrm$Gene_Symbol,
                                                       Pvlaue.Tri.Di= di.tri$p.value,
                                                       Diff.Gain= Diff.Tri.Di, 
                                                       Pvlaue.Di.Mono= Di.Mono$p.value,
                                                       Diff.Di_Mono= Diff.Di.Mono, 
                                                       Pvlaue.Tri.Mono= Tri.Mono$p.value,
                                                       Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_prot_category_LowPloidy<-t.test_prot_category_LowPloidy
    }
  }
}

t.test_prot_category_LowPloidy<-distinct(t.test_prot_category_LowPloidy) 

t.test_prot_category_LowPloidy<-t.test_prot_category_LowPloidy[order(t.test_prot_category_LowPloidy$Pvlaue.Tri.Mono),]# *** genes

# show ribosomes still buffered
# compare low aneuploid and high aneuploid cell lines difference expressions are correlated. 



### sub-step 3: get RNA difference in cells with low aneuploidy
# Low ploidy analysis: RNA expression difference: 
# prep data frame: 
colnames(RNA.Expression.filtered)<-sub("[..].*", "", as.character(colnames(RNA.Expression.filtered)))#12755 genes
RNA.Expression.LowPloidy<-RNA.Expression.filtered[,-c(1)]

RNA.Expression.LowPloidy<- subset(RNA.Expression.LowPloidy, Cell_line %in% high.aneuploid.cells)

# Run loop to get difference in experssion per gene in only low aneuploid cells
t.test_RNA_category_LowPloidy<- data.frame(RNA_ID=character(),
                                          RNA_Name=character(),
                                          Pvalue.Gain= numeric(),
                                          Diff.Gain= numeric(), 
                                          Pvalue.Loss= numeric(),
                                          Diff.Loss= numeric(),
                                          Pvalue.GainvLoss= numeric(),
                                          Diff.GainvLoss= numeric())

for (i in 2:length(RNA.Expression.LowPloidy)){
  
  TestRNA_name=sub("[..].*", "", as.character(colnames(RNA.Expression.LowPloidy[i])))
  TestRNA=colnames(RNA.Expression.LowPloidy[i])
  
  testRNAChrm <- filter(RNA_Info3, RNA_Info3$"Approved Symbol"==TestRNA_name) #get test RNA data
  if (length(testRNAChrm$Chromosome)!=0){ #Make sure RNA data is in RNA_info3, else it crashes
    TestArm <- testRNAChrm$arm #Find test RNA chromosome arm
    TestChrm<- as.numeric(testRNAChrm$Chromosome) #Find test RNA chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test RNA location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    RNA_Expression_TestRNA <- RNA.Expression.LowPloidy %>% select(Cell_line, all_of(TestRNA))
    
    #Combine RNA data with aneuploidy data for test RNA location
    RNA.effect.Chrm<-merge(y= RNA_Expression_TestRNA, x= testchrm.percell, 
                           by.y="Cell_line", by.x="DepMap_ID", 
                           sort = TRUE)# 368 cells, 12757 RNAs 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(RNA.effect.Chrm, arm_call==1) #all cells triploid for testRNA chrm
    Chrm.di<- filter(RNA.effect.Chrm, arm_call==0) #all cells diploid for testRNA chrm
    Chrm.mono<- filter(RNA.effect.Chrm, arm_call==-1) #all cells mono for testRNA chrm
    
    ##Make data frame with t-test info about RNA expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=3 &  
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=3 & 
        colSums(!is.na(Chrm.mono[8]))>=3) { 
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with RNA data
                     Chrm.di[,8], # [,8] because that is column with RNA data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Mono.Di<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                      Chrm.di[,8], # [,8] because that is column with RNA data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Mono.Di<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                       Chrm.tri[,8], # [,8] because that is column with RNA data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_RNA_category_LowPloidy<- rbind(t.test_RNA_category_LowPloidy, data.frame(
        RNA_ID=colnames(RNA.effect.Chrm[8]),
        RNA_Name=TestRNA_name,
        Pvalue.Gain= di.tri$p.value,
        Diff.Gain= Diff.Tri.Di, 
        Pvalue.Loss= Mono.Di$p.value,
        Diff.Loss= Diff.Mono.Di, 
        Pvalue.GainvLoss= Tri.Mono$p.value,
        Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_RNA_category_LowPloidy<-t.test_RNA_category_LowPloidy
    }
  }
}

t.test_RNA_category_LowPloidy<-distinct(t.test_RNA_category_LowPloidy) 
t.test_RNA_category_LowPloidy<-t.test_RNA_category_LowPloidy[order(t.test_RNA_category_LowPloidy$Pvalue.GainvLoss),]


### Sub step 4: combine low ploidy RNA and protein data
CN.Diff.RNA.Prot_HighPloidy<-merge(x=t.test_RNA_category_LowPloidy, 
                                 y=t.test_prot_category_LowPloidy, 
                                 by.x="RNA_Name", by.y="Protein_Name") #4147

colnames(CN.Diff.RNA.Prot_HighPloidy)<-c("RNA_Name", "RNA_ID", 
                                       "RNA.Pvalue.Gain", "RNA.Diff.Gain", 
                                       "RNA.Pvalue.Loss", "RNA.Diff.Loss", 
                                       "RNA.Pvalue.GainvLoss", "RNA.Diff.GainvLoss", 
                                       "Protein_ID", "Protein.Pvalue.Gain", 
                                       "Protein.Diff.Gain", "Protein.Pvalue.Loss", 
                                       "Protein.Diff.Loss", "Protein.Pvalue.GainvLoss", 
                                       "Protein.Diff.GainvLoss")

## add categories based on cutoffs.  (-Inf,-0.1,0.25,Inf)
CN.Diff.RNA.Prot_HighPloidy$Three.RNA.Gain<- cut(CN.Diff.RNA.Prot_HighPloidy$RNA.Diff.Gain,
                                                    breaks=c(-Inf,-0.1,0.25,Inf),
                                                    include.lowest=TRUE,
                                                    labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_HighPloidy$Three.RNA.Loss<- cut(CN.Diff.RNA.Prot_HighPloidy$RNA.Diff.Loss,
                                                    breaks=c(-Inf,-0.25,0.1,Inf),
                                                    include.lowest=TRUE,
                                                    labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.RNA.Prot_HighPloidy$Three.Protein.Gain<- cut(CN.Diff.RNA.Prot_HighPloidy$Protein.Diff.Gain,
                                                        breaks=c(-Inf,-0.1,0.25,Inf),
                                                        include.lowest=TRUE,
                                                        labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_HighPloidy$Three.Protein.Loss<- cut(CN.Diff.RNA.Prot_HighPloidy$Protein.Diff.Loss,
                                                        breaks=c(-Inf,-0.25,0.1,Inf),
                                                        include.lowest=TRUE,
                                                        labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.RNA.Prot_HighPloidy$Three.Protein.Loss<-factor(CN.Diff.RNA.Prot_HighPloidy$Three.Protein.Loss, 
                                                          levels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_HighPloidy$Three.RNA.Loss<-factor(CN.Diff.RNA.Prot_HighPloidy$Three.RNA.Loss, 
                                                      levels=c("Anti-Scaling","Buffering","Scaling"))


#setwd(DataFileLocation)
write.csv(CN.Diff.RNA.Prot_LowPloidy, #change name as needed
          file =paste("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_LowPloidy.min3cells.csv", sep=','), 
          row.names = TRUE)

#CN.Diff.RNA.Prot_LowPloidy<- read.csv("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_LowPloidy.min10cells.csv")
#CN.Diff.RNA.Prot_LowPloidy<-CN.Diff.RNA.Prot_LowPloidy[,-1]

### Repeat above code with "High ploidy cells" to get difference upon gain in the top 1/4th high aneuploid cells 
write.csv(CN.Diff.RNA.Prot_HighPloidy, #change name as needed
          file =paste("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_HighPloidy.min3cells.csv", sep=','), 
          row.names = TRUE)

#write.csv(Protein.effect.Chrm$DepMap_ID, #save a list of the cell lines used in analysis. 
#          file =paste("CellLines_ProdeinDosageCompensationManuscript.csv", sep=','), 
#          row.names = TRUE)

###### Control dataset: . NSCLC lung cancer cells only ####
# the dataset this section generates is available in supplementary data 3, sheet 8 "NSCLC_Only" 


### Protein difference & p-value in lung cancer cells only 
## make list of pvalue and difference for lung cancer cells per gene 
## for each gene in Protein data. using only cells that have both RNA and Protein data

### sub-step 1: Get protein and RNA expression data set for only lung cancer cells: 
# top cancer types: NSCLC (158), melanoma (102), glioma (90), colorectal_adenocarcinoma (81)

Cell_Line_Info <- read_delim(file="sample_info.csv", 
                             delim=",")
lung.cells <- subset(Cell_Line_Info, lineage_subtype=="NSCLC")

### Make quick pie graph of which cell line lineages
Cell_Lines_Used<- subset(Cell_Line_Info, DepMap_ID %in% Protein.Expression.filtered_min10Cells$Cell_line)
Cell_Lines_Used<- subset(Cell_Lines_Used, DepMap_ID %in% aneuploid$DepMap_ID)

Table1<- as.data.frame(table(Cell_Lines_Used$lineage_subtype))
Table1<- Table1[order(-Table1$Freq),]
Table2<- Table1[1:6,]
levels(Table2$Var1) <- c(levels(Table2$Var1), "Other")
Table2[7,]<- c("Other", 192)
Table2$Freq<-as.numeric(Table2$Freq)
Table2$Var1 <- factor(Table2$Var1, levels = Table2$Var1)

ggplot(Table2, aes(x="", y=Freq, fill=Var1))+
  geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0)+
  theme_void()
# plot.pie.CellLine_Subtypes
# 6x4

# get protein data only for proteins with minimum 10 cells with protein data
Protein.Expression.filtered_min10Cells<- Protein.Expression.filtered %>% select(one_of(CN.Diff.xRNA.yProt.ThreeGroups$Protein_ID))
# Now add Cell_lines back in. 
Protein.Expression.filtered_min10Cells$Cell_line<- Protein.Expression.filtered$Cell_line

Protein.Expression_Lung<- subset(Protein.Expression.filtered_min10Cells, Cell_line %in% lung.cells$DepMap_ID)
# number of cells= 64
# number of genes 9414

### sub-step 2: get protein difference in lung cells 
Protein.ID<-data.frame(Protein_ID=
                         colnames(Protein.Expression_Lung[3:length(Protein.Expression_Lung)]))

# Now merge depmap column name info with Protein info. 
# First add Uniprot ID and gene symbol info to depmap data. 
Depmap.Protein.info2<-merge(x= Protein.ID, y= Protein_ProID, 
                            by.x="Protein_ID", by.y="Protein_Id", 
                            sort = TRUE)# 3 collumns, 12755. all Proteins given Uniprot IDs.
# Now combine all genes with same uniprot ID: 
Depmap.Protein.info3<-merge(x= Depmap.Protein.info2, y= Protein_Info4, 
                            by.x="Uniprot_Acc", by.y="UniProt ID(supplied by UniProt)", 
                            sort = TRUE)# 10594 Proteins
# Find those genes without uniprot id: 
No_UniprotID <- anti_join(Depmap.Protein.info2, Protein_Info4, #finding genes_Symbol without match
                          by = c("Uniprot_Acc" = "UniProt ID(supplied by UniProt)"))
#...find the genes with matching gene symbols
Depmap.Protein.info4<-merge(x= No_UniprotID, y= Protein_Info4, 
                            by.x="Gene_Symbol", by.y="Approved Symbol", 
                            sort = TRUE)# 10 collumns, 253 genes, only no gene-symbol genes with uniprot_ID

# Merge genes with gene symbol and genes with only uniprot ID. 
Depmap.Protein.info5<-merge(x= Depmap.Protein.info3, y= Depmap.Protein.info4, 
                            all=TRUE)# 11 collumns, 12100 genes


t.test_prot_category_Lung<- data.frame(Protein_ID=character(),
                                            Protein_Name=character(),
                                            Pvlaue.Tri.Di= numeric(),
                                            Diff.Gain= numeric(), 
                                            Pvlaue.Di.Mono= numeric(),
                                            Diff.Di_Mono= numeric(),
                                            Pvlaue.Tri.Mono= numeric(),
                                            Diff.GainvLoss= numeric())


for (i in 1:(length(Protein.Expression_Lung)-1)){
  TestProt = colnames(Protein.Expression_Lung[i])
  
  testProtChrm <- filter(Depmap.Protein.info5, Protein_ID==TestProt) #get test protein data
  
  
  #use If statement to check that Protein Info has Protein of interest. 
  if (length(testProtChrm[,1])!=0){
    TestArm <- testProtChrm$arm #Find test protein chromosome arm
    TestChrm<- as.numeric(testProtChrm$Chromosome) #Find test protein chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test protein location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    Protein_Expression_TestProt <- Protein.Expression_Lung %>% select(Cell_line, all_of(TestProt))
    
    #Combine protein data with aneuploidy data for test protein location
    Protein.effect.Chrm<-merge(y= Protein_Expression_TestProt, x= testchrm.percell, 
                               by.y="Cell_line", by.x="DepMap_ID", 
                               sort = TRUE)# 368 cells, 12757 proteins 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(Protein.effect.Chrm, arm_call==1) #all cells triploid for testProt chrm
    Chrm.di<- filter(Protein.effect.Chrm, arm_call==0) #all cells diploid for testProt chrm
    Chrm.mono<- filter(Protein.effect.Chrm, arm_call==-1) #all cells mono for testProt chrm
    
    ##Make data frame with t-test info about protein expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 & #[8] is column with protein expression data. check it has values. 
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 & #
        colSums(!is.na(Chrm.mono[8]))>=10) { #
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with Protein data
                     Chrm.di[,8], # [,8] because that is column with Protein data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Di.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                      Chrm.di[,8], # [,8] because that is column with Protein data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Di.Mono<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE)  #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                       Chrm.tri[,8], # [,8] because that is column with Protein data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_prot_category_Lung<- rbind(t.test_prot_category_Lung, 
                                             data.frame(Protein_ID=TestProt,
                                                        Protein_Name=testProtChrm$Gene_Symbol,
                                                        Pvlaue.Tri.Di= di.tri$p.value,
                                                        Diff.Gain= Diff.Tri.Di, 
                                                        Pvlaue.Di.Mono= Di.Mono$p.value,
                                                        Diff.Di_Mono= Diff.Di.Mono, 
                                                        Pvlaue.Tri.Mono= Tri.Mono$p.value,
                                                        Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_prot_category_Lung<-t.test_prot_category_Lung
    }
  }
}

t.test_prot_category_Lung<-distinct(t.test_prot_category_Lung) 

t.test_prot_category_Lung<-t.test_prot_category_Lung[order(t.test_prot_category_Lung$Pvlaue.Tri.Mono),]# *** genes
# length (number of genes) = 1684



### sub-step 3: get RNA difference in lung cells
# Lung cancer analysis: RNA expression difference: 
# prep data frame: 
colnames(RNA.Expression.filtered)<-sub("[..].*", "", as.character(colnames(RNA.Expression.filtered)))#12755 genes
RNA.Expression.Lung<-RNA.Expression.filtered[,-c(1)]

RNA.Expression.Lung<- subset(RNA.Expression.Lung, Cell_line %in% lung.cells$DepMap_ID)
#19145 genes, 64 cells

# Run loop to get difference in experssion per gene in only low aneuploid cells
t.test_RNA_category_Lung<- data.frame(RNA_ID=character(),
                                           RNA_Name=character(),
                                           Pvalue.Gain= numeric(),
                                           Diff.Gain= numeric(), 
                                           Pvalue.Loss= numeric(),
                                           Diff.Loss= numeric(),
                                           Pvalue.GainvLoss= numeric(),
                                           Diff.GainvLoss= numeric())

for (i in 2:length(RNA.Expression.Lung)){
  
  TestRNA_name=sub("[..].*", "", as.character(colnames(RNA.Expression.Lung[i])))
  TestRNA=colnames(RNA.Expression.Lung[i])
  
  testRNAChrm <- filter(RNA_Info3, RNA_Info3$"Approved Symbol"==TestRNA) #get test RNA data
  if (length(testRNAChrm$Chromosome)!=0){ #Make sure RNA data is in RNA_info3, else it crashes
    TestArm <- testRNAChrm$arm #Find test RNA chromosome arm
    TestChrm<- as.numeric(testRNAChrm$Chromosome) #Find test RNA chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test RNA location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    RNA_Expression_TestRNA <- RNA.Expression.Lung %>% select(Cell_line, all_of(TestRNA))
    
    #Combine RNA data with aneuploidy data for test RNA location
    RNA.effect.Chrm<-merge(y= RNA_Expression_TestRNA, x= testchrm.percell, 
                           by.y="Cell_line", by.x="DepMap_ID", 
                           sort = TRUE)# 368 cells, 12757 RNAs 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(RNA.effect.Chrm, arm_call==1) #all cells triploid for testRNA chrm
    Chrm.di<- filter(RNA.effect.Chrm, arm_call==0) #all cells diploid for testRNA chrm
    Chrm.mono<- filter(RNA.effect.Chrm, arm_call==-1) #all cells mono for testRNA chrm
    
    ##Make data frame with t-test info about RNA expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 &  
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 & 
        colSums(!is.na(Chrm.mono[8]))>=10) { 
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with RNA data
                     Chrm.di[,8], # [,8] because that is column with RNA data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Mono.Di<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                      Chrm.di[,8], # [,8] because that is column with RNA data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Mono.Di<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                       Chrm.tri[,8], # [,8] because that is column with RNA data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_RNA_category_Lung<- rbind(t.test_RNA_category_Lung, data.frame(
        RNA_ID=colnames(RNA.effect.Chrm[8]),
        RNA_Name=TestRNA_name,
        Pvalue.Gain= di.tri$p.value,
        Diff.Gain= Diff.Tri.Di, 
        Pvalue.Loss= Mono.Di$p.value,
        Diff.Loss= Diff.Mono.Di, 
        Pvalue.GainvLoss= Tri.Mono$p.value,
        Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_RNA_category_Lung<-t.test_RNA_category_Lung
    }
  }
}

t.test_RNA_category_Lung<-distinct(t.test_RNA_category_Lung) 
t.test_RNA_category_Lung<-t.test_RNA_category_Lung[order(t.test_RNA_category_Lung$Pvalue.GainvLoss),]
# length (# of RNA genes) : 4815

### Sub step 4: combine low ploidy RNA and protein data
CN.Diff.RNA.Prot_Lung<-merge(x=t.test_RNA_category_Lung, 
                                   y=t.test_prot_category_Lung, 
                                   by.x="RNA_Name", by.y="Protein_Name") #7874 genes

colnames(CN.Diff.RNA.Prot_Lung)<-c("RNA_Name", "RNA_ID", 
                                         "RNA.Pvalue.Gain", "RNA.Diff.Gain", 
                                         "RNA.Pvalue.Loss", "RNA.Diff.Loss", 
                                         "RNA.Pvalue.GainvLoss", "RNA.Diff.GainvLoss", 
                                         "Protein_ID", "Protein.Pvalue.Gain", 
                                         "Protein.Diff.Gain", "Protein.Pvalue.Loss", 
                                         "Protein.Diff.Loss", "Protein.Pvalue.GainvLoss", 
                                         "Protein.Diff.GainvLoss")

## add categories based on cutoffs.  (-Inf,-0.1,0.25,Inf)
CN.Diff.RNA.Prot_Lung$Three.RNA.Gain<- cut(CN.Diff.RNA.Prot_Lung$RNA.Diff.Gain,
                                                 breaks=c(-Inf,-0.1,0.25,Inf),
                                                 include.lowest=TRUE,
                                                 labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_Lung$Three.RNA.Loss<- cut(CN.Diff.RNA.Prot_Lung$RNA.Diff.Loss,
                                                 breaks=c(-Inf,-0.25,0.1,Inf),
                                                 include.lowest=TRUE,
                                                 labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.RNA.Prot_Lung$Three.Protein.Gain<- cut(CN.Diff.RNA.Prot_Lung$Protein.Diff.Gain,
                                                     breaks=c(-Inf,-0.1,0.25,Inf),
                                                     include.lowest=TRUE,
                                                     labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_Lung$Three.Protein.Loss<- cut(CN.Diff.RNA.Prot_Lung$Protein.Diff.Loss,
                                                     breaks=c(-Inf,-0.25,0.1,Inf),
                                                     include.lowest=TRUE,
                                                     labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.RNA.Prot_Lung$Three.Protein.Loss<-factor(CN.Diff.RNA.Prot_Lung$Three.Protein.Loss, 
                                                       levels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_Lung$Three.RNA.Loss<-factor(CN.Diff.RNA.Prot_Lung$Three.RNA.Loss, 
                                                   levels=c("Anti-Scaling","Buffering","Scaling"))


###Bar graph of RNA and Protein quantiles Scaling/Buffered

#format data so I can plot all 4 groups: RNA/Prot gain/loss
dat.m_lung <- melt(CN.Diff.RNA.Prot_Lung, id.vars='RNA_ID', 
                   measure.vars=c('Three.RNA.Gain','Three.Protein.Gain', 
                                  'Three.RNA.Loss', 'Three.Protein.Loss'))

## Barplot of scaling/buffering/anti-scaling categories
pdf(height=4, width=4, "plot.Protein.RNA.Bargraph.ThreeCat.Gain_Lung.min10.pdf")
ggplot(data= dat.m_lung, aes(x=variable, fill=value)) + 
  geom_bar(position="fill")+ #dodge (next to eachother)
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()
#4x4
# plot.Protein.RNA.Bargraph.ThreeCat.Gain_Lung.min10.pdf

# percent of genes in each group
# RNA gain = 9.5, 37.9 and 52.6
table(subset(dat.m_lung,variable=="Three.RNA.Gain")$value)/length(subset(dat.m_lung,variable=="Three.RNA.Gain")$value)

# protein gain = 11.6, 52.5, 35.9
table(subset(dat.m_lung,variable=="Three.Protein.Gain")$value)/length(subset(dat.m_lung,variable=="Three.RNA.Gain")$value)
# RNA loss = 3.7, 26.8, 69.4
table(subset(dat.m_lung,variable=="Three.RNA.Loss")$value)/length(subset(dat.m_lung,variable=="Three.RNA.Gain")$value)
# Protein loss = 12.2, 56.2, 31.6
table(subset(dat.m_lung,variable=="Three.Protein.Loss")$value)/length(subset(dat.m_lung,variable=="Three.RNA.Gain")$value)


write.csv(CN.Diff.RNA.Prot_Lung, #change name as needed
          file =paste("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_Lung.min10cells.csv", sep=','), 
          row.names = TRUE)

write.csv(lung.cells$DepMap_ID, #save a list of the cell lines used in analysis. 
          file =paste("CellLines_LungCancer", sep=','), 
          row.names = TRUE)

#CN.Diff.RNA.Prot_Lung<- read.csv("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_Lung.min10cells.csv")


###### Control dataset: . No mutations #####
# the dataset this section generates is available in supplementary data 3, sheet 3 "No_Mutations" 


# get protein data only for proteins with minimum 10 cells with protein data
Protein.Expression.filtered_min10Cells<- Protein.Expression.filtered %>% select(one_of(CN.Diff.xRNA.yProt.ThreeGroups$Protein_ID))
# Now add Cell_lines back in. 
Protein.Expression.filtered_min10Cells$Cell_line<- Protein.Expression.filtered$Cell_line

#colnames(RNA.Expression.filtered)<-sub("[..].*", "", as.character(colnames(RNA.Expression.filtered)))#12755 genes
#RNA.Expression.filtered<-RNA.Expression.filtered[,-c(1)]

### remove mutated genes 
# Get list of genes and mutations per cell line
# filter mutation list to get cell lines in our dataset
# add protein ID to each gene mutation info
# for-loop through mutation list and for each cell line, remove all mutated genes

# Mutational score
# get all mutations in all cell lines
# subset only cells in our analysis
# cell line list from Protein_RNA_filtered_CellLine.R
# count table for gene mutations per gene
setwd(DataFileLocation)
mutation<-read_csv("CCLE_mutations.csv")
mutations1<-mutation[,c(1,2,16)] # 
mutations1<- subset(mutations1, DepMap_ID %in% Protein.Expression.filtered_min10Cells$Cell_line) #length = 50607 mutations

mutations.info<-merge(mutations1, Protein_Info4, by.x="Entrez_Gene_Id", by.y="Entrez Gene ID")#50260
mutations.info2<-anti_join(mutations1, Protein_Info4, by=c("Entrez_Gene_Id" = "Entrez Gene ID") ) #347
mutations.info3<-merge(mutations.info2, Protein_Info4, by.x="Hugo_Symbol", by.y="Approved Symbol") #320
#mutations.info4<-anti_join(mutations.info2, Protein_Info4, by=c("Hugo_Symbol" = "Approved Symbol") )
# 27 genes not found back in the protein info dataset. no known proteins, perhaps. Most at orfs
# example 9 of these genes are CXorf22, and three are C10orf12. 
mutations.info5<- merge(x= mutations.info, y= mutations.info3, 
                        all=TRUE) #50580

mutations.info6<- merge(x= mutations.info5, y= Protein_ProID, 
                        by.x="Hugo_Symbol", 
                        by.y= "Gene_Symbol") #35926
mutations.info7<-anti_join(mutations.info5, Protein_ProID, by=c("Hugo_Symbol" = "Gene_Symbol") ) #
mutations.info8<-merge(mutations.info7, Protein_ProID, 
                       by.x="UniProt ID(supplied by UniProt)", by.y="Uniprot_Acc") #93
mutations.info9<- merge(x= mutations.info6, y= mutations.info8, 
                        all=TRUE) #36019 mutations (many are same genes mutated in diff cells)

# now go through mutations list and delete corresponding mutated protein in that cell line
Protein.Expression_noMut<-Protein.Expression.filtered_min10Cells
RNA.Expression_noMut<- RNA.Expression.filtered
RNA_Col<-colnames(RNA.Expression_noMut)
Prot_Col<-colnames(Protein.Expression_noMut)

for (i in 1:length(mutations.info9$Entrez_Gene_Id)){
  MutCell=mutations.info9$DepMap_ID[i] 
  MutProt=mutations.info9$Protein_Id[i]
  MutGene=mutations.info9$Hugo_Symbol[i]
  #Protein: get row name for row with cell line for protein expression
  MutRow<-rownames(Protein.Expression_noMut[Protein.Expression_noMut$Cell_Line==MutCell,])
  # Replace protein expression data that had mutated protein name in mutated cell line with NA:  
  if (MutProt %in% Prot_Col){
    Protein.Expression_noMut[MutRow, MutProt]<-NA 
  }
  if (MutGene %in% RNA_Col){
    # Replace RNA expression data that had mutated RNA name in mutated cell line with NA:  
    RNA.Expression_noMut[MutRow, MutGene]<-NA 
  }
}
# proteins: 9414 genes
# RNA: 19146 genes
Protein.Expression_noMut
RNA.Expression_noMut



## Step 2: get difference data
## Now that I have dataset with no mutations
## find difference in expression upon chrm gain and loss
## then find mean RNA/Protein difference in expression upon chrm gain/loss 
## also find number of genes

### sub-step 2: get protein difference in no mutation genes 
Protein.ID<-data.frame(Protein_ID=
                         colnames(Protein.Expression_noMut[3:length(Protein.Expression_noMut)]))



t.test_prot_category_noMut<- data.frame(Protein_ID=character(),
                                       Protein_Name=character(),
                                       Pvlaue.Tri.Di= numeric(),
                                       Diff.Gain= numeric(), 
                                       Pvlaue.Di.Mono= numeric(),
                                       Diff.Di_Mono= numeric(),
                                       Pvlaue.Tri.Mono= numeric(),
                                       Diff.GainvLoss= numeric())


for (i in 1:(length(Protein.Expression_noMut)-1)){
  TestProt = colnames(Protein.Expression_noMut[i])
  
  testProtChrm <- filter(Depmap.Protein.info5, Protein_ID==TestProt) #get test protein data
  
  
  #use If statement to check that Protein Info has Protein of interest. 
  if (length(testProtChrm[,1])!=0){
    TestArm <- testProtChrm$arm #Find test protein chromosome arm
    TestChrm<- as.numeric(testProtChrm$Chromosome) #Find test protein chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test protein location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    Protein_Expression_TestProt <- Protein.Expression_noMut %>% select(Cell_line, all_of(TestProt))
    
    #Combine protein data with aneuploidy data for test protein location
    Protein.effect.Chrm<-merge(y= Protein_Expression_TestProt, x= testchrm.percell, 
                               by.y="Cell_line", by.x="DepMap_ID", 
                               sort = TRUE)# 368 cells, 12757 proteins 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(Protein.effect.Chrm, arm_call==1) #all cells triploid for testProt chrm
    Chrm.di<- filter(Protein.effect.Chrm, arm_call==0) #all cells diploid for testProt chrm
    Chrm.mono<- filter(Protein.effect.Chrm, arm_call==-1) #all cells mono for testProt chrm
    
    ##Make data frame with t-test info about protein expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 & #[8] is column with protein expression data. check it has values. 
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 & #usually 10, now 3**
        colSums(!is.na(Chrm.mono[8]))>=10) { #usually 10, now 3**
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with Protein data
                     Chrm.di[,8], # [,8] because that is column with Protein data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Di.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                      Chrm.di[,8], # [,8] because that is column with Protein data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Di.Mono<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE)  #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                       Chrm.tri[,8], # [,8] because that is column with Protein data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_prot_category_noMut<- rbind(t.test_prot_category_noMut, 
                                        data.frame(Protein_ID=TestProt,
                                                   Protein_Name=testProtChrm$Gene_Symbol,
                                                   Pvlaue.Tri.Di= di.tri$p.value,
                                                   Diff.Gain= Diff.Tri.Di, 
                                                   Pvlaue.Di.Mono= Di.Mono$p.value,
                                                   Diff.Di_Mono= Diff.Di.Mono, 
                                                   Pvlaue.Tri.Mono= Tri.Mono$p.value,
                                                   Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_prot_category_noMut<-t.test_prot_category_noMut
    }
  }
}

t.test_prot_category_noMut<-distinct(t.test_prot_category_noMut) 

t.test_prot_category_noMut<-t.test_prot_category_noMut[order(t.test_prot_category_noMut$Pvlaue.Tri.Mono),]# *** genes
# length (number of genes) = 9443  



### sub-step 3: get RNA difference in lung cells
# Lung cancer analysis: RNA expression difference: 
# prep data frame: 

# Run loop to get difference in experssion per gene in only low aneuploid cells
t.test_RNA_category_noMut<- data.frame(RNA_ID=character(),
                                      RNA_Name=character(),
                                      Pvalue.Gain= numeric(),
                                      Diff.Gain= numeric(), 
                                      Pvalue.Loss= numeric(),
                                      Diff.Loss= numeric(),
                                      Pvalue.GainvLoss= numeric(),
                                      Diff.GainvLoss= numeric())

for (i in 2:length(RNA.Expression_noMut)){
  
  TestRNA_name=sub("[..].*", "", as.character(colnames(RNA.Expression_noMut[i])))
  TestRNA=colnames(RNA.Expression_noMut[i])
  
  testRNAChrm <- filter(RNA_Info3, RNA_Info3$"Approved Symbol"==TestRNA) #get test RNA data
  if (length(testRNAChrm$Chromosome)!=0){ #Make sure RNA data is in RNA_info3, else it crashes
    TestArm <- testRNAChrm$arm #Find test RNA chromosome arm
    TestChrm<- as.numeric(testRNAChrm$Chromosome) #Find test RNA chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test RNA location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    RNA_Expression_TestRNA <- RNA.Expression_noMut %>% select(Cell_line, all_of(TestRNA))
    
    #Combine RNA data with aneuploidy data for test RNA location
    RNA.effect.Chrm<-merge(y= RNA_Expression_TestRNA, x= testchrm.percell, 
                           by.y="Cell_line", by.x="DepMap_ID", 
                           sort = TRUE)# 368 cells, 12757 RNAs 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(RNA.effect.Chrm, arm_call==1) #all cells triploid for testRNA chrm
    Chrm.di<- filter(RNA.effect.Chrm, arm_call==0) #all cells diploid for testRNA chrm
    Chrm.mono<- filter(RNA.effect.Chrm, arm_call==-1) #all cells mono for testRNA chrm
    
    ##Make data frame with t-test info about RNA expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 &  
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 & 
        colSums(!is.na(Chrm.mono[8]))>=10) { 
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with RNA data
                     Chrm.di[,8], # [,8] because that is column with RNA data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Mono.Di<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                      Chrm.di[,8], # [,8] because that is column with RNA data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Mono.Di<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                       Chrm.tri[,8], # [,8] because that is column with RNA data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_RNA_category_noMut<- rbind(t.test_RNA_category_noMut, data.frame(
        RNA_ID=colnames(RNA.effect.Chrm[8]),
        RNA_Name=TestRNA_name,
        Pvalue.Gain= di.tri$p.value,
        Diff.Gain= Diff.Tri.Di, 
        Pvalue.Loss= Mono.Di$p.value,
        Diff.Loss= Diff.Mono.Di, 
        Pvalue.GainvLoss= Tri.Mono$p.value,
        Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_RNA_category_noMut<-t.test_RNA_category_noMut
    }
  }
}

t.test_RNA_category_noMut<-distinct(t.test_RNA_category_noMut) 
t.test_RNA_category_noMut<-t.test_RNA_category_noMut[order(t.test_RNA_category_noMut$Pvalue.GainvLoss),]
# length (# of RNA genes) : 

### Sub step 4: combine low ploidy RNA and protein data
CN.Diff.RNA.Prot_noMut<-merge(x=t.test_RNA_category_noMut, 
                             y=t.test_prot_category_noMut, 
                             by.x="RNA_Name", by.y="Protein_Name") #7874 genes

colnames(CN.Diff.RNA.Prot_noMut)<-c("RNA_Name", "RNA_ID", 
                                   "RNA.Pvalue.Gain", "RNA.Diff.Gain", 
                                   "RNA.Pvalue.Loss", "RNA.Diff.Loss", 
                                   "RNA.Pvalue.GainvLoss", "RNA.Diff.GainvLoss", 
                                   "Protein_ID", "Protein.Pvalue.Gain", 
                                   "Protein.Diff.Gain", "Protein.Pvalue.Loss", 
                                   "Protein.Diff.Loss", "Protein.Pvalue.GainvLoss", 
                                   "Protein.Diff.GainvLoss")

## add categories based on cutoffs.  (-Inf,-0.1,0.25,Inf)
CN.Diff.RNA.Prot_noMut$Three.RNA.Gain<- cut(CN.Diff.RNA.Prot_noMut$RNA.Diff.Gain,
                                           breaks=c(-Inf,-0.1,0.25,Inf),
                                           include.lowest=TRUE,
                                           labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_noMut$Three.RNA.Loss<- cut(CN.Diff.RNA.Prot_noMut$RNA.Diff.Loss,
                                           breaks=c(-Inf,-0.25,0.1,Inf),
                                           include.lowest=TRUE,
                                           labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.RNA.Prot_noMut$Three.Protein.Gain<- cut(CN.Diff.RNA.Prot_noMut$Protein.Diff.Gain,
                                               breaks=c(-Inf,-0.1,0.25,Inf),
                                               include.lowest=TRUE,
                                               labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_noMut$Three.Protein.Loss<- cut(CN.Diff.RNA.Prot_noMut$Protein.Diff.Loss,
                                               breaks=c(-Inf,-0.25,0.1,Inf),
                                               include.lowest=TRUE,
                                               labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.RNA.Prot_noMut$Three.Protein.Loss<-factor(CN.Diff.RNA.Prot_noMut$Three.Protein.Loss, 
                                                 levels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_noMut$Three.RNA.Loss<-factor(CN.Diff.RNA.Prot_noMut$Three.RNA.Loss, 
                                             levels=c("Anti-Scaling","Buffering","Scaling"))


#setwd(DataFileLocation)
write.csv(CN.Diff.RNA.Prot_noMut, #change name as needed
          file =paste("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_noMut_Min10Cells.csv", sep=','), 
          row.names = TRUE)
#CN.Diff.RNA.Prot_noMut<- read.csv("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_noMut_Min10Cells.csv")


###### Control dataset: . Remove low RNA expression genes #####
# the dataset this section generates is available in supplementary data 3, sheet 5 "No_Low_Expression" 

## remove genes that have lowest 20% of RNA expression
## analyze genes with top 80% of expression

#make dataframe with mean RNA expression per gene
Mean_RNA_perGene<-data.frame(RNA_Name=as.character(), 
                             MeanRNAExpress=as.numeric())

for (i in 1:(length(RNA.Expression.filtered)-1)){
  Mean_RNA_perGene<- rbind(Mean_RNA_perGene, 
                           data.frame(RNA_Name=colnames(RNA.Expression.filtered[i]), 
                                      MeanRNAExpress=mean(RNA.Expression.filtered[,i], na.rm=TRUE)))
}

Mean_RNA_perGene<-Mean_RNA_perGene[-c(1),]# remove cell_lines
Mean_RNA_perGene<- Mean_RNA_perGene %>%
  mutate(quantile = ntile(MeanRNAExpress, 10))
Mean_RNA_perGene$RNA_Name[1:10]
table(Mean_RNA_perGene$quantile) #19143

# add protein_ID per RNA gene to get protein data/gene (not all RNA have Protein_ID attached)
# Protein_ProID has 12755 proteins
Mean_RNA_perGene2<- merge(x= Mean_RNA_perGene, y= Protein_ProID, 
                          by.x="RNA_Name", by.y= "Gene_Symbol") # 12105


# then split into high-low expression categories
LowRNAExpress<- subset(Mean_RNA_perGene2, quantile <= 2) #get all genes with lowest 20% RNA expression levels
HighRNAExpress<- subset(Mean_RNA_perGene2, quantile > 2) #get all genes with higher 80% RNA expression levels


# Filter RNA & protein data to get higher expression genes only
Protein.Expression_highRNA<- Protein.Expression.filtered_min10Cells[, colnames(Protein.Expression.filtered_min10Cells) %in% HighRNAExpress$Protein_Id]
Protein.Expression_highRNA$Cell_Line<- Protein.Expression.filtered_min10Cells$Cell_line #also add back cell line info
# number of cells= 371,  genes 9383

RNA.Expression_highRNA<- RNA.Expression.filtered[, colnames(RNA.Expression.filtered) %in% HighRNAExpress$RNA_Name]
RNA.Expression_highRNA$Cell_Line<- RNA.Expression.filtered$Cell_line
# number of cells= 371,  genes 11510


## Step 2: get difference data
## Now that I have dataset with no mutations
## find difference in expression upon chrm gain and loss
## then find mean RNA/Protein difference in expression upon chrm gain/loss 
## also find number of genes

### sub-step 2: get protein difference in no mutation genes 
Protein.ID<-data.frame(Protein_ID=
                         colnames(Protein.Expression_highRNA[3:length(Protein.Expression_highRNA)]))



t.test_prot_category_highRNA<- data.frame(Protein_ID=character(),
                                        Protein_Name=character(),
                                        Pvlaue.Tri.Di= numeric(),
                                        Diff.Gain= numeric(), 
                                        Pvlaue.Di.Mono= numeric(),
                                        Diff.Di_Mono= numeric(),
                                        Pvlaue.Tri.Mono= numeric(),
                                        Diff.GainvLoss= numeric())


for (i in 1:(length(Protein.Expression_highRNA)-1)){
  TestProt = colnames(Protein.Expression_highRNA[i])
  
  testProtChrm <- filter(Depmap.Protein.info5, Protein_ID==TestProt) #get test protein data
  
  
  #use If statement to check that Protein Info has Protein of interest. 
  if (length(testProtChrm[,1])!=0){
    TestArm <- testProtChrm$arm #Find test protein chromosome arm
    TestChrm<- as.numeric(testProtChrm$Chromosome) #Find test protein chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test protein location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    Protein_Expression_TestProt <- Protein.Expression_highRNA %>% select(Cell_Line, all_of(TestProt))
    
    #Combine protein data with aneuploidy data for test protein location
    Protein.effect.Chrm<-merge(y= Protein_Expression_TestProt, x= testchrm.percell, 
                               by.y="Cell_Line", by.x="DepMap_ID", 
                               sort = TRUE)# 368 cells, 12757 proteins 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(Protein.effect.Chrm, arm_call==1) #all cells triploid for testProt chrm
    Chrm.di<- filter(Protein.effect.Chrm, arm_call==0) #all cells diploid for testProt chrm
    Chrm.mono<- filter(Protein.effect.Chrm, arm_call==-1) #all cells mono for testProt chrm
    
    ##Make data frame with t-test info about protein expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 & #[8] is column with protein expression data. check it has values. 
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 & #usually 10, now 3**
        colSums(!is.na(Chrm.mono[8]))>=10) { #usually 10, now 3**
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with Protein data
                     Chrm.di[,8], # [,8] because that is column with Protein data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Di.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                      Chrm.di[,8], # [,8] because that is column with Protein data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Di.Mono<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE)  #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                       Chrm.tri[,8], # [,8] because that is column with Protein data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_prot_category_highRNA<- rbind(t.test_prot_category_highRNA, 
                                         data.frame(Protein_ID=TestProt,
                                                    Protein_Name=testProtChrm$Gene_Symbol,
                                                    Pvlaue.Tri.Di= di.tri$p.value,
                                                    Diff.Gain= Diff.Tri.Di, 
                                                    Pvlaue.Di.Mono= Di.Mono$p.value,
                                                    Diff.Di_Mono= Diff.Di.Mono, 
                                                    Pvlaue.Tri.Mono= Tri.Mono$p.value,
                                                    Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_prot_category_highRNA<-t.test_prot_category_highRNA
    }
  }
}

t.test_prot_category_highRNA<-distinct(t.test_prot_category_highRNA) 

t.test_prot_category_highRNA<-t.test_prot_category_highRNA[order(t.test_prot_category_highRNA$Pvlaue.Tri.Mono),]# *** genes
# length (number of genes) = 9257



### sub-step 3: get RNA difference in lung cells
# Lung cancer analysis: RNA expression difference: 
# prep data frame: 


# Run loop to get difference in experssion per gene in only low aneuploid cells
t.test_RNA_category_highRNA<- data.frame(RNA_ID=character(),
                                       RNA_Name=character(),
                                       Pvalue.Gain= numeric(),
                                       Diff.Gain= numeric(), 
                                       Pvalue.Loss= numeric(),
                                       Diff.Loss= numeric(),
                                       Pvalue.GainvLoss= numeric(),
                                       Diff.GainvLoss= numeric())

for (i in 2:length(RNA.Expression_highRNA)){
  
  TestRNA_name=sub("[..].*", "", as.character(colnames(RNA.Expression_highRNA[i])))
  TestRNA=colnames(RNA.Expression_highRNA[i])
  
  testRNAChrm <- filter(RNA_Info3, RNA_Info3$"Approved Symbol"==TestRNA) #get test RNA data
  if (length(testRNAChrm$Chromosome)!=0){ #Make sure RNA data is in RNA_info3, else it crashes
    TestArm <- testRNAChrm$arm #Find test RNA chromosome arm
    TestChrm<- as.numeric(testRNAChrm$Chromosome) #Find test RNA chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test RNA location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    RNA_Expression_TestRNA <- RNA.Expression_highRNA %>% select(Cell_Line, all_of(TestRNA))
    
    #Combine RNA data with aneuploidy data for test RNA location
    RNA.effect.Chrm<-merge(y= RNA_Expression_TestRNA, x= testchrm.percell, 
                           by.y="Cell_Line", by.x="DepMap_ID", 
                           sort = TRUE)# 368 cells, 12757 RNAs 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(RNA.effect.Chrm, arm_call==1) #all cells triploid for testRNA chrm
    Chrm.di<- filter(RNA.effect.Chrm, arm_call==0) #all cells diploid for testRNA chrm
    Chrm.mono<- filter(RNA.effect.Chrm, arm_call==-1) #all cells mono for testRNA chrm
    
    ##Make data frame with t-test info about RNA expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 &  
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 & 
        colSums(!is.na(Chrm.mono[8]))>=10) { 
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with RNA data
                     Chrm.di[,8], # [,8] because that is column with RNA data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Mono.Di<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                      Chrm.di[,8], # [,8] because that is column with RNA data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Mono.Di<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                       Chrm.tri[,8], # [,8] because that is column with RNA data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_RNA_category_highRNA<- rbind(t.test_RNA_category_highRNA, data.frame(
        RNA_ID=colnames(RNA.effect.Chrm[8]),
        RNA_Name=TestRNA_name,
        Pvalue.Gain= di.tri$p.value,
        Diff.Gain= Diff.Tri.Di, 
        Pvalue.Loss= Mono.Di$p.value,
        Diff.Loss= Diff.Mono.Di, 
        Pvalue.GainvLoss= Tri.Mono$p.value,
        Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_RNA_category_highRNA<-t.test_RNA_category_highRNA
    }
  }
}
t.test_RNA_category_highRNA<-distinct(t.test_RNA_category_highRNA) 
t.test_RNA_category_highRNA<-t.test_RNA_category_highRNA[order(t.test_RNA_category_highRNA$Pvalue.GainvLoss),]
# length (# of RNA genes) : 10659

### Sub step 4: combine low ploidy RNA and protein data
CN.Diff.RNA.Prot_highRNA<-merge(x=t.test_RNA_category_highRNA, 
                              y=t.test_prot_category_highRNA, 
                              by.x="RNA_Name", by.y="Protein_Name") #7874 genes

colnames(CN.Diff.RNA.Prot_highRNA)<-c("RNA_Name", "RNA_ID", 
                                    "RNA.Pvalue.Gain", "RNA.Diff.Gain", 
                                    "RNA.Pvalue.Loss", "RNA.Diff.Loss", 
                                    "RNA.Pvalue.GainvLoss", "RNA.Diff.GainvLoss", 
                                    "Protein_ID", "Protein.Pvalue.Gain", 
                                    "Protein.Diff.Gain", "Protein.Pvalue.Loss", 
                                    "Protein.Diff.Loss", "Protein.Pvalue.GainvLoss", 
                                    "Protein.Diff.GainvLoss")

## add categories based on cutoffs.  (-Inf,-0.1,0.25,Inf)
CN.Diff.RNA.Prot_highRNA$Three.RNA.Gain<- cut(CN.Diff.RNA.Prot_highRNA$RNA.Diff.Gain,
                                            breaks=c(-Inf,-0.1,0.25,Inf),
                                            include.lowest=TRUE,
                                            labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_highRNA$Three.RNA.Loss<- cut(CN.Diff.RNA.Prot_highRNA$RNA.Diff.Loss,
                                            breaks=c(-Inf,-0.25,0.1,Inf),
                                            include.lowest=TRUE,
                                            labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.RNA.Prot_highRNA$Three.Protein.Gain<- cut(CN.Diff.RNA.Prot_highRNA$Protein.Diff.Gain,
                                                breaks=c(-Inf,-0.1,0.25,Inf),
                                                include.lowest=TRUE,
                                                labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_highRNA$Three.Protein.Loss<- cut(CN.Diff.RNA.Prot_highRNA$Protein.Diff.Loss,
                                                breaks=c(-Inf,-0.25,0.1,Inf),
                                                include.lowest=TRUE,
                                                labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.RNA.Prot_highRNA$Three.Protein.Loss<-factor(CN.Diff.RNA.Prot_highRNA$Three.Protein.Loss, 
                                                  levels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_highRNA$Three.RNA.Loss<-factor(CN.Diff.RNA.Prot_highRNA$Three.RNA.Loss, 
                                              levels=c("Anti-Scaling","Buffering","Scaling"))


#setwd(DataFileLocation)

#write.csv(CN.Diff.RNA.Prot_highRNA, #change name as needed
#          file =paste("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_highRNA_Min10Cells.csv", sep=','), 
#          row.names = TRUE) #9289 genes

CN.Diff.RNA.Prot_highRNA<-read.csv("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_highRNA_Min10Cells.csv")


###### Control dataset: . Remove flipped Gene CN and Chrm CN genes/cell #####
# the dataset this section generates is available in supplementary data 3, sheet 4 "No_Flipped" 

## Sub-step 1: filter protein/RNA expression data and get rid of "flipped" genes
## Flipped: genes whose gene copy number increases if located on lost chromosome
## or whose copy number decreases if located on gained chromosome

## Remove genes whose gene copy number was gained upon cells with coresponding chrmosome loss
## and remove genes whose gene copy number was lost uponc chrmosome gain
## "Gene CN flipped" dataset generated in TSG.OG_CNV_Difference_v2.R 
List.GeneCN.flipped.ChrmCN<- read.csv("List.Flipped.GeneCN.ChrmCN.csv")

### remove flipped genes 
# Get list of "flipped" genes
# add protein ID to each flipped gene info
# for-loop through flipped list for each cell line, remove all flipped genes

Flipped.info<-merge(List.GeneCN.flipped.ChrmCN, Protein_Info4, by.x="Gene_ID", by.y="Approved Symbol")#

# now go through flipped list and delete corresponding flipped protein in that cell line
Protein.Expression_noFlip<-Protein.Expression.filtered_min10Cells
RNA.Expression_noFlip<- RNA.Expression.filtered
RNA_Col<-colnames(RNA.Expression_noFlip)

#for proteins: get gene name from protein ID
#get protein expression with gene name labels
Prot_Col<- sub(".*[.]", "", as.character(colnames(Protein.Expression_noFlip)))
Prot_Col<- sub("_.*", "", Prot_Col)
Protein.Expression_noFlip2<-Protein.Expression_noFlip
colnames(Protein.Expression_noFlip2)<-Prot_Col

#now replace mutated genes/cell with NA in both protein and RNA datasets
for (i in 1:length(Flipped.info$Gene_ID)){
  FlipCell=as.character(Flipped.info$Cell_line[i])
  FlipGene=as.character(Flipped.info$Gene_ID[i])
  #Protein: get row name for row with cell line for protein expression
  FlipRow<-rownames(RNA.Expression_noFlip[RNA.Expression_noFlip$Cell_line==FlipCell,])
  # Replace protein expression data that had mutated protein name in mutated cell line with NA:  
  if (FlipGene %in% Prot_Col){
    Protein.Expression_noFlip2[FlipRow, FlipGene]<-NA 
  }
  if (FlipGene %in% RNA_Col){
    # Replace RNA expression data that had mutated RNA name in mutated cell line with NA:  
    RNA.Expression_noFlip[FlipRow, FlipGene]<-NA 
  }
}
# proteins: 9383 genes, 367 cell lines 
# RNA: 11510 genes, 367 cell lines 
# This takes a few hours to run, FYI

#now add protein_ID's back to the protein dataframe
colnames(Protein.Expression_noFlip2)<-colnames(Protein.Expression_noFlip)
RNA.Expression_noFlip

# cell lines = 371
# number of genes = 9414

setwd(DataFileLocation) 
write.csv(RNA.Expression_noFlip, file="RNA.expression.data_noFlipped.GeneCN.ChrmCN.csv")
write.csv(Protein.Expression_noFlip2, file="Protein.expression.data_noFlipped.GeneCN.ChrmCN.csv")


#### Step 2: get difference data
## Now that I have dataset with no mutations
## find difference in expression upon chrm gain and loss
## then find mean RNA/Protein difference in expression upon chrm gain/loss 
## also find number of genes

### sub-step 2: get protein difference in no mutation genes 
Protein.ID<-data.frame(Protein_ID=
                         colnames(Protein.Expression_noFlip2[3:length(Protein.Expression_noFlip2)]))


t.test_prot_category_noFlip<- data.frame(Protein_ID=character(),
                                          Protein_Name=character(),
                                          Pvlaue.Tri.Di= numeric(),
                                          Diff.Gain= numeric(), 
                                          Pvlaue.Di.Mono= numeric(),
                                          Diff.Di_Mono= numeric(),
                                          Pvlaue.Tri.Mono= numeric(),
                                          Diff.GainvLoss= numeric())


for (i in 1:(length(Protein.Expression_noFlip2)-1)){
  TestProt = colnames(Protein.Expression_noFlip2[i])
  
  testProtChrm <- filter(Depmap.Protein.info5, Protein_ID==TestProt) #get test protein data
  
  
  #use If statement to check that Protein Info has Protein of interest. 
  if (length(testProtChrm[,1])!=0){
    TestArm <- testProtChrm$arm #Find test protein chromosome arm
    TestChrm<- as.numeric(testProtChrm$Chromosome) #Find test protein chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test protein location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    Protein_Expression_TestProt <- Protein.Expression_noFlip2 %>% select(Cell_line, all_of(TestProt))
    
    #Combine protein data with aneuploidy data for test protein location
    Protein.effect.Chrm<-merge(y= Protein_Expression_TestProt, x= testchrm.percell, 
                               by.y="Cell_line", by.x="DepMap_ID", 
                               sort = TRUE)# 368 cells, 12757 proteins 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(Protein.effect.Chrm, arm_call==1) #all cells triploid for testProt chrm
    Chrm.di<- filter(Protein.effect.Chrm, arm_call==0) #all cells diploid for testProt chrm
    Chrm.mono<- filter(Protein.effect.Chrm, arm_call==-1) #all cells mono for testProt chrm
    
    ##Make data frame with t-test info about protein expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 & #[8] is column with protein expression data. check it has values. 
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 & #usually 10, now 3**
        colSums(!is.na(Chrm.mono[8]))>=10) { #usually 10, now 3**
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with Protein data
                     Chrm.di[,8], # [,8] because that is column with Protein data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Di.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                      Chrm.di[,8], # [,8] because that is column with Protein data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Di.Mono<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE)  #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                       Chrm.tri[,8], # [,8] because that is column with Protein data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_prot_category_noFlip<- rbind(t.test_prot_category_noFlip, 
                                           data.frame(Protein_ID=TestProt,
                                                      Protein_Name=testProtChrm$Gene_Symbol,
                                                      Pvlaue.Tri.Di= di.tri$p.value,
                                                      Diff.Gain= Diff.Tri.Di, 
                                                      Pvlaue.Di.Mono= Di.Mono$p.value,
                                                      Diff.Di_Mono= Diff.Di.Mono, 
                                                      Pvlaue.Tri.Mono= Tri.Mono$p.value,
                                                      Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_prot_category_noFlip<-t.test_prot_category_noFlip
    }
  }
}

t.test_prot_category_noFlip<-distinct(t.test_prot_category_noFlip) 

t.test_prot_category_noFlip<-t.test_prot_category_noFlip[order(t.test_prot_category_noFlip$Pvlaue.Tri.Mono),]# *** genes
# length (number of genes) = 9421



### sub-step 3: get RNA difference in lung cells
# Lung cancer analysis: RNA expression difference: 
# prep data frame: 

# Run loop to get difference in experssion per gene in only low aneuploid cells
t.test_RNA_category_noFlip<- data.frame(RNA_ID=character(),
                                         RNA_Name=character(),
                                         Pvalue.Gain= numeric(),
                                         Diff.Gain= numeric(), 
                                         Pvalue.Loss= numeric(),
                                         Diff.Loss= numeric(),
                                         Pvalue.GainvLoss= numeric(),
                                         Diff.GainvLoss= numeric())

for (i in 2:length(RNA.Expression_noFlip)){
  
  TestRNA_name=sub("[..].*", "", as.character(colnames(RNA.Expression_noFlip[i])))
  TestRNA=colnames(RNA.Expression_noFlip[i])
  
  testRNAChrm <- filter(RNA_Info3, RNA_Info3$"Approved Symbol"==TestRNA) #get test RNA data
  if (length(testRNAChrm$Chromosome)!=0){ #Make sure RNA data is in RNA_info3, else it crashes
    TestArm <- testRNAChrm$arm #Find test RNA chromosome arm
    TestChrm<- as.numeric(testRNAChrm$Chromosome) #Find test RNA chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test RNA location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    RNA_Expression_TestRNA <- RNA.Expression_noFlip %>% select(Cell_line, all_of(TestRNA))
    
    #Combine RNA data with aneuploidy data for test RNA location
    RNA.effect.Chrm<-merge(y= RNA_Expression_TestRNA, x= testchrm.percell, 
                           by.y="Cell_line", by.x="DepMap_ID", 
                           sort = TRUE)# 368 cells, 12757 RNAs 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(RNA.effect.Chrm, arm_call==1) #all cells triploid for testRNA chrm
    Chrm.di<- filter(RNA.effect.Chrm, arm_call==0) #all cells diploid for testRNA chrm
    Chrm.mono<- filter(RNA.effect.Chrm, arm_call==-1) #all cells mono for testRNA chrm
    
    ##Make data frame with t-test info about RNA expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 &  
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 & 
        colSums(!is.na(Chrm.mono[8]))>=10) { 
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with RNA data
                     Chrm.di[,8], # [,8] because that is column with RNA data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Mono.Di<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                      Chrm.di[,8], # [,8] because that is column with RNA data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Mono.Di<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                       Chrm.tri[,8], # [,8] because that is column with RNA data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_RNA_category_noFlip<- rbind(t.test_RNA_category_noFlip, data.frame(
        RNA_ID=colnames(RNA.effect.Chrm[8]),
        RNA_Name=TestRNA_name,
        Pvalue.Gain= di.tri$p.value,
        Diff.Gain= Diff.Tri.Di, 
        Pvalue.Loss= Mono.Di$p.value,
        Diff.Loss= Diff.Mono.Di, 
        Pvalue.GainvLoss= Tri.Mono$p.value,
        Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_RNA_category_noFlip<-t.test_RNA_category_noFlip
    }
  }
}

t.test_RNA_category_noFlip<-distinct(t.test_RNA_category_noFlip) 
t.test_RNA_category_noFlip<-t.test_RNA_category_noFlip[order(t.test_RNA_category_noFlip$Pvalue.GainvLoss),]
# length (# of RNA genes) : 17774


### Sub step 4: combine low ploidy RNA and protein data
CN.Diff.RNA.Prot_noFlip<-merge(x=t.test_RNA_category_noFlip, 
                                y=t.test_prot_category_noFlip, 
                                by.x="RNA_Name", by.y="Protein_Name") #7874 genes

colnames(CN.Diff.RNA.Prot_noFlip)<-c("RNA_Name", "RNA_ID", 
                                      "RNA.Pvalue.Gain", "RNA.Diff.Gain", 
                                      "RNA.Pvalue.Loss", "RNA.Diff.Loss", 
                                      "RNA.Pvalue.GainvLoss", "RNA.Diff.GainvLoss", 
                                      "Protein_ID", "Protein.Pvalue.Gain", 
                                      "Protein.Diff.Gain", "Protein.Pvalue.Loss", 
                                      "Protein.Diff.Loss", "Protein.Pvalue.GainvLoss", 
                                      "Protein.Diff.GainvLoss")

## add categories based on cutoffs.  (-Inf,-0.1,0.25,Inf)
CN.Diff.RNA.Prot_noFlip$Three.RNA.Gain<- cut(CN.Diff.RNA.Prot_noFlip$RNA.Diff.Gain,
                                              breaks=c(-Inf,-0.1,0.25,Inf),
                                              include.lowest=TRUE,
                                              labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_noFlip$Three.RNA.Loss<- cut(CN.Diff.RNA.Prot_noFlip$RNA.Diff.Loss,
                                              breaks=c(-Inf,-0.25,0.1,Inf),
                                              include.lowest=TRUE,
                                              labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.RNA.Prot_noFlip$Three.Protein.Gain<- cut(CN.Diff.RNA.Prot_noFlip$Protein.Diff.Gain,
                                                  breaks=c(-Inf,-0.1,0.25,Inf),
                                                  include.lowest=TRUE,
                                                  labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_noFlip$Three.Protein.Loss<- cut(CN.Diff.RNA.Prot_noFlip$Protein.Diff.Loss,
                                                  breaks=c(-Inf,-0.25,0.1,Inf),
                                                  include.lowest=TRUE,
                                                  labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.RNA.Prot_noFlip$Three.Protein.Loss<-factor(CN.Diff.RNA.Prot_noFlip$Three.Protein.Loss, 
                                                    levels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_noFlip$Three.RNA.Loss<-factor(CN.Diff.RNA.Prot_noFlip$Three.RNA.Loss, 
                                                levels=c("Anti-Scaling","Buffering","Scaling"))


#setwd(DataFileLocation)

write.csv(CN.Diff.RNA.Prot_noFlip, #change name as needed
          file =paste("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_noFlip_Min10Cells.csv", sep=','), 
          row.names = TRUE)
CN.Diff.RNA.Prot_noFlip<- read.csv("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_noFlip_Min10Cells.csv")


###### Control dataset: . Remove 20% of the "least reliable" (low reproducibility) genes from Upadhya & Ryan BioRxiv 2021 ####
# the dataset this section generates is available in supplementary data 3, sheet 6 "No_Low_Reproducibility" 

### Substep 1: Filter data for genes that are "more reliable" according to
## Upadhya & Ryan, bioRxiv 2021
## Experimental reproducibility limits the correlation between mRNA and protein abundances in tumour proteomic profiles
## data from supplemental table 2
## Aggregated Protein Reliability rank
## discard the 20% of genes with the "least reliable" scores
## Reliability score based on better protein expression reproducibility 
##       between TCGA and CCLE sample repeat expression similarity

## First correlate our protein & RNA expression dataset
## with Upadhya & Ryan, bioRxiv 2021 reliability scores
## then cut out the bottom 20% of genes 
## then find difference in expression for genes with top 80% reliability

Reliability_Scores<-read.xlsx(file="Upadhya.Ryan_BioRxiv_2021_S2_Protein_reliability.xlsx", 
                              sheetName = "B. Protein reproducibility rank", header=TRUE)
colnames(Reliability_Scores)<-c("Gene_Symbol", "Ovarian.Reproducibility.Rank", 
                                "Colon.Reproducibility.Rank", "CCLE.Reproducibility.Rank", 
                                "Aggregated.Reproducibility.Rank")

## remove genes that have lowest 20% of Reliability scores
## analyze genes with top 80% of Reliability scores

# Protein_ProID has 12755 proteins
Reliability_Scores.Info<- merge(x= Reliability_Scores, y= Protein_ProID, 
                          by.x="Gene_Symbol", by.y= "Gene_Symbol") # 12105
Reliability_Scores.Info<- Reliability_Scores.Info %>%
  mutate(quantile = ntile(Aggregated.Reproducibility.Rank, 10)) #split reliability scores into 10


# then split into high-low expression categories
LowReliability<- subset(Reliability_Scores.Info, quantile <= 2) #get all genes with lowest 20% RNA expression levels
HighReliability<- subset(Reliability_Scores.Info, quantile > 2) #get all genes with higher 80% RNA expression levels


# Filter RNA & protein data to get higher expression genes only
Protein.Expression_Reliable<- Protein.Expression.filtered_min10Cells[, colnames(Protein.Expression.filtered_min10Cells) %in% HighReliability$Protein_Id]
Protein.Expression_Reliable$Cell_Line<- Protein.Expression.filtered_min10Cells$Cell_line #also add back cell line info
# number of cells= ,  genes 

RNA.Expression_Reliable<- RNA.Expression.filtered[, colnames(RNA.Expression.filtered) %in% HighReliability$Gene_Symbol]
RNA.Expression_Reliable$Cell_Line<- RNA.Expression.filtered$Cell_line
# number of cells= ,  genes 
length(Protein.Expression_Reliable) #4008 genes
length(RNA.Expression_Reliable) # 3982




#### Step 2: Get difference data
## Now that I have dataset with no "unreliable" genes
## find difference in expression upon chrm gain and loss
## then find mean RNA/Protein difference in expression upon chrm gain/loss 
## also find number of genes

### sub-step 2: get protein difference in no mutation genes 
Protein.ID<-data.frame(Protein_ID=
                         colnames(Protein.Expression_Reliable[3:length(Protein.Expression_Reliable)]))



t.test_prot_category_Reliable<- data.frame(Protein_ID=character(),
                                         Protein_Name=character(),
                                         Pvlaue.Tri.Di= numeric(),
                                         Diff.Gain= numeric(), 
                                         Pvlaue.Di.Mono= numeric(),
                                         Diff.Di_Mono= numeric(),
                                         Pvlaue.Tri.Mono= numeric(),
                                         Diff.GainvLoss= numeric())


for (i in 1:(length(Protein.Expression_Reliable)-1)){
  TestProt = colnames(Protein.Expression_Reliable[i])
  
  testProtChrm <- filter(Depmap.Protein.info5, Protein_ID==TestProt) #get test protein data
  
  
  #use If statement to check that Protein Info has Protein of interest. 
  if (length(testProtChrm[,1])!=0){
    TestArm <- testProtChrm$arm #Find test protein chromosome arm
    TestChrm<- as.numeric(testProtChrm$Chromosome) #Find test protein chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test protein location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    Protein_Expression_TestProt <- Protein.Expression_Reliable %>% select(Cell_Line, all_of(TestProt))
    
    #Combine protein data with aneuploidy data for test protein location
    Protein.effect.Chrm<-merge(y= Protein_Expression_TestProt, x= testchrm.percell, 
                               by.y="Cell_Line", by.x="DepMap_ID", 
                               sort = TRUE)# 368 cells, 12757 proteins 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(Protein.effect.Chrm, arm_call==1) #all cells triploid for testProt chrm
    Chrm.di<- filter(Protein.effect.Chrm, arm_call==0) #all cells diploid for testProt chrm
    Chrm.mono<- filter(Protein.effect.Chrm, arm_call==-1) #all cells mono for testProt chrm
    
    ##Make data frame with t-test info about protein expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 & #[8] is column with protein expression data. check it has values. 
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 & #usually 10, now 3**
        colSums(!is.na(Chrm.mono[8]))>=10) { #usually 10, now 3**
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with Protein data
                     Chrm.di[,8], # [,8] because that is column with Protein data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Di.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                      Chrm.di[,8], # [,8] because that is column with Protein data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Di.Mono<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE)  #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                       Chrm.tri[,8], # [,8] because that is column with Protein data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_prot_category_Reliable<- rbind(t.test_prot_category_Reliable, 
                                          data.frame(Protein_ID=TestProt,
                                                     Protein_Name=testProtChrm$Gene_Symbol,
                                                     Pvlaue.Tri.Di= di.tri$p.value,
                                                     Diff.Gain= Diff.Tri.Di, 
                                                     Pvlaue.Di.Mono= Di.Mono$p.value,
                                                     Diff.Di_Mono= Diff.Di.Mono, 
                                                     Pvlaue.Tri.Mono= Tri.Mono$p.value,
                                                     Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_prot_category_Reliable<-t.test_prot_category_Reliable
    }
  }
}

t.test_prot_category_Reliable<-distinct(t.test_prot_category_Reliable) 

t.test_prot_category_Reliable<-t.test_prot_category_Reliable[order(t.test_prot_category_Reliable$Pvlaue.Tri.Mono),]# *** genes
# length (number of genes) = 4006 genes



### sub-step 3: get RNA difference in lung cells
# Lung cancer analysis: RNA expression difference: 
# prep data frame: 

# Run loop to get difference in experssion per gene in only low aneuploid cells
t.test_RNA_category_Reliable<- data.frame(RNA_ID=character(),
                                        RNA_Name=character(),
                                        Pvalue.Gain= numeric(),
                                        Diff.Gain= numeric(), 
                                        Pvalue.Loss= numeric(),
                                        Diff.Loss= numeric(),
                                        Pvalue.GainvLoss= numeric(),
                                        Diff.GainvLoss= numeric())

for (i in 2:length(RNA.Expression_Reliable)){
  
  TestRNA_name=sub("[..].*", "", as.character(colnames(RNA.Expression_Reliable[i])))
  TestRNA=colnames(RNA.Expression_Reliable[i])
  
  testRNAChrm <- filter(RNA_Info3, RNA_Info3$"Approved Symbol"==TestRNA) #get test RNA data
  if (length(testRNAChrm$Chromosome)!=0){ #Make sure RNA data is in RNA_info3, else it crashes
    TestArm <- testRNAChrm$arm #Find test RNA chromosome arm
    TestChrm<- as.numeric(testRNAChrm$Chromosome) #Find test RNA chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test RNA location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    RNA_Expression_TestRNA <- RNA.Expression_Reliable %>% select(Cell_Line, all_of(TestRNA))
    
    #Combine RNA data with aneuploidy data for test RNA location
    RNA.effect.Chrm<-merge(y= RNA_Expression_TestRNA, x= testchrm.percell, 
                           by.y="Cell_Line", by.x="DepMap_ID", 
                           sort = TRUE)# 368 cells, 12757 RNAs 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(RNA.effect.Chrm, arm_call==1) #all cells triploid for testRNA chrm
    Chrm.di<- filter(RNA.effect.Chrm, arm_call==0) #all cells diploid for testRNA chrm
    Chrm.mono<- filter(RNA.effect.Chrm, arm_call==-1) #all cells mono for testRNA chrm
    
    ##Make data frame with t-test info about RNA expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 &  
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 & 
        colSums(!is.na(Chrm.mono[8]))>=10) { 
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with RNA data
                     Chrm.di[,8], # [,8] because that is column with RNA data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Mono.Di<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                      Chrm.di[,8], # [,8] because that is column with RNA data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Mono.Di<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                       Chrm.tri[,8], # [,8] because that is column with RNA data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_RNA_category_Reliable<- rbind(t.test_RNA_category_Reliable, data.frame(
        RNA_ID=colnames(RNA.effect.Chrm[8]),
        RNA_Name=TestRNA_name,
        Pvalue.Gain= di.tri$p.value,
        Diff.Gain= Diff.Tri.Di, 
        Pvalue.Loss= Mono.Di$p.value,
        Diff.Loss= Diff.Mono.Di, 
        Pvalue.GainvLoss= Tri.Mono$p.value,
        Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_RNA_category_Reliable<-t.test_RNA_category_Reliable
    }
  }
}

t.test_RNA_category_Reliable<-distinct(t.test_RNA_category_Reliable) 
t.test_RNA_category_Reliable<-t.test_RNA_category_Reliable[order(t.test_RNA_category_Reliable$Pvalue.GainvLoss),]
# length (# of RNA genes) : 

### Sub step 4: combine low ploidy RNA and protein data
CN.Diff.RNA.Prot_Reliable<-merge(x=t.test_RNA_category_Reliable, 
                               y=t.test_prot_category_Reliable, 
                               by.x="RNA_Name", by.y="Protein_Name") #7874 genes

colnames(CN.Diff.RNA.Prot_Reliable)<-c("RNA_Name", "RNA_ID", 
                                     "RNA.Pvalue.Gain", "RNA.Diff.Gain", 
                                     "RNA.Pvalue.Loss", "RNA.Diff.Loss", 
                                     "RNA.Pvalue.GainvLoss", "RNA.Diff.GainvLoss", 
                                     "Protein_ID", "Protein.Pvalue.Gain", 
                                     "Protein.Diff.Gain", "Protein.Pvalue.Loss", 
                                     "Protein.Diff.Loss", "Protein.Pvalue.GainvLoss", 
                                     "Protein.Diff.GainvLoss")

## add categories based on cutoffs.  (-Inf,-0.1,0.25,Inf)
CN.Diff.RNA.Prot_Reliable$Three.RNA.Gain<- cut(CN.Diff.RNA.Prot_Reliable$RNA.Diff.Gain,
                                             breaks=c(-Inf,-0.1,0.25,Inf),
                                             include.lowest=TRUE,
                                             labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_Reliable$Three.RNA.Loss<- cut(CN.Diff.RNA.Prot_Reliable$RNA.Diff.Loss,
                                             breaks=c(-Inf,-0.25,0.1,Inf),
                                             include.lowest=TRUE,
                                             labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.RNA.Prot_Reliable$Three.Protein.Gain<- cut(CN.Diff.RNA.Prot_Reliable$Protein.Diff.Gain,
                                                 breaks=c(-Inf,-0.1,0.25,Inf),
                                                 include.lowest=TRUE,
                                                 labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_Reliable$Three.Protein.Loss<- cut(CN.Diff.RNA.Prot_Reliable$Protein.Diff.Loss,
                                                 breaks=c(-Inf,-0.25,0.1,Inf),
                                                 include.lowest=TRUE,
                                                 labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.RNA.Prot_Reliable$Three.Protein.Loss<-factor(CN.Diff.RNA.Prot_Reliable$Three.Protein.Loss, 
                                                   levels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_Reliable$Three.RNA.Loss<-factor(CN.Diff.RNA.Prot_Reliable$Three.RNA.Loss, 
                                               levels=c("Anti-Scaling","Buffering","Scaling"))

###Bar graph of RNA and Protein quantiles Scaling/Buffered

#format data so I can plot all 4 groups: RNA/Prot gain/loss
dat.m_reliable <- melt(CN.Diff.RNA.Prot_Reliable, id.vars='RNA_ID', 
                   measure.vars=c('Three.RNA.Gain','Three.Protein.Gain', 
                                  'Three.RNA.Loss', 'Three.Protein.Loss'))

## Barplot of scaling/buffering/anti-scaling categories
pdf(height=4, width=4, "plot.Protein.RNA.Bargraph.ThreeCat.Gain_Reliable.min10.pdf")
ggplot(data= dat.m_reliable, aes(x=variable, fill=value)) + 
  geom_bar(position="fill")+ #dodge (next to eachother)
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()
#4x4
# plot.Protein.RNA.Bargraph.ThreeCat.Gain_Lung.min10.pdf

#setwd(DataFileLocation)

write.csv(CN.Diff.RNA.Prot_Reliable, #change name as needed
          file =paste("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_Reliable_Min10Cells.csv", sep=','), 
          row.names = TRUE)

#CN.Diff.RNA.Prot_Reliable<- read.csv("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_Reliable_Min10Cells.csv")

###### Control dataset: . Merged: no mutations, flipped, low RNA & low reliability genes ####
# the dataset this section generates is available in supplementary data 3, sheet 7 "Merge_Control" 

### Remove mutated genes/cell and then remove flipped genes/cell
## and then remove genes in bottom 20% RNA expression and bottom 20% reliability genes

### sub-step 1: Get protein and RNA expression data set 
# Prep protein data
Protein.ID<-data.frame(Protein_ID=
                         colnames(Protein.Expression.filtered[3:length(Protein.Expression.filtered)]))

# Now merge depmap column name info with Protein info. 
# First add Uniprot ID and gene symbol info to depmap data. 
Depmap.Protein.info2<-merge(x= Protein.ID, y= Protein_ProID, 
                            by.x="Protein_ID", by.y="Protein_Id", 
                            sort = TRUE)# 3 collumns, 12755. all Proteins given Uniprot IDs.
# Now combine all genes with same uniprot ID: 
Depmap.Protein.info3<-merge(x= Depmap.Protein.info2, y= Protein_Info4, 
                            by.x="Uniprot_Acc", by.y="UniProt ID(supplied by UniProt)", 
                            sort = TRUE)# 10594 Proteins
# Find those genes without uniprot id: 
No_UniprotID <- anti_join(Depmap.Protein.info2, Protein_Info4, #finding genes_Symbol without match
                          by = c("Uniprot_Acc" = "UniProt ID(supplied by UniProt)"))
#...find the genes with matching gene symbols
Depmap.Protein.info4<-merge(x= No_UniprotID, y= Protein_Info4, 
                            by.x="Gene_Symbol", by.y="Approved Symbol", 
                            sort = TRUE)# 10 collumns, 253 genes, only no gene-symbol genes with uniprot_ID

# Merge genes with gene symbol and genes with only uniprot ID. 
Depmap.Protein.info5<-merge(x= Depmap.Protein.info3, y= Depmap.Protein.info4, 
                            all=TRUE)# 11 collumns, 12100 genes

setwd(DataFileLocation)
CN.Diff.xRNA.yProt.ThreeGroups<- read.csv("RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points_3cat.csv")

# Now: get protein data only for proteins with minimum 10 cells with protein data
Protein.Expression.filtered_min10Cells<- Protein.Expression.filtered %>% select(one_of(CN.Diff.xRNA.yProt.ThreeGroups$Protein_ID))
# Now add Cell_lines back in. 
Protein.Expression.filtered_min10Cells$Cell_line<- Protein.Expression.filtered$Cell_line


colnames(RNA.Expression.filtered)<-sub("[..].*", "", as.character(colnames(RNA.Expression.filtered)))#12755 genes


#### 1) First remove genes with lowest 20% reliability score
### See above code to get list of 80% reliable genes
# Filter RNA & protein data to get higher expression genes only
Protein.Expression_noMut_noFlip_noLowRNA_nolowRely<- Protein.Expression.filtered_min10Cells[, colnames(Protein.Expression.filtered_min10Cells) %in% HighReliability$Protein_Id]
Protein.Expression_noMut_noFlip_noLowRNA_nolowRely$Cell_Line<- Protein.Expression.filtered_min10Cells$Cell_line #also add back cell line info
# genes = 4008

RNA.Expression_noMut_noFlip_noLowRNA_nolowRely<- RNA.Expression.filtered[, colnames(RNA.Expression.filtered) %in% HighReliability$Gene_Symbol]
RNA.Expression_noMut_noFlip_noLowRNA_nolowRely$Cell_Line<- RNA.Expression.filtered$Cell_line
# genes = 3982

#### 2) next remove genes with the 20% of lowest RNA expression levels
### See above code to get list of 80% top RNA expressing genes
Protein.Expression_noMut_noFlip_noLowRNA_nolowRely2<- Protein.Expression_noMut_noFlip_noLowRNA_nolowRely[, colnames(Protein.Expression_noMut_noFlip_noLowRNA_nolowRely) %in% HighRNAExpress$Protein_Id]
Protein.Expression_noMut_noFlip_noLowRNA_nolowRely2$Cell_Line<- Protein.Expression_noMut_noFlip_noLowRNA_nolowRely$Cell_line #also add back cell line info
# number of cells= 371,  genes 3984

RNA.Expression_noMut_noFlip_noLowRNA_nolowRely2<- RNA.Expression_noMut_noFlip_noLowRNA_nolowRely[, colnames(RNA.Expression_noMut_noFlip_noLowRNA_nolowRely) %in% HighRNAExpress$RNA_Name]
RNA.Expression_noMut_noFlip_noLowRNA_nolowRely2$Cell_Line<- RNA.Expression_noMut_noFlip_noLowRNA_nolowRely$Cell_line
# Genes 3943


#### 3) Now remove mutated genes/cell #
# now go through mutations list and delete corresponding mutated protein in that cell line
Protein.Expression_noMut_noFlip_noLowRNA_nolowRely3<-Protein.Expression_noMut_noFlip_noLowRNA_nolowRely2
RNA.Expression_noMut_noFlip_noLowRNA_nolowRely3<- RNA.Expression_noMut_noFlip_noLowRNA_nolowRely2
RNA_Col<-colnames(RNA.Expression_noMut_noFlip_noLowRNA_nolowRely3)
Prot_Col<-colnames(Protein.Expression_noMut_noFlip_noLowRNA_nolowRely3)

for (i in 1:length(mutations.info9$Entrez_Gene_Id)){
  MutCell=mutations.info9$DepMap_ID[i] 
  MutProt=mutations.info9$Protein_Id[i]
  MutGene=mutations.info9$Hugo_Symbol[i]
  #Protein: get row name for row with cell line for protein expression
  MutRow<-rownames(Protein.Expression_noMut_noFlip_noLowRNA_nolowRely3[Protein.Expression_noMut_noFlip_noLowRNA_nolowRely3$Cell_Line==MutCell,])
  # Replace protein expression data that had mutated protein name in mutated cell line with NA:  
  if (MutProt %in% Prot_Col){
    Protein.Expression_noMut_noFlip_noLowRNA_nolowRely3[MutRow, MutProt]<-NA 
  }
  if (MutGene %in% RNA_Col){
    # Replace RNA expression data that had mutated RNA name in mutated cell line with NA:  
    RNA.Expression_noMut_noFlip_noLowRNA_nolowRely3[MutRow, MutGene]<-NA 
  }
}
# proteins: 3985 genes, 371 cell lines 
# RNA: 3948 genes, 371 cell lines 
length(Protein.Expression_noMut_noFlip_noLowRNA_nolowRely3)
length(RNA.Expression_noMut_noFlip_noLowRNA_nolowRely3)

#***#
#### 4) Now remove flipped genes/cell ##
# now go through flipped list and delete corresponding flipped protein in that cell line
setwd(DataFileLocation)
RNA.Expression_noFlip<- read.csv("RNA.expression.data_noFlipped.GeneCN.ChrmCN.csv")
Protein.Expression_noFlip2<- read.csv("Protein.expression.data_noFlipped.GeneCN.ChrmCN.csv")

Protein.Expression_noMut_noFlip_noLowRNA_nolowRely4<-Protein.Expression_noMut_noFlip_noLowRNA_nolowRely3
RNA.Expression_noMut_noFlip_noLowRNA_nolowRely4<- RNA.Expression_noMut_noFlip_noLowRNA_nolowRely3
RNA_Col<-colnames(RNA.Expression_noMut_noFlip_noLowRNA_nolowRely4)

#for proteins: get gene name from protein ID
#get protein expression with gene name labels
Prot_Col<- sub(".*[.]", "", as.character(colnames(Protein.Expression_noMut_noFlip_noLowRNA_nolowRely4)))
Prot_Col<- sub("_.*", "", Prot_Col)
Protein.Expression_noMut_noFlip_noLowRNA_nolowRely5<-Protein.Expression_noMut_noFlip_noLowRNA_nolowRely4
colnames(Protein.Expression_noMut_noFlip_noLowRNA_nolowRely5)<-Prot_Col

#now replace mutated genes/cell with NA in both protein and RNA datasets
for (i in 1:length(Flipped.info$Gene_ID)){
  FlipCell=as.character(Flipped.info$Cell_line[i])
  FlipGene=as.character(Flipped.info$Gene_ID[i])
  #Protein: get row name for row with cell line for protein expression
  FlipRow<-rownames(RNA.Expression_noFlip[RNA.Expression_noFlip$Cell_line==FlipCell,])
  # Replace protein expression data that had mutated protein name in mutated cell line with NA:  
  if (FlipGene %in% Prot_Col){
    Protein.Expression_noMut_noFlip_noLowRNA_nolowRely5[FlipRow, FlipGene]<-NA 
  }
  if (FlipGene %in% RNA_Col){
    # Replace RNA expression data that had mutated RNA name in mutated cell line with NA:  
    RNA.Expression_noMut_noFlip_noLowRNA_nolowRely4[FlipRow, FlipGene]<-NA 
  }
}
# proteins: 3985 genes
# RNA: 3948 genes
# This takes a few hours to run, FYI

#now add protein_ID's back to the protein dataframe
colnames(Protein.Expression_noMut_noFlip_noLowRNA_nolowRely5)<-colnames(Protein.Expression_noMut_noFlip_noLowRNA_nolowRely4)
RNA.Expression_noMut_noFlip_noLowRNA_nolowRely4



#setwd(DataFileLocation)
write.csv(RNA.Expression_noMut_noFlip_noLowRNA_nolowRely4, file="RNA.expression.data_noMut_noFlip_noLowRNA_noLowReliable.csv")
write.csv(Protein.Expression_noMut_noFlip_noLowRNA_nolowRely5, file="Protein.expression.data_noMut_noFlip_noLowRNA_noLowReliable.csv")

# Only analyze NSCLC cell lines?? (did not end up doing this)
# did not end up using this because it reduced the number of genes too much. 
# Protein.Expression_noMut_noFlip_noLowRNA_nolowRely5<-read.csv("Protein.expression.data_noMut_noFlip_noLowRNA_noLowReliable.csv")
# Protein.Expression_noMut_noFlip_noLowRNA_nolowRely5$Cell_line<-Protein.Expression.filtered$Cell_line
# Protein.Expression_noMut_noFlip_noLowRNA_nolowRely_NSCLC<- subset(Protein.Expression_noMut_noFlip_noLowRNA_nolowRely5, Protein.Expression_noMut_noFlip_noLowRNA_nolowRely5$Cell_line %in% lung.cells$DepMap_ID)
# only 927 proteins had 10+ datapoints per gain/neutral/loss condition. 
# less than a thousand genes, not enough datapoints to draw strong conclusions


### Sub-step 2: get protein difference in merges control set


t.test_prot_category_noMut_noFlip_noLowRNA_nolowRely<- data.frame(Protein_ID=character(),
                                                                  Protein_Name=character(),
                                                                  Pvlaue.Tri.Di= numeric(),
                                                                  Diff.Gain= numeric(), 
                                                                  Pvlaue.Di.Mono= numeric(),
                                                                  Diff.Di_Mono= numeric(),
                                                                  Pvlaue.Tri.Mono= numeric(),
                                                                  Diff.GainvLoss= numeric())


for (i in 1:(length(Protein.Expression_noMut_noFlip_noLowRNA_nolowRely)-1)){
  TestProt = colnames(Protein.Expression_noMut_noFlip_noLowRNA_nolowRely[i])
  
  testProtChrm <- filter(Depmap.Protein.info5, Protein_ID==TestProt) #get test protein data
  
  
  #use If statement to check that Protein Info has Protein of interest. 
  if (length(testProtChrm[,1])!=0){
    TestArm <- testProtChrm$arm #Find test protein chromosome arm
    TestChrm<- as.numeric(testProtChrm$Chromosome) #Find test protein chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test protein location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    Protein_Expression_TestProt <- Protein.Expression_noMut_noFlip_noLowRNA_nolowRely %>% select(Cell_Line, all_of(TestProt))
    
    #Combine protein data with aneuploidy data for test protein location
    Protein.effect.Chrm<-merge(y= Protein_Expression_TestProt, x= testchrm.percell, 
                               by.y="Cell_Line", by.x="DepMap_ID", 
                               sort = TRUE)# 368 cells, 12757 proteins 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(Protein.effect.Chrm, arm_call==1) #all cells triploid for testProt chrm
    Chrm.di<- filter(Protein.effect.Chrm, arm_call==0) #all cells diploid for testProt chrm
    Chrm.mono<- filter(Protein.effect.Chrm, arm_call==-1) #all cells mono for testProt chrm
    
    ##Make data frame with t-test info about protein expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 & #[8] is column with protein expression data. check it has values. 
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 & # min 10 datapoints/category
        colSums(!is.na(Chrm.mono[8]))>=10) { #
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with Protein data
                     Chrm.di[,8], # [,8] because that is column with Protein data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Di.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                      Chrm.di[,8], # [,8] because that is column with Protein data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Di.Mono<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE)  #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with Protein data
                       Chrm.tri[,8], # [,8] because that is column with Protein data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_prot_category_noMut_noFlip_noLowRNA_nolowRely<- rbind(t.test_prot_category_noMut_noFlip_noLowRNA_nolowRely, 
                                                                   data.frame(Protein_ID=TestProt,
                                                                              Protein_Name=testProtChrm$Gene_Symbol,
                                                                              Pvlaue.Tri.Di= di.tri$p.value,
                                                                              Diff.Gain= Diff.Tri.Di, 
                                                                              Pvlaue.Di.Mono= Di.Mono$p.value,
                                                                              Diff.Di_Mono= Diff.Di.Mono, 
                                                                              Pvlaue.Tri.Mono= Tri.Mono$p.value,
                                                                              Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_prot_category_noMut_noFlip_noLowRNA_nolowRely<-t.test_prot_category_noMut_noFlip_noLowRNA_nolowRely
    }
  }
}

length(t.test_prot_category_noMut_noFlip_noLowRNA_nolowRely)
# length (number of genes) = 4010
# minimum of 10 datapoints per category



### sub-step 3: get RNA difference
#RNA expression difference: 

#format RNA data collumn names
RNA.Expression_noMut_noFlip_noLowRNA_nolowRely2<-RNA.Expression_noMut_noFlip_noLowRNA_nolowRely
#colnames(RNA.Expression_noMut_noFlip_noLowRNA_nolowRely2)<-sub(".*[..]", "", as.character(colnames(RNA.Expression_noMut_noFlip_noLowRNA_nolowRely2)))
#colnames(RNA.Expression_noMut_noFlip_noLowRNA_nolowRely2)<-sub("_.*", "", as.character(colnames(RNA.Expression_noMut_noFlip_noLowRNA_nolowRely2)))


# Run loop to get difference in experssion per gene in only low aneuploid cells
#minimum 10 datapoints per category
t.test_RNA_category_noMut_noFlip_noLowRNA_nolowRely<- data.frame(RNA_ID=character(),
                                                                 RNA_Name=character(),
                                                                 Pvalue.Gain= numeric(),
                                                                 Diff.Gain= numeric(), 
                                                                 Pvalue.Loss= numeric(),
                                                                 Diff.Loss= numeric(),
                                                                 Pvalue.GainvLoss= numeric(),
                                                                 Diff.GainvLoss= numeric())

for (i in 1:(length(RNA.Expression_noMut_noFlip_noLowRNA_nolowRely2)-1)){
  
  TestRNA_name=colnames(RNA.Expression_noMut_noFlip_noLowRNA_nolowRely2[i])
  TestRNA=colnames(RNA.Expression_noMut_noFlip_noLowRNA_nolowRely2[i])
  
  testRNAChrm <- filter(RNA_Info3, RNA_Info3$"Approved Symbol"==TestRNA_name) #get test RNA data
  if (length(testRNAChrm$Chromosome)!=0){ #Make sure RNA data is in RNA_info3, else it crashes
    TestArm <- testRNAChrm$arm #Find test RNA chromosome arm
    TestChrm<- as.numeric(testRNAChrm$Chromosome) #Find test RNA chromosome
    #filter cell data and get only aneuploidy scores for chrm arm of test RNA location
    testchrm.percell <- filter(aneuploid, chrom==TestChrm) 
    testchrm.percell <- filter(testchrm.percell, arm==TestArm)
    
    RNA_Expression_TestRNA <- RNA.Expression_noMut_noFlip_noLowRNA_nolowRely2 %>% select(Cell_Line, all_of(TestRNA))
    
    #Combine RNA data with aneuploidy data for test RNA location
    RNA.effect.Chrm<-merge(y= RNA_Expression_TestRNA, x= testchrm.percell, 
                           by.y="Cell_Line", by.x="DepMap_ID", 
                           sort = TRUE)# 368 cells, 12757 RNAs 
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(RNA.effect.Chrm, arm_call==1) #all cells triploid for testRNA chrm
    Chrm.di<- filter(RNA.effect.Chrm, arm_call==0) #all cells diploid for testRNA chrm
    Chrm.mono<- filter(RNA.effect.Chrm, arm_call==-1) #all cells mono for testRNA chrm
    
    ##Make data frame with t-test info about RNA expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 &  
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 & #min 10 datapoints
        colSums(!is.na(Chrm.mono[8]))>=10) { 
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with RNA data
                     Chrm.di[,8], # [,8] because that is column with RNA data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Mono.Di<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                      Chrm.di[,8], # [,8] because that is column with RNA data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Mono.Di<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      # Tri vs monosomy
      Tri.Mono<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                       Chrm.tri[,8], # [,8] because that is column with RNA data
                       mu = 0, 
                       alt = "two.sided",
                       conf.level = 0.99) #get p-value
      Diff.Tri.Mono<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.mono[,8], na.rm=TRUE) #get difference
      
      t.test_RNA_category_noMut_noFlip_noLowRNA_nolowRely<- rbind(t.test_RNA_category_noMut_noFlip_noLowRNA_nolowRely, data.frame(
        RNA_ID=colnames(RNA.Expression_noMut_noFlip_noLowRNA_nolowRely[i]),
        RNA_Name=TestRNA_name,
        Pvalue.Gain= di.tri$p.value,
        Diff.Gain= Diff.Tri.Di, 
        Pvalue.Loss= Mono.Di$p.value,
        Diff.Loss= Diff.Mono.Di, 
        Pvalue.GainvLoss= Tri.Mono$p.value,
        Diff.GainvLoss= Diff.Tri.Mono))
    } else {
      t.test_RNA_category_noMut_noFlip_noLowRNA_nolowRely<-t.test_RNA_category_noMut_noFlip_noLowRNA_nolowRely
    }
  }
}

length(t.test_RNA_category_noMut_noFlip_noLowRNA_nolowRely$RNA_ID)
# length (# of RNA genes) :  3825
# Note: protein length = 4010



### Sub step 4: combine low ploidy RNA and protein data
CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely<-merge(x=t.test_RNA_category_noMut_noFlip_noLowRNA_nolowRely, 
                                                        y=t.test_prot_category_noMut_noFlip_noLowRNA_nolowRely, 
                                                        by.x="RNA_Name", by.y="Protein_Name") # genes

colnames(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely)<-c("RNA_Name", "RNA_ID", 
                                                              "RNA.Pvalue.Gain", "RNA.Diff.Gain", 
                                                              "RNA.Pvalue.Loss", "RNA.Diff.Loss", 
                                                              "RNA.Pvalue.GainvLoss", "RNA.Diff.GainvLoss", 
                                                              "Protein_ID", "Protein.Pvalue.Gain", 
                                                              "Protein.Diff.Gain", "Protein.Pvalue.Loss", 
                                                              "Protein.Diff.Loss", "Protein.Pvalue.GainvLoss", 
                                                              "Protein.Diff.GainvLoss")

## add categories based on cutoffs.  (-Inf,-0.1,0.25,Inf)
CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Three.RNA.Gain<- cut(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA.Diff.Gain,
                                                                      breaks=c(-Inf,-0.1,0.25,Inf),
                                                                      include.lowest=TRUE,
                                                                      labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Three.RNA.Loss<- cut(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA.Diff.Loss,
                                                                      breaks=c(-Inf,-0.25,0.1,Inf),
                                                                      include.lowest=TRUE,
                                                                      labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Three.Protein.Gain<- cut(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Protein.Diff.Gain,
                                                                          breaks=c(-Inf,-0.1,0.25,Inf),
                                                                          include.lowest=TRUE,
                                                                          labels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Three.Protein.Loss<- cut(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Protein.Diff.Loss,
                                                                          breaks=c(-Inf,-0.25,0.1,Inf),
                                                                          include.lowest=TRUE,
                                                                          labels=c("Scaling", "Buffering","Anti-Scaling"))
CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Three.Protein.Loss<-factor(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Three.Protein.Loss, 
                                                                            levels=c("Anti-Scaling","Buffering","Scaling"))
CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Three.RNA.Loss<-factor(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Three.RNA.Loss, 
                                                                        levels=c("Anti-Scaling","Buffering","Scaling"))


#setwd(DataFileLocation)

write.csv(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely, #change name as needed
          file =paste("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_noMut_noFlip_noLowRNA_nolowRep.min10cells.csv", sep=','), 
          row.names = TRUE)

write.csv(Protein.Expression_noMut_noFlip_noLowRNA_nolowRely, #change name as needed
          file =paste("Protein_Expression_noMut_noFlip_noLowRNA_nolowRely.csv", sep=','), 
          row.names = TRUE)
write.csv(RNA.Expression_noMut_noFlip_noLowRNA_nolowRely, #change name as needed
          file =paste("RNA_Expression_noMut_noFlip_noLowRNA_nolowRely.csv", sep=','), 
          row.names = TRUE)

Protein.Expression_noMut_noFlip_noLowRNA_nolowRely<-read.csv("Protein_Expression_noMut_noFlip_noLowRNA_nolowRely.csv")
RNA.Expression_noMut_noFlip_noLowRNA_nolowRely<-read.csv("RNA_Expression_noMut_noFlip_noLowRNA_nolowRely.csv")
CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely<-read.csv("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_noMut_noFlip_noLowRNA_nolowRep.min10cells.csv")

# mean difference and mean %change upon gain loss in rna/protein
# protein Gain: FC= 0.183, %change= 13.5
# protein loss: FC= -0.1355691, %change= -8.97%
# RNA Gain: FC= 0.3065611, %change= 23.7%
# RNA loss: FC= -0.2412671, %change= -15.4%


### For comparison, bulk analysis is: 
# CN.Diff.xRNA.yProt.ThreeGroups
# protein gain: FC= 0.1603291, %change= 11.75% 
# protein loss: FC = -0.1268382, %change = -8.42%


### Plot scatterplot and barplots of difference categories
##Scatterplots with categories colored in: 
#RNA Chrm gain graphs
ggplot(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely, 
       aes(x=RNA.Diff.Gain, y=-log2(RNA.Pvalue.Gain), color=Three.RNA.Gain))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in RNA expression")+
  ylab("log2(p-value)")+
  labs(color = "Quartiles:")+ #legend title
  theme_classic()+
  coord_cartesian(xlim=c(-3, 3), ylim=c(0,100))+
  ggtitle("RNA Gain scatterplot quartiles: NSCL noMut HighRNA")

#RNA Loss graphs
ggplot(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely, 
       aes(x=RNA.Diff.Loss, y=-log2(RNA.Pvalue.Loss), color=Three.RNA.Loss))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in RNA expression")+
  ylab("log2(p-value)")+
  labs(color = "Quartiles:")+ #legend title
  theme_classic()+
  coord_cartesian(xlim=c(-3, 3), ylim=c(0,100))+
  ggtitle("RNA Loss scatterplot quartiles: NSCL noMut HighRNA")

#Protein Chrm gain graphs
ggplot(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely, 
       aes(x=Protein.Diff.Gain, y=-log2(Protein.Pvalue.Gain), color=Three.Protein.Gain))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in protein expression")+
  ylab("log2(p-value)")+
  labs(color = "Quartiles:")+ #legend title
  theme_classic()+
  coord_cartesian(xlim=c(-3, 3), ylim=c(0,100))+
  ggtitle("Protein Gain scatterplot quartiles: NSCL noMut HighRNA")

#Protein Loss graphs
ggplot(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely, 
       aes(x=Protein.Diff.Loss, y=-log2(Protein.Pvalue.Loss), color=Three.Protein.Loss))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in protein expression")+
  ylab("log2(p-value)")+
  labs(color = "Quartiles:")+ #legend title
  theme_classic()+
  coord_cartesian(xlim=c(-3, 3), ylim=c(0,100))+
  ggtitle("protein Loss scatterplot quartiles: NSCL noMut HighRNA")


###Bar graph of RNA and Protein quantiles Scaling/Buffered

#format data so I can plot all 4 groups: RNA/Prot gain/loss
dat.m_merge <- melt(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely, id.vars='RNA_ID', 
                    measure.vars=c('Three.RNA.Gain','Three.Protein.Gain', 
                                   'Three.RNA.Loss', 'Three.Protein.Loss'))

## RNA & Protein Gain
ggplot(data= dat.m_merge, aes(x=variable, fill=value)) + 
  geom_bar(position="fill")+ #dodge (next to eachother)
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("RNA Difference upon chrm arm gain: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  #coord_cartesian(ylim=c(0,6000))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Relationship between RNA and Protein Difference\nupon chromosome arm gain: per Category")
# 4x4
# plot.Protein.RNA.Bargraph.ThreeCat.Gain_noMut_noFlip_noLowRNA_nolowRely

# number of genes per group: AS, Buffering, Scaling (total 4010)
table(subset(dat.m_merge, variable=="Three.RNA.Gain")$value) #178, 1309, 2523
# percent RNA gain: 4.4, 32.6, 62.9
table(subset(dat.m_merge, variable=="Three.Protein.Gain")$value) #178, 1309, 2523
# percent RNA gain: 6.4, 58.9, 34.7
table(subset(dat.m_merge, variable=="Three.RNA.Loss")$value) 
# percent RNA loss: 7.7, 36.4, 55.8
table(subset(dat.m_merge, variable=="Three.Protein.Loss")$value) #178, 1309, 2523
# percent RNA gain: 7.7, 67.5, 24.7


### RNA and Proteien expression correlation
## merged control sets
## Density plot of change in RNA & Protein upon chrm gain/Loss 
## using only data with 10+ datapoints per condition. 

CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely
# analyzed 7737 genes, 64 cell lines
# CHROMOSOME GAIN: correlate RNA and Protein expression difference
# Mimimum of 3 datapoints per condition (gain, no aneu, loss)
ggplot(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely, aes(x=CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA.Diff.Gain, 
                                                             y=CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Protein.Diff.Gain))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  xlab("Difference in NSCL RNA expression")+
  ylab("Difference in NSCL protein expression")+
  theme_classic()+
  geom_hline(yintercept=0.0)+
  geom_vline(xintercept=0.0)+
  geom_smooth(method="lm", color="Red")+
  coord_cartesian(xlim=c(-0.8, 0.8), ylim=c(-0.8, 0.8))+
  ggtitle("Chromosome gain")
#4x4
# plot.Protein.RNA.expression.density.Gain_noMut_noFlip_noLowRNA_nolowRely
Chrm_Gain_noMut_noFlip_noLowRNA_nolowRely<-cor.test(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Protein.Diff.Gain, 
                                                    CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA.Diff.Gain, 
                                                    method="pearson")
#cor=0.629
#P < 2E-16


# CHOMOSOME LOSS: correlate RNA and Protein expression difference
ggplot(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely, aes(x=RNA.Diff.Loss, y=Protein.Diff.Loss))+
  stat_density_2d(alpha=0.3, geom= "polygon", color="black",
                  aes() )+
  xlab("Difference in NSCL RNA expression")+
  ylab("Difference in NSCL protein expression")+
  theme_classic()+
  geom_hline(yintercept=0.0)+
  geom_vline(xintercept=0.0)+
  geom_smooth(method="lm", color="Red")+
  coord_cartesian(xlim=c(-0.7, 0.7), ylim=c(-0.7,0.7))+
  ggtitle("Chromosome loss")
#4x4
#plot.Protein.RNA.expression.density_noMut_noFlip_noLowRNA_nolowRely
Chrm_Loss_Merge<-cor.test(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA.Diff.Loss, 
                          CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Protein.Diff.Loss, 
                          method="pearson")
Chrm_Loss_Merge
#cor=0.654
#P<2E-16

## RNA & Protein Gain
ggplot(data= dat.m_merge, aes(x=variable, fill=value)) + 
  geom_bar(position="fill")+ #dodge (next to eachother)
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("RNA Difference upon chrm arm gain: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  #coord_cartesian(ylim=c(0,6000))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Relationship between RNA and Protein Difference\nupon chromosome arm gain: per Category")
#4x4
# plot.Protein.RNA.Bargraph.ThreeCat.Gain_noMut_noFlip_noLowRNA_nolowRely



###### Plot & overview of expression difference per data control #####
# data generated in this section is available in supplementary figure 3, sheet 2. 

## I want to have an overview of average RNA/Protein expression changes upon gain/loss
## from all the different above filtered datasets
## minimum 10 datapoints per category (gain loss neutral) for all
## datasets I will use: 
# CN.Diff.RNA.Prot_Reliable
# CN.Diff.RNA.Prot_noMut
# CN.Diff.RNA.Prot_noFlip
# CN.Diff.RNA.Prot_highRNA
# CN.Diff.RNA.Prot_Lung
# CN.Diff.xRNA.yProt
# CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely

# Note that log2 fold change (difference score) is not the percent difference!!! 
## Percent changes are available in supplementary figure 3, sheet 2. 


## Plot 1: Heatmap of mean expression per dataset & gene difference
Mean.Diff.Controls<-data.frame(Data=as.character(), 
                               Num.Genes=as.numeric(), 
                               Num.Cells=as.numeric(), 
                               Diff.Gain.RNA=as.numeric(), 
                               Diff.Gain.Protein=as.numeric(),
                               Diff.Loss.RNA=as.numeric(), 
                               Diff.Loss.Protein=as.numeric())
# All data
Mean.Diff.Controls<-rbind(Mean.Diff.Controls, 
                          data.frame(Data= "All",
                            Num.Genes=length(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name), 
                            Num.Cells=as.numeric(367), 
                            Diff.Gain.RNA=mean(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain), 
                            Diff.Gain.Protein=mean(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain), 
                            Diff.Loss.RNA=mean(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss), 
                            Diff.Loss.Protein=mean(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss)  ))

# No Mutations 
Mean.Diff.Controls<-rbind(Mean.Diff.Controls, 
                          data.frame(Data= "No mutations",
                                     Num.Genes=length(CN.Diff.RNA.Prot_noMut$RNA_Name), 
                                     Num.Cells=as.numeric(367), 
                                     Diff.Gain.RNA=mean(CN.Diff.RNA.Prot_noMut$RNA.Diff.Gain, na.rm=TRUE), 
                                     Diff.Gain.Protein=mean(CN.Diff.RNA.Prot_noMut$Protein.Diff.Gain, na.rm=TRUE), 
                                     Diff.Loss.RNA=mean(CN.Diff.RNA.Prot_noMut$RNA.Diff.Loss , na.rm=TRUE), 
                                     Diff.Loss.Protein=mean(CN.Diff.RNA.Prot_noMut$Protein.Diff.Loss, na.rm=TRUE)  ))

# Flipped Gene CV 
Mean.Diff.Controls<-rbind(Mean.Diff.Controls, 
                          data.frame(Data= "No flipped genes",
                                     Num.Genes= length(CN.Diff.RNA.Prot_noFlip$RNA_Name), 
                                     Num.Cells= as.numeric(367), 
                                     Diff.Gain.RNA= mean(CN.Diff.RNA.Prot_noFlip$RNA.Diff.Gain), 
                                     Diff.Gain.Protein= mean(CN.Diff.RNA.Prot_noFlip$Protein.Diff.Gain), 
                                     Diff.Loss.RNA= mean(CN.Diff.RNA.Prot_noFlip$RNA.Diff.Loss), 
                                     Diff.Loss.Protein= mean(CN.Diff.RNA.Prot_noFlip$Protein.Diff.Loss)  ))

# no low RNA expression 
Mean.Diff.Controls<-rbind(Mean.Diff.Controls, 
                          data.frame(Data= "No low RNA expression",
                                     Num.Genes=length(CN.Diff.RNA.Prot_highRNA$RNA_Name), 
                                     Num.Cells=as.numeric(367), 
                                     Diff.Gain.RNA=mean(CN.Diff.RNA.Prot_highRNA$RNA.Diff.Gain), 
                                     Diff.Gain.Protein=mean(CN.Diff.RNA.Prot_highRNA$Protein.Diff.Gain), 
                                     Diff.Loss.RNA=mean(CN.Diff.RNA.Prot_highRNA$RNA.Diff.Loss), 
                                     Diff.Loss.Protein=mean(CN.Diff.RNA.Prot_highRNA$Protein.Diff.Loss)  ))

# No least reliable genes 
Mean.Diff.Controls<-rbind(Mean.Diff.Controls, 
                          data.frame(Data= "No low reproducibility genes",
                                     Num.Genes=length(CN.Diff.RNA.Prot_Reliable$RNA_Name), 
                                     Num.Cells= as.numeric(367), 
                                     Diff.Gain.RNA=mean(CN.Diff.RNA.Prot_Reliable$RNA.Diff.Gain), 
                                     Diff.Gain.Protein=mean(CN.Diff.RNA.Prot_Reliable$Protein.Diff.Gain), 
                                     Diff.Loss.RNA=mean(CN.Diff.RNA.Prot_Reliable$RNA.Diff.Loss), 
                                     Diff.Loss.Protein=mean(CN.Diff.RNA.Prot_Reliable$Protein.Diff.Loss)  ))

# Merge
Mean.Diff.Controls<-rbind(Mean.Diff.Controls, 
                          data.frame(Data= "Merge: no mutations, flipped, \nlow RNA or low reproducibility genes",
                                     Num.Genes=length(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA_Name), 
                                     Num.Cells= as.numeric(367), 
                                     Diff.Gain.RNA=mean(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA.Diff.Gain), 
                                     Diff.Gain.Protein=mean(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Protein.Diff.Gain), 
                                     Diff.Loss.RNA=mean(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA.Diff.Loss), 
                                     Diff.Loss.Protein=mean(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Protein.Diff.Loss)  ))



# NSCLC 
Mean.Diff.Controls<-rbind(Mean.Diff.Controls, 
                          data.frame(Data= "NSCLC",
                                     Num.Genes=length(CN.Diff.RNA.Prot_Lung$RNA_Name), 
                                     Num.Cells=as.numeric(64), 
                                     Diff.Gain.RNA=mean(CN.Diff.RNA.Prot_Lung$RNA.Diff.Gain), 
                                     Diff.Gain.Protein=mean(CN.Diff.RNA.Prot_Lung$Protein.Diff.Gain), 
                                     Diff.Loss.RNA=mean(CN.Diff.RNA.Prot_Lung$RNA.Diff.Loss), 
                                     Diff.Loss.Protein=mean(CN.Diff.RNA.Prot_Lung$Protein.Diff.Loss)  ))

write.csv(Mean.Diff.Controls, "Table.Overview.NumGenes.meanDiff.perSingleControl_withMerge.csv")
# Mean.Diff.Controls<- read.csv("Table.Overview.NumGenes.meanDiff.perSingleControl_withMerge.csv")

# Format data for a heatmao plot
Differences.control.m<- melt(Mean.Diff.Controls[c(1,4:7)])
Differences.control.m$variable<-factor(Differences.control.m$variable, 
                                          levels=c("Diff.Loss.Protein", "Diff.Loss.RNA", 
                                                   "Diff.Gain.Protein", "Diff.Gain.RNA"))


ggplot(Differences.control.m, aes(x=Data, y=variable))+ 
  geom_raster(aes(fill = value), hjust=0.5, vjust=0.5, interpolate=FALSE)+
  xlab("")+
  ylab("Mean Difference in expression")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  scale_fill_gradient2(low="dodgerblue3", mid="white", high="Red", midpoint=0, 
                       name="Difference", limits=c(-0.4, 0.4)) +
  #coord_flip()+
  ggtitle("Difference in gene expression upon \n chromosome gain or loss \n by filtered dataset")
#5x5 
# Plot.heatmap.diff.byControl

## Plot 2: barplot of difference: 
Differences.control.m$variable<-factor(Differences.control.m$variable, 
                                       levels= c("Diff.Gain.RNA", "Diff.Gain.Protein", 
                                                 "Diff.Loss.RNA", "Diff.Loss.Protein"))

pdf(height=4, width=6, file="plot.bar.MeanDiff.perControl_withMerge.pdf")
ggplot(Differences.control.m, aes(x=factor(variable), y=value))+ 
  geom_bar(position="dodge", stat = "identity", aes(fill = Data))+
  xlab("")+
  ylab("Mean Difference in expression")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  scale_fill_manual(name="Dataset", 
                    values=c("grey12", "grey25", "grey37", "grey50", "grey67", "grey75", "grey87"))+
  ggtitle("Difference in gene expression \n by control dataset")
dev.off()


##Plot 3: 
## boxplot of delta difference relative to all datapoints dataset
Mean.Diff.Controls2<- Mean.Diff.Controls
Mean.Diff.Controls2$Diff.Gain.RNA<- Mean.Diff.Controls$Diff.Gain.RNA/Mean.Diff.Controls$Diff.Gain.RNA[1]
Mean.Diff.Controls2$Diff.Gain.Protein<- Mean.Diff.Controls$Diff.Gain.Protein/Mean.Diff.Controls$Diff.Gain.Protein[1]
Mean.Diff.Controls2$Diff.Loss.RNA<- Mean.Diff.Controls$Diff.Loss.RNA/Mean.Diff.Controls$Diff.Loss.RNA[1]
Mean.Diff.Controls2$Diff.Loss.Protein<- Mean.Diff.Controls$Diff.Loss.Protein/Mean.Diff.Controls$Diff.Loss.Protein[1]


Delta.control.m<- melt(Mean.Diff.Controls2[c(1,4:7)])
Delta.control.m$variable<-factor(Delta.control.m$variable, 
                                       levels=c("Diff.Gain.RNA", "Diff.Gain.Protein", 
                                                "Diff.Loss.RNA", "Diff.Loss.Protein"))
pdf(height=4, width=6, file="plot.bar.Delta.MeanDiff.perControl_withMerge.pdf")
ggplot(Delta.control.m, aes(x=factor(variable), y=value))+ 
  geom_bar(position="dodge", stat = "identity", aes(fill = Data))+
  xlab("")+
  ylab("Delta (Difference in expression)")+
  theme_classic()+
  coord_cartesian(ylim=c(0.85,1.5))+
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  scale_fill_manual(name="Dataset", 
                    values=c("grey12", "grey25", "grey37", "grey50", "grey67", "grey75", "grey87"))+
  ggtitle("Difference in gene expression \n by control dataset")
dev.off()

### Plot 4: boxplot of all datapoints
# get dataframe of dataset, gene difference, and datapoints 
# get data from all control datasets, for all gene differences (RNA/Protein gain/loss)
## Note: datapoints are not merged by gene name, all datapoints used

## datasets I will use: 
# CN.Diff.xRNA.yProt.ThreeGroups
# CN.Diff.RNA.Prot_noMut <- read.csv()
# CN.Diff.RNA.Prot_noFlip <- read.csv()
# CN.Diff.RNA.Prot_highRNA <- read.csv()
# CN.Diff.RNA.Prot_Reliable <- read.csv()
# CN.Diff.RNA.Prot_Lung <- read.csv()
# CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely 


Diff.Data.PerControl<- data.frame(Dataset= as.character(), 
                                  DifferenceType= as.character(), 
                                  Difference = as.numeric())
# All # CN.Diff.xRNA.yProt.ThreeGroups
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="All", 
                                        DifferenceType= "RNA.Gain", 
                                        Difference= CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="All", 
                                        DifferenceType= "Protein.Gain", 
                                        Difference= CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="All", 
                                        DifferenceType= "RNA.Loss", 
                                        Difference= CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="All", 
                                        DifferenceType= "Protein.Loss", 
                                        Difference= CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss ))
# no mutations # CN.Diff.RNA.Prot_noMut
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="No mutations", 
                                        DifferenceType= "RNA.Gain", 
                                        Difference= CN.Diff.RNA.Prot_noMut$RNA.Diff.Gain ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="No mutations", 
                                        DifferenceType= "Protein.Gain", 
                                        Difference= CN.Diff.RNA.Prot_noMut$Protein.Diff.Gain ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="No mutations", 
                                        DifferenceType= "RNA.Loss", 
                                        Difference= CN.Diff.RNA.Prot_noMut$RNA.Diff.Loss ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="No mutations", 
                                        DifferenceType= "Protein.Loss", 
                                        Difference= CN.Diff.RNA.Prot_noMut$Protein.Diff.Loss ))
# no Flipped # CN.Diff.RNA.Prot_noFlip
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="No flipped genes", 
                                        DifferenceType= "RNA.Gain", 
                                        Difference= CN.Diff.RNA.Prot_noFlip$RNA.Diff.Gain ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="No flipped genes", 
                                        DifferenceType= "Protein.Gain", 
                                        Difference= CN.Diff.RNA.Prot_noFlip$Protein.Diff.Gain ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="No flipped genes", 
                                        DifferenceType= "RNA.Loss", 
                                        Difference= CN.Diff.RNA.Prot_noFlip$RNA.Diff.Loss ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="No flipped genes", 
                                        DifferenceType= "Protein.Loss", 
                                        Difference= CN.Diff.RNA.Prot_noFlip$Protein.Diff.Loss ))

# No low RNA only # CN.Diff.RNA.Prot_highRNA
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="High RNA only", 
                                        DifferenceType= "RNA.Gain", 
                                        Difference= CN.Diff.RNA.Prot_highRNA$RNA.Diff.Gain ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="High RNA only", 
                                        DifferenceType= "Protein.Gain", 
                                        Difference= CN.Diff.RNA.Prot_highRNA$Protein.Diff.Gain ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="High RNA only", 
                                        DifferenceType= "RNA.Loss", 
                                        Difference= CN.Diff.RNA.Prot_highRNA$RNA.Diff.Loss ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="High RNA only", 
                                        DifferenceType= "Protein.Loss", 
                                        Difference= CN.Diff.RNA.Prot_highRNA$Protein.Diff.Loss ))

# Reliable Genes only # CN.Diff.RNA.Prot_Reliable
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="High reproducibility genes", 
                                        DifferenceType= "RNA.Gain", 
                                        Difference= CN.Diff.RNA.Prot_Reliable$RNA.Diff.Gain ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="High reproducibility genes", 
                                        DifferenceType= "Protein.Gain", 
                                        Difference= CN.Diff.RNA.Prot_Reliable$Protein.Diff.Gain ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="High reproducibility genes", 
                                        DifferenceType= "RNA.Loss", 
                                        Difference= CN.Diff.RNA.Prot_Reliable$RNA.Diff.Loss ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="High reproducibility genes", 
                                        DifferenceType= "Protein.Loss", 
                                        Difference= CN.Diff.RNA.Prot_Reliable$Protein.Diff.Loss ))

# Merge no mutations, no flipped, no low RNA, no low reproducible # CN.Diff.RNA.Prot_Reliable
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="Merge", 
                                        DifferenceType= "RNA.Gain", 
                                        Difference= CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA.Diff.Gain ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="Merge", 
                                        DifferenceType= "Protein.Gain", 
                                        Difference= CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Protein.Diff.Gain ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="Merge", 
                                        DifferenceType= "RNA.Loss", 
                                        Difference= CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA.Diff.Loss ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="Merge", 
                                        DifferenceType= "Protein.Loss", 
                                        Difference= CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Protein.Diff.Loss ))

# NSCLC    # CN.Diff.RNA.Prot_Lung
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="NSCLC", 
                                        DifferenceType= "RNA.Gain", 
                                        Difference= CN.Diff.RNA.Prot_Lung$RNA.Diff.Gain ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="NSCLC", 
                                        DifferenceType= "Protein.Gain", 
                                        Difference= CN.Diff.RNA.Prot_Lung$Protein.Diff.Gain ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="NSCLC", 
                                        DifferenceType= "RNA.Loss", 
                                        Difference= CN.Diff.RNA.Prot_Lung$RNA.Diff.Loss ))
Diff.Data.PerControl<- rbind(Diff.Data.PerControl, 
                             data.frame(Dataset="NSCLC", 
                                        DifferenceType= "Protein.Loss", 
                                        Difference= CN.Diff.RNA.Prot_Lung$Protein.Diff.Loss ))

## Plot 5: boxplots with errorbars for difference in gene expression
# Diff.Data.PerControl (length 172644 datapoints) (haha, that's so many genes) 
#     it's a lot because it's all the datasets on top of eachother
# column names "Dataset"        "DifferenceType" "Difference"

pdf(height=4, width=6, file="plot.bax.MeanDiff.perControl.NoOutlier_withMerge.pdf")
ggplot(Diff.Data.PerControl, aes(x=factor(DifferenceType), y=Difference))+ 
  geom_boxplot(aes(fill = Dataset), outlier.shape=NA)+ # no outliers
  xlab("")+
  ylab("Difference in expression")+
  theme_classic()+
  coord_cartesian(ylim=c(-1,1))+
  geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  scale_fill_manual(name="Dataset", 
                    values=c("grey12", "grey25", "grey37", "grey50", "grey62", "grey75", "grey90"))+
  ggtitle("Difference in gene expression by control dataset")
dev.off()
# 4x6

# with outliers
pdf(height=4, width=6, file="plot.bax.MeanDiff.perControl.Outlier.pdf")
ggplot(Diff.Data.PerControl, aes(x=factor(DifferenceType), y=Difference))+ 
  geom_boxplot(aes(fill = Dataset), outlier.size = 0.5)+ # no outliers
  xlab("")+
  ylab("Difference in expression")+
  theme_classic()+
  #coord_cartesian(ylim=c(-1,1))+
  geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  scale_fill_manual(name="Dataset", 
                    values=c("grey15", "grey31", "grey47", "grey62", "grey78", "grey90"))+
  ggtitle("Difference in gene expression by control dataset")
dev.off()


## Do t-tests to see if significant differences
## sig if less than 0.01 <-- ! Note! 
## corrected using bonferroni correction (24 tests,so multiply by 24) 
## control vs NSCLC 
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain, 
       CN.Diff.RNA.Prot_Lung$RNA.Diff.Gain) #P 0.001269 --> 0.02538  NS 
t.test(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain, 
       CN.Diff.RNA.Prot_Lung$Protein.Diff.Gain) #P = 0.007042 --> NS
x<-t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss, 
       CN.Diff.RNA.Prot_Lung$RNA.Diff.Loss) #P < 2E-16 --> corrected < 2E-16 ***
x<-t.test(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss, 
       CN.Diff.RNA.Prot_Lung$Protein.Diff.Loss) #P 8E-06  --> corrected 2E-04 **

## control vs High reproducibility genes 
x<-t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain, 
       CN.Diff.RNA.Prot_Reliable$RNA.Diff.Gain) #P = 2.438E-06 --> 6E-05 ***
x<-t.test(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain, 
       CN.Diff.RNA.Prot_Reliable$Protein.Diff.Gain) #P = 8.48E-08 --> 2E-06 ***
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss, 
       CN.Diff.RNA.Prot_Reliable$RNA.Diff.Loss) #P = NS
t.test(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss, 
       CN.Diff.RNA.Prot_Reliable$Protein.Diff.Loss) #P = 0.017 --> corrected NS

## control vs high RNA express genes 
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain, 
       CN.Diff.RNA.Prot_highRNA$RNA.Diff.Gain) #P = NS
t.test(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain, 
       CN.Diff.RNA.Prot_highRNA$Protein.Diff.Gain) #P = NS
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss, 
       CN.Diff.RNA.Prot_highRNA$RNA.Diff.Loss) #P = NS
t.test(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss, 
       CN.Diff.RNA.Prot_highRNA$Protein.Diff.Loss) #P = NS

## control vs NO FLIPPED genes 
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain, 
       CN.Diff.RNA.Prot_noFlip$RNA.Diff.Gain) #P = 0.01787 --> corrected NS
t.test(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain, 
       CN.Diff.RNA.Prot_noFlip$Protein.Diff.Gain) #P = NS
x<- t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss, 
       CN.Diff.RNA.Prot_noFlip$RNA.Diff.Loss) #P = 2.766E-13--> corrected 7E-12 ***
t.test(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss, 
       CN.Diff.RNA.Prot_noFlip$Protein.Diff.Loss) #P = 0.02624  --> corrected NS

## control vs no mutated genes 
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain, 
       CN.Diff.RNA.Prot_noMut$RNA.Diff.Gain) #P = NS
t.test(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain, 
       CN.Diff.RNA.Prot_noMut$Protein.Diff.Gain) #P = NS
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss, 
       CN.Diff.RNA.Prot_noMut$RNA.Diff.Loss) #P = NS
t.test(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss, 
       CN.Diff.RNA.Prot_noMut$Protein.Diff.Loss) #P = NS

## control vs merge (no flip, no mut, no low rna, no low reproducibility) genes 
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain, 
       CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA.Diff.Gain) #P = 3.1E-06 --> 7E-05
t.test(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain, 
       CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Protein.Diff.Gain) #P = 9.8E-08 --> 2E-06
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss, 
       CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA.Diff.Loss) #P = NS
t.test(CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss, 
       CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Protein.Diff.Loss) #P = 0.01753 --> NS


###### Calculate and plot RNA-Protein expression delta per control dataset #####
## datasets I will use: 
# CN.Diff.RNA.Prot_Reliable
# CN.Diff.RNA.Prot_noMut
# CN.Diff.RNA.Prot_noFlip
# CN.Diff.RNA.Prot_highRNA
# CN.Diff.RNA.Prot_Lung
# CN.Diff.xRNA.yProt


Delta.Data.PerControl<- data.frame(Dataset= as.character(), 
                                  DifferenceType= as.character(), 
                                  Delta = as.numeric())
# All # CN.Diff.xRNA.yProt.ThreeGroups
Delta.Data.PerControl<- rbind(Delta.Data.PerControl, 
                             data.frame(Dataset="All", 
                                        DifferenceType= "Gain", 
                                        Delta= CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain - CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain))
Delta.Data.PerControl<- rbind(Delta.Data.PerControl, 
                             data.frame(Dataset="All", 
                                        DifferenceType= "Loss", 
                                        Delta= CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss - CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss))

# no mutations # CN.Diff.RNA.Prot_noMut
Delta.Data.PerControl<- rbind(Delta.Data.PerControl, 
                             data.frame(Dataset="No mutations", 
                                        DifferenceType= "Gain", 
                                        Delta= CN.Diff.RNA.Prot_noMut$RNA.Diff.Gain - CN.Diff.RNA.Prot_noMut$Protein.Diff.Gain))
Delta.Data.PerControl<- rbind(Delta.Data.PerControl, 
                             data.frame(Dataset="No mutations", 
                                        DifferenceType= "Loss", 
                                        Delta= CN.Diff.RNA.Prot_noMut$RNA.Diff.Loss - CN.Diff.RNA.Prot_noMut$Protein.Diff.Loss))

# no Flipped # CN.Diff.RNA.Prot_noFlip
Delta.Data.PerControl<- rbind(Delta.Data.PerControl, 
                             data.frame(Dataset="No flipped genes", 
                                        DifferenceType= "Gain", 
                                        Delta= CN.Diff.RNA.Prot_noFlip$RNA.Diff.Gain - CN.Diff.RNA.Prot_noFlip$Protein.Diff.Gain))
Delta.Data.PerControl<- rbind(Delta.Data.PerControl, 
                             data.frame(Dataset="No flipped genes", 
                                        DifferenceType= "Loss", 
                                        Delta= CN.Diff.RNA.Prot_noFlip$RNA.Diff.Loss - CN.Diff.RNA.Prot_noFlip$Protein.Diff.Loss))

# High RNA only # CN.Diff.RNA.Prot_highRNA
Delta.Data.PerControl<- rbind(Delta.Data.PerControl, 
                             data.frame(Dataset="No low RNA genes", 
                                        DifferenceType= "Gain", 
                                        Delta= CN.Diff.RNA.Prot_highRNA$RNA.Diff.Gain - CN.Diff.RNA.Prot_highRNA$Protein.Diff.Gain))
Delta.Data.PerControl<- rbind(Delta.Data.PerControl, 
                             data.frame(Dataset="No low RNA genes", 
                                        DifferenceType= "Loss", 
                                        Delta= CN.Diff.RNA.Prot_highRNA$RNA.Diff.Loss - CN.Diff.RNA.Prot_highRNA$Protein.Diff.Loss))

# Reliable Genes only # CN.Diff.RNA.Prot_Reliable
Delta.Data.PerControl<- rbind(Delta.Data.PerControl, 
                             data.frame(Dataset="No low reproducibility genes", 
                                        DifferenceType= "Gain", 
                                        Delta= CN.Diff.RNA.Prot_Reliable$RNA.Diff.Gain - CN.Diff.RNA.Prot_Reliable$Protein.Diff.Gain))
Delta.Data.PerControl<- rbind(Delta.Data.PerControl, 
                             data.frame(Dataset="No low reproducibility genes", 
                                        DifferenceType= "Loss", 
                                        Delta= CN.Diff.RNA.Prot_Reliable$RNA.Diff.Loss - CN.Diff.RNA.Prot_Reliable$Protein.Diff.Loss))

# Merge no mutation, no flip, no low RNA, no low reproduce # CN.Diff.RNA.Prot_Reliable
Delta.Data.PerControl<- rbind(Delta.Data.PerControl, 
                              data.frame(Dataset="Merge", 
                                         DifferenceType= "Gain", 
                                         Delta= CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA.Diff.Gain - CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Protein.Diff.Gain))
Delta.Data.PerControl<- rbind(Delta.Data.PerControl, 
                              data.frame(Dataset="Merge", 
                                         DifferenceType= "Loss", 
                                         Delta= CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA.Diff.Loss - CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Protein.Diff.Loss))

# NSCLC    # CN.Diff.RNA.Prot_Lung
Delta.Data.PerControl<- rbind(Delta.Data.PerControl, 
                             data.frame(Dataset="NSCLC", 
                                        DifferenceType= "Gain", 
                                        Delta= CN.Diff.RNA.Prot_Lung$RNA.Diff.Gain - CN.Diff.RNA.Prot_Lung$Protein.Diff.Gain))
Delta.Data.PerControl<- rbind(Delta.Data.PerControl, 
                             data.frame(Dataset="NSCLC", 
                                        DifferenceType= "Loss", 
                                        Delta= CN.Diff.RNA.Prot_Lung$RNA.Diff.Loss - CN.Diff.RNA.Prot_Lung$Protein.Diff.Loss))


## Plot 5: boxplots with errorbars for difference in gene expression
# Diff.Data.PerControl (length 172,644 datapoints) (haha, that's so many genes) 
#     it's a lot because it's all the datasets on top of eachother

# with outliers
pdf(height=4, width=5, file="plot.bax.Delta.RNAProt.perControl_withMerge.pdf")
ggplot(Delta.Data.PerControl, aes(x=factor(DifferenceType), y=Delta))+ 
  geom_boxplot(aes(fill = Dataset), outlier.shape = NA)+ # no outliers
  xlab("")+
  ylab("RNA difference - Protein difference\n upon aneuploidy")+
  theme_classic()+
  coord_cartesian(ylim=c(-1,1))+
  geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, color="black"), 
        axis.text.y = element_text(color="black"))+
  scale_fill_manual(name="Dataset", 
                    values=c("grey12", "grey25", "grey37", "grey50", "grey62", "grey75", "grey90"))+
  ggtitle("Delta in RNA-Protein expression differnce by control dataset")
dev.off()



### test for significance
## Do t-tests to see if significant differences
## sig if less than 0.01 <-- ! Note! 
## corrected using bonferroni correction (12 tests,so multiply by 12) 
## control vs NSCLC 
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain - CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain, 
       CN.Diff.RNA.Prot_Lung$RNA.Diff.Gain - CN.Diff.RNA.Prot_Lung$Protein.Diff.Gain) 
#P = 3.3E-09
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss - CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss, 
       CN.Diff.RNA.Prot_Lung$RNA.Diff.Loss - CN.Diff.RNA.Prot_Lung$Protein.Diff.Loss) 
#P < 2.2E-16 

## control vs High reproducibility genes 
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain - CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain, 
          CN.Diff.RNA.Prot_Reliable$RNA.Diff.Gain - CN.Diff.RNA.Prot_Reliable$Protein.Diff.Gain) 
#P = NS
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss - CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss, 
          CN.Diff.RNA.Prot_Reliable$RNA.Diff.Loss - CN.Diff.RNA.Prot_Reliable$Protein.Diff.Loss) 
#P = NS

## control vs high RNA express genes 
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain - CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain, 
       CN.Diff.RNA.Prot_highRNA$RNA.Diff.Gain - CN.Diff.RNA.Prot_highRNA$Protein.Diff.Gain) 
#P = NS
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss - CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss, 
       CN.Diff.RNA.Prot_highRNA$RNA.Diff.Loss - CN.Diff.RNA.Prot_highRNA$Protein.Diff.Loss) 
#P = NS

## control vs NO FLIPPED genes 
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain - CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain, 
       CN.Diff.RNA.Prot_noFlip$RNA.Diff.Gain - CN.Diff.RNA.Prot_noFlip$Protein.Diff.Gain) 
#P = 0.047 --> NS
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss - CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss, 
       CN.Diff.RNA.Prot_noFlip$RNA.Diff.Loss - CN.Diff.RNA.Prot_noFlip$Protein.Diff.Loss) 
#P = 1.12E-08--> 1.3E-07

## control vs no mutated genes 
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain - CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain, 
       CN.Diff.RNA.Prot_noMut$RNA.Diff.Gain - CN.Diff.RNA.Prot_noMut$Protein.Diff.Gain) 
#P = NS
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss - CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss, 
       CN.Diff.RNA.Prot_noMut$RNA.Diff.Loss - CN.Diff.RNA.Prot_noMut$Protein.Diff.Loss) 
#P = NS

## control vs merge (no flip, no mut, no low rna, no low reproducibility) genes 
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Gain - CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Gain, 
       CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA.Diff.Gain - CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Protein.Diff.Gain) 
#P = NS
t.test(CN.Diff.xRNA.yProt.ThreeGroups$RNA.Diff.Loss - CN.Diff.xRNA.yProt.ThreeGroups$Protein.Diff.Loss, 
       CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA.Diff.Loss - CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$Protein.Diff.Loss) 
#P = NS


###### Compare coding and non-coding RNA buffering #####
## Reviewer asked us to look at non-coding RNA (pseudogenes), and see if there is a 
# difference in mean RNA regulation between coding and non-coding RNA. 

### Step 1: get non-coding RNA expression data
### If you use this data about non-coding RNA (ncRNA), cite: 
## Wright MW et al. Human genomics 2011 Jan;5(2)90-98 PMID 21296742
#  Naming 'junk': human non-protein coding RNA (ncRNA) gene nomenclature
setwd("/Users/user/Dropbox (Sheltzer Lab)/Sheltzer Lab Shared/Klaske/RNA_Protein_expression_Paper")
NonCodingRNA.names<-read_tsv("NoncodingRNA.WrightMWetal2011.txt")

### also get ebi CCLE RNA-seq data with from the expression atlas 
##  Barretina J, et al. (2012) the cancer cell line encyclopedia enables predictive modelling of anticancer drug sensitivity
##  expression values across all genes (TPM) 
setwd("/Volumes/Schukken_SSD/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/")
RNA_Expression_with_ncRNA<- read_tsv("E-MTAB-2770-query-results.tpms.tsv")
# 56443 genes 


##  step 2:   get list of cells in dataset
# get only cell lines with aneuploidy data

celllines_ncRNA_info<-colnames(RNA_Expression_with_ncRNA)
celllines_ncRNA_info<-celllines_ncRNA_info[-1]
celllines_ncRNA_info<-celllines_ncRNA_info[-1]
celllines_ncRNA_info<-sub(",.*", "", celllines_ncRNA_info)

setwd("/Volumes/Schukken_SSD/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/cleaned up code")
aneuploid<- read.csv("arm_calls_data_bendavid.csv", header=TRUE) #had depmap ID for all cell lines with known aneuploidy

setwd("/Volumes/Schukken_SSD/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison")
Cell_Line_Info <- read_delim(file="sample_info.csv", 
                             delim=",") # had cell line info including depmap id and cell name
Cell_Line_Info_ncRNA<-subset(Cell_Line_Info, cell_line_name %in% celllines_ncRNA_info)
# 951 cell lines with aneuploidy info in this dataset
x<- Cell_Line_Info_ncRNA$cell_line_name
x<- append(x, c("Gene Name", "Gene ID"))

RNA_Exp_anCell_ncRNA<-RNA_Expression_with_ncRNA[, sub(",.*", "", colnames(RNA_Expression_with_ncRNA)) %in%  x] 
length(RNA_Exp_anCell_ncRNA$`Gene ID`)
length(RNA_Exp_anCell_ncRNA)
# 56,443 genes
# 951 cell lines + gene name and ID


### Step 3: get location data for each non-coding gene
# Get ensembl gene ID info, chromosome data for each ensemble ID: 
setwd("/Volumes/Schukken_SSD/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison")
Ensemble_Gene_Info<- read_tsv("Gene_Ensembl_ID_location.txt")
RNA_Exp_anCell_ncRNA2<- subset(RNA_Exp_anCell_ncRNA, `Gene ID` %in% Ensemble_Gene_Info$`Gene stable ID`)
#56,118 genes remain
Ensemble_Gene_Info$arm<- substr(Ensemble_Gene_Info$`Karyotype band`,1,1) # take first character of band for arm
Ensemble_Gene_Info$chrmArm<-paste0(Ensemble_Gene_Info$`Chromosome/scaffold name`, Ensemble_Gene_Info$arm)


### Step 4: quantile normalize RNA expression data
### quantile normalization from https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/ 
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

## Normalize RPE1 RNA data
### NOTE THAT NORMALIZED DATA IS NOT LOG2 TRANSFORMED!!! 
rownames(RNA_Exp_anCell_ncRNA2)<-RNA_Exp_anCell_ncRNA2$`Gene ID`
RNA_Exp_anCell_ncRNA3<-RNA_Exp_anCell_ncRNA2
RNA_Exp_anCell_ncRNA3[is.na(RNA_Exp_anCell_ncRNA3)] <- 0

RNA_Exp_anCell_ncRNA3<- quantile_normalisation(as.data.frame(RNA_Exp_anCell_ncRNA3[,3:length(RNA_Exp_anCell_ncRNA3)]))
RNA_Exp_anCell_ncRNA3<- data.frame(RNA_Exp_anCell_ncRNA3)
RNA_Exp_anCell_ncRNA3$'Gene ID'<-rownames(RNA_Exp_anCell_ncRNA2)

## each cell line now has an average RNA expression of 17.65 


### Step 4: Find average fold change 
#### AND DO LOG2 TRANSFORMATION!!!! 

# First flip the dataset so rows are cell lines and collumns are genes
RNA_Exp_anCell_ncRNA4<-t(RNA_Exp_anCell_ncRNA3) 
RNA_Exp_anCell_ncRNA4<-as.data.frame(RNA_Exp_anCell_ncRNA4)
rownames(RNA_Exp_anCell_ncRNA4)<- colnames(RNA_Exp_anCell_ncRNA3)
colnames(RNA_Exp_anCell_ncRNA4)<- RNA_Exp_anCell_ncRNA3$`Gene ID`

# now make sure cell line names are in same format. in our dataset all "-" and spaces are "." 
RNA_Exp_anCell_ncRNA4$Cell_Line_name<- sub("\\.\\..*", "", colnames(RNA_Exp_anCell_ncRNA3) ) 
Cell_Line_Info_ncRNA$cell_line_dots<- sub("[(]", ".", Cell_Line_Info_ncRNA$cell_line_name)
Cell_Line_Info_ncRNA$cell_line_dots<- sub("[)]", ".", Cell_Line_Info_ncRNA$cell_line_dots)
Cell_Line_Info_ncRNA$cell_line_dots<- sub("-", ".", Cell_Line_Info_ncRNA$cell_line_dots)
Cell_Line_Info_ncRNA$cell_line_dots<- sub("-", ".", Cell_Line_Info_ncRNA$cell_line_dots)
Cell_Line_Info_ncRNA$cell_line_dots<- sub("-", ".", Cell_Line_Info_ncRNA$cell_line_dots)
Cell_Line_Info_ncRNA$cell_line_dots<- sub("-", ".", Cell_Line_Info_ncRNA$cell_line_dots)
Cell_Line_Info_ncRNA$cell_line_dots<- sub("[ ]", ".", Cell_Line_Info_ncRNA$cell_line_dots)
Cell_Line_Info_ncRNA$cell_line_dots<- sub("[ ]", ".", Cell_Line_Info_ncRNA$cell_line_dots)
Cell_Line_Info_ncRNA$cell_line_dots<- sub("[ ]", ".", Cell_Line_Info_ncRNA$cell_line_dots)
Cell_Line_Info_ncRNA$cell_line_dots<- sub("[ ]", ".", Cell_Line_Info_ncRNA$cell_line_dots)
Cell_Line_Info_ncRNA$cell_line_dots<- sub("[/]", ".", Cell_Line_Info_ncRNA$cell_line_dots)
Cell_Line_Info_ncRNA$cell_line_dots<- sub("[/]", ".", Cell_Line_Info_ncRNA$cell_line_dots)
Cell_Line_Info_ncRNA$cell_line_dots<- sub("[;]", ".", Cell_Line_Info_ncRNA$cell_line_dots)
Cell_Line_Info_ncRNA$cell_line_dots<- sub("\\.\\..*", "", Cell_Line_Info_ncRNA$cell_line_dots)
Cell_Line_Info_ncRNA$cell_line_dots<- sub("SK.N.BE.2.", "SK.N.BE.2", Cell_Line_Info_ncRNA$cell_line_dots)
Cell_Line_Info_ncRNA$cell_line_dots<- ifelse(substring(Cell_Line_Info_ncRNA$cell_line_dots,1,1) %in% c(0:9), paste0("X", Cell_Line_Info_ncRNA$cell_line_dots), Cell_Line_Info_ncRNA$cell_line_dots)

# merge rna epression data with cell line info (including depmap I)
RNA_Exp_anCell_ncRNA5<- merge(RNA_Exp_anCell_ncRNA4, Cell_Line_Info_ncRNA[,c(1, 27)], 
                              by.x="Cell_Line_name", by.y="cell_line_dots") #add cell line depmap info
length(RNA_Exp_anCell_ncRNA5$Cell_Line_name)

#  log2 transform the normalized data: *****
RNA_Exp_anCell_ncRNA6<- RNA_Exp_anCell_ncRNA5 
for (i in 2:length(RNA_Exp_anCell_ncRNA6)){
  RNA_Exp_anCell_ncRNA6[,i]<- log(as.numeric(as.character(RNA_Exp_anCell_ncRNA5[,i])),2)
}

### now replace -Inf and Inf with NaN
RNA_Exp_anCell_ncRNA6 <- do.call(data.frame, # Replace Inf in data by NA
                   lapply(RNA_Exp_anCell_ncRNA6,
                          function(x) replace(x, is.infinite(x), NA)))
RNA_Exp_anCell_ncRNA6$DepMap_ID<- RNA_Exp_anCell_ncRNA5$DepMap_ID

### write this log2 transformed normalized RNA expression data
setwd("Users/user/Documents")
write.csv(RNA_Exp_anCell_ncRNA6, "E-MTAB-2770-RNA.ncRNA.log2.norm.csv")

RNA_Exp_anCell_ncRNA6<- read_csv("E-MTAB-2770-RNA.ncRNA.log2.norm.csv")

RNA.Diff.ncRNA.log<- data.frame(Ensembl_ID=character(),
                            ChrmArm=character(),
                             Diff.Gain=as.numeric(), 
                             Pvalue.Gain=as.numeric(), 
                             Diff.Loss=as.numeric(), 
                             Pvalue.Loss=as.numeric() ) 

for (i in 1:length(RNA_Exp_anCell_ncRNA3$'Gene ID')){
  
  Ensembl.ID=RNA_Exp_anCell_ncRNA3$'Gene ID'[i]
  
  testRNAChrm <- subset(Ensemble_Gene_Info, `Gene stable ID`==Ensembl.ID) #get test RNA data
  
  if (testRNAChrm$`Chromosome/scaffold name`[1] %in% c(1:22)){ #Make sure RNA data is in RNA_info3, else it crashes
    #filter cell data and get only aneuploidy scores for chrm arm of test RNA location
    testchrm.percell <- filter(aneuploid, chrom==testRNAChrm$`Chromosome/scaffold name`[1]) 
    testchrm.percell <- filter(testchrm.percell, arm==testRNAChrm$arm[1])
    
    
    RNA_Expression_TestRNA <- RNA_Exp_anCell_ncRNA6 %>% select(DepMap_ID, all_of(Ensembl.ID))
    
    
    #Combine RNA data with aneuploidy data for test RNA location
    RNA.effect.Chrm<-merge(y= RNA_Expression_TestRNA, x= testchrm.percell, 
                           by.y="DepMap_ID", by.x="DepMap_ID", 
                           sort = TRUE)# 368 cells, 12757 RNAs 
    RNA.effect.Chrm[,8]<-(RNA.effect.Chrm[,8]) #formatting
    
    ## put data into categories based on chrm arm number
    Chrm.tri<- filter(RNA.effect.Chrm, arm_call==1) #all cells with chrm gain for  testRNA chrm
    Chrm.di<- filter(RNA.effect.Chrm, arm_call==0) #all cells with chrm neutral ploidy for  testRNA chrm
    Chrm.mono<- filter(RNA.effect.Chrm, arm_call==-1) #all cells with chrm loss for  testRNA chrm
    
    ##Make data frame with t-test info about RNA expression per arm number category.
    ## Return this data at the end of the function. 
    if (colSums(!is.na(Chrm.tri[8]))>=10 &  #minimum 10 cells in "tri" category
        #Added if statement so only genes with 2+ values per condition are analyzed
        #if I don't do this, the t-test crashes and I get no values. 
        colSums(!is.na(Chrm.di[8]))>=10 &
        colSums(!is.na(Chrm.mono[8]))>=10) {
      # Trisomy vs disomy
      di.tri<-t.test(Chrm.tri[,8], # [,8] because that is column with RNA data
                     Chrm.di[,8], # [,8] because that is column with RNA data
                     mu = 0, 
                     alt = "two.sided",
                     conf.level = 0.99) #get p-value
      Diff.Tri.Di<-mean(Chrm.tri[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference
      #Disomy vs monosome
      Mono.Di<-t.test(Chrm.mono[,8], # [,8] because that is column with RNA data
                      Chrm.di[,8], # [,8] because that is column with RNA data
                      mu = 0, 
                      alt = "two.sided",
                      conf.level = 0.99) #get p-value
      Diff.Mono.Di<- mean(Chrm.mono[,8], na.rm=TRUE) - mean(Chrm.di[,8], na.rm=TRUE) #get difference

      RNA.Diff.ncRNA.log<- rbind(RNA.Diff.ncRNA.log, data.frame(
        Ensembl_ID=Ensembl.ID, 
        ChrmArm=testRNAChrm$chrmArm[1], #can have multiple genes because multiple rna's and proteins per gene (same chrm data though, so only get top row)
        Pvalue.Gain= di.tri$p.value,
        Diff.Gain= Diff.Tri.Di, 
        Pvalue.Loss= Mono.Di$p.value,
        Diff.Loss= Diff.Mono.Di
        ))
    }
  }
}
length(RNA.Diff.ncRNA.log$Ensembl_ID) # 53344 genes had 10+ points/caegory (of total 56,118 genes)
# after log2 transform and eliminating -Inf values==> 37983 genes had diff values

### some tests could not be done becasue data was "essentially constant" added P-value == NA

setwd("/Users/user/Documents")
#write.csv(RNA.Diff.ncRNA, 
#          file =paste("RNA_ncRNA_diff_pvalue_all.csv", sep=','), 
#          row.names = TRUE)

setwd("/Volumes/Schukken_SSD/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison")
#write.csv(RNA.Diff.ncRNA.log, 
#          file =paste("RNA_ncRNA_diff_pvalue_all_log2.csv", sep=','), 
#          row.names = TRUE)


RNA.Diff.ncRNA.log<-read.csv("RNA_ncRNA_diff_pvalue_all_log2.csv", 
                                 dec=".", header = TRUE, sep=",")
RNA.Diff.ncRNA.log<- RNA.Diff.ncRNA.log[,-1]



### Step 5: Compare non-coding to coding RNA

### add buffering/scaling/AS categories 
RNA.Diff.ncRNA.ThreeGroups<- RNA.Diff.ncRNA.log
RNA.Diff.ncRNA.ThreeGroups$RNA.Gain.Category<- cut(RNA.Diff.ncRNA.ThreeGroups$Diff.Gain,
                                                    breaks=c(-Inf,-0.1,0.25,Inf),
                                                    include.lowest=TRUE,
                                                    labels=c("Anti-Scaling","Buffering","Scaling"))
RNA.Diff.ncRNA.ThreeGroups$RNA.Loss.Category<- cut(RNA.Diff.ncRNA.ThreeGroups$Diff.Loss,
                                                    breaks=c(-Inf,-0.25,0.1,Inf),
                                                    include.lowest=TRUE,
                                                    labels=c("Scaling", "Buffering","Anti-Scaling"))
RNA.Diff.ncRNA.ThreeGroups$RNA.Loss.Category<-factor(RNA.Diff.ncRNA.ThreeGroups$RNA.Loss.Category, 
                                                      levels=c("Anti-Scaling","Buffering","Scaling"))

#### seperate into coding and non-coding categories
RNA_diff_ncRNA_only<- subset(RNA.Diff.ncRNA.ThreeGroups, `Ensembl_ID` %in% NonCodingRNA.names$`Ensembl gene ID`)
# 4,760 non-coding RNA 

RNA_diff_NOncRNA<- subset(RNA.Diff.ncRNA.ThreeGroups, ! `Ensembl_ID` %in% NonCodingRNA.names$`Ensembl gene ID`)
# 33,223 coding genes with location data

## mean difference in expression upon chromosome gain or loss
mean(RNA_diff_ncRNA_only$Diff.Gain) #non-coding mean diff upon gain = 0.167
mean(RNA_diff_NOncRNA$Diff.Gain)  #coding mean diff upon gain = 0.188
mean(RNA_diff_ncRNA_only$Diff.Loss)  #non-coding mean diff upon loss = -0.295 
mean(RNA_diff_NOncRNA$Diff.Loss)  # coding mean diff upon loss = -0.2877

## Diff gain: 
t.test(RNA_diff_ncRNA_only$Diff.Gain, RNA_diff_NOncRNA$Diff.Gain)
# 0.003987  ncRNA more buffered

## Diff loss: 
t.test(RNA_diff_ncRNA_only$Diff.Loss, RNA_diff_NOncRNA$Diff.Loss)
# NS 


#### now also look at specific types of non-coding RNAs
RNA_diff_ncRNA_only2<- merge(RNA_diff_ncRNA_only, NonCodingRNA.names[,c(10,13)], 
                             by.x="Ensembl_ID", by.y="Ensembl gene ID")
RNA_diff_ncRNA_only2$`Group name`<- sub("MicroRNA.*", "MicroRNAs", RNA_diff_ncRNA_only2$`Group name`)
RNA_diff_ncRNA_only2$`Group name`<- sub("Small nucleolar RNA.*", "Small nucleolar RNA", RNA_diff_ncRNA_only2$`Group name`)
RNA_diff_ncRNA_only2$`Group name`<- sub("Long non-coding.*", "Long non-coding RNAs", RNA_diff_ncRNA_only2$`Group name`)
RNA_diff_ncRNA_only2$`Group name`<- factor(RNA_diff_ncRNA_only2$`Group name`, 
                                             level=c("Antisense RNAs", "Divergent transcripts", "Intronic transcripts", 
                                                      "Long intergenic non-protein coding RNAs", "Long non-coding RNAs", 
                                                      "MicroRNAs", "Small nucleolar RNA",
                                                      "5S ribosomal RNAs", "Nuclear-encoded mitochondrial transfer RNAs",
                                                      "Overlapping transcripts", " RNAs, Ro60-associated Y", "Small Cajal body-specific RNAs", 
                                                      "Small nuclear RNAs", "Variant U1 small nuclear RNAs", "Vault RNAs"
                                                      ))
table(RNA_diff_ncRNA_only2$`Group name`)


RNA_diff_lincRNA<- subset(RNA_diff_ncRNA_only2, `Group name` == "Long intergenic non-protein coding RNAs")
mean(RNA_diff_lincRNA$Diff.Gain) #0.1392154
mean(RNA_diff_lincRNA$Diff.Loss) #-0.2491367
# 2171 Long intergenic non-protein coding RNAs


RNA_diff_microRNA<- subset(RNA_diff_ncRNA_only2, `Group name` == "MicroRNAs")
mean(RNA_diff_microRNA$Diff.Gain) #0.08508667 
mean(RNA_diff_microRNA$Diff.Loss) #-0.2191795
# 1394 MicroRNAs


RNA_diff_asRNA<- subset(RNA_diff_ncRNA_only2, `Group name` == "Antisense RNAs")
mean(RNA_diff_asRNA$Diff.Gain) #0.1949068 
mean(RNA_diff_asRNA$Diff.Loss) #-0.3533054
# 1698 Antisense RNAs


RNA_diff_dtRNA<- subset(RNA_diff_ncRNA_only2, `Group name` == "Divergent transcripts")
mean(RNA_diff_dtRNA$Diff.Gain) #0.2977186 
mean(RNA_diff_dtRNA$Diff.Loss) #-0.3556969
# 566 Divergent transcripts


RNA_diff_lncRNA<- subset(RNA_diff_ncRNA_only2, `Group name` == "Long non-coding RNAs")
mean(RNA_diff_lncRNA$Diff.Gain) #0.1405558 
mean(RNA_diff_lncRNA$Diff.Loss) #-0.2815256
# 387 Long non-coding RNAs

RNA_diff_snRNA<- subset(RNA_diff_ncRNA_only2, `Group name` == "Small nucleolar RNA")
mean(RNA_diff_snRNA$Diff.Gain) #0.1679094 
mean(RNA_diff_snRNA$Diff.Loss) #-0.3002559
# 318 Small nucleolar RNA

RNA_diff_itRNA<- subset(RNA_diff_ncRNA_only2, `Group name` == "Intronic transcripts")
mean(RNA_diff_itRNA$Diff.Gain) #0.1240303 
mean(RNA_diff_itRNA$Diff.Loss) #-0.3594183
# 122 Intronic transcripts


### get mean and SD per ncRNA type and plot as boxplot
## there's probably a cleaner way to do this, but eh... 


RNA.Diff.ncRNA.ThreeGroups2<- merge(RNA.Diff.ncRNA.ThreeGroups, NonCodingRNA.names[,c(10,13)], 
                             by.x="Ensembl_ID", by.y="Ensembl gene ID", all.x=TRUE)
RNA.Diff.ncRNA.ThreeGroups2$`Group name`[is.na(RNA.Diff.ncRNA.ThreeGroups2$`Group name`)]<-"Coding"
RNA.Diff.ncRNA.ThreeGroups2$`Group name`<- sub("MicroRNA.*", "MicroRNAs", RNA.Diff.ncRNA.ThreeGroups2$`Group name`)
RNA.Diff.ncRNA.ThreeGroups2$`Group name`<- sub("Small nucleolar RNA.*", "Small nucleolar RNA", RNA.Diff.ncRNA.ThreeGroups2$`Group name`)
RNA.Diff.ncRNA.ThreeGroups2$`Group name`<- sub("Long non-coding.*", "Long non-coding RNAs", RNA.Diff.ncRNA.ThreeGroups2$`Group name`)
RNA.Diff.ncRNA.ThreeGroups2$`Group name`<- factor(RNA.Diff.ncRNA.ThreeGroups2$`Group name`, 
                                           level=c("Coding","Antisense RNAs", "Divergent transcripts",
                                                   "Long intergenic non-protein coding RNAs", "Long non-coding RNAs", 
                                                   "MicroRNAs", "Small nucleolar RNA",
                                                   "Small nuclear RNAs"
                                           ))

ggplot(RNA.Diff.ncRNA.ThreeGroups2, aes(x=`Group name`, y=Diff.Gain, fill=`Group name`))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("grey50","Red2", "orange2", "Gold2", 
                             "green2", "cyan2", "blue", "purple2",
                             "grey10", "grey20","grey30", "grey40", 
                             "grey60", "grey70", "grey80"))+
  theme_classic()+
  ggtitle("Diff upon gain")+
  #stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(-1.5,1.5))+
  xlab("")+
  ylab("Difference upon gain")
# plot.box.ncRNA.mean.Gain_log2
# 8x5

ggplot(RNA.Diff.ncRNA.ThreeGroups2, aes(x=`Group name`, y=Diff.Loss, fill=`Group name`))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("grey50","Red2", "orange2", "Gold2", 
                             "green2", "cyan2", "blue", "purple2",
                             "grey10", "grey20","grey30", "grey40", 
                             "grey60", "grey70", "grey80"))+
  theme_classic()+
  ggtitle("Diff upon Loss")+
  #stat_summary(fun=mean, geom="point", shape=20, size=2, color="black", fill="black")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(-1.5,1.5))+
  xlab("")+
  ylab("Difference upon loss")
# plot.box.ncRNA.mean.Loss_log2
# 8x5


## t.test for significance upon GAIN
t.test(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Coding")$Diff.Gain, 
       subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Antisense RNAs")$Diff.Gain )
# Antisense RNAs  p-value = NS

t.test(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Coding")$Diff.Gain, 
       subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Divergent transcripts")$Diff.Gain )
# Divergent transcripts  p-value = 3.729e-10, divergent scaling --> ***

t.test(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Coding")$Diff.Gain, 
       subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Long intergenic non-protein coding RNAs")$Diff.Gain )
# Long intergenic non-protein coding RNAs  p-value= 0.003518, ncRNA buffered --> NS after correcting (0.0633)

t.test(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Coding")$Diff.Gain, 
       subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Long non-coding RNAs")$Diff.Gain )
# Long non-coding RNAs  p-value = NS, ncRNA buffered

t.test(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Coding")$Diff.Gain, 
       subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="MicroRNAs")$Diff.Gain )
# MicroRNAs  p-value = 2.3e-11, ncRNA buffered --> ***

t.test(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Coding")$Diff.Gain, 
       subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Small nucleolar RNA")$Diff.Gain )
# Small nucleolar RNA  p-value = NS 

t.test(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Coding")$Diff.Gain, 
       subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Small nuclear RNAs")$Diff.Gain )
# Small nuclear RNAs  p-value = NS  ncRNA buffered

t.test(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Coding")$Diff.Gain, 
       subset(RNA.Diff.ncRNA.ThreeGroups2, is.na(`Group name`))$Diff.Gain )
# Other RNAs  p-value = 0.0286, ncRNA buffered ("OTHER" RNAs)--> after correcting NS



## t.test for significance upon LOSS
t.test(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Coding")$Diff.Loss, 
       subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Antisense RNAs")$Diff.Loss )
# Antisense RNAs  p-value 5.466e-10, ncRNA scaling --> ***

t.test(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Coding")$Diff.Loss, 
       subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Divergent transcripts")$Diff.Loss )
# Divergent transcripts  p-value = 3.278e-07, ncRNA scaling --> ***

t.test(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Coding")$Diff.Loss, 
       subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Long intergenic non-protein coding RNAs")$Diff.Loss )
# Long intergenic non-protein coding RNAs  p-value = 0.006299, ncRNA buffered --> after correcting == 0.1133 NS

t.test(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Coding")$Diff.Loss, 
       subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Long non-coding RNAs")$Diff.Loss )
# Long non-coding RNAs  p-value = NS, ncRNA 

t.test(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Coding")$Diff.Loss, 
       subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="MicroRNAs")$Diff.Loss )
# MicroRNAs  p-value = 2.425e-08, ncRNA buffered --> ***

t.test(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Coding")$Diff.Loss, 
       subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Small nucleolar RNA")$Diff.Loss )
# Small nucleolar RNA  p-value = NS 

t.test(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Coding")$Diff.Loss, 
       subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Small nuclear RNAs")$Diff.Loss )
# Small nuclear RNAs  p-value = 0.03756 ncRNA buffered--> NS after correcting

t.test(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name`=="Coding")$Diff.Loss, 
       subset(RNA.Diff.ncRNA.ThreeGroups2, is.na(`Group name`))$Diff.Loss )
# Other RNAs  p-value = NS, ncRNA scale


## to correct for multiple tests, p-values were multipled by 18
## OR, significance < 2.78E-03  (*)
## ** = 2.8E-04
## *** = 2.8 E-05 or less


##Scatterplots with categories colored in: 
## scatterplot were not used in paper.
#ncRNA Chrm gain graphs
colnames(RNA_diff_ncRNA_only)
ggplot(RNA_diff_ncRNA_only, aes(x=Diff.Gain, y=-log2(Pvalue.Gain), color=RNA.Gain.Category))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in ncRNA expression")+
  ylab("log2(p-value)")+
  labs(color = "Category")+ #legend title
  theme_classic()+
  #coord_cartesian(xlim=c(-4,4), ylim=c(0,200))+
  ggtitle("ncRNA Gain scatterplot")
#plot.scatter.ncRNA.category.Gain_log2

#ncRNA Loss graphs
ggplot(RNA_diff_ncRNA_only, aes(x=Diff.Loss, y=-log2(Pvalue.Loss), color=RNA.Loss.Category))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in ncRNA expression")+
  ylab("log2(p-value)")+
  labs(color = "Category:")+ #legend title
  theme_classic()+
  #coord_cartesian(xlim=c(-4, 4), ylim=c(0,200))+
  ggtitle("ncRNA Loss scatterplot")


#coding RNA Chrm gain graphs
colnames(RNA_diff_NOncRNA)
ggplot(RNA_diff_NOncRNA, aes(x=Diff.Gain, y=-log2(Pvalue.Gain), color=RNA.Gain.Category))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in coding-RNA expression")+
  ylab("log2(p-value)")+
  labs(color = "Category")+ #legend title
  theme_classic()+
  #coord_cartesian(xlim=c(-4, 4), ylim=c(0,250))+
  ggtitle("coding-RNA Gain scatterplot")

#coding-RNA Loss graphs
ggplot(RNA_diff_NOncRNA, aes(x=Diff.Loss, y=-log2(Pvalue.Loss), color=RNA.Loss.Category))+
  geom_point(size=2)+
  scale_color_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Difference in coding-RNA expression")+
  ylab("log2(p-value)")+
  labs(color = "Category:")+ #legend title
  theme_classic()+
  #coord_cartesian(xlim=c(-4, 4), ylim=c(0,250))+
  ggtitle("coding-RNA Loss scatterplot")


#scatterplot based on RNA type
ggplot(RNA_diff_ncRNA_only2, aes(x=Diff.Gain, y=-log2(Pvalue.Gain), color=`Group name`))+ #color by rna type
  geom_point(size=2)+
  scale_color_manual(values=c("Red2", "orange2", "Gold2", 
                              "green2", "cyan2", "blue", "Purple2",
                              "grey10", "grey20","grey30", "grey40", 
                              "grey50", "grey60", "grey70", "grey80"))+
  xlab("Difference in ncRNA expression")+
  ylab("log2(p-value)")+
  labs(color = "Category:")+ #legend title
  theme_classic()+
  #coord_cartesian(xlim=c(-5, 5), ylim=c(0,100))+ #zoomed in
  #coord_cartesian(xlim=c(-100, 100), ylim=c(0,200))+  #not zoomed in
  ggtitle("ncRNA Gain scatterplot")
# zoom = xlim=c(-5, 5), ylim=c(0,50)
# No zoom = xlim=c(-100, 100), ylim=c(0,200)

#chrm loss (long intergenic non-protein coding RNAs AS & buffered most)
ggplot(RNA_diff_ncRNA_only2, aes(x=Diff.Loss, y=-log2(Pvalue.Loss), color=RNA_diff_ncRNA_only2$`Group name`))+ #color by rna type
  geom_point(size=2)+
  scale_color_manual(values=c("grey50","Red2", "orange2", "Gold2", 
                              "green2", "cyan2", "blue", "purple2",
                              "grey10", "grey20","grey30", "grey40", 
                              "grey50", "grey60", "grey70", "grey80"))+
  xlab("Difference in ncRNA expression")+
  ylab("log2(p-value)")+
  labs(color = "Category:")+ #legend title
  theme_classic()+
  #coord_cartesian(xlim=c(-5, 5), ylim=c(0,100))+ #zoomed in
  #coord_cartesian(xlim=c(-100, 100), ylim=c(0,200))+  #not zoomed in
  ggtitle("ncRNA Loss scatterplot")
# zoom = xlim=c(-5, 5), ylim=c(0,50)
# No zoom = xlim=c(-100, 100), ylim=c(0,200)



### Bar graph of RNA and Protein quantiles Scaling/Buffered

## coding RNA bargraph
dat.m <- melt(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name` == "Coding"), 
              id.vars='Ensembl_ID', 
              measure.vars=c('RNA.Gain.Category','RNA.Loss.Category'))
ggplot(data= dat.m, aes(x=variable, fill=value)) + 
  geom_bar(position="fill")+ #dodge (next to eachother)
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("coding RNA Difference upon aneuploidy: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  #coord_cartesian(ylim=c(0,6000))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("boxplot coding RNA category upon aneuploidy")
#4x4
# plot.codingRNA.Bargraph.ThreeCat_log_coding



## ncRNA bargraph: microRNAs
dat.m <- melt(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name` == "MicroRNAs"), 
              id.vars='Ensembl_ID', measure.vars=c('RNA.Gain.Category','RNA.Loss.Category'))
ggplot(data= dat.m, aes(x=variable, fill=value)) + 
  geom_bar(position="fill")+ #dodge (next to eachother)
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("MicroRNA Difference upon aneuploidy: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  #coord_cartesian(ylim=c(0,6000))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("boxplot MicroRNA category upon aneuploidy")
# 4x4
# plot.ncRNA.Bargraph.ThreeCat_log_MicroRNA


## ncRNA bargraph: Long intergenic non-protein coding RNAs
dat.m <- melt(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name` == "Long intergenic non-protein coding RNAs"), 
              id.vars='Ensembl_ID', measure.vars=c('RNA.Gain.Category','RNA.Loss.Category'))
ggplot(data= dat.m, aes(x=variable, fill=value)) + 
  geom_bar(position="fill")+ #dodge (next to eachother)
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("lincRNA Difference upon aneuploidy: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  #coord_cartesian(ylim=c(0,6000))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("boxplot lincRNA category upon aneuploidy")
# 4x4
# plot.ncRNA.Bargraph.ThreeCat_log_lincRNA

## ncRNA bargraph: Long  non-coding RNAs
dat.m <- melt(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name` == "Long non-coding RNAs"), 
              id.vars='Ensembl_ID', measure.vars=c('RNA.Gain.Category','RNA.Loss.Category'))
ggplot(data= dat.m, aes(x=variable, fill=value)) + 
  geom_bar(position="fill")+ #dodge (next to eachother)
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("lncRNA Difference upon aneuploidy: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  #coord_cartesian(ylim=c(0,6000))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("boxplot lncRNA category upon aneuploidy")
# 4x4
# plot.ncRNA.Bargraph.ThreeCat_log_lncRNA

## ncRNA bargraph: divergent transcripts RNAs
dat.m <- melt(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name` == "Divergent transcripts"), 
              id.vars='Ensembl_ID', measure.vars=c('RNA.Gain.Category','RNA.Loss.Category'))
ggplot(data= dat.m, aes(x=variable, fill=value)) + 
  geom_bar(position="fill")+ #dodge (next to eachother)
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("MicroRNA Difference upon aneuploidy: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  #coord_cartesian(ylim=c(0,6000))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("boxplot divergent transc. category upon aneuploidy")
# 4x4
# plot.ncRNA.Bargraph.ThreeCat_log_divergent

## ncRNA bargraph: sn RNAs
dat.m <- melt(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name` == "Small nuclear RNAs"), 
              id.vars='Ensembl_ID', measure.vars=c('RNA.Gain.Category','RNA.Loss.Category'))
ggplot(data= dat.m, aes(x=variable, fill=value)) + 
  geom_bar(position="fill")+ #dodge (next to eachother)
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("snRNA Difference upon aneuploidy: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  #coord_cartesian(ylim=c(0,6000))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("boxplot snRNA category upon aneuploidy")
# 4x4
# plot.ncRNA.Bargraph.ThreeCat_log_snRNA


## ncRNA bargraph: sn RNAs
dat.m <- melt(subset(RNA.Diff.ncRNA.ThreeGroups2, `Group name` == "Antisense RNAs"), 
              id.vars='Ensembl_ID', measure.vars=c('RNA.Gain.Category','RNA.Loss.Category'))
ggplot(data= dat.m, aes(x=variable, fill=value)) + 
  geom_bar(position="fill")+ #dodge (next to eachother)
  scale_fill_manual(values=c("salmon2", "Gold2", "palegreen3"))+
  xlab("Antisense RNA Difference upon aneuploidy: category")+
  ylab("Percent of genes")+
  labs(fill = "Protein Difference:\nCategory")+ #legend title
  theme_classic()+
  #coord_cartesian(ylim=c(0,6000))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("boxplot Antisense RNA category upon aneuploidy")
# 4x4
# plot.ncRNA.Bargraph.ThreeCat_log_antisenseRNA


  