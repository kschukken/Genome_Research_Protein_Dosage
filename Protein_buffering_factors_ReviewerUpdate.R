##### Factors affecting for protein buffering  ####
## Factors that influence gene attenuation at Protein/RNA level: 3 cat BUFFERING
## plotted as AUC ROC
## Both for CCLE data and ovarian tumor data (TCGA)

## look up ROC Area under the curve
## 210429: Buffering factors for 3-category distribution: anti-scale, scale and buffering
## 210330
## By Klaske M. Schukken


## Protein expression data: 
##  https://www.cell.com/cell/fulltext/S0092-8674(19)31385-6#secsectitle0190
## Nusinow et al. Cell, 2020, Quantitative proteomics of the cancer cell line encyclopedia
## Chromosome arm data: Uri Ben-David paper 
## RNA expression data:  CCLE depmap data.


# Protein expression for genes with minimum 10 cells per category (gain & loss) (RNA & Protein)
# and only cells with RNA & DNA expression data. 
# from Protein_RNA_Corr_min10.Filtered.R

library('ggplot2')
library('tidyverse')
library('readxl')
library('reshape2')
library('BBmisc')
library('dplyr')
library("ROCR")
library("pROC") #for auc() calculations, and roc() curves
library("readODS")
library("ggVennDiagram")

# Steps: 
# 1) Get filtered data (RNA & Protein & aneuploidy)
#       only cells that have aneuploudy, RNA and Protein data
#       only genes with minimum of 10 cells per category. 
#       get data about location of each gene
#       get data about difference upon gain or loss, per gene
# 2) Get factor data from various sources 
#       extract MobiBD Disorder score per gene, etc.
# 3) Merge factor data with difference per gene
#       for gain and loss, for RNA and Protein
#       get ROC AUC and correlations (?) 
# 4) Plot difference (gain and loss) vs. factor, and ROC AUC

##### Step 1: Get filtered RNA and Protein data. ####
## Step 1: Get filtered RNA and Protein data. 

#### Protein expression data 
# filtered to only get data from cells that have RNA and Protein expression
## Protein info: 12755 genes
setwd()
Protein_Info3<-read_csv2("Protein_location_info.csv")
Protein_ProID<-read_csv("Protein_ID_info.csv")
#Protein_ProID<-Protein_ProID[,-c(1)]# Remove "X" collumn

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


###  Aneuploidy data:
aneuploid<- read.csv("arm_calls_data_bendavid.csv", header=TRUE)


# Protein / RNA expression difference upon chrm arm gain/loss, per gene 
CN.Diff.xRNA.yProt.ThreeGroups <- read.csv("RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points_3cat.csv", header=TRUE) #9,414 genes
##!! or get this dataframe from supplementary data *! 

# Note: loss of about 2k genes (protein ID)-- 
# due to needing both RNA and Protein expression/gene
# and due to getting 10+ cell line with chrm gain/loss for chrm arm of gene location.

#add info about protein ID, uniprot ID, etc. 
CN.Diff.xRNA.yProt.ThreeGroups2<- merge(x=CN.Diff.xRNA.yProt.ThreeGroups, y=Protein_ProID, by.x="Protein_ID", by.y="Protein_Id")

# Make new category based on -0.25 to 0.25
CN.Diff.xRNA.yProt.ThreeGroups3<- CN.Diff.xRNA.yProt.ThreeGroups2

CN.Diff.xRNA.yProt.ThreeGroups3$Protein.Gain.Category2<-NA
for (i in 1:length(CN.Diff.xRNA.yProt.ThreeGroups3$Protein.Gain.Category2)){
if (CN.Diff.xRNA.yProt.ThreeGroups3$Protein.Diff.Gain[i]<  -0.25) {
  CN.Diff.xRNA.yProt.ThreeGroups3$Protein.Gain.Category2[i]<- "Anti-Scaling"
} else if (CN.Diff.xRNA.yProt.ThreeGroups3$Protein.Diff.Gain[i]>  0.25) {
  CN.Diff.xRNA.yProt.ThreeGroups3$Protein.Gain.Category2[i]<- "Scaling"
} else {
  CN.Diff.xRNA.yProt.ThreeGroups3$Protein.Gain.Category2[i]<- "Buffered"
}
}


CN.Diff.xRNA.yProt.ThreeGroups3$Protein.Loss.Category2<-NA
for (i in 1:length(CN.Diff.xRNA.yProt.ThreeGroups3$Protein.Loss.Category2)){
  if (CN.Diff.xRNA.yProt.ThreeGroups3$Protein.Diff.Loss[i]>  0.25) {
    CN.Diff.xRNA.yProt.ThreeGroups3$Protein.Loss.Category2[i]<- "Anti-Scaling"
  } else if (CN.Diff.xRNA.yProt.ThreeGroups3$Protein.Diff.Loss[i]<  -0.25) {
    CN.Diff.xRNA.yProt.ThreeGroups3$Protein.Loss.Category2[i]<- "Scaling"
  } else {
    CN.Diff.xRNA.yProt.ThreeGroups3$Protein.Loss.Category2[i]<- "Buffered"
  }
}

###Optional: 
#download the factor results & scores, so you don't need to calculcate them yourself
# If you upload these datasets to R, you can go down to STEP 4: Plot ROC bargraph, mean expression boxplot
All.Factors.Diff<- read.csv(file= "Protein.AllFactors.csv")
Corr.protein.ROCAUC<- read.csv(file="Factor.ROCAUC.correlationscore.csv")



##### Step 2: Get potential buffering factor info ####
## Step 2: Get potential buffering factor info
## Protein Intrinsic Disorder data from: MobiBD. Version: 4.0 - Release: 2020_09
#  “mobidb_result.tsv”
#  https://mobidb.bio.unipd.it
#get Disorder data from MobiBD

#setwd()#set working directory to where you downloaded MobiBD data
Intrinsic_Disorder_All<-read.delim2("mobidb_result.tsv", 
                                         dec=".", header = TRUE, sep="\t")

Intrinsic_Disorder_MobiBD<- filter(Intrinsic_Disorder_All, feature=="prediction-disorder-mobidb_lite")#NS difference
#Intrinsic_Disorder_2<- filter(Intrinsic_Disorder_All, feature=="curated-disorder-merge")# NS
#Intrinsic_Disorder_3<- filter(Intrinsic_Disorder_All, feature=="prediction-disorder-vsl")# NS
#Intrinsic_Disorder_4<- filter(Intrinsic_Disorder_All, feature=="prediction-disorder-glo")# NS (RNA loss correlates, but don't trust, low corr)
#Intrinsic_Disorder_5<- filter(Intrinsic_Disorder_All, feature=="prediction-disorder-espN")# NS
lip_anchor<- filter(Intrinsic_Disorder_All, feature=="prediction-lip-anchor")#NS difference
Intrinsic_Disorder_binding<- filter(Intrinsic_Disorder_All, feature=="derived-binding_mode_disorder_to_disorder-mobi")#NS difference
Intrinsic_Disorder_bindingdo<- filter(Intrinsic_Disorder_All, feature=="derived-binding_mode_disorder_to_order-mobi")# Yes, but don't trust it
#Intrinsic_Disorder_transmembrane<- filter(Intrinsic_Disorder_All, feature=="prediction-transmembrane-uniprot")#NS
#Intrinsic_Disorder_signalPep<- filter(Intrinsic_Disorder_All, feature=="prediction-signal_peptide-uniprot")#NS
Intrinsic_Disorder_CoilCoil<- filter(Intrinsic_Disorder_All, feature=="prediction-coiled_coil-uniprot")#
Intrinsic_Disorder_LowComplex<- filter(Intrinsic_Disorder_All, feature=="prediction-low_complexity-merge")#
Intrinsic_Disorder_Polarity<- filter(Intrinsic_Disorder_All, feature=="prediction-polar-mobidb_lite_sub")#
Intrinsic_Disorder_Polyampholyte<- filter(Intrinsic_Disorder_All, feature=="prediction-polyampholyte-mobidb_lite_sub")#
Intrinsic_Disorder_homology<- filter(Intrinsic_Disorder_All, feature=="homology-domain-merge")#technically RNA corrolates, but I don't buy it, very small corr


# Note: acc is Uniprot_Acc

####
### Get mRNA decay rates 
## data from: Yang et al. 2003, Genome Res. Decay Rates of Human mRNAs: Correlation With Functional Characteristics and Sequence Attributes
## Sup Table 9. 
#setwd()#set working directory to where you downloaded this data
mRNA_Decay<-read_xlsx("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/Yang.2003.mRNADecay.Rates.xlsx")
mRNA_Decay$Rate_1<- gsub(".\\s", "", mRNA_Decay$Rate_1)
mRNA_Decay$Rate_2<- gsub(".\\s", "", mRNA_Decay$Rate_2)
mRNA_Decay$Rate_3<- mRNA_Decay$Rate_1 #combine data from both rate columns
for (i in 1:length(mRNA_Decay$Accession)) {
  if (is.na(mRNA_Decay$Rate_3[i])){
    mRNA_Decay$Rate_3[i]<-mRNA_Decay$Rate_2[i] 
  }
}

headerRows<-c()# find rows with header titles instead of data
for (i in 1:length(mRNA_Decay$Accession)) {
  if (mRNA_Decay$Rate_3[i]=="Rate"){
    headerRows<-append(headerRows, c(i))
  }
}
mRNA_Decay2<-mRNA_Decay[-c(headerRows), ]#remove header rows
mRNA_Decay2<-mRNA_Decay2[,-c(3,4) ]#remove rate_1 and _2, already have _3
mRNA_Decay2$StdDev<-as.numeric(mRNA_Decay2$StdDev)
mRNA_Decay2$Rate_3<-as.numeric(mRNA_Decay2$Rate_3)
mRNA_Decay2$Accession<-as.factor(mRNA_Decay2$Accession)

#write.table(mRNA_Decay2$Accession, "mRNA.Decay_RefSeq.csv")

mRNADecay_info<- read_delim("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/uniprot_mRNA.Decay_info.tab", 
                                delim="\t") 

mRNA_Decay3<-merge(x=mRNA_Decay2, y=mRNADecay_info, 
                      by.x="Accession", by.y="yourlist:M202104065C475328CEF75220C360D524E9D456CE1A9041H") # genes

#write_csv(mRNA_Decay3, "mRNA_Decay_Mean_Yang.etal.csv")
#mRNA_Decay3<-read_csv(file="mRNA_Decay_Mean_Yang.etal.csv")



###
## NCBI genome browser data. Human RefSeq data. 
## NOTE:! multiple 5'UTR and/or 3'UTR per gene are possible. 
##      keep UTR data seperate from other factor lists, otherwise it skews data for other genes 
## Downloaded 210409
## 5'UTR and 3'UTR length data
## https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1081368649_vYvfytavAH0WudzGAUbcIpcRFbCA

#bin            Indexing field to speed chromosome range queries.
#name	          Name of gene (usually transcript_id from GTF)
#chrom	  	    Reference sequence chromosome or scaffold
#strand   	    + or - for strand
#txStart	    	Transcription start position (or end position for minus strand item)
#txEnd	        Transcription end position (or start position for minus strand item)
#cdsStart	      Coding region start (or end position for minus strand item)
#cdsEnd	        Coding region end (or start position for minus strand item)
#exonCount	    Number of exons
#exonStarts	    Exon start positions (or end positions for minus strand item)
#exonEnds	 	    Exon end positions (or start positions for minus strand item)
#score		      score
#name2		      Alternate name (e.g. gene_id from GTF)
#cdsStartStat	  Status of CDS start annotation (none, unknown, incomplete, or complete)
#cdsEndStat	    Status of CDS end annotation (none, unknown, incomplete, or complete)
#exonFrames		  Exon frame {0,1,2}, or -1 if no frame for exon

#172766 genes. from NCBI.
NCBI.genedata<- read_tsv("NCBI.Human.RefSeq.tsv")
#NCBI.genedata$UTR5<-NA
#NCBI.genedata$UTR3<-NA

for (i in 1:length(NCBI.genedata$name)){
  if (NCBI.genedata$strand[i]=="+"){ #if positive strand
    NCBI.genedata$UTR5[i]<-abs(NCBI.genedata$cdsStart[i]-NCBI.genedata$txStart[i])
    NCBI.genedata$UTR3[i]<-abs(NCBI.genedata$txEnd[i]-NCBI.genedata$cdsEnd[i])
  } else { #if negaitve strand
    NCBI.genedata$UTR5[i]<-abs(NCBI.genedata$txEnd[i]-NCBI.genedata$cdsEnd[i])
    NCBI.genedata$UTR3[i]<-abs(NCBI.genedata$cdsStart[i]-NCBI.genedata$txStart[i])
  }
}

#write_csv(NCBI.genedata, "NCBI.GeneData.3UTR.5UTR.csv")
NCBI.genedata<- read_csv("NCBI.GeneData.3UTR.5UTR.csv")

NCBI.genedata.NR<-subset(NCBI.genedata, startsWith(NCBI.genedata$name, "NR_")) #only munually curated, non coding genes
NCBI.genedata.NM<-subset(NCBI.genedata, startsWith(NCBI.genedata$name, "NM_")) #only munually curated, mRNA genes
NCBI.genedata.NM.5UTR<-unique(NCBI.genedata.NM[,c(13,17)]) #5' UTR data, unique entries, 33k
NCBI.genedata.NM.3UTR<-unique(NCBI.genedata.NM[,c(13,18)])#3'UTR data, unique entries, 26k




### Protein half-life data from: Methieson et al. 2018 Systematic analysis of protein turnover in primary cells
# 41467_2018_3106_MOESM5_ESM.xlsx
#setwd()#set working directory to where you downloaded this data
ProteinHL<-read_xlsx("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/41467_2018_3106_MOESM5_ESM.xlsx")
ProteinHL$Mean_HalfLife<-NA

for (i in 1:length(ProteinHL$Mean_HalfLife)) {
ProteinHL$Mean_HalfLife[i]<-mean(c(ProteinHL$`Bcells replicate 1 half_life`[i], 
                                 ProteinHL$`Bcells replicate 2 half_life`[i], 
                              ProteinHL$`NK cells replicate 1 half_life`[i], 
                              ProteinHL$`NK cells replicate 2 half_life`[i],
                              ProteinHL$`Hepatocytes replicate 1 half_life`[i], 
                              ProteinHL$`Hepatocytes replicate 2 half_life`[i], 
                              ProteinHL$`Monocytes replicate 1 half_life`[i], 
                              ProteinHL$`Monocytes replicate 2 half_life`[i], 
                              ProteinHL$`Mouse Neurons, replicate 3 half_life`[i], 
                              ProteinHL$`Mouse Neurons, replicate 4 half_life`[i]), 
                              na.rm=TRUE)
}

#write_csv(ProteinHL, "ProteinHalfLife.Methieson.csv")

ProteinHL<- read_csv("ProteinHalfLife.Methieson.csv")


#### Human transcription and translation rates, calculated. 
## from: Hausser, Mayo, Keren, &Alon. 2019. 
# Central dogma rates and the trade-off between precision and economy 
# in gene expression. Nature communications, 10 (68)
# trancription rate is equal to mRNA abundance* rate of mRNA decay, directly proportional.
# just do mRNA abundance.  mRNA abundance= transcription * decay

Central.Dogma.info<- read_ods("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/human.Transcription.Translation.Rates.ods") 
# variable                                                description                unit
# RNA.RPKM  log10 RPKM in RNAseq experiment of Eichhorn et al. (HeLa) reads / nucleotides
# RPF.RPKM      log10 RPKM in RP experiment of Eichhorn et al. (HeLa) reads / nucleotides
#        m       log10 mRNA abundance (estimated from Eichorn et al.)     copies per cell
#       bm  log10 transcription rate (estimated from Eichhorn et al.)              mRNA/h
#       bp    log10 translation rate (estimated from Eichhorn et al.)    protein/(mRNA.h)
#       lp                        log10 protein length (from UniProt)          aminoacids
#       lm log10 premRNA length (from RefSeq chromosomal coordinates)         nucleotides
#        p           log10 protein abundance estimated from bm and bp     proteins / cell

Central.Dogma<- read_ods("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/human.Transcription.Translation.Rates.ods", 
                         sheet=2) 
colnames(Central.Dogma)<- c("EnsembleID", "RNA.RPKM", "RPF.RPKM", "m", "bm", "bp", "lp", "lm", "p")
# added ensemble_ID column
# now to add Uniprot_accession ID, and/or gene symbol to the Ensemble ID data
#write.table(Central.Dogma$EnsembleID, "/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/CentralDogma_EnsembleID.csv")
#write.table(Central.Dogma$EnsembleID[5200:8449], "/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/CentralDogma_EnsembleID2.csv")

#Ran emsemble_ID's through BioDB.net https://biodbnet-abcc.ncifcrf.gov/db/db2dbRes.php 
# got Gene IDs and Gene symbols and uniprot IDs for each ensemble ID. 

CD_EnsembleID_info<- read_delim("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/bioDBnet_EnsembleID_Info.txt", 
                                delim="\t") #only 5200 ENSEMBLE ID's found back. from 8k genes in Central Dogma
CD_EnsembleID_info2<- read_delim("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/bioDBnet_EnsembleID_Info2.txt", 
                                delim="\t") #only 5200 ENSEMBLE ID's found back. from 8k genes in Central Dogma
CD_EnsembleID_info3<- rbind(CD_EnsembleID_info, CD_EnsembleID_info2)

Central.Dogma2<-merge(x=Central.Dogma, y=CD_EnsembleID_info3, 
                      by.x="EnsembleID", by.y="Ensembl Gene ID") #8449 genes

####
#Intrinsic Protein variance, and intrinsic RNA variance
# calculated in Variance_aneuploidy.R
# calculated from filtered Depmap data (cells with both RNA & Protein data only)
# from cells not aneuploid for chrm arm gene is located on
# from Variance_aneuploidy.R
#setwd()#set working directory to where you downloaded this data
Protein.var.sd<-read.delim2("Protein.Variance.SD.VarMean.csv", 
                            dec=".", header = TRUE, sep=",")
RNA.var.sd<-read.delim2("RNA.Variance.SD.VarMean.csv", 
                        dec=".", header = TRUE, sep=",")

Protein.var.sd$log2.CoeffVar<- -log2(Protein.var.sd$Variance.Mean.noAn)
RNA.var.sd$log2.CoeffVar<- -log2(RNA.var.sd$Variance.Mean.noAn)

####
#Phosphosite plus 
#https://www.phosphosite.org/staticDownloads
# downloaded April 27, 2021
#Datasets: Acetylation, methylation, phosphorylation, regulatory, ubiquitination, sumoylation
#setwd()#set working directory to where you downloaded this data
#Acetyl
Acetylation<-read.csv("Acetylation_site_dataset2.csv")
Acetylation2<-subset(Acetylation, Acetylation$ORGANISM=="human")
Acetylation3<-table(Acetylation2$ACC_ID) #9505 proteins, 14880 ACC_IDs

Acetylation4<-data.frame(ACC_ID=rownames(Acetylation3), 
                         Acetylation=Acetylation3) #14880 genes, 7629 of which are 0
#Methyl
Methylation<-read_tsv("Methylation_site_dataset2.tsv")
Methylation2<-subset(Methylation, Methylation$ORGANISM=="human")
Methylation3<-table(Methylation2$ACC_ID) #ACC_IDs

Methylation4<-data.frame(ACC_ID=rownames(Methylation3), 
                         Methylation=Methylation3)#5691 genes, none are 0. 

#Phosphorylation
Phosphorylation<-read_tsv("Phosphorylation_site_dataset2.tsv")
Phosphorylation2<-subset(Phosphorylation, Phosphorylation$ORGANISM=="human")
Phosphorylation3<-table(Phosphorylation2$ACC_ID) #ACC_IDs

Phosphorylation4<-data.frame(ACC_ID=rownames(Phosphorylation3), 
                            Phosphorylation=Phosphorylation3)#19833 genes, of which none 0

#ubiquitination
ubiquitination<-read_tsv("ubiquitination_site_dataset2.tsv")
ubiquitination2<-subset(ubiquitination, ubiquitination$ORGANISM=="human")
ubiquitination3<-table(ubiquitination2$ACC_ID) #ACC_IDs

ubiquitination4<-data.frame(ACC_ID=rownames(ubiquitination3), 
                            ubiquitination=ubiquitination3)#12435 genes, of which none 0

#Sumoylation
Sumoylation<-read_tsv("Sumoylation_site_dataset2.tsv")
Sumoylation2<-subset(Sumoylation, Sumoylation$ORGANISM=="human")
Sumoylation3<-table(Sumoylation2$ACC_ID) #ACC_IDs

Sumoylation4<-data.frame(ACC_ID=rownames(Sumoylation3), 
                         Sumoylation=Sumoylation3)#2669 genes, of which none 0

#regulatory
regulatory<-read_tsv("Regulatory_sites2.tsv")
regulatory2<-subset(regulatory, regulatory$ORGANISM=="human")
regulatory3<-table(regulatory2$ACC_ID) #ACC_IDs

regulatory4<-data.frame(ACC_ID=rownames(regulatory3), 
                        regulatory=regulatory3)#3110 genes, of which none 0

# did not end up using phosphosite plus data on galactose acytylation because not enough datapoints
# O-Gal-N-Ac

#O-GlcNAc


### Gene amplification percent ###
## From TSG.OG_CNV_Difference
## Percent of cells with gene amplifications(>2 DNA copy number)
## file from TSG.OG.CNV.Difference_v2.R
GeneAmplification<-read.csv(file="AmplificationRatio_perGene1.75.csv")

GeneAmplification$AmplificationRatio

### number of complexes gene is in: frequency in CORUM (all complexes dataset)
## Downloaded 210512
CORUM.Complex.all<-read_tsv(file="/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/CORUMallComplexes.tsv")

CORUM.Complex.all<-subset(CORUM.Complex.all, Organism=="Human")
CORUM.Complex.subunits<-as.vector(CORUM.Complex.all$`subunits(UniProt IDs)`)
CORUM.Complex.subunits<-strsplit(CORUM.Complex.subunits, ";")
CORUM.complex.subunits2 = c()
for (i in 1:length(CORUM.Complex.subunits)){
  CORUM.complex.subunits2<-append(CORUM.complex.subunits2, CORUM.Complex.subunits[[i]])
}
CORUM.complex.subunits3<-table(CORUM.complex.subunits2) #number of times genes occur in CORUM dataset
CORUM.complex.subunits4<-data.frame(ACC_ID=rownames(CORUM.complex.subunits3), 
                     gene=CORUM.complex.subunits3)#3110 genes, of which none 0


### number of Protein-Protein interactions
### HIPPIE version 2.2 updated 2/14/2009. downloaded 210513
## http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php
## Get interactions with >0.6 confidence (0=none, 1=complete confidence)
HIPPIE.PPI<-read_tsv(file="HIPPIE-current.mitab.tsv")
HIPPIE.PPI$"Confidence Value"<-as.numeric(HIPPIE.PPI$"Confidence Value")
HIPPIE.PPI.2<- subset(HIPPIE.PPI, HIPPIE.PPI$`Confidence Value`>=0.6 ) #only get PPI with >0.6 confidence

HIPPIE.PPI.3<-table(HIPPIE.PPI.2$`Gene Name Interactor A`) #number of times genes occur in CORUM dataset
HIPPIE.PPI.3.2<-data.frame(PPI=HIPPIE.PPI.3)#3110 genes, of which none 0
HIPPIE.PPI.4<-table(HIPPIE.PPI.2$`Gene Name Interactor B`) #number of times genes occur in CORUM dataset
HIPPIE.PPI.4.2<-data.frame(PPI=HIPPIE.PPI.4)#3110 genes, of which none 0
HIPPIE.PPI.all<-merge(HIPPIE.PPI.3.2, HIPPIE.PPI.4.2, 
                      by.x="PPI.Var1", by.y="PPI.Var1", all=TRUE)
HIPPIE.PPI.all$PPI.Freq.x[is.na(HIPPIE.PPI.all$PPI.Freq.x)]<-0
HIPPIE.PPI.all$PPI.Freq.y[is.na(HIPPIE.PPI.all$PPI.Freq.y)]<-0
HIPPIE.PPI.all$PPI.Freq<-(as.numeric(HIPPIE.PPI.all$PPI.Freq.x)+ as.numeric(HIPPIE.PPI.all$PPI.Freq.y))

###
#Non-exponential decay
#https://www.sciencedirect.com/science/article/pii/S009286741631248X
#Kinetic Analysis of Protein Stability Reveals Age-Dependent Degradation
#Science direct
NonExponentialDecay<-read_xlsx("Non-exponential decay.xlsx")
NonExponentialDecay2<-NonExponentialDecay[,-c(5,6)] # Δ-score (Non-exponential decay delta), higher= non exponential decay
NonExponentialDecay2$`Δ-score (Non-exponential decay delta)` #1569 genes

###
#Dependency score, depmap
Dependency<-read_xlsx("Dependency scores.xlsx")
Dependency2<-Dependency[,-c(3)] # Dependency Score
Dependency2$`Dependency.Score` #8964



###
# Aggregation score
# https://www.sciencedirect.com/science/article/pii/S2211124713005664
# Widespread Aggregation and Neurodegenerative Diseases Are Associated with Supersaturated Proteins

Aggregation<-read_xlsx("Aggregation score.xlsx")
Aggregation2<-Aggregation[,-c(4,5)] # Aggregation score
Aggregation2$`Aggregation score`

###
# Mutational score
# get all mutations in all cell lines
# subset only cells in our analysis
# cell line list from Protein_RNA_filtered_CellLine.R
# count table for gene mutations per gene
RNA_Protein_CellLines<- read_delim(file="RNA_Protein_Shared_cells.csv", 
                                   delim=";")

RNA_Protein_Cell_Info<-merge(y= RNA_Protein_CellLines, x= aneuploid, 
                             by.y="Cell_line", by.x="DepMap_ID", 
                             sort = TRUE, all=FALSE) 

mutation<-read_csv("CCLE_mutations.csv")
mutations1<-mutation[,c(1,2,16)] # 
mutations1<- subset(mutations1, DepMap_ID %in% RNA_Protein_Cell_Info$DepMap_ID)
mutations2<- table(mutations1$Hugo_Symbol)
mutations2<-data.frame(Hugo_Symbol=rownames(mutations2), 
                        Mutations=mutations2)#18586 genes

table(mutation$'Variant_Classification')
mutations.subset<- subset(mutation, Variant_Classification %in% c("Nonsense_Mutation", "Missense_Mutation"))
mutations.subset1<-mutations.subset[,c(1,2,16)] # 
mutations.subset1<- subset(mutations.subset1, DepMap_ID %in% RNA_Protein_Cell_Info$DepMap_ID)
mutations.subset2<- table(mutations.subset1$Hugo_Symbol)
mutations.nonsense.missense<-data.frame(Hugo_Symbol=rownames(mutations.subset2), 
                       mutations.nonsense.missense=mutations.subset2)#18586 genes

mutations.subset<- subset(mutation, Variant_Classification %in% c("Frame_Shift_Ins", "Frame_Shift_Del"))
mutations.subset3<-mutations.subset[,c(1,2,16)] # 
mutations.subset3<- subset(mutations.subset3, DepMap_ID %in% RNA_Protein_Cell_Info$DepMap_ID)
mutations.subset4<- table(mutations.subset3$Hugo_Symbol)
mutations.FrameShift<-data.frame(Hugo_Symbol=rownames(mutations.subset4), 
                                        mutations.FrameShift=mutations.subset4)#18586 genes

###
# Protein Reproducibility scores 
# Upadhya & Ryan BioRxiv 2021
# DOI:10.1101/2021.09.22.461108
# Corpus ID: 237637972
# Experimental reproducibility limits the correlation between mRNA and protein abundances in tumour proteomic profiles
# Supplementary table 2 
# 
Reproducibility<-read_xlsx("Upadhya.Ryan_BioRxiv_2021_S2_Protein_reliability.xlsx", sheet = "B. Protein reproducibility rank")
Reproducibility2<-Reproducibility[,c(1,5)] # aggregated reproducibility rank
colnames(Reproducibility2)<- c("Gene_ID", "Aggregated Reproducibility Rank")
Reproducibility2$`Aggregated Reproducibility Rank`

##### Step 3: Combine all factors into one dataset ####

# Step 1: Combine all factors into one dataset, add NA to non-info rows
#        Trim dataset: only keep useful collumns
#        Save dataset. 

#Factor.corr #9414 genes

colnames(Intrinsic_Disorder_MobiBD)<- c("acc", "feature", "start..end", "content_fraction", "Intrinsic_Disorder_MobiBD", "length")
colnames(lip_anchor)<- c("acc", "feature", "start..end", "content_fraction", "Loops.in.Protein(ANCHOR)", "length")
colnames(Intrinsic_Disorder_LowComplex)<- c("acc", "feature", "start..end", "content_fraction", "Intrinsic_Disorder_LowComplexity", "length")
colnames(Intrinsic_Disorder_Polarity)<- c("acc", "feature", "start..end", "content_fraction", "Intrinsic_Disorder_Polarity", "length")
colnames(Intrinsic_Disorder_Polyampholyte)<- c("acc", "feature", "start..end", "content_fraction", "Intrinsic_Disorder_Polyampholyte", "length")
colnames(Intrinsic_Disorder_homology)<- c("acc", "feature", "start..end", "content_fraction", "Intrinsic_Disorder_homology", "length")

#ProteinHL$Mean_HalfLife # 7125 data points
#Central.Dogma2$m # log10 abundance of mRNA, ~8k datapoints, 7029
#Central.Dogma2$p # log10 abundance of protein, ~8k datapoints
#Central.Dogma2$bm # log10 transcription rate, ~8k datapoints
#Central.Dogma2$bp # log10 translation rate, ~8k datapoints
#Central.Dogma2$lm # log10 pre-mRNA length, ~8k datapoints
#Central.Dogma2$lp # log10 protein length, ~8k datapoints
#mRNA_Decay3$Rate_3 # mRNA_Decay_rate, 2664 genes
#Protein.var.sd$log2.Protein.CoeffVar<-log2(Protein.var.sd$Variance.Mean.noAn) #intrinsic protein variance, 9256
#RNA.var.sd$log2.RNA.CoeffVar<-log2(RNA.var.sd$Variance.Mean.noAn) #instrinsic RNA variance, 9256
#GeneAmplification$AmplificationRatio #9401
#Acetylation4$Acetylation.Freq #4079, where not in database =0
#Methylation4$Methylation.Freq #3550, where not in database =0
#Phosphorylation4$Phosphorylation.Freq #8100, where not in database =0
#ubiquitination4$ubiquitination.Freq # 7035, where not in database =0
#Sumoylation4$Sumoylation.Freq # 1560, where not in database =0
#regulatory4$regulatory.Freq # 1754, where not in database =0
#CORUM.complex.subunits4$gene.Freq #2443, where not in database =0
#HIPPIE.PPI.all$PPI.Freq #8137 (confidence >0.6), where not in database =0
#Aggregation2$`Aggregation score`
#Dependency2$`Dependency.Score`  #8964
#NonExponentialDecay2$`Δ-score (Non-exponential decay delta)` #1569
#mutations2$Hugo_Symbol
#mutations.FrameShift$Hugo_Symbol
#mutations.nonsense.missense$Hugo_Symbol

### !! Do not add 5'UTR and 3'UTR, they have multiple datapoints per gene and can mess up counts, correlations, etc. 


## Make dataframe
All.Factors.Diff<-CN.Diff.xRNA.yProt.ThreeGroups2

## EDIT!!!:  add all factors one by one. merge via Gene_Symbol or Uniprot ID.
All.Factors.Diff<-merge(x=All.Factors.Diff, y= Reproducibility2, 
                        by.x="Gene_Symbol", by.y="Gene_ID", 
                        sort = TRUE, all.x = TRUE)

## some factors only include genes if they have X, so all other genes are set to 0. 
All.Factors.Diff$gene.Freq[is.na(All.Factors.Diff$gene.Freq)]<-0 #all genes not in in CORUM dataset ==0 freq
All.Factors.Diff$Mutations.Freq[is.na(All.Factors.Diff$Mutations.Freq)]<-0 # all genes with no mutations, have mutation==0
All.Factors.Diff$mutations.FrameShift.Freq[is.na(All.Factors.Diff$mutations.FrameShift.Var1)]<-0 
All.Factors.Diff$mutations.nonsense.missense.Freq[is.na(All.Factors.Diff$mutations.nonsense.missense.Freq)]<-0 
All.Factors.Diff$PPI.Freq[is.na(All.Factors.Diff$PPI.Freq)]<-0 
All.Factors.Diff$ubiquitination.Freq[is.na(All.Factors.Diff$ubiquitination.Freq)]<-0 
All.Factors.Diff$Acetylation.Freq[is.na(All.Factors.Diff$Acetylation.Freq)]<-0 
All.Factors.Diff$Methylation.Freq[is.na(All.Factors.Diff$Methylation.Freq)]<-0 
All.Factors.Diff$Phosphorylation.Freq[is.na(All.Factors.Diff$Phosphorylation.Freq)]<-0 
All.Factors.Diff$Sumoylation.Freq[is.na(All.Factors.Diff$Sumoylation.Freq)]<-0 
All.Factors.Diff$regulatory.Freq[is.na(All.Factors.Diff$regulatory.Freq)]<-0 

All.Factors.Diff$Aggregation.score<- as.numeric(as.character(All.Factors.Diff$Aggregation.score)) #correct format

colnames(CN.Diff.xRNA.yProt.ThreeGroups2)
# write.csv(All.Factors.Diff, file= "Protein.AllFactors_v2.csv")
## !! Or, you can download All factor values per gene from supplemental data: 
# All.Factors.Diff<- read.csv(file= "Protein.AllFactors_v2.csv")

All.Factors.Diff.noAS.Gain<-subset(All.Factors.Diff, Protein.Diff.Gain> -0.1) # Buffer vs scaling
All.Factors.Diff.noAS.Loss<-subset(All.Factors.Diff, Protein.Diff.Loss< 0.1) # Buffer vs scaling

All.Factors.Diff.noS.Gain<-subset(All.Factors.Diff, Protein.Diff.Gain< 0.25) # Buffer vs anti-scaling
All.Factors.Diff.noS.Loss<-subset(All.Factors.Diff, Protein.Diff.Loss> -0.25) # Buffer vs anti-scaling

#Protein_ID: sp.Q9NQ94.A1CF_HUMAN  sp.P01023.A2MG_HUMAN  sp.A8K2U0.A2ML1_HUMAN
#RNA_ID: A1CF..29974.   A2M..2.        A2ML1..144568.
#RNA_Name: A1CF  A2M   A2ML1
#Uniprot_Acc: "A0AV96"   "A0AVF1"   "A0AVT1"   "A0FGR8-2"


##### Step 4: RUN ROC analysis as FUNCTION: calculate ROC AUC ####

## Define function for getting ROC data per dataset, per category: 


# make list of column names and what they are titles: 
NameTestCollumns= c("PPI.Freq", "gene.Freq", "regulatory.Freq", "Sumoylation.Freq", 
                    "ubiquitination.Freq", "Phosphorylation.Freq", "Methylation.Freq", "Acetylation.Freq", 
                    "AmplificationRatio", "log2.CoeffVar.RNA", "log2.Protein.CoeffVar", "Rate_3", 
                    "m", "p", "lp", "lm", "bm", "bp", 
                    "Mean_HalfLife", "Loops.in.Protein.ANCHOR.", "Intrinsic_Disorder_MobiBD", "Intrinsic_Disorder_LowComplexity", 
                    "Intrinsic_Disorder_Polarity", "Intrinsic_Disorder_Polyampholyte", "Intrinsic_Disorder_homology", 
                    "Dependency.Score", "Aggregation.score", "Non.exponential.decay.delta", "Mutations.Freq", 
                    "mutations.FrameShift.Freq", "mutations.nonsense.missense.Freq", "Aggregated Reproducibility Rank")

testdatanames= c("Protein-protein interaction", "Protein complex (CORUM)", "Protein regulatory sites", "Sumoylation sites", 
                 "Ubiquitination sites", "Phosphorylation sites", "Methylation sites", "Acetylation sites", 
                 "Percent gene amplification", "RNA neutral variance", "Protein neutral variance", "mRNA decay rate", 
                 "mRNA abundance", "Protein abundance", "Protein length", "mRNA length", "Transcription rate", "Translation rate", 
                 "Protein half life", "Loops in protein score", "Intrinsic protein disorder", "Low complexity score", 
                 "Protein polarity", "Protein polyampholyte score", "Homology score", 
                 "Dependency score", "Aggregation score", "Non-exponential decay delta", "Mutation count (all)", 
                 "Mutation count (frame shifts)", "Mutation count (nonsence and missense)", "Reproducibility Rank")


dataset= All.Factors.Diff


FindROCAUC<- function(dataset, NameTestCollumns, testdatanames) {
  Corr.protein.loop<-data.frame(Data=as.character(), 
                                Chrm_CN=character(), 
                                Gene_Type=character(), 
                                Category=character(),
                                p_value=as.numeric(), 
                                Corr_coef=as.numeric(), 
                                significant=character(), 
                                ROCauc=as.numeric())
  
  
  for (i in 1:length(NameTestCollumns)){
    
    TestDataCol= NameTestCollumns[i]
    name= testdatanames[i]
    TestData<- as.numeric(as.character(dataset %>% pull(TestDataCol)))
    
    # Buffering, Gain: 
    Gain.Protein.corr<- cor.test(TestData, dataset$Protein.Diff.Gain, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Gain", 
                                        Gene_Type="Protein", 
                                        Category="Buffering",
                                        p_value=Gain.Protein.corr$p.value, 
                                        Corr_coef=Gain.Protein.corr$estimate, 
                                        significant=Gain.Protein.corr$p.value<0.05, 
                                        ROCauc=auc(dataset$Three.Protein.Gain=="Buffering", 
                                                   TestData, na.rm=TRUE ))) 
    
    # Buffering, Loss: 
    Loss.Protein.corr<-cor.test(TestData, dataset$Protein.Diff.Loss, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Loss", 
                                        Gene_Type="Protein", 
                                        Category="Buffering",
                                        p_value=Loss.Protein.corr$p.value, 
                                        Corr_coef=Loss.Protein.corr$estimate, 
                                        significant=Loss.Protein.corr$p.value<0.01, 
                                        ROCauc=auc(dataset$Three.Protein.Loss=="Buffering", 
                                                   TestData, na.rm=TRUE )))
    
    # AS, Gain: 
    Gain.Protein.corr<- cor.test(TestData, dataset$Protein.Diff.Gain, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Gain", 
                                        Gene_Type="Protein", 
                                        Category="Anti-Scaling",
                                        p_value=Gain.Protein.corr$p.value, 
                                        Corr_coef=Gain.Protein.corr$estimate, 
                                        significant=Gain.Protein.corr$p.value<0.05, 
                                        ROCauc=auc(dataset$Three.Protein.Gain=="Anti-Scaling", 
                                                   TestData, na.rm=TRUE ))) 
    
    # AS, Loss: 
    Loss.Protein.corr<-cor.test(TestData, dataset$Protein.Diff.Loss, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Loss", 
                                        Gene_Type="Protein", 
                                        Category="Anti-Scaling",
                                        p_value=Loss.Protein.corr$p.value, 
                                        Corr_coef=Loss.Protein.corr$estimate, 
                                        significant=Loss.Protein.corr$p.value<0.01, 
                                        ROCauc=auc(dataset$Three.Protein.Loss=="Anti-Scaling", 
                                                   TestData, na.rm=TRUE ))) 
    
    # Scaling, Gain: 
    Gain.Protein.corr<- cor.test(TestData, dataset$Protein.Diff.Gain, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Gain", 
                                        Gene_Type="Protein", 
                                        Category="Scaling",
                                        p_value=Gain.Protein.corr$p.value, 
                                        Corr_coef=Gain.Protein.corr$estimate, 
                                        significant=Gain.Protein.corr$p.value<0.05, 
                                        ROCauc=auc(dataset$Three.Protein.Gain=="Scaling", 
                                                   TestData, na.rm=TRUE ))) 
    
    # Scaling, Loss: 
    Loss.Protein.corr<-cor.test(TestData, dataset$Protein.Diff.Loss, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Loss", 
                                        Gene_Type="Protein", 
                                        Category="Scaling",
                                        p_value=Loss.Protein.corr$p.value, 
                                        Corr_coef=Loss.Protein.corr$estimate, 
                                        significant=Loss.Protein.corr$p.value<0.01, 
                                        ROCauc=auc(dataset$Three.Protein.Loss=="Scaling", 
                                                   TestData, na.rm=TRUE ))) 
    
  }
  return(Corr.protein.loop)
}
Corr.protein.ROCAUC<- FindROCAUC(dataset, NameTestCollumns, testdatanames)


# Below: Do same ROC analysis for +/-0.25 cutoffs, instead of -0.1 to 0.25 cutoffs.
# This is for supplementat data. 
# Replace "Corr.protein.ROCAUC" with "Corr.protein.ROCAUC_0.25AS" in code below 
# to plot this alternative cutoff analysis
 FindROCAUC_0.25AS<- function(dataset, NameTestCollumns, testdatanames) {
  Corr.protein.loop<-data.frame(Data=as.character(), 
                                Chrm_CN=character(), 
                                Gene_Type=character(), 
                                Category=character(),
                                p_value=as.numeric(), 
                                Corr_coef=as.numeric(), 
                                significant=character(), 
                                ROCauc=as.numeric())
  
  
  for (i in 1:length(NameTestCollumns)){
    
    TestDataCol= NameTestCollumns[i]
    name= testdatanames[i]
    TestData<- as.numeric(as.character(dataset %>% pull(TestDataCol)))
    
    # Buffering, Gain: 
    Gain.Protein.corr<- cor.test(TestData, dataset$Protein.Diff.Gain, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Gain", 
                                        Gene_Type="Protein", 
                                        Category="Buffering",
                                        p_value=Gain.Protein.corr$p.value, 
                                        Corr_coef=Gain.Protein.corr$estimate, 
                                        significant=Gain.Protein.corr$p.value<0.05, 
                                        ROCauc=auc(dataset$Protein.Diff.Gain>-0.25 & dataset$Protein.Diff.Gain< 0.25, 
                                                   TestData, na.rm=TRUE ))) 
    
    # Buffering, Loss: 
    Loss.Protein.corr<-cor.test(TestData, dataset$Protein.Diff.Loss, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Loss", 
                                        Gene_Type="Protein", 
                                        Category="Buffering",
                                        p_value=Loss.Protein.corr$p.value, 
                                        Corr_coef=Loss.Protein.corr$estimate, 
                                        significant=Loss.Protein.corr$p.value<0.05, 
                                        ROCauc=auc(dataset$Protein.Diff.Loss>-0.25 & dataset$Protein.Diff.Loss< 0.25, 
                                                   TestData, na.rm=TRUE )))
    
    # AS, Gain: 
    Gain.Protein.corr<- cor.test(TestData, dataset$Protein.Diff.Gain, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Gain", 
                                        Gene_Type="Protein", 
                                        Category="Anti-Scaling",
                                        p_value=Gain.Protein.corr$p.value, 
                                        Corr_coef=Gain.Protein.corr$estimate, 
                                        significant=Gain.Protein.corr$p.value<0.05, 
                                        ROCauc=auc(dataset$Protein.Diff.Gain<  -0.25 , 
                                                   TestData, na.rm=TRUE ))) 
    
    # AS, Loss: 
    Loss.Protein.corr<-cor.test(TestData, dataset$Protein.Diff.Loss, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Loss", 
                                        Gene_Type="Protein", 
                                        Category="Anti-Scaling",
                                        p_value=Loss.Protein.corr$p.value, 
                                        Corr_coef=Loss.Protein.corr$estimate, 
                                        significant=Loss.Protein.corr$p.value<0.05, 
                                        ROCauc=auc(dataset$Protein.Diff.Loss> 0.25, 
                                                   TestData, na.rm=TRUE ))) 
    
    # Scaling, Gain: 
    Gain.Protein.corr<- cor.test(TestData, dataset$Protein.Diff.Gain, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Gain", 
                                        Gene_Type="Protein", 
                                        Category="Scaling",
                                        p_value=Gain.Protein.corr$p.value, 
                                        Corr_coef=Gain.Protein.corr$estimate, 
                                        significant=Gain.Protein.corr$p.value<0.05, 
                                        ROCauc=auc(dataset$Three.Protein.Gain=="Scaling", 
                                                   TestData, na.rm=TRUE ))) 
    
    # Scaling, Loss: 
    Loss.Protein.corr<-cor.test(TestData, dataset$Protein.Diff.Loss, method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Loss", 
                                        Gene_Type="Protein", 
                                        Category="Scaling",
                                        p_value=Loss.Protein.corr$p.value, 
                                        Corr_coef=Loss.Protein.corr$estimate, 
                                        significant=Loss.Protein.corr$p.value<0.05, 
                                        ROCauc=auc(dataset$Three.Protein.Loss=="Scaling", 
                                                   TestData, na.rm=TRUE ))) 
    
  }
  return(Corr.protein.loop)
}
# Corr.protein.ROCAUC_0.25AS<- FindROCAUC_0.25AS(dataset, NameTestCollumns, testdatanames)

 FindCorAbs<- function(dataset, NameTestCollumns, testdatanames) {
  Corr.protein.loop<-data.frame(Data=as.character(), 
                                Chrm_CN=character(), 
                                p_value=as.numeric(), 
                                Corr_coef=as.numeric(), 
                                significant=character())
  
  
  for (i in 1:length(NameTestCollumns)){
    
    TestDataCol= NameTestCollumns[i]
    name= testdatanames[i]
    TestData<- as.numeric(as.character(dataset %>% pull(TestDataCol)))
    
    # Gain: 
    Gain.Protein.corr.2<- cor.test(TestData, abs(dataset$Protein.Diff.Gain), method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Gain", 
                                        Gene_Type="Protein",
                                        p_value=Gain.Protein.corr.2$p.value, 
                                        Corr_coef=Gain.Protein.corr.2$estimate, 
                                        significant=Gain.Protein.corr.2$p.value*62<0.05)) 
    
    # Loss: 
    Loss.Protein.corr.2<-cor.test(TestData, abs(dataset$Protein.Diff.Loss), method="pearson")
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Loss", 
                                        Gene_Type="Protein", 
                                        p_value=Loss.Protein.corr.2$p.value, 
                                        Corr_coef=Loss.Protein.corr.2$estimate, 
                                        significant=Loss.Protein.corr.2$p.value*62<0.05))
    
  }
  return(Corr.protein.loop)
}
# Corr.protein.abs<- FindCorAbs(dataset, NameTestCollumns, testdatanames)
# Pearson correlation between absolute protein difference upon aneuploidy and factor values
# Significance is true if p< 0.05 after correcting for multiple tests (62 tests)
# p-values are pearson p value (uncorrected, multiply by 62 to correct)



## Now look at mean factor score per category: 
FindMeanFactorPerCategory<- function(dataset, NameTestCollumns, testdatanames) {
  Mean.Factor<-data.frame(Factor=as.character(),
                          Gain.AntiScaling= as.numeric(), 
                          Gain.Buffering= as.numeric(),  
                          Gain.Scaling= as.numeric(), 
                          Loss.AntiScaling= as.numeric(), 
                          Loss.Buffering= as.numeric(), 
                          Loss.Scaling= as.numeric() )
  
  for (i in 1:length(NameTestCollumns)){
    TestColName= NameTestCollumns[i]
    TestColData= dataset %>% pull(TestColName)
    name= testdatanames[i]
    # get mean expression of variable per group, 
    # normalize to 0-1 scale: (value-minimum)/max-minimum
    
    table.Gain <- data.frame( tapply(TestColData, dataset$Three.Protein.Gain, mean, na.rm=TRUE) )
    table.Gain$StandardData<- table.Gain[,1] 
    
    table.Loss <- data.frame( tapply(TestColData, dataset$Three.Protein.Loss, mean, na.rm=TRUE) )
    table.Loss$StandardData<- table.Loss[,1]
    
    Mean.Factor<-rbind(Mean.Factor, 
                       data.frame(Factor=name, 
                                  Gain.AntiScaling= table.Gain[1,2], 
                                  Gain.Buffering= table.Gain[2,2], 
                                  Gain.Scaling= table.Gain[3,2], 
                                  Loss.AntiScaling= table.Loss[1,2], 
                                  Loss.Buffering= table.Loss[2,2], 
                                  Loss.Scaling= table.Loss[3,2] ) )
  }
  return(Mean.Factor)
} # not standardized, just give values
FindMeanFactorPerCategory2<- function(dataset, NameTestCollumns, testdatanames) {
  Mean.Factor2<-data.frame(Factor=character(), 
                           Gain.AntiScaling= as.numeric(), 
                           Gain.Buffering= as.numeric(),  
                           Gain.Scaling= as.numeric(), 
                           Loss.AntiScaling= as.numeric(), 
                           Loss.Buffering= as.numeric(), 
                           Loss.Scaling= as.numeric() )
  
  for (i in 1:length(NameTestCollumns)){
    TestColName= NameTestCollumns[i]
    TestColData= dataset %>% pull(TestColName)
    name= testdatanames[i]
    # get mean expression of variable per group, 
    # normalize to 0-1 scale: (value-minimum)/max-minimum
    
    table.Gain <- data.frame( tapply(TestColData, dataset$Three.Protein.Gain, mean, na.rm=TRUE) )
    table.Gain$StandardData<- (table.Gain[,1] - mean(TestColData, na.rm=TRUE) )/ sd(TestColData, na.rm=TRUE)
    
    table.Loss <- data.frame( tapply(TestColData, dataset$Three.Protein.Loss, mean, na.rm=TRUE) )
    table.Loss$StandardData<- (table.Loss[,1] - mean(TestColData, na.rm=TRUE) )/ sd(TestColData, na.rm=TRUE)
    
    Mean.Factor2<-rbind(Mean.Factor2, 
                        data.frame(Factor=name, 
                                   Gain.AntiScaling= table.Gain[1,2], 
                                   Gain.Buffering= table.Gain[2,2], 
                                   Gain.Scaling= table.Gain[3,2], 
                                   Loss.AntiScaling= table.Loss[1,2], 
                                   Loss.Buffering= table.Loss[2,2], 
                                   Loss.Scaling= table.Loss[3,2] ) )
  }
  return(Mean.Factor2)
} # Standardized: x minus mean, divided by standard dev. 
##
Mean.Factor.PerCategory.notNorm<- FindMeanFactorPerCategory(dataset, NameTestCollumns, testdatanames)
#Mean.Factor.PerCategory<- FindMeanFactorPerCategory2(dataset, NameTestCollumns, testdatanames)


#setwd()#set working directory to where you downloaded this data
# write.csv(Mean.Factor.PerCategory.notNorm, file="Mean.Factor.PerCategory.notStandardized")
#write.csv(Mean.Factor.PerCategory, file="Mean.Factor.PerCategory.Standardized")

# Mean.Factor.PerCategory.notNorm<- read.csv("Mean.Factor.PerCategory.notStandardized")

#####         Plot ROC bargraph, mean expression boxplot ####

### boxplots of mean factor score/value in Buffered, scaling and antiscaling genes upon either gain and loss

y.lab="Ubiquitination" #Name of factor you are looking at
testcollumn<- All.Factors.Diff$ubiquitination.Freq #get the collumn of factor you want to look at


## Chrm gain categories
ggplot(All.Factors.Diff, aes(x=Three.Protein.Gain, y=testcollumn, fill=Three.Protein.Gain))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("salmon", "Gold2", "palegreen3"))+
  theme_classic()+
  ggtitle("Gain")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(-0.2,30))+
  #stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="black", fill="black")+
  xlab("")+
  ylab(y.lab)
# 5x4
# plot.mean.factor.Gain.boxplot.Ubiquitin

#Chrm loss categories
ggplot(All.Factors.Diff, aes(x=Three.Protein.Loss, y=testcollumn, fill=Three.Protein.Loss))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("salmon", "Gold2", "palegreen3"))+
  theme_classic()+
  ggtitle("Loss")+
  geom_hline(yintercept=0)+
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="black", fill="black")+
  coord_cartesian(ylim=c(-0.2, 30))+
  xlab("")+
  ylab(y.lab)
# 5x4
# plot.mean.factor.Loss.boxplot.Ubiquitin_mean
# colors= salmon", "Gold2", "palegreen3
# colors2= "darkorange2", "Gold2", "green3"

G.AS<- subset(All.Factors.Diff, Three.Protein.Gain=="Anti-Scaling") #chrm gain, category anti scaling
G.B<- subset(All.Factors.Diff, Three.Protein.Gain=="Buffering") #chrm gain, category buffered
G.S<- subset(All.Factors.Diff, Three.Protein.Gain=="Scaling") #chrm gain, category scaling
L.AS<- subset(All.Factors.Diff, Three.Protein.Loss=="Anti-Scaling")
L.B<- subset(All.Factors.Diff, Three.Protein.Loss=="Buffering")
L.S<- subset(All.Factors.Diff, Three.Protein.Loss=="Scaling")

## Test significance between groups
# change the collumn to the factor you want to test. 
t.test(G.AS$Phosphorylation.Freq, G.B$Phosphorylation.Freq)
t.test(G.S$Phosphorylation.Freq,  G.B$Phosphorylation.Freq)
t.test(L.AS$Phosphorylation.Freq, L.B$Phosphorylation.Freq)
t.test(L.S$Phosphorylation.Freq,  L.B$Phosphorylation.Freq)



### Plot ROC bargraphs:  ###

## But first Add 5' and 3' UTR  AUC to data ***
NCBI.genedata.NM.5UTR$UTR5 #5' UTR data, unique entries, 33k
TestData<- NCBI.genedata.NM.5UTR
name<- "5' UTR length" 

Diff.test<-merge(x=TestData, y= CN.Diff.xRNA.yProt.ThreeGroups2, 
                 by.x="name2", by.y="Gene_Symbol", 
                 sort = TRUE)
TestDataCol<-Diff.test$UTR5

Gain.Protein.corr<- cor.test(TestDataCol, Diff.test$Protein.Diff.Gain, method="pearson")
Corr.protein.ROCAUC_0.25AS<-rbind(Corr.protein.ROCAUC_0.25AS, 
                                   data.frame(data=name, 
                                              Chrm_CN="Gain", 
                                              Gene_Type="Protein", 
                                              Category= "Buffering",
                                              p_value=Gain.Protein.corr$p.value, 
                                              Corr_coef=Gain.Protein.corr$estimate, 
                                              significant=Gain.Protein.corr$p.value<0.05, 
                                              ROCauc=auc(Diff.test$Three.Protein.Gain=="Buffering", 
                                                         TestDataCol ))) 
Loss.Protein.corr<- cor.test(TestDataCol, Diff.test$Protein.Diff.Loss, method="pearson")
Corr.protein.ROCAUC_0.25AS<-rbind(Corr.protein.ROCAUC_0.25AS, 
                                   data.frame(data=name, 
                                              Chrm_CN="Loss", 
                                              Gene_Type="Protein", 
                                              Category= "Buffering",
                                              p_value=Loss.Protein.corr$p.value, 
                                              Corr_coef=Loss.Protein.corr$estimate, 
                                              significant=Loss.Protein.corr$p.value<0.01, 
                                              ROCauc=auc(Diff.test$Three.Protein.Loss=="Buffering", 
                                                         TestDataCol ))) 


## Now add 3' UTR 
NCBI.genedata.NM.3UTR$UTR3#3'UTR data, unique entries, 26k
TestData<- NCBI.genedata.NM.3UTR
name<- "3' UTR length" 

Diff.test<-merge(x=TestData, y= CN.Diff.xRNA.yProt.ThreeGroups2, 
                 by.x="name2", by.y="Gene_Symbol", 
                 sort = TRUE)
TestDataCol<-Diff.test$UTR3

Gain.Protein.corr<- cor.test(TestDataCol, Diff.test$Protein.Diff.Gain, method="pearson")
Corr.protein.ROCAUC_0.25AS<-rbind(Corr.protein.ROCAUC_0.25AS, 
                                   data.frame(data=name, 
                                              Chrm_CN="Gain", 
                                              Gene_Type="Protein", 
                                              Category= "Buffering",
                                              p_value=Gain.Protein.corr$p.value, 
                                              Corr_coef=Gain.Protein.corr$estimate, 
                                              significant=Gain.Protein.corr$p.value<0.05, 
                                              ROCauc=auc(Diff.test$Three.Protein.Gain=="Buffering", 
                                                         TestDataCol ))) 
Loss.Protein.corr<- cor.test(TestDataCol, Diff.test$Protein.Diff.Loss, method="pearson")
Corr.protein.ROCAUC_0.25AS<-rbind(Corr.protein.ROCAUC_0.25AS, 
                                   data.frame(data=name, 
                                              Chrm_CN="Loss", 
                                              Gene_Type="Protein", 
                                              Category= "Buffering",
                                              p_value=Loss.Protein.corr$p.value, 
                                              Corr_coef=Loss.Protein.corr$estimate, 
                                              significant=Loss.Protein.corr$p.value<0.01, 
                                              ROCauc=auc(Diff.test$Three.Protein.Loss=="Buffering", 
                                                         TestDataCol ))) 

write.csv(Corr.protein.ROCAUC_0.25AS, file="Factor.ROCAUC.correlationscore_0.25AS.csv")
## ! OR you can download the file from supplementary: 
# Corr.protein.ROCAUC<- read.csv(file="Factor.ROCAUC.correlationscore.csv")

## now that 5' and 3 UTR have been added, remove dataset specific factors: 
Corr.protein.noHighRNAVar2<-Corr.protein.ROCAUC_0.25AS
Corr.protein.noHighRNAVar2<-Corr.protein.noHighRNAVar2[-c(49:66, 151:156, 169:186),] #remove RNA/Protein neutral variance, & gene amp ratio


#Subset data by type: 
Corr.protein.noHighRNAVar2.Gain<- subset(Corr.protein.noHighRNAVar2, Chrm_CN=="Gain" & Gene_Type=="Protein" & Category=="Buffering")
Corr.protein.noHighRNAVar2.Gain<- Corr.protein.noHighRNAVar2.Gain[order(Corr.protein.noHighRNAVar2.Gain$ROCauc),]

Corr.protein.noHighRNAVar2.Loss<- subset(Corr.protein.noHighRNAVar2, Chrm_CN=="Loss" & Gene_Type=="Protein" & Category=="Buffering")
Corr.protein.noHighRNAVar2.Loss<- Corr.protein.noHighRNAVar2.Loss[order(Corr.protein.noHighRNAVar2.Loss$ROCauc),]


## Buffering ROC AUC bargraph
Corr.protein.noHighRNAVar2.Gain$data<- factor(Corr.protein.noHighRNAVar2.Gain$data, level=Corr.protein.noHighRNAVar2.Gain$data)
Corr.protein.noHighRNAVar2.Loss$data<- factor(Corr.protein.noHighRNAVar2.Loss$data, level=Corr.protein.noHighRNAVar2.Gain$data)

ggplot(Corr.protein.noHighRNAVar2.Gain, aes(x=data, y=ROCauc))+ #5x5 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.45, 0.6))+
  ggtitle("Buffering upon chrm gain, factors")
# 5x4
# plot.ROCauc.Protein.gain.Buffer_AS0.25_v1

ggplot(Corr.protein.noHighRNAVar2.Loss, aes(x=data, y=ROCauc))+ #5x5 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.45, 0.6))+
  ggtitle("Buffering upon chrm loss factors")
# 5x4
# plot.ROCauc.Protein.loss.Buffer_AS0.25_v1

## Test of AUC significantly different or not: 
roc.test(auc(All.Factors.Diff$Three.Protein.Gain=="Buffering", 
             All.Factors.Diff$Non.exponential.decay.delta), 
         auc(All.Factors.Diff$Three.Protein.Loss=="Buffering", 
             All.Factors.Diff$Non.exponential.decay.delta), 
         reuse.auc=FALSE, paired=FALSE, method="venkatraman")
# p-value = 0.0035 for original NED ROCs
# p-value = 0.004 for updated -0.25 to 0.25 cutoff NED ROCs
# test type: Venkatraman comparison of unpaired ROC curves


#### plot dataset specific factors##
Corr.protein.ROCAUC.dataset<-Corr.protein.ROCAUC
Corr.protein.ROCAUC.dataset<-Corr.protein.ROCAUC.dataset[c(55:66, 151:155, 169:186),] #get RNA/Protein neutral variance, & gene amp ratio

#Subset data by type: 
Corr.protein.ROCAUC.dataset.Gain<- subset(Corr.protein.ROCAUC.dataset, Chrm_CN=="Gain" & Gene_Type=="Protein" & Category=="Buffering")
Corr.protein.ROCAUC.dataset.Gain<- Corr.protein.ROCAUC.dataset.Gain[order(Corr.protein.ROCAUC.dataset.Gain$ROCauc),]

Corr.protein.ROCAUC.dataset.Loss<- subset(Corr.protein.ROCAUC.dataset, Chrm_CN=="Loss" & Gene_Type=="Protein" & Category=="Buffering")
Corr.protein.ROCAUC.dataset.Loss<- Corr.protein.ROCAUC.dataset.Loss[order(Corr.protein.ROCAUC.dataset.Loss$ROCauc),]

## Buffering ROC AUC bargraph
Corr.protein.ROCAUC.dataset.Gain$data<- factor(Corr.protein.ROCAUC.dataset.Gain$data, level=Corr.protein.ROCAUC.dataset.Gain$data)
Corr.protein.ROCAUC.dataset.Loss$data<- factor(Corr.protein.ROCAUC.dataset.Loss$data, level=Corr.protein.ROCAUC.dataset.Gain$data)

ggplot(Corr.protein.ROCAUC.dataset.Gain, aes(x=data, y=ROCauc))+ #5x4 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.45, 0.65))+
  ggtitle("Buffering upon chrm gain, dataset specific factors")
# 5x4
# plot.ROCauc.Protein.gain.Buffer.datasetspecific

ggplot(Corr.protein.ROCAUC.dataset.Loss, aes(x=data, y=ROCauc))+ #5x4 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.45, 0.65))+
  ggtitle("Buffering upon chrm loss, dataset specific factors")
# 5x4
# plot.ROCauc.Protein.loss.Buffer.datasetspecific






### Anti Scaling ##
#Corr.protein.noHighRNAVar2<-Corr.protein.noHighRNAVar

Corr.protein.AS.Gain<- subset(Corr.protein.noHighRNAVar2, Chrm_CN=="Gain" & Gene_Type=="Protein" & Category=="Anti-Scaling")
Corr.protein.AS.Gain<- Corr.protein.AS.Gain[order(Corr.protein.AS.Gain$ROCauc),]

Corr.protein.AS.Loss<- subset(Corr.protein.noHighRNAVar2, Chrm_CN=="Loss" & Gene_Type=="Protein" & Category=="Anti-Scaling")
Corr.protein.AS.Loss<- Corr.protein.AS.Loss[order(Corr.protein.AS.Loss$ROCauc),]

ggplot(Corr.protein.AS.Gain, aes(x=reorder(data, ROCauc), y=ROCauc))+ #5x5 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.40, 0.65))+
  ggtitle("Anti-Scaling upon chrm gain, factors")
# 5x4
# plot.ROCauc.Protein.gain.AS_v1

ggplot(Corr.protein.AS.Loss, aes(x=reorder(data, ROCauc), y=ROCauc))+ #5x5 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.40, 0.65))+
  ggtitle("Anti-Scaling upon chrm loss factors")
# 5x4
# plot.ROCauc.Protein.loss.AS_v1




###Scaling ###
#Corr.protein.noHighRNAVar2<-Corr.protein.noHighRNAVar

Corr.protein.S.Gain<- subset(Corr.protein.noHighRNAVar2, Chrm_CN=="Gain" & Gene_Type=="Protein" & Category=="Scaling")
Corr.protein.S.Gain<- Corr.protein.S.Gain[order(Corr.protein.S.Gain$ROCauc),]

Corr.protein.S.Loss<- subset(Corr.protein.noHighRNAVar2, Chrm_CN=="Loss" & Gene_Type=="Protein" & Category=="Scaling")
Corr.protein.S.Loss<- Corr.protein.S.Loss[order(Corr.protein.S.Loss$ROCauc),]

ggplot(Corr.protein.S.Gain, aes(x=reorder(data, ROCauc), y=ROCauc))+ #5x5 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.44, 0.6))+
  ggtitle("Scaling upon chrm gain, factors")
# 5x4
# plot.ROCauc.Protein.gain.Scaling_v1

ggplot(Corr.protein.S.Loss, aes(x=reorder(data, ROCauc), y=ROCauc))+ #5x5 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.44, 0.6))+
  ggtitle("Scaling upon chrm loss factors")
# 5x4
# plot.ROCauc.Protein.loss.Scaling_v1



#####         PLOT ROC AUC curve for : buffer, Anti-Scaling, and Scale ####
# automatically plot all data. Put pdfs in below folder: 
#setwd()#set working directory to where you downloaded this data

### EDIT! change name and test collumn you want to test as needed: 
## Get ROC AUC for specific factor of interest. 

name <- "Reproducibility rank" #name of factor
TestDataCol <- All.Factors.Diff$`Aggregated Reproducibility Rank` #factor you want to see ROC curve of
# Category== Buffering
# Automatically make all the plots. 
# 3x3
pdf(file = paste0("Protein.Gain.ROC.Buffer.", name, ".pdf"),
    width = 3, 
    height = 3)
roc(All.Factors.Diff$Three.Protein.Gain=="Buffering", TestDataCol, 
    smoothed = TRUE,
    plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
    print.auc=TRUE, show.thres=TRUE)
dev.off()


pdf(file = paste0("Protein.Loss.ROC.Buffer.", name, ".pdf"),
    width = 3, 
    height = 3)
roc(All.Factors.Diff$Three.Protein.Loss=="Buffering", TestDataCol,
    smoothed = TRUE,
    plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
    print.auc=TRUE, show.thres=TRUE)
dev.off()




#Anti scaling
pdf(file = paste0("Protein.Gain.ROC.AS.", name, ".pdf"),
    width = 3, 
    height = 3)
roc(All.Factors.Diff$Three.Protein.Gain=="Anti-Scaling", TestDataCol, 
    smoothed = TRUE,
    plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
    print.auc=TRUE, show.thres=TRUE)
dev.off()


pdf(file = paste0("Protein.Loss.ROC.AS.", name, ".pdf"),
    width = 3, 
    height = 3)
roc(All.Factors.Diff$Three.Protein.Loss=="Anti-Scaling", TestDataCol,
    smoothed = TRUE,
    plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
    print.auc=TRUE, show.thres=TRUE)
dev.off()

# Scaling
pdf(file = paste0("Protein.Gain.ROC.Scaling.", name, ".pdf"),
    width = 3, 
    height = 3)
roc(All.Factors.Diff$Three.Protein.Gain=="Scaling", TestDataCol, 
    smoothed = TRUE,
    plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
    print.auc=TRUE, show.thres=TRUE)
dev.off()


pdf(file = paste0("Protein.Loss.ROC.Scaling.", name, ".pdf"),
    width = 3, 
    height = 3)
roc(All.Factors.Diff$Three.Protein.Loss=="Scaling", TestDataCol,
    smoothed = TRUE,
    plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
    print.auc=TRUE, show.thres=TRUE)
dev.off()

#####         PLOT correlation coefficient, color significance #####
# Corr.protein.abs

# seperate data by data upon chrm gain, and loss
Corr.factor.gain<- subset(Corr.protein.abs, Chrm_CN=="Gain") 
Corr.factor.loss<- subset(Corr.protein.abs, Chrm_CN=="Loss")
Corr.factor.gain<-Corr.factor.gain[-c(9:11, 26, 29:31),] #get only dataset specific factors
Corr.factor.loss<-Corr.factor.loss[-c(9:11, 26, 29:31),] #get only dataset specific factors

# order data by largest correlation upon chrm gain
Corr.factor.gain<- Corr.factor.gain[order(Corr.factor.gain$Corr_coef),]
Corr.factor.loss<- Corr.factor.loss[order(Corr.factor.loss$Corr_coef),]

#factor data so you can plot it accordingly
Corr.factor.gain$data<- factor(Corr.factor.gain$data, level=Corr.factor.gain$data)
Corr.factor.loss$data<- factor(Corr.factor.loss$data, level=Corr.factor.loss$data)

ggplot(Corr.factor.gain, aes(x=data, y=Corr_coef, fill= significant))+ 
  geom_bar(stat="identity")+
  ylab("Correlation coefficient")+
  xlab("")+ 
  scale_fill_manual(values=c("black", "Gold2"))+
  theme_classic()+
  geom_hline(yintercept=0, color="black")+
  coord_flip(ylim=c(-0.2,0.2))+
  ggtitle("Factors for chrm gain, corr")
# plot.Difference.Factors.Corr.Gain
# 5x5

ggplot(Corr.factor.loss, aes(x=data, y=Corr_coef, fill= significant))+ #5x5 figure
  geom_bar(stat="identity")+
  ylab("Correlation coefficient")+
  scale_fill_manual(values=c("black", "Gold2"))+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0, color="black")+
  coord_flip(ylim=c(-0.2,0.2))+
  ggtitle("Factors for chrm loss, corr")
# plot.Difference.Factors.Corr.Loss
# 5x5



### Correlation coefficint graphs for dataset specific factors

# seperate data by data upon chrm gain, and loss
Corr.factor.gain.dataset<- subset(Corr.protein.abs, Chrm_CN=="Gain") 
Corr.factor.loss.dataset<- subset(Corr.protein.abs, Chrm_CN=="Loss")
Corr.factor.gain.dataset<-Corr.factor.gain.dataset[-c(1:9, 12:25, 27:28),] #get only dataset specific factors
Corr.factor.loss.dataset<-Corr.factor.loss.dataset[-c(1:9, 12:25, 27:28),] #get only dataset specific factors


# order data by largest correlation upon chrm gain
Corr.factor.gain.dataset<- Corr.factor.gain.dataset[order(Corr.factor.gain.dataset$Corr_coef),]
Corr.factor.loss.dataset<- Corr.factor.loss.dataset[order(Corr.factor.loss.dataset$Corr_coef),]

#factor data so you can plot it accordingly
Corr.factor.gain.dataset$data<- factor(Corr.factor.gain.dataset$data, level=Corr.factor.gain.dataset$data)
Corr.factor.loss.dataset$data<- factor(Corr.factor.loss.dataset$data, level=Corr.factor.loss.dataset$data)

ggplot(Corr.factor.gain.dataset, aes(x=data, y=Corr_coef, fill= significant))+ 
  geom_bar(stat="identity")+
  ylab("Correlation coefficient")+
  xlab("")+ 
  scale_fill_manual(values=c("black", "Gold2"))+
  theme_classic()+
  geom_hline(yintercept=0, color="black")+
  coord_flip(ylim=c(-0.2,0.2))+
  ggtitle("Factors for chrm gain, corr")
# plot.Difference.Factors.Corr.Gain.dataset
# 5x3

ggplot(Corr.factor.loss.dataset, aes(x=data, y=Corr_coef, fill= significant))+ #5x5 figure
  geom_bar(stat="identity")+
  ylab("Correlation coefficient")+
  scale_fill_manual(values=c("black", "Gold2"))+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0, color="black")+
  coord_flip(ylim=c(-0.2,0.2))+
  ggtitle("Factors for chrm loss, corr")
# plot.Difference.Factors.Corr.Loss.dataset
# 5x3


##### Ovarian cancer ROC analysis #####
## get ovarian cancer data from TCGA_ProteinDiff_Ovarian.R 

t.test_TCGA_ovarian.ThreeGroups<- read.csv("Protein_TCGA_Difference_Pvalue_min10points_3cat.csv")
All.Factors.Diff.Ovarian<-All.Factors.Diff[-c(1,3,4, 6:21)] #Remove mutations, too specific to CCLE data

All.Factors.Diff.Ovarian2<-merge(x=t.test_TCGA_ovarian.ThreeGroups, y=All.Factors.Diff.Ovarian, 
                                 by.x= "Gene", by.y= "Gene_Symbol")

NameTestCollumns= c("PPI.Freq", "gene.Freq", "regulatory.Freq", "Sumoylation.Freq", 
                    "ubiquitination.Freq", "Phosphorylation.Freq", "Methylation.Freq", "Acetylation.Freq", 
                    "AmplificationRatio", "log2.CoeffVar.RNA", "log2.Protein.CoeffVar", "Rate_3", 
                    "m", "p", "lp", "lm", "bm", "bp", 
                    "Mean_HalfLife", "Loops.in.Protein.ANCHOR.", "Intrinsic_Disorder_MobiBD", "Intrinsic_Disorder_LowComplexity", 
                    "Intrinsic_Disorder_Polarity", "Intrinsic_Disorder_Polyampholyte", "Intrinsic_Disorder_homology", 
                    "Dependency.Score", "Aggregation.score", "Non.exponential.decay.delta", "Mutations.Freq", "mutations.FrameShift.Freq", "mutations.nonsense.missense.Freq")

testdatanames= c("Protein-protein interaction", "Protein complex (CORUM)", "Protein regulatory sites", "Sumoylation sites", 
                 "Ubiquitination sites", "Phosphorylation sites", "Methylation sites", "Acetylation sites", 
                 "Percent gene amplification", "RNA neutral variance", "Protein neutral variance", "mRNA decay rate", 
                 "mRNA abundance", "Protein abundance", "Protein length", "mRNA length", "Transcription rate", "Translation rate", 
                 "Protein half life", "Loops in protein score", "Intrinsic protein disorder", "Low complexity score", 
                 "Protein polarity", "Protein polyampholyte score", "Homology score", 
                 "Dependency score", "Aggregation score", "Non-exponential decay delta", "Mutation count (all)", "Mutation count (frame shifts)", "Mutation count (nonsence and missense)")

FindROCAU.ovarian2<- function(dataset, NameTestCollumns, testdatanames) {
  Corr.protein.loop<-data.frame(Data=as.character(), 
                                Chrm_CN=character(), 
                                Gene_Type=character(), 
                                Category=character(),
                                p_value=as.numeric(), 
                                Corr_coef=as.numeric(), 
                                significant=character(), 
                                ROCauc=as.numeric())
  
  
  for (i in 1:length(NameTestCollumns)){
    
    TestDataCol= NameTestCollumns[i]
    name= testdatanames[i]
    TestData<- as.numeric(as.character(dataset %>% pull(TestDataCol)))
    
    Gain.Protein.corr<- cor.test(TestData, dataset$Diff.Gain, method="pearson")
    Loss.Protein.corr<-cor.test(TestData, dataset$Diff.Loss, method="pearson")
    # Buffering, Gain: 
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Gain", 
                                        Gene_Type="Protein", 
                                        Category="Buffering",
                                        p_value=Gain.Protein.corr$p.value, 
                                        Corr_coef=Gain.Protein.corr$estimate, 
                                        significant=Gain.Protein.corr$p.value<0.05, 
                                        ROCauc=auc(dataset$Three.Protein.Gain.x=="Buffering", 
                                                   TestData, na.rm=TRUE ))) 
    
    # Buffering, Loss: 
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Loss", 
                                        Gene_Type="Protein", 
                                        Category="Buffering",
                                        p_value=Loss.Protein.corr$p.value, 
                                        Corr_coef=Loss.Protein.corr$estimate, 
                                        significant=Loss.Protein.corr$p.value<0.01, 
                                        ROCauc=auc(dataset$Three.Protein.Loss.x=="Buffering", 
                                                   TestData, na.rm=TRUE )))
    
    # AS, Gain: 
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Gain", 
                                        Gene_Type="Protein", 
                                        Category="Anti-Scaling",
                                        p_value=Gain.Protein.corr$p.value, 
                                        Corr_coef=Gain.Protein.corr$estimate, 
                                        significant=Gain.Protein.corr$p.value<0.05, 
                                        ROCauc=auc(dataset$Three.Protein.Gain.x=="Anti-Scaling", 
                                                   TestData, na.rm=TRUE ))) 
    
    # AS, Loss: 
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Loss", 
                                        Gene_Type="Protein", 
                                        Category="Anti-Scaling",
                                        p_value=Loss.Protein.corr$p.value, 
                                        Corr_coef=Loss.Protein.corr$estimate, 
                                        significant=Loss.Protein.corr$p.value<0.01, 
                                        ROCauc=auc(dataset$Three.Protein.Loss.x=="Anti-Scaling", 
                                                   TestData, na.rm=TRUE ))) 
    
    # Scaling, Gain: 
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Gain", 
                                        Gene_Type="Protein", 
                                        Category="Scaling",
                                        p_value=Gain.Protein.corr$p.value, 
                                        Corr_coef=Gain.Protein.corr$estimate, 
                                        significant=Gain.Protein.corr$p.value<0.05, 
                                        ROCauc=auc(dataset$Three.Protein.Gain.x=="Scaling", 
                                                   TestData, na.rm=TRUE ))) 
    
    # Scaling, Loss: 
    Corr.protein.loop<-rbind(Corr.protein.loop, 
                             data.frame(data=name, 
                                        Chrm_CN="Loss", 
                                        Gene_Type="Protein", 
                                        Category="Scaling",
                                        p_value=Loss.Protein.corr$p.value, 
                                        Corr_coef=Loss.Protein.corr$estimate, 
                                        significant=Loss.Protein.corr$p.value<0.01, 
                                        ROCauc=auc(dataset$Three.Protein.Loss.x=="Scaling", 
                                                   TestData, na.rm=TRUE ))) 
    
  }
  return(Corr.protein.loop)
}

Ovarian.protein.ROCAUC<- FindROCAU.ovarian2(All.Factors.Diff.Ovarian2, NameTestCollumns, testdatanames)



###    Plot ROC AUC 
## But first Add 5' and 3' UTR  AUC to data 
NCBI.genedata.NM.5UTR$UTR5 #5' UTR data, unique entries, 33k
TestData<- NCBI.genedata.NM.5UTR
name<- "5' UTR length" 

Diff.test<-merge(x=TestData, y= All.Factors.Diff.Ovarian2, 
                 by.x="name2", by.y="Gene", 
                 sort = TRUE)
TestDataCol<-Diff.test$UTR5

Gain.Protein.corr<- cor.test(TestDataCol, Diff.test$Diff.Gain, method="pearson")
Ovarian.protein.ROCAUC<-rbind(Ovarian.protein.ROCAUC, 
                                  data.frame(data=name, 
                                             Chrm_CN="Gain", 
                                             Gene_Type="Protein", 
                                             Category= "Buffering",
                                             p_value=Gain.Protein.corr$p.value, 
                                             Corr_coef=Gain.Protein.corr$estimate, 
                                             significant=Gain.Protein.corr$p.value<0.05, 
                                             ROCauc=auc(Diff.test$Three.Protein.Gain.x=="Buffering", 
                                                        TestDataCol ))) 
Loss.Protein.corr<- cor.test(TestDataCol, Diff.test$Diff.Loss, method="pearson")
Ovarian.protein.ROCAUC<-rbind(Ovarian.protein.ROCAUC, 
                                  data.frame(data=name, 
                                             Chrm_CN="Loss", 
                                             Gene_Type="Protein", 
                                             Category= "Buffering",
                                             p_value=Loss.Protein.corr$p.value, 
                                             Corr_coef=Loss.Protein.corr$estimate, 
                                             significant=Loss.Protein.corr$p.value<0.01, 
                                             ROCauc=auc(Diff.test$Three.Protein.Loss.x=="Buffering", 
                                                        TestDataCol ))) 


## Now add 3' UTR 
NCBI.genedata.NM.3UTR$UTR3#3'UTR data, unique entries, 26k
TestData<- NCBI.genedata.NM.3UTR
name<- "3' UTR length" 

Diff.test<-merge(x=TestData, y= All.Factors.Diff.Ovarian2, 
                 by.x="name2", by.y="Gene", 
                 sort = TRUE)
TestDataCol<-Diff.test$UTR3

Gain.Protein.corr<- cor.test(TestDataCol, Diff.test$Diff.Gain, method="pearson")
Ovarian.protein.ROCAUC<-rbind(Ovarian.protein.ROCAUC, 
                                  data.frame(data=name, 
                                             Chrm_CN="Gain", 
                                             Gene_Type="Protein", 
                                             Category= "Buffering",
                                             p_value=Gain.Protein.corr$p.value, 
                                             Corr_coef=Gain.Protein.corr$estimate, 
                                             significant=Gain.Protein.corr$p.value<0.05, 
                                             ROCauc=auc(Diff.test$Three.Protein.Gain.x=="Buffering", 
                                                        TestDataCol ))) 
Loss.Protein.corr<- cor.test(TestDataCol, Diff.test$Diff.Loss, method="pearson")
Ovarian.protein.ROCAUC<-rbind(Ovarian.protein.ROCAUC, 
                                  data.frame(data=name, 
                                             Chrm_CN="Loss", 
                                             Gene_Type="Protein", 
                                             Category= "Buffering",
                                             p_value=Loss.Protein.corr$p.value, 
                                             Corr_coef=Loss.Protein.corr$estimate, 
                                             significant=Loss.Protein.corr$p.value<0.01, 
                                             ROCauc=auc(Diff.test$Three.Protein.Loss.x=="Buffering", 
                                                        TestDataCol ))) 

write.csv(Ovarian.protein.ROCAUC, file="Factor.ROCAUC.correlationscore_Ovarian.csv")
## ! OR you can download the file from supplementary: 
Ovarian.protein.ROCAUC<- read.csv(file="Factor.ROCAUC.correlationscore_Ovarian.csv")


## now that 5' and 3' UTR have been added, remove dataset specific factors: 
Ovarian.protein.ROCAUC2<-Ovarian.protein.ROCAUC
Ovarian.protein.ROCAUC2<-Ovarian.protein.ROCAUC2[-c(49:66, 151:156, 169:186),] #remove RNA/Protein neutral variance, & gene amp ratio


#Subset data by type: 
Ovarian.protein.ROCAUC2.Gain<- subset(Ovarian.protein.ROCAUC2, Chrm_CN=="Gain" & Gene_Type=="Protein" & Category=="Buffering")
Ovarian.protein.ROCAUC2.Gain<- Ovarian.protein.ROCAUC2.Gain[order(Ovarian.protein.ROCAUC2.Gain$ROCauc),]

Ovarian.protein.ROCAUC2.Loss<- subset(Ovarian.protein.ROCAUC2, Chrm_CN=="Loss" & Gene_Type=="Protein" & Category=="Buffering")
Ovarian.protein.ROCAUC2.Loss<- Ovarian.protein.ROCAUC2.Loss[order(Ovarian.protein.ROCAUC2.Loss$ROCauc),]


## Buffering ROC AUC bargraph
Ovarian.protein.ROCAUC2.Gain$data<- factor(Ovarian.protein.ROCAUC2.Gain$data, level=Ovarian.protein.ROCAUC2.Gain$data)
Ovarian.protein.ROCAUC2.Loss$data<- factor(Ovarian.protein.ROCAUC2.Loss$data, level=Ovarian.protein.ROCAUC2.Loss$data)

ggplot(Ovarian.protein.ROCAUC2.Gain, aes(x=data, y=ROCauc))+ #5x5 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.45, 0.6))+
  ggtitle("Buffering upon chrm gain, factors")
# 5x4
# plot.ROCauc.Protein.gain.Buffer_Ovarian

ggplot(Ovarian.protein.ROCAUC2.Loss, aes(x=data, y=ROCauc))+ #5x5 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.45, 0.6))+
  ggtitle("Buffering upon chrm loss factors")
# 5x4
# plot.ROCauc.Protein.loss.Buffer_Ovarian



#### plot dataset specific factors ##
## No mutations because mutation data is very dataset specific
Ovarian.protein.ROCAUC.dataset<-Ovarian.protein.ROCAUC
Ovarian.protein.ROCAUC.dataset<-Ovarian.protein.ROCAUC.dataset[c(55:66, 151:155),] #get RNA/Protein neutral variance, & gene amp ratio

#Subset data by type: 
Ovarian.protein.ROCAUC.dataset.Gain<- subset(Ovarian.protein.ROCAUC.dataset, Chrm_CN=="Gain" & Gene_Type=="Protein" & Category=="Buffering")
Ovarian.protein.ROCAUC.dataset.Gain<- Ovarian.protein.ROCAUC.dataset.Gain[order(Ovarian.protein.ROCAUC.dataset.Gain$ROCauc),]

Ovarian.protein.ROCAUC.dataset.Loss<- subset(Ovarian.protein.ROCAUC.dataset, Chrm_CN=="Loss" & Gene_Type=="Protein" & Category=="Buffering")
Ovarian.protein.ROCAUC.dataset.Loss<- Ovarian.protein.ROCAUC.dataset.Loss[order(Ovarian.protein.ROCAUC.dataset.Loss$ROCauc),]

## Buffering ROC AUC bargraph
Ovarian.protein.ROCAUC.dataset.Gain$data<- factor(Ovarian.protein.ROCAUC.dataset.Gain$data, level=Ovarian.protein.ROCAUC.dataset.Gain$data)
Ovarian.protein.ROCAUC.dataset.Loss$data<- factor(Ovarian.protein.ROCAUC.dataset.Loss$data, level=Ovarian.protein.ROCAUC.dataset.Loss$data)

ggplot(Ovarian.protein.ROCAUC.dataset.Gain, aes(x=data, y=ROCauc))+ #5x4 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.45, 0.65))+
  ggtitle("Ovarian buffering upon chrm gain, dataset specific factors")
# 5x4
# plot.ROCauc.Protein.gain.Buffer.datasetspecific_Ovarian

ggplot(Ovarian.protein.ROCAUC.dataset.Loss, aes(x=data, y=ROCauc))+ #5x4 figure
  geom_bar(stat="identity")+
  ylab("ROC AUC")+
  xlab("")+ 
  theme_classic()+
  geom_hline(yintercept=0.5, color="black")+
  coord_flip(ylim=c(0.45, 0.65))+
  ggtitle("Ovarian buffering upon chrm loss, dataset specific factors")
# 5x4
# plot.ROCauc.Protein.loss.Buffer.datasetspecific_Ovarian



##### Random Permutation test to add significance to Buffering Factor ROC AUC #####
# 210908
# We want to do random permutations, >1000 permutations per buffering factor
# then we test what the ROC values are for 1000x permutations per factor
# then we see if the actual factor ROC is significantly greater than 95% of random permutation. 
# which ones are greater than random? 

All.Factors.Diff

NameTestCollumns= c("PPI.Freq", "gene.Freq", "regulatory.Freq", "Sumoylation.Freq", 
                    "ubiquitination.Freq", "Phosphorylation.Freq", "Methylation.Freq", "Acetylation.Freq", 
                    "AmplificationRatio", "log2.CoeffVar.RNA", "log2.Protein.CoeffVar", "Rate_3", 
                    "m", "p", "lp", "lm", "bm", "bp", 
                    "Mean_HalfLife", "Loops.in.Protein.ANCHOR.", "Intrinsic_Disorder_MobiBD", "Intrinsic_Disorder_LowComplexity", 
                    "Intrinsic_Disorder_Polarity", "Intrinsic_Disorder_Polyampholyte", "Intrinsic_Disorder_homology", 
                    "Dependency.Score", "Aggregation.score", "Non.exponential.decay.delta", "Mutations.Freq", "mutations.FrameShift.Freq", "mutations.nonsense.missense.Freq")

testdatanames= c("Protein-protein interaction", "Protein complex (CORUM)", "Protein regulatory sites", "Sumoylation sites", 
                 "Ubiquitination sites", "Phosphorylation sites", "Methylation sites", "Acetylation sites", 
                 "Percent gene amplification", "RNA neutral variance", "Protein neutral variance", "mRNA decay rate", 
                 "mRNA abundance", "Protein abundance", "Protein length", "mRNA length", "Transcription rate", "Translation rate", 
                 "Protein half life", "Loops in protein score", "Intrinsic protein disorder", "Low complexity score", 
                 "Protein polarity", "Protein polyampholyte score", "Homology score", 
                 "Dependency score", "Aggregation score", "Non-exponential decay delta", "Mutation count (all)", "Mutation count (frame shifts)", "Mutation count (nonsence and missense)")


dataset= All.Factors.Diff


## Basic Permutation test Function
# where N = number of repeats for test
# where i = one for each factor
N= 10000

PermuteFunction <- function(y=All.Factors.Diff[,i], x=All.Factors.Diff$Three.Protein.Gain){
  TestDataCol= NameTestCollumns[i]
  name= testdatanames[i]
  TestData<- as.numeric(as.character(dataset %>% pull(TestDataCol)))
  
  
  model.ROC= auc(sample(x=="Buffering", replace=F), 
                 TestData, na.rm=TRUE )
  # permutate category, then runs model
  return(model.ROC)
}
PermuteFunction.ChrmLoss <- function(y=All.Factors.Diff[,i], x=All.Factors.Diff$Three.Protein.Loss){
  TestDataCol= NameTestCollumns[i]
  name= testdatanames[i]
  TestData<- as.numeric(as.character(dataset %>% pull(TestDataCol)))
  
  
  model.ROC= auc(sample(x=="Buffering", replace=F), 
                 TestData, na.rm=TRUE )
  # permutate category, then runs model
  return(model.ROC)
}

#### For chromosome GAIN: 
# to perform the iterations of te permutate function we can use replicate()
# i=1
# Permute_N<- replicate(N, PermuteFunction())

# now to perform iterations of random permutations for each factor.
Random.Permute.Gain<-data.frame(Factor=as.character(), 
                                Chrm_CN=as.character(),
                                RandomPermutation=c(), 
                                Top5percent=numeric(), 
                                Bottom5percent=numeric())

for (i in 1:length(testdatanames)){
  Permute_N<- replicate(N, PermuteFunction())
  Random.Permute.Gain<- rbind (Random.Permute.Gain, 
                                 data.frame(Factor= testdatanames[i], 
                                            Chrm_CN="Gain",
                                            RandomPermutation= I(list(Permute_N)), 
                                            Top5percent=sort(Permute_N, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                            Bottom5percent=sort(Permute_N)[500]) ) # bottom 500th value is lower 5% cutoff for 10000 values
  
}

# Now see if our ROC value was greater than random, 95% of the time. 

Gain.ROC<- subset(Corr.protein.ROCAUC, Chrm_CN=="Gain" & Category=="Buffering")

Random.Permute.Gain2<- merge(Random.Permute.Gain, Gain.ROC, 
                             by.x= "Factor", by.y= "data")

# now get significance cutoffs: 
Random.Permute.Gain2$Significance<-NA
for (w in 1:length(Random.Permute.Gain2$Factor)){
  cutoffs<-quantile(Random.Permute.Gain2$RandomPermutation[[w]], c(.95, .995, .9995))
  if (Random.Permute.Gain2$ROCauc[w]> cutoffs[3]){
    Random.Permute.Gain2$Significance[w]<-"***"
  } else if (Random.Permute.Gain2$ROCauc[w]> cutoffs[2]){
    Random.Permute.Gain2$Significance[w]<-"**"
  } else if (Random.Permute.Gain2$ROCauc[w]> cutoffs[1]){
    Random.Permute.Gain2$Significance[w]<-"*"
  } else {
    Random.Permute.Gain2$Significance[w]<-"NS"
  }
}

#write.csv(Random.Permute.Gain2, file="/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/ROC.Factors.RandomPermutationTest.ChrmGain.csv")
#Where "significant is TRUE/FALSE value for whether or not correlation coefficient p-value is <0.05
# while significance is the threshold for random permutation analtsis with 10k random permutations
# and whtehr the ROC value is >5% of random ROC permutations (*), >0.5% of random permutations (**), or >0.05% of random permutations
# This does not open well in R. So I manually removed the list of random permutation ROC values so that the rest of the data is accessable.
# renamed file: ROC.Factors.RandomPermutationTest.ChrmGain_noPermutationList.csv
Random.Permute.Gain3<-read.csv("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/ROC.Factors.RandomPermutationTest.ChrmGain_noPermutationList.csv")

## plot Protein-Protein Interaction histogram: 
ggplot(data.frame(Random.Permute.Gain2$RandomPermutation[[26]]), aes(x=Random.Permute.Gain2.RandomPermutation..26..))+
  geom_histogram(bins=50)+
  theme_classic()+
  xlab("Random Permutation ROC AUC for\nProtein-Protein Interactions, chromosome gain")+
  ylab("Frequency")+
  geom_vline(xintercept= Random.Permute.Gain2$ROCauc[26], linetype="dotted")+
  geom_vline(xintercept= 0.5)+
  scale_x_continuous(limits=c(0.45, 0.60))
## plot.RandomPermutation.Factor.ROC.PPI.Gain
# 4x4

## plot Ubiquitination histogram: 
ggplot(data.frame(Random.Permute.Gain2$RandomPermutation[[31]]), aes(x=Random.Permute.Gain2.RandomPermutation..31..))+
  geom_histogram(bins=50)+
  theme_classic()+
  xlab("Random Permutation ROC AUC for\nUbiquitination sites, chromosome gain")+
  ylab("Frequency")+
  geom_vline(xintercept= Random.Permute.Gain2$ROCauc[31], linetype="dotted")+
  geom_vline(xintercept= 0.5)+
  scale_x_continuous(limits=c(0.45, 0.60))
## plot.RandomPermutation.Factor.ROC.Ubiquitination.Gain
# 4x4



### For Chromosome LOSS: 
PermuteFunction.ChrmLoss <- function(y=All.Factors.Diff[,i], x=All.Factors.Diff$Three.Protein.Loss){
  TestDataCol= NameTestCollumns[i]
  name= testdatanames[i]
  TestData<- as.numeric(as.character(dataset %>% pull(TestDataCol)))
  
  
  model.ROC= auc(sample(x=="Buffering", replace=F), 
                 TestData, na.rm=TRUE )
  # permutate category, then runs model
  return(model.ROC)
}

# to perform the iterations of te permutate function we can use replicate()
# i=1
# Permute_N<- replicate(N, PermuteFunction())

# now to perform iterations of random permutations for each factor.
Random.Permute.Loss<-data.frame(Factor=as.character(), 
                                Chrm_CN=as.character(),
                                RandomPermutation=c(), 
                                Top5percent=numeric(), 
                                Bottom5percent=numeric())

for (i in 1:length(testdatanames)){
  Permute_N<- replicate(N, PermuteFunction())
  Random.Permute.Loss<- rbind (Random.Permute.Loss, 
                               data.frame(Factor= testdatanames[i], 
                                          Chrm_CN="Loss",
                                          RandomPermutation= I(list(Permute_N)), 
                                          Top0.5percent=sort(Permute_N, decreasing=TRUE)[50], # top 500th value is top 5% cutoff for 10000 values
                                          Bottom0.5percent=sort(Permute_N)[50]) ) # bottom 500th value is lower 5% cutoff for 10000 values
  
}

# Now see if our ROC value was greater than random, 95% of the time. 

Loss.ROC<- subset(Corr.protein.ROCAUC, Chrm_CN=="Loss" & Category=="Buffering")

Random.Permute.Loss2<- merge(Random.Permute.Loss, Loss.ROC, 
                             by.x= "Factor", by.y= "data")

# now get significance cutoffs: 
Random.Permute.Loss2$Significance<-NA
for (w in 1:length(Random.Permute.Loss2$Factor)){
  cutoffs<-quantile(Random.Permute.Loss2$RandomPermutation[[w]], c(.95, .995, .9995))
  if (Random.Permute.Loss2$ROCauc[w]> cutoffs[3]){
    Random.Permute.Loss2$Significance[w]<-"***"
  } else if (Random.Permute.Loss2$ROCauc[w]> cutoffs[2]){
    Random.Permute.Loss2$Significance[w]<-"**"
  } else if (Random.Permute.Loss2$ROCauc[w]> cutoffs[1]){
    Random.Permute.Loss2$Significance[w]<-"*"
  } else {
    Random.Permute.Loss2$Significance[w]<-"NS"
  }
}

# write.csv(Random.Permute.Loss2, 
#          file="/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/ROC.Factors.RandomPermutationTest.ChrmLoss.csv")
# Where "significant is TRUE/FALSE value for whether or not correlation coefficient p-value is <0.05
# while significance is the threshold for random permutation analtsis with 10k random permutations
# and whtehr the ROC value is >5% of random ROC permutations (*), >0.5% of random permutations (**), or >0.05% of random permutations
# This does not open well in R. So I manually removed the list of random permutation ROC values so that the rest of the data is accessable.
# renamed file: ROC.Factors.RandomPermutationTest.ChrmLoss_noPermutationList.csv
Random.Permute.Loss3<-read.csv("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/ROC.Factors.RandomPermutationTest.ChrmLoss_noPermutationList.csv")

## plot Protein-Protein Interaction histogram: 
ggplot(data.frame(Random.Permute.Loss2$RandomPermutation[[26]]), aes(x=Random.Permute.Loss2.RandomPermutation..26..))+
  geom_histogram(bins=50)+
  theme_classic()+
  xlab("Random Permutation ROC AUC for\nProtein-Protein Interactions, chromosome loss")+
  ylab("Frequency")+
  geom_vline(xintercept= Random.Permute.Loss2$ROCauc[26], linetype="dotted")+
  geom_vline(xintercept= 0.5)+
  scale_x_continuous(limits=c(0.45, 0.60))
## plot.RandomPermutation.Factor.ROC.PPI.Loss
# 4x4

## plot Ubiquitination histogram: 
ggplot(data.frame(Random.Permute.Loss2$RandomPermutation[[31]]), aes(x=Random.Permute.Loss2.RandomPermutation..31..))+
  geom_histogram(bins=50)+
  theme_classic()+
  xlab("Random Permutation ROC AUC for\nUbiquitination sites, chromosome loss")+
  ylab("Frequency")+
  geom_vline(xintercept= Random.Permute.Loss2$ROCauc[31], linetype="dotted")+
  geom_vline(xintercept= 0.5)+
  scale_x_continuous(limits=c(0.45, 0.60))
## plot.RandomPermutation.Factor.ROC.Ubiquitination.Loss
# 4x4

##### Random Permutation test for the 3' and 5' UTR regions ROC AUC #####
# 210928
# I also want to test 5' UTR and 3' UTR data, but there are multiple 5' and 3' datapoints per gene, 
# so I do that seperately so as not to mess with other data values/frequencies. 

### 5' UTR data random permutation test ###
# Step one: get 5' UTR data (from above)

dataset= All.Factors.Diff

## add 5' UTR data
TestData<- NCBI.genedata.NM.5UTR
name<- "5' UTR length" 

Diff.Factor.5UTR<-merge(x=TestData, y= All.Factors.Diff, 
                 by.x="name2", by.y="Gene_Symbol", 
                 sort = TRUE)


## Write the basic Permutation test Function
# where N = number of repeats for test
# where i = one for each factor
N= 10000

PermuteFunction.ChrmGain.5UTR <- function(y=Diff.Factor.5UTR$UTR5, x=Diff.Factor.5UTR$Three.Protein.Gain){
  TestDataCol= "UTR5"
  name= "5' UTR"
  TestData<- as.numeric(as.character(Diff.Factor.5UTR %>% pull(TestDataCol)))
  
  
  model.ROC= auc(sample(x=="Buffering", replace=F), 
                 TestData, na.rm=TRUE )
  # permutate category, then runs model
  return(model.ROC)
}
PermuteFunction.ChrmLoss.5UTR <- function(y=Diff.Factor.5UTR$UTR5, x=Diff.Factor.5UTR$Three.Protein.Loss){
  TestDataCol= "UTR5"
  name= "5' UTR"
  TestData<- as.numeric(as.character(Diff.Factor.5UTR %>% pull(TestDataCol)))
  
  
  model.ROC= auc(sample(x=="Buffering", replace=F), 
                 TestData, na.rm=TRUE )
  # permutate category, then runs model
  return(model.ROC)
}

#### For 5' UTR chromosome GAIN: 
# to perform the iterations of te permutate function we can use replicate()
# i=1
# Permute_N<- replicate(N, PermuteFunction())

Permute_N.5UTRGain<- replicate(N, PermuteFunction.ChrmGain.5UTR())
write.csv(Permute_N.5UTRGain, "RandomPermute.ROC.ChrmGain.Buffering.5UTR.csv")
Random.Permute.Gain.5UTR<- data.frame(Factor= "5' UTR length", 
                                          Chrm_CN="Gain",
                                          Top5percent=sort(Permute_N.5UTRGain, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                          Bottom5percent=sort(Permute_N.5UTRGain)[500])  # bottom 500th value is lower 5% cutoff for 10000 values


# Now see if our ROC value was greater than random, 95% of the time. 
setwd() #set location
Corr.protein.ROCAUC<- read.csv(file="Factor.ROCAUC.correlationscore.csv")

Gain.ROC.5UTR<- subset(Corr.protein.ROCAUC, Chrm_CN=="Gain" & Category=="Buffering")

Random.Permute.Gain.5UTR2<- merge(Random.Permute.Gain.5UTR, Gain.ROC.5UTR, 
                             by.x= "Factor", by.y= "data")

# now get significance cutoffs: 
Random.Permute.Gain.5UTR2$Significance<-NA
for (w in 1:length(Random.Permute.Gain.5UTR2$Factor)){
  cutoffs<-quantile(Permute_N.5UTRGain, c(.95, .995, .9995))
  if (Random.Permute.Gain.5UTR2$ROCauc[w]> cutoffs[3]){
    Random.Permute.Gain.5UTR2$Significance[w]<-"***"
  } else if (Random.Permute.Gain.5UTR2$ROCauc[w]> cutoffs[2]){
    Random.Permute.Gain.5UTR2$Significance[w]<-"**"
  } else if (Random.Permute.Gain.5UTR2$ROCauc[w]> cutoffs[1]){
    Random.Permute.Gain.5UTR2$Significance[w]<-"*"
  } else {
    Random.Permute.Gain.5UTR2$Significance[w]<-"NS"
  }
}
setwd() #set location you want to write data to. 
write.csv(Random.Permute.Gain.5UTR2, file="/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/ROC.Factors.RandomPermutationTest.ChrmGain.5UTR.csv")
# Where "significant is TRUE/FALSE value for whether or not correlation coefficient p-value is <0.05
# while significance is the threshold for random permutation analtsis with 10k random permutations
# and whether the ROC value is >5% of random ROC permutations (*), >0.5% of random permutations (**), or >0.05% of random permutations
# I updated the overview file with data for each factor's permutation analysis. 


### For 5' UTR Chromosome LOSS: 
Permute_N_Loss.5UTR<- replicate(N, PermuteFunction.ChrmLoss.5UTR())
write.csv(Permute_N_Loss.5UTR, "RandomPermute.ROC.ChrmLoss.Buffering.5UTR.csv")
Random.Permute.Loss.5UTR<- data.frame(Factor= "5' UTR length", 
                                          Chrm_CN="Loss",
                                          Top0.5percent=sort(Permute_N_Loss.5UTR, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                          Bottom0.5percent=sort(Permute_N_Loss.5UTR)[500]) # bottom 500th value is lower 5% cutoff for 10000 values


# Now see if our ROC value was greater than random, 95% of the time. 

Loss.ROC.5UTR<- subset(Corr.protein.ROCAUC, Chrm_CN=="Loss" & Category=="Buffering")

Random.Permute.Loss.5UTR2<- merge(Random.Permute.Loss.5UTR, Loss.ROC.5UTR, 
                             by.x= "Factor", by.y= "data")

# now get significance cutoffs: 
Random.Permute.Loss.5UTR2$Significance<-NA
for (w in 1:length(Random.Permute.Loss.5UTR2$Factor)){
  cutoffs<-quantile(Permute_N_Loss.5UTR, c(.95, .995, .9995))
  if (Random.Permute.Loss.5UTR2$ROCauc[w]> cutoffs[3]){
    Random.Permute.Loss.5UTR2$Significance[w]<-"***"
  } else if (Random.Permute.Loss.5UTR2$ROCauc[w]> cutoffs[2]){
    Random.Permute.Loss.5UTR2$Significance[w]<-"**"
  } else if (Random.Permute.Loss.5UTR2$ROCauc[w]> cutoffs[1]){
    Random.Permute.Loss.5UTR2$Significance[w]<-"*"
  } else {
    Random.Permute.Loss.5UTR2$Significance[w]<-"NS"
  }
}

write.csv(Random.Permute.Loss.5UTR2, 
          file="/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/ROC.Factors.RandomPermutationTest.ChrmLoss.5UTR.csv")
# Where "significant is TRUE/FALSE value for whether or not correlation coefficient p-value is <0.05
# while significance is the threshold for random permutation analtsis with 10k random permutations
# and whtehr the ROC value is >5% of random ROC permutations (*), >0.5% of random permutations (**), or >0.05% of random permutations
# I also manually updated list for overview of factor significance. 




### 3' UTR data random permutation test ###
# get 3' UTR data (from above)

dataset= All.Factors.Diff

## add 3' UTR data
TestData<- NCBI.genedata.NM.3UTR
name<- "3' UTR length" 

Diff.Factor.3UTR<-merge(x=TestData, y= All.Factors.Diff, 
                        by.x="name2", by.y="Gene_Symbol", 
                        sort = TRUE)


## Write the basic Permutation test Function
# where N = number of repeats for test
# where i = one for each factor
N= 10000

PermuteFunction.ChrmGain.3UTR <- function(y=Diff.Factor.3UTR$UTR3, x=Diff.Factor.3UTR$Three.Protein.Gain){
  TestDataCol= "UTR3"
  name= "3' UTR"
  TestData<- as.numeric(as.character(Diff.Factor.3UTR %>% pull(TestDataCol)))
  
  
  model.ROC= auc(sample(x=="Buffering", replace=F), 
                 TestData, na.rm=TRUE )
  # permutate category, then runs model
  return(model.ROC)
}
PermuteFunction.ChrmLoss.3UTR <- function(y=Diff.Factor.3UTR$UTR3, x=Diff.Factor.3UTR$Three.Protein.Loss){
  TestDataCol= "UTR3"
  name= "3' UTR"
  TestData<- as.numeric(as.character(Diff.Factor.3UTR %>% pull(TestDataCol)))
  
  
  model.ROC= auc(sample(x=="Buffering", replace=F), 
                 TestData, na.rm=TRUE )
  # permutate category, then runs model
  return(model.ROC)
}

#### For 3' UTR chromosome GAIN: 
# to perform the iterations of te permutate function we can use replicate()
# i=1
# Permute_N<- replicate(N, PermuteFunction())

Permute_N.3UTRGain<- replicate(N, PermuteFunction.ChrmGain.3UTR())
write.csv(Permute_N.3UTRGain, "RandomPermute.ROC.ChrmGain.Buffering.3UTR.csv")
Random.Permute.Gain.3UTR<- data.frame(Factor= "3' UTR length", 
                                      Chrm_CN="Gain",
                                      Top5percent=sort(Permute_N.3UTRGain, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                      Bottom5percent=sort(Permute_N.3UTRGain)[500])  # bottom 500th value is lower 5% cutoff for 10000 values


# Now see if our ROC value was greater than random, 95% of the time. 
# setwd()
#Corr.protein.ROCAUC<- read.csv(file="Factor.ROCAUC.correlationscore.csv")

Gain.ROC.3UTR<- subset(Corr.protein.ROCAUC, Chrm_CN=="Gain" & Category=="Buffering")

Random.Permute.Gain.3UTR2<- merge(Random.Permute.Gain.3UTR, Gain.ROC.3UTR, 
                                  by.x= "Factor", by.y= "data")

# now get significance cutoffs: 
Random.Permute.Gain.3UTR2$Significance<-NA
for (w in 1:length(Random.Permute.Gain.3UTR2$Factor)){
  cutoffs<-quantile(Permute_N.3UTRGain, c(.95, .995, .9995))
  if (Random.Permute.Gain.3UTR2$ROCauc[w]> cutoffs[3]){
    Random.Permute.Gain.3UTR2$Significance[w]<-"***"
  } else if (Random.Permute.Gain.3UTR2$ROCauc[w]> cutoffs[2]){
    Random.Permute.Gain.3UTR2$Significance[w]<-"**"
  } else if (Random.Permute.Gain.3UTR2$ROCauc[w]> cutoffs[1]){
    Random.Permute.Gain.3UTR2$Significance[w]<-"*"
  } else {
    Random.Permute.Gain.3UTR2$Significance[w]<-"NS"
  }
}
setwd()
write.csv(Random.Permute.Gain.3UTR2, file="/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/ROC.Factors.RandomPermutationTest.ChrmGain.3UTR.csv")
# Where "significant is TRUE/FALSE value for whether or not correlation coefficient p-value is <0.05
# while significance is the threshold for random permutation analtsis with 10k random permutations
# and whtehr the ROC value is >5% of random ROC permutations (*), >0.5% of random permutations (**), or >0.05% of random permutations
# I updated the overview file with data for each factor's permutation analysis. 


### For 3' UTR Chromosome LOSS: 
Permute_N_Loss.3UTR<- replicate(N, PermuteFunction.ChrmLoss.3UTR())
write.csv(Permute_N_Loss.3UTR, "RandomPermute.ROC.ChrmLoss.Buffering.3UTR.csv")
Random.Permute.Loss.3UTR<- data.frame(Factor= "3' UTR length", 
                                      Chrm_CN="Loss",
                                      Top0.5percent=sort(Permute_N_Loss.3UTR, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                      Bottom0.5percent=sort(Permute_N_Loss.3UTR)[500]) # bottom 500th value is lower 5% cutoff for 10000 values


# Now see if our ROC value was greater than random, 95% of the time. 

Loss.ROC.3UTR<- subset(Corr.protein.ROCAUC, Chrm_CN=="Loss" & Category=="Buffering")

Random.Permute.Loss.3UTR2<- merge(Random.Permute.Loss.3UTR, Loss.ROC.3UTR, 
                                  by.x= "Factor", by.y= "data")

# now get significance cutoffs: 
Random.Permute.Loss.3UTR2$Significance<-NA
for (w in 1:length(Random.Permute.Loss.3UTR2$Factor)){
  cutoffs<-quantile(Permute_N_Loss.3UTR, c(.95, .995, .9995))
  if (Random.Permute.Loss.3UTR2$ROCauc[w]> cutoffs[3]){
    Random.Permute.Loss.3UTR2$Significance[w]<-"***"
  } else if (Random.Permute.Loss.3UTR2$ROCauc[w]> cutoffs[2]){
    Random.Permute.Loss.3UTR2$Significance[w]<-"**"
  } else if (Random.Permute.Loss.3UTR2$ROCauc[w]> cutoffs[1]){
    Random.Permute.Loss.3UTR2$Significance[w]<-"*"
  } else {
    Random.Permute.Loss.3UTR2$Significance[w]<-"NS"
  }
}

write.csv(Random.Permute.Loss.3UTR2, 
          file="/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/ROC.Factors.RandomPermutationTest.ChrmLoss.3UTR.csv")
# Where "significant is TRUE/FALSE value for whether or not correlation coefficient p-value is <0.05
# while significance is the threshold for random permutation analtsis with 10k random permutations
# and whtehr the ROC value is >5% of random ROC permutations (*), >0.5% of random permutations (**), or >0.05% of random permutations
# I also manually updated list for overview of factor significance. 

##### Random Permutation test : Ovarian #####
# 211124
# for Ovarian TCGA data
# We want to do random permutations, >1000 permutations per buffering factor
# then we test what the ROC values are for 1000x permutations per factor
# then we see if the actual factor ROC is significantly greater than 95% of random permutation. 
# which ones are greater than random? 

# All.Factors.Diff.Ovarian2

NameTestCollumns= c("PPI.Freq", "gene.Freq", "regulatory.Freq", "Sumoylation.Freq", 
                    "ubiquitination.Freq", "Phosphorylation.Freq", "Methylation.Freq", "Acetylation.Freq", 
                    "AmplificationRatio", "log2.CoeffVar.RNA", "log2.Protein.CoeffVar", "Rate_3", 
                    "m", "p", "lp", "lm", "bm", "bp", 
                    "Mean_HalfLife", "Loops.in.Protein.ANCHOR.", "Intrinsic_Disorder_MobiBD", "Intrinsic_Disorder_LowComplexity", 
                    "Intrinsic_Disorder_Polarity", "Intrinsic_Disorder_Polyampholyte", "Intrinsic_Disorder_homology", 
                    "Dependency.Score", "Aggregation.score", "Non.exponential.decay.delta", "Mutations.Freq", "mutations.FrameShift.Freq", "mutations.nonsense.missense.Freq")

testdatanames= c("Protein-protein interaction", "Protein complex (CORUM)", "Protein regulatory sites", "Sumoylation sites", 
                 "Ubiquitination sites", "Phosphorylation sites", "Methylation sites", "Acetylation sites", 
                 "Percent gene amplification", "RNA neutral variance", "Protein neutral variance", "mRNA decay rate", 
                 "mRNA abundance", "Protein abundance", "Protein length", "mRNA length", "Transcription rate", "Translation rate", 
                 "Protein half life", "Loops in protein score", "Intrinsic protein disorder", "Low complexity score", 
                 "Protein polarity", "Protein polyampholyte score", "Homology score", 
                 "Dependency score", "Aggregation score", "Non-exponential decay delta", "Mutation count (all)", "Mutation count (frame shifts)", "Mutation count (nonsence and missense)")



## Basic Permutation test Function
# where N = number of repeats for test
# where i = one for each factor
N= 10000

PermuteFunction.ChrmGain_Ovarian <- function(y=All.Factors.Diff[,i], x=All.Factors.Diff.Ovarian2$Three.Protein.Gain.x){
  TestDataCol= NameTestCollumns[i]
  name= testdatanames[i]
  TestData<- as.numeric(as.character(dataset %>% pull(TestDataCol)))
  
  
  model.ROC= auc(sample(x=="Buffering", replace=F), 
                 TestData, na.rm=TRUE )
  # permutate category, then runs model
  return(model.ROC)
}
PermuteFunction.ChrmLoss_Ovarian <- function(y=All.Factors.Diff[,i], x=All.Factors.Diff.Ovarian2$Three.Protein.Loss.x){
  TestDataCol= NameTestCollumns[i]
  name= testdatanames[i]
  TestData<- as.numeric(as.character(dataset %>% pull(TestDataCol)))
  
  
  model.ROC= auc(sample(x=="Buffering", replace=F), 
                 TestData, na.rm=TRUE )
  # permutate category, then runs model
  return(model.ROC)
}

#### For chromosome GAIN: 
# to perform the iterations of te permutate function we can use replicate()
# i=1
# Permute_N<- replicate(N, PermuteFunction.ChrmGain_Ovarian())

# now to perform iterations of random permutations for each factor.
Random.Permute.Gain_Ovarian<-data.frame(Factor=as.character(), 
                                Chrm_CN=as.character(),
                                RandomPermutation=c(), 
                                Top5percent=numeric(), 
                                Bottom5percent=numeric())

for (i in 1:length(testdatanames)){
  Permute_N<- replicate(N, PermuteFunction.ChrmGain_Ovarian())
  Random.Permute.Gain_Ovarian<- rbind (Random.Permute.Gain_Ovarian, 
                               data.frame(Factor= testdatanames[i], 
                                          Chrm_CN="Gain",
                                          RandomPermutation= I(list(Permute_N)), 
                                          Top5percent=sort(Permute_N, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                          Bottom5percent=sort(Permute_N)[500]) ) # bottom 500th value is lower 5% cutoff for 10000 values
  
}

# Now see if our ROC value was greater than random, 95% of the time. 

Gain.ROC<- subset(Ovarian.protein.ROCAUC, Chrm_CN=="Gain" & Category=="Buffering")

Random.Permute.Gain_Ovarian2<- merge(Random.Permute.Gain_Ovarian, Gain.ROC, 
                             by.x= "Factor", by.y= "data")

# now get significance cutoffs: 
Random.Permute.Gain_Ovarian2$Significance<-NA
for (w in 1:length(Random.Permute.Gain_Ovarian2$Factor)){
  cutoffs<-quantile(Random.Permute.Gain_Ovarian2$RandomPermutation[[w]], c(.95, .995, .9995))
  if (Random.Permute.Gain_Ovarian2$ROCauc[w]> cutoffs[3]){
    Random.Permute.Gain_Ovarian2$Significance[w]<-"***"
  } else if (Random.Permute.Gain_Ovarian2$ROCauc[w]> cutoffs[2]){
    Random.Permute.Gain_Ovarian2$Significance[w]<-"**"
  } else if (Random.Permute.Gain_Ovarian2$ROCauc[w]> cutoffs[1]){
    Random.Permute.Gain_Ovarian2$Significance[w]<-"*"
  } else {
    Random.Permute.Gain_Ovarian2$Significance[w]<-"NS"
  }
}

Random.Permute.Gain_Ovarian3<-Random.Permute.Gain_Ovarian2[,-3]
write.csv(Random.Permute.Gain_Ovarian3, file="/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/ROC.Factors.RandomPermutationTest.ChrmGain_ovarian.csv")
# Where "significant is TRUE/FALSE value for whether or not correlation coefficient p-value is <0.05
# while significance is the threshold for random permutation analtsis with 10k random permutations
# and whether the ROC value is >5% of random ROC permutations (*), >0.5% of random permutations (**), or >0.05% of random permutations

# Random.Permute.Gain_Ovarian2<-read.csv("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/ROC.Factors.RandomPermutationTest.ChrmGain_ovarian.csv")

## plot Protein-Protein Interaction histogram: 
ggplot(data.frame(Random.Permute.Gain_Ovarian2$RandomPermutation[[26]]), aes(x=Random.Permute.Gain_Ovarian2.RandomPermutation..26..))+
  geom_histogram(bins=50)+
  theme_classic()+
  xlab("Random Permutation ROC AUC for\nProtein-Protein Interactions, chromosome gain")+
  ylab("Frequency")+
  geom_vline(xintercept= Random.Permute.Gain_Ovarian2$ROCauc[26], linetype="dotted")+
  geom_vline(xintercept= 0.5)+
  scale_x_continuous(limits=c(0.45, 0.60))
## plot.RandomPermutation.Factor.ROC.PPI.Gain_Ovarian
# 4x4



### For Chromosome LOSS: 
PermuteFunction.ChrmLoss_Ovarian <- function(y=All.Factors.Diff[,i], x=All.Factors.Diff.Ovarian2$Three.Protein.Loss.x){
  TestDataCol= NameTestCollumns[i]
  name= testdatanames[i]
  TestData<- as.numeric(as.character(dataset %>% pull(TestDataCol)))
  
  
  model.ROC= auc(sample(x=="Buffering", replace=F), 
                 TestData, na.rm=TRUE )
  # permutate category, then runs model
  return(model.ROC)
}

# to perform the iterations of te permutate function we can use replicate()
# i=1
# Permute_N<- replicate(N, PermuteFunction.ChrmLoss_Ovarian())

# now to perform iterations of random permutations for each factor.
Random.Permute.Loss_Ovarian<-data.frame(Factor=as.character(), 
                                Chrm_CN=as.character(),
                                RandomPermutation=c(), 
                                Top5percent=numeric(), 
                                Bottom5percent=numeric())

for (i in 1:length(testdatanames)){
  Permute_N<- replicate(N, PermuteFunction.ChrmLoss_Ovarian())
  Random.Permute.Loss_Ovarian<- rbind (Random.Permute.Loss_Ovarian, 
                               data.frame(Factor= testdatanames[i], 
                                          Chrm_CN="Loss",
                                          RandomPermutation= I(list(Permute_N)), 
                                          Top0.5percent=sort(Permute_N, decreasing=TRUE)[50], # top 500th value is top 5% cutoff for 10000 values
                                          Bottom0.5percent=sort(Permute_N)[50]) ) # bottom 500th value is lower 5% cutoff for 10000 values
  
}

# Now see if our ROC value was greater than random, 95% of the time. 

Loss.ROC<- subset(Ovarian.protein.ROCAUC, Chrm_CN=="Loss" & Category=="Buffering")

Random.Permute.Loss_Ovarian2<- merge(Random.Permute.Loss_Ovarian, Loss.ROC, 
                             by.x= "Factor", by.y= "data")

# now get significance cutoffs: 
Random.Permute.Loss_Ovarian2$Significance<-NA
for (w in 1:length(Random.Permute.Loss_Ovarian2$Factor)){
  cutoffs<-quantile(Random.Permute.Loss_Ovarian2$RandomPermutation[[w]], c(.95, .995, .9995))
  if (Random.Permute.Loss_Ovarian2$ROCauc[w]> cutoffs[3]){
    Random.Permute.Loss_Ovarian2$Significance[w]<-"***"
  } else if (Random.Permute.Loss_Ovarian2$ROCauc[w]> cutoffs[2]){
    Random.Permute.Loss_Ovarian2$Significance[w]<-"**"
  } else if (Random.Permute.Loss_Ovarian2$ROCauc[w]> cutoffs[1]){
    Random.Permute.Loss_Ovarian2$Significance[w]<-"*"
  } else {
    Random.Permute.Loss_Ovarian2$Significance[w]<-"NS"
  }
}
Random.Permute.Loss_Ovarian3<- Random.Permute.Loss_Ovarian2[-3,]

write.csv(Random.Permute.Loss_Ovarian3, 
          file="/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/ROC.Factors.RandomPermutationTest.ChrmLoss_Ovarian.csv")
# Where "significant is TRUE/FALSE value for whether or not correlation coefficient p-value is <0.05
# while significance is the threshold for random permutation analtsis with 10k random permutations
# and whether the ROC value is >5% of random ROC permutations (*), >0.5% of random permutations (**), or >0.05% of random permutations
#Random.Permute.Loss3<-read.csv("/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/ROC.Factors.RandomPermutationTest.ChrmLoss_Ovarian.csv")

## plot Protein-Protein Interaction histogram: 
ggplot(data.frame(Random.Permute.Loss_Ovarian3$RandomPermutation[[26]]), aes(x=Random.Permute.Loss_Ovarian3.RandomPermutation..26..))+
  geom_histogram(bins=50)+
  theme_classic()+
  xlab("Random Permutation ROC AUC for\nProtein-Protein Interactions, chromosome loss")+
  ylab("Frequency")+
  geom_vline(xintercept= Random.Permute.Loss_Ovarian3$ROCauc[26], linetype="dotted")+
  geom_vline(xintercept= 0.5)+
  scale_x_continuous(limits=c(0.45, 0.60))
## plot.RandomPermutation.Factor.ROC.PPI.Loss
# 4x4



##### Random Permutation test for the 3' and 5' UTR regions 
# 211124
# I also want to test 5' UTR and 3' UTR data, but there are multiple 5' and 3' datapoints per gene, 
# so I do that seperately so as not to mess with other data values/frequencies. 

### 5' UTR data random permutation test ###
# Step one: get 5' UTR data (from above)

dataset= All.Factors.Diff.Ovarian2

## add 5' UTR data 
TestData<- NCBI.genedata.NM.5UTR
name<- "5' UTR length" 

Diff.Factor.5UTR<-merge(x=TestData, y= All.Factors.Diff, 
                        by.x="name2", by.y="Gene_Symbol", 
                        sort = TRUE)


## Write the basic Permutation test Function
# where N = number of repeats for test
# where i = one for each factor
N= 10000

PermuteFunction.ChrmGain.5UTR <- function(y=Diff.Factor.5UTR$UTR5, x=Diff.Factor.5UTR$Three.Protein.Gain){
  TestDataCol= "UTR5"
  name= "5' UTR"
  TestData<- as.numeric(as.character(Diff.Factor.5UTR %>% pull(TestDataCol)))
  
  
  model.ROC= auc(sample(x=="Buffering", replace=F), 
                 TestData, na.rm=TRUE ) 
  # Permutate category, then runs model 
  return(model.ROC)
}
PermuteFunction.ChrmLoss.5UTR <- function(y=Diff.Factor.5UTR$UTR5, x=Diff.Factor.5UTR$Three.Protein.Loss){
  TestDataCol= "UTR5"
  name= "5' UTR"
  TestData<- as.numeric(as.character(Diff.Factor.5UTR %>% pull(TestDataCol)))
  
  
  model.ROC= auc(sample(x=="Buffering", replace=F), 
                 TestData, na.rm=TRUE )
  # permutate category, then runs model
  return(model.ROC)
}

#### For 5' UTR chromosome GAIN: 
# to perform the iterations of te permutate function we can use replicate()
# i=1
# Permute_N<- replicate(N, PermuteFunction())

Permute_N.5UTRGain<- replicate(N, PermuteFunction.ChrmGain.5UTR())
write.csv(Permute_N.5UTRGain, "RandomPermute.ROC.ChrmGain.Buffering.5UTR_Ovarian.csv")
Random.Permute.Gain.5UTR<- data.frame(Factor= "5' UTR length", 
                                      Chrm_CN="Gain",
                                      Top5percent=sort(Permute_N.5UTRGain, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                      Bottom5percent=sort(Permute_N.5UTRGain)[500])  # bottom 500th value is lower 5% cutoff for 10000 values


# Now see if our ROC value was greater than random, 95% of the time. 
# setwd()
Ovarian.protein.ROCAUC<- read.csv(file="Factor.ROCAUC.correlationscore.csv")

Gain.ROC.5UTR_Ovarian<- subset(Ovarian.protein.ROCAUC, Chrm_CN=="Gain" & Category=="Buffering")

Random.Permute.Gain.5UTR2<- merge(Random.Permute.Gain.5UTR, Gain.ROC.5UTR_Ovarian, 
                                  by.x= "Factor", by.y= "data")

# now get significance cutoffs: 
Random.Permute.Gain.5UTR2$Significance<-NA
for (w in 1:length(Random.Permute.Gain.5UTR2$Factor)){
  cutoffs<-quantile(Permute_N.5UTRGain, c(.95, .995, .9995))
  if (Random.Permute.Gain.5UTR2$ROCauc[w]> cutoffs[3]){
    Random.Permute.Gain.5UTR2$Significance[w]<-"***"
  } else if (Random.Permute.Gain.5UTR2$ROCauc[w]> cutoffs[2]){
    Random.Permute.Gain.5UTR2$Significance[w]<-"**"
  } else if (Random.Permute.Gain.5UTR2$ROCauc[w]> cutoffs[1]){
    Random.Permute.Gain.5UTR2$Significance[w]<-"*"
  } else {
    Random.Permute.Gain.5UTR2$Significance[w]<-"NS"
  }
}
# setwd()
write.csv(Random.Permute.Gain.5UTR2, file="/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/ROC.Factors.RandomPermutationTest.ChrmGain.5UTR_Ovarian.csv")
# Where "significant is TRUE/FALSE value for whether or not correlation coefficient p-value is <0.05
# while significance is the threshold for random permutation analtsis with 10k random permutations
# and whether the ROC value is >5% of random ROC permutations (*), >0.5% of random permutations (**), or >0.05% of random permutations
# I updated the overview file with data for each factor's permutation analysis. 



### For 5' UTR Chromosome LOSS: 
Permute_N_Loss.5UTR<- replicate(N, PermuteFunction.ChrmLoss.5UTR())
write.csv(Permute_N_Loss.5UTR, "RandomPermute.ROC.ChrmLoss.Buffering.5UTR_Ovarian.csv")
Random.Permute.Loss.5UTR<- data.frame(Factor= "5' UTR length", 
                                      Chrm_CN="Loss",
                                      Top0.5percent=sort(Permute_N_Loss.5UTR, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                      Bottom0.5percent=sort(Permute_N_Loss.5UTR)[500]) # bottom 500th value is lower 5% cutoff for 10000 values


# Now see if our ROC value was greater than random, 95% of the time. 

Loss.ROC.5UTR_Ovarian<- subset(Ovarian.protein.ROCAUC, Chrm_CN=="Loss" & Category=="Buffering")

Random.Permute.Loss.5UTR2<- merge(Random.Permute.Loss.5UTR, Loss.ROC.5UTR_Ovarian, 
                                  by.x= "Factor", by.y= "data")

# now get significance cutoffs: 
Random.Permute.Loss.5UTR2$Significance<-NA
for (w in 1:length(Random.Permute.Loss.5UTR2$Factor)){
  cutoffs<-quantile(Permute_N_Loss.5UTR, c(.95, .995, .9995))
  if (Random.Permute.Loss.5UTR2$ROCauc[w]> cutoffs[3]){
    Random.Permute.Loss.5UTR2$Significance[w]<-"***"
  } else if (Random.Permute.Loss.5UTR2$ROCauc[w]> cutoffs[2]){
    Random.Permute.Loss.5UTR2$Significance[w]<-"**"
  } else if (Random.Permute.Loss.5UTR2$ROCauc[w]> cutoffs[1]){
    Random.Permute.Loss.5UTR2$Significance[w]<-"*"
  } else {
    Random.Permute.Loss.5UTR2$Significance[w]<-"NS"
  }
}

write.csv(Random.Permute.Loss.5UTR2, 
          file="/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/ROC.Factors.RandomPermutationTest.ChrmLoss.5UTR_Ovarian.csv")
# Where "significant is TRUE/FALSE value for whether or not correlation coefficient p-value is <0.05
# while significance is the threshold for random permutation analtsis with 10k random permutations
# and whtehr the ROC value is >5% of random ROC permutations (*), >0.5% of random permutations (**), or >0.05% of random permutations
# I also manually updated list for overview of factor significance. 




### 3' UTR data random permutation test ###
# get 3' UTR data (from above)

dataset= All.Factors.Diff.Ovarian2

## add 3' UTR data
TestData<- NCBI.genedata.NM.3UTR
name<- "3' UTR length" 

Diff.Factor.3UTR<-merge(x=TestData, y= All.Factors.Diff.Ovarian2, 
                        by.x="name2", by.y="Gene", 
                        sort = TRUE)


## Write the basic Permutation test Function
# where N = number of repeats for test
# where i = one for each factor
N= 10000

PermuteFunction.ChrmGain.3UTR <- function(y=Diff.Factor.3UTR$UTR3, x=Diff.Factor.3UTR$Three.Protein.Gain.x){
  TestDataCol= "UTR3"
  name= "3' UTR"
  TestData<- as.numeric(as.character(Diff.Factor.3UTR %>% pull(TestDataCol)))
  
  
  model.ROC= auc(sample(x=="Buffering", replace=F), 
                 TestData, na.rm=TRUE )
  # permutate category, then runs model
  return(model.ROC)
}
PermuteFunction.ChrmLoss.3UTR <- function(y=Diff.Factor.3UTR$UTR3, x=Diff.Factor.3UTR$Three.Protein.Loss.x){
  TestDataCol= "UTR3"
  name= "3' UTR"
  TestData<- as.numeric(as.character(Diff.Factor.3UTR %>% pull(TestDataCol)))
  
  
  model.ROC= auc(sample(x=="Buffering", replace=F), 
                 TestData, na.rm=TRUE )
  # permutate category, then runs model
  return(model.ROC)
}

#### For 3' UTR chromosome GAIN: 
# to perform the iterations of te permutate function we can use replicate()
# i=1
# Permute_N<- replicate(N, PermuteFunction())

Permute_N.3UTRGain<- replicate(N, PermuteFunction.ChrmGain.3UTR())
write.csv(Permute_N.3UTRGain, "RandomPermute.ROC.ChrmGain.Buffering.3UTR_Ovarian.csv")
Random.Permute.Gain.3UTR<- data.frame(Factor= "3' UTR length", 
                                      Chrm_CN="Gain",
                                      Top5percent=sort(Permute_N.3UTRGain, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                      Bottom5percent=sort(Permute_N.3UTRGain)[500])  # bottom 500th value is lower 5% cutoff for 10000 values


# Now see if our ROC value was greater than random, 95% of the time. 
# setwd()
#Corr.protein.ROCAUC<- read.csv(file="Factor.ROCAUC.correlationscore.csv")

Gain.ROC.3UTR<- subset(Ovarian.protein.ROCAUC, Chrm_CN=="Gain" & Category=="Buffering")

Random.Permute.Gain.3UTR2<- merge(Random.Permute.Gain.3UTR, Gain.ROC.3UTR, 
                                  by.x= "Factor", by.y= "data")

# now get significance cutoffs: 
Random.Permute.Gain.3UTR2$Significance<-NA
for (w in 1:length(Random.Permute.Gain.3UTR2$Factor)){
  cutoffs<-quantile(Permute_N.3UTRGain, c(.95, .995, .9995))
  if (Random.Permute.Gain.3UTR2$ROCauc[w]> cutoffs[3]){
    Random.Permute.Gain.3UTR2$Significance[w]<-"***"
  } else if (Random.Permute.Gain.3UTR2$ROCauc[w]> cutoffs[2]){
    Random.Permute.Gain.3UTR2$Significance[w]<-"**"
  } else if (Random.Permute.Gain.3UTR2$ROCauc[w]> cutoffs[1]){
    Random.Permute.Gain.3UTR2$Significance[w]<-"*"
  } else {
    Random.Permute.Gain.3UTR2$Significance[w]<-"NS"
  }
}
# setwd()
write.csv(Random.Permute.Gain.3UTR2, file="/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/ROC.Factors.RandomPermutationTest.ChrmGain.3UTR_Ovarian.csv")
# Where "significant is TRUE/FALSE value for whether or not correlation coefficient p-value is <0.05
# while significance is the threshold for random permutation analtsis with 10k random permutations
# and whtehr the ROC value is >5% of random ROC permutations (*), >0.5% of random permutations (**), or >0.05% of random permutations
# I updated the overview file with data for each factor's permutation analysis. 


### For 3' UTR Chromosome LOSS: 
Permute_N_Loss.3UTR<- replicate(N, PermuteFunction.ChrmLoss.3UTR())
write.csv(Permute_N_Loss.3UTR, "RandomPermute.ROC.ChrmLoss.Buffering.3UTR_ovarian.csv")
Random.Permute.Loss.3UTR<- data.frame(Factor= "3' UTR length", 
                                      Chrm_CN="Loss",
                                      Top0.5percent=sort(Permute_N_Loss.3UTR, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                      Bottom0.5percent=sort(Permute_N_Loss.3UTR)[500]) # bottom 500th value is lower 5% cutoff for 10000 values


# Now see if our ROC value was greater than random, 95% of the time. 

Loss.ROC.3UTR<- subset(Ovarian.protein.ROCAUC, Chrm_CN=="Loss" & Category=="Buffering")

Random.Permute.Loss.3UTR2<- merge(Random.Permute.Loss.3UTR, Loss.ROC.3UTR, 
                                  by.x= "Factor", by.y= "data")

# now get significance cutoffs: 
Random.Permute.Loss.3UTR2$Significance<-NA
for (w in 1:length(Random.Permute.Loss.3UTR2$Factor)){
  cutoffs<-quantile(Permute_N_Loss.3UTR, c(.95, .995, .9995))
  if (Random.Permute.Loss.3UTR2$ROCauc[w]> cutoffs[3]){
    Random.Permute.Loss.3UTR2$Significance[w]<-"***"
  } else if (Random.Permute.Loss.3UTR2$ROCauc[w]> cutoffs[2]){
    Random.Permute.Loss.3UTR2$Significance[w]<-"**"
  } else if (Random.Permute.Loss.3UTR2$ROCauc[w]> cutoffs[1]){
    Random.Permute.Loss.3UTR2$Significance[w]<-"*"
  } else {
    Random.Permute.Loss.3UTR2$Significance[w]<-"NS"
  }
}

write.csv(Random.Permute.Loss.3UTR2, 
          file="/Users/user/Documents/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/ROC.Factors.RandomPermutationTest.ChrmLoss.3UTR.csv")
# Where "significant is TRUE/FALSE value for whether or not correlation coefficient p-value is <0.05
# while significance is the threshold for random permutation analtsis with 10k random permutations
# and whtehr the ROC value is >5% of random ROC permutations (*), >0.5% of random permutations (**), or >0.05% of random permutations
# I also manually updated list for overview of factor significance. 


##### Random Permutation test : 5D- buffering factors boxplot #####
# 210908
# We want to do random permutations, >1000 permutations per buffering factor boxplot
# to see whether the t-test p-values are maintained. 
# which ones are greater than random? 

All.Factors.Diff

NameTestCollumns= c("PPI.Freq", "gene.Freq", "regulatory.Freq", "Sumoylation.Freq", 
                    "ubiquitination.Freq", "Phosphorylation.Freq", "Methylation.Freq", "Acetylation.Freq", 
                    "AmplificationRatio", "log2.CoeffVar.RNA", "log2.Protein.CoeffVar", "Rate_3", 
                    "m", "p", "lp", "lm", "bm", "bp", 
                    "Mean_HalfLife", "Loops.in.Protein.ANCHOR.", "Intrinsic_Disorder_MobiBD", "Intrinsic_Disorder_LowComplexity", 
                    "Intrinsic_Disorder_Polarity", "Intrinsic_Disorder_Polyampholyte", "Intrinsic_Disorder_homology", 
                    "Dependency.Score", "Aggregation.score", "Non.exponential.decay.delta", "Mutations.Freq", "mutations.FrameShift.Freq", "mutations.nonsense.missense.Freq")

testdatanames= c("Protein-protein interaction", "Protein complex (CORUM)", "Protein regulatory sites", "Sumoylation sites", 
                 "Ubiquitination sites", "Phosphorylation sites", "Methylation sites", "Acetylation sites", 
                 "Percent gene amplification", "RNA neutral variance", "Protein neutral variance", "mRNA decay rate", 
                 "mRNA abundance", "Protein abundance", "Protein length", "mRNA length", "Transcription rate", "Translation rate", 
                 "Protein half life", "Loops in protein score", "Intrinsic protein disorder", "Low complexity score", 
                 "Protein polarity", "Protein polyampholyte score", "Homology score", 
                 "Dependency score", "Aggregation score", "Non-exponential decay delta", "Mutation count (all)", "Mutation count (frame shifts)", "Mutation count (nonsence and missense)")

dataset= All.Factors.Diff

# All.Factors.Diff, aes(x=Three.Protein.Gain, y=testcollumn)


## Basic Permutation test Function
# where N = number of repeats for test
# where i = one for each factor
N= 10000

PermuteFunction.ChrmGain.Buff.Scale <- function(y=All.Factors.Diff[,i], x=All.Factors.Diff$Three.Protein.Gain){
  TestDataCol= NameTestCollumns[i]
  name= testdatanames[i]
  TestData<- as.numeric(as.character(dataset %>% pull(TestDataCol)))
  
  
  Random.ttest= t.test(TestData[sample(x=="Buffering", replace=F)], 
                       TestData[sample(x=="Scaling", replace=F)])
  # permutate category, then runs model
  return(Random.ttest$p.value)
}
PermuteFunction.ChrmGain.Buff.AS <- function(y=All.Factors.Diff[,i], x=All.Factors.Diff$Three.Protein.Gain){
  TestDataCol= NameTestCollumns[i]
  name= testdatanames[i]
  TestData<- as.numeric(as.character(dataset %>% pull(TestDataCol)))
  
  
  Random.ttest= t.test(TestData[sample(x=="Buffering", replace=F)], 
                       TestData[sample(x=="Anti-Scaling", replace=F)])
  # permutate category, then runs model
  return(Random.ttest$p.value)
}
PermuteFunction.ChrmLoss.Buff.Scale <- function(y=All.Factors.Diff[,i], x=All.Factors.Diff$Three.Protein.Loss){
  TestDataCol= NameTestCollumns[i]
  name= testdatanames[i]
  TestData<- as.numeric(as.character(dataset %>% pull(TestDataCol)))
  
  
  Random.ttest= t.test(TestData[sample(x=="Buffering", replace=F)], 
                       TestData[sample(x=="Scaling", replace=F)])
  # permutate category, then runs model
  return(Random.ttest$p.value)
}
PermuteFunction.ChrmLoss.Buff.AS <- function(y=All.Factors.Diff[,i], x=All.Factors.Diff$Three.Protein.Loss){
  TestDataCol= NameTestCollumns[i]
  name= testdatanames[i]
  TestData<- as.numeric(as.character(dataset %>% pull(TestDataCol)))
  
  
  Random.ttest= t.test(TestData[sample(x=="Buffering", replace=F)], 
                       TestData[sample(x=="Anti-Scaling", replace=F)])
  # permutate category, then runs model
  return(Random.ttest$p.value)
}

#### For chromosome GAIN: 
# to perform the iterations of te permutate function we can use replicate()
# i=1
# Permute_N<- replicate(N, PermuteFunction())

# Fig5Dtest is list of i values corresponding to boxplot terms in fig 5D 
Fig5Dtest<-c(28, 10, 26, 5,1,2)

# now to perform iterations of random permutations for each factor.
Random.Permute.Gain.category<-data.frame(Factor=as.character(), 
                                Chrm_CN=as.character(),
                                ttest.pvalues.Buff.Scale=c(), 
                                ttest.pvalues.Buff.AS=c(), 
                                Top5percent.Scale=numeric(), 
                                Bottom5percent.Scale=numeric(),
                                Top5percent.AS=numeric(), 
                                Bottom5percent.AS=numeric())

for (i in c(28, 10, 26, 5,1,2) ){
  Permute_N_Scale<- replicate(N, PermuteFunction.ChrmGain.Buff.Scale())
  Permute_N_AS<- replicate(N, PermuteFunction.ChrmGain.Buff.AS())
  
  Random.Permute.Gain.category<- rbind (Random.Permute.Gain.category, 
                               data.frame(Factor= testdatanames[i], 
                                          Chrm_CN="Gain",
                                          ttest.pvalues.Buff.Scale= I(list(Permute_N_Scale)), 
                                          ttest.pvalues.Buff.AS= I(list(Permute_N_AS)), 
                                          Top5percent.Scale=sort(Permute_N_Scale, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                          Bottom5percent.Scale=sort(Permute_N_Scale)[500],  
                                          Bottom.Scale=min(Permute_N_Scale), 
                               Top5percent.AS=sort(Permute_N_AS, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                               Bottom5percent.AS=sort(Permute_N_AS)[500]), 
                               Bottom.AS=min(Permute_N_AS)
                               ) # bottom 500th value is lower 5% cutoff for 10000 values
  
}

# Now see if our t-test p-value  was less than random, 95% of the time. 
# look at p-values for 
x<-c()
y<-c()
for (i in 1:6){
  x<- append(x, min(Random.Permute.Gain.category$ttest.pvalues.Buff.Scale[[i]]) )
  y<- append(y, min(Random.Permute.Gain.category$ttest.pvalues.Buff.AS[[i]]) )
}

Random.Permute.Gain.category$Bottom.Scale<-x
Random.Permute.Gain.category$Bottom.AS<- y

write.csv(Random.Permute.Gain.category, file="/Volumes/Schukken_SSD/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/Factors.RandomPermutationTest.Category.ChrmGain.csv")
## for chromosome gain all Buffer vs scaing and all buffered vs Anti-Scaling data
## had their actual t-test between categories (Figure 5D) per factor lower than the 5% random t-test cutoff
## showing that these p-values were more significant than random chance in this dataset. 

#                           Factor Chrm_CN ttest.pvalues.Buff.Scale ttest.pvalues.Buff.AS
#1     Non-exponential decay delta    Gain             0.887486....          0.447465....
#2            RNA neutral variance    Gain             0.293626....          0.527408....
#3                Dependency score    Gain             0.561493....          0.943156....
#4            Ubiquitination sites    Gain             0.218870....          0.393000....
#5     Protein-protein interaction    Gain             0.626907....          0.942442....
#6         Protein complex (CORUM)    Gain             0.568433....          0.483447....
#7 Ovarian_Protein complex (CORUM)    Gain             0.306814....          0.264448....
#8        Ovarian_Dependency score    Gain             0.709295....          0.271740....
#Top5percent.Scale Bottom5percent.Scale Top5percent.AS Bottom5percent.AS Botton.Scale    Botton.AS
#1         0.9606339           0.13052747      0.9570715        0.07308991 1.407699e-03 3.214599e-05
#2         0.9601355           0.13694001      0.9550431        0.07191588 3.797800e-03 1.099387e-04
#3         0.9645035           0.13637815      0.9529093        0.07129145 1.475204e-03 9.404772e-05
#4         0.9614003           0.13048316      0.9546440        0.07475670 9.932620e-04 8.311637e-07
#5         0.9590462           0.13157108      0.9533592        0.06424552 3.828220e-04 2.772443e-05
#6         0.9598603           0.13037247      0.9525739        0.06497395 2.184883e-04 7.606768e-06
#7         0.9505399           0.06271098      0.9459577        0.04819380 7.653789e-05 4.899720e-05
#8         0.9543176           0.06329547      0.9514679        0.05313208 3.526028e-06 1.566212e-05


#### For chromosome LOSS: 
# to perform the iterations of te permutate function we can use replicate()
# i=1
# Permute_N<- replicate(N, PermuteFunction())

# Fig5Dtest is list of i values corresponding to boxplot terms in fig 5D 
Fig5Dtest<-c(28, 10, 26, 5,1,2)

# now to perform iterations of random permutations for each factor.
Random.Permute.Loss.category<-data.frame(Factor=as.character(), 
                                         Chrm_CN=as.character(),
                                         ttest.pvalues.Buff.Scale=c(), 
                                         ttest.pvalues.Buff.AS=c(), 
                                         Top5percent.Scale=numeric(), 
                                         Bottom5percent.Scale=numeric(),
                                         Top5percent.AS=numeric(), 
                                         Bottom5percent.AS=numeric())

for (i in c(28, 10, 26, 5,1,2) ){
  Permute_N_Scale<- replicate(N, PermuteFunction.ChrmLoss.Buff.Scale())
  Permute_N_AS<- replicate(N, PermuteFunction.ChrmLoss.Buff.AS())
  
  Random.Permute.Loss.category<- rbind (Random.Permute.Loss.category, 
                                        data.frame(Factor= testdatanames[i], 
                                                   Chrm_CN="Loss",
                                                   ttest.pvalues.Buff.Scale= I(list(Permute_N_Scale)), 
                                                   ttest.pvalues.Buff.AS= I(list(Permute_N_AS)), 
                                                   Top5percent.Scale=sort(Permute_N_Scale, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                                   Bottom5percent.Scale=sort(Permute_N_Scale)[500],  
                                                   Bottom.Scale=min(Permute_N_Scale), 
                                                   Top5percent.AS=sort(Permute_N_AS, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                                   Bottom5percent.AS=sort(Permute_N_AS)[500], 
                                                   Bottom.AS=min(Permute_N_AS))
  ) # bottom 500th value is lower 5% cutoff for 10,000 values
  
}

# look at p-values for each categorical 
x<-c()
y<-c()
for (i in 1:6){
  x<- append(x, min(Random.Permute.Loss.category$ttest.pvalues.Buff.Scale[[i]]) )
  y<- append(y, min(Random.Permute.Loss.category$ttest.pvalues.Buff.AS[[i]]) )
}

Random.Permute.Loss.category$Botton.Scale<-x
Random.Permute.Loss.category$Botton.AS<- y


write.csv(Random.Permute.Loss.category, file="/Volumes/Schukken_SSD/Depmap Aneuploidy/Protein.Quantile/RNA_Protein_comparison/Buffering_factors/RandomPermutation/Factors.RandomPermutationTest.Category.ChrmLoss.csv")
# Now see if our t-test p-value  was less than random, 95% of the time. by comparing to actual factor categorical cutoffs
## for chromosome loss all Buffer vs scaing and all buffered vs Anti-Scaling data
## had their actual t-test between categories (Figure 5D) per factor lower than the 5% random t-test cutoff
## showing that these p-values were more significant than random chance in this dataset. 

#                           Factor Chrm_CN ttest.pvalues.Buff.Scale ttest.pvalues.Buff.AS
#1     Non-exponential decay delta    Gain             0.740107....          0.063469....
#2            RNA neutral variance    Gain             0.833588....          0.592496....
#3                Dependency score    Gain             0.858945....          0.933694....
#4            Ubiquitination sites    Gain             0.099239....          0.925212....
#5     Protein-protein interaction    Gain             0.387154....          0.586649....
#6         Protein complex (CORUM)    Gain             0.325300....          0.611650....
#7 Ovarian_Protein complex (CORUM)    Gain             0.879048....          0.335962....
#8        Ovarian_Dependency score    Gain             0.869087....          0.990797....

#Top5percent.Scale Bottom5percent.Scale Top5percent.AS Bottom5percent.AS Botton.Scale    Botton.AS
#1         0.9586612           0.10606949      0.9563292        0.07187615 2.498501e-03 4.833214e-05
#2         0.9575406           0.11634452      0.9553399        0.07780728 3.555554e-03 2.315410e-04
#3         0.9595724           0.12396965      0.9568119        0.07100033 2.991096e-03 1.101464e-04
#4         0.9601459           0.10933270      0.9539134        0.07308793 1.287353e-03 2.233219e-05
#5         0.9613373           0.11468925      0.9577826        0.06594083 8.430026e-04 6.279277e-06
#6         0.9567850           0.11182942      0.9521905        0.06326073 5.660182e-04 1.421211e-06
#7         0.9459275           0.06042190      0.9514083        0.05201056 9.775394e-06 5.940676e-06
#8         0.9562617           0.06205504      0.9488484        0.05299193 8.692719e-06 1.309959e-05

#### Ovarian
All.Factors.Diff.Ovarian2
#test  protein complex (#2) and dependency score (# 26)
N=10000
#test ovarian CORUM boxplot data
# make 10,000 random permutations and do t-test with those random permutations
# see if our data is less than the bottom 5% of random tests
# see if actual data is less than the lowest of the 10,000 random tests

## chrm gain: 
## test Protein complexes: 
Permute_N_Scale<- replicate(N, PermuteFunction.ChrmGain.Buff.Scale(y=All.Factors.Diff.Ovarian2[,2], 
                                                                   x=All.Factors.Diff.Ovarian2$Three.Protein.Gain.x
                                                                   ))
## test Protein complexes AS: 
Permute_N_AS<- replicate(N, PermuteFunction.ChrmGain.Buff.AS(y=All.Factors.Diff.Ovarian2[,2], 
                                                             x=All.Factors.Diff.Ovarian2$Three.Protein.Gain.x
))

Random.Permute.Gain.category<- rbind (Random.Permute.Gain.category, 
                                      data.frame(Factor= paste0("Ovarian_",testdatanames[2]), 
                                                 Chrm_CN="Gain",
                                                 ttest.pvalues.Buff.Scale= I(list(Permute_N_Scale)), 
                                                 ttest.pvalues.Buff.AS= I(list(Permute_N_AS)), 
                                                 Top5percent.Scale=sort(Permute_N_Scale, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                                 Bottom5percent.Scale=sort(Permute_N_Scale)[500],  
                                                 Botton.Scale=min(Permute_N_Scale), 
                                                 Top5percent.AS=sort(Permute_N_AS, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                                 Bottom5percent.AS=sort(Permute_N_AS)[500], 
                                                 Botton.AS=min(Permute_N_AS) ) ) 

### Now test dependency ovarian TCGA data
Permute_N_Scale<- replicate(N, PermuteFunction.ChrmGain.Buff.Scale(y=All.Factors.Diff.Ovarian2[,26], 
                                                                   x=All.Factors.Diff.Ovarian2$Three.Protein.Gain.x
))
## test Protein complexes AS: 
Permute_N_AS<- replicate(N, PermuteFunction.ChrmGain.Buff.AS(y=All.Factors.Diff.Ovarian2[,26], 
                                                             x=All.Factors.Diff.Ovarian2$Three.Protein.Gain.x
))

Random.Permute.Gain.category<- rbind (Random.Permute.Gain.category, 
                                      data.frame(Factor= paste0("Ovarian_",testdatanames[26]), 
                                                 Chrm_CN="Gain",
                                                 ttest.pvalues.Buff.Scale= I(list(Permute_N_Scale)), 
                                                 ttest.pvalues.Buff.AS= I(list(Permute_N_AS)), 
                                                 Top5percent.Scale=sort(Permute_N_Scale, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                                 Bottom5percent.Scale=sort(Permute_N_Scale)[500],  
                                                 Botton.Scale=min(Permute_N_Scale), 
                                                 Top5percent.AS=sort(Permute_N_AS, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                                 Bottom5percent.AS=sort(Permute_N_AS)[500], 
                                                 Botton.AS=min(Permute_N_AS) ) ) 


# chrm loss Ovarian
## test Protein complexes: 
Permute_N_Scale<- replicate(N, PermuteFunction.ChrmLoss.Buff.Scale(y=All.Factors.Diff.Ovarian2[,2], 
                                                                   x=All.Factors.Diff.Ovarian2$Three.Protein.Gain.x
))
Permute_N_AS<- replicate(N, PermuteFunction.ChrmLoss.Buff.AS(y=All.Factors.Diff.Ovarian2[,2], 
                                                             x=All.Factors.Diff.Ovarian2$Three.Protein.Gain.x
))

Random.Permute.Loss.category<- rbind (Random.Permute.Loss.category, 
                                      data.frame(Factor= paste0("Ovarian_",testdatanames[2]), 
                                                 Chrm_CN="Gain",
                                                 ttest.pvalues.Buff.Scale= I(list(Permute_N_Scale)), 
                                                 ttest.pvalues.Buff.AS= I(list(Permute_N_AS)), 
                                                 Top5percent.Scale=sort(Permute_N_Scale, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                                 Bottom5percent.Scale=sort(Permute_N_Scale)[500],  
                                                 Botton.Scale=min(Permute_N_Scale), 
                                                 Top5percent.AS=sort(Permute_N_AS, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                                 Bottom5percent.AS=sort(Permute_N_AS)[500], 
                                                 Botton.AS=min(Permute_N_AS) ) ) 

### Now test dependency ovarian TCGA data
Permute_N_Scale<- replicate(N, PermuteFunction.ChrmLoss.Buff.Scale(y=All.Factors.Diff.Ovarian2[,26], 
                                                                   x=All.Factors.Diff.Ovarian2$Three.Protein.Gain.x
))

Permute_N_AS<- replicate(N, PermuteFunction.ChrmLoss.Buff.AS(y=All.Factors.Diff.Ovarian2[,26], 
                                                             x=All.Factors.Diff.Ovarian2$Three.Protein.Gain.x
))

Random.Permute.Loss.category<- rbind (Random.Permute.Loss.category, 
                                      data.frame(Factor= paste0("Ovarian_",testdatanames[26]), 
                                                 Chrm_CN="Gain",
                                                 ttest.pvalues.Buff.Scale= I(list(Permute_N_Scale)), 
                                                 ttest.pvalues.Buff.AS= I(list(Permute_N_AS)), 
                                                 Top5percent.Scale=sort(Permute_N_Scale, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                                 Bottom5percent.Scale=sort(Permute_N_Scale)[500],  
                                                 Botton.Scale=min(Permute_N_Scale), 
                                                 Top5percent.AS=sort(Permute_N_AS, decreasing=TRUE)[500], # top 500th value is top 5% cutoff for 10000 values
                                                 Bottom5percent.AS=sort(Permute_N_AS)[500], 
                                                 Botton.AS=min(Permute_N_AS) ) ) 




