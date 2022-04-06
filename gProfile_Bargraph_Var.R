###### g:profile Gene enrichemnt analysis and plots ####
## 210303
## g:profile analysis 
## Categorical (Buffer/scaling/antiscaling, down /up-regulated) Protein expression 
## upon chromosome gain or loss. 

## Klaske Schukken

library(ggplot2)
library(tidyverse) 
library(readxl)
library(ggpubr)
library("cowplot")
library(plyr)
library('gprofiler2')


# Notes: RNA/Protein expression difference upon chrm gain/loss (difference, not corr!)
# Difference in expression can be + or -. 
# with filtered cell lines only (only cells with RNA & Protein data)

# To do: 
# per category I am looking at, run samples through g:profiler 
# either as categorical groups or as ordered list. 
# Sub-Step 1: Get g:profiler data for RNA/Protein Chrm Gain/Loss 
# Sub-Step 2: filter for top 5-10 GO terms 
# SubStep 3: bargraph -log2() of p-value for top terms for all conditions. 

###### Step 1: Get data #####

setwd() #set location 
#make in Protein_RNA_expression.PerCell_v2.R
# also available from supplementary data ** sheet 2.
CN.Diff.xRNA.yProt.ThreeGroups<- read.csv("RNA.Protein_Mono.Di.Tri_Difference_Pvalue_min10points_3category.csv")
CN.Diff.xRNA.yProt.ThreeGroups<- CN.Diff.xRNA.yProt.ThreeGroups[,-c(1)]


### These datasets are from g:profiler, the online interface, where I submitted RNA and proteins 
# as ordered lists (highest to lowest, and lowest to highest) based on their correlation with cellular aneuploidy
# then I uploaded these lists to this script. 

## available from Supplementary Data **
RNA.negSig<-read.delim2("gProfiler_RNA.negCor.AnScore.expression.BHSig.csv", 
                        dec=".", header = TRUE, sep=",") 
RNA.posSig<-read.delim2("gProfiler_RNA.posCor.AnScore.expression.BHSig.csv", 
                        dec=".", header = TRUE, sep=",") 

Protein.negSig<-read.delim2("gProfiler_Protein.NegCor.AneuploidyScore.expression.BHSig.csv", 
                            dec=".", header = TRUE, sep=",") 
Protein.posSig<-read.delim2("gProfiler_Protein.PosCorr.AneuploidyScore.BHSig.csv", 
                            dec=".", header = TRUE, sep=",") 



###### Figure 4: gene expression correlates with aneuploidy score #####
## Figure 4: gene expression correlates with aneuploidy score
## Look at RNA and Proteins positively or negatively correlated with aneuploidy score

#make bargraph of key terms/hits RNA
# RNA positively and negatively correlated with aneuploidy
R.KeyNegSig<-RNA.negSig[order(RNA.negSig$adjusted_p_value),]
R.KeyNegSig1<- subset(R.KeyNegSig, source %in%  c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
R.KeyNegSig2<- R.KeyNegSig1[c(1,9,10,12,38, 13,25,26,32,37),]#ribosome (2,4,7,), translation, protease
R.KeyNegSig2$term_name<- factor(R.KeyNegSig2$term_name, levels= R.KeyNegSig2$term_name)
R.KeyNegSig2$Termtype<-c("Ribosome","Ribosome","Ribosome","Ribosome","Ribosome",
                       "RNA processing","RNA processing","RNA processing", "RNA processing","RNA processing")

R.KeyPosSig<-RNA.posSig[order(RNA.posSig$adjusted_p_value),]
R.KeyPosSig1<- subset(R.KeyPosSig, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
R.KeyPosSig2<-R.KeyPosSig1[c(2,8,19,91,174,201,  412,413),]#protein complexes, DNA damage response, heat shock
R.KeyPosSig2$term_name<- factor(R.KeyPosSig2$term_name, levels= R.KeyPosSig2$term_name)
R.KeyPosSig2$Termtype<-c("Adhesion", "Adhesion", "Adhesion", "Adhesion", "Adhesion", "Adhesion", 
                         "Heat acclimation", "Heat acclimation")

#subset(R.KeyPosSig1, grepl( "heat", R.KeyPosSig1$term_name, fixed = TRUE))
# look for heat shock terms.

plot.KeyTermsRNA<- ggplot(R.KeyPosSig2, 
                       aes(x=term_name, y=-log2(adjusted_p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key RNA terms")+
  ggtitle ("Positively correlated with aneuploidy")+
  theme_classic()+
  ylim(0,140)+
  coord_flip()
plot.KeyTermsRNA#7x5
# plot.RNA.Pos.AneuScoreCorr.KeyTerms_3


#make bargraph of key terms/hits PROTEIN positively and negatively regulated with aneuploidy score
KeyNegSig<-Protein.negSig[order(Protein.negSig$adjusted_p_value),]
KeyNegSig2<-KeyNegSig[c(6,8,9,15, 2,24,27),]#ribosome, translation, protease
KeyNegSig2$term_name<- factor(KeyNegSig2$term_name, levels= KeyNegSig2$term_name)
KeyNegSig2$Termtype<-c("Ribosome","Ribosome","Ribosome","Ribosome",
                       "RNA processing", "RNA processing", "RNA processing")

KeyPosSig<-Protein.posSig[order(Protein.posSig$adjusted_p_value),]
#KeyPosSig2<-KeyPosSig[c(1,6,9,12,13,15,16,21,23,30, 14,18,20,26,28,53, 45,49),]#protein complexes, DNA damage response, heat shock
#KeyPosSig2<-KeyPosSig[c(1,20,26, 45,49),]#protein complexes, DNA damage response, heat shock. (1,6,9,10,14,18,20,26,28,53, 45,49)
KeyPosSig2<-KeyPosSig[c(22,29,50, 45,49),]#membranes, heat shock.
KeyPosSig2$Termtype<-c("Membrane","Membrane", "Membrane",
                       "Heat shock response", "Heat shock response")
KeyPosSig2$term_name<- factor(KeyPosSig2$term_name, levels= KeyPosSig2$term_name)


plot.KeyTerms<- ggplot(KeyPosSig2, 
                          aes(x=term_name, y=-log2(adjusted_p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Biological terms")+
  ggtitle ("Positively correlated with aneuploidy")+
  theme_classic()+
  coord_flip()+
  ylim(0,13)
plot.KeyTerms
# 6x4
# plot.PosCorr.AneuScore.terms



##
##origin of terms pos and negatively correlated with aneuploidy score
# this did not end up in the paper, but interesting
Protein.AnScore.Source<- rbind.fill.matrix(t(table(Protein.posSig$source)), t(table(Protein.negSig$source)))
row.names(Protein.AnScore.Source)<- c("Positive", "Negative") 
Protein.AnScore.Source2<-t(Protein.AnScore.Source)
GOMF<- c(0,0)
MIRNA<-c(0,0)
TF<-c(0,0)
Protein.AnScore.Source3 <- rbind( Protein.AnScore.Source2[1:3,], GOMF, Protein.AnScore.Source2[ 4:6,], MIRNA, Protein.AnScore.Source2[7,], TF ,Protein.AnScore.Source2[8,] )
row.names(Protein.AnScore.Source3)<- c("CORUM", "GO:BP", "GO:CC", "GO:MF", "HP", "HPA", "KEGG", "MIRNA", "REAC", "TF", "WP") 
Protein.AnScore.Source4<-Protein.AnScore.Source3[c(1,9),]
Protein.AnScore.Source_melt <- melt(Protein.AnScore.Source2, id = colnames)

RNA.AnScore.Source<- rbind(table(RNA.posSig$source), table(RNA.negSig$source))
row.names(RNA.AnScore.Source)<- c("Positive", "Negative") #Weird error: 
RNA.AnScore.Source2<-t(RNA.AnScore.Source)
RNA.AnScore.Source3<-RNA.AnScore.Source2[c(1,9),]
RNA.AnScore.Source_melt <- melt(RNA.AnScore.Source3, id = colnames)


plot.CVterms<- ggplot(Protein.AnScore.Source_melt, 
                      aes(x=Var1, y= value, fill=Var2))+
  geom_col(position = "dodge")+
  ylab("Count") +
  xlab("Enriched terms types")+
  scale_fill_discrete(name = "Correlated with \naneuploidy score")+
  ggtitle ("Origin of protein terms correlated with aneuploidy")+
  theme_classic()

plot.CVterms
# Lots of CORUM terms (protein complexes) pos correlated with aneuploidy
# lots of REAC terms (reactosome genes, signalling pathways) negatively correlated with aneuploidy. 



###### PROTEIN AntiScale/buffer/scale key terms ####
# Figure C/D: Buffering/scaling terms
## Three categories based on difference-- PROTEIN
# categorized by scaling, anti-scaling and buffering
# data generated in Protein_RNA_Corr_min10.Filtered.R

### Run through g:profiler
## !! Run g:profiler below, or upload results from supplementary data **

### Find proteins buffered upon gain: 
Prot.Gain.Buffering<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.Protein.Gain=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Prot.Gain.Buffering.result<-Prot.Gain.Buffering$result

## Proteins buffered upon loss:
Prot.Loss.Buffering<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.Protein.Loss=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Prot.Loss.Buffering.result<-Prot.Loss.Buffering$result

### Find proteins scaling upon gain: 
Prot.Gain.Scale<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.Protein.Gain=="Scaling")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Prot.Gain.Scale.result<-Prot.Gain.Scale$result

## Proteins scaling upon loss:
Prot.Loss.Scale<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.Protein.Loss=="Scaling")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Prot.Loss.Scale.result<-Prot.Loss.Scale$result


### Find proteins anti-scaling upon gain: 
Prot.Gain.AntiScale<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.Protein.Gain=="Anti-Scaling")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Prot.Gain.AntiScale.result<-Prot.Gain.AntiScale$result

## Proteins anti-scaling upon loss:
Prot.Loss.AntiScale<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.Protein.Loss=="Anti-Scaling")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Prot.Loss.AntiScale.result<-Prot.Loss.AntiScale$result

##Buffering terms
P.3Cat.Gain.Buff<-Prot.Gain.Buffering.result[order(Prot.Gain.Buffering.result$p_value),]
P.3Cat.Gain.Buff1<- subset(P.3Cat.Gain.Buff, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.3Cat.Gain.Buff2<-P.3Cat.Gain.Buff1[c(1, 9,13,60, 28,35,72),]#168, CORUM root, Ribosomal 32, cell cycle, splisosome
P.3Cat.Gain.Buff2$term_name<- factor(P.3Cat.Gain.Buff2$term_name, levels= P.3Cat.Gain.Buff2$term_name)
P.3Cat.Gain.Buff2$Termtype<-c("CORUM root", "Ribosomal", "Ribosomal", "Ribosomal", "Cell cycle", "Cell cycle", "Cell cycle")

P.3Cat.Loss.Buff<-Prot.Loss.Buffering.result[order(Prot.Loss.Buffering.result$p_value),]
P.3Cat.Loss.Buff1<- subset(P.3Cat.Loss.Buff, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.3Cat.Loss.Buff2<-P.3Cat.Loss.Buff1[c(1, 7,9,75, 91,105,113),]#180, CORUM root,  Ribosomal 24, cell cycle, splicing,
P.3Cat.Loss.Buff2$term_name<- factor(P.3Cat.Loss.Buff2$term_name, levels= P.3Cat.Loss.Buff2$term_name)
P.3Cat.Loss.Buff2$Termtype<-c("CORUM root", "Ribosomal", "Ribosomal", "Ribosomal", "Cell cycle", "Cell cycle", "Cell cycle")

## Scaling terms
P.3Cat.Gain.Scale<-Prot.Gain.Scale.result[order(Prot.Gain.Scale.result$p_value),]
#P.3Cat.Gain.Scale1<- subset(P.3Cat.Gain.Scale, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.3Cat.Gain.Scale2<-P.3Cat.Gain.Scale[c(2,3,4,5,14),]# 
P.3Cat.Gain.Scale2$term_name<- factor(P.3Cat.Gain.Scale2$term_name, levels= P.3Cat.Gain.Scale2$term_name)
P.3Cat.Gain.Scale2$Termtype<-c("Membrane", "Membrane", "Membrane", "Membrane", "Authophagy")

P.3Cat.Loss.Scale<-Prot.Loss.Scale.result[order(Prot.Loss.Scale.result$p_value),]
P.3Cat.Loss.Scale1<- subset(P.3Cat.Loss.Scale, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.3Cat.Loss.Scale2<-P.3Cat.Loss.Scale1[c(1,3,6,10,13),]#36, metabolism, 
P.3Cat.Loss.Scale2$term_name<- factor(P.3Cat.Loss.Scale2$term_name, levels= P.3Cat.Loss.Scale2$term_name)
P.3Cat.Loss.Scale2$Termtype<-c("Metabolism","Metabolism","Metabolism","Metabolism","Metabolism")

## anti-scaling terms
P.3Cat.Gain.AntiScale<-Prot.Gain.AntiScale.result[order(Prot.Gain.AntiScale.result$p_value),]
P.3Cat.Gain.AntiScale1<- subset(P.3Cat.Gain.AntiScale, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.3Cat.Gain.AntiScale2<-P.3Cat.Gain.AntiScale1[c(2,4,26,35, 3,5,16),]#58 total,  extracellular(adhesion, matrix, cell movement), ribosome, 
P.3Cat.Gain.AntiScale2$term_name<- factor(P.3Cat.Gain.AntiScale2$term_name, levels= P.3Cat.Gain.AntiScale2$term_name)
P.3Cat.Gain.AntiScale2$Termtype<-c("Extracellular","Extracellular","Extracellular","Extracellular", "Ribosomal","Ribosomal","Ribosomal")

P.3Cat.Loss.AntiScale<-Prot.Loss.AntiScale.result[order(Prot.Loss.AntiScale.result$p_value),]
P.3Cat.Loss.AntiScale1<- subset(P.3Cat.Loss.AntiScale, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.3Cat.Loss.AntiScale2<-P.3Cat.Loss.AntiScale1[c(2,4,11,26),]#135 total, extracellular (adhesion, matrix, cell movement), 
P.3Cat.Loss.AntiScale2$term_name<- factor(P.3Cat.Loss.AntiScale2$term_name, levels= P.3Cat.Loss.AntiScale2$term_name)
P.3Cat.Loss.AntiScale2$Termtype<-c("Extracellular","Extracellular","Extracellular","Extracellular")


ggplot(P.3Cat.Gain.Scale2, 
                          aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key RNA terms")+
  ggtitle ("Biological terms enriched: Protein gain Scale")+
  theme_classic()+
  ylim(0,160)+
  coord_flip()
#7x5
# plot.gprofiler.3Term.Protein.Gain.Scale_v3

# Write those terms as csv files
#setwd()
#write.csv(Prot.Loss.AntiScale.result[,1:13], 
#          file =paste("gprofiler.Protein.Loss.3cat.AntiScale.csv",
#                      sep=','), row.names = TRUE)

#P.3Cat.Gain.Scale<- read.csv("gprofiler.Protein.Gain.3cat.Scale.csv")
#P.3Cat.Loss.Scale<- read.csv("gprofiler.Protein.Loss.3cat.Scale.csv")



###### RNA     AntiScale/buffer/scale key terms ####
## Figure C/D: Buffering/scaling terms
## Three categories based on difference-- RNA
# categorized by scaling, anti-scaling and buffering
# data generated in Protein_RNA_Corr_min10.Filtered.R

### Run through g:profiler
## !! Run g:profiler below, or upload results from supplementary data **

### Find RNAs buffered upon gain: 
RNA.Gain.Buffering<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.RNA.Gain=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
RNA.Gain.Buffering.result<-RNA.Gain.Buffering$result

## RNAs buffered upon loss:
RNA.Loss.Buffering<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.RNA.Loss=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
RNA.Loss.Buffering.result<-RNA.Loss.Buffering$result

### Find RNAs scaling upon gain: 
RNA.Gain.Scale<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.RNA.Gain=="Scaling")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
RNA.Gain.Scale.result<-RNA.Gain.Scale$result

## RNAs scaling upon loss:
RNA.Loss.Scale<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.RNA.Loss=="Scaling")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
RNA.Loss.Scale.result<-RNA.Loss.Scale$result


### Find RNAs anti-scaling upon gain: 
RNA.Gain.AntiScale<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.RNA.Gain=="Anti-Scaling")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
RNA.Gain.AntiScale.result<-RNA.Gain.AntiScale$result

## RNAs anti-scaling upon loss:
RNA.Loss.AntiScale<- gost(
  subset(CN.Diff.xRNA.yProt.ThreeGroups, Three.RNA.Loss=="Anti-Scaling")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.xRNA.yProt.ThreeGroups$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
RNA.Loss.AntiScale.result<-RNA.Loss.AntiScale$result

##Buffering terms RNA
R.3cat.Gain.Buff<-RNA.Gain.Buffering.result[order(RNA.Gain.Buffering.result$p_value),]
R.3cat.Gain.Buff1<- subset(R.3cat.Gain.Buff, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
R.3cat.Gain.Buff2<-R.3cat.Gain.Buff1[c(1,2),]#168, CORUM root, Ribosomal 32, cell cycle, splisosome
R.3cat.Gain.Buff2$term_name<- factor(R.3cat.Gain.Buff2$term_name, levels= R.3cat.Gain.Buff2$term_name)
R.3cat.Gain.Buff2$Termtype<-c("mRNA processing", "mRNA processing")

R.3cat.Loss.Buff<-RNA.Loss.Buffering.result[order(RNA.Loss.Buffering.result$p_value),]
R.3cat.Loss.Buff1<- subset(R.3cat.Loss.Buff, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
R.3cat.Loss.Buff2<-R.3cat.Loss.Buff1[c(1),]#180, CORUM root,  Ribosomal 24, cell cycle, splicing,
R.3cat.Loss.Buff2$term_name<- factor(R.3cat.Loss.Buff2$term_name, levels= R.3cat.Loss.Buff2$term_name)
R.3cat.Loss.Buff2$Termtype<-c("Extracellular")

## Scaling terms RNA
R.3cat.Gain.Scale<-RNA.Gain.Scale.result[order(RNA.Gain.Scale.result$p_value),]
R.3cat.Gain.Scale1<- subset(R.3cat.Gain.Scale, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
R.3cat.Gain.Scale2<-R.3cat.Gain.Scale1[c(1, 5),]#4 total. signaling pathways?
R.3cat.Gain.Scale2$term_name<- factor(R.3cat.Gain.Scale2$term_name, levels= R.3cat.Gain.Scale2$term_name)
R.3cat.Gain.Scale2$Termtype<-c("CORUM root", "Senescence")

R.3cat.Loss.Scale<-RNA.Loss.Scale.result[order(RNA.Loss.Scale.result$p_value),]
R.3cat.Loss.Scale1<- subset(R.3cat.Loss.Scale, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
R.3cat.Loss.Scale2<-R.3cat.Loss.Scale1[c(1,7,5, 9, 65),]#36, ribosome, rna process 11,12,58, 16,26,43,
R.3cat.Loss.Scale2$term_name<- factor(R.3cat.Loss.Scale2$term_name, levels= R.3cat.Loss.Scale2$term_name)
R.3cat.Loss.Scale2$Termtype<-c("metabolic process","metabolic process","metabolic process",
                               "CORUM root", "cell cycle")

## anti-scaling terms RNA
R.3cat.Gain.AntiScale<-RNA.Gain.AntiScale.result[order(RNA.Gain.AntiScale.result$p_value),]
R.3cat.Gain.AntiScale1<- subset(R.3cat.Gain.AntiScale, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
R.3cat.Gain.AntiScale2<-R.3cat.Gain.AntiScale1[c(1,3,5,6,17),]#58 total,  extracellular(adhesion, matrix, cell movement), ribosome, 
R.3cat.Gain.AntiScale2$term_name<- factor(R.3cat.Gain.AntiScale2$term_name, levels= R.3cat.Gain.AntiScale2$term_name)
R.3cat.Gain.AntiScale2$Termtype<-c("Extracellular","Extracellular","Extracellular","Extracellular","Extracellular")

R.3cat.Loss.AntiScale<-RNA.Loss.AntiScale.result[order(RNA.Loss.AntiScale.result$p_value),]
R.3cat.Loss.AntiScale1<- subset(R.3cat.Loss.AntiScale, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
R.3cat.Loss.AntiScale2<-R.3cat.Loss.AntiScale1[c(1,4,7,14,21),]#135 total, extracellular (adhesion, matrix, cell movement), 
R.3cat.Loss.AntiScale2$term_name<- factor(R.3cat.Loss.AntiScale2$term_name, levels= R.3cat.Loss.AntiScale2$term_name)
R.3cat.Loss.AntiScale2$Termtype<-c("Extracellular","Extracellular","Extracellular","Extracellular","Extracellular")


ggplot(R.3cat.Loss.Scale2, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key RNA terms")+
  ggtitle ("Biological terms enriched: RNA loss scale")+
  theme_classic()+
  ylim(0,160)+
  coord_flip()
#7x5
# plot.gprofiler.3Term.RNA.Loss.Scale

# Write those terms as csv files
#setwd()
#write.csv(RNA.Gain.Scale.result[,1:13], 
#          file =paste("gprofiler.RNA.Gain.3cat.Scale.csv",
#                      sep=','), row.names = TRUE)



###### High aneuploid cells vs low aneuploid cells, PROTEIN buffering key terms ####
# Sup Fig 5
#  PROTEIN Buffering 
# data generated in Protein_RNA_expression.PerCell.R and Protein_expression_data_GeneScore.R

### Run through g:profiler
## !! Run g:profiler below, or upload results from supplementary data **


### Find proteins buffered upon gain: 
HighAneu.Prot.Gain.Buffering.result<- gost(
  subset(CN.Diff.RNA.Prot_HighPloidy, Three.Protein.Gain=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.RNA.Prot_HighPloidy$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
HighAneu.Prot.Gain.Buffering.result<-HighAneu.Prot.Gain.Buffering.result$result


HighAneu.Prot.Loss.Buffering.result<- gost(
  subset(CN.Diff.RNA.Prot_HighPloidy, Three.Protein.Loss=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.RNA.Prot_HighPloidy$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
HighAneu.Prot.Loss.Buffering.result<-HighAneu.Prot.Loss.Buffering.result$result


## Proteins buffered upon gain, low aneuploid cells:
LowAneu.Prot.Gain.Buffering.result<- gost(
  subset(CN.Diff.RNA.Prot_LowPloidy, Three.Protein.Gain=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.RNA.Prot_LowPloidy$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
LowAneu.Prot.Gain.Buffering.result<-LowAneu.Prot.Gain.Buffering.result$result


## Proteins buffered upon loss:
LowAneu.Prot.Loss.Buffering.result<- gost(
  subset(CN.Diff.RNA.Prot_LowPloidy, Three.Protein.Loss=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.RNA.Prot_LowPloidy$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
LowAneu.Prot.Loss.Buffering.result<-LowAneu.Prot.Loss.Buffering.result$result

##Buffering terms HIGH ANEUPLOID CELL GROUP: 
P.Gain.Buff.HighAn<-HighAneu.Prot.Gain.Buffering.result[order(HighAneu.Prot.Gain.Buffering.result$p_value),]
P.Gain.Buff.HighAn1<- subset(P.Gain.Buff.HighAn, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.Gain.Buff.HighAn2<-P.Gain.Buff.HighAn1[c(1, 2,4,6, 14,17,54),]#168, CORUM root, Ribosomal 32, cell cycle, splisosome
P.Gain.Buff.HighAn2$term_name<- factor(P.Gain.Buff.HighAn2$term_name, levels= P.Gain.Buff.HighAn2$term_name)
P.Gain.Buff.HighAn2$Termtype<-c("CORUM root", "RNA processing","RNA processing","RNA processing","Ribosomal", "Ribosomal", "Ribosomal")

P.Loss.Buff.HighAn<-HighAneu.Prot.Loss.Buffering.result[order(HighAneu.Prot.Loss.Buffering.result$p_value),]
P.Loss.Buff.HighAn1<- subset(P.Loss.Buff.HighAn, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.Loss.Buff.HighAn2<-P.Loss.Buff.HighAn1[c(1, 2,4,10, 7,12,17),]#168, CORUM root, Ribosomal 32, cell cycle, splisosome
P.Loss.Buff.HighAn2$term_name<- factor(P.Loss.Buff.HighAn2$term_name, levels= P.Loss.Buff.HighAn2$term_name)
P.Loss.Buff.HighAn2$Termtype<-c("CORUM root", "RNA processing","RNA processing","RNA processing","Ribosomal", "Ribosomal", "Ribosomal" )

##Buffering terms LOW ANEUPLOID CELL GROUP: 
P.Gain.Buff.LowAn<-LowAneu.Prot.Gain.Buffering.result[order(LowAneu.Prot.Gain.Buffering.result$p_value),]
#P.Gain.Buff.LowAn1<- subset(P.Gain.Buff.LowAn, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.Gain.Buff.LowAn2<-P.Gain.Buff.LowAn[c(1,3, 2,6,23, 9),]#168, CORUM root, Ribosomal 32, cell cycle, splisosome
P.Gain.Buff.LowAn2$term_name<- factor(P.Gain.Buff.LowAn2$term_name, levels= P.Gain.Buff.LowAn2$term_name)
P.Gain.Buff.LowAn2$Termtype<-c("CORUM root", "CORUM root", "RNA processing","RNA processing","RNA processing", "Ribosomal")


P.Loss.Buff.LowAn<-LowAneu.Prot.Loss.Buffering.result[order(LowAneu.Prot.Loss.Buffering.result$p_value),]
#P.Loss.Buff.LowAn1<- subset(P.Loss.Buff.LowAn, source %in% c("CORUM" ,"GO:BP", "GO:MF", "WP", "KEGG", "REAC"))
P.Loss.Buff.LowAn2<-P.Loss.Buff.LowAn[c(1, 2, 3,4),]#168, CORUM root, Ribosomal 32, cell cycle, splisosome
P.Loss.Buff.LowAn2$term_name<- factor(P.Loss.Buff.LowAn2$term_name, levels= P.Loss.Buff.LowAn2$term_name)
P.Loss.Buff.LowAn2$Termtype<-c("RNA processing", "Ribosomal", "CORUM root", "CORUM root")

ggplot(P.Gain.Buff.HighAn2, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key RNA terms")+
  ggtitle ("Biological terms enriched: Protein gain Buffered")+
  theme_classic()+
  ylim(0,170)+
  coord_flip()
#7x5
# plot.gprofiler.HighAneuploid.Gain.Protein.Buffered

# Write those terms as csv files
setwd()
write.csv(P.Loss.Buff.LowAn[,1:13], 
          file =paste("gprofiler.LowAneu.Protein.Loss.Buffer.csv",
                      sep=','), row.names = TRUE)


###### TCGA Ovarian terms enriched in buffered genes ####

### TCGA Ovarian terms enriched: 
### Run through g:profiler
Ovarian.Gain.result<- gost(
  subset(t.test_TCGA_ovarian.ThreeGroups, Three.Protein.Gain=="Buffering" )$Gene,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(t.test_TCGA_ovarian.ThreeGroups$Gene),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Ovarian.Gain.result<-Ovarian.Gain.result$result
# CORUM root

Ovarian.Loss.result<- gost(
  subset(t.test_TCGA_ovarian.ThreeGroups, Three.Protein.Loss=="Buffering" )$Gene,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(t.test_TCGA_ovarian.ThreeGroups$Gene),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Ovarian.Loss.result<-Ovarian.Loss.result$result
# CORUM root and non-membrane-bounded organelle (ribosomal)

R.Ovarian.Gain.result<- gost(
  subset(t.test_TCGA_ovarian.ThreeGroups, Three.RNA.Gain=="Buffering" )$Gene,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(t.test_TCGA_ovarian.ThreeGroups$Gene),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
R.Ovarian.Gain.result<- R.Ovarian.Gain.result$result
# None

R.Ovarian.Loss.result<- gost(
  subset(t.test_TCGA_ovarian.ThreeGroups, Three.RNA.Loss=="Buffering" )$Gene,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(t.test_TCGA_ovarian.ThreeGroups$Gene),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
R.Ovarian.Loss.result<-R.Ovarian.Loss.result$result
# three terms: all extracellular terms

write.csv(Ovarian.Loss.result, file="gprofiler.Protein.Loss.TCGA.Ovarian.Buffer.csv")
write.csv(Ovarian.Gain.result, file="gprofiler.Protein.Gain.TCGA.Ovarian.Buffer.csv")
write.csv(R.Ovarian.Loss.result, file="gprofiler.RNA.Loss.TCGA.Ovarian.Buffer.csv")



###### Low/med/high expression variance (in gain/loss condition)  (Did not use this) ####
### Find terms enriched in buffered proteins based on low-medium-high variance
## Did not end up using this in the paper. 

## Proteins low variance upon loss:
CCLA.LowVar.Loss.result<- gost(
  subset(t.test_prot_category, -log(abs(Var.Loss)) > 0)$Protein_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(t.test_prot_category$Protein_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
CCLA.LowVar.Loss.result<-CCLA.LowVar.Loss.result$result
CCLA.LowVar.Loss.result<-CCLA.LowVar.Loss.result[order(CCLA.LowVar.Loss.result$p_value),]
CCLA.LowVar.Loss.result$term_name
# low variance: Membrane, extracellular matrix, motility/adhesion, signalling
# high variance:  metabolism, RNA processing, HPA terms 
# medium variance (-1 to 0): CORUM ROOT, cytosol, DNA repair, cell cycle, TF, HPA terms

## Proteins low variance upon gain:
CCLA.LowVar.Gain.result<- gost(
  subset(t.test_prot_category, -log(abs(Var.Loss)) > -1 & -log(abs(Var.Loss)) < 0)$Protein_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(t.test_prot_category$Protein_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
CCLA.LowVar.Gain.result<-CCLA.LowVar.Gain.result$result
CCLA.LowVar.Gain.result<-CCLA.LowVar.Gain.result[order(CCLA.LowVar.Gain.result$p_value),]
CCLA.LowVar.Gain.result$term_name
# low var: membrane associated, extracellular matrix, migration/adhesion, ribosome
# high var: RNA processing, metabolism/catabolism, HPA terms,  
# Medium variance (-1 to 0): CORUM ROOT, cytosol, DNA repair, cell cycle, TF, HPA terms

### Conclusion: 
## Buffered genes tend to have lower variation. 
## Genes with low variability in cells aneuploid for that chromosome are: membrane and extracellular and Ribosome
## Genes with high variability in cells aneuploid for that chromosome are: RNA processing, metabolism, and HPA terms (tissue specific)



###### Double buffered genes: both CCLE and (yeast, down syndrone, Stingele, and Ovarian tumor) buffered #####
## only yeach & CCLE double buffered terms had significant term enrichment
## for other datasets, there just were't enough genes to reach significance. 

## got datasetw from: 
# Protein_DownSyndrome.R
# Protein_filter_Stingele.R
# TCGA_ProteinDiff_Ovarian.R

## Sub-Step 1: get data from different data analysis files, or download: 
## got datasetw from:  Yeast_prteomics.R 
setwd() #! set this to correct location! 
Yeast.Prot<- read_xlsx("yeast-human aneuploidy.xlsx") 
# CCLE y=`Human Tri vs. Di - protein`, Yeastx=`Yeast TMT - protein`)

#get data: Protein_filter_Stingele.R
## Stingle RPE1 and HCT116 Ts5 vs control 
setwd()
Stingele.DepmapProt.Chrm5.3<- read.csv("Stingele.DepmapProt.Chrm5.csv")

## got datasetw from: 
# Protein_DownSyndrome.R
## (meandiff = stingele chrm 5genes, cells ts5)
setwd()
Ts21.Prot<- read_xlsx("trisomy 21 comparison.xlsx")
Ts21.Prot ## y=`Ts vs. Ds` (CCLE), x=`Ts21 vs. WT...4`) (down syndrome, Ts21)

## got datasetw from: 
# TCGA_ProteinDiff_Ovarian.R
setwd()
t.test_TCGA_ovarian.ThreeGroups<- read.csv("Protein_TCGA_Difference_Pvalue_min10points_3cat.csv")
Diff.pvalue.TCGAovarianx.CCLEy<-merge(x=t.test_TCGA_ovarian.ThreeGroups, 
                                      y=CN.Diff.xRNA.yProt.ThreeGroups, 
                                      by.x="Gene", 
                                      by.y= "RNA_Name") #2991 genes
Diff.pvalue.TCGAovarianx.CCLEy # TCGA x=Diff.Gain, CCLE y=Protein.Diff.Gain


# Sub-Step 2: Run g:profiler for all 4 datasets
#Stingele (stable aneuploid) & CCLE buffered
Ts5.Gain.Buffering<- gost(
  subset(Stingele.DepmapProt.Chrm5.3, Stingele.DepmapProt.Chrm5.3$meanDiff<0.25&Stingele.DepmapProt.Chrm5.3$Protein.Diff.Tri_Di<0.25)$Gene_Name1,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(Stingele.DepmapProt.Chrm5.3$Gene_Name1),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Ts5.Gain.Buffering.result<-Ts5.Gain.Buffering$result
# no enrichment
# 72 gene double buffered from 372 total


# Down Syndrome & CCLE buffered
Ts21.Prot ## y=`Ts vs. Ds` (CCLE), x=`Ts21 vs. WT...4`) 
DS.Gain.Buffering<- gost(
  subset(Ts21.Prot, `Ts vs. Ds`<0.25&`Ts21 vs. WT...4`<0.25)$Gene_name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(Ts21.Prot$Gene_name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
DS.Gain.Buffering.result<-DS.Gain.Buffering$result
# No significant enrichment
# 15 genes double buffered from 42 total

# Yeast & CCLE buffered
Yeast.Gain.Buffering<- gost(
  subset(Yeast.Prot, `Human Tri vs. Di - protein`<0.25&`Yeast TMT - protein`<0.25)$`Human Gene`,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(Yeast.Prot$`Human Gene`),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Yeast.Gain.Buffering.result<-Yeast.Gain.Buffering$result
# 98 double buffered from 738 total
# Ribosomes, Nonsence-mediated decay and maybe protein complexes

# Ovarian gain and CCLE gain buffer
# Diff.pvalue.TCGAovarianx.CCLEy # TCGA x=Diff.Gain, CCLE y=Protein.Diff.Gain
Ovarian.Gain.Buffering<- gost(
  subset(Diff.pvalue.TCGAovarianx.CCLEy, Diff.Gain<0.25&Protein.Diff.Gain<0.25)$Gene,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(Diff.pvalue.TCGAovarianx.CCLEy$Gene),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Ovarian.Gain.Buffering.result<-Ovarian.Gain.Buffering$result
# 906 genes double buffered from 1313 total
# CORUM root, RNA processing (slicing), Ribosome

# Ovarian loss
# Diff.pvalue.TCGAovarianx.CCLEy # TCGA x=Diff.Gain, CCLE y=Protein.Diff.Gain
Ovarian.Loss.Buffering<- gost(
  subset(Diff.pvalue.TCGAovarianx.CCLEy, Diff.Loss> -0.25&Protein.Diff.Loss> -0.25)$Gene,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(Diff.pvalue.TCGAovarianx.CCLEy$Gene),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Ovarian.Loss.Buffering.result<-Ovarian.Loss.Buffering$result
# CORUM root, Ribosome
#937 double buffered from 1313 total


## Plot terms: 
## Yeast
P.Yeast.CCLE.Buffer<-Yeast.Gain.Buffering.result[order(Yeast.Gain.Buffering.result$p_value),]
P.Yeast.CCLE.Buffer2<-P.Yeast.CCLE.Buffer[c(2,6,17, 5, 28,31,33),]#135 total, extracellular (adhesion, matrix, cell movement), 
P.Yeast.CCLE.Buffer2$term_name<- factor(P.Yeast.CCLE.Buffer2$term_name, levels= P.Yeast.CCLE.Buffer2$term_name)
P.Yeast.CCLE.Buffer2$Termtype<-c("Ribosome","Ribosome","Ribosome", "Protein Complex", 
                                 "RNA processing","RNA processing","RNA processing")


ggplot(P.Yeast.CCLE.Buffer2, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key RNA terms")+
  ggtitle ("Biological terms enriched: Yeast & cancer cell buffered")+
  theme_classic()+
  ylim(0,30)+
  coord_flip()
# 7x5
# plot.gprofiler.Yeast.CCLE.Buffered

write.csv(P.Yeast.CCLE.Buffer[,1:13], file="gprofiler.CCLE.Yeast.Gain.Buffered.csv")

## Ovarian gain: 
P.Ovarian.Gain.CCLE.Buffer<-Ovarian.Gain.Buffering.result[order(Ovarian.Gain.Buffering.result$p_value),]
P.Ovarian.Gain.CCLE.Buffer2<-P.Ovarian.Gain.CCLE.Buffer[c(1,2, 3,5,9, 4,7,12),]#135 total, extracellular (adhesion, matrix, cell movement), 
P.Ovarian.Gain.CCLE.Buffer2$term_name<- factor(P.Ovarian.Gain.CCLE.Buffer2$term_name, levels= P.Ovarian.Gain.CCLE.Buffer2$term_name)
P.Ovarian.Gain.CCLE.Buffer2$Termtype<-c("Protein Complex","Protein Complex", "RNA processing","RNA processing","RNA processing",
                                        "Ribosome","Ribosome","Ribosome")


ggplot(P.Ovarian.Gain.CCLE.Buffer2, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key RNA terms")+
  ggtitle ("Biological terms enriched upon chrm gain: Ovarian  & cancer cell buffered")+
  theme_classic()+
  ylim(0,30)+
  coord_flip()
# 7x5
# plot.gprofiler.Gain.Ovarian.CCLE.Buffered


## Ovarian loss: 
P.Ovarian.Loss.CCLE.Buffer<-Ovarian.Loss.Buffering.result[order(Ovarian.Loss.Buffering.result$p_value),]
P.Ovarian.Loss.CCLE.Buffer2<-P.Ovarian.Loss.CCLE.Buffer[c(1,4, 5,7,14, 6,9),]#135 total, extracellular (adhesion, matrix, cell movement), 
P.Ovarian.Loss.CCLE.Buffer2$term_name<- factor(P.Ovarian.Loss.CCLE.Buffer2$term_name, levels= P.Ovarian.Loss.CCLE.Buffer2$term_name)
P.Ovarian.Loss.CCLE.Buffer2$Termtype<-c("Protein Complex","Protein Complex",
                                        "Ribosome","Ribosome","Ribosome", 
                                        "RNA processing","RNA processing")

ggplot(P.Ovarian.Loss.CCLE.Buffer2, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key RNA terms")+
  ggtitle ("Biological terms enriched upon chrm loss: Ovarian  & cancer cell buffered")+
  theme_classic()+
  ylim(0,30)+
  coord_flip()
# 7x5
# plot.gprofiler.Loss.Ovarian.CCLE.Buffered

setwd() #set to location you want to write data results to
write.csv(P.Ovarian.Loss.CCLE.Buffer[,1:13], file="gProfiler.CCLE.and.Ovarian.Loss.Buffering.csv")
write.csv(P.Ovarian.Gain.CCLE.Buffer[,1:13], file="gProfiler.CCLE.and.Ovarian.Gain.Buffering.csv")
write.csv(P.Yeast.CCLE.Buffer[,1:13], file="gProfiler.CCLE.and.Yeast.Buffering.csv")
write.csv(DS.Gain.Buffering.result[,1:13], file="gProfiler.CCLE.and.DownSyndrome.Buffering.csv")
write.csv(Ts5.Gain.Buffering.result[,1:13], file="gProfiler.CCLE.and.Ts5.Buffering.csv")


###### g:profiler for lung cancer cells buffering #####
## data from Protein_RNA_expression.PerCell_v2.R, difference in expression
## looking at only non-small cell lung cancer
# CN.Diff.RNA.Prot_Reliable
CN.Diff.RNA.Prot_Lung<- read.csv("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_Lung.min10cells.csv")

#NSCL non-small cell lung cancer cells
NSCL.Gain.Buffering<- gost(
  subset(CN.Diff.RNA.Prot_Lung, Three.Protein.Gain=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.RNA.Prot_Lung$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
NSCL.Gain.Buffering.Result<-NSCL.Gain.Buffering$result
# CORUM, RNA processing, Ribosome

## Lung Gain: 
P.Lung.Gain.Buffer<-NSCL.Gain.Buffering.Result[order(NSCL.Gain.Buffering.Result$p_value),]
P.Lung.Gain.Buffer2<-P.Lung.Gain.Buffer[c(1,2,3,4,5),]#135 total, extracellular (adhesion, matrix, cell movement), 
P.Lung.Gain.Buffer2$term_name<- factor(P.Lung.Gain.Buffer2$term_name, levels= P.Lung.Gain.Buffer2$term_name)
P.Lung.Gain.Buffer2$Termtype<-c("Protein Complex","Protein Complex",
                                "Ribosome",
                                "RNA processing","RNA processing")
pdf(height=4, width=5, file="plot.gprofiler.NSCL.Gain.Buffered")
ggplot(P.Lung.Gain.Buffer2, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key Protein terms")+
  ggtitle ("Biological terms enriched upon chrm gain: NSCL buffered")+
  theme_classic()+
  ylim(0,75)+
  coord_flip()
dev.off()
# 7x5
# plot.gprofiler.NSCL.Gain.Buffered

write.csv(NSCL.Gain.Buffering.Result[,1:13], file="gProfiler.NSCL.Gain.Buffering.min10.csv")



## NSCL loss

#NSCL non-small cell lung cancer cells
NSCL.Loss.Buffering<- gost(
  subset(CN.Diff.RNA.Prot_Lung, Three.Protein.Loss=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.RNA.Prot_Lung$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
NSCL.Loss.Buffering.Result<-NSCL.Loss.Buffering$result
# CORUM

## Lung Gain: 
P.Lung.Loss.Buffer<-NSCL.Loss.Buffering.Result[order(NSCL.Loss.Buffering.Result$p_value),]
P.Lung.Loss.Buffer2<-P.Lung.Loss.Buffer[c(1,2),]#135 total, extracellular (adhesion, matrix, cell movement), 
P.Lung.Loss.Buffer2$term_name<- factor(P.Lung.Loss.Buffer2$term_name, levels= P.Lung.Loss.Buffer2$term_name)
P.Lung.Loss.Buffer2$Termtype<-c("Protein Complex","Protein Complex")

pdf(height=4, width=5, file="plot.gprofiler.NSCL.Loss.Buffered")
ggplot(P.Lung.Loss.Buffer2, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key Protein terms")+
  ggtitle ("Biological terms enriched upon chrm loss: NSCL buffered")+
  theme_classic()+
  ylim(0,75)+
  coord_flip()
# 7x5
# plot.gprofiler.NSCL.Loss.Buffered
dev.off()

write.csv(NSCL.Loss.Buffering.Result[,1:13], file="gProfiler.NSCL.Loss.Buffering.min10.csv")



###### g:profiler for no low reproducible gene buffering #####
## data from Protein_RNA_expression.PerCell_v2.R, difference in expression
## looking at only genes that have a higher reliability score
# CN.Diff.RNA.Prot_Reliable
CN.Diff.RNA.Prot_Reliable<- read.csv("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_Reliable_Min10Cells.csv")

#NSCL non-small cell lung cancer cells
Reliable.Gain.Buffering<- gost(
  subset(CN.Diff.RNA.Prot_Reliable, Three.Protein.Gain=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.RNA.Prot_Reliable$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Reliable.Gain.Buffering.Result<-Reliable.Gain.Buffering$result
# CORUM, RNA processing, Ribosome

## reliable Gain: 
P.reliable.Gain.Buffer<-Reliable.Gain.Buffering.Result[order(Reliable.Gain.Buffering.Result$p_value),]
P.reliable.Gain.Buffer2<-P.reliable.Gain.Buffer[c(1,11,21, 3,8,9, 10,32,36),]#135 total, extracellular (adhesion, matrix, cell movement), 
P.reliable.Gain.Buffer2$term_name<- factor(P.reliable.Gain.Buffer2$term_name, levels= P.reliable.Gain.Buffer2$term_name)
P.reliable.Gain.Buffer2$Termtype<-c("RNA processing","RNA processing","RNA processing", 
                                    "Protein Complex","Protein Complex", "Protein Complex",
                                    "Ribosome","Ribosome","Ribosome")

pdf(height=4, width=5, file="plot.gprofiler.Reliable.Gain.Buffered")
ggplot(P.reliable.Gain.Buffer2, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key Protein terms")+
  ggtitle ("Biological terms enriched upon chrm gain: NSCL buffered")+
  theme_classic()+
  ylim(0,75)+
  coord_flip()
# 7x5
# plot.gprofiler.Reliable.Gain.Buffered
dev.off()

write.csv(NSCL.Gain.Buffering.Result[,1:13], file="gProfiler.NSCL.Gain.Buffering.min10.csv")



## NSCL loss

#NSCL non-small cell reliable cancer cells
Reliable.Loss.Buffering<- gost(
  subset(CN.Diff.RNA.Prot_Reliable, Three.Protein.Loss=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.RNA.Prot_Reliable$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Reliable.Loss.Buffering.Result<-Reliable.Loss.Buffering$result
# CORUM

## reliable Gain: 
P.reliable.Loss.Buffer<-Reliable.Loss.Buffering.Result[order(Reliable.Loss.Buffering.Result$p_value),]
P.reliable.Loss.Buffer2<-P.reliable.Loss.Buffer[c(1,2,6, 4,11,12, 5,7,10),]#135 total, extracellular (adhesion, matrix, cell movement), 
P.reliable.Loss.Buffer2$term_name<- factor(P.reliable.Loss.Buffer2$term_name, levels= P.reliable.Loss.Buffer2$term_name)
P.reliable.Loss.Buffer2$Termtype<-c("Protein Complex","Protein Complex", "Protein Complex", 
                                    "Ribosome","Ribosome","Ribosome",
                                    "RNA processing","RNA processing","RNA processing")

pdf(height=4, width=5, file="plot.gprofiler.Reliable.Loss.Buffered")
ggplot(P.reliable.Loss.Buffer2, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key Protein terms")+
  ggtitle ("Biological terms enriched upon chrm loss: Reliable buffered")+
  theme_classic()+
  ylim(0,75)+
  coord_flip()
# 7x5
# plot.gprofiler.Reliable.Loss.Buffered
dev.off()

write.csv(Reliable.Loss.Buffering.Result[,1:13], file="gProfiler.Reliable.Loss.Buffering.min10.csv")



###### g:profiler for genes no mutations, no flipping, no low RNA, no low reproducible, minimum 10 cells #####
## data from Protein_RNA_expression.PerCell_v2.R, difference in expression
## looking at only genes with no mutations, no flipped genes, no low RNA genes, and no low reliability genes
CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely<- read.csv("RNA.Protein_Loss.Neutral.Gain_Difference_Pvalue_noMut_noFlip_noLowRNA_nolowRep.min10cells.csv")

#Merge gain
Gain.Buffering.Merge<- gost(
  subset(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely, Three.Protein.Gain=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Gain.Buffering.Merge2<-Gain.Buffering.Merge$result
# 

## Lung Gain: 
P.Merge.Gain.Buffer<-Gain.Buffering.Merge2[order(Gain.Buffering.Merge2$p_value),]
P.Merge.Gain.Buffer2<-P.Merge.Gain.Buffer[c(1,7,11, 3,8, 10,67,78),]#135 total, extracellular (adhesion, matrix, cell movement), 
P.Merge.Gain.Buffer2$term_name<- factor(P.Merge.Gain.Buffer2$term_name, levels= P.Merge.Gain.Buffer2$term_name)
P.Merge.Gain.Buffer2$Termtype<-c("RNA processing","RNA processing", "RNA processing", 
                                 "Protein Complex","Protein Complex",
                                 "Ribosome", "Ribosome", "Ribosome")

ggplot(P.Merge.Gain.Buffer2, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key Protein terms")+
  ggtitle ("Biological terms enriched upon chrm gain: Merge buffered")+
  theme_classic()+
  ylim(0,75)+
  coord_flip()
# 7x5
# plot.gprofiler.Merge.min10.Gain.Buffered.pdf

write.csv(Gain.Buffering.Merge2[,1:13], file="gProfiler.Merge.min10.Gain.Buffering.csv")



## Merge loss

Loss.Buffering.Merge<- gost(
  subset(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely, Three.Protein.Loss=="Buffering")$RNA_Name,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(CN.Diff.RNA.Prot_noMut_noFlip_noLowRNA_nolowRely$RNA_Name),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Loss.Buffering.Merge2<-Loss.Buffering.Merge$result
# 

## Merge loss: 
P.Merge.Loss.Buffer<-Loss.Buffering.Merge2[order(Loss.Buffering.Merge2$p_value),]
P.Merge.Loss.Buffer2<-P.Merge.Loss.Buffer[c(1,2,3, 4,12,24, 5,7,10),]#135 total, extracellular (adhesion, matrix, cell movement), 
P.Merge.Loss.Buffer2$term_name<- factor(P.Merge.Loss.Buffer2$term_name, levels= P.Merge.Loss.Buffer2$term_name)
P.Merge.Loss.Buffer2$Termtype<-c("Protein Complex","Protein Complex","Protein Complex",
                                 "Ribosome", "Ribosome", "Ribosome", 
                                 "RNA processing","RNA processing", "RNA processing")

ggplot(P.Merge.Loss.Buffer2, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key Protein terms")+
  ggtitle ("Biological terms enriched upon chrm loss: Merge buffered")+
  theme_classic()+
  ylim(0,75)+
  coord_flip()
# 7x5
# plot.gprofiler.Merge.min10.Loss.Buffered.pdf

write.csv(Loss.Buffering.Merge2[,1:13], file="gProfiler.Merge.min10.Loss.Buffering.csv")



###### g:profiler double buffered RPPA +/- Mass Spec genes (None enriched) #####
## data from Protein_RNA_expression.PerCell_v2.R, difference in expression
## looking at genes that have antibody-based abundance data (RPPA) and mass spec data
RPPA_data_MassSpec<- read.csv("RPPA.MassSpec.difference.csv")

# I examined both "double buffered" RPPA and Mass spec, as well as only buffered in RPPA data
# both yielded no significant results. too few datapoints. 
#RPPA gain
Gain.Buffering.RPPA.MS<- gost(
  #subset(RPPA_data_MassSpec, Protein.Diff.Gain.RPPA< 0.25 & Protein.Diff.Gain<0.25)$Gene_Symbol, #double buffered RPPA and mass spec
  subset(RPPA_data_MassSpec, Protein.Diff.Gain.RPPA< 0.25)$Gene_Symbol, #RPPA buffered
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(RPPA_data_MassSpec$Gene_Symbol),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Gain.Buffering.RPPA.MS2<-Gain.Buffering.RPPA.MS$result
# Not enough data

## RPPA loss
Loss.Buffering.RPPA.MS<- gost(
  subset(RPPA_data_MassSpec, Protein.Diff.Loss.RPPA> -0.25)$Gene_Symbol, #RPPA buffered
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(RPPA_data_MassSpec$Gene_Symbol),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
Gain.Buffering.RPPA.MS2<-Loss.Buffering.RPPA.MS$result
# Not enough data to get significant results. 