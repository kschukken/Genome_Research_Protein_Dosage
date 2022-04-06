# Genome_Research_Protein_Dosage
Depository of all R files for Extensive protein dosage compensation in aneuploid human cancers

This is a repository for the 2022 Genome Research article "Extensive protein dosage compensation in aneuploid human cancers" by Klaske M Schukken and Jason M Sheltzer

The READ_ME_R.xlsx is an overview of which R file generates which datasets and graphs, also listed below.

Many of the code files need additional datasets retrieved from outside sources-- The source of the datasets are often listed in the ## Get data ## section of the R files. If data source is not listed in the R code itself, it will be credited in the "Extensive protein dosage compensation in aneuploid human cancers" article. Several intermediate datasets are made available for ease of use.

All code was build on R version 3.6.3 GUI 1.70 El Capitan build, and R studio version 1.1.383.

In order to run the code, open the .R files, download the indicated datasets, set the proper working directory, and run through the code.

Updates relative to the previous bioRxiv iteration include an aditional non-coding RNA analysis, additional secondary statistical controls for our buffering factor analysis, updated cell line information, and an update to figure names and numbers. 

ORGANIZATION

R CODE:

Protein_RNA_filtered_CellLine.R

•	All figures

•	Filter protein expression and RNA expression data to isolate cell lines with RNA expression, protein expression, and arm call data
Protein_filtered_analysis_3cat.R

•	Figure 1, S1

•	Calculate protein expression difference for all genes upon chrm arm gain or loss (ex. gain chrm 5, protein expression diff in chrm 1,2,3,4,5,6, etc. )


RNA_filtered_analysis_3cat.R

•	Figure 1, S1

•	Calculates RNA expression difference for all genes upon chrm arm gain or loss (ex. gain chrm 5, RNA expression diff in chrm 1,2,3,4,5,6, etc. )


ChrmArm.Scatterplot.Protein.RNA.R

•	Figure 1

•	Protein and RNA expression difference per chrm arm scatterplot (ex. Gain chromosome 5q: protein expression difference on chrm 5q, and on other 
chromosomes)

DNA_Protein_RNA_filtered.R

•	Figure 1, S1

•	Make graphs of chrm arm difference from neutral-ploidy, for specific cells (ex. MCF7)

Protein_RNA_Heatmap_filtered_min10percategory.R

•	Figure 1

•	Generate heatmaps: protein/RNA expression difference upon chrm arm gain/loss, and additional summary heatmap

Protein_RNA_expression.PerCell_v3.R

•	Figure 2,3,5,6, S2, S3, S4, S5, S8

•	Difference upon chrm gain/loss dataset(s), categorize (scaling, buffering, anti-scaling). Boxplots of RNA or Protein expression vs DNA copy number 
category (gain, neutral, loss). Near diploid vs near triploid difference analysis. High aneuploidy score and low aneuploidy score difference analysis. 
Non-coding RNA (microRNA) analysis.

Protein_RNA_1d.plot_v2

•	Figure 3, 6, S2, S6

•	Plot 1dimentional density plot of gene expression changes per group (Protein and RNA). Also look at low vs high aneuploidy cells only

gProfile_Quantile_Bargraph_v2.R

•	Figure 3, 4, S5, S6, S7, S8

•	Make graphs of key Gene Ontology enrichment terms (for multiple analysis')

Yeast_proteomics_v2.R

•	Figure S6

•	Aneuploid yeast analysis (Protein & RNA)

DownSyndrome_analysis_v2.R

•	Figure S6

•	Down syndrome fibroblast analysis (Protein and RNA)

Stingele_Protein_filter_v2.R

•	Figure S6

•	Stable aneuploid cell line analysis (Protein and RNA)

Aneuploid_Score.R

•	Figure 4, S7

•	Generate cellular aneuploidy scores

Protein_expression_data_GeneScore.R

•	Figure 4, S7, S8

•	Protein expression correlate with cellular aneuploidy, and high/low aneuploid cell only analysis

CCLE_RNA_Analysis_geneAn.R

•	Figure S7

•	RNA aneuploidy score correlation analysis

Variance_aneuploidy.R

•	Figure 5

•	Calculate RNA and Protein neutral-ploidy variance

Protein_buffering_factors_v2.R

•	Figure 5

•	Protein buffering factors. Calculate ROC AUC and make boxplots of scores/values

TSG.OG_CNV_Difference_v2.R

•	Figure 6

•	Oncogene and tumor suppressor gene difference upon chromosome copy number changes and gene copy number changes

TCGA_ProteinDiff_Ovarian.R

•	Figure 7, S9

•	Analysis of ovarian tumor data from TCGA. 



ADDITIONAL FILES:

arm_calls_data_bendavid.csv

•	Aneuploid arm calls per chromosome arm per cell line analyzed. Used to generate aneuploidy score per cell line, and used to identify which cell lines 
have chromosome arm gains and losses per chromosome arm. Data from Cohen-Sharir et al. 2021 Nature, from the lab of Dr. Uri Ben-David.

bp_per_arm.csv

•	number of basepairs per chromosome arm. Used to generate a basepair-based aneuploidy score

Score2.aneuploid.cell_line.csv

•	Cellular aneuploidy scores per cell line. Multupiple forms of the aneuploidy score are calculated, the gene_ploidy_score was used in the manuscript 
unless otherwise indicated.

RNA_Protein_Shared_cells.csv

•	Cell lines with both RNA expression data and protein expression data. Not all of these cells are also present in the Cohen-Sharir arm call data.

Protein_location_info.csv

•	Information about the chromosome arm location per gene.

Protein_ID_info.csv

•	Protein_IDs matched with gene symbol and uniprot accesion

RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points.csv

•	Dataset of gene expression difference on cells with chromosome gain or loss, relative to cells with a neutral ploidy for chromosome arm a gene is 
located on.

RNA.Protein_loss.neutral.gain_Difference_Pvalue_min10points_3cat.csv

•	Dataset of gene expression difference on cells with chromosome gain or loss, relative to cells with a neutral ploidy for chromosome arm a gene is 
located on. with gene classifications as either "Anti-Scaling", "Buffering" or "Scaling" upon chromosome gain and loss

RNA.Protein_loss.neutral.gain_Difference_Pvalue_*****.csv

•	Dataset of gene expression difference, as above, for a subset of cells and/or genes corresponding to the indicated control dataset. genes 
classifications as either "Anti-Scaling", "Buffering" or "Scaling" upon chromosome gain and loss. 

Trisomy 21 comparison.xlsx

•	Protein expression difference between Down Syndrome (Ts21) fibroblasts and control fibroblasts, set relative to cancer cell line protein expression 
difference. Yeast data from: Liu et al. 2017, Nature Communications. Systematic proteome and proteostasis profiling in human Trisomy 21 fibroblast cells

yeast-human aneuploidy.xlsx

•	Yeast protein expression difference upon chromosome gain. Data originated from: Dephoure et al. 2014, eLife; 3:e03023 Quantitative proteomic analysis 
reveals posttranslational responses to aneuploidy in yeast

AmplificationRatio_perGene1.75.csv

•	Calculation of the percent of cells in the filtered dataset with gene amplifications per gene. Used in initial factor analysis, but not included in 
the final manuscript. Data from CCLE gene copy number data.

Protein.AllFactors.csv

•	A list of all genes, difference upon chromosome gain or loss, and factor scores/values per gene. 5' UTR length and 3' UTR length were not included, 
because there are frequently multiple 5' or 3' UTR lengths per gene, and this duplication can skew data when investigating other factors. See 
Protein_buffering_factors_v2.R, to incorporate 5' UTR length and 3'UTR length if desired.

Protein.AllFactors_Ovarian.csv

•	A list of all genes, difference upon chromosome gain or loss, and factor scores/values per gene for the ovarian tumor dataset. 

Factor.ROCAUC.correlationscore.csv

•	A list of all investigated factors, their AUC ROC values with regard to buffering, scaling and anti-scaling genes, and the correlation coefficient 
with RNA and Protein expression difference.

Go term gene list

•	A list of gene names or protein IDs associated with each gene group examined in the paper. GO terms can also be found on the gene ontology website. 

