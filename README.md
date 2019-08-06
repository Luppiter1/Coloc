# coloc-utils
R tool to perform colocalisation analysis on genetic variants and perform downstream analysis. Utilises the R package coloc (https://cran.r-project.org/web/packages/coloc/index.html). GWAS results for a phenotype of interest are compared with eqtl data (e.g. downloaded from GTEx https://gtexportal.org/) to test for colocalisation of SNPs associated with both the phenotype of interest and the expression level of a nearby gene. More information about the colocalisation method can be found here: https://academic.oup.com/hmg/article/24/12/3305/621728

##############################################################################################################################

Required R packages - tested versions:

ggplot2 >= 3.1.0; purrr >= 0.3.0; tibble >= 1.4.2; dplyr >= 0.7.7; tidyr >= 0.8.2; stringr >= 1.3.1; readr >= 1.3.1; forcats >= 0.3.0; coloc >= 3.1; data.table >= 1.11.8; tidyverse >= 1.2.1;

To run:
        Rscript gwas-file.csv tissue_list.txt base_folder_path bp_window gwas_size eqtl_size
        
Options

bp_window = The range in base pairs surrounding the lead snp which is used to test for colocalisation (recommended value = 200000)

gwas_size = Number of samples in GWAS dataset 

eqtl_size = Number of samples in eQTL dataset

base_folder_path should contain a folder called gtex_eqtl, containing the gtex eqtl results files (ee below for details)


##############################################################################################################################

OUTPUT:

colocABFtest-results.csv contains results summary of colocalisation tests performed on each region of interest in different tissues.

nsnp:           Number of SNPs within window

PP.H0.abf:      Probability that neither trait has a genetic association in the region

PP.H1.abf:      Probability only trait 1 (GWAS-phenotype) has a genetic association in the region

PP.H2.abf:      Probability only trait 2 (eQTL-expression) has a genetic association in the region

PP.H3.abf:      Probability both traits are associated, but with different causal variants

PP.H4.abf:      Probability both traits are associated and share a single causal variant

tissue:         Tissue for eQTL gene expression association

Lead_SNP        Lead SNP RS number

SNP_LOC         Lead SNP genome coordinate

##############################################################################################################################

INPUT GWAS: 

gwas_file.csv:- GWAS results file in CSV format must contain at least the following columns:

SNP:      RS number for the SNP being tested

chr:      Chromosome number

BP:       Genome coordinate of SNP

BETA:     Regression Coefficient

SE:       Standard Error

namesnp:  SNP ID (Chr_Position_Allele1_Allele2)

ALLELE0:  Reference Allele 

ALLELE1:  Minor Allele

A1FREQ:   Frequency of Allele1

leadsnp:  RS number for lead SNP

##############################################################################################################################

INPUT eQTL:

base_folder_path:- path containing GTEx eqtle data folder
Create a folder in this directory called "gtex_eqtl". eQTL association data must be placed inside this folder, and must be in CSV format. Each file should correspond to a different tissue-chromosome combination (i.e. 23 chromosome files per tissue), each of which must contain at least the following columns:

chr:      Chromosome number

position: Genome coordinate of SNP

slope:    Regression Coefficient

slope_se: Standard Error
  
newid:    SNP ID (Chr_Position_Allele1_Allele2)

ALLELE1:  Minor Allele

ALLELE2:  Reference Allele

maf:      Frequency of Allele1

leadsnp:  RS number for lead SNP

Please note that other columns apart from those listed will be permited, but only those listed will be used for colocalisation analysis.
Files within gtex_folder_path should be named with an underscore seperating the tissue name and chromosome number like so: [tissue]\_[chromosome].txt

##############################################################################################################################

tissue_list.txt:- List of tissues names to be analysed

Names witin this file should be the tissue name prefixes of the filenames in gtex_folder_path. So for the filename "Artery_Coronary_11.txt" the tissue name contained in tissue_list.txt should be "Artery_Coronary"

##############################################################################################################################

Questions about this script or suggestions: blakeleyp@gmail.com
