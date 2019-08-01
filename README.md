# coloc-utils
R tool to perform colocalisation analysis on genetic variants and downstream analysis. Utilises the R package coloc (https://cran.r-project.org/web/packages/coloc/index.html). GWAS results for a phenotype of interest are compared with eqtl data (e.g. downloaded from GTEx https://gtexportal.org/) to test for colocalisation of SNPs associated with both the phenotype of interest and the expression level of a nearby gene. More information about the colocalisation method can be found here: https://academic.oup.com/hmg/article/24/12/3305/621728

##############################################################################################################################
Required R packages - tested versions:

ggplot2 >= 3.1.0; purrr >= 0.3.0; tibble >= 1.4.2; dplyr >= 0.7.7; tidyr >= 0.8.2; stringr >= 1.3.1; readr >= 1.3.1; forcats >= 0.3.0; coloc >= 3.1; data.table >= 1.11.8; tidyverse >= 1.2.1;

To run:
        Rscript gwas-file.csv tissue_list.txt gtex_folder_path
       
##############################################################################################################################

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

gtex_folder_path:- GTEx folder path containing eQTL data
eQTL association data must be a folder containing CSV files corresponding to different tissues, each of which must contain at least the following columns:

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
