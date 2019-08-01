# Coloc
R tool to perform colocalisation analysis on genetic variants. GWAS results for a phenotype of interest are compared with eqtl data (e.g. downloaded from GTEx https://gtexportal.org/) to test for colocalisation of SNPs associated with both the phenotype of interest and the expression level of a nearby gene. More information about the colocalisation method can be found here: https://academic.oup.com/hmg/article/24/12/3305/621728

Required R packages - tested versions:

ggplot2   3.1.0     purrr     0.3.0
tibble    1.4.2     dplyr       0.7.7
tidyr     0.8.2     stringr     1.3.1
readr     1.3.1     forcats     0.3.0
coloc     3.1       data.table  1.11.8
tidyverse 1.2.1

Input Files:
1) GWAS results

GWAS results file in CSV format must contain at least the following columns:

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





2) eQTL data

eQTL association data must be a CSV file or multiple CSV files corresponding to different tissues, which must contain at least the following columns:

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

Questions about this script or suggestions: blakeleyp@gmail.com
