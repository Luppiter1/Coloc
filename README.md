# Coloc
R tool to perform colocalisation analysis on genetic variants. GWAS results for a phenotype of interest are compared with eqtl data (e.g. downloaded from GTEx https://gtexportal.org/) to test for colocalisation of SNPs associated with both the phenotype of interest and the expression level of a nearby gene. More information about the colocalisation method can be found here: https://academic.oup.com/hmg/article/24/12/3305/621728

Input Files:

1) GWAS results

GWAS results file in CSV format containing the following columns:

SNP:      RS number for the SNP being tested
CHR:      Chromosome number
BP:       Genome coordinate od SNP
BETA:     Regression Coefficient
SE:       Standard Error
namesnp:  SNP ID (Chr_Position_Allele1_Allele2)
ALLELE0:  Allele0
ALLELE1:  Allele1
A1FREQ:   Frequency of Allele1
leadsnp:  RS number for lead SNP

2) eQTL data

eQTL association data must be a CSV file or multiple CSV files corresponding to different tissues, which contain the following columns:

