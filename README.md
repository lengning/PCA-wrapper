# PCA-wrapper
PCA codes and an example data set (random numbers generated from Poisson(1000))

Example commands to run the script:
- Rscript PCA.R PCA_example.csv 4 10 T
- Rscript PCA.R PCA_example2.txt 3 0 T
- Rscript PCA.R PCA_example3.tab 5 50 F

The 3rd term should be the name of the input data set. 
Currently the program takes csv files or tab delimited file.
The input file will be treated as a tab delimited file if the suffix is not '.csv'.
Rows are genes and columns are samples. Row names and column names are required in the input file.

The 4th term defines number of PCs to output (define it as k, default k =5)

The 5th term defines the lower limit of detection threshold (default is 0). Genes with max expression below this threshold will be removed.

The 6th term should be T or F. If T is specified, median-by-ratio normalization will be performed prior to PCA analysis (default is F).

To simplify your input command and use the default values, you may run

Rscript PCA.R PCA_example.csv

or

Rscript PCA.R PCA_example.csv 4


The code will first rescale the data (after filtering out low expressers
and normalization, if specified) to gene specific z scores, then perform PCA.
The output files are:

- prefix_PC_pairs.pdf:
pairwise plots of the transformed data; k PCs will be shown
- prefix_loading.csv
: gene loading for the top k PCs
- prefix_sort_by_absloading.csv
: For each of the top k PCs, genes are sorted by their absolute loadings in each PC
- prefix_perc_sdev.csv
: Percentage of SD explained by each PC

The ‘prefix’ is defined as the filename of the input file (the string before the suffix, like the string before ‘.csv’) 


I’m currently using a quick and dirty way to keep the R open when browsing the results - the code will terminate the R process after 10^30 sec.
