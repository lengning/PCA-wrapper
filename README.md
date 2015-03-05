# PCA-wrapper
PCA codes and an example data set (random numbers generated from Poisson(1000))

An example command to run the script:
Rscript PCA.R PCA_example.csv 4 10 T

You’ll see the 3rd input is name of the input data set. Currently I require the data set to be in csv format. Rows are genes and columns are samples. Row names and column names are required in the input file.

The 4th input is number of PCs to output (define it as k, default k =5)

The 5th input is the lower limit of detection threshold (default is 0). Genes with max expression below this threshold will be removed.

The 6th input should be T or F. If T is specified, median-by-ratio normalization will be performed prior to PCA analysis (default is F).

To simplify your input command and use the default values, you may run
Rscript PCA.R PCA_example.csv

or

Rscript PCA.R PCA_example.csv 4


The code will first rescale the data to gene specific z scores, then perform PCA.
The output files are:

- prefix_PC_pairs.pdf:
pairwise plots of the transformed data; k PCs will be shown
- prefix_loading.csv
: gene loading for the top k PCs
- prefix_sort_by_absloading.csv
: For each of the top k PCs, genes are sorted by their absolute loadings in each PC
- prefix_perc_sdev.csv
: Percentage of SD explained by each PC

The ‘prefix’ is defined as the filename of the input file (the string before ‘.csv’) 


I’m currently using a quick and dirty way to keep the R open when browsing the results - the code will terminate the R process after 10^30 sec.
