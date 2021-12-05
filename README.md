All of the analysis are performed on Linux system (Ubuntu 18.04.6 LTS).

Before performing any analysis, one needs to download "bedtools" and "bigWigToBedGraph" for the analysis.

The version of bedtools is v2.26.0.

bigWigToBedGraph can be downloaded by using the command line below:

conda install -c bioconda ucsc-bigwigtobedgraph



All the R codes run on R Studio version 1.1.463 with R version 3.6.3. 

In total, there are 6 R codes used for performing the analysis published in the paper,

and each code can be operated independently. 

Besides, there are 5 shell scripts used to produce necessary files for running R codes,

and we documented the usage of these shell scripts and a few command lines based on tools under Linux system in R codes. 


