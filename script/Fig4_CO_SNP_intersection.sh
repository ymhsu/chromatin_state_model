#!/bin/bash
# Program: bedtools intersect
# make the intersection between CO data of Rowan and Ian with initial 100-kb bed files 
#SBATCH --job-name=bedtools

set -euo pipefail




#create the intersection between Rowans' COs and the 100-bin bed file
bedtools intersect -a ../data/den_table_100k_bed -b ../data/Fig1/R_CO_final_mid_bed -wb | awk '{print$1"\t"$2"\t"$3"\t"$4}' > ../data/Fig4/den_table_100k_RCO_raw_mid_bed

#create the intersection between Blackwells' COs and the 100-bin bed file
bedtools intersect -a ../data/den_table_100k_bed -b ../data/Fig4/Ian_CO_mid_5dataset -wb | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$9}' > ../data/Fig4/den_table_100k_ICO_raw_mid_bed

#create the intersection between Blackwells' SNPs and the 100-bin bed file
bedtools intersect -a ../data/Fig4/Ian_pop_passed_SNP_bed -b ../data/den_table_100k_bed -wb | awk '{print$8"\t"$4}' > ../data/Fig4/den_table_100k_ISNP_intersect



