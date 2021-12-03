#!/bin/bash
# Program:
# Bedtools intersect to take the information of 9 features including genomic/epigenomic features into the initial bed files using different bin sizes
#SBATCH --job-name=bedtools

#Bedtools intersect
#gene/TEs data and density table
set -euo pipefail




for bin in 0_5k 1k 2k 3k 4k 5k 10k 15k 20k 50k 100k 200k 500k 1000k
do

echo ${bin}

bedtools intersect -a ../data/Fig1/den_table_${bin}_bed -b ../data/Fig1/R_CO_final_bed -wb | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > ../data/Fig1/den_table_${bin}_RCO_raw_bed
bedtools intersect -a ../data/Fig1/den_table_${bin}_bed -b ../data/Fig1/R_CO_final_mid_bed -wb | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > ../data/Fig1/den_table_${bin}_RCO_raw_mid_bed
bedtools intersect -a ../data/Fig1/epimark_data/GSM3674617_H3K27me3_leaf_R1.bedgraph.gz -b ../data/Fig1/den_table_${bin}_bed -wb | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$8}' > ../data/Fig1/den_table_${bin}_H3K27me3_SV_intersect
bedtools intersect -a ../data/Fig1/epimark_data/GSM3674621_H3K4me1_leaf_R1.bedgraph.gz -b ../data/Fig1/den_table_${bin}_bed -wb | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$8}' > ../data/Fig1/den_table_${bin}_H3K4me1_SV_intersect
bedtools intersect -a ../data/Fig1/epimark_data/GSM3674620_H3K4me3_leaf_R1.bedgraph.gz -b ../data/Fig1/den_table_${bin}_bed -wb | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$8}' > ../data/Fig1/den_table_${bin}_H3K4me3_SV_intersect
bedtools intersect -a ../data/Fig1/epimark_data/GSM4734580_H3K9me2_leaf_R1.bedgraph.gz -b ../data/Fig1/den_table_${bin}_bed -wb | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$8}' > ../data/Fig1/den_table_${bin}_H3K9me2_SV_intersect
bedtools intersect -a ../data/Fig1/epimark_data/GSM1289358_DNase_7d_seedling_R1.bedgraph.gz -b ../data/Fig1/den_table_${bin}_bed -wb | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$8}' > ../data/Fig1/den_table_${bin}_DNase_SV_intersect
bedtools intersect -a ../data/Fig1/TAIR10_protein_coding_genes_bed_sorted_merge_noname -b ../data/Fig1/den_table_${bin}_bed -wb | awk '{print$1"\t"$2"\t"$3"\t"$7}' > ../data/Fig1/den_table_${bin}_gene_SV_intersect
bedtools intersect -a ../data/Fig1/TAIR10_protein_coding_genes_TSS_bed_sorted -b ../data/Fig1/den_table_${bin}_bed -wb > ../data/Fig1/den_table_${bin}_TSS_SV_raw_bed
bedtools intersect -a ../data/Fig1/TAIR10_transposable_elements_sorted_merge -b ../data/Fig1/den_table_${bin}_bed -wb | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$8}' > ../data/Fig1/den_table_${bin}_TE_SV_intersect
bedtools intersect -a ../data/Fig1/epimark_data/GSM3674715_ATAC_leaf_R1.bedgraph.gz -b ../data/Fig1/den_table_${bin}_bed -wb | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$8}' > ../data/Fig1/den_table_${bin}_ATAC_SV_intersect

#used for testing SNP effect
bedtools intersect -a ../data/Fig5/den_table_100k_bed -b ../data/Fig4/Ian_CO_mid_5dataset -wb | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$9}' > ../data/Fig4/den_table_${bin}_ICO_raw_mid_bed

done

for i in $(cat ../data/Fig1/epimark_data/SP_Fig1_intersection_list)
do

echo ${i}

bedtools intersect -a ../data/Fig1/epimark_data/${i} -b ../data/Fig1/den_table_100k_bed -wb > ../data/Fig1/den_table_100k_$(echo ${i} | sed -E 's/\..*\.gz//g')_SV_intersect

done
