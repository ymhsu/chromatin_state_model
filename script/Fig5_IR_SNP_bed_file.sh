#!/bin/bash
# Program: bedtools
# generate modified 10-state segments and SNP density and IR table for the CO rate prediction (density table)
#SBATCH --job-name=bedtools

set -euo pipefail

#extract syntenic regions between Col and Ler in Arabidopsis genome
bedtools subtract -a ../data/Fig2/Ara_genome_bed -b ../data/Fig2/SV_raw > ../data/Fig5/Ara_genome_subtract_SV_bed

#assign the information of 9 states to syntenic regions
bedtools intersect -a ../data/Fig5/state_9_total_modified -b ../data/Fig5/Ara_genome_subtract_SV_bed > ../data/Fig5/Ara_genome_subtract_SV_state_bed

#create the modified version of 10-state segments
cat ../data/Fig2/SV_raw ../data/Fig5/Ara_genome_subtract_SV_state_bed | bedtools sort | awk -v OFS="\t" '{print$1, $2, $3, $4, NR}' > ../data/Fig5/state_10_raw_state8_m_order

#get the intersection of IR and 10-state segment, add the information of IR (size and transcription direction) on segments intersected with IR
bedtools intersect -a ../data/Fig5/state_10_raw_state8_m_order -b ../data/Fig3/TAIR10_protein_coding_genes_IR_bed_trimmed_f -wb | awk -v OFS="\t" '{print$1, $2, $3, $4, $8-$7, $11}' > ../data/Fig5/state_10_Col_Ler_IR_intersect

#add the information of IR on segments not intersected with IR
bedtools subtract -a ../data/Fig5/state_10_raw_state8_m_order -b ../data/Fig5/state_10_Col_Ler_IR_intersect | awk -v OFS="\t" '{print$1, $2, $3, $4, 0, "no"}' > ../data/Fig5/state_10_Col_Ler_IR_subtract

#merge two 10-state segments into the final bed file for the following analysis
cat ../data/Fig5/state_10_Col_Ler_IR_intersect ../data/Fig5/state_10_Col_Ler_IR_subtract | bedtools sort > ../data/Fig5/state_10_Col_Ler_IR_total

#let state_10_Col_Ler_IR_total intersect with the initial bed file, then intersect with SNP information 

for bin in 0_5k 1k 2k 3k 4k 5k 10k 15k 20k 50k 100k 200k 500k 1000k
do

echo ${bin}
bedtools intersect -a ../data/Fig5/state_10_Col_Ler_IR_total -b ../data/den_table_${bin}_bed -wb | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8"\t"$9"\t"$10}' > ../data/Fig5/den_table_${bin}_state_10_Ler_IR_intersect

bedtools intersect -a ../data/Fig4/Ian_pop_passed_SNP_bed -b ../data/Fig5/den_table_${bin}_state_10_Ler_IR_intersect -wb | awk '{print$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$4}' | grep "SNP_Col_Ler" > ../data/Fig5/den_table_${bin}_state_10_Ler_IR_ISNP_intersect

bedtools intersect -a ../data/den_table_${bin}_bed -b ../data/Fig1/R_CO_final_mid_bed -wb | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > ../data/Fig1/den_table_${bin}_RCO_raw_mid_bed

done

