#!/bin/bash


set -euo pipefail

#download the GFF3 file containing genes and TEs
wget -P ../data/Fig2/ https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes_transposons.gff

#remove the third column that are transposable_element or transposon_fragment and create raw_gene_annotation_v1
grep -v "ID=AT[1-5]TE\|Parent=AT[1-5]TE" ../data/Fig2/TAIR10_GFF3_genes_transposons.gff | grep -vE "chromosome|ChrC|ChrM" > ../data/Fig2/raw_gene_annotation_v1

#create the list of transposable element genes
grep "transposable_element_gene" ../data/Fig2/raw_gene_annotation_v1 | awk '{print$9}' | awk '{print substr($1,4,9)}' > ../data/Fig2/TE_gene_list



#create the file only with the information of transposable element genes
for i in $(cat ../data/Fig2/TE_gene_list)
do

echo ${i}

grep ${i} ../data/Fig2/raw_gene_annotation_v1 >> ../data/Fig2/TE_gene_info

done




#remove the information of transposable element genes from raw_gene_annotation_v1 to create TAIR10_GFF3_genes_noTEgenes_nochromosome.gff
grep -vf ../data/Fig2/TE_gene_info ../data/Fig2/raw_gene_annotation_v1 > ../data/Fig2/TAIR10_GFF3_genes_noTEgenes_nochromosome.gff

#extract protein coding genes from TAIR10_GFF3_genes_noTEgenes_nochromosome.gff
cat ../data/Fig2/TAIR10_GFF3_genes_noTEgenes_nochromosome.gff | grep "gene" | grep "protein" | awk '{print$1"\t"$4-1"\t"$5"\t""protein_coding_gene""\t"$7"\t"substr($9, 4, 9)}' > ../data/Fig2/TAIR10_protein_coding_genes_bed


#use bedtools sort to sort the file 
cat ../data/Fig2/TAIR10_protein_coding_genes_bed | bedtools sort  > ../data/Fig2/TAIR10_protein_coding_genes_bed_sorted


#use bedtools merge to produce the boundary of each non-overlapping or overlapping genes
bedtools merge -i ../data/Fig2/TAIR10_protein_coding_genes_bed_sorted -c 6 -o collapse -delim "|" > ../data/Fig2/TAIR10_protein_coding_genes_bed_sorted_merge




##The code below are used for producing data for Figure 1
#create TSS file
cat ../data/Fig2/TAIR10_protein_coding_genes_bed_sorted | awk -v OFS="\t" '{if($5 == "+") print $1,$2,$2+1,"TSS",$5,$6; else print $1,$3-1,$3,"TSS",$5,$6}' > ../data/Fig1/TAIR10_protein_coding_genes_TSS_bed_sorted

#create the file of protein coding genes without the ID name 
cat ../data/Fig2/TAIR10_protein_coding_genes_bed_sorted_merge | awk '{print$1"\t"$2"\t"$3}' > ../data/Fig1/TAIR10_protein_coding_genes_bed_sorted_merge_noname

#create the file of Transposable elements
cat -n ../data/Fig2/TAIR10_GFF3_genes_transposons.gff | grep "ID=AT[1-5]TE\|Parent=AT[1-5]TE" | grep "transposable_element" | awk '{print$2"\t"$5-1"\t"$6"\t""transposable_element"}' | bedtools sort > ../data/Fig1/TAIR10_transposable_elements_sorted
bedtools merge -i ../data/Fig1/TAIR10_transposable_elements_sorted | awk '{print$1"\t"$2"\t"$3"\t""TE"}' > ../data/Fig1/TAIR10_transposable_elements_sorted_merge

