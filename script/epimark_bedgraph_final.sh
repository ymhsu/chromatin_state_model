#!/bin/bash


# wget to download histone modification marks (except for two data from E-MTAB-7370)
#histone modification data
set -euo pipefail

for i in $(cat ../data/Fig1/epimark_data/epimark_download_list)
do

echo ${i}
wget -P ../data/Fig1/epimark_data/ ${i}


done


#rename the downloaded file
for n in $(cat ../data/Fig1/epimark_data/epimark_rename_list)

do

    first=$(echo $n | cut -d "|" -f 1)

    second=$(echo $n | cut -d "|" -f 2)

    mv ../data/Fig1/epimark_data/${first} ../data/Fig1/epimark_data/${second}


done

#create the list of bw file for the following transformation
ls ../data/Fig1/epimark_data/ | grep "\.bw" > ../data/Fig1/epimark_data/epimark_bw_list


#use bigWigToBedGraph to transform bw file to raw bedgraph files

for i in $(cat ../data/Fig1/epimark_data/epimark_bw_list)
do

echo ${i}
echo ${i} | sed 's/\.bw//g' > ../data/Fig1/epimark_data/fqin
cat ../data/Fig1/epimark_data/fqin

#Now we need to use bigWigToBedGraph to transform bw file to bedgraph files
bigWigToBedGraph ../data/Fig1/epimark_data/${i} >(gzip > ../data/Fig1/epimark_data/$(cat ../data/Fig1/epimark_data/fqin).bedgraph.gz)

rm ../data/Fig1/epimark_data/${i}

done

rm ../data/Fig1/epimark_data/fqin
ls ../data/Fig1/epimark_data/ | grep "raw\.bedgraph\.gz" > ../data/Fig1/epimark_data/epimark_raw_bedgraph_list


#normalize these raw bedgraph file into the same bed format

for i in $(cat ../data/Fig1/epimark_data/epimark_raw_bedgraph_list)

do

echo ${i}
gunzip -c ../data/Fig1/epimark_data/${i} | awk '{if(length($4) == 0) print$1"\t"$2"\t"$3"\t"1; else print$1"\t"$2"\t"$3"\t"$4}' | awk '{if(length($1)==1) print"Chr"$1"\t"$2"\t"$3"\t"$4; else print"Chr"substr($1,4)"\t"$2"\t"$3"\t"$4}' | awk '/Chr[1-5]/ {print$0}' | bedtools sort | gzip > ../data/Fig1/epimark_data/$(echo ${i} | sed 's/_raw//g')
rm ../data/Fig1/epimark_data/${i}

done


#download two E-MTAB-7370 data and create their bedgraph file (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7370/files/processed/?ref=E-MTAB-7370)
#Two Chip-Seq data of unopened flower 
for i in 11 13 14
do

echo ${i}
wget -P ../data/Fig1/epimark_data/ https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7370/E-MTAB-7370.processed.${i}.zip

done

for i in 11 13 14
do

echo ${i}
unzip ../data/Fig1/epimark_data/E-MTAB-7370.processed.${i}.zip -d ../data/Fig1/epimark_data
rm ../data/Fig1/epimark_data/E-MTAB-7370.processed.${i}.zip

done



for n in $(cat ../data/Fig1/epimark_data/E-MTAB-7370_list)

do

first=$(echo $n | cut -d "|" -f 1)

second=$(echo $n | cut -d "|" -f 2)


##create bedgraph file with NR for the latter use of join
cat -n ../data/Fig1/epimark_data/${first} | grep "Chr" | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$5}' > ../data/Fig1/epimark_data/${first}_NR

##create input file with NR for the latter use of join (I guess the input file named as H3K9me2, so this file can be used for other histone modification marks)
cat -n ../data/Fig1/epimark_data/wt_H3K9me2_Rep1_input_libsize_norm_cov.bedGraph | grep "Chr" | awk '{print$1"\t"$5}' > ../data/Fig1/epimark_data/wt_H3K9me2_Rep1_input_libsize_norm_cov.bedGraph_NR


##use join to two abovementioned files, and calculate the ratio between Chip and input
join -1 1 -2 1 ../data/Fig1/epimark_data/${first}_NR ../data/Fig1/epimark_data/wt_H3K9me2_Rep1_input_libsize_norm_cov.bedGraph_NR | awk '{if($6 == 0) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$5; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$5/$6}' | awk '{print$2"\t"$3"\t"$4"\t"$7}' | gzip > ../data/Fig1/epimark_data/${second}

done

rm ../data/Fig1/epimark_data/wt_*_Rep1_*
