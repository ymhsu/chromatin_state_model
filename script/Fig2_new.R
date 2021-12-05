#using "p_load" from the package "pacman" to install and load necessary packages
install.packages("pacman")
library(pacman)

Packages <- c("scales", "tidyverse", "ggrepel", "ggsci", "ggpubr", "doMC", "doParallel", "foreach", "slider", "cowplot", "combinat")
p_load(Packages, character.only = TRUE)

#lapply(Packages, library, character.only = TRUE)

#change the directory "chromatin_state_model" as the working directory (the link below is an example)
setwd("/data/projects/thesis/INRA_project/Ara_TE_task/R_markdown/Model_1st/chromatin_state_model/")

#######The analysis of Fig 2A#######
#pie chart for the comparison between states and CO
#import the state data
chromatin_state_total_SV <- read_delim("./data/Fig2/chromatin_state_total_SV", delim = "\t", col_names = c("Chr", "str", "end", "state"))

#sum of bp of states in the whole genome
state_bp_sum <- chromatin_state_total_SV %>%  
  group_by(state) %>%
  summarise(sum_bp = sum(end-str)) %>%
  ungroup() %>%
  mutate(sum_total = sum(sum_bp), frac = sum_bp/sum_total)

#import CO in 9 states
chromatin_state_total_R_CO_raw_noSV <- read_delim("./data/Fig2/chromatin_state_total_noSV_R_CO", delim = "\t", col_names = c("Chr", "str", "end", "state", "Chr_CO","str_CO", "end_CO", "Sel_420", "CO_l"))

#calculate sum of CO in 9 states
state9_CO_sum <- chromatin_state_total_R_CO_raw_noSV %>%
  mutate(CO_n = (end-str)/(end_CO-str_CO)) %>%
  group_by(state) %>%
  summarise(sum_CO = sum(CO_n)) 

#import CO intersected with SVs
SV_raw_R_CO <- read_delim("./data/Fig2/SV_raw_R_CO", delim = "\t", col_names = c("Chr", "str", "end", "state", "Chr_CO","str_CO", "end_CO", "Sel_420", "CO_l"))

#calculate sum of CO located in SV
SV_CO_sum <- SV_raw_R_CO %>%
  mutate(CO_n = (end-str)/(end_CO-str_CO)) %>%
  group_by(state) %>%
  summarise(sum_CO = sum(CO_n)) 

#create the genome-wide and CO fraction in 10 states
state10_CO <- tibble(
  state = c(state9_CO_sum$state, "SV"),
  sum_CO = c(state9_CO_sum$sum_CO, SV_CO_sum$sum_CO),
  sum_bp = c(state_bp_sum$sum_bp), 
  ratio_bp = state_bp_sum$frac
) %>%
  mutate(CO_frac = sum_CO/sum(sum_CO)) %>%
  mutate(the_genome_fraction = ratio_bp*100, CO_fraction = CO_frac*100) %>%
  select(state, the_genome_fraction, CO_fraction) %>%
  gather(key = type, value = fraction, c("the_genome_fraction", "CO_fraction")) %>%
  mutate(percent_fraction = percent(fraction/100, accuracy = 0.2)) %>%
  group_by(type) %>%
  mutate(cul_frac = rev(cumsum(rev(fraction)))) %>%
  mutate(half_frac = fraction/2) %>%
  mutate(posi_y = cul_frac-fraction+half_frac) %>%
  select(-cul_frac, -half_frac) %>%
  ungroup() %>%
  mutate(label_pie = str_c(state, ": ", percent_fraction)) 

#assign factors for two kinds of fractions for the following figure order
state10_CO$type <- factor(state10_CO$type, levels = c("the_genome_fraction", "CO_fraction"))

state10_CO <- state10_CO %>%
  split(.$type)

#https://stackoverflow.com/questions/46277894/ggplot2-pie-plot-with-geom-text-repel (assign the position of legend in bar plot first, then assign them in pie chart)
#https://datavizpyr.com/remove-border-facet-in-ggplot2/ (remove spacing between grouped figures)
#https://github.com/kassambara/ggpubr/issues/235 (assign margin of the figure before using merging them)
#https://ggplot2.tidyverse.org/reference/element.html

Fig2_2A <- vector(mode = "list", length = 2)

for (i in seq_along(1:2)) {
  Fig2_2A[[i]] <- state10_CO[[i]] %>%
    ggplot(aes(x="", y=fraction, fill=state))+
    geom_bar(width = 1, stat = "identity", show.legend = FALSE) + 
    coord_polar("y", start=0) +
    theme(axis.text.x=element_blank()) +
    geom_label_repel(aes(label = percent_fraction, y = posi_y), fontface = "bold", size=6, show.legend = F, min.segment.length = 0.1, seed = 42, box.padding = 0.5, nudge_x = 0, nudge_y = 0) +
    scale_fill_manual(values=c(pal_npg("nrc", alpha = 0.7)(10))) +
    facet_wrap(~type, nrow = 1) +
    labs(x = NULL, y = NULL, fill = NULL) + theme_bw() + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), strip.text.x = element_text(colour = "black", face = "bold", size = 20), legend.title = element_blank(), plot.margin = margin(0,0,0,3, "cm"))
  
}

#######The analysis for Fig 2B#######
#The objective of this analysis is to calculate recombination rate in gene bodies and flanking regions of genes
#For flanking regions, we only focused on regions not being polluted by any genes.
#To simply our anysis, we just extracted non-overlapping genes located in syntenic regions.
#However, we still need to extract all protein coding genes to precisely locate the 3-kb flanking regions before removing genes in which we are not interested
#Open the terminal, run the shell script "genomic_features.sh" in the directory "script/" to create the bed file of protein coding genes.
TAIR10_protein_coding_genes_bed_sorted <- read_delim("./data/Fig2/TAIR10_protein_coding_genes_bed_sorted", delim = "\t", col_names = c("Chr", "str", "end", "type", "strand", "gene_name"))
TAIR10_protein_coding_genes_bed_sorted_merge <- read_delim("./data/Fig2/TAIR10_protein_coding_genes_bed_sorted_merge", delim = "\t", col_names = c("Chr", "str", "end", "gene_name"))
Ara_genome <- read_delim("./data/Fig2/Ara_genome", delim = "\t", col_names = c("Chr", "chr_end"))

#produce the strand file of each genes (for later defining TSS and TTS)
strand <- TAIR10_protein_coding_genes_bed_sorted %>%
  select(gene_name, strand)

#define 3-kb border on two sides of all protein coding genes (overlapping genes already merged)
TAIR10_protein_coding_genes_bed_light_merge_b_raw <- TAIR10_protein_coding_genes_bed_sorted_merge %>%
  left_join(strand) %>%
  #overlapping genes will not be joined
  replace_na(list(strand = "no")) %>%
  arrange(Chr, str) %>%
  group_by(Chr) %>%
  mutate(pre_end = lag(end)) %>%
  replace_na(list(pre_end=0)) %>%
  #filter str - pre_end means the overlapped genes do not have borders
  filter(str - pre_end > 0) %>%
  mutate(exten_l = if_else(str - pre_end > 3000, str - 3000, str - round((str - pre_end)/2))) %>%
  mutate(exten_r = lead(exten_l)) %>%
  mutate(exten_r = if_else(is.na(exten_r) == TRUE, end + 3000, exten_r)) %>%
  mutate(exten_r = if_else(exten_r - end > 3000, end + 3000, exten_r)) %>%
  select(-pre_end) %>%
  select(Chr, exten_l_str = exten_l, exten_l_end = str, exten_r_str = end, exten_r_end = exten_r, gene_name, strand) %>%
  gather(c(exten_l_str, exten_r_str), key = border_type_str, value = str) %>%
  gather(c(2, 3), key = border_type_end, value = end) %>%
  filter(str_sub(border_type_str, 7, 7) == str_sub(border_type_end, 7, 7)) %>%
  arrange(Chr, str) %>%
  ungroup() %>%
  group_by(gene_name) %>%
  mutate(type = if_else(strand == "no", "flanking", if_else(strand == "+" & str == min(str) | strand == "-" & end == max(end), "TSS_3kb", "TTS_3kb"))) %>%
  ungroup() %>%
  filter(str != end) %>%
  select(Chr, str, end, gene_name, strand, type) %>% left_join(Ara_genome) %>%
  filter(str <= chr_end) %>%
  mutate(end = if_else(end>chr_end, chr_end, end)) %>%
  select(-chr_end) %>%
  #only non-overlapping genes are kept
  filter(strand != "no") 
  
#create non-overlapping genes for identifying genes which are not overlapping SVs
TAIR10_protein_coding_genes_bed_non_overlapping <- TAIR10_protein_coding_genes_bed_sorted_merge %>%
  left_join(strand) %>%
  drop_na() %>%
  mutate(type = "protein_coding_genes")

write_delim(TAIR10_protein_coding_genes_bed_non_overlapping, "./data/Fig2/TAIR10_protein_coding_genes_bed_non_overlapping", delim = "\t", col_names = FALSE)

#Open the terminal, run the command below in the directory "data/Fig2/" to identify protein-coding genes not overlapping SVs. 
#bedtools subtract -a TAIR10_protein_coding_genes_bed_non_overlapping -b SV_raw > TAIR10_protein_coding_genes_bed_non_overlapping_noSV

#import non_overlapping genes that already subtract SVs 
TAIR10_protein_coding_genes_bed_non_overlapping_noSV <- read_delim("./data/Fig2/TAIR10_protein_coding_genes_bed_non_overlapping_noSV", delim = "\t", col_names = c("Chr", "str", "end", "gene_name", "strand", "type")) %>%
  mutate(status = "noSV")

#keep intact non-overlapping genes without being polluted by SVs (larger than 100 bps)
TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact <- TAIR10_protein_coding_genes_bed_non_overlapping %>%
  left_join(TAIR10_protein_coding_genes_bed_non_overlapping_noSV) %>%
  drop_na() %>%
  filter(end - str >= 100)

#create the list of no-SV intact genes for extracting the same gene list of flanking regions
TAIR10_protein_coding_genes_noSV_intact_filter <- TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact %>%
  select(gene_name, status)

#create the bed file of 3-kb flanking regions
TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_b_raw <- TAIR10_protein_coding_genes_bed_light_merge_b_raw %>%
  left_join(TAIR10_protein_coding_genes_noSV_intact_filter) %>%
  drop_na()

#create the bed file of flanking regions and gene bodies of genes for removing SVs of flanking regions later
TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_integral <- bind_rows(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact, TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_b_raw) %>%
  arrange(Chr, str)

write_delim(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_integral, "./data/Fig2/TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_integral", col_names = FALSE, delim = "\t")


#Open the terminal, run the command below in the directory "data/Fig2/" to remove SVs of flanking regions 
#bedtools subtract -a TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_integral -b SV_raw > TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_integral_flanking_noSV
#Then, create the final integral list of syntenic genes with gene bodies and flanking regions

#create the gene list of gene bodies and flanking regions with the old coordinates
pc_gene_list <- TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_integral %>%
  select(gene_name, type, str_old = str, end_old = end)

#load the gene list that already subtracted SVs
TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_integral_flanking_noSV <- read_delim("./data/Fig2/TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_integral_flanking_noSV", delim = "\t", col_names = c("Chr", "str", "end", "gene_name", "strand", "type", "status"))

#flanking regions which are not broken by SVs at the two boundaries will be kept (the if_else part below)
#and the coordinates of kept regions will be modified if necessary 
TAIR10_gene_noSV_non_overlapping_3kb_integral_temp <- TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_integral_flanking_noSV %>%
  left_join(pc_gene_list) %>%
  mutate(filter_label = if_else(type == "protein_coding_genes", T, 
                                if_else(strand == "+" & type == "TSS_3kb" & end == end_old, T, 
                                        if_else(strand == "+" & type == "TTS_3kb" & str == str_old, T,
                                                if_else(strand == "-" & type == "TTS_3kb" & end == end_old, T,
                                                        if_else(strand == "-" & type == "TSS_3kb" & str == str_old, T, F)))))) %>%
  mutate(SV_label = if_else(filter_label == FALSE, 1, 0)) %>%
  filter(SV_label != 1) %>%
  select(gene_name, strand, type, str_integral = str, end_integral = end)


#The abovementioned part is to produce the borders of each syntenic intact genes
#Afterward, we create bins of flanking regions and gene bodies of concerned genes. We forced gene bodies to have 100 bins, and divided each flanking regions into bins by every 100 bps.
#The strategy below is to naively create 3-kb borders for each no-SV intact genes at the beginning, then use the integral file above to correct the binning file 


#create the list of the concerned genes
TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l <- TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact %>%
  split(.$gene_name)

#create border of each gene (TSS/TTS) (the whole gene set)
TAIR10_GFF3_genes_intact_l_b <- vector("list", length = length(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l))

for (i in seq_along(TAIR10_GFF3_genes_intact_l_b)) {
  posi_b_l <- seq(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l[[i]]$str-3000, TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l[[i]]$str, length.out = 31)
  posi_b_r <- seq(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l[[i]]$end, TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l[[i]]$end + 3000, length.out = 31)
  
  TAIR10_GFF3_genes_intact_l_b[[i]] <- tibble(
    Chr = rep(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l[[i]]$Chr, 150),
    str = c(posi_b_l[11:30], posi_b_r[11:30], posi_b_l[6:30], posi_b_r[6:30], posi_b_l[1:30], posi_b_r[1:30]),
    end = c(posi_b_l[12:31], posi_b_r[12:31], posi_b_l[7:31], posi_b_r[7:31], posi_b_l[2:31], posi_b_r[2:31]),
    gene_name = rep(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l[[i]]$gene_name, 150),
    strand = rep(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l[[i]]$strand, 150)
  ) %>%
    mutate(type = if_else(strand == "+", c(rep("TSS_2kb", 20), rep("TTS_2kb", 20), rep("TSS_2.5kb", 25), rep("TTS_2.5kb", 25), rep("TSS_3kb", 30), rep("TTS_3kb", 30)), c(rep("TTS_2kb", 20), rep("TSS_2kb", 20), rep("TTS_2.5kb", 25), rep("TSS_2.5kb", 25), rep("TTS_3kb", 30), rep("TSS_3kb", 30)))) %>%
    mutate(order = if_else(strand == "+", c(rep(seq(1:20), 2), rep(seq(1:25), 2), rep(seq(1:30), 2)), c(rep(rev(1:20), 2), rep(rev(1:25), 2), rep(rev(1:30), 2)))) %>%
    select(Chr, str, end, gene_name, type, strand, order) %>%
    mutate(type_size = str_sub(type, 5, 10)) %>%
    split(.$type_size) %>%
    map(. %>% select(-type_size))
}


#generate the list of gene body
TAIR10_GFF3_genes_intact_l_g <- vector("list", length = length(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l))

for (i in seq_along(TAIR10_GFF3_genes_intact_l_b)) {
  posi_g  <- floor(seq(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l[[i]]$str, TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l[[i]]$end, length.out = 101))
  
  TAIR10_GFF3_genes_intact_l_g[[i]] <-tibble(
    Chr = rep(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l[[i]]$Chr, 100),
    str = posi_g[1:100],
    end = posi_g[2:101],
    gene_name = rep(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l[[i]]$gene_name, 100),
    type = "protein_coding_genes",
    strand = rep(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l[[i]]$strand, 100)
  ) %>%
    mutate(order = if_else(strand == "+", seq(1:100), rev(1:100)))
  }


#combine the gene body and the extension together
TAIR10_GFF3_genes_intact_noSV_body_border_2kb <- vector("list", length = length(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l))
TAIR10_GFF3_genes_intact_noSV_body_border_2_5kb <- vector("list", length = length(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l))
TAIR10_GFF3_genes_intact_noSV_body_border_3kb <- vector("list", length = length(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l))

for (i in seq_along(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_l)) {
  TAIR10_GFF3_genes_intact_noSV_body_border_2kb[[i]] <- bind_rows(TAIR10_GFF3_genes_intact_l_b[[i]]$`2kb`, TAIR10_GFF3_genes_intact_l_g[[i]]) %>%
    arrange(Chr, str) 
  
  TAIR10_GFF3_genes_intact_noSV_body_border_2_5kb[[i]] <- bind_rows(TAIR10_GFF3_genes_intact_l_b[[i]]$`2.5kb`, TAIR10_GFF3_genes_intact_l_g[[i]]) %>%
    arrange(Chr, str) 
  
  TAIR10_GFF3_genes_intact_noSV_body_border_3kb[[i]] <- bind_rows(TAIR10_GFF3_genes_intact_l_b[[i]]$`3kb`, TAIR10_GFF3_genes_intact_l_g[[i]]) %>%
    arrange(Chr, str)
}


#modify the last bin using the chromosome length of Arabidopsis
TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total <- list(bind_rows(TAIR10_GFF3_genes_intact_noSV_body_border_2kb), bind_rows(TAIR10_GFF3_genes_intact_noSV_body_border_2_5kb), bind_rows(TAIR10_GFF3_genes_intact_noSV_body_border_3kb)) %>%
  map(. %>% filter(str >= 0 & end >= 0)) %>%
  map(. %>% left_join(Ara_genome)) %>%
  map(. %>% filter(str <= chr_end)) %>%
  map(. %>% mutate(end = if_else(end>chr_end, chr_end, end))) %>%
  map(. %>% select(-chr_end))


#removed bins not in the precise boundary of flanking regions, and modified the coordinates of bins if necessary
TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total_f <- TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total[[3]] %>%
  left_join(TAIR10_gene_noSV_non_overlapping_3kb_integral_temp) %>%
  drop_na() %>%
  filter(str < end_integral) %>%
  filter(end > str_integral) %>%
  mutate(str_n = if_else(str < str_integral, str_integral, str), end_n = if_else(end > end_integral, end_integral, end)) %>%
  select(-str, -end) %>%
  select(Chr, str = str_n, end = end_n, gene_name, strand, type, order) 



write_delim(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total_f, "./data/Fig2/TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total_f", delim = "\t", col_names = FALSE)

#Open the terminal, run the command below in the directory "data/Fig2/" to create the file of COs intersecting the final binned gene file
#bedtools intersect -a TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total_f -b R_CO_final_bed -wb > TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total_RCO

#load the file of "TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total_RCO", and calculate sum of COs
TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total_RCO_sum <- read_delim("./data/Fig2/TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total_RCO", delim = "\t", col_names = c("Chr", "str", "end", "gene_name", "strand", "type", "order", "Chr_CO", "str_CO", "end_CO", "sel_420", "CO_label")) %>%
  mutate(CO_n = (end-str)/(end_CO-str_CO)) %>%
  group_by(gene_name, type, order) %>%
  summarise(sum_CO = sum(CO_n))

#attach the previous files with all summed CO to the list of binning file to calculate summed CO of all bins
gene_sum_CO_3kb_non_overlapped <- TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total_f %>%
  left_join(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total_RCO_sum) %>%
  replace_na(list(sum_CO = 0)) %>%
  mutate(type = if_else(type == "protein_coding_genes", type, str_sub(type, 1, 3))) %>%
  ungroup() %>%
  group_by(type, order) %>%
  summarise(sum_size = sum(end - str), sum_CO = sum(sum_CO)) %>%
  split(.$type)

#create the raw figure of CO of genes and their 3-kb flanking regions
row_n = nrow(gene_sum_CO_3kb_non_overlapped$TSS)

gene_figure_raw_3kb_non_overlapped <- tibble(
  type = c(gene_sum_CO_3kb_non_overlapped$TSS$type, gene_sum_CO_3kb_non_overlapped$protein_coding_genes$type, gene_sum_CO_3kb_non_overlapped$TTS$type), 
  bin = c(gene_sum_CO_3kb_non_overlapped$TSS$order, gene_sum_CO_3kb_non_overlapped$protein_coding_genes$order, gene_sum_CO_3kb_non_overlapped$TTS$order),
  sum_size = c(gene_sum_CO_3kb_non_overlapped$TSS$sum_size, gene_sum_CO_3kb_non_overlapped$protein_coding_genes$sum_size, gene_sum_CO_3kb_non_overlapped$TTS$sum_size),
  sum_CO = c(gene_sum_CO_3kb_non_overlapped$TSS$sum_CO, gene_sum_CO_3kb_non_overlapped$protein_coding_genes$sum_CO, gene_sum_CO_3kb_non_overlapped$TTS$sum_CO)
) %>%
  mutate(Recrate = sum_CO/2182/2/sum_size*1000000*100) %>%
  mutate(mid_posi_bin = if_else(type == "TSS", (bin + bin - 1)/2, if_else(type == "TTS", (bin + bin - 1)/2 + 100 + row_n, (bin + bin - 1)/2 + row_n))) %>%
  #group_by(type) %>%
  #summarise(Recrate = sum(sum_CO)/sum(sum_size)/2182/2*10^8)
  #View()
  ggplot(mapping = aes(mid_posi_bin, Recrate)) +
  geom_col(width = 1)+
  labs(x = "bin of genes", y = "Recombination rate (cM/Mb)")

gene_type_l_3kb_non_overlapped <- tibble(
  mid_posi_bin = c(0, row_n, row_n, row_n+100, row_n+100, 2*row_n+100),
  Recrate = rep(0, 6), 
  gene_type = c("Flanking TSS", "Flanking TSS", "gene_body", "gene_body", "Flanking TTS", "Flanking TTS")
)

#Open the terminal, run the command below in the directory "data/Fig2/" to produce state files with bins of genes
#bedtools intersect -a TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total_f -b state_10_raw -wb | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11}' > TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total_state
gene_sum_state_3kb_non_overlapped <- read_delim("./data/Fig2/TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total_state", delim = "\t", col_names = c("Chr", "str", "end", "gene_name", "strand", "type", "order", "state")) %>%
  mutate(type = if_else(type == "protein_coding_genes", type, str_sub(type, 1, 3))) %>%
  group_by(type, order, state) %>%
  summarise(sum_bp_state = sum(end - str)) %>%
  ungroup() %>%
  group_by(type, order) %>%
  mutate(sum_bp_total = sum(sum_bp_state)) %>%
  mutate(state_fraction = sum_bp_state/sum_bp_total) %>%
  mutate(bin = order) %>%
  ungroup() %>%
  mutate(state_m = str_sub(state, 6, 6)) %>%
  select(-order) %>%
  mutate(mid_posi_bin = if_else(type == "TSS", (bin + bin - 1)/2, if_else(type == "TTS", (bin + bin - 1)/2 + 100 + row_n, (bin + bin - 1)/2 + row_n)))

#using genome-wide CO rate of 9 state (no SV) to predict CO rate of each bin of gene bodies and flanking regions
predicted_bins_non_overlapped <- gene_sum_state_3kb_non_overlapped %>%
  arrange(mid_posi_bin) %>%
  left_join(experimental_rec_state) %>%
  mutate(experimental_rec_part = experimental_rec * state_fraction) %>%
  group_by(mid_posi_bin) %>%
  summarise(experimental_rec_f = sum(experimental_rec_part)) 

###produce Fig 2B###
#Fig2B top
Fig2_gene_n_top <- ggplot(gene_sum_state_3kb_non_overlapped) + geom_col(aes(x = mid_posi_bin, y = state_fraction, fill = state_m), width = 1) +
  theme_bw() +
  scale_fill_manual(values=c(pal_npg("nrc", alpha = 0.7)(10))) + 
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.position = "bottom", 
        axis.text = element_text(size = 18, face = "bold"), plot.margin = margin(0,0,0,3, "cm")) + scale_x_continuous(name = c(""), breaks=c(0, 30, 80, 130, 160), labels=c("-3 kb", "TSS", "gene bodies", "TTS", "3 kb")) +
  scale_y_continuous(name = c("The genome fraction of chromatin states"), breaks = c(0, 0.5, 1)) + guides(fill = guide_legend(nrow = 1))


#Fig2B bottom
Fig2_gene_n_bottom <- gene_figure_raw_3kb_non_overlapped +
  #geom_line(data = predicted_100kb_bins_non_overlapped, aes(mid_posi_bin, predicted_rec_f), size = 1, color = "red") + 
  geom_line(data = predicted_bins_non_overlapped, aes(mid_posi_bin, experimental_rec_f), size = 1, color = "blue") + 
  geom_line(data = gene_type_l_3kb_non_overlapped, aes(mid_posi_bin, Recrate, color = gene_type), size = 4) + theme_bw() + 
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), 
        axis.text = element_text(size = 18, face = "bold"), legend.position = c(0.5, 0.9), legend.box.background = element_rect(colour = "black"), plot.margin = margin(0,0,0,3, "cm")) + 
  scale_x_continuous(name = "", breaks=c(0, 30, 80, 130, 160), labels=c("-3 kb", "TSS", "gene bodies", "TTS", "3 kb")) +
  scale_y_continuous(breaks = seq(0, 6, 0.5)) 

######The analysis of Fig 2C######
#Open the terminal, run the command below in the directory "data/Fig2/" to create the table with all information of protein coding genes and intergenic regions
#bedtools subtract -a Ara_genome_bed -b TAIR10_protein_coding_genes_bed_sorted | awk '{print$1"\t"$2"\t"$3"\t""IR""\t""no""\t""IR"NR}' > TAIR10_protein_coding_genes_IR_only_temp_bed
#cat TAIR10_protein_coding_genes_bed_sorted TAIR10_protein_coding_genes_IR_only_temp_bed | bedtools sort > TAIR10_protein_coding_genes_IR_bed 
#load the table containing protein coding genes and intergenic regions 
TAIR10_protein_coding_genes_IR_bed <- read_delim("./data/Fig2/TAIR10_protein_coding_genes_IR_bed", delim = "\t", col_names = c("Chr", "str", "end", "type", "strand", "ID"))

#filter out IR events
TAIR10_protein_coding_genes_IR_bed_trimmed <- TAIR10_protein_coding_genes_IR_bed %>%
  group_by(Chr) %>%
  filter(str != min(str) & end != max(end)) %>%
  mutate(pre_strand = lag(strand), aft_strand = lead(strand)) %>%
  filter(type == "IR") %>%
  mutate(IR_type = if_else(pre_strand == aft_strand, "parallel", if_else(pre_strand == "+" & aft_strand == "-", "convergent", "divergent")))

write_delim(TAIR10_protein_coding_genes_IR_bed_trimmed, "./data/Fig2/TAIR10_protein_coding_genes_IR_bed_trimmed", delim = "\t", col_names = FALSE)  

#Open the terminal, run the command below in the directory "data/Fig2/" to subtract SV, then extracted intact IR in syntenic regions
#bedtools subtract -a TAIR10_protein_coding_genes_IR_bed_trimmed -b SV_raw > TAIR10_protein_coding_genes_IR_bed_trimmed_no_SV
#build up the list of IR event in syntenic regions and syntenic IR size for each IR
#we also extracted IR with SV to investigate whether "spreading" occurs in SV nearby regions

TAIR10_protein_coding_genes_IR_bed_trimmed_no_SV_sum_IR_size <- read_delim("./data/Fig2/TAIR10_protein_coding_genes_IR_bed_trimmed_no_SV", delim = "\t", col_names = c("Chr", "str", "end", "type", "strand", "ID", "pre_strand", "aft_strand", "IR_type")) %>%
  group_by(ID) %>%
  mutate(syn_IR_size = sum(end-str)) %>%
  select(ID, Chr, str, end, type, IR_type, syn_IR_size) %>%
  ungroup() %>%
  select(-ID)

write_delim(TAIR10_protein_coding_genes_IR_bed_trimmed_no_SV_sum_IR_size, "./analysis_output/TAIR10_protein_coding_genes_IR_bed_trimmed_no_SV_sum_IR_size", delim = "\t", col_names = FALSE)

#extract IR without any SV
TAIR10_protein_coding_genes_IR_bed_trimmed_no_SV_sum_IR_size %>%
  #group_by(ID) %>%
  filter(end-str == syn_IR_size)

#load IR without SV
TAIR10_protein_coding_genes_IR_bed_trimmed_no_SV <- read_delim("./data/Fig2/TAIR10_protein_coding_genes_IR_bed_trimmed_no_SV", delim = "\t", col_names = c("Chr", "str", "end", "type", "strand", "ID", "pre_strand", "aft_strand", "IR_type")) %>%
  group_by(ID) %>%
  mutate(count = n()) %>%
  filter(count == 1) %>%
  mutate(syn_size = end-str) %>%
  select(ID, syn_size, count)

#extract syntenic IR which are intact events 
TAIR10_protein_coding_genes_IR_bed_trimmed_no_SV_intact <- TAIR10_protein_coding_genes_IR_bed_trimmed %>%
  left_join(TAIR10_protein_coding_genes_IR_bed_trimmed_no_SV) %>%
  filter(end - str == syn_size)


#extracted IR with the size at least 100 bp 
TAIR10_protein_coding_genes_IR_bed_trimmed_100bp <- TAIR10_protein_coding_genes_IR_bed_trimmed %>%
  filter(end-str >= 100) %>%
  split(.$ID)


#generate bins of IR (the function)
produce_IR_100bins <- function(data) {
  posi_g  <- floor(seq(data$str, data$end, length.out = 101))
  
  x <-tibble(
    Chr = rep(data$Chr, 100),
    str = posi_g[1:100],
    end = posi_g[2:101],
    IR_name = rep(data$ID, 100),
    type = "IR",
    strand = rep(data$strand, 100),
    IR_type = rep(data$IR_type, 100),
    pre_strand = rep(data$pre_strand, 100),
    aft_strand = rep(data$aft_strand, 100)
  ) %>% 
    mutate(bin = if_else(pre_strand == "-" & aft_strand == "-", rev(c(1:n())), c(1:n())))
  x
}

##combine the information of 100 bins of every IR together, then get their intersection with Rowan's intervals and 9 states (including IR with or without SVs)
TAIR10_protein_coding_genes_IR_100bins_integral <- TAIR10_protein_coding_genes_IR_bed_trimmed_100bp %>%
  map(. %>% produce_IR_100bins()) %>%
  bind_rows() %>%
  arrange(Chr, str)

write_delim(TAIR10_protein_coding_genes_IR_100bins_integral, "./data/Fig2/TAIR10_protein_coding_genes_IR_100bins_integral", delim = "\t", col_names = FALSE)

#Open the terminal, run the command below in the directory "data/Fig2/" to get intersected CO info. 
#bedtools intersect -a TAIR10_protein_coding_genes_IR_100bins_integral -b R_CO_final_bed -wb > TAIR10_protein_coding_genes_IR_100bins_integral_RCO
TAIR10_protein_coding_genes_IR_100bins_integral_RCO_sum <- read_delim("./data/Fig2/TAIR10_protein_coding_genes_IR_100bins_integral_RCO", delim = "\t", col_names = c("Chr", "str", "end", "IR_name", "type", "strand", "IR_type", "pre_strand", "aft_strand", "bin", "Chr_CO", "str_CO", "end_CO", "sel", "label")) %>%
  mutate(CO_n_f = (end-str)/(end_CO-str_CO)) %>%
  group_by(IR_name, bin) %>%
  summarise(sum_CO = sum(CO_n_f))

#Open the terminal, run the command below in the directory "data/Fig2/" to get intersected state info. 
#bedtools intersect -a TAIR10_protein_coding_genes_IR_100bins_integral -b state_10_raw -wb | awk -v OFS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14}' | bedtools sort > TAIR10_protein_coding_genes_IR_100bins_integral_state
TAIR10_protein_coding_genes_IR_100bins_integral_state_sum <- read_delim("./data/Fig2/TAIR10_protein_coding_genes_IR_100bins_integral_state", delim = "\t", col_names = c("Chr", "str", "end", "IR_name", "type", "strand", "IR_type", "pre_strand", "aft_strand", "bin", "state")) %>%
  group_by(IR_name, bin, state) %>%
  summarise(sum_state_bp = sum(end-str))

#generate the list of 10 states of each bin of each IR for the following computation of fractions of 10 states
bin_state_list <- tibble(bin = rep(c(1:100), each = 10), 
                         state = rep(c(str_c("state", c(1:9)), "SV"), 100))


TAIR10_protein_coding_genes_IR_100bins_integral_state <-TAIR10_protein_coding_genes_IR_100bins_integral %>%
  left_join(bin_state_list)

#divide the IR data into two parts (with or without SVs)
TAIR10_protein_coding_genes_IR_bed_trimmed_no_SV_intact_list <- TAIR10_protein_coding_genes_IR_bed_trimmed_no_SV_intact %>%
  ungroup() %>%
  select(IR_name = ID) %>%
  mutate(SV_status = "no_SV")

#create the raw file for following calculation of predicted CO rate for all events or 3 types of transcriptions (IR with or without SV are separated)
TAIR10_protein_coding_genes_IR_100bins_integral_state_frac <- TAIR10_protein_coding_genes_IR_100bins_integral_state %>%
  left_join(TAIR10_protein_coding_genes_IR_bed_trimmed_no_SV_intact_list) %>%
  left_join(TAIR10_protein_coding_genes_IR_100bins_integral_state_sum) %>%
  replace_na(list(sum_state_bp = 0, SV_status = "with_SV")) %>%
  group_by(SV_status, IR_type, bin, state) %>%
  summarise(sum_state_bp = sum(sum_state_bp)) %>%
  group_by(SV_status, IR_type, bin) %>%
  mutate(sum_bin_bp = sum(sum_state_bp)) %>%
  mutate(state_fraction = sum_state_bp/sum_bin_bp) %>%
  mutate(mid_posi_bin = (bin + bin - 1)/2)


#based on the fraction of 10 states, we can calculate the theoretical recrate of each bin of integrated IR by using observed genome-wide recrate of 10 state
#calculate genome-wide recrate of 10 states using the data produced in the analysis of Fig 2A.
#the table of the genome-wide recrate of 10 states
experimental_rec_state <- state_bp_sum %>%
  left_join(state9_CO_sum) %>%
  replace_na(list(sum_CO = 0)) %>%
  mutate(sum_CO = if_else(state == "SV", 17077 - sum(sum_CO), sum_CO)) %>%
  mutate(experimental_rec = sum_CO/sum_bp/2182/2*10^8) %>%
  select(state, experimental_rec)

#the table of theoretical recrate of each bin of different IR_type based on genome-wide recrate
TAIR10_protein_coding_genes_IR_100bins_integral_state_frac_rec <- TAIR10_protein_coding_genes_IR_100bins_integral_state_frac %>%
  arrange(mid_posi_bin) %>%
  left_join(experimental_rec_state) %>%
  mutate(experimental_rec_part = experimental_rec * state_fraction) %>%
  group_by(SV_status, IR_type, mid_posi_bin) %>%
  summarise(experimental_rec_f = sum(experimental_rec_part))

#the table of theoretical recrate of each bin of different IR_type based on genome-wide recrate(with and without SV and pooling different transcription types)
TAIR10_protein_coding_genes_IR_100bins_integral_state_frac_rec_pooled <- TAIR10_protein_coding_genes_IR_100bins_integral_state_frac %>%
  #filter(bin == 1 & state == "state1")
  group_by(bin, mid_posi_bin, state) %>%
  summarise(sum_state_bp = sum(sum_state_bp)) %>%
  group_by(bin) %>%
  mutate(sum_bin_bp = sum(sum_state_bp), state_fraction = sum_state_bp/sum_bin_bp) %>%
  left_join(experimental_rec_state) %>%
  mutate(experimental_rec_part = experimental_rec * state_fraction) %>%
  group_by(bin, mid_posi_bin) %>%
  summarise(experimental_rec_f = sum(experimental_rec_part))  


##the first version of this figure##
#the comparison between the fraction of 10 states and recrate in each bin
#then we also compare the experimental recrate with theoretical recrate produced by genome-wide recrate of 10 states

#prepare CO rate data for plotting 
IR_rec_plot_d <- TAIR10_protein_coding_genes_IR_100bins_integral %>%
  left_join(TAIR10_protein_coding_genes_IR_bed_trimmed_no_SV_intact_list) %>%
  left_join(TAIR10_protein_coding_genes_IR_100bins_integral_RCO_sum) %>%
  replace_na(list(sum_CO = 0, SV_status = "with_SV")) %>%
  group_by(SV_status, IR_type, bin) %>%
  summarise(sum_CO = sum(sum_CO), sum_bp = sum(end-str)) %>%
  mutate(Recrate = sum_CO/sum_bp/2182/2*10^8, mid_posi_bin = (bin + bin - 1)/2) %>%
  left_join(TAIR10_protein_coding_genes_IR_100bins_integral_state_frac_rec) 
#filter(IR_type == "parallel") %>%

#the function of plotting rec rate
IR_Rec_function <-  function(data) {data %>% ggplot() +
    geom_col(mapping = aes(mid_posi_bin, Recrate), width = 1)+
    geom_line(aes(mid_posi_bin, experimental_rec_f), size = 1, color = "blue") + 
    theme_bw() +
    labs(x = "bin of IR", y = "Recombination rate (cM/Mb)") + 
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 14, face = "bold"),
          legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), 
          axis.text = element_text(size = 18, face = "bold"), legend.position = c(0.85, 0.9), legend.box.background = element_rect(colour = "black"), plot.margin = margin(0,0,0,0, "cm")) + 
    scale_x_continuous(name = "") +
    scale_y_continuous(breaks = seq(0, 6, 0.5)) 
}

#the function of plotting state
IR_state_function <-  function(data) {data %>%
    ggplot() + geom_col(aes(x = mid_posi_bin, y = state_fraction, fill = state_m), width = 1) +
    theme_bw() +
    scale_fill_manual(values=c(pal_npg("nrc", alpha = 0.7)(10))) + 
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 14, face = "bold"),
          legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.position = "bottom", 
          axis.text = element_text(size = 18, face = "bold"), plot.margin = margin(0,0,0,0, "cm")) + scale_x_continuous(name = c("")) +
    scale_y_continuous(name = c("The genome fraction of chromatin states"), breaks = c(0, 0.5, 1)) + guides(fill = guide_legend(nrow = 1)) 
}


#plot the CO rate and the distribution of all IR events
IR_figure_bottom <- IR_rec_plot_d %>%
  group_by(bin) %>%
  summarise(sum_CO = sum(sum_CO), sum_bp = sum(sum_bp)) %>%
  mutate(Recrate = sum_CO/sum_bp/2182/2/10^-8) %>%
  left_join(TAIR10_protein_coding_genes_IR_100bins_integral_state_frac_rec_pooled) %>%
  #filter(SV_status == "no_SV") %>%
  IR_Rec_function() +
  theme(plot.margin = margin(0,0.5,0,3, "cm"))

IR_figure_top <- TAIR10_protein_coding_genes_IR_100bins_integral_state_frac %>%
  group_by(bin, mid_posi_bin, state) %>%
  summarise(sum_state_bp = sum(sum_state_bp)) %>%
  group_by(bin) %>%
  mutate(sum_bin_bp = sum(sum_state_bp), state_fraction = sum_state_bp/sum_bin_bp) %>%
  mutate(state_m = if_else(state == "SV", "SV", str_sub(state, 6, 7))) %>%
  #filter(SV_status == "no_SV") %>%
  #filter(state != "SV") %>% 
  IR_state_function() +
  theme(plot.margin = margin(0,0.5,0,3, "cm"))

IR_figure <- ggarrange(IR_figure_top, IR_figure_bottom, ncol = 1, labels = c("C", ""), font.label = list(size = 36, color = "black", face = "bold"), label.x = 0.05)

ggsave("./analysis_output/IR_figure.jpeg", IR_figure, width = 300, height = 400, units = c("mm"), dpi = 320)


#####merge Fig 2A/B/C ########
Fig2_test_right <- ggarrange(Fig2_gene_n_top, Fig2_gene_n_bottom, ncol = 1, labels = c("B", ""), font.label = list(size = 36, color = "black", face = "bold"), label.x = 0.05)
Fig2_test_left <- ggarrange(Fig2_2A[[1]], Fig2_2A[[2]], ncol = 1, labels = c("A", ""), font.label = list(size = 36, color = "black", face = "bold"), label.x = 0.05)

Fig2_test <- ggarrange(Fig2_test_left, Fig2_test_right, IR_figure, ncol = 3, nrow = 1, widths = c(1,1.1,1.1)) + 
  theme(plot.background = element_rect(fill="white", color = NA))

ggsave("./analysis_output/Fig2_test_protein_coding_gene_v3.jpeg", Fig2_test, width = 640, height = 400, units = c("mm"), dpi = 320)



########SP figure for the chromatin state and CO rate profile in different groups of gene size####### (#SP Fig S3)
#plot the gene figure using different quantiles of gene sizes
#label protein coding genes with their categories of size
gene_size_fct <- TAIR10_gene_noSV_non_overlapping_3kb_integral_temp %>%
  filter(type == "protein_coding_genes") %>%
  mutate(size = end_integral - str_integral) %>%
  mutate(quantile_size = if_else(size <= quantile(size, 0.25), "genes ≤ 25%", 
                                 if_else(size > quantile(size, 0.25) & size <= quantile(size, 0.5), "25% < genes ≤ 50%",
                                         if_else(size > quantile(size, 0.5) & size <= quantile(size, 0.75), "50% < genes ≤ 75%", "75% < genes ≤ 100%")))) %>%
  select(gene_name, quantile_size) %>%
  mutate(quantile_size = factor(quantile_size, levels = c("genes ≤ 25%", "25% < genes ≤ 50%", "50% < genes ≤ 75%", "75% < genes ≤ 100%")))

#This is to know the median and maximum size of 4 categories of genes
TAIR10_gene_noSV_non_overlapping_3kb_integral_temp %>%
  filter(type == "protein_coding_genes") %>%
  mutate(size = end_integral - str_integral) %>%
  mutate(quantile_size = if_else(size <= quantile(size, 0.25), "genes ≤ 25%", 
                                 if_else(size > quantile(size, 0.25) & size <= quantile(size, 0.5), "25% < genes ≤ 50%",
                                         if_else(size > quantile(size, 0.5) & size <= quantile(size, 0.75), "50% < genes ≤ 75%", "75% < genes ≤ 100%")))) %>%
  group_by(quantile_size) %>%
  summarise(median_size = median(size), max_size = max(size))

#create the table with bins containing the information of states and categories of gene sizes
gene_sum_state_3kb_non_overlapped_quantile_gene_size <- 
read_delim("./data/Fig2/TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total_state", delim = "\t", col_names = c("Chr", "str", "end", "gene_name", "strand", "type", "order", "state")) %>%
  mutate(type = if_else(type == "protein_coding_genes", type, str_sub(type, 1, 3))) %>%
  left_join(gene_size_fct) %>%
  group_by(quantile_size, type, order, state) %>%
  summarise(sum_bp_state = sum(end - str)) %>%
  ungroup() %>%
  group_by(quantile_size, type, order) %>%
  mutate(sum_bp_total = sum(sum_bp_state)) %>%
  mutate(state_fraction = sum_bp_state/sum_bp_total) %>%
  mutate(bin = order) %>%
  ungroup() %>%
  select(-order) %>%
  mutate(mid_posi_bin = if_else(type == "TSS", (bin + bin - 1)/2, if_else(type == "TTS", (bin + bin - 1)/2 + 100 + row_n, (bin + bin - 1)/2 + row_n)))

#plot fractions of states of 4 types of genes
gene_size_top_state <- ggplot(gene_sum_state_3kb_non_overlapped_quantile_gene_size) + geom_col(aes(x = mid_posi_bin, y = state_fraction, fill = state), width = 1) +
  theme_bw() +
  scale_fill_manual(values=c(pal_npg("nrc", alpha = 0.7)(10))) + 
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.position = "bottom", 
        axis.text = element_text(size = 18, face = "bold"), plot.margin = margin(0,0.2,1,0.2, "cm")) + scale_x_continuous(name = c(""), breaks=c(0, 30, 80, 130, 160), labels=c("-3 kb", "TSS", "gene bodies", "TTS", "3 kb")) +
  scale_y_continuous(name = c("The genome fraction of chromatin states"), breaks = c(0, 0.5, 1)) + guides(fill = guide_legend(nrow = 1)) + facet_wrap(~quantile_size, nrow = 1)


#calculate the predicted CO rate for 4 types of gene sizes
predicted_100kb_bins_non_overlapped_quantile_gene_size <- gene_sum_state_3kb_non_overlapped_quantile_gene_size %>%
  arrange(mid_posi_bin) %>%
  left_join(experimental_rec_state) %>%
  mutate(experimental_rec_part = experimental_rec * state_fraction) %>%
  group_by(quantile_size, mid_posi_bin) %>%
  summarise(experimental_rec_f = sum(experimental_rec_part)) 


#create the table for plots of CO rate of different sizes of genes
gene_sum_CO_3kb_non_overlapped_quantile_gene_size <- TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total_f %>%
  left_join(TAIR10_protein_coding_genes_bed_non_overlapping_noSV_intact_3kb_total_RCO_sum) %>%
  replace_na(list(sum_CO = 0)) %>%
  left_join(gene_size_fct) %>%
  mutate(type = if_else(type == "protein_coding_genes", type, str_sub(type, 1, 3))) %>%
  ungroup() %>%
  group_by(quantile_size, type, order) %>%
  summarise(sum_size = sum(end - str), sum_CO = sum(sum_CO)) %>%
  mutate(Recrate = sum_CO/2182/2/sum_size*1000000*100) %>%
  mutate(mid_posi_bin = if_else(type == "TSS", (order + order - 1)/2, if_else(type == "TTS", (order + order - 1)/2 + 100 + row_n, (order + order - 1)/2 + row_n)))



#plot CO rate of 4 types of genes
gene_size_bottom_CO <- gene_sum_CO_3kb_non_overlapped_quantile_gene_size %>%
  ggplot(mapping = aes(mid_posi_bin, Recrate)) +
  geom_col(width = 1)+
  labs(x = "bin of genes", y = "Recombination rate (cM/Mb)") + facet_wrap(~quantile_size, nrow = 1) +
  #geom_line(data = predicted_100kb_bins_non_overlapped, aes(mid_posi_bin, predicted_rec_f), size = 1, color = "red") + 
  geom_line(data = predicted_100kb_bins_non_overlapped_quantile_gene_size, aes(mid_posi_bin, experimental_rec_f), size = 1, color = "blue") + 
  geom_line(data = gene_type_l_3kb_non_overlapped, aes(mid_posi_bin, Recrate, color = gene_type), size = 4) + theme_bw() + 
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), 
        axis.text = element_text(size = 18, face = "bold"), legend.position = c(0.85, 0.9), legend.box.background = element_rect(colour = "black"), plot.margin = margin(1,0.2,0,0.2, "cm")) + 
  scale_x_continuous(name = "", breaks=c(0, 30, 80, 130, 160), labels=c("-3 kb", "TSS", "gene bodies", "TTS", "3 kb")) +
  scale_y_continuous(breaks = seq(0, 8, 0.5)) 



#merge state and CO rate of 4 sizes of genes
gene_size_merged <- ggarrange(gene_size_top_state, gene_size_bottom_CO, nrow = 2, ncol = 1)
ggsave("./analysis_output/gene_size_merged.jpeg", gene_size_merged, width = 480, height = 360, units = c("mm"), dpi = 320)


######SP Fig S4-S6 (IRs are classified into three categories based on the transcription direction of genes on their both sides)#######
#categorize IRs into 4 quantiles (SV_status polled)
IR_size_4quantiles_list <- TAIR10_protein_coding_genes_IR_bed_trimmed %>%
  mutate(IR_name = ID) %>%
  left_join(TAIR10_protein_coding_genes_IR_bed_trimmed_no_SV_intact_list) %>%
  replace_na(list(SV_status = "with_SV")) %>%
  #group_by(SV_status, IR_type) %>%
  group_by(IR_type) %>%
  mutate(size = end-str) %>%
  mutate(quantile_size = if_else(size <= quantile(size, 0.25), "IR ≤ 25%", 
                                 if_else(size > quantile(size, 0.25) & size <= quantile(size, 0.5), "25% < IR ≤ 50%",
                                         if_else(size > quantile(size, 0.5) & size <= quantile(size, 0.75), "50% < IR ≤ 75%", "75% < IR ≤ 100%")))) %>%
  ungroup() %>%
  mutate(quantile_size = factor(quantile_size, levels = c("IR ≤ 25%", "25% < IR ≤ 50%",  "50% < IR ≤ 75%", "75% < IR ≤ 100%"))) %>%
  #group_by(IR_type, quantile_size) %>%
  #summarise(median_size = median(size), count = n(), min_size = min(size), max_size = max(size))
  select(IR_name, SV_status, IR_type, quantile_size)


#investigate the size of IR in 4 quantiles in terms of different transcription categories
IR_info <- TAIR10_protein_coding_genes_IR_bed_trimmed %>%
  mutate(IR_name = ID) %>%
  left_join(TAIR10_protein_coding_genes_IR_bed_trimmed_no_SV_intact_list) %>%
  replace_na(list(SV_status = "with_SV")) %>%
  group_by(SV_status, IR_type) %>%
  mutate(size = end-str) %>%
  mutate(quantile_size = if_else(size <= quantile(size, 0.25), "IR ≤ 25%", 
                                 if_else(size > quantile(size, 0.25) & size <= quantile(size, 0.5), "25% < IR ≤ 50%",
                                         if_else(size > quantile(size, 0.5) & size <= quantile(size, 0.75), "50% < IR ≤ 75%", "75% < IR ≤ 100%")))) %>%
  mutate(qunatile_size = factor(quantile_size, levels = c("IR ≤ 25%", "25% < IR ≤ 50%",  "50% < IR ≤ 75%", "75% < IR ≤ 100%"))) %>%
  group_by(SV_status, IR_type, quantile_size) %>%
  summarise(median_size = median(size), count = n(), min_size = min(size), max_size = max(size)) 

write_csv(IR_info, "./analysis_output/IR_size_info.csv", col_names = TRUE)

#show the info of mininum/maximum/median size, counts of each quantile of different IR categories
TAIR10_protein_coding_genes_IR_bed_trimmed_no_SV_intact %>%
  group_by(IR_type) %>%
  mutate(size = end-str) %>%
  mutate(quantile_size = if_else(size <= quantile(size, 0.25), "≤ 25%", 
                                 if_else(size > quantile(size, 0.25) & size <= quantile(size, 0.5), "25% < IR ≤ 50%",
                                         if_else(size > quantile(size, 0.5) & size <= quantile(size, 0.75), "50% < IR ≤ 75%", "75% < IR ≤ 100%")))) %>%
  mutate(qunatile_size = factor(quantile_size, levels = c("≤ 25%", "25% < IR ≤ 50%",  "50% < IR ≤ 75%", "75% < IR ≤ 100%"))) %>%
  group_by(IR_type, quantile_size) %>%
  summarise(median_size = median(size), count = n(), min_size = min(size), max_size = max(size)) 


#create the information of 10 states of each bin in 4 categories of IR sizes in terms of 3 IR types (all IR, with or without SVs)
TAIR10_protein_coding_genes_IR_100bins_integral_4quantiles_state_frac <- TAIR10_protein_coding_genes_IR_100bins_integral_state %>%
  left_join(TAIR10_protein_coding_genes_IR_100bins_integral_state_sum) %>%
  replace_na(list(sum_state_bp = 0)) %>%
  left_join(IR_size_4quantiles_list) %>%
  #filter(is.na(quantile_size) == TRUE)
  #group_by(SV_status, IR_type, quantile_size, bin, state) %>%
  group_by(IR_type, quantile_size, bin, state) %>%
  summarise(sum_state_bp = sum(sum_state_bp)) %>%
  #group_by(SV_status, IR_type, quantile_size, bin) %>%
  group_by(IR_type, quantile_size, bin) %>%
  mutate(sum_bin_bp = sum(sum_state_bp)) %>%
  mutate(state_fraction = sum_state_bp/sum_bin_bp) %>%
  mutate(mid_posi_bin = (bin + bin - 1)/2)

#transform the previous table into a list based on 3 IR types
IR_4q_state_frac_l <- TAIR10_protein_coding_genes_IR_100bins_integral_4quantiles_state_frac %>%
  #split(.$SV_status) %>%
  split(.$IR_type)



#create the information of CO rate of each bin in 4 categories of IR sizes in terms of 3 IR types
TAIR10_protein_coding_genes_IR_100bins_integral_4quantiles_state_frac_rec <- TAIR10_protein_coding_genes_IR_100bins_integral_4quantiles_state_frac %>%
  arrange(mid_posi_bin) %>%
  left_join(experimental_rec_state) %>%
  mutate(experimental_rec_part = experimental_rec * state_fraction) %>%
  #group_by(SV_status, IR_type, quantile_size, sum_bin_bp, mid_posi_bin) %>%
  group_by(IR_type, quantile_size, sum_bin_bp, mid_posi_bin) %>%
  #mutate(experimental_rec_f = sum(experimental_rec_part)) %>%
  summarise(experimental_rec_bin = sum(experimental_rec_part)) %>%
  #group_by(SV_status, IR_type, quantile_size) %>%
  group_by(IR_type, quantile_size) %>%
  mutate(experimental_rec_mean = sum(experimental_rec_bin*sum_bin_bp)/sum(sum_bin_bp))

#transform the previous table into a list based on 3 IR types
IR_4q_CO_frac_l <- TAIR10_protein_coding_genes_IR_100bins_integral %>%
  left_join(TAIR10_protein_coding_genes_IR_100bins_integral_RCO_sum) %>%
  replace_na(list(sum_CO = 0)) %>%
  left_join(IR_size_4quantiles_list) %>%
  #group_by(SV_status, IR_type, quantile_size, bin) %>%
  group_by(IR_type, quantile_size, bin) %>%
  summarise(sum_CO = sum(sum_CO), sum_bp = sum(end-str)) %>%
  mutate(Recrate = sum_CO/sum_bp/2182/2*10^8, mid_posi_bin = (bin + bin - 1)/2) %>%
  #group_by(SV_status, IR_type, quantile_size) %>%
  group_by(IR_type, quantile_size) %>%
  mutate(sum_CO_total = sum(sum_CO), sum_bp_total = sum(sum_bp)) %>%
  mutate(Recrate_mean = sum_CO_total/sum_bp_total/2182/2*10^8) %>%
  left_join(TAIR10_protein_coding_genes_IR_100bins_integral_4quantiles_state_frac_rec) %>%
  #split(.$SV_status) %>%
  split(.$IR_type)


#create functions for making figures
#the function for plotting the information of state distribution
IR_figure_state_size_f <- function(data){  
  data %>%
    #filter(IR_type == "parallel") %>%
    ggplot() + geom_col(aes(x = mid_posi_bin, y = state_fraction, fill = state), width = 1) +
    theme_bw() +
    scale_fill_manual(values=c(pal_npg("nrc", alpha = 0.7)(10))) + 
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 14, face = "bold"),
          legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.position = "bottom", 
          axis.text = element_text(size = 18, face = "bold"), plot.margin = margin(0,0.2,1,0.2, "cm")) + scale_x_continuous(name = c("")) +
    scale_y_continuous(name = c("The genome fraction of chromatin states"), breaks = c(0, 0.5, 1)) + guides(fill = guide_legend(nrow = 1)) +
    facet_wrap(~quantile_size, nrow = 1)
}

#the function for plotting CO rate
IR_figure_CO_size_f <- function(data) {  
  data %>%
    ggplot() +
    geom_col(mapping = aes(mid_posi_bin, Recrate), width = 1)+
    geom_line(aes(mid_posi_bin, experimental_rec_bin), size = 1, color = "blue") + 
    labs(x = "bin of IR", y = "Recombination rate (cM/Mb)") + 
    theme_bw() +
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 14, face = "bold"),
          legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), 
          axis.text = element_text(size = 18, face = "bold"), legend.position = c(0.85, 0.9), legend.box.background = element_rect(colour = "black"), plot.margin = margin(1,0.2,0,0.2, "cm")) + 
    scale_x_continuous(name = "") +
    scale_y_continuous(breaks = seq(0, 8, 1.5)) +
    facet_wrap(~quantile_size, nrow = 1)
}

#create the plots for state distribution of 3 IR types
IR_figure_state_size_l_all <- IR_4q_state_frac_l %>%
  map(. %>% IR_figure_state_size_f)

#create plots for CO rate of 3 IR types
IR_figure_CO_size_l_all <- vector("list", length = 3)

for (i in seq_along(IR_figure_CO_size_l_all)) {
  IR_figure_CO_size_l_all[[i]] <-  IR_4q_CO_frac_l[[i]] %>%
    IR_figure_CO_size_f() +
    geom_line(aes(mid_posi_bin, Recrate_mean), linetype = "dashed") +
    geom_line(aes(mid_posi_bin, experimental_rec_mean), linetype = "dashed", color = "blue")
}

#combine plots of 10-state distribution and CO rate
IR_figure_state_CO_l_all <- list(IR_figure_state_size_l_all, IR_figure_CO_size_l_all) %>%
  pmap(ggarrange, nrow = 2)

paths_IR_CO_rate_all <- str_c("./analysis_output/IR_all_state_CO_size_", c("convergent", "divergent", "parallel"),".jpeg")

#export SP Fig S4-S6
pwalk(list(paths_IR_CO_rate_all, IR_figure_state_CO_l_all), ggsave, width = 480, height = 360, units = c("mm"), dpi = 320)


#Fig S2
#use the experimental average recombination rate of 10 states to predict CO rates using 100-kb bins
#create the initial bed files using different bins for producing bed files intersecting features
Ara_Chr_label <- vector(mode = "list", length = 5)
bin_size <- c(500, 1000, 2000, 3000, 4000, 5000, 10000, 15000, 20000, 50000, 100000, 200000, 500000, 1000000)

Chr_label = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

Ara_genome_bed <-
  read_delim(
    "data/Fig2/Ara_genome_bed",
    delim = "\t",
    col_names = c("Chr", "str", "end")
  )

df_new = vector(mode = "list", length = 5)

for (chr in seq_along(Ara_Chr_label)) {
  for (size in seq_along(bin_size)) {
    Ara_Chr_label[[chr]][[size]] <-
      c(
        seq(Ara_genome_bed$str[[chr]], Ara_genome_bed$end[[chr]], by = bin_size[[size]]),
        Ara_genome_bed$end[[chr]]
      )
    
    df_new[[chr]][[size]] <- tibble(
      Chr = rep(Chr_label[[chr]], length(Ara_Chr_label[[chr]][[size]])),
      str = c(Ara_Chr_label[[chr]][[size]]),
      size_l = bin_size[[size]]
    ) %>%
      mutate(end = lead(str)) %>%
      drop_na() %>%
      mutate(bin = seq(1:n())) %>%
      mutate(chr_bin = str_c(Chr, "_", bin)) %>%
      select(Chr, str, end, chr_bin, size_l)
  }
}

den_table_list <-
  bind_rows(df_new[[1]], df_new[[2]], df_new[[3]], df_new[[4]], df_new[[5]]) %>%
  split(.$size_l)

write_delim(den_table_list[[11]], "./data/den_table_100k_bed", delim = "\t", col_names = FALSE)

#create the list of states for being joined by the information of sum of state 
den_table_list_state <- bind_rows(den_table_list) %>%
  mutate(state1 = "state1", state2 = "state2", state3 = "state3", state4 = "state4", state5 = "state5", state6 = "state6", state7 = "state7", state8 = "state8", state9 = "state9", SV = "SV") %>%
  gather(key = "state", value = "state_repeat", c(str_c("state", c(1:9)), "SV")) %>%
  select(-state_repeat) %>%
  split(.$size_l)


#make the intersection between 10 states and the 100-bin bed file, and import 10-state information and calculate the sum of base pairs of 10 states
#Open the terminal, run the command below in the directory "data/Fig2/"
#bedtools intersect -a state_10_raw -b ../den_table_100k_bed -wb | bedtools sort | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8}' > den_table_100k_state_10_intersect
#bedtools intersect -a ../den_table_100k_bed -b R_CO_final_bed -wb | awk -v OFS="\t" '{print$1, $2, $3, $4, $6, $7, $8, $9, $10}' > den_table_100k_RCO_raw_bed

#calculate the sum of base pairs of each state in 100-kb bins
den_table_state_sum_100k <- read_delim(
  "./data/Fig2/den_table_100k_state_10_intersect",
  delim = "\t",
  col_names = c("Chr", "str", "end", "state", "str_bin", "end_bin", "chr_bin")
) %>%
  group_by(chr_bin, state) %>%
  summarise(sum_state_bp = sum(end-str), .groups = "drop")

#calculate the sum of COs in 100-kb bins
den_table_RCO_sum_100k <- read_delim(
  "./data/Fig2/den_table_100k_RCO_raw_bed",
  delim = "\t",
  col_names = c(
    "Chr",
    "str",
    "end",
    "chr_bin",
    "Chr_CO",
    "str_CO",
    "end_CO",
    "sel_420",
    "CO_l"
  )
) %>%
  mutate(CO_n = (end - str) / (end_CO - str_CO)) %>%
  group_by(chr_bin) %>%
  summarise(sum_CO = sum(CO_n))


#create the fraction of 10 states in each bin, then use 10 states from "experimental_rec_state" to produce predicted CO rate
#and calculate experimental CO rate of each bin
den_table_state_f_100k_plot <- den_table_list_state[[11]] %>%
  left_join(den_table_state_sum_100k) %>%
  replace_na(list(sum_state_bp = 0)) %>%
  mutate(den_state = sum_state_bp/(end-str)) %>%
  select(-sum_state_bp) %>%
  spread(key = state, value = den_state) %>%
  gather(key = "state", value = "state_frac", c(6:15)) %>%
  left_join(experimental_rec_state) %>%
  mutate(rec_part = state_frac*experimental_rec) %>%
  group_by(Chr, str, end, chr_bin) %>%
  summarise(rec_based_on_state = sum(rec_part)) %>%
  left_join(den_table_RCO_sum_100k) %>%
  replace_na(list(sum_CO = 0)) %>%
  mutate(Recrate_Rowan = (sum_CO) / 2 / 2182 / (end - str) * 10 ^ 8) 


#calculate explained variance of this prediction
explained_var_test_plot <- round(1-sum((den_table_state_f_100k_plot$Recrate_Rowan - den_table_state_f_100k_plot$rec_based_on_state)^2)/sum((den_table_state_f_100k_plot$Recrate_Rowan - mean(den_table_state_f_100k_plot$Recrate_Rowan))^2),2)

#produce Fig S2
#define pericentromeric regions based on methylation levels (based on Underwood's research, https://genome.cshlp.org/content/early/2018/03/09/gr.227116.117)
Ara_peri_posi <- tibble(
  Chr = str_c("Chr", 1:5),
  str_peri = c(11420001, 910001, 10390001, 1070001, 8890001),
  end_peri = c(18270000, 7320000, 16730000, 6630000, 15550000)
)

SP_Figure_S2 <- den_table_state_f_100k_plot %>%
  left_join(Ara_peri_posi, by = "Chr") %>%
  mutate(status = if_else(
    end < str_peri |
      str > end_peri,
    "arms",
    if_else(
      str > str_peri &
        end < end_peri,
      "pericentromeric_regions",
      "overlapped"
    )
  )) %>% ggscatter("rec_based_on_state", "Recrate_Rowan", color = "status", palette = c(pal_npg("nrc", alpha = 1)(6))[c(1,3,4)]) +
  geom_abline(slope = 1, size = 1) +
  theme_bw() +
  labs(y = "Recombination rate (cM/Mb)", x = "Predicted recombination rate (cM/Mb)") +
  theme(
    #axis.ticks = element_blank(),
    strip.text.x = element_text(
      colour = "black",
      face = "bold",
      size = 20
    ),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    legend.position="bottom"
  ) + scale_y_continuous(limits =  c(0,12), breaks=seq(0, 12,2)) + geom_text(x =  1.5, y = 10, label = str_c("R^2 == ", 0.24), parse = TRUE, size = 6)


ggsave("analysis_output/SP_Figure_S2.jpeg", SP_Figure_S2, width = 300, height = 225, units = c("mm"))
