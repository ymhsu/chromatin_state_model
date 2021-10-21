#import necessary packages
Packages <- c("scales", "tidyverse", "ggrepel", "ggsci", "ggpubr", "doMC", "doParallel", "foreach", "slider", "cowplot", "ggpmisc", "ggformula")
lapply(Packages, library, character.only = TRUE)

##Figure 1##
#change the current directory as the working directory
setwd("/data/projects/thesis/INRA_project/Ara_TE_task/R_markdown/Model_1st/chromatin_state_model/")

#create the initial bed files using different bins for producing bed files intersecting features
Ara_Chr_label <- vector(mode = "list", length = 5)
bin_size <- c(500, 1000, 2000, 3000, 4000, 5000, 10000, 15000, 20000, 50000, 100000, 200000, 500000, 1000000)

Chr_label = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

Ara_genome_bed <-
  read_delim(
    "data/Fig1/Ara_genome_bed",
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


bin_name <-
  c("0_5k", "1k", "2k", "3k", "4k", "5k", "10k", "15k", "20k", "50k", "100k", "200k", "500k", "1000k")

paths_den_table <-
  str_c(
    "data/Fig1/",
    "den_table_",
    bin_name,
    "_bed"
  )

pwalk(list(
  den_table_list,
  paths_den_table,
  delim = "\t",
  col_names = FALSE
),
write_delim)

#there are 9 features including, gene, TE, TSS, H3K4me1, H3K4me3, H3K9me2, H3K27me3, ATAC, DNase, their information are shown below,
 
#TSS
#cat TAIR10_protein_coding_genes_bed_sorted | awk -v OFS="\t" '{if($5 == "+") print $1,$2,$2+1,"TSS",$5,$6; else print $1,$3-1,$3,"TSS",$5,$6}' > TAIR10_protein_coding_genes_TSS_bed_sorted

#produce the protein-coding genes without the names of genes (4th column was deleted from gene: TAIR10_protein_coding_genes_bed_sorted_merge, the production of this file refers to Fig2_new.R)
#cat TAIR10_protein_coding_genes_bed_sorted_merge | awk '{print$1"\t"$2"\t"$3}' > TAIR10_protein_coding_genes_bed_sorted_merge_noname

#TE (TAIR10_GFF3_genes_transposons.gff was downloaded from TAIR10 website, https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release%2FTAIR10_gff3)
#cat -n TAIR10_GFF3_genes_transposons.gff | awk '{if ($1 >= 590265) print$0}' | grep "transposable_element" | awk '{print$2"\t"$5-1"\t"$6"\t""transposable_element"}' | bedtools sort > ./TAIR10_transposable_elements_sorted
#bedtools merge -i TAIR10_transposable_elements_sorted | awk '{print$1"\t"$2"\t"$3"\t""TE"}' > TAIR10_transposable_elements_sorted_merge

#H3K4me1
#GSM3674621_H3K4me1_leaf_R1.bedgraph.gz

#H3K4me3
#GSM3674620_H3K4me3_leaf_R1.bedgraph.gz

#H3K9me2
#GSM4734580_H3K9me2_leaf_R1.bedgraph.gz

#H3K27me3
#GSM3674617_H3K27me3_leaf_R1.bedgraph.gz

#ATAC
#GSM3674715_ATAC_leaf_R1.bedgraph.gz

#DNase
#GSM1289358_DNase_7d_seedling_R1.bedgraph.gz

#using "Fig1_9features_intersect_14bins_f.sh" to get the bed (SVs intersecting bins or features were not removed)

#read CO file and calculate recombination rate
paths_den_table_RCO <-
  str_c(
    "data/Fig1/",
    "den_table_",
    bin_name,
    "_RCO_raw_bed"
  )

den_table_CO_sum = vector("list", length = length(bin_size))

for (size in seq_along(bin_size)) {
  den_table_CO_sum[[size]] <-
    read_delim(
      paths_den_table_RCO[[size]],
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
}


#define pericentromeric regions based on methylation levels (based on Underwood's research, https://genome.cshlp.org/content/early/2018/03/09/gr.227116.117)
Ara_peri_posi <- tibble(
  Chr = str_c("Chr", 1:5),
  str_peri = c(11420001, 910001, 10390001, 1070001, 8890001),
  end_peri = c(18270000, 7320000, 16730000, 6630000, 15550000)
)



##calculate sum of features
#set paths for importing 9 features
paths_den_table_gene_SV <- str_c("data/Fig1/", "den_table_", bin_name, "_gene_SV_intersect")
paths_den_table_TE_SV <- str_c("data/Fig1/", "den_table_", bin_name, "_TE_SV_intersect")
paths_den_table_H3K27me3_SV <- str_c("data/Fig1/", "den_table_", bin_name,"_H3K27me3_SV_intersect")
paths_den_table_H3K4me3_SV <- str_c("data/Fig1/", "den_table_", bin_name, "_H3K4me3_SV_intersect")
paths_den_table_H3K4me1_SV <- str_c("data/Fig1/", "den_table_", bin_name, "_H3K4me1_SV_intersect")
paths_den_table_H3K9me2_SV <- str_c("data/Fig1/", "den_table_", bin_name, "_H3K9me2_SV_intersect")
paths_den_table_DNase_SV <- str_c("data/Fig1/", "den_table_", bin_name, "_DNase_SV_intersect")
paths_den_table_ATAC_SV <- str_c("data/Fig1/", "den_table_", bin_name, "_ATAC_SV_intersect")
paths_den_table_TSS_SV <- str_c("data/Fig1/", "den_table_", bin_name, "_TSS_SV_raw_bed")

#create lists for 9 features
den_table_gene_sum_SV = vector("list", length = length(bin_size))
den_table_TE_sum_SV = vector("list", length = length(bin_size))
den_table_H3K27me3_sum_SV = vector("list", length = length(bin_size))
den_table_H3K4me3_sum_SV = vector("list", length = length(bin_size))
den_table_H3K4me1_sum_SV = vector("list", length = length(bin_size))
den_table_H3K9me2_sum_SV = vector("list", length = length(bin_size))
den_table_TSS_sum_SV = vector("list", length = length(bin_name))
den_table_DNase_sum_SV = vector("list", length = length(bin_name))
den_table_ATAC_sum_SV = vector("list", length = length(bin_name))



#calculation of sum of signals of 9 features in each bin
for (size in seq_along(bin_size)) {
  den_table_gene_sum_SV[[size]] <-
    read_delim(
      paths_den_table_gene_SV[[size]],
      delim = "\t",
      col_names = c("Chr", "str", "end", "chr_bin")
    ) %>%
    mutate(size = end - str) %>%
    group_by(chr_bin) %>%
    summarise(gene_sum_size = sum(size))
  
  den_table_TE_sum_SV[[size]] <-
    read_delim(
      paths_den_table_TE_SV[[size]],
      delim = "\t",
      col_names = c("Chr", "str", "end", "S_family", "chr_bin")
    ) %>%
    mutate(size = end - str) %>%
    group_by(chr_bin) %>%
    summarise(TE_sum_size = sum(size))
  
  den_table_H3K27me3_sum_SV[[size]] <-
    read_delim(
      paths_den_table_H3K27me3_SV[[size]],
      delim = "\t",
      col_names = c("Chr", "str", "end", "intensity", "chr_bin")
    ) %>%
    group_by(chr_bin) %>%
    summarise(sum_H3K27me3 = sum((end - str) * intensity))
  
  
  den_table_H3K4me1_sum_SV[[size]] <-
    read_delim(
      paths_den_table_H3K4me1_SV[[size]],
      delim = "\t",
      col_names = c("Chr", "str", "end", "intensity", "chr_bin")
    ) %>%
    group_by(chr_bin) %>%
    summarise(sum_H3K4me1 = sum((end - str) * intensity))
  
  den_table_H3K4me3_sum_SV[[size]] <-
    read_delim(
      paths_den_table_H3K4me3_SV[[size]],
      delim = "\t",
      col_names = c("Chr", "str", "end", "intensity", "chr_bin")
    ) %>%
    group_by(chr_bin) %>%
    summarise(sum_H3K4me3 = sum((end - str) * intensity))
  
  
  den_table_H3K9me2_sum_SV[[size]] <-
    read_delim(
      paths_den_table_H3K9me2_SV[[size]],
      delim = "\t",
      col_names = c("Chr", "str", "end", "intensity", "chr_bin")
    ) %>%
    group_by(chr_bin) %>%
    summarise(sum_H3K9me2 = sum((end - str) * intensity))
  
  
  den_table_TSS_sum_SV[[size]] <- read_delim(paths_den_table_TSS_SV[[size]], delim = "\t", col_names = c("Chr", "str", "end", "feature", "strand", "gene_name", "Chr_bin", "str_bin", "end_bin", "chr_bin", "size_l")) %>%
    group_by(chr_bin) %>%
    summarise(sum_TSS = sum(end-str))
  
  den_table_DNase_sum_SV[[size]] <-
    read_delim(
      paths_den_table_DNase_SV[[size]],
      delim = "\t",
      col_names = c("Chr", "str", "end", "intensity", "chr_bin")
    ) %>%
    group_by(chr_bin) %>%
    summarise(sum_DNase = sum((end - str) * intensity))
  
  den_table_ATAC_sum_SV[[size]] <-
    read_delim(
      paths_den_table_ATAC_SV[[size]],
      delim = "\t",
      col_names = c("Chr", "str", "end", "intensity", "chr_bin")
    ) %>%
    group_by(chr_bin) %>%
    summarise(sum_ATAC = sum((end - str) * intensity))
}

#join the sum of 9 features and COs to the raw table of bins for the following analysis
den_table_f_cor_SV = vector("list", length = length(bin_size))

for (size in seq_along(bin_size)) {
  den_table_f_cor_SV[[size]] <- den_table_list[[size]] %>%
    left_join(den_table_H3K27me3_sum_SV[[size]], by = "chr_bin") %>%
    left_join(den_table_H3K4me1_sum_SV[[size]], by = "chr_bin") %>%
    left_join(den_table_H3K4me3_sum_SV[[size]], by = "chr_bin") %>%
    left_join(den_table_H3K9me2_sum_SV[[size]], by = "chr_bin") %>%
    left_join(den_table_gene_sum_SV[[size]], by = "chr_bin") %>%
    left_join(den_table_TE_sum_SV[[size]], by = "chr_bin") %>%
    left_join(den_table_TSS_sum_SV[[size]], by = "chr_bin") %>%
    left_join(den_table_ATAC_sum_SV[[size]], by = "chr_bin") %>%
    left_join(den_table_DNase_sum_SV[[size]], by = "chr_bin") %>%
    left_join(den_table_CO_sum[[size]], by = c("chr_bin")) %>%
    replace_na(
      list(
        TE_sum_size = 0,
        sum_H3K27me3 = 0,
        sum_H3K4me3 = 0,
        sum_H3K4me1 = 0,
        sum_H3K9me2 = 0,
        gene_sum_size = 0,
        sum_TSS = 0,
        sum_ATAC = 0,
        sum_DNase = 0,
        sum_CO = 0
      )
    ) %>%
    mutate(
      den_gene = gene_sum_size / (end - str),
      den_TE = TE_sum_size / (end - str),
      den_H3K27me3 = sum_H3K27me3 / (end - str),
      den_H3K4me3 = sum_H3K4me3 / (end - str),
      den_H3K4me1 = sum_H3K4me1 / (end - str),
      den_H3K9me2 = sum_H3K9me2 / (end - str),
      den_TSS = sum_TSS*10000/(end-str),
      den_ATAC = sum_ATAC/(end-str),
      den_DNase = sum_DNase/(end-str)
    ) %>%
    mutate(RecRate = (sum_CO) / 2 / 2182 / (end - str) * 10 ^ 8) %>%
    mutate(
      den_DNase = if_else(
        den_DNase == "NaN" |
          den_DNase == 0 |
          den_DNase == "Inf",
        10 ^ -10,
        den_DNase
      )
    ) %>%
    mutate(
      den_H3K27me3 = if_else(
        den_H3K27me3 == "NaN" |
          den_H3K27me3 == 0 |
          den_H3K27me3 == "Inf",
        10 ^ -10,
        den_H3K27me3
      )
    ) %>%
    mutate(
      den_H3K4me3 = if_else(
        den_H3K4me3 == "NaN" |
          den_H3K4me3 == 0 | den_H3K4me3 == "Inf",
        10 ^ -10,
        den_H3K4me3
      )
    ) %>%
    mutate(
      den_H3K4me1 = if_else(
        den_H3K4me1 == "NaN" |
          den_H3K4me1 == 0 | den_H3K4me1 == "Inf",
        10 ^ -10,
        den_H3K4me1
      )
    ) %>% 
    mutate(
      den_H3K9me2 = if_else(
        den_H3K9me2 == "NaN" |
          den_H3K9me2 == 0 | den_H3K9me2 == "Inf",
        10 ^ -10,
        den_H3K9me2
      )
    ) %>%
    select(
      -sum_H3K27me3,
      -sum_H3K4me3,
      -sum_H3K4me1,
      -sum_H3K9me2)
}


#extract 100 kb for producing figures, and divided them into separate table for ploting figure later
data_100kb_cor_feature_SV <- den_table_f_cor_SV[[11]] %>%
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
  )) %>%
  gather(key = "all_feature", value = "den_feature", c("den_gene", "den_TE", "den_H3K4me1", "den_H3K4me3", "den_H3K9me2", "den_H3K27me3", "den_TSS", "den_ATAC", "den_DNase")) %>%
  mutate(all_feature = str_remove(all_feature, "den_")) %>%
  mutate(all_feature = if_else(all_feature == "H3K4me1" | all_feature == "H3K4me3" | all_feature == "H3K9me2" | all_feature == "H3K27me3" | all_feature == "ATAC", str_c(all_feature, "_leaves"), 
                               if_else(all_feature == "DNase", str_c(all_feature, "_seedlings"), all_feature))) %>%
  mutate(all_feature = factor(all_feature, levels = c("gene", "TE", "TSS", "H3K4me1_leaves", "H3K4me3_leaves", "H3K9me2_leaves", "H3K27me3_leaves", "ATAC_leaves", "DNase_seedlings"))) %>%
  split(.$all_feature)


#produce the function for calculating explained variance
explained_var_predicted_single_feature_smaller_spar <- function(data){
  prediction <- smooth.spline(data$den_feature, data$RecRate, spar = 1)
  predicted_y <- predict(prediction, data$den_feature)$y
  1-sum((data$RecRate-predicted_y)^2)/sum((data$RecRate-mean(data$RecRate))^2)
}

#produce the function for ploting figures
Fig_1_100kb_cor_feature_func_smaller_spar <- function(data){
  data %>% ggscatter("den_feature", "RecRate", color = "status", palette = c(pal_npg("nrc", alpha = 1)(6))[c(1,3,4)]) +
    geom_spline(
      aes(den_feature, RecRate),
      color = "black",
      size  = 1,
      spar = 1
    ) +
    facet_wrap( ~ all_feature, nrow = 1, scales = "free") +
    theme_bw() +
    labs(y = "Recombination rate(cM/Mb)", x = "The density of features") +
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
    )
}


#produce explained variance of 9 features in a vector
explained_var_vector_SV <- round(unlist(data_100kb_cor_feature_SV %>%
                                          map(. %>% explained_var_predicted_single_feature_smaller_spar())), 2)



#produce the list of 9 figures from 9 features
isolate_feature_figure_SV <- vector("list", length(data_100kb_cor_feature_SV))

for (i in seq_along(isolate_feature_figure_SV)) {
  x_location = 0.23*(quantile(data_100kb_cor_feature_SV[[i]]$den_feature, 0.975) - quantile(data_100kb_cor_feature_SV[[i]]$den_feature, 0.025)) + quantile(data_100kb_cor_feature_SV[[i]]$den_feature, 0.025)
  
  x_limits = c(quantile(data_100kb_cor_feature_SV[[i]]$den_feature, 0.025), quantile(data_100kb_cor_feature_SV[[i]]$den_feature, 0.975))
  
  isolate_feature_figure_SV[[i]] <- data_100kb_cor_feature_SV[[i]] %>%
    Fig_1_100kb_cor_feature_func_smaller_spar() + 
    theme(axis.title.y = element_text(color = "white"), axis.title.x = element_text(color = "white"),  legend.position = "", plot.margin = margin(0,0.5,0,0.2, "cm")) +
    geom_text(x =  x_location, y = 15, label = str_c("R^2 == ", explained_var_vector_SV[[i]]), parse = TRUE, size = 6) +
    scale_x_continuous(limits = c(x_limits[[1]], x_limits[[2]]))
}


#add the label of y-axis for the combined figure
isolate_feature_figure_SV[[4]] <- isolate_feature_figure_SV[[4]] +
  theme(axis.title.y = element_text(color = "black"))

#add the label of x-axis for the combined figure
isolate_feature_figure_SV[[8]] <- isolate_feature_figure_SV[[8]] +
  theme(axis.title.x = element_text(color = "black"))

Fig1_test_SV <- ggarrange(isolate_feature_figure_SV[[1]], isolate_feature_figure_SV[[2]], isolate_feature_figure_SV[[3]], 
                          isolate_feature_figure_SV[[4]], isolate_feature_figure_SV[[5]], isolate_feature_figure_SV[[6]],
                          isolate_feature_figure_SV[[7]], isolate_feature_figure_SV[[8]], isolate_feature_figure_SV[[9]], nrow = 3, ncol = 3, align = "v")

ggsave("analysis/Fig1.jpeg", Fig1_test_SV, width = 400, height = 330, units = c("mm"))


##Table S2##
#perform the additive model and the model with interactions terms using data with SV
summary_linear <- vector("list", length = length(bin_size))
summary_interaction <- vector("list", length = length(bin_size))

for (i in seq_along(bin_size)) {
  lm_linear <- lm(RecRate ~ den_gene + den_TE + den_TSS + den_H3K4me1 + den_H3K4me3 + den_H3K9me2 + den_H3K27me3 + den_ATAC + den_DNase, data = den_table_f_cor_SV[[i]])
  lm_interaction <- lm(RecRate ~ (den_gene + den_TE + den_H3K4me1 + den_H3K4me3 + den_H3K9me2 + den_H3K27me3 + den_TSS + den_DNase + den_ATAC)^2, data = den_table_f_cor_SV[[i]])
  summary_linear[[i]] <- summary(lm_linear)
  summary_interaction[[i]] <- summary(lm_interaction)
}

Table_S2 <- tibble(
  para = c("Intercept", "gene", "TE", "TSS", "H3K4me1", "H3K4me3", "H3K9me2", "H3K27me3", "ATAC", "DNase", "R_square"),
  add_50_k = round(c(summary_linear[[10]]$coefficients[,1], summary_linear[[10]]$r.squared), 2),
  add_100_k = round(c(summary_linear[[11]]$coefficients[,1], summary_linear[[11]]$r.squared), 2),
  add_200_k = round(c(summary_linear[[12]]$coefficients[,1], summary_linear[[12]]$r.squared), 2),
  add_500_k = round(c(summary_linear[[13]]$coefficients[,1], summary_linear[[13]]$r.squared), 2)
)

write_delim(Table_S2, "./analysis/Table_S2", delim = "\t", col_names = TRUE)


#To investigate the predictive power of the additive one and the one with interaction term
#We do the fit in one chromosome, and predict CO rate in 4 other chromosomes
den_table_f_cor_SV_Chr <- den_table_f_cor_SV %>%
  map(. %>% split(.$Chr))

#perform fitting in one chromosome
summary_linear_Chr_100k <- vector("list", length = length(c(1:5)))
summary_interaction_Chr_100k <- vector("list", length = length(c(1:5)))

for (i in seq_along(c(1:5))) {
  lm_linear <- lm(RecRate ~ den_gene + den_TE + den_TSS + den_H3K4me1 + den_H3K4me3 + den_H3K9me2 + den_H3K27me3 + den_ATAC + den_DNase, data = den_table_f_cor_SV_Chr[[11]][[i]])
  lm_interaction <- lm(RecRate ~ (den_gene + den_TE + den_H3K4me1 + den_H3K4me3 + den_H3K9me2 + den_H3K27me3 + den_TSS + den_DNase + den_ATAC)^2, data = den_table_f_cor_SV_Chr[[11]][[i]])
  summary_linear_Chr_100k[[i]] <- summary(lm_linear)
  summary_interaction_Chr_100k[[i]] <- summary(lm_interaction)
}


#create the function to get R square
R_square_lm_Chr <- function(data, data2){
  table_for_R <- data %>%
    mutate(intercept_para = data2$coefficients[,1][[1]]) %>%
    mutate(den_gene_para = data2$coefficients[,1][[2]]) %>%
    mutate(den_TE_para = data2$coefficients[,1][[3]]) %>%
    mutate(den_TSS_para = data2$coefficients[,1][[4]]) %>%
    mutate(den_H3K4me1_para = data2$coefficients[,1][[5]]) %>%
    mutate(den_H3K4me3_para = data2$coefficients[,1][[6]]) %>%
    mutate(den_H3K9me2_para = data2$coefficients[,1][[7]]) %>%
    mutate(den_H3K27me3_para = data2$coefficients[,1][[8]]) %>%
    mutate(den_ATAC_para = data2$coefficients[,1][[9]]) %>%
    mutate(den_DNase_para = data2$coefficients[,1][[10]]) %>%
    mutate(Recrate_pre = intercept_para + den_gene_para*den_gene + den_TE_para*den_TE + den_TSS_para*den_TSS + den_H3K4me1_para*den_H3K4me1 + den_H3K4me3_para*den_H3K4me3+den_H3K9me2_para*den_H3K9me2+den_ATAC_para*den_ATAC+den_DNase_para*den_DNase + den_H3K27me3_para*den_H3K27me3)
  
  1-sum((table_for_R$Recrate_pre-table_for_R$RecRate)^2)/sum((mean(table_for_R$RecRate)-table_for_R$RecRate)^2)
}

#calculating R square for the additive model
R_square_lm_result <- vector("list", length = 5)


for (i in seq_along(R_square_lm_result)) {
  for (j in seq_along(R_square_lm_result)) {
    R_square_lm_result[[i]][[j]] <- R_square_lm_Chr(den_table_f_cor_SV_Chr[[11]][[j]], summary_linear_Chr_100k[[i]])
  }
}

R_square_based_on_lm_Chr_table <- tibble(
  Chr = str_c("Chr", c(1:5), "_predict"),
  Chr1_fit = round(R_square_lm_result[[1]], 3),
  Chr2_fit = round(R_square_lm_result[[2]], 3),
  Chr3_fit = round(R_square_lm_result[[3]], 3),
  Chr4_fit = round(R_square_lm_result[[4]], 3),
  Chr5_fit = round(R_square_lm_result[[5]], 3)
)

write_delim(R_square_based_on_lm_Chr_table, "./analysis/R_square_based_on_lm_Chr_table", delim = "\t", col_names = TRUE)



#calculating R square for the model with interaction terms
#create the function to get R square
R_square_lm_interaction_Chr <- function(data, data2){
  table_for_R <- data %>%
    mutate(intercept_para = data2$coefficients[,1][[1]]) %>%
    mutate(den_gene_para = data2$coefficients[,1][[2]]) %>%
    mutate(den_TE_para = data2$coefficients[,1][[3]]) %>%
    mutate(den_H3K4me1_para = data2$coefficients[,1][[4]]) %>%
    mutate(den_H3K4me3_para = data2$coefficients[,1][[5]]) %>%
    mutate(den_H3K9me2_para = data2$coefficients[,1][[6]]) %>%
    mutate(den_H3K27me3_para = data2$coefficients[,1][[7]]) %>%
    mutate(den_TSS_para = data2$coefficients[,1][[8]]) %>%
    mutate(den_DNase_para = data2$coefficients[,1][[9]]) %>%
    mutate(den_ATAC_para = data2$coefficients[,1][[10]]) %>%
    mutate(den_gene_TE_para = data2$coefficients[,1][[11]]) %>%
    mutate(den_gene_H3K4me1_para = data2$coefficients[,1][[12]]) %>%
    mutate(den_gene_H3K4me3_para = data2$coefficients[,1][[13]]) %>%
    mutate(den_gene_H3K9me2_para = data2$coefficients[,1][[14]]) %>%
    mutate(den_gene_H3K27me3_para = data2$coefficients[,1][[15]]) %>%
    mutate(den_gene_TSS_para = data2$coefficients[,1][[16]]) %>%
    mutate(den_gene_DNase_para = data2$coefficients[,1][[17]]) %>%
    mutate(den_gene_ATAC_para = data2$coefficients[,1][[18]]) %>%
    mutate(den_TE_H3K4me1_para = data2$coefficients[,1][[19]]) %>%
    mutate(den_TE_H3K4me3_para = data2$coefficients[,1][[20]]) %>%
    mutate(den_TE_H3K9me2_para = data2$coefficients[,1][[21]]) %>%
    mutate(den_TE_H3K27me3_para = data2$coefficients[,1][[22]]) %>%
    mutate(den_TE_TSS_para = data2$coefficients[,1][[23]]) %>%
    mutate(den_TE_DNase_para = data2$coefficients[,1][[24]]) %>%
    mutate(den_TE_ATAC_para = data2$coefficients[,1][[25]]) %>%
    mutate(den_H3K4me1_H3K4me3_para = data2$coefficients[,1][[26]]) %>%
    mutate(den_H3K4me1_H3K9me2_para = data2$coefficients[,1][[27]]) %>%
    mutate(den_H3K4me1_H3K27me3_para = data2$coefficients[,1][[28]]) %>%
    mutate(den_H3K4me1_TSS_para = data2$coefficients[,1][[29]]) %>%
    mutate(den_H3K4me1_DNase_para = data2$coefficients[,1][[30]]) %>%
    mutate(den_H3K4me1_ATAC_para = data2$coefficients[,1][[31]]) %>%
    mutate(den_H3K4me3_H3K9me2_para = data2$coefficients[,1][[32]]) %>%
    mutate(den_H3K4me3_H3K27me3_para = data2$coefficients[,1][[33]]) %>%
    mutate(den_H3K4me3_TSS_para = data2$coefficients[,1][[34]]) %>%
    mutate(den_H3K4me3_DNase_para = data2$coefficients[,1][[35]]) %>%
    mutate(den_H3K4me3_ATAC_para = data2$coefficients[,1][[36]]) %>%
    mutate(den_H3K9me2_H3K27me3_para = data2$coefficients[,1][[37]]) %>%
    mutate(den_H3K9me2_TSS_para = data2$coefficients[,1][[38]]) %>%
    mutate(den_H3K9me2_DNase_para = data2$coefficients[,1][[39]]) %>%
    mutate(den_H3K9me2_ATAC_para = data2$coefficients[,1][[40]]) %>%
    mutate(den_H3K27me3_TSS_para = data2$coefficients[,1][[41]]) %>%
    mutate(den_H3K27me3_DNase_para = data2$coefficients[,1][[42]]) %>%
    mutate(den_H3K27me3_ATAC_para = data2$coefficients[,1][[43]]) %>%
    mutate(den_TSS_DNase_para = data2$coefficients[,1][[44]]) %>%
    mutate(den_TSS_ATAC_para = data2$coefficients[,1][[45]]) %>%
    mutate(den_DNase_ATAC_para = data2$coefficients[,1][[46]]) %>%
    mutate(Recrate_pre = intercept_para + den_gene_para*den_gene + den_TE_para*den_TE + den_TSS_para*den_TSS + 
             den_H3K4me1_para*den_H3K4me1 + den_H3K4me3_para*den_H3K4me3 + den_H3K9me2_para*den_H3K9me2 + 
             den_ATAC_para*den_ATAC + den_DNase_para*den_DNase + den_H3K27me3_para*den_H3K27me3+
             den_gene_TE_para*den_gene*den_TE+
             den_gene_H3K4me1_para*den_gene*den_H3K4me1+
             den_gene_H3K4me3_para*den_gene*den_H3K4me3+
             den_gene_H3K9me2_para*den_gene*den_H3K9me2+
             den_gene_H3K27me3_para*den_gene*den_H3K27me3+
             den_gene_TSS_para*den_gene*den_TSS+
             den_gene_DNase_para*den_gene*den_DNase+
             den_gene_ATAC_para*den_gene*den_ATAC+
             den_TE_H3K4me1_para*den_TE*den_H3K4me1+
             den_TE_H3K4me3_para*den_TE*den_H3K4me3+
             den_TE_H3K9me2_para*den_TE*den_H3K9me2+
             den_TE_H3K27me3_para*den_TE*den_H3K27me3+
             den_TE_TSS_para*den_TE*den_TSS+
             den_TE_DNase_para*den_TE*den_DNase+
             den_TE_ATAC_para*den_TE*den_ATAC+
             den_H3K4me1_H3K4me3_para*den_H3K4me1*den_H3K4me3+
             den_H3K4me1_H3K9me2_para*den_H3K4me1*den_H3K9me2+
             den_H3K4me1_H3K27me3_para*den_H3K4me1*den_H3K27me3+
             den_H3K4me1_TSS_para*den_H3K4me1*den_TSS+
             den_H3K4me1_DNase_para*den_H3K4me1*den_DNase+
             den_H3K4me1_ATAC_para*den_H3K4me1*den_ATAC+
             den_H3K4me3_H3K9me2_para*den_H3K4me3*den_H3K9me2+
             den_H3K4me3_H3K27me3_para*den_H3K4me3*den_H3K27me3+
             den_H3K4me3_TSS_para*den_H3K4me3*den_TSS+
             den_H3K4me3_DNase_para*den_H3K4me3*den_DNase+
             den_H3K4me3_ATAC_para*den_H3K4me3*den_ATAC+
             den_H3K9me2_H3K27me3_para*den_H3K9me2*den_H3K27me3+
             den_H3K9me2_TSS_para*den_H3K9me2*den_TSS+
             den_H3K9me2_DNase_para*den_H3K9me2*den_DNase+
             den_H3K9me2_ATAC_para*den_H3K9me2*den_ATAC+
             den_H3K27me3_TSS_para*den_H3K27me3*den_TSS+
             den_H3K27me3_DNase_para*den_H3K27me3*den_DNase+
             den_H3K27me3_ATAC_para*den_H3K27me3*den_ATAC+
             den_TSS_DNase_para*den_TSS*den_DNase+
             den_TSS_ATAC_para*den_TSS*den_ATAC+
             den_DNase_ATAC_para*den_DNase*den_ATAC)
  
  1-sum((table_for_R$Recrate_pre-table_for_R$RecRate)^2)/sum((mean(table_for_R$RecRate)-table_for_R$RecRate)^2)
}

R_square_lm_result_interaction <- vector("list", length = 5)


for (i in seq_along(R_square_lm_result_interaction)) {
  for (j in seq_along(R_square_lm_result_interaction)) {
    R_square_lm_result_interaction[[i]][[j]] <- R_square_lm_interaction_Chr(den_table_f_cor_SV_Chr[[11]][[j]], summary_interaction_Chr_100k[[i]])
  }
}

R_square_lm_interaction_Chr(den_table_f_cor_SV_Chr[[11]][[5]], summary_interaction_Chr_100k[[1]])

R_square_based_on_lm_interaction_Chr_table <- tibble(
  Chr = str_c("Chr", c(1:5), "_predict"),
  Chr1_fit = round(R_square_lm_result_interaction[[1]], 3),
  Chr2_fit = round(R_square_lm_result_interaction[[2]], 3),
  Chr3_fit = round(R_square_lm_result_interaction[[3]], 3),
  Chr4_fit = round(R_square_lm_result_interaction[[4]], 3),
  Chr5_fit = round(R_square_lm_result_interaction[[5]], 3)
)

write_delim(R_square_based_on_lm_interaction_Chr_table, "./analysis/R_square_based_on_lm_interaction_Chr_table", delim = "\t", col_names = TRUE)


##Figure S1##
#In order to indentify the signal of epigenetic marks in different tissue
#we downloaded features from different tissues from ncbi or EMBL database
#import the list of epigenetic marks
epi_mark_intersect_file_list <- read_delim("./data/Fig1/epigenomic_mark_intersect_file_list", delim = "\t", col_names = "name") %>%
  mutate(name_v2 = str_replace(name, "den_table_100k_", "")) %>%
  mutate(name_v2 = str_replace(name_v2, "_SV_intersect", "")) %>%
  mutate(name_v3 = sub("[A-Z0-9a-z-]*_[HDA][A-Z0-9a-z]*_(.*)_R[0-9]", "\\1", name_v2)) %>%
  mutate(label = sub("([A-Z0-9a-z-]*_[HDA][A-Z0-9a-z]*)_.*_R[0-9]", "\\1", name_v2)) %>%
  mutate(name_v3 = str_replace(name_v3, "_", " ")) %>%
  mutate(name_v3 = str_replace(name_v3, "_", " ")) %>%
  mutate(name_v3 = str_replace(name_v3, "_", " ")) %>%
  mutate(name_v3 = str_replace(name_v3, "roots non hair cells", "roots (non hair cells)")) %>%
  mutate(name_v3 = str_replace(name_v3, "4week", "4-week")) %>%
  mutate(name_v3 = str_replace(name_v3, "open flowers", "opened flowers")) %>%
  mutate(name_v3 = str_replace(name_v3, "unopen flowers", "unopened flowers")) %>%
  mutate(name_v3 = str_replace(name_v3, "microspore", "microspores")) %>%
  mutate(label_n = c(1:n()))

#calculate the sum of signals of epigenetic marks
epi_mark_related_replicates <- vector("list", length = nrow(epi_mark_intersect_file_list))

for (i in seq_along(epi_mark_related_replicates)) {
  epi_mark_related_replicates[[i]] <- read_delim(str_c("./data/Fig1/", epi_mark_intersect_file_list$name[[i]]), delim = "\t", col_names = c("Chr", "str", "end", "intensity", "Chr_100k", "str_100k", "end_100k", "chr_bin", "bin_size")) %>%
    mutate(name = epi_mark_intersect_file_list$name_v2[[i]]) %>%
    group_by(name, chr_bin) %>%
    summarise(sum_intensity = sum((end - str) * intensity))
}

#add the name used in figures later
sum_epi_mark_related_replicates <- bind_rows(epi_mark_related_replicates) %>%
  mutate(name_v2 = name) %>%
  ungroup() %>%
  select(-name) %>%
  left_join(epi_mark_intersect_file_list) %>%
  select(chr_bin, sum_intensity, name_light = name_v3, label_n) %>%
  split(.$label_n)

#create the function for producing figures
explained_var_predicted_single_feature_SP_Fig_1 <- function(data, data2){
  data_x <- den_table_f_cor_SV[[11]] %>%
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
    )) %>%
    left_join(data) %>%
    replace_na(list(name_light = data2, sum_intensity = 0)) %>%
    mutate(den_feature = sum_intensity/(end-str))
  
  x_location = 0.23*(quantile(data_x$den_feature, 0.975) - quantile(data_x$den_feature, 0.025)) + quantile(data_x$den_feature, 0.025)
  
  x_limits = c(quantile(data_x$den_feature, 0.025), quantile(data_x$den_feature, 0.975))
  
  
  prediction <- smooth.spline(data_x$den_feature, data_x$RecRate, spar = 1.2)
  predicted_y <- predict(prediction, data_x$den_feature)$y
  explained_var_vector <- round(1-sum((data_x$RecRate-predicted_y)^2)/sum((data_x$RecRate-mean(data_x$RecRate))^2), 2)
  
  data_x %>% ggscatter("den_feature", "RecRate", color = "status", palette = c(pal_npg("nrc", alpha = 1)(6))[c(1,3,4)]) +
    geom_spline(
      aes(den_feature, RecRate),
      color = "black",
      size  = 1,
      spar = 1
    ) +
    facet_wrap( ~ name_light, nrow = 1, scales = "free") +
    theme_bw() +
    labs(y = "Recombination rate (cM/Mb)", x = "The density of features") +
    theme(
      #axis.ticks = element_blank(),
      strip.text.x = element_text(
        colour = "black",
        face = "bold",
        size = 16
      ),
      legend.text = element_text(size = 12, face = "bold"),
      legend.title = element_blank(),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 14),
      legend.position="bottom"
    ) + 
    theme(axis.title.y = element_text(color = "white"), axis.title.x = element_text(color = "white"),  legend.position = "", plot.margin = margin(0,0.5,0,0.5, "cm")) +
    geom_text(x =  x_location, y = 15, label = str_c("R^2 == ", explained_var_vector), parse = TRUE, size = 6) +
    scale_x_continuous(limits = c(x_limits[[1]], x_limits[[2]]))
}

#create the list for producing separate plots of 6 epigenetic features
SP1_fig_list <- vector("list", length(epi_mark_intersect_file_list$name_v3))

for (i in seq_along(epi_mark_intersect_file_list$name_v3)) {
  SP1_fig_list[[i]] <- explained_var_predicted_single_feature_SP_Fig_1(sum_epi_mark_related_replicates[[i]], epi_mark_intersect_file_list$name_v3[[i]])
}


#create SP1 by combining plots (4 plots (sources) for each feature)
SP1_f_raw_H3K4me1 <- list(SP1_fig_list[[1]], SP1_fig_list[[2]], SP1_fig_list[[3]], SP1_fig_list[[4]])
SP1_f_raw_H3K4me3 <- list(SP1_fig_list[[5]], SP1_fig_list[[6]], SP1_fig_list[[7]], SP1_fig_list[[8]])
SP1_f_raw_H3K9me2 <- list(SP1_fig_list[[9]], SP1_fig_list[[10]], SP1_fig_list[[11]], SP1_fig_list[[12]])
SP1_f_raw_H3K27me3 <- list(SP1_fig_list[[13]], SP1_fig_list[[14]], SP1_fig_list[[15]], SP1_fig_list[[16]])
SP1_f_raw_ATAC <- list(SP1_fig_list[[17]], SP1_fig_list[[18]], SP1_fig_list[[19]], SP1_fig_list[[20]])
SP1_f_raw_DNase <- list(SP1_fig_list[[21]], SP1_fig_list[[22]], SP1_fig_list[[23]], SP1_fig_list[[24]])
epi_marks_labels <- c("H3K4me1", "H3K4me3", "H3K9me2", "H3K27me3", "ATAC", "DNase")
SP1_f_raw_epi_feature <- list(SP1_f_raw_H3K4me1, SP1_f_raw_H3K4me3, SP1_f_raw_H3K9me2, SP1_f_raw_H3K27me3, SP1_f_raw_ATAC, SP1_f_raw_DNase)
SP_Fig1_sub_lab <- c("A", "B", "C", "D", "E", "F")

SP_Fig1_part_epi_marks <- vector("list", length = length(epi_marks_labels))


for (i in seq_along(epi_marks_labels)) {
  
  raw <- ggarrange(SP1_f_raw_epi_feature[[i]][[1]], SP1_f_raw_epi_feature[[i]][[2]], SP1_f_raw_epi_feature[[i]][[3]], SP1_f_raw_epi_feature[[i]][[4]], nrow = 2, ncol = 2) +
    theme(plot.margin = margin(0.5,0,0,0, "cm"), plot.background = element_rect(fill = "white", color = "white"))  
  
  SP_Fig1_part_epi_marks[[i]] <- annotate_figure(raw, top = text_grob(epi_marks_labels[[i]],  color = "black", face = "bold", size = 26),
                                                 bottom = text_grob("The density of features", color = "black", size = 16, face = "bold"),
                                                 left = text_grob("Recombination rate (cM/Mb)", color = "black", rot = 90, size = 16, face = "bold"),
                                                 fig.lab = SP_Fig1_sub_lab[[i]], fig.lab.face = "bold", fig.lab.size = 16)  +
    theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"), plot.background = element_rect(fill = "white", color = "white")) 
}

SP_Fig1 <- ggarrange(SP_Fig1_part_epi_marks[[1]], SP_Fig1_part_epi_marks[[2]], SP_Fig1_part_epi_marks[[3]], SP_Fig1_part_epi_marks[[4]], SP_Fig1_part_epi_marks[[5]], SP_Fig1_part_epi_marks[[6]],nrow = 3, ncol = 2)

ggsave("./analysis/SP_Fig1.jpeg", SP_Fig1, width = 480, height = 640, units = c("mm"))
