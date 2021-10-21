#trick
#one wants to split up one figure into two in a quick way
#https://pinetools.com/split-image

setwd("/data/projects/thesis/INRA_project/Ara_TE_task/R_markdown/Model_1st/chromatin_state_model/")

#import necessary packages
Packages <- c("scales", "tidyverse", "ggrepel", "ggsci", "ggpubr", "doMC", "doParallel", "foreach", "slider", "cowplot", "ggpmisc", "combinat")
lapply(Packages, library, character.only = TRUE)

#this analysis is to check whether different sizes of bins will give us different optimized recrate
Ara_Chr_label <- vector(mode = "list", length = 5)
bin_size <- c(500, 1000, 2000, 3000, 4000, 5000, 10000, 15000, 20000, 50000, 100000, 200000, 500000, 1000000)
bin_name <-
  c("0_5k", "1k", "2k", "3k", "4k", "5k", "10k", "15k", "20k", "50k", "100k", "200k", "500k", "1000k")
Chr_label = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

Ara_genome_bed <-
  read_delim(
    "./data/Fig4/Ara_genome_bed",
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
      size_l = rep(bin_size[[size]], length(Ara_Chr_label[[chr]][[size]]))
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

path_den_table <- str_c("/data/projects/thesis/INRA_project/Ara_TE_task/R_markdown/Model_1st/",
                        "den_table_",
                        bin_name,
                        "_bed")

pwalk(list(den_table_list, path_den_table), write_delim, delim = "\t", col_names = FALSE)


#read CO file and calculate recombination rate
paths_den_table_RCO <-
  str_c(
    "./data/Fig4/",
    "den_table_",
    bin_name,
    "_RCO_raw_bed"
  )

#produce the file of the sum of CO (Rowan)
den_table_CO_sum = vector("list", length = length(bin_name))

for (i in seq_along(bin_name)) {
  den_table_CO_sum[[i]] <-
    read_delim(
      paths_den_table_RCO[[i]],
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
    summarise(sum_CO = sum(CO_n), .groups = "drop")  
}



#produce the file of the sum of CO (Ian)
cross_combination <- str_c(c("Col_Bur", "Col_Ws", "Col_Clc", "Col_Ct", "Col_Ler"), "_raw")

paths_den_table_ICO_combined <-
  str_c(
    "./data/Fig4/",
    "den_table_",
    bin_name,
    "_ICO_raw_bed"
  )

den_table_ICO_sum = vector("list", length = length(bin_name))

for (i in seq_along(bin_name)) {
  den_table_ICO_sum[[i]] <-
    read_delim(
      paths_den_table_ICO_combined[[i]],
      delim = "\t",
      col_names = c(
        "Chr",
        "str",
        "end",
        "chr_bin",
        "Chr_CO",
        "str_CO",
        "end_CO",
        "mid_CO",
        "size",
        "cross"
      )
    ) %>%
    mutate(CO_n = (end - str) / (end_CO - str_CO)) %>%
    mutate(cross = factor(.$cross, levels = cross_combination)) %>%
    group_by(chr_bin, cross) %>%
    summarise(sum_ICO = sum(CO_n), .groups = "drop") %>%
    split(.$cross)
}

den_table_f_all_CO = vector("list", length(bin_name))

for (i in seq_along(bin_name)) {
  den_table_f_all_CO[[i]] <- den_table_list[[i]] %>%
    left_join(den_table_CO_sum[[i]]) %>%
    left_join(den_table_ICO_sum[[i]]$Col_Bur_raw) %>%
    mutate(sum_ICO_Bur = sum_ICO) %>%
    select(-cross, -sum_ICO) %>%
    left_join(den_table_ICO_sum[[i]]$Col_Ct_raw) %>%
    mutate(sum_ICO_Ct = sum_ICO) %>%
    select(-cross, -sum_ICO) %>%
    left_join(den_table_ICO_sum[[i]]$Col_Ler_raw) %>%
    mutate(sum_ICO_Ler = sum_ICO) %>%
    select(-cross, -sum_ICO) %>%
    left_join(den_table_ICO_sum[[i]]$Col_Ws_raw) %>%
    mutate(sum_ICO_Ws = sum_ICO) %>%
    select(-cross, -sum_ICO) %>%
    left_join(den_table_ICO_sum[[i]]$Col_Clc_raw) %>%
    mutate(sum_ICO_Clc = sum_ICO) %>%
    select(-cross, -sum_ICO) %>%
    replace_na(list(sum_ICO_Ler = 0, sum_ICO_Bur = 0, sum_ICO_Ws = 0, sum_ICO_Ct = 0, sum_ICO_Clc = 0)) %>%
    mutate(Recrate_Ler = sum_ICO_Ler/2/245/(end-str)*10^8, Recrate_Clc = sum_ICO_Clc/2/189/(end-str)*10^8, Recrate_Ct = sum_ICO_Ct/2/305/(end-str)*10^8, Recrate_Ws = sum_ICO_Ws/2/188/(end-str)*10^8, Recrate_Bur = sum_ICO_Bur/2/180/(end-str)*10^8) %>%
    mutate(Recrate_I_equal = (Recrate_Clc+Recrate_Ws+Recrate_Ct+Recrate_Bur+Recrate_Ler)/5, Recrate_I_noLer = (Recrate_Clc+Recrate_Ws+Recrate_Ct+Recrate_Bur)/4) %>%
    #select(sum_ICO_Bur, sum_ICO_Clc, sum_ICO_Ct, sum_ICO_Ler, sum_ICO_Ws) %>%
    replace_na(list(sum_CO = 0, sum_TSS = 0, sum_5_prime_UTR = 0)) %>%
    mutate(Recrate_Rowan = sum_CO / 2 / 2182 / (end - str) * 10 ^ 8, size_bin = end-str) %>%
    group_by(Chr) %>%
    mutate(
      genetic_length = sum(Recrate_Rowan * (end - str) * 10 ^ -6),
      genetic_length_Clc = sum(Recrate_Clc*size_bin/1000000),
      physical_length = (end - str) * 10 ^ -6
    ) %>%
    ungroup() %>%
    mutate(var_rec_Rowan = Recrate_Rowan/2182/2/(end-str)*10^8, var_rec_Ler = Recrate_Ler/245/2/(end-str)*10^8, var_rec_Clc = Recrate_Clc/189/2/(end-str)*10^8, var_rec_Ct = Recrate_Ct/305/2/(end-str)*10^8, var_rec_Ws = Recrate_Ws/188/2/(end-str)*10^8, var_rec_Bur = Recrate_Bur/180/2/(end-str)*10^8) %>%
    mutate(var_rec_Ian = (var_rec_Clc+var_rec_Ct+var_rec_Ws+var_rec_Bur)/16) %>%
    mutate(var_rec_df_Rowan_Ian = var_rec_Rowan+var_rec_Ian) %>%
    mutate(sd_rec_df_Rowan_Ian = sqrt(var_rec_df_Rowan_Ian))
}

#To make a table with ICO
den_table_f_CO_Ian = vector("list", length(bin_name))

for (i in seq_along(bin_name)) {
  den_table_f_CO_Ian[[i]] <- den_table_f_all_CO[[i]] %>%
    select(Chr, str, end, chr_bin, Recrate_Ler_Rowan = Recrate_Rowan, Recrate_Ler, Recrate_Clc, Recrate_Ct, Recrate_Bur, Recrate_Ws, Recrate_I_equal) %>%
    gather(key = "cross", value = "Recrate", c("Recrate_Ler_Rowan","Recrate_Ler", "Recrate_Clc", "Recrate_Ct", "Recrate_Bur", "Recrate_Ws")) %>%
    mutate(cross = str_replace(cross, "Recrate", "Col")) %>%
    mutate(source = if_else(cross == "Col_Ler_Rowan", "Rowan", "Ian")) %>%
    mutate(cross = if_else(cross == "Col_Ler_Rowan", "Col_Ler", cross))
}


paths_den_table_ISNP_state <-
  str_c(
    "./data/Fig4/",
    "den_table_",
    bin_name,
    "_state_10_ISNP_intersect"
  )

paths_den_table_ISNP_s_state <-
  str_c(
    "./data/Fig4/",
    "den_table_",
    bin_name,
    "_state_10_ISNP_s_intersect"
  )

den_table_ISNP_sum = vector(mode = "list", length = length(bin_name))
den_table_ISNP_s_sum = vector(mode = "list", length = length(bin_name))

for (size in seq_along(bin_name)) {
  den_table_ISNP_sum[[size]] <- read_delim(
    paths_den_table_ISNP_state[[size]],
    delim = "\t",
    col_names = c("Chr", "str", "end", "feature", "label", "str_bin", "end_bin", "chr_bin", "SNP_type")
  ) %>%
    group_by(chr_bin, SNP_type) %>%
    summarise(sum_SNP = n()) %>%
    mutate(SNP_type = str_replace(SNP_type, "SNP_", "")) %>%
    select(chr_bin, cross = SNP_type, sum_SNP)
  
  den_table_ISNP_s_sum[[size]] <- read_delim(
    paths_den_table_ISNP_s_state[[size]],
    delim = "\t",
    col_names = c("Chr", "str", "end", "feature", "label", "str_bin", "end_bin", "chr_bin", "SNP_type")
  ) %>%
    group_by(chr_bin, SNP_type) %>%
    summarise(sum_SNP_s = n()) %>%
    mutate(SNP_type = str_replace(SNP_type, "SNP_", ""))%>%
    select(chr_bin, cross = SNP_type, sum_SNP_s)
}


den_table_f_CO_SNP_Ian_100kb <- den_table_f_CO_Ian[[11]] %>%
  left_join(den_table_ISNP_sum[[11]]) %>%
  left_join(den_table_ISNP_s_sum[[11]]) %>%
  replace_na(list(sum_SNP = 0, sum_SNP_s = 0)) %>%
  mutate(den_SNP_bin = 1000*sum_SNP/(end - str), den_SNP_s_bin = 1000*sum_SNP_s/(end-str)) %>%
  mutate(diff_Recrate = Recrate - Recrate_I_equal) %>%
  ungroup() %>%
  mutate(cross_f = str_c(cross, "_", source)) %>%
  split(.$cross_f)


#The result obviously showed that both SNP dataset sources showed the similar quadratic effect in the changing CO rate

Ara_cen_posi <- tibble(
  Chr = str_c("Chr", 1:5),
  str_cen = c(13920001, 2950001, 12680001, 3390001, 10950001),
  end_cen = c(15970000, 4750000, 14750000, 4820000, 13240000)
)
#the masked regions are the Col genetic background in Clc 
masked_region <- tibble(
  Chr = c("Chr2", "Chr2", "Chr4"),
  masked_str = c(7000000, 16500000, 12500000),
  masked_end = c(10000000, 18500000, 18500000)
)

den_table_f_CO_SNP_Ian_100kb_combined_masked_region <- bind_rows(den_table_f_CO_SNP_Ian_100kb) %>%
  left_join(masked_region) %>%
  drop_na() %>%
  mutate(overlapped_masked = if_else(end <= masked_str | str >= masked_end, "no", "yes")) %>%
  group_by(Chr, str, end, chr_bin) %>%
  mutate(overlapped_masked_value = if_else(overlapped_masked == "yes", 1, 0)) %>%
  summarise(sum_overlapped_masked_value = sum(overlapped_masked_value)) %>%
  arrange(Chr, str, end) %>%
  mutate(overlapped_masked = if_else(sum_overlapped_masked_value != 0, "yes", "no")) %>%
  select(Chr, str, end, chr_bin, overlapped_masked)

Ian_pop_num <- tibble(cross_f = c(str_c(c("Col_Bur", "Col_Ler", "Col_Ws", "Col_Ct", "Col_Clc"),"_Ian"), "Col_Ler_Rowan"),
                      pop_n = c(180, 245, 188, 305, 189, 2182))

den_table_f_CO_SNP_Ian_100kb_combined <- bind_rows(den_table_f_CO_SNP_Ian_100kb) %>%
  left_join(Ian_pop_num) %>%
  mutate(size_bin = (end-str)) %>%
  left_join(den_table_f_CO_SNP_Ian_100kb_combined_masked_region) %>%
  replace_na(list(overlapped_masked = "no")) %>%
  filter(overlapped_masked == "no") %>%
  group_by(cross_f) %>%
  mutate(quantile_den_SNP = if_else(den_SNP_bin <= quantile(den_SNP_bin, 0.5), 1, 2)) %>%
  ungroup()


den_table_f_CO_SNP_Ian_100kb_combined_noCO_removed_list <- den_table_f_CO_SNP_Ian_100kb_combined %>%
  filter(Recrate !=0) %>%
  split(.$cross_f) 

names(den_table_f_CO_SNP_Ian_100kb_combined_noCO_removed_list)

nls_model_Ian <- vector("list", length = length(Ian_pop_num$cross_f))


for (i in seq_along(nls_model_Ian)) {
  nls_model_Ian[[i]] <- summary(nls(
    Recrate ~ (intercept+linear*den_SNP_bin)*exp(-ex*den_SNP_bin),
    data = den_table_f_CO_SNP_Ian_100kb_combined_noCO_removed_list[[i]],
    start = list(intercept = 1, linear = 0.1, ex = 0.02)))
}


nls_model_estimate <- bind_rows(nls_model_Ian[[1]]$coefficients[,1], nls_model_Ian[[2]]$coefficients[,1], nls_model_Ian[[3]]$coefficients[,1], nls_model_Ian[[4]]$coefficients[,1], nls_model_Ian[[5]]$coefficients[,1], nls_model_Ian[[6]]$coefficients[,1]) %>%
  mutate(cross_f = names(den_table_f_CO_SNP_Ian_100kb_combined_noCO_removed_list)) 

nls_model_f <- tibble(den_SNP_bin = rep(seq(0,25,0.5),6),
                      cross_f = rep(c(str_c(c("Col_Bur", "Col_Ler", "Col_Ws", "Col_Ct", "Col_Clc"),"_Ian"), "Col_Ler_Rowan"), each = 51)) %>%
  left_join(nls_model_estimate) %>%
  mutate(fit_final = (intercept+linear*den_SNP_bin)*exp(-ex*den_SNP_bin)) %>%
  mutate(cross_f = str_replace(cross_f, "_Rowan", " (Rowan et al.)")) %>%
  mutate(cross_f = str_replace(cross_f, "_Ian", " (Blackwell et al.)"))

den_table_f_CO_SNP_Ian_100kb_combined_noCO_removed_list_f <- bind_rows(den_table_f_CO_SNP_Ian_100kb_combined_noCO_removed_list) %>%
  mutate(cross_f = str_replace(cross_f, "_Rowan", " (Rowan et al.)")) %>%
  mutate(cross_f = str_replace(cross_f, "_Ian", " (Blackwell et al.)"))

den_table_f_CO_SNP_Ian_100kb_combined_noCO_removed_list_f$cross_f <- factor(den_table_f_CO_SNP_Ian_100kb_combined_noCO_removed_list_f$cross_f, levels = c("Col_Ler (Rowan et al.)", str_c("Col_", c("Ler", "Bur", "Clc", "Ws", "Ct"), " (Blackwell et al.)")))

nls_model_f$cross_f <- factor(nls_model_f$cross_f, levels = c("Col_Ler (Rowan et al.)", str_c("Col_", c("Ler", "Bur", "Clc", "Ws", "Ct"), " (Blackwell et al.)")))

Fig4_top <- den_table_f_CO_SNP_Ian_100kb_combined_noCO_removed_list_f %>%
  ggplot() +
  geom_point(aes(den_SNP_bin, Recrate)) +
  geom_line(aes(den_SNP_bin, fit_final), color = "red", data = nls_model_f) +
  facet_wrap(~ cross_f, nrow = 2, scales = "free")+
  theme_bw() +
  labs(y = "recombination rate (cM/Mb)", x = "The density of SNPs (counts/kb)") +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), 
        axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 16, face = "bold"))



ISNP_total <- den_table_f_CO_SNP_Ian_100kb_combined %>%
  select(chr_bin, cross, den_SNP_bin, cross_f) %>%
  filter(cross_f != "Col_Ler_Rowan") %>%
  select(-cross_f) %>%
  mutate(cross = str_c("den_SNP_bin_", cross)) %>%
  spread(key = "cross", value = den_SNP_bin)



#based on combn and the vector I produced, I made 1024 kinds of combination for the calculation of chi-square values
Bur_cross <- c("1_Clc", "1_Ct", "1_Ler", "1_Ws")
Clc_cross <- c("2_Bur", "2_Ct", "2_Ler", "2_Ws")
Ws_cross <- c("3_Bur", "3_Ct", "3_Ler", "3_Clc")
Ct_cross <- c("4_Bur", "4_Ws", "4_Ler", "4_Clc")
Ler_cross <- c("5_Bur", "5_Ct", "5_Ws", "5_Clc")

#test "combn" using 2 levels of vectors
as_tibble(t(combn(c(Bur_cross, Clc_cross, Ws_cross, Ct_cross, Ler_cross), 2))) %>%
  View()


#produce 1024 combinations of shuffling SNP from 5 crosses 
shuffling_cross_combn <- as_tibble(t(combn(c(Bur_cross, Clc_cross, Ws_cross, Ct_cross, Ler_cross), 5))) %>%
  mutate(order_1 = as.double(str_sub(V1, 1, 1)), order_2 = as.double(str_sub(V2, 1, 1)), order_3 = as.double(str_sub(V3, 1, 1)), order_4 = as.double(str_sub(V4, 1, 1)), order_5 = as.double(str_sub(V5, 1, 1))) %>%
  filter(order_1 == 1 & order_2 == 2 & order_3 == 3 & order_4 == 4 & order_5 == 5)


#produce the list of the combination for the reshuffling test
shuffling_cross_combn_list <- shuffling_cross_combn %>%
  mutate(V1 = str_replace(V1, "1", "Col"), V2 = str_replace(V2, "2", "Col"), V3 = str_replace(V3, "3", "Col"), V4 = str_replace(V4, "4", "Col"), V5 = str_replace(V5, "5", "Col")) %>%
  select(order_1 = V1, order_2 = V2, order_3 = V3, order_4 = V4, order_5 = V5) %>%
  mutate(label = rep(c(1:(n()/8)), each = 8)) %>%
  mutate(label_2 = rep(c(1:8), times = n()/8)) %>%
  gather(key = "order", value = "chi_square_cross", c(str_c("order_", seq(1:5)))) %>%
  mutate(order = str_remove(order, "order_")) %>%
  mutate(order = as.double(order)) %>%
  arrange(label, order) %>%
  select(chi_square_cross, order, label, label_2) %>%
  split(.$label) %>%
  map(. %>% split(.$label_2))




#produce the chi-square value using the first 50% of shuffling SNP density. (exponential) (For the chi-square of reshuffled ones)
#(ex. chi_square_Bur means this value is calculated through the application of the SNP density of Bur)
sum_chi_square_SNP_effect_first_q_50 <- function(a, data1, data2){
  x <- data1 %>%
    left_join(ISNP_total) %>%
    mutate(constant = 1) %>%
    mutate(linear = a[1]) %>%
    mutate(exponential = a[2]) %>%
    mutate(SNP_effect_total_Bur = (linear*den_SNP_bin_Col_Bur+constant)*exp(-exponential*den_SNP_bin_Col_Bur)) %>%
    mutate(SNP_effect_total_Ct = (linear*den_SNP_bin_Col_Ct+constant)*exp(-exponential*den_SNP_bin_Col_Ct)) %>%
    mutate(SNP_effect_total_Ws = (linear*den_SNP_bin_Col_Ws+constant)*exp(-exponential*den_SNP_bin_Col_Ws)) %>%
    mutate(SNP_effect_total_Clc = (linear*den_SNP_bin_Col_Clc+constant)*exp(-exponential*den_SNP_bin_Col_Clc)) %>%
    mutate(SNP_effect_total_Ler = (linear*den_SNP_bin_Col_Ler+constant)*exp(-exponential*den_SNP_bin_Col_Ler)) %>%
    group_by(chr_bin) %>%
    ungroup() %>%
    select(Chr, str, end, chr_bin, size_bin, pop_n, cross, Recrate, SNP_effect_total_Bur, SNP_effect_total_Ct, SNP_effect_total_Clc, SNP_effect_total_Ws, SNP_effect_total_Ler) %>%
    gather(key = "SNP_effect_total_cross", value = "SNP_effect_total_value", str_c("SNP_effect_total_", c("Bur", "Ct", "Ws", "Clc", "Ler"))) %>%
    mutate(chi_square_cross = str_c("Col_", str_sub(SNP_effect_total_cross, 18, 21))) %>%
    filter(cross != SNP_effect_total_cross) %>%
    mutate(order = if_else(cross == "Col_Bur", 1, if_else(cross == "Col_Clc", 2, if_else(cross == "Col_Ws", 3, if_else(cross == "Col_Ct", 4, 5))))) %>%
    left_join(den_table_f_CO_SNP_Ian_100kb_combined_masked_region) %>%
    replace_na(list(overlapped_masked = "no")) %>%
    filter(overlapped_masked == "no") %>%
    left_join(data2) %>%
    drop_na() %>%
    group_by(chr_bin) %>%
    arrange(Chr, chr_bin) %>%
    mutate(base_r_parameter = sqrt(sum(Recrate^2)/sum(SNP_effect_total_value^2))) %>%
    ungroup() %>%
    mutate(chi_square = (Recrate-base_r_parameter*SNP_effect_total_value)^2/(10^2*base_r_parameter*SNP_effect_total_value/(size_bin*pop_n*2*10^-6))) %>%
    mutate(chi_square = if_else(base_r_parameter == 0, 0, chi_square)) %>%
    summarise(sum_chi_square = sum(chi_square))
  
  x$sum_chi_square
}


#only used 5 pops of Ian for the optimization (events with larger 50 % of SNP density are removed)
den_table_f_CO_SNP_Ian_100kb_combined_v2 <- den_table_f_CO_SNP_Ian_100kb_combined %>%
  filter(cross_f != "Col_Ler_Rowan") %>%
  filter(quantile_den_SNP == 1)


#function 50% quantile of SNP density new (exponential) (for the reference chi-square)
ref_SNP_effect_ex_v2_first_50_q <- function(a, data){
  x <- data %>%
    mutate(constant = 1) %>%
    mutate(linear = a[1]) %>%
    mutate(exponential = a[2]) %>%
    mutate(SNP_effect_total = (linear*den_SNP_bin+constant)*exp(-exponential*den_SNP_bin)) %>%
    mutate(SNP_effect_total = if_else(SNP_effect_total <= 10^-4, 10^-4, SNP_effect_total)) %>%
    group_by(chr_bin) %>%
    mutate(base_r_parameter = sqrt(sum(Recrate^2)/sum(SNP_effect_total^2))) %>%
    mutate(Recrate_sum = sum(Recrate)) %>%
    ungroup() %>%
    #filter(Recrate_sum!=0) %>%
    mutate(chi_square = (Recrate-base_r_parameter*SNP_effect_total)^2/(10^2*base_r_parameter*SNP_effect_total/(size_bin*pop_n*2*10^-6))) %>%
    mutate(chi_square = if_else(base_r_parameter == 0, 0, chi_square)) %>%
    summarise(sum_chi_square = sum(chi_square))
  
  x$sum_chi_square
}

#perform optimization for the reference chi-square
op_ISNP_effect_ex_v4_first_50_q <- optim(
  rep(1,2),
  ref_SNP_effect_ex_v2_first_50_q,
  data = den_table_f_CO_SNP_Ian_100kb_combined_v2,
  control = list(maxit=800, REPORT=1, trace=6),
  method = c("L-BFGS-B"),
  lower = rep(10^-10, 2)
)


#the optimization for 1024 combinations of shuffling SNP density (the first 50% quantile)
#exponential
registerDoParallel(cores = 8)
getDoParWorkers()

op_sum_chi_square_SNP_effect_shffuling_first_50q_threads_v2 <- 
  foreach(i=c(1:128), .packages = c("tidyverse")) %:%
  foreach(j=1:8, .packages = c("tidyverse")) %dopar% {
    optim(
      c(1, 1),
      sum_chi_square_SNP_effect_first_q_50,
      data1 = den_table_f_CO_SNP_Ian_100kb_combined_v2,
      data2 = shuffling_cross_combn_list[[i]][[j]],
      control = list(maxit=800, REPORT=1, trace=6),
      method = c("L-BFGS-B"),
      lower = rep(10^-10, 2)
    )
  }

#check whether upper limit can give us different result
op_sum_chi_square_SNP_effect_shffuling_first_50q_threads_v3 <- 
  foreach(i=c(1:2), .packages = c("tidyverse")) %:%
  foreach(j=1:8, .packages = c("tidyverse")) %dopar% {
    optim(
      c(10^-9, 10^-9),
      sum_chi_square_SNP_effect_first_q_50,
      data1 = den_table_f_CO_SNP_Ian_100kb_combined_v2,
      data2 = shuffling_cross_combn_list[[i]][[j]],
      control = list(maxit=800, REPORT=1, trace=6),
      method = c("L-BFGS-B"),
      lower = rep(10^-10, 2),
      upper = rep(10^-8, 2)
    )
  }


##the final chi-square result based on the final "sum_chi_square_SNP_effect_first_q_50" with three parameters
sum_chi_square_list_first_50q_v2 <- vector("list", length = length(op_sum_chi_square_SNP_effect_shffuling_first_50q_threads_v2))

for (i in seq_along(sum_chi_square_list_first_50q_v2)) {
  for (j in seq_along(c(1:8))) {
    sum_chi_square_list_first_50q_v2[[i]][[j]] <- tibble(sum_chi_square = op_sum_chi_square_SNP_effect_shffuling_first_50q_threads_v2[[i]][[j]]$value, 
                                                         label = i,
                                                         label_2 = j,
                                                         constant = 1,
                                                         linear = op_sum_chi_square_SNP_effect_shffuling_first_50q_threads_v2[[i]][[j]]$par[[1]],
                                                         exponential = op_sum_chi_square_SNP_effect_shffuling_first_50q_threads_v2[[i]][[j]]$par[[2]],
                                                         convergence = op_sum_chi_square_SNP_effect_shffuling_first_50q_threads_v2[[i]][[j]]$convergence) 
  }
}

sum_chi_square_list_first_50q_base_r_v2_t <- bind_rows(sum_chi_square_list_first_50q_v2)
write_delim(sum_chi_square_list_first_50q_base_r_v2_t, "analysis/sum_chi_square_list_first_50q_base_r_v2_t", delim = "\t", col_names = TRUE)


sum_chi_square_fg <- function(data1, data2){
  x <- bind_rows(data1) %>%
    mutate(label = seq(1:n())) %>%
    mutate(optim_result = data2$value) %>% 
    #filter(sum_chi_square < optim_result) 
    ggplot() + 
    geom_histogram(aes(sum_chi_square), bins = 30) +
    geom_vline(xintercept = data2$value, color = "red") +
    theme_bw() +
    theme(
      strip.text.x = element_text(
        colour = "black",
        face = "bold",
        size = 20
      ),
      legend.text = element_text(size = 8, face = "bold"),
      legend.title = element_blank(),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 12)
    ) +
    labs(x = "chi-square", y = "counts")
  
  x
}

#produce the bottom of Fig 4
Fig4_bottom <- sum_chi_square_fg(sum_chi_square_list_first_50q_base_r_v2_t, op_ISNP_effect_ex_v4_first_50_q) +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), 
        axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 16, face = "bold"))




Fig4_f <- ggarrange(Fig4_top, Fig4_bottom, nrow = 2)

ggsave("./analysis/Fig4_final.jpeg", Fig4_f, width = 330, height = 300, units = c("mm"), dpi = 320)


#SNP effect in the segments of states
#----------------------------------------------------------------------------------------------------#
state_10_raw_order_new_final <- read_delim("./data/Fig4/state_10_raw_order_new_final", delim = "\t", col_names = c("Chr", "str", "end", "state", "label")) 

state_10_raw_order_new_final_sum_SNP <- read_delim("./data/Fig4/state_10_raw_order_new_final_ISNP_Col_Ler", delim = "\t", col_names = c("Chr", "str", "end", "state", "label")) %>%
  group_by(label) %>%
  summarise(sum_SNP = n())

state_10_raw_order_new_final_sum_CO <- read_delim("./data/state_10_raw_order_new_final_RCO", delim = "\t", col_names = c("Chr", "str", "end", "state", "label", "Chr_CO", "str_CO", "end_CO", "sel_420", "CO_l")) %>%
  mutate(CO_n = (end-str)/(end_CO-str_CO)) %>%
  group_by(label) %>%
  summarise(sum_CO = sum(CO_n))


library(ggformula)




scatterplot_state_pooled <- state_10_raw_order_new_final %>%
  left_join(state_10_raw_order_new_final_sum_SNP) %>%
  left_join(state_10_raw_order_new_final_sum_CO) %>%
  replace_na(list(sum_SNP = 0, sum_CO = 0)) %>%
  mutate(den_SNP_kb = sum_SNP/(end-str)*1000, CO_rate = sum_CO/(end-str)/2182/2*10^8) %>%
  mutate(state = if_else(state == "state1" | state == "state3" | state == "state6" | state == "state7", "state1367", if_else(state == "state2" | state == "state4" | state == "state5", "state245", state))) %>%
  filter(state != "SV" & state != "state9")

scatterplot_state_pooled_l <- scatterplot_state_pooled %>%
  split(.$state)

scatterplot_state_pooled_l$state8

nls_scatterplot_state_pooled <- vector("list", length = length(seq_along(scatterplot_state_pooled_l)))

for (i in seq_along(nls_scatterplot_state_pooled)) {
  nls_scatterplot_state_pooled[[i]] <- summary(nls(
    CO_rate ~ (intercept+linear*den_SNP_kb)*exp(-ex*den_SNP_kb),
    data = scatterplot_state_pooled_l[[i]],
    start = list(intercept = 0, linear = 0.1, ex = 0.02)))
}

nls_state_model_estimate_SNP <- bind_rows(nls_scatterplot_state_pooled[[1]]$coefficients[,1], nls_scatterplot_state_pooled[[2]]$coefficients[,1], nls_scatterplot_state_pooled[[3]]$coefficients[,1]) %>%
  mutate(state = c("state1367", "state245", "state8")) 


nls_state_model_f_SNP <- tibble(den_SNP_kb = rep(seq(0,20,0.5),3),
                                state = rep(c("state1367", "state245", "state8"), each = 41)) %>%
  left_join(nls_state_model_estimate_SNP) %>%
  mutate(fit_final = (intercept+linear*den_SNP_kb)*exp(-ex*den_SNP_kb)) %>%
  split(.$state)

state8_scatter_SNP_CO <- scatterplot_state_pooled %>%
  #filter(den_SNP_kb != 0) %>%
  filter(state == "state8") %>%
  filter(den_SNP_kb <= 8.89) %>%
  ggplot() +
  geom_point(aes(den_SNP_kb, CO_rate)) +
  geom_line(aes(den_SNP_kb, fit_final), color = "red", data = nls_state_model_f_SNP$state8) +
  xlim(0, 9)+
  ylim(0, 10)+
  facet_wrap(~state, scales = "free", nrow = 5) +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), panel.spacing.x = unit(0,"line"), 
        axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 18, hjust = 0.5, vjust = 1, face = "bold")) +
  xlab("The SNP density (counts/kb)") +
  ylab("Recombination rate (cM/Mb)")


state1367_scatter_SNP_CO <- scatterplot_state_pooled %>%
  #filter(den_SNP_kb != 0) %>%
  filter(state == "state1367") %>%
  filter(den_SNP_kb <= 1.09) %>%
  ggplot() +
  geom_point(aes(den_SNP_kb, CO_rate)) +
  geom_line(aes(den_SNP_kb, fit_final), color = "red", data = nls_state_model_f_SNP$state1367) +
  xlim(0, 1.5)+
  ylim(0, 10)+
  facet_wrap(~state, scales = "free", nrow = 5)+
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), panel.spacing.x = unit(0,"line"), 
        axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 18, hjust = 0.5, vjust = 1, face = "bold")) +
  xlab("The SNP density (counts/kb)") +
  ylab("Recombination rate (cM/Mb)")

state245_scatter_SNP_CO <- scatterplot_state_pooled %>%
  #filter(den_SNP_kb != 0) %>%
  filter(state == "state245") %>%
  filter(den_SNP_kb <= 2.46) %>%
  ggplot() +
  geom_point(aes(den_SNP_kb, CO_rate)) +
  geom_line(aes(den_SNP_kb, fit_final), color = "red", data = nls_state_model_f_SNP$state245) +
  xlim(0, 3)+
  ylim(0, 10)+
  facet_wrap(~state, scales = "free", nrow = 5) +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), panel.spacing.x = unit(0,"line"), 
        axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 18, hjust = 0.5, vjust = 1, face = "bold")) +
  xlab("The SNP density (counts/kb)") +
  ylab("Recombination rate (cM/Mb)")

scatterplot_state_pooled %>%
  #filter(den_SNP_kb != 0) %>%
  group_by(state) %>%
  summarise(SNP_q50 = quantile(den_SNP_kb, 0.5))

scatterplot_state_pooled_SNP_fitted <- ggarrange(state8_scatter_SNP_CO, state1367_scatter_SNP_CO, state245_scatter_SNP_CO, ncol = 3)
ggsave("./analysis/scatterplot_state_pooled_SNP_fitted.jpeg", scatterplot_state_pooled_SNP_fitted, width = 500, height = 250, units = c("mm"), dpi = 320)