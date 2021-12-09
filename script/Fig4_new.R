#change the directory "chromatin_state_model" as the working directory (the link below is an example)
#setwd("/data/projects/thesis/INRA_project/Ara_TE_task/R_markdown/Model_1st/chromatin_state_model/")

#using "p_load" from the package "pacman" to install and load necessary packages
install.packages("pacman", repos = "https://mirror.ibcp.fr/pub/CRAN/")
library(pacman)

Packages <- c("scales", "tidyverse", "ggrepel", "ggsci", "ggpubr", "doMC", "doParallel", "foreach", "slider", "cowplot", "combinat")
p_load(Packages, character.only = TRUE)

#lapply(Packages, library, character.only = TRUE)

#create the initial bed files using different bins for producing bed files intersecting features
Ara_Chr_label <- vector(mode = "list", length = 5)
bin_size <- c(500, 1000, 2000, 3000, 4000, 5000, 10000, 15000, 20000, 50000, 100000, 200000, 500000, 1000000)
bin_name <-
  c("0_5k", "1k", "2k", "3k", "4k", "5k", "10k", "15k", "20k", "50k", "100k", "200k", "500k", "1000k")
Chr_label = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

Ara_genome_bed <-
  read_delim(
    "./data/Fig2/Ara_genome_bed",
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

write_delim(den_table_list[[11]], "./data/den_table_100k_bed", delim = "\t", col_names = FALSE)


#produce the table of SNP from 5 F2 populations from Blakewell's research
SNP_tag_Ian <- tibble(
  cross_raw = c(1:5),
  cross = str_c("SNP_Col_", c("Ct", "Ws", "Bur", "Clc", "Ler"))
)

#run the command below in the directory "./data/Fig4" to procude the decompressed bed file
system(paste("cd ./data/Fig4", "&& gunzip -c Ian_pop_passed_SNP_bed_raw.gz > Ian_pop_passed_SNP_bed_raw", sep = " "))

Ian_pop_passed_SNP_bed <- read_delim("./data/Fig4/Ian_pop_passed_SNP_bed_raw", col_names = c("Chr", "str", "end", "cross_raw"), delim = "\t") %>%
  mutate(Chr = str_c("Chr", Chr)) %>%
  left_join(SNP_tag_Ian) %>%
  select(-cross_raw)

write_delim(Ian_pop_passed_SNP_bed, "./data/Fig4/Ian_pop_passed_SNP_bed", col_names = FALSE, delim = "\t")


#run the shell script "Fig4_CO_SNP_intersection.sh" in the directory "script/" to procude the bed file for the intersection between CO/SNP and 100-kb bins.
system(paste("cd ./script", "&& bash Fig4_CO_SNP_intersection.sh", sep = " "))

#produce the file of the sum of Rowans' COs
den_table_100k_RCO_sum <- read_delim(
  "./data/Fig4/den_table_100k_RCO_raw_mid_bed",
  delim = "\t",
  col_names = c(
    "Chr",
    "str",
    "end",
    "chr_bin"
  )
) %>%
  mutate(CO_n = (end - str)) %>%
  group_by(chr_bin) %>%
  summarise(sum_CO = sum(CO_n), .groups = "drop")

#produce the file of the sum of Ians' COs
cross_combination <- str_c(c("Col_Bur", "Col_Ws", "Col_Clc", "Col_Ct", "Col_Ler"), "_raw")

den_table_100k_ICO_sum <-
    read_delim(
      "./data/Fig4/den_table_100k_ICO_raw_mid_bed",
      delim = "\t",
      col_names = c(
        "Chr",
        "str",
        "end",
        "chr_bin",
        "cross"
      )
    ) %>%
  mutate(CO_n = (end - str)) %>%
  mutate(cross = factor(.$cross, levels = cross_combination)) %>%
  group_by(chr_bin, cross) %>%
  summarise(sum_ICO = sum(CO_n), .groups = "drop") %>%
  split(.$cross)

#calculate CO rate from 6 populations (5 are from Blackwell's research, and 1 is from Rowan's data)
den_table_100k_f_all_CO <- den_table_list[[11]] %>%
    left_join(den_table_100k_RCO_sum) %>%
    left_join(den_table_100k_ICO_sum$Col_Bur_raw) %>%
    mutate(sum_ICO_Bur = sum_ICO) %>%
    select(-cross, -sum_ICO) %>%
    left_join(den_table_100k_ICO_sum$Col_Ct_raw) %>%
    mutate(sum_ICO_Ct = sum_ICO) %>%
    select(-cross, -sum_ICO) %>%
    left_join(den_table_100k_ICO_sum$Col_Ler_raw) %>%
    mutate(sum_ICO_Ler = sum_ICO) %>%
    select(-cross, -sum_ICO) %>%
    left_join(den_table_100k_ICO_sum$Col_Ws_raw) %>%
    mutate(sum_ICO_Ws = sum_ICO) %>%
    select(-cross, -sum_ICO) %>%
    left_join(den_table_100k_ICO_sum$Col_Clc_raw) %>%
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


#transform the previous table with CO information into another format for the following analysis
den_table_100k_f_ICO <- den_table_100k_f_all_CO %>%
    select(Chr, str, end, chr_bin, Recrate_Ler_Rowan = Recrate_Rowan, Recrate_Ler, Recrate_Clc, Recrate_Ct, Recrate_Bur, Recrate_Ws, Recrate_I_equal) %>%
    gather(key = "cross", value = "Recrate", c("Recrate_Ler_Rowan","Recrate_Ler", "Recrate_Clc", "Recrate_Ct", "Recrate_Bur", "Recrate_Ws")) %>%
    mutate(cross = str_replace(cross, "Recrate", "Col")) %>%
    mutate(source = if_else(cross == "Col_Ler_Rowan", "Rowan", "Ian")) %>%
    mutate(cross = if_else(cross == "Col_Ler_Rowan", "Col_Ler", cross))

#calculate the sum of SNP of 5 F2 populations in 100-kb bins
den_table_100k_ISNP_sum <- read_delim(
  "./data/Fig4/den_table_100k_ISNP_intersect",
  delim = "\t",
  col_names = c("chr_bin", "SNP_type")
) %>%
  group_by(chr_bin, SNP_type) %>%
  summarise(sum_SNP = n(), .groups = "drop") %>%
  mutate(SNP_type = str_replace(SNP_type, "SNP_", "")) %>%
  select(chr_bin, cross = SNP_type, sum_SNP)


#create the table with CO rate and SNP density of the corresponding F2 population in 100-kb bins
den_table_100kb_f_ICO_ISNP <- den_table_100k_f_ICO %>%
  left_join(den_table_100k_ISNP_sum) %>%
  replace_na(list(sum_SNP = 0)) %>%
  mutate(den_SNP_bin = 1000*sum_SNP/(end - str)) %>%
  mutate(diff_Recrate = Recrate - Recrate_I_equal) %>%
  ungroup() %>%
  mutate(cross_f = str_c(cross, "_", source)) %>%
  split(.$cross_f)


#define the centromeric regions of Arabidopsis
Ara_cen_posi <- tibble(
  Chr = str_c("Chr", 1:5),
  str_cen = c(13920001, 2950001, 12680001, 3390001, 10950001),
  end_cen = c(15970000, 4750000, 14750000, 4820000, 13240000)
)

#the masked regions are the Col genetic background in Clc (part of regions in Chr2 and Chr4)
masked_region <- tibble(
  Chr = c("Chr2", "Chr2", "Chr4"),
  masked_str = c(7000000, 16500000, 12500000),
  masked_end = c(10000000, 18500000, 18500000)
)

den_table_100kb_f_ICO_ISNP_combined_masked_region <- bind_rows(den_table_100kb_f_ICO_ISNP) %>%
  left_join(masked_region) %>%
  drop_na() %>%
  mutate(overlapped_masked = if_else(end <= masked_str | str >= masked_end, "no", "yes")) %>%
  group_by(Chr, str, end, chr_bin) %>%
  mutate(overlapped_masked_value = if_else(overlapped_masked == "yes", 1, 0)) %>%
  summarise(sum_overlapped_masked_value = sum(overlapped_masked_value)) %>%
  arrange(Chr, str, end) %>%
  mutate(overlapped_masked = if_else(sum_overlapped_masked_value != 0, "yes", "no")) %>%
  select(Chr, str, end, chr_bin, overlapped_masked)

#create the table for the following optimization by maximizing log likelihood
Ian_pop_num <- tibble(cross_f = c(str_c(c("Col_Bur", "Col_Ler", "Col_Ws", "Col_Ct", "Col_Clc"),"_Ian"), "Col_Ler_Rowan"),
                      pop_n = c(180, 245, 188, 305, 189, 2182))

den_table_100kb_f_ICO_ISNP_combined <- bind_rows(den_table_100kb_f_ICO_ISNP) %>%
  left_join(Ian_pop_num) %>%
  mutate(size_bin = (end-str)) %>%
  left_join(den_table_100kb_f_ICO_ISNP_combined_masked_region) %>%
  replace_na(list(overlapped_masked = "no")) %>%
  filter(overlapped_masked == "no") %>%
  group_by(cross_f) %>%
  mutate(quantile_den_SNP = if_else(den_SNP_bin <= quantile(den_SNP_bin, 0.5), 1, 2)) %>%
  ungroup()


den_table_100kb_f_ICO_ISNP_combined_list <- den_table_100kb_f_ICO_ISNP_combined %>%
  split(.$cross_f) 


#create the function for producing the log likelihood using 3 parameters for SNP effect
SNP_effect_ex_f <- function(a, data){
  x <- data %>%
    mutate(intercept = a[1], linear = a[2], ex = a[3], physical_length = (end-str)/10^6, sum_CO = Recrate*physical_length*2*pop_n/100) %>%
    mutate(predicted_rec = (intercept+linear*den_SNP_bin)*exp(-ex*den_SNP_bin)) %>% 
    mutate(p_CO = if_else(predicted_rec*physical_length/100 < 1, predicted_rec*physical_length/100, 1-10^-10), p_noCO = 1-p_CO) %>%
    mutate(log_likelihood = sum_CO*log(p_CO) + (pop_n*2-sum_CO)*log(p_noCO)) 
  
  sum(x$log_likelihood)  
  
}

#perform optimization using 3 parameters for 6 pops 
registerDoParallel(cores = 6)
getDoParWorkers()

op_SNP_effect_ex <-
  foreach(j=1:6, .packages = c("tidyverse")) %dopar% {
    optim(
      c(1, 0.1, 0.02),
      SNP_effect_ex_f,
      data = den_table_100kb_f_ICO_ISNP_combined_list[[j]],
      control = list(fnscale=-1, maxit=800, REPORT=1, trace=6),
      method = c("L-BFGS-B"),
      #upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.9999, 0.9999, 0.9999, 0.9999, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
      lower = rep(10^-8, 3))
  }

#create the function for the optimization using 2 parameters for SNP effect (the linear term excluded)
SNP_effect_H0_f <- function(a, data){
  x <- data %>%
    mutate(intercept = a[1], ex = a[2], physical_length = (end-str)/10^6, sum_CO = Recrate*physical_length*2*pop_n/100) %>%
    mutate(predicted_rec = intercept*exp(-ex*den_SNP_bin)) %>% 
    mutate(p_CO = if_else(predicted_rec*physical_length/100 < 1, predicted_rec*physical_length/100, 1-10^-10), p_noCO = 1-p_CO) %>%
    mutate(log_likelihood = sum_CO*log(p_CO) + (pop_n*2-sum_CO)*log(p_noCO)) 
  
  sum(x$log_likelihood)  
  
}

#perform optimization using 2 parameters for 6 pops
op_SNP_effect_H0 <-
  foreach(j=1:6, .packages = c("tidyverse")) %dopar% {
    optim(
      c(1, 0.02),
      SNP_effect_H0_f,
      data = den_table_100kb_f_ICO_ISNP_combined_list[[j]],
      control = list(fnscale=-1, maxit=800, REPORT=1, trace=6),
      method = c("L-BFGS-B"),
      #upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.9999, 0.9999, 0.9999, 0.9999, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
      lower = rep(10^-8, 3))
  }


#compare two models using likelihood ratio test
ll_test_result <- vector("list", length = length(op_SNP_effect_H0))
p_value_ll_test <- vector("list", length = length(op_SNP_effect_H0))
SNP_effect_ex_tibble <- vector("list", length = length(op_SNP_effect_H0))

for (i in seq_along(ll_test_result)) {
  #create the ratio of likelihood ratio test
  ll_test_result[[i]] <- -2*(op_SNP_effect_H0[[i]]$value-op_SNP_effect_ex[[i]]$value)
  
  #use the critical value of chi square to get the p-value
  p_value_ll_test[[i]] <- signif(pchisq(ll_test_result[[i]], 1, lower.tail = FALSE), 2)
  
  #merge the table including the information of parameters, p-value and cross information 
  SNP_effect_ex_tibble[[i]] <- tibble(
    para_type = c("intercept", "linear", "ex"),
    para = op_SNP_effect_ex[[i]]$par) %>%
    spread(key = para_type, value = para) %>%
    mutate(cross_f = names(den_table_100kb_f_ICO_ISNP_combined_list)[[i]], 
           p_value = if_else(substr(p_value_ll_test[[i]], 2, 2) == ".", as.character(p_value_ll_test[[i]]), str_c(str_sub(p_value_ll_test[[i]], 1, 1), ".0", str_sub(p_value_ll_test[[i]], 2, nchar(p_value_ll_test[[i]]))))) 
  
}

#bind the table for creating another table for plotting fitted line
SNP_ll_model_estimate <- bind_rows(SNP_effect_ex_tibble) 


SNP_ll_model_f <- tibble(den_SNP_bin = rep(seq(0,25,0.5),6),
                         cross_f = rep(c(str_c(c("Col_Bur", "Col_Ler", "Col_Ws", "Col_Ct", "Col_Clc"),"_Ian"), "Col_Ler_Rowan"), each = 51)) %>%
  left_join(SNP_ll_model_estimate) %>%
  mutate(fit_final = (intercept+linear*den_SNP_bin)*exp(-ex*den_SNP_bin)) %>%
  mutate(cross_f = str_replace(cross_f, "_Rowan", " (Rowan et al.)")) %>%
  mutate(cross_f = str_replace(cross_f, "_Ian", " (Blackwell et al.)")) %>%
  mutate(cross_f = factor(cross_f, levels = c("Col_Ler (Rowan et al.)", str_c("Col_", c("Ler", "Bur", "Clc", "Ws", "Ct"), " (Blackwell et al.)")))) %>%
  split(.$cross_f)

#create the table for the raw scatterplot of 6 pops
den_table_100kb_f_ICO_ISNP_combined_list_f <- bind_rows(den_table_100kb_f_ICO_ISNP_combined_list) %>%
  mutate(cross_f = str_replace(cross_f, "_Rowan", " (Rowan et al.)")) %>%
  mutate(cross_f = str_replace(cross_f, "_Ian", " (Blackwell et al.)")) %>%
  mutate(cross_f = factor(cross_f, levels = c("Col_Ler (Rowan et al.)", str_c("Col_", c("Ler", "Bur", "Clc", "Ws", "Ct"), " (Blackwell et al.)")))) %>%
  split(.$cross_f)

#create Fig4
Fig4_top_ll_raw = vector("list", length = length(den_table_100kb_f_ICO_ISNP_combined_list_f))

for (i in seq_along(Fig4_top_ll_raw)) {
  #set the location for p-value from likelihood ratio test
  x_location_Fig4 = 0.7*(max(den_table_100kb_f_ICO_ISNP_combined_list_f[[i]]$end)-min(den_table_100kb_f_ICO_ISNP_combined_list_f[[i]]$str))/10^6
  y_location_Fig4 = 0.8*(max(den_table_100kb_f_ICO_ISNP_combined_list_f[[i]]$Recrate) %/% 1 + 1)
  
  Fig4_top_ll_raw[[i]] <- den_table_100kb_f_ICO_ISNP_combined_list_f[[i]] %>%
    ggplot() +
    geom_point(aes(den_SNP_bin, Recrate)) +
    geom_line(aes(den_SNP_bin, fit_final), color = "red", data = SNP_ll_model_f[[i]]) +
    facet_wrap(~ cross_f, nrow = 1, scales = "free")+
    geom_text(x =  x_location_Fig4, y = y_location_Fig4, label = str_c("p = ", SNP_ll_model_f[[i]]$p_value[[1]]), parse = FALSE, size = 5) +
    theme_bw() +
    labs(y = "recombination rate (cM/Mb)", x = "The density of SNPs (counts/kb)") +
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
          legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
          axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 16, face = "bold"))
  
  
}


Fig4_top_ll_merged <- ggarrange(Fig4_top_ll_raw[[1]], Fig4_top_ll_raw[[2]], Fig4_top_ll_raw[[3]], Fig4_top_ll_raw[[4]], Fig4_top_ll_raw[[5]], Fig4_top_ll_raw[[6]], nrow = 2, ncol = 3)

Fig4_top_ll_merged_f <- annotate_figure(Fig4_top_ll_merged, bottom = text_grob("The density of SNPs (counts/kb)", face = "bold", size = 18), left = text_grob("Recombination rate (cM/Mb)", rot = 90, face = "bold", size = 18)) +
  theme(plot.margin = margin(0,0.5,0.5,0, "cm"), plot.background = element_rect(fill = "white", color = "white")) 

ggsave("./analysis_output/Fig4_top_ll_merged_f.jpeg", Fig4_top_ll_merged_f, width = 330, height = 200, units = c("mm"), dpi = 320)

###reshuffling test for SP figure S7###
#prepare the table with the SNP density of 5 pops
ISNP_total <- den_table_100kb_f_ICO_ISNP_combined %>%
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
    left_join(den_table_100kb_f_ICO_ISNP_combined_masked_region) %>%
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


#only used 5 pops of Ian for the optimization (bins with the SNP density larger than the 50th percentile are removed)
den_table_100kb_f_ICO_ISNP_combined_v2 <- den_table_100kb_f_ICO_ISNP_combined %>%
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
  data = den_table_100kb_f_ICO_ISNP_combined_v2,
  control = list(maxit=800, REPORT=1, trace=6),
  method = c("L-BFGS-B"),
  lower = rep(10^-10, 2)
)


#the optimization for 1024 combinations of shuffling SNP density (the 50th percentile)
#exponential
registerDoParallel(cores = 8)
getDoParWorkers()

op_sum_chi_square_SNP_effect_shffuling_first_50q_threads_v2 <- 
  foreach(i=c(1:128), .packages = c("tidyverse")) %:%
  foreach(j=1:8, .packages = c("tidyverse")) %dopar% {
    optim(
      c(1, 1),
      sum_chi_square_SNP_effect_first_q_50,
      data1 = den_table_100kb_f_ICO_ISNP_combined_v2,
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
      data1 = den_table_100kb_f_ICO_ISNP_combined_v2,
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

#create and store the table with the optimization result
sum_chi_square_list_first_50q_base_r_v2_t <- bind_rows(sum_chi_square_list_first_50q_v2)
write_delim(sum_chi_square_list_first_50q_base_r_v2_t, "analysis_output/sum_chi_square_list_first_50q_base_r_mid_CO_t", delim = "\t", col_names = TRUE)

#create the function for producing SP figure S7
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

#create SP figure S7
SP_FigS7 <- sum_chi_square_fg(sum_chi_square_list_first_50q_base_r_v2_t, op_ISNP_effect_ex_v4_first_50_q) +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), 
        axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 16, face = "bold"))


ggsave("./analysis_output/SP_FigS7.jpeg", SP_FigS7, width = 330, height = 200, units = c("mm"), dpi = 320)
