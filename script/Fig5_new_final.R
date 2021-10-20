#In the final version of the model, we stayed with 10 states without pooling some states together as intra or intergenic states
#and we also let IR and SNP effect of different states be the same

#import necessary packages
Packages <- c("scales", "tidyverse", "ggrepel", "ggsci", "ggpubr", "doMC", "doParallel", "foreach", "slider", "cowplot", "ggpmisc")
lapply(Packages, library, character.only = TRUE)

#import data (the segment of state data and CO data)
#state data in df binning table
#Every state segment with intersected IR info
paths_den_table_IR_state_Ler_IR_seg <-
  str_c(
    "/data/projects/thesis/INRA_project/Ara_TE_task/R_markdown/Model_1st/paper_draft_analysis/data/modeling_IR_SNP_all_IR_event/",
    "den_table_",
    bin_name,
    "_state_10_Ler_IR_intersect"
  )

#Every state segment with SNP info (state 8 already modified, 9:8:9 are switched into 9:9:9)
paths_den_table_IR_state_Ler_IR_ISNP_seg <-
  str_c(
    "/data/projects/thesis/INRA_project/Ara_TE_task/R_markdown/Model_1st/paper_draft_analysis/data/modeling_IR_SNP_all_IR_event/",
    "den_table_",
    bin_name,
    "_state_10_Ler_IR_ISNP_intersect"
  )

den_table_ISNP_IR_state_Ler_sep_seg <- vector("list", length = length(bin_size))

for (size in seq_along(bin_name)) {
  
  x <- read_delim(paths_den_table_IR_state_Ler_IR_ISNP_seg[[size]], delim = "\t", col_names = c("Chr", "str", "end", "feature", "size_IR", "IR_type", "str_bin", "end_bin", "chr_bin", "ISNP")) %>%
    group_by(Chr, str, end, feature) %>%
    summarise(sum_ISNP = n()) %>%
    mutate(den_ISNP_seg = 1000*sum_ISNP/(end - str)) %>%
    select(-sum_ISNP)
  
  
  
  den_table_ISNP_IR_state_Ler_sep_seg[[size]] <- read_delim(paths_den_table_IR_state_Ler_IR_seg[[size]], delim = "\t", col_names = c("Chr", "str", "end", "feature", "size_IR", "IR_type", "str_bin", "end_bin", "chr_bin")) %>%
    left_join(x, by = c("Chr", "str", "end", "feature")) %>%
    replace_na(list(den_ISNP_seg = 0)) %>%
    mutate(size = (end - str)) %>%
    #mutate(den_SNP_seg_0_1k = den_SNP_seg*0.1, den_SNP_seg_0_3k = den_SNP_seg*0.3, den_SNP_seg_0_5k = den_SNP_seg*0.5, den_SNP_seg_3k = den_SNP_seg*3, den_SNP_seg_5k = den_SNP_seg*5, den_SNP_seg_10k = den_SNP_seg*10) %>%
    group_by(chr_bin) %>%
    mutate(size_bin = sum(end-str)) %>%
    arrange(Chr, str) %>%
    ungroup() 
}

#calculate CO rate in df binning table
Ara_Chr_label <- vector(mode = "list", length = 5)

bin_size <- c(500, 1000, 2000, 3000, 4000, 5000, 10000, 15000, 20000, 50000, 100000, 200000, 500000, 1000000)

bin_name <-
  c("0_5k", "1k", "2k", "3k", "4k", "5k", "10k", "15k", "20k", "50k", "100k", "200k", "500k", "1000k")
Chr_label = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

Ara_genome_bed <-
  read_delim(
    "/data/projects/thesis/INRA_project/Ara_TE_task/R_markdown/Model_1st/paper_draft_analysis/data/Fig4/Ara_genome_bed",
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
    "/data/projects/thesis/INRA_project/Ara_TE_task/R_markdown/Model_1st/paper_draft_analysis/data/Fig4/",
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
    "/data/projects/thesis/INRA_project/Ara_TE_task/R_markdown/Model_1st/paper_draft_analysis/data/Fig4/",
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

######modeling (adding SNP/IR effect in state9/SV)######
#####make the function for optimization
####the function for running SNP/IR effect (10 states, same effect for these two effects)###
state_9_same_SNP_IR_para_modu_v1_consistent <-  function(a, data, data2) {
  para_d <- tibble(
    feature = c(str_c("state", c(1:9)), "SV"),
    raw_rec = c(a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10]),
    SNP_para1 = rep(a[11], 10),
    SNP_para2 = rep(a[12], 10),
    IR_para1 = rep(a[13], 10),
    IR_para2 = rep(a[14], 10),
    IR_para3 = rep(a[15], 10)
  )
  
  x = data %>%
    left_join(para_d, by = c("feature")) %>%
    #mutate(SNP_para1 = as.double(SNP_para1), SNP_para2 = as.double(SNP_para2)) %>%
    mutate(raw_recrate_v2 = if_else(size_IR == 0, raw_rec*(1 + SNP_para1*den_ISNP_seg)*exp(-SNP_para2*den_ISNP_seg)*size/size_bin, raw_rec*(1 + SNP_para1*den_ISNP_seg)*exp(-SNP_para2*den_ISNP_seg)/(IR_para1 + IR_para2*exp(-IR_para3*size_IR/1000))*size/size_bin)) %>%
    group_by(chr_bin) %>%
    summarise(recrate_non_rescalling = sum(raw_recrate_v2)) %>%
    mutate(Chr = str_sub(chr_bin, 1, 4), label = str_sub(chr_bin, 6, 20)) %>%
    mutate(label = as.double(label)) %>%
    arrange(Chr, label)
  
  x2 <- data2 %>%
    left_join(x, by = c("Chr", "chr_bin")) %>%
    drop_na() %>%
    group_by(Chr) %>%
    mutate(pre_genetic_length = sum(recrate_non_rescalling*physical_length)) %>%
    mutate(rescaled_model_rec_rates = recrate_non_rescalling*genetic_length/pre_genetic_length) 
  
  diff <- x2$Recrate_Rowan - x2$rescaled_model_rec_rates
  sum(diff^2)
}

###run optimizations
#create initial parameters
state_9_same_SNP_3_IR_seg_rob_initial_par <- vector(mode = "list", length = 8)

for (i in c(1:8)) {
  state_9_same_SNP_3_IR_seg_rob_initial_par[[i]] <- c(1 * (0.5 + runif(15)))
}

#perform optimization (50k ~ 500k)
registerDoParallel(cores = 8)
getDoParWorkers()

test_foreach_modeling_same_SNP_3_IR_10_states_10_500k_all_IR_consistent <-
  foreach(i=c(10:13), .packages = c("tidyverse")) %:%
  foreach(j=1:8, .packages = c("tidyverse")) %dopar% {
    optim(
      state_9_same_SNP_3_IR_seg_rob_initial_par[[j]],
      state_9_same_SNP_IR_para_modu_v1_consistent,
      data = den_table_ISNP_IR_state_Ler_sep_seg[[i]],
      data2 = den_table_f_all_CO[[i]],
      control = list(maxit=800, REPORT=1, trace=6),
      method = c("L-BFGS-B"),
      #upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.9999, 0.9999, 0.9999, 0.9999, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
      lower = rep(10^-8, 15))
  }
test_foreach_modeling_same_SNP_3_IR_10_states_10_500k_all_IR_consistent[[4]][[1]]$value

#rescaling optimized CO rate of 10 states of 6 bins
rescaling_state_9_seg_same_SNP_3_IR_modu_v2_test_rec <- function(a, data, data2) {
  para_d <- tibble(
    feature = c(str_c("state", c(1:9)), "SV"),
    raw_rec = c(a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10]),
    SNP_para1 = rep(a[11], 10),
    SNP_para2 = rep(a[12], 10),
    IR_para1 = rep(a[13], 10),
    IR_para2 = rep(a[14], 10),
    IR_para3 = rep(a[15], 10)
  )
  
  x = data %>%
    left_join(para_d, by = c("feature")) %>%
    mutate(raw_recrate_v2 = if_else(size_IR == 0, raw_rec*(1 + SNP_para1*den_ISNP_seg)*exp(-SNP_para2*den_ISNP_seg)*size/size_bin, raw_rec*(1 + SNP_para1*den_ISNP_seg)*exp(-SNP_para2*den_ISNP_seg)/(IR_para1 + IR_para2*exp(-IR_para3*size_IR/1000))*size/size_bin)) %>%
    group_by(chr_bin) %>%
    summarise(recrate_non_rescalling = sum(raw_recrate_v2)) %>%
    mutate(Chr = str_sub(chr_bin, 1, 4), label = str_sub(chr_bin, 6, 20)) %>%
    mutate(label = as.double(label)) %>%
    arrange(Chr, label)
  
  x2 <- data2 %>%
    left_join(x, by = c("Chr", "chr_bin")) %>%
    drop_na() %>%
    group_by(Chr) %>%
    mutate(pre_genetic_length = sum(recrate_non_rescalling*physical_length)) %>%
    mutate(rescaled_model_rec_rates = recrate_non_rescalling*genetic_length/pre_genetic_length) 
  
  x3 <- x2 %>%
    group_by(Chr) %>%
    mutate(pred_genetic_length = sum(recrate_non_rescalling*size_bin/1000000)) %>%
    select(Chr, genetic_length, pred_genetic_length, physical_length) %>%
    group_by(Chr, genetic_length, pred_genetic_length) %>%
    mutate(physical_length_Chr = sum(physical_length), rescaling_factor = genetic_length/pred_genetic_length) %>%
    select(-physical_length) %>%
    distinct() %>%
    ungroup() %>%
    mutate(weighted_rescaling_factor = rescaling_factor*physical_length_Chr/sum(physical_length_Chr)) %>%
    mutate(sum_weighted_rescaling_factor = sum(weighted_rescaling_factor))
  
  
  c(a[1:10]*x3$sum_weighted_rescaling_factor[[1]], a[11:length(a)])
  
}

final_op_para_modified_9_state_3IR_SNP_consistent_first_trial <- vector("list", length = 4)

for (i in seq_along(1:4)) {
  for (j in seq_along(1:8)) {
    final_op_para_modified_9_state_3IR_SNP_consistent_first_trial[[i]][[j]] <- rescaling_state_9_seg_same_SNP_3_IR_modu_v2_test_rec(test_foreach_modeling_same_SNP_3_IR_10_states_10_500k_all_IR_consistent[[i]][[j]]$par, den_table_ISNP_IR_state_Ler_sep_seg[[i+9]], den_table_f_all_CO[[i+9]])
  }
}

transform_para_tibble_10_states <- function(data){
  x = tibble(
    para = c(str_c("state", c(1:9)), "SV", "SNP_para1", "SNP_para2", str_c("IR_para", c(1:3))),
    r1 = data[[1]],
    r2 = data[[2]], 
    r3 = data[[3]], 
    r4 = data[[4]],
    r5 = data[[5]],
    r6 = data[[6]], 
    r7 = data[[7]],
    r8 = data[[8]]
  )
  
  x
}

final_op_para_modified_same_3IR_SNP_10k_500k_10_states_consistent_rep1 <- final_op_para_modified_9_state_3IR_SNP_consistent_first_trial %>%
  map(. %>% transform_para_tibble_10_states()) %>%
  bind_rows() %>%
  mutate(bin = rep(bin_size[10:13], each = 15))


write_delim(final_op_para_modified_same_3IR_SNP_10k_500k_10_states_consistent_rep1, "./analysis/final_op_para_modified_same_3IR_SNP_10k_500k_10_states_consistent_rep1", delim = "\t", col_names = TRUE)

final_op_para_modified_same_3IR_SNP_10k_500k_10_states_consistent_rep1_l <- final_op_para_modified_same_3IR_SNP_10k_500k_10_states_consistent_rep1 %>%
  split(.$bin)

#####Fig 5: make new landscape plot based on 100-kb bins
#create the modified function that can show rescaled CO rate (least square 10 states/same SNP and IR effect, including state9/SV)
state_9_same_SNP_IR_para_modu_v2_rec <-  function(a, data, data2) {
  para_d <- tibble(
    feature = c(str_c("state", c(1:9)), "SV"),
    raw_rec = c(a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10]),
    SNP_para1 = rep(a[11], 10),
    SNP_para2 = rep(a[12], 10),
    IR_para1 = rep(a[13], 10),
    IR_para2 = rep(a[14], 10),
    IR_para3 = rep(a[15], 10)
  )
  
  x = data %>%
    left_join(para_d, by = c("feature")) %>%
    #mutate(SNP_para1 = as.double(SNP_para1), SNP_para2 = as.double(SNP_para2)) %>%
    mutate(raw_recrate_v2 = if_else(size_IR == 0, raw_rec*(1 + SNP_para1*den_ISNP_seg)*exp(-SNP_para2*den_ISNP_seg)*size/size_bin, raw_rec*(1 + SNP_para1*den_ISNP_seg)*exp(-SNP_para2*den_ISNP_seg)/(IR_para1 + IR_para2*exp(-IR_para3*size_IR/1000))*size/size_bin)) %>%
    group_by(chr_bin) %>%
    summarise(recrate_non_rescalling = sum(raw_recrate_v2)) %>%
    mutate(Chr = str_sub(chr_bin, 1, 4), label = str_sub(chr_bin, 6, 20)) %>%
    mutate(label = as.double(label)) %>%
    arrange(Chr, label)
  
  x2 <- data2 %>%
    left_join(x, by = c("Chr", "chr_bin")) %>%
    drop_na() %>%
    group_by(Chr) %>%
    mutate(pre_genetic_length = sum(recrate_non_rescalling*physical_length)) %>%
    mutate(rescaled_model_rec_rates = recrate_non_rescalling*genetic_length/pre_genetic_length) 
  
  x2$rescaled_model_rec_rates
}

#create the table for later using the middle of arms to produce the inset of landscape
Ara_peri_posi <- tibble(
  Chr = str_c("Chr", 1:5),
  str_peri = c(11420001, 910001, 10390001, 1070001, 8890001),
  end_peri = c(18270000, 7320000, 16730000, 6630000, 15550000)
)

Fig5_tag <- Ara_genome %>% 
  left_join(Ara_peri_posi) %>%
  mutate(mid_arm = (chr_end+end_peri)/2000000)

#create the list of main and inset figures (Chr separated) (using r5 in 100kb of final_op_para_modified_same_3IR_SNP_10k_500k_10_states_consistent_rep1)
Fig5_main_v2 <- vector("list", length = length(Fig5_tag$Chr))
Fig5_inset_v2 <- vector("list", length = length(Fig5_tag$Chr))
Fig5_f_v2 <- vector("list", length = length(Fig5_tag$Chr))

test_table <- den_table_f_all_CO[[11]] %>%
  mutate(predicted_Recrate = state_9_same_SNP_IR_para_modu_v2_rec(final_op_para_modified_same_3IR_SNP_10k_500k_10_states_consistent_rep1_l$`1e+05`$r5, den_table_ISNP_IR_state_Ler_sep_seg[[11]], den_table_f_all_CO[[11]]))

R_square_Fig5 <- 1-sum((test_table$predicted_Recrate-test_table$Recrate_Rowan)^2)/sum((test_table$Recrate_Rowan-mean(test_table$Recrate_Rowan))^2)

#https://stackoverflow.com/questions/5219671/it-is-possible-to-create-inset-graphs
Fig5_tag$Chr[[1]]
for (i in seq_along(Fig5_tag$Chr)) {
  
  data_Fig5 <- den_table_f_all_CO[[11]] %>%
    mutate(predicted_Recrate = state_9_same_SNP_IR_para_modu_v1_rec(final_op_para_modified_same_3IR_SNP_10k_500k_10_states_consistent_rep1_l$`1e+05`$r5, den_table_ISNP_IR_state_Ler_sep_seg[[11]], den_table_f_all_CO[[11]])) %>%
    #filter(Chr == "Chr1" | Chr == "Chr3" | Chr == "Chr5") %>%
    filter(Chr == Fig5_tag$Chr[[i]])
  
  R_square_Fig5 <- round(1-sum((data_Fig5$predicted_Recrate-data_Fig5$Recrate_Rowan)^2)/sum((data_Fig5$Recrate_Rowan-mean(data_Fig5$Recrate_Rowan))^2), 2)
  x_location_Fig5 = 0.1*(max(data_Fig5$end)-min(data_Fig5$str))
  y_location_Fig5 = 0.8*(max(data_Fig5$Recrate_Rowan) %/% 1 + 1) 
  
  
  Fig5_main_v2[[i]] <-  data_Fig5 %>%
    gather(key = "Recrate_type", value = "Recrate", c("Recrate_Rowan", "predicted_Recrate")) %>%
    mutate(Recrate_type = if_else(Recrate_type == "Recrate_Rowan", "Experimental", "Predicted")) %>%
    ggplot() +
    geom_line(aes((str+end)/2000000, Recrate, color = Recrate_type), size = 1) +
    facet_wrap(~Chr, nrow = 5) +
    theme_bw() +
    labs(y = "Recombination rate (cM/Mb)", x = "Mb") +
    scale_color_manual(values=c(pal_npg("nrc", alpha = 1)(2))) +
    geom_text(x =  x_location_Fig5/1000000, y =y_location_Fig5, label = str_c("R^2 == ", R_square_Fig5), parse = TRUE, size = 6) +
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 16, face = "bold"),
          legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), panel.spacing.x = unit(0,"line"), 
          axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 18, face = "bold"), legend.position = c(0.1, 0.92), legend.box.background = element_rect(colour = "black")) 
  
  Fig5_inset_v2[[i]] <- data_Fig5 %>%
    gather(key = "Recrate_type", value = "Recrate", c("Recrate_Rowan", "predicted_Recrate")) %>%
    mutate(Recrate_type = if_else(Recrate_type == "Recrate_Rowan", "Experimental", "Predicted_state_SNP_IR")) %>%
    ggplot() +
    geom_line(aes((str+end)/2000000, Recrate, color = Recrate_type), size = 0.8) +
    #facet_wrap(~Chr, nrow = 5) +
    scale_color_manual(values=c(pal_npg("nrc", alpha = 1)(2))) +
    xlim(c(Fig5_tag$mid_arm[[i]]-3, Fig5_tag$mid_arm[[i]]+3))+
    theme_bw() +
    labs(y = "", x = "") +
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
          legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), panel.spacing.x = unit(0,"line"), 
          axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 18, face = "bold"), legend.position = "") 
  
  
  Fig5_f_v2[[i]] <- ggdraw() +
    draw_plot(Fig5_main_v2[[i]]) +
    draw_plot(Fig5_inset_v2[[i]], x = 0.62, y = .55, width = .35, height = .35)
}

path_Fig5_v2 <- str_c("Fig5_Chr", seq(1,5), "_new_model_v2.jpeg")

pwalk(list(path_Fig5_v2, Fig5_f_v2), ggsave, path = "./analysis/", width = 320, height = 192, units = c("mm"), dpi = 320)

Fig5_f_v2_gra_abs_test <- Fig5_main_v2[[1]] +
  theme(legend.position = c(0.8, 0.85))

ggsave("./analysis/Fig5_Chr1_new_model_gra_abs_v1.jpeg", Fig5_f_v2_gra_abs_test, width = 150, height = 100, units = c("mm"), dpi = 320)


#select the best optimization result for producing Fig S7
SP_table_op_para_combined_v2 <- tibble(
  para = final_op_para_modified_same_3IR_SNP_10k_500k_10_states_rep2_list[[1]]$para,
  para_50k = final_op_para_modified_9_state_3IR_SNP_consistent_first_trial[[1]][[4]],
  para_100k = final_op_para_modified_9_state_3IR_SNP_consistent_first_trial[[2]][[5]],
  para_200k = final_op_para_modified_9_state_3IR_SNP_consistent_first_trial[[3]][[7]],
  para_500k = final_op_para_modified_9_state_3IR_SNP_consistent_first_trial[[4]][[5]]
)


write_csv(SP_table_op_para_combined_v2, "./analysis/SP_table_op_para_combined_v2.csv", col_names = TRUE)

####produce the scatter plot for figure S7 (state9/SV SNP/IR effect added)####
den_table_f_scatter_plot_raw_rec = vector("list", length = length(c(1:4)))

#using the best case from the optimization
optimized_par_50_500k <- list(final_op_para_modified_9_state_3IR_SNP_consistent_first_trial[[1]][[4]], final_op_para_modified_9_state_3IR_SNP_consistent_first_trial[[2]][[5]], final_op_para_modified_9_state_3IR_SNP_consistent_first_trial[[3]][[7]], final_op_para_modified_9_state_3IR_SNP_consistent_first_trial[[4]][[5]])

#create the table for plots
for (i in c(1:4)) {
  den_table_f_scatter_plot_raw_rec[[i]] <- den_table_f_all_CO[[i+9]] %>%
    left_join(Ara_peri_posi, by = "Chr") %>%
    mutate(status = if_else(end < str_peri | str > end_peri, "arms", if_else(str > str_peri & end < end_peri, "pericentromeric_regions", "overlapped"))) %>%
    mutate(predicted_Recrate =state_9_same_SNP_IR_para_modu_v2_rec(optimized_par_50_500k[[i]], den_table_ISNP_IR_state_Ler_sep_seg[[i+9]], den_table_f_all_CO[[i+9]]))
}

#ploting formula
den_table_f_scatter_plot_raw_rec <- den_table_f_scatter_plot_raw_rec %>%  
  map(. %>% mutate(size_l_text = str_c(size_l/1000, "kb"))) 


scatter_plots_pre_ex <- function(data){  
  y <- data$predicted_Recrate
  x <- data$Recrate_Rowan
  explained_var_vector = round(1 - sum((x-y)^2)/sum((mean(x)-x)^2), 2)
  x_location = 0.2*(max(data$Recrate_Rowan) %/% 1 + 1) 
  y_location = 0.9*(max(data$Recrate_Rowan) %/% 1 + 1) 
  
  figure <- data %>%
    #ggscatter("RecRate", "predicted_Recrate", color = "status") +
    ggscatter("predicted_Recrate", "Recrate_Rowan", color = "status", palette = c(pal_npg("nrc", alpha = 1)(6))[c(1,3,4)]) +
    #stat_cor(method = "pearson", label.x =10,  label.y = 1, cor.coef.name = "R", size = 10, digits = 3) +
    scale_x_continuous(limits = c(0, max(data$Recrate_Rowan) %/% 1 + 1), breaks = seq(0, 5 * (max(data$Recrate_Rowan) %/% 5 + 1), 5)) +
    scale_y_continuous(limits = c(0, max(data$Recrate_Rowan) %/% 1 + 1), breaks = seq(0, 5 * (max(data$Recrate_Rowan) %/% 5 + 1), 5)) +
    geom_text(x =  x_location, y =y_location, label = str_c("R^2 == ", explained_var_vector), parse = TRUE, size = 6) +
    geom_abline(slope = 1, size = 1) +
    facet_wrap(~size_l_text, nrow = 2, scales = "free_x") +
    theme_bw() +
    labs(x = "Predicted recombination rate (cM/Mb)", y = "Experimental recombination rate (cM/Mb)") +
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
          legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), panel.spacing.x = unit(0,"line"), 
          axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 18, face = "bold"), legend.position = c(0.8, 0.15), legend.box.background = element_rect(colour = "black")) 
  
  figure
}


Fig_SX <- ggarrange(scatter_plots_pre_ex(den_table_f_scatter_plot_raw_rec[[1]]), scatter_plots_pre_ex(den_table_f_scatter_plot_raw_rec[[2]]), scatter_plots_pre_ex(den_table_f_scatter_plot_raw_rec[[3]]), scatter_plots_pre_ex(den_table_f_scatter_plot_raw_rec[[4]]), nrow = 2, ncol = 2)

Fig_SX_v2 <- annotate_figure(Fig_SX, 
                             bottom = text_grob("Predicted recombination rate (cM/Mb)", color = "black", size = 18, face = "bold"),
                             left = text_grob("Experimental recombination rate (cM/Mb)", color = "black", rot = 90, size = 18, face = "bold")) +
  theme(plot.margin = margin(0,0.5,0.5,0, "cm"), plot.background = element_rect(fill = "white", color = "white")) 

ggsave("./analysis/Fig_SX_new.png", Fig_SX_v2, width = 360, height = 280, units = c("mm"))


#perform op from Chr1 to Chr5 separately, and predict other chromosomes using 100-kb bins
#perform optimization (50k ~ 500k)
registerDoParallel(cores = 8)
getDoParWorkers()

#separate 5 Chrs
den_table_ISNP_IR_state_Ler_sep_seg_Chr_100k <- den_table_ISNP_IR_state_Ler_sep_seg[[11]] %>%
  split(.$Chr)

den_table_f_all_CO_Chr_100k <- den_table_f_all_CO[[11]] %>%
  split(.$Chr)

den_table_f_all_CO_Chr_100k[[1]]

test_foreach_modeling_same_SNP_3_IR_10_states_10_500k_all_IR_consistent_Chr_100k <-
  foreach(i=c(1:5), .packages = c("tidyverse")) %:%
  foreach(j=1:8, .packages = c("tidyverse")) %dopar% {
    optim(
      state_9_same_SNP_3_IR_seg_rob_initial_par[[j]],
      state_9_same_SNP_IR_para_modu_v1_consistent,
      data = den_table_ISNP_IR_state_Ler_sep_seg_Chr_100k[[i]],
      data2 = den_table_f_all_CO_Chr_100k[[i]],
      control = list(maxit=800, REPORT=1, trace=6),
      method = c("L-BFGS-B"),
      #upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.9999, 0.9999, 0.9999, 0.9999, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
      lower = rep(10^-8, 15))
  }

state_9_same_SNP_IR_para_modu_v1_consistent_rec <- function(a, data, data2) {
  para_d <- tibble(
    feature = c(str_c("state", c(1:9)), "SV"),
    raw_rec = c(a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10]),
    SNP_para1 = rep(a[11], 10),
    SNP_para2 = rep(a[12], 10),
    IR_para1 = rep(a[13], 10),
    IR_para2 = rep(a[14], 10),
    IR_para3 = rep(a[15], 10)
  )
  
  x = data %>%
    left_join(para_d, by = c("feature")) %>%
    #mutate(SNP_para1 = as.double(SNP_para1), SNP_para2 = as.double(SNP_para2)) %>%
    mutate(raw_recrate_v2 = if_else(size_IR == 0, raw_rec*(1 + SNP_para1*den_ISNP_seg)*exp(-SNP_para2*den_ISNP_seg)*size/size_bin, raw_rec*(1 + SNP_para1*den_ISNP_seg)*exp(-SNP_para2*den_ISNP_seg)/(IR_para1 + IR_para2*exp(-IR_para3*size_IR/1000))*size/size_bin)) %>%
    group_by(chr_bin) %>%
    summarise(recrate_non_rescalling = sum(raw_recrate_v2)) %>%
    mutate(Chr = str_sub(chr_bin, 1, 4), label = str_sub(chr_bin, 6, 20)) %>%
    mutate(label = as.double(label)) %>%
    arrange(Chr, label)
  
  x2 <- data2 %>%
    left_join(x, by = c("Chr", "chr_bin")) %>%
    drop_na() %>%
    group_by(Chr) %>%
    mutate(pre_genetic_length = sum(recrate_non_rescalling*physical_length)) %>%
    mutate(rescaled_model_rec_rates = recrate_non_rescalling*genetic_length/pre_genetic_length) 
  
  x2$rescaled_model_rec_rates
}


op_para_Chr_list <- list(
  test_foreach_modeling_same_SNP_3_IR_10_states_10_500k_all_IR_consistent_Chr_100k[[1]][[3]]$par,
  test_foreach_modeling_same_SNP_3_IR_10_states_10_500k_all_IR_consistent_Chr_100k[[2]][[4]]$par,
  test_foreach_modeling_same_SNP_3_IR_10_states_10_500k_all_IR_consistent_Chr_100k[[3]][[1]]$par,
  test_foreach_modeling_same_SNP_3_IR_10_states_10_500k_all_IR_consistent_Chr_100k[[4]][[7]]$par,
  test_foreach_modeling_same_SNP_3_IR_10_states_10_500k_all_IR_consistent_Chr_100k[[5]][[5]]$par
)

R_square_Chr_list <- vector("list", length = length(op_para_Chr_list))

for (i in seq_along(R_square_Chr_list)) {
  for (j in seq_along(R_square_Chr_list)) {
    predicted_rec <- state_9_same_SNP_IR_para_modu_v1_consistent_rec(op_para_Chr_list[[i]], den_table_ISNP_IR_state_Ler_sep_seg_Chr_100k[[j]], den_table_f_all_CO_Chr_100k[[j]])
    R_square_Chr_list[[i]][[j]] <- 1 - sum((den_table_f_all_CO_Chr_100k[[j]]$Recrate_Rowan - predicted_rec)^2)/sum((den_table_f_all_CO_Chr_100k[[j]]$Recrate_Rowan - mean(den_table_f_all_CO_Chr_100k[[j]]$Recrate_Rowan))^2)
  }
}
R_square_Chr_list[[5]]

R_square_based_on_op_Chr_table <- tibble(
  Chr = str_c("Chr", c(1:5), "_predict"),
  Chr1_fit = round(R_square_Chr_list[[1]], 3),
  Chr2_fit = round(R_square_Chr_list[[2]], 3),
  Chr3_fit = round(R_square_Chr_list[[3]], 3),
  Chr4_fit = round(R_square_Chr_list[[4]], 3),
  Chr5_fit = round(R_square_Chr_list[[5]], 3)
)

write_delim(R_square_based_on_op_Chr_table, "./analysis/R_square_based_on_op_Chr_table_v2", delim = "\t", col_names = TRUE)
