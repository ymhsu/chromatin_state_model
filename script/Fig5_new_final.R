#This script intends to create all the figures and tables related to our final model with 15 parameters.
#We create Figure 5, Figure S7, Figure S8, Table S3, Table S4, and SP file2 here.

#In the full model, we have 10 parameters for each state (10 parameters), 2 parameters of SNP effect and 3 parameters for IR effect.
#note that IR/SNP effect are the same for each segment of 10 states

#change the directory "chromatin_state_model" as the working directory
setwd("/data/projects/thesis/INRA_project/Ara_TE_task/R_markdown/Model_1st/chromatin_state_model/")

#using "p_load" from the package "pacman" to install and load necessary packages
install.packages("pacman")
library(pacman)

Packages <- c("scales", "tidyverse", "ggrepel", "ggsci", "ggpubr", "doMC", "doParallel", "foreach", "slider", "cowplot", "combinat")
p_load(Packages, character.only = TRUE)

#lapply(Packages, library, character.only = TRUE)

#Set given labels for the following analysis
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

bin_name <-
  c("0_5k", "1k", "2k", "3k", "4k", "5k", "10k", "15k", "20k", "50k", "100k", "200k", "500k", "1000k")

paths_den_table <-
  str_c(
    "data/",
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

#use the original 9-chromatin state file to modify state 8 sandwiched by two state 9 segments into state 9
state_9_total_modified <- read_delim("./data/Fig5/chromatin_state_total_bed", col_names = c("Chr", "str", "end", "state"), delim = "\t") %>%
  arrange(Chr, str) %>%
  group_by(Chr) %>%
  mutate(pre_state = lag(state), aft_state = lead(state)) %>%
  mutate(state = if_else(state == "state8" & pre_state == "state9" & aft_state == "state9", "state9", state)) %>%
  select(-pre_state, -aft_state) %>%
  ungroup() 

write_delim(state_9_total_modified, "./data/Fig5/state_9_total_modified", delim = "\t", col_names = FALSE)


#Open the terminal, run the command below in the directory "./data/Fig4" to procude the decompressed bed file
#gunzip -c Ian_pop_passed_SNP_bed_raw.gz > Ian_pop_passed_SNP_bed_raw
SNP_tag_Ian <- tibble(
  cross_raw = c(1:5),
  cross = str_c("SNP_Col_", c("Ct", "Ws", "Bur", "Clc", "Ler"))
)

Ian_pop_passed_SNP_bed <- read_delim("./data/Fig4/Ian_pop_passed_SNP_bed_raw", col_names = c("Chr", "str", "end", "cross_raw"), delim = "\t") %>%
  mutate(Chr = str_c("Chr", Chr)) %>%
  left_join(SNP_tag_Ian) %>%
  select(-cross_raw)

write_delim(Ian_pop_passed_SNP_bed, "./data/Fig4/Ian_pop_passed_SNP_bed", col_names = FALSE, delim = "\t")

#Open the terminal, run the shell script "Fig5_IR_SNP_bed_file.sh" below in the directory "./script" for generating the modified 10-state segments
#Then create the 10-state segments with the information of SNP density and IR using different sizes of bins 

#Import data (the segment of state data and CO data)
#State data in tables with different size of bins
#Every state segment with intersected IR info
paths_den_table_IR_state_Ler_IR_seg <-
  str_c(
    "./data/Fig5/",
    "den_table_",
    bin_name,
    "_state_10_Ler_IR_intersect"
  )

#Every state segment with SNP info (state 8 already modified, 9:8:9 are switched into 9:9:9)
paths_den_table_IR_state_Ler_IR_ISNP_seg <-
  str_c(
    "./data/Fig5/",
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
    group_by(chr_bin) %>%
    mutate(size_bin = sum(end-str)) %>%
    arrange(Chr, str) %>%
    ungroup() 
}

#read CO file and calculate recombination rate
paths_den_table_RCO_mid <-
  str_c(
    "./data/Fig1/",
    "den_table_",
    bin_name,
    "_RCO_raw_mid_bed"
  )

#produce the file of all bins (50-500 k) with the sum of CO (Rowan)
#calculate the sum of CO
den_table_CO_mid_sum = vector("list", length = length(bin_name))

for (i in seq_along(bin_name)) {
  den_table_CO_mid_sum[[i]] <-
    read_delim(
      paths_den_table_RCO_mid[[i]],
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

#the list of different bins with recombination rate based on CO intervals of Rowan
den_table_f_all_RCO = vector("list", length(bin_name))

for (i in seq_along(bin_name)) {
  den_table_f_all_RCO[[i]] <- den_table_list[[i]] %>%
    left_join(den_table_CO_mid_sum[[i]]) %>%
    replace_na(list(sum_CO = 0)) %>%
    mutate(Recrate_Rowan = sum_CO / 2 / 2182 / (end - str) * 10 ^ 8, size_bin = end-str) %>%
    group_by(Chr) %>%
    mutate(
      genetic_length = sum(Recrate_Rowan * (end - str) * 10 ^ -6),
      physical_length = (end - str) * 10 ^ -6
    ) %>%
    ungroup() %>%
    mutate(var_rec_Rowan = Recrate_Rowan/2182/2/(end-str)*10^8) 
}

######modeling (adding SNP/IR effect in state9/SV)######
#####make the function for optimization (To maximize log likelihood score)
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
    #calculation for the full model (10 states + IR effect + SNP effect), for rescaling see 14 lines below
    mutate(raw_recrate_v2 = if_else(size_IR == 0, raw_rec*(1 + SNP_para1*den_ISNP_seg)*exp(-SNP_para2*den_ISNP_seg)*size/size_bin, raw_rec*(1 + SNP_para1*den_ISNP_seg)*exp(-SNP_para2*den_ISNP_seg)/(IR_para1 + IR_para2*exp(-IR_para3*size_IR/1000))*size/size_bin)) %>%
    #calculation for the model (10 states + SNP effect)
    #mutate(raw_recrate_v2 = raw_rec*(1 + SNP_para1*den_ISNP_seg)*exp(-SNP_para2*den_ISNP_seg)*size/size_bin) %>%
    #calculation for the model (10 states + IR effect)
    #mutate(raw_recrate_v2 = if_else(size_IR == 0, raw_rec*size/size_bin, raw_rec/(IR_para1 + IR_para2*exp(-IR_para3*size_IR/1000))*size/size_bin)) %>%
    #calculation for the model (10 states)
    #mutate(raw_recrate_v2 = raw_rec*size/size_bin) %>%
    group_by(chr_bin) %>%
    summarise(recrate_non_rescaling = sum(raw_recrate_v2)) %>%
    mutate(Chr = str_sub(chr_bin, 1, 4), label = str_sub(chr_bin, 6, 20)) %>%
    mutate(label = as.double(label)) %>%
    arrange(Chr, label)
  
  x2 <- data2 %>%
    left_join(x, by = c("Chr", "chr_bin")) %>%
    drop_na() %>%
    group_by(Chr) %>%
    mutate(pre_genetic_length = sum(recrate_non_rescaling*physical_length)) %>%
    mutate(rescaled_model_rec_rates = recrate_non_rescaling*genetic_length/pre_genetic_length) %>%
    mutate(p_CO = if_else(rescaled_model_rec_rates*physical_length/100 < 1, rescaled_model_rec_rates*physical_length/100, 1-10^-10), p_noCO = 1-p_CO) %>%
    #use the following line to bypass the genetic length rescaling
    #mutate(p_CO = if_else(recrate_non_rescaling*physical_length/100 < 1, recrate_non_rescaling*physical_length/100, 1-10^-10), p_noCO = 1-p_CO) %>%
    mutate(log_likelihood = sum_CO*log(p_CO) + (2182*2-sum_CO)*log(p_noCO))
  
  sum(x2$log_likelihood)
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

test_foreach_modeling_same_SNP_3_IR_10_states_50_500k_all_IR_consistent <-
  foreach(i=c(10:13), .packages = c("tidyverse")) %:%
  foreach(j=1:8, .packages = c("tidyverse")) %dopar% {
    optim(
      state_9_same_SNP_3_IR_seg_rob_initial_par[[j]],
      state_9_same_SNP_IR_para_modu_v1_consistent,
      data = den_table_ISNP_IR_state_Ler_sep_seg[[i]],
      data2 = den_table_f_all_RCO[[i]],
      control = list(fnscale=-1, maxit=800, REPORT=1, trace=6),
      method = c("L-BFGS-B"),
      #upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.9999, 0.9999, 0.9999, 0.9999, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
      lower = rep(10^-8, 15))
  }

#among replicates of the optimization result, their recombination rates of each state create different predicted genetic length.
#Thus we need another rescaling process to make these predicted recombination rate comparable. 
#rescaling optimized CO rate of 10 states of 4 bins
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
    final_op_para_modified_9_state_3IR_SNP_consistent_first_trial[[i]][[j]] <- rescaling_state_9_seg_same_SNP_3_IR_modu_v2_test_rec(test_foreach_modeling_same_SNP_3_IR_10_states_50_500k_all_IR_consistent[[i]][[j]]$par, den_table_ISNP_IR_state_Ler_sep_seg[[i+9]], den_table_f_all_RCO[[i+9]])
  }
}

#the function to produce table with rescaled recombination rate and other parameters
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

final_op_para_modified_same_3IR_SNP_50k_500k_10_states_consistent_rep1 <- final_op_para_modified_9_state_3IR_SNP_consistent_first_trial %>%
  map(. %>% transform_para_tibble_10_states()) %>%
  bind_rows() %>%
  mutate(bin = rep(bin_size[10:13], each = 15))


write_delim(final_op_para_modified_same_3IR_SNP_50k_500k_10_states_consistent_rep1, "./analysis_output/final_op_para_modified_same_3IR_SNP_50k_500k_10_states_consistent_rep1_ll", delim = "\t", col_names = TRUE)

final_op_para_modified_same_3IR_SNP_50k_500k_10_states_consistent_rep1_l <- final_op_para_modified_same_3IR_SNP_50k_500k_10_states_consistent_rep1 %>%
  split(.$bin)

#create the modified function that can show rescaled CO rate (10 states/same SNP and IR effect, including state9/SV)
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

#extract the best fit of each bin for the following analysis
table_para_ll <- vector("list", length = length(final_op_para_modified_9_state_3IR_SNP_consistent_first_trial))

for (i in seq_along(table_para_ll)) {
  for (j in seq_along(final_op_para_modified_9_state_3IR_SNP_consistent_first_trial[[1]])) {
    x <- state_9_same_SNP_IR_para_modu_v1_consistent(final_op_para_modified_9_state_3IR_SNP_consistent_first_trial[[i]][[j]], den_table_ISNP_IR_state_Ler_sep_seg[[i+9]], data2 = den_table_f_all_RCO[[i+9]])
    test_table <- den_table_f_all_RCO[[i+9]] %>%
      mutate(predicted_Recrate = state_9_same_SNP_IR_para_modu_v2_rec(final_op_para_modified_9_state_3IR_SNP_consistent_first_trial[[i]][[j]], den_table_ISNP_IR_state_Ler_sep_seg[[i+9]], den_table_f_all_RCO[[i+9]]))
    
    R_square <- 1-sum((test_table$predicted_Recrate-test_table$Recrate_Rowan)^2)/sum((test_table$Recrate_Rowan-mean(test_table$Recrate_Rowan))^2)
    
    table_para_ll[[i]][[j]] <- tibble(
      bin_size = bin_size[i+9],
      rep = j,
      para = final_op_para_modified_9_state_3IR_SNP_consistent_first_trial[[i]][[j]],
      ll = rep(x, length(final_op_para_modified_9_state_3IR_SNP_consistent_first_trial[[i]][[j]])), 
      R_square = R_square
    ) 
  }
}

#only keep the best fit with largest likelihood and make the table as a list
table_para_ll_max <- bind_rows(table_para_ll) %>%
  group_by(bin_size) %>%
  filter(ll == max(ll)) %>%
  split(.$bin_size)

#####Fig 5: make new landscape plot based on 100-kb bins
#create the table for later using the middle of arms to produce the inset of landscape
Ara_peri_posi <- tibble(
  Chr = str_c("Chr", 1:5),
  str_peri = c(11420001, 910001, 10390001, 1070001, 8890001),
  end_peri = c(18270000, 7320000, 16730000, 6630000, 15550000)
)

Fig5_tag <- Ara_genome_bed %>% 
  left_join(Ara_peri_posi) %>%
  mutate(mid_arm = (end+end_peri)/2000000)

#create the list of main and inset figures (Chr separated) 
Fig5_main_v2 <- vector("list", length = length(Fig5_tag$Chr))
Fig5_inset_v2 <- vector("list", length = length(Fig5_tag$Chr))
Fig5_f_v2 <- vector("list", length = length(Fig5_tag$Chr))


#https://stackoverflow.com/questions/5219671/it-is-possible-to-create-inset-graphs
Fig5_tag$Chr[[1]]
for (i in seq_along(Fig5_tag$Chr)) {
  
  data_Fig5 <- den_table_f_all_RCO[[11]] %>%
    mutate(predicted_Recrate = state_9_same_SNP_IR_para_modu_v2_rec(table_para_ll_max[[2]]$para, den_table_ISNP_IR_state_Ler_sep_seg[[11]], den_table_f_all_RCO[[11]])) %>%
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

path_Fig5_v2 <- str_c("Fig5_Chr", seq(1,5), "_new_model_v3.jpeg")

pwalk(list(path_Fig5_v2, Fig5_f_v2), ggsave, path = "./analysis_output/", width = 320, height = 192, units = c("mm"), dpi = 320)

Fig5_f_v2_gra_abs_test <- Fig5_main_v2[[1]] +
  theme(legend.position = c(0.8, 0.85))

ggsave("./analysis_output/Fig5_Chr1_new_model_gra_abs_v1.jpeg", Fig5_f_v2_gra_abs_test, width = 150, height = 100, units = c("mm"), dpi = 320)


#select the best optimization result for producing Fig S7
#one may have different best result from the optimization
SP_table_op_para_combined_v2 <- tibble(
  para = c(final_op_para_modified_same_3IR_SNP_50k_500k_10_states_consistent_rep1$para[1:15], "R_square"),
  para_50k = c(table_para_ll_max[[1]]$para, table_para_ll_max[[1]]$R_square[[1]]),
  para_100k = c(table_para_ll_max[[2]]$para, table_para_ll_max[[2]]$R_square[[1]]),
  para_200k = c(table_para_ll_max[[3]]$para, table_para_ll_max[[3]]$R_square[[1]]),
  para_500k = c(table_para_ll_max[[4]]$para, table_para_ll_max[[4]]$R_square[[1]])
) %>%
  mutate(para_50k = if_else(para_50k >= 0.001, round(para_50k, 3), signif(para_50k, 3))) %>%
  mutate(para_100k = if_else(para_100k >= 0.001, round(para_100k, 3), signif(para_100k, 3))) %>%
  mutate(para_200k = if_else(para_200k >= 0.001, round(para_200k, 3), signif(para_200k, 3))) %>%
  mutate(para_500k = if_else(para_500k >= 0.001, round(para_500k, 3), signif(para_500k, 3))) 

signif

write_csv(SP_table_op_para_combined_v2, "./analysis_output/SP_table_op_para_combined_v2_ll.csv", col_names = TRUE)

####produce the scatter plot for figure S7 (state9/SV SNP/IR effect added)####
den_table_f_scatter_plot_raw_rec = vector("list", length = length(c(1:4)))

#using the best case from the optimization
optimized_par_50_500k <- list(table_para_ll_max[[1]]$para, table_para_ll_max[[2]]$para, table_para_ll_max[[3]]$para, table_para_ll_max[[4]]$para)

#create the table for plots
for (i in c(1:4)) {
  den_table_f_scatter_plot_raw_rec[[i]] <- den_table_f_all_RCO[[i+9]] %>%
    left_join(Ara_peri_posi, by = "Chr") %>%
    mutate(status = if_else(end < str_peri | str > end_peri, "arms", if_else(str > str_peri & end < end_peri, "pericentromeric_regions", "overlapped"))) %>%
    mutate(predicted_Recrate =state_9_same_SNP_IR_para_modu_v2_rec(optimized_par_50_500k[[i]], den_table_ISNP_IR_state_Ler_sep_seg[[i+9]], den_table_f_all_RCO[[i+9]]))
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

ggsave("./analysis_output/Fig_SX_new_ll.png", Fig_SX_v2, width = 360, height = 280, units = c("mm"))


#perform op from Chr1 to Chr5 separately, and predict other chromosomes using 100-kb bins
#perform optimization (50k ~ 500k)
registerDoParallel(cores = 8)
getDoParWorkers()

#separate 5 Chrs
den_table_ISNP_IR_state_Ler_sep_seg_Chr_100k <- den_table_ISNP_IR_state_Ler_sep_seg[[11]] %>%
  split(.$Chr)

den_table_f_all_RCO_Chr_100k <- den_table_f_all_RCO[[11]] %>%
  split(.$Chr)


test_foreach_modeling_same_SNP_3_IR_10_states_10_500k_all_IR_consistent_Chr_100k <-
  foreach(i=c(1:5), .packages = c("tidyverse")) %:%
  foreach(j=1:8, .packages = c("tidyverse")) %dopar% {
    optim(
      state_9_same_SNP_3_IR_seg_rob_initial_par[[j]],
      state_9_same_SNP_IR_para_modu_v1_consistent,
      data = den_table_ISNP_IR_state_Ler_sep_seg_Chr_100k[[i]],
      data2 = den_table_f_all_RCO_Chr_100k[[i]],
      control = list(fnscale=-1, maxit=800, REPORT=1, trace=6),
      method = c("L-BFGS-B"),
      #upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.9999, 0.9999, 0.9999, 0.9999, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
      lower = rep(10^-8, 15))
  }


#create the function to extract which replicate has the best result
ll_gather_best_f <- function(data){
  x = length(data)
  y = length(data[[1]])
  
  ll_list <- vector("list", length = x)
  for (i in c(1:x)) {
    for (j in c(1:y)) {
      ll_list[[i]] <- append(ll_list[[i]], data[[i]][[j]]$value)
    }
  }
  
  tibble(Chr = rep(str_c("Chr", c(1:x)), each = y),
         rep = rep(c(1:y), x),
         ll = unlist(ll_list)) %>%
    group_by(Chr) %>%
    filter(ll == max(ll))
}

best_ll_Chr <- ll_gather_best_f(test_foreach_modeling_same_SNP_3_IR_10_states_10_500k_all_IR_consistent_Chr_100k)
test_foreach_modeling_same_SNP_3_IR_10_states_10_500k_all_IR_consistent_Chr_100k[[1]][[5]]$par

#extract the information of parameters from each replicate
op_para_Chr_list_raw <- vector("list", length = length(best_ll_Chr$Chr))

for (i in seq_along(op_para_Chr_list)) {
  for (j in seq_along(test_foreach_modeling_same_SNP_3_IR_10_states_10_500k_all_IR_consistent_Chr_100k[[1]])) {
    op_para_Chr_list_raw[[i]][[j]] <- tibble(
      Chr = rep(best_ll_Chr$Chr[[i]], 15),
      rep = rep(j, 15),
      para = test_foreach_modeling_same_SNP_3_IR_10_states_10_500k_all_IR_consistent_Chr_100k[[i]][[j]]$par
    ) 
  }
}

#keep the best fit of each chromosome for the following analysis
op_para_Chr_list <- bind_rows(op_para_Chr_list_raw) %>%
  left_join(best_ll_Chr) %>%
  drop_na() %>%
  mutate(Chr_light = as.double(str_sub(Chr, 4, 4))) %>%
  split(.$Chr_light)


#calculae R square of the prediction of the best fit
R_square_Chr_list <- vector("list", length = length(op_para_Chr_list))

for (i in seq_along(R_square_Chr_list)) {
  for (j in seq_along(R_square_Chr_list)) {
    predicted_rec <- state_9_same_SNP_IR_para_modu_v2_rec(op_para_Chr_list[[i]]$para, den_table_ISNP_IR_state_Ler_sep_seg_Chr_100k[[j]], den_table_f_all_RCO_Chr_100k[[j]])
    R_square_Chr_list[[i]][[j]] <- 1 - sum((den_table_f_all_RCO_Chr_100k[[j]]$Recrate_Rowan - predicted_rec)^2)/sum((den_table_f_all_RCO_Chr_100k[[j]]$Recrate_Rowan - mean(den_table_f_all_RCO_Chr_100k[[j]]$Recrate_Rowan))^2)
  }
}


R_square_based_on_op_Chr_table <- tibble(
  Chr = str_c("Chr", c(1:5), "_predict"),
  Chr1_fit = round(R_square_Chr_list[[1]], 3),
  Chr2_fit = round(R_square_Chr_list[[2]], 3),
  Chr3_fit = round(R_square_Chr_list[[3]], 3),
  Chr4_fit = round(R_square_Chr_list[[4]], 3),
  Chr5_fit = round(R_square_Chr_list[[5]], 3)
)

write_delim(R_square_based_on_op_Chr_table, "./analysis_output/R_square_based_on_op_Chr_table_v3", delim = "\t", col_names = TRUE)


