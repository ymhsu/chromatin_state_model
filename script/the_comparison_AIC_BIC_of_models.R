#change the directory "chromatin_state_model" as the working directory
setwd("/data/projects/thesis/INRA_project/Ara_TE_task/R_markdown/Model_1st/chromatin_state_model/")

#using "p_load" from the package "pacman" to install and load necessary packages
install.packages("pacman")
library(pacman)

Packages <- c("scales", "tidyverse", "ggrepel", "ggsci", "ggpubr", "doMC", "doParallel", "foreach", "slider", "cowplot", "combinat")
p_load(Packages, character.only = TRUE)
#lapply(Packages, library, character.only = TRUE)

#For validating the tradeoff between the fit and complexity of our models, we compared AIC and BIC of four stages of our model.
#These four stages are 10 states (non rescaled), 10 states with the IR effect (non rescaled), 10 states with IR/SNP effects (non rescaled) and 10 states with IR/SNP effects (rescaled).

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
paths_den_table_RCO <-
  str_c(
    "./data/Fig1/",
    "den_table_",
    bin_name,
    "_RCO_raw_mid_bed"
  )

#produce the file of all bins (50-500 k) with the sum of CO (Rowan)
#calculate the sum of CO
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

#the list of different bins with recombination rate based on CO intervals of Rowan
den_table_f_all_RCO = vector("list", length(bin_name))

for (i in seq_along(bin_name)) {
  den_table_f_all_RCO[[i]] <- den_table_list[[i]] %>%
    left_join(den_table_CO_sum[[i]]) %>%
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

######modeling######
###model 1: 10 states without rescaling
#perform the optimization based on only 10 states (log likelihood) 
#the function of producing log likelihood
op_f_stage1 <-  function(a, data, data2) {
  para_d <- tibble(
    feature = c(str_c("state", c(1:9)), "SV"),
    raw_rec = c(a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10])
  )
  
  x = data %>%
    left_join(para_d, by = c("feature")) %>%
    mutate(raw_recrate_v2 = raw_rec*size/size_bin) %>%
    group_by(chr_bin) %>%
    summarise(recrate_non_rescaling = sum(raw_recrate_v2)) %>%
    mutate(Chr = str_sub(chr_bin, 1, 4), label = str_sub(chr_bin, 6, 20)) %>%
    mutate(label = as.double(label)) %>%
    arrange(Chr, label)
  
  x2 <- data2 %>%
    left_join(x, by = c("Chr", "chr_bin")) %>%
    drop_na() %>%
    mutate(p_CO = if_else(recrate_non_rescaling*physical_length/100 < 1, recrate_non_rescaling*physical_length/100, 1-10^-10), p_noCO = 1-p_CO) %>%
    mutate(log_likelihood = sum_CO*log(p_CO) + (2182*2-sum_CO)*log(p_noCO))
  
  sum(x2$log_likelihood)
}

#run optimizations
#create initial parameters
op_stage1_initial_par <- vector(mode = "list", length = 8)

for (i in c(1:8)) {
  op_stage1_initial_par[[i]] <- c(1 * (0.5 + runif(10)))
}

#perform optimization (50k ~ 500k)
registerDoParallel(cores = 8)
getDoParWorkers()

op_stage1_result_50_500k <-
  foreach(i=c(10:13), .packages = c("tidyverse")) %:%
  foreach(j=1:8, .packages = c("tidyverse")) %dopar% {
    optim(
      op_stage1_initial_par[[j]],
      op_f_stage1,
      data = den_table_ISNP_IR_state_Ler_sep_seg[[i]],
      data2 = den_table_f_all_RCO[[i]],
      control = list(fnscale=-1, maxit=800, REPORT=1, trace=6),
      method = c("L-BFGS-B"),
      #upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.9999, 0.9999, 0.9999, 0.9999, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
      lower = rep(10^-8, 10))
  }

###model 2: 10 states with the IR effect (non rescaled)
#the function of producing log likelihood 
op_f_stage2 <-  function(a, data, data2) {
  para_d <- tibble(
    feature = c(str_c("state", c(1:9)), "SV"),
    raw_rec = c(a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10]),
    IR_para1 = rep(a[11], 10),
    IR_para2 = rep(a[12], 10),
    IR_para3 = rep(a[13], 10)
  )
  
  x = data %>%
    left_join(para_d, by = c("feature")) %>%
    mutate(raw_recrate_v2 = if_else(size_IR == 0, raw_rec*size/size_bin, raw_rec/(IR_para1 + IR_para2*exp(-IR_para3*size_IR/1000))*size/size_bin)) %>%
    group_by(chr_bin) %>%
    summarise(recrate_non_rescaling = sum(raw_recrate_v2)) %>%
    mutate(Chr = str_sub(chr_bin, 1, 4), label = str_sub(chr_bin, 6, 20)) %>%
    mutate(label = as.double(label)) %>%
    arrange(Chr, label)
  
  x2 <- data2 %>%
    left_join(x, by = c("Chr", "chr_bin")) %>%
    drop_na() %>%
    mutate(p_CO = if_else(recrate_non_rescaling*physical_length/100 < 1, recrate_non_rescaling*physical_length/100, 1-10^-10), p_noCO = 1-p_CO) %>%
    mutate(log_likelihood = sum_CO*log(p_CO) + (2182*2-sum_CO)*log(p_noCO))   
  
  sum(x2$log_likelihood)
}

#run optimizations
#create initial parameters
op_stage2_initial_par <- vector(mode = "list", length = 8)

for (i in c(1:8)) {
  op_stage2_initial_par[[i]] <- c(1 * (0.5 + runif(13)))
}

#perform optimization (50k ~ 500k)
registerDoParallel(cores = 12)
getDoParWorkers()

op_stage2_result_50_500k <-
  foreach(i=c(10:13), .packages = c("tidyverse")) %:%
  foreach(j=1:8, .packages = c("tidyverse")) %dopar% {
    optim(
      op_stage2_initial_par[[j]],
      op_f_stage2,
      data = den_table_ISNP_IR_state_Ler_sep_seg[[i]],
      data2 = den_table_f_all_RCO[[i]],
      control = list(fnscale=-1, maxit=800, REPORT=1, trace=6),
      method = c("L-BFGS-B"),
      #upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.9999, 0.9999, 0.9999, 0.9999, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
      lower = rep(10^-8, 13))
  }


###model 3: 10 states with the IR/SNP effect (non rescaled)
#the function of producing log likelihood 
op_f_stage3 <-  function(a, data, data2) {
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
    summarise(recrate_non_rescaling = sum(raw_recrate_v2)) %>%
    mutate(Chr = str_sub(chr_bin, 1, 4), label = str_sub(chr_bin, 6, 20)) %>%
    mutate(label = as.double(label)) %>%
    arrange(Chr, label)
  
  x2 <- data2 %>%
    left_join(x, by = c("Chr", "chr_bin")) %>%
    drop_na() %>%
    mutate(p_CO = if_else(recrate_non_rescaling*physical_length/100 < 1, recrate_non_rescaling*physical_length/100, 1-10^-10), p_noCO = 1-p_CO) %>%
    mutate(log_likelihood = sum_CO*log(p_CO) + (2182*2-sum_CO)*log(p_noCO)) 
  
  sum(x2$log_likelihood)
}

#run optimizations
#create initial parameters
op_stage3_initial_par <- vector(mode = "list", length = 8)

for (i in c(1:8)) {
  op_stage3_initial_par[[i]] <- c(1 * (0.5 + runif(15)))
}

#perform optimization (50k ~ 500k)
registerDoParallel(cores = 8)
getDoParWorkers()

op_stage3_result_50_500k <-
  foreach(i=c(10:13), .packages = c("tidyverse")) %:%
  foreach(j=1:8, .packages = c("tidyverse")) %dopar% {
    optim(
      op_stage3_initial_par[[j]],
      op_f_stage3,
      data = den_table_ISNP_IR_state_Ler_sep_seg[[i]],
      data2 = den_table_f_all_RCO[[i]],
      control = list(fnscale=-1, maxit=800, REPORT=1, trace=6),
      method = c("L-BFGS-B"),
      #upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.9999, 0.9999, 0.9999, 0.9999, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
      lower = rep(10^-8, 15))
  }

op_stage2_result_50_500k[[4]]
optim(
  op_stage3_initial_par[[5]],
  op_f_stage3,
  data = den_table_ISNP_IR_state_Ler_sep_seg[[13]],
  data2 = den_table_f_all_RCO[[13]],
  control = list(fnscale=-1, maxit=800, REPORT=1, trace=6),
  method = c("L-BFGS-B"),
  #upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.9999, 0.9999, 0.9999, 0.9999, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
  lower = rep(10^-8, 15))
op_f_stage3(c(2.87667, 1.88043, 2.51343, 2.30418, 2.03941, 2.32579, 2.48509, 0.887931, 1.50436, 2.83146, 3.73607, 1e-08, 1e-08, 0.210961, 1.63419), den_table_ISNP_IR_state_Ler_sep_seg[[13]], den_table_f_all_RCO[[13]]) %>%
  View()

###model 4: 10 states with the IR/SNP effect (rescaled)
#the function of producing log likelihood 
op_f_stage4 <-  function(a, data, data2) {
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
    mutate(rescaled_model_rec_rates = recrate_non_rescalling*genetic_length/pre_genetic_length) %>%
    mutate(p_CO = if_else(rescaled_model_rec_rates*physical_length/100 < 1, rescaled_model_rec_rates*physical_length/100, 1-10^-10), p_noCO = 1-p_CO) %>%
    mutate(log_likelihood = sum_CO*log(p_CO) + (2182*2-sum_CO)*log(p_noCO))
  
  sum(x2$log_likelihood)
}

#run optimizations
#create initial parameters
op_stage4_initial_par <- vector(mode = "list", length = 8)

for (i in c(1:8)) {
  op_stage4_initial_par[[i]] <- c(1 * (0.5 + runif(15)))
}

#perform optimization (50k ~ 500k)
registerDoParallel(cores = 8)
getDoParWorkers()

op_stage4_result_50_500k <-
  foreach(i=c(10:13), .packages = c("tidyverse")) %:%
  foreach(j=1:8, .packages = c("tidyverse")) %dopar% {
    optim(
      op_stage4_initial_par[[j]],
      op_f_stage4,
      data = den_table_ISNP_IR_state_Ler_sep_seg[[i]],
      data2 = den_table_f_all_RCO[[i]],
      control = list(fnscale=-1, maxit=800, REPORT=1, trace=6),
      method = c("L-BFGS-B"),
      #upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 0.9999, 0.9999, 0.9999, 0.9999, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
      lower = rep(10^-8, 15))
  }

#create the function for calculating explained variance of 4 models
#10 states
op_f_var_stage1 <-  function(a, data, data2) {
  para_d <- tibble(
    feature = c(str_c("state", c(1:9)), "SV"),
    raw_rec = c(a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10])
  )
  
  x = data %>%
    left_join(para_d, by = c("feature")) %>%
    mutate(raw_recrate_v2 = raw_rec*size/size_bin) %>%
    group_by(chr_bin) %>%
    summarise(recrate_non_rescaling = sum(raw_recrate_v2)) %>%
    mutate(Chr = str_sub(chr_bin, 1, 4), label = str_sub(chr_bin, 6, 20)) %>%
    mutate(label = as.double(label)) %>%
    arrange(Chr, label)
  
  x2 <- data2 %>%
    left_join(x, by = c("Chr", "chr_bin")) %>%
    drop_na() 
  
  ex_rec <- x2$Recrate_Rowan
  pred_rec <- x2$recrate_non_rescaling
  explained_var_vector = round(1 - sum((ex_rec-pred_rec)^2)/sum((mean(ex_rec)-ex_rec)^2), 2)
  explained_var_vector
}

#10 states + IR
op_f_var_stage2 <-  function(a, data, data2) {
  para_d <- tibble(
    feature = c(str_c("state", c(1:9)), "SV"),
    raw_rec = c(a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10]),
    IR_para1 = rep(a[11], 10),
    IR_para2 = rep(a[12], 10),
    IR_para3 = rep(a[13], 10)
  )
  
  x = data %>%
    left_join(para_d, by = c("feature")) %>%
    mutate(raw_recrate_v2 = if_else(size_IR == 0, raw_rec*size/size_bin, raw_rec/(IR_para1 + IR_para2*exp(-IR_para3*size_IR/1000))*size/size_bin)) %>%
    group_by(chr_bin) %>%
    summarise(recrate_non_rescaling = sum(raw_recrate_v2)) %>%
    mutate(Chr = str_sub(chr_bin, 1, 4), label = str_sub(chr_bin, 6, 20)) %>%
    mutate(label = as.double(label)) %>%
    arrange(Chr, label)
  
  x2 <- data2 %>%
    left_join(x, by = c("Chr", "chr_bin")) %>%
    drop_na()    
  
  ex_rec <- x2$Recrate_Rowan
  pred_rec <- x2$recrate_non_rescaling
  explained_var_vector = round(1 - sum((ex_rec-pred_rec)^2)/sum((mean(ex_rec)-ex_rec)^2), 2)
  explained_var_vector
}

#10 states + IR + SNP
op_f_var_stage3 <-  function(a, data, data2) {
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
    summarise(recrate_non_rescaling = sum(raw_recrate_v2)) %>%
    mutate(Chr = str_sub(chr_bin, 1, 4), label = str_sub(chr_bin, 6, 20)) %>%
    mutate(label = as.double(label)) %>%
    arrange(Chr, label)
  
  x2 <- data2 %>%
    left_join(x, by = c("Chr", "chr_bin")) %>%
    drop_na() 
  
  ex_rec <- x2$Recrate_Rowan
  pred_rec <- x2$recrate_non_rescaling
  explained_var_vector = round(1 - sum((ex_rec-pred_rec)^2)/sum((mean(ex_rec)-ex_rec)^2), 2)
  explained_var_vector
  
}

#10 states + IR + SNP + rescaling
op_f_var_stage4 <-  function(a, data, data2) {
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
  
  ex_rec <- x2$Recrate_Rowan
  pred_rec <- x2$rescaled_model_rec_rates
  explained_var_vector = round(1 - sum((ex_rec-pred_rec)^2)/sum((mean(ex_rec)-ex_rec)^2), 2)
  explained_var_vector
}

#create the function for calculating AIC/BIC of four models, and include the information of R square
ll_gather_f <- function(data, data2, data3, data4, data5){
  x = length(data)
  y = length(data[[1]])
  
  bin_number <- vector("character", length = 0)
  for (i in c(1:x)) {
    bin_number <- append(bin_number, length(data2[[i]]$Chr))
  }
  
  ll_list <- vector("list", length = x)
  for (i in c(1:x)) {
    for (j in c(1:y)) {
      ll_list[[i]] <- append(ll_list[[i]], data[[i]][[j]]$value)
    }
  }
  
  R_square_list <- vector("list", length = x)
  
  for (i in c(1:x)) {
    for (j in c(1:y)) {
      x <- data3(data[[i]][[j]]$par, data = data4[[i]], data2 = data2[[i]])
      R_square_list[[i]] <- append(R_square_list[[i]], x)
    }
  }
  
  tibble(bin_size = rep(c(50, 100, 200, 500), each = 8),
         rep = rep(c(1:8), 4),
         bin_number = as.double(rep(bin_number, each = 8)), 
         ll = unlist(ll_list),
         R_square = unlist(R_square_list)) %>%
    group_by(bin_size) %>%
    filter(ll == max(ll)) %>%
    mutate(num_p = data5, AIC = 2*num_p-2*ll, BIC = num_p*log(bin_number)-2*ll)
}


#below are the data necessary for the function "ll_gather_f"
#data2
den_table_f_all_RCO_50_500kb_l <- list(den_table_f_all_RCO[[10]], den_table_f_all_RCO[[11]], den_table_f_all_RCO[[12]], den_table_f_all_RCO[[13]])
#data3
op_f_var_list <- list(op_f_var_stage1, op_f_var_stage2, op_f_var_stage3, op_f_var_stage4)
#data4
den_table_ISNP_IR_state_Ler_sep_seg_l <- list(den_table_ISNP_IR_state_Ler_sep_seg[[10]], den_table_ISNP_IR_state_Ler_sep_seg[[11]], den_table_ISNP_IR_state_Ler_sep_seg[[12]], den_table_ISNP_IR_state_Ler_sep_seg[[13]])

#produce the final table with AIC/BIC information
sp_data_AIC_BIC_4models <- ll_gather_f(op_stage1_result_50_500k, den_table_f_all_RCO_50_500kb_l, op_f_var_list[[1]], den_table_ISNP_IR_state_Ler_sep_seg_l, 10) %>%
  bind_rows(ll_gather_f(op_stage2_result_50_500k, den_table_f_all_RCO_50_500kb_l, op_f_var_list[[2]], den_table_ISNP_IR_state_Ler_sep_seg_l, 13)) %>%
  bind_rows(ll_gather_f(op_stage3_result_50_500k, den_table_f_all_RCO_50_500kb_l, op_f_var_list[[3]], den_table_ISNP_IR_state_Ler_sep_seg_l, 15)) %>%
  bind_rows(ll_gather_f(op_stage4_result_50_500k, den_table_f_all_RCO_50_500kb_l, op_f_var_list[[4]], den_table_ISNP_IR_state_Ler_sep_seg_l, 15)) %>%
  ungroup() %>%
  mutate(stage = rep(c("10 states", "10 states + IR", "10_states + IR + SNP", "10_states + IR + SNP + rescaling"), each = 4)) %>%
  arrange(bin_size) %>%
  select(bin_size, AIC, BIC, R_square, stage) %>%
  mutate(AIC = round(AIC, 2), BIC = round(BIC, 2))

write_delim(sp_data_AIC_BIC_4models, "./analysis_output/sp_data_AIC_BIC_4models", delim = "\t", col_names = TRUE)
