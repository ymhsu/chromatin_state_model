Packages <- c("scales", "tidyverse", "ggrepel", "ggsci", "ggpubr", "doMC", "doParallel", "foreach", "slider")
lapply(Packages, library, character.only = TRUE)


#change the current directory as the working directory
setwd("/data/projects/thesis/INRA_project/Ara_TE_task/R_markdown/Model_1st/chromatin_state_model/")

#Using all intergenic regions between protein coding genes to perform the analysis

#load all IR 
Intergenic_region_transcription_bed <- read_delim("./data/Fig3/TAIR10_protein_coding_genes_IR_bed_trimmed_f", delim = "\t", col_names = c("Chr", "str", "end", "feature", "group_label", "transcription"))

#load IR with intersected Rowan's CO interval
Intergenic_region_transcription_RCO_bed <- read_delim("./data/Fig3/TAIR10_protein_coding_genes_IR_bed_trimmed_RCO", delim = "\t", col_names = c("Chr", "str", "end", "feature", "group_label", "transcription", "Chr_CO", "str_CO", "end_CO", "sel_420", "CO_l"))

#calculate sum of CO for each IR that intersects Rowan's COs
RCO_sum_transcription <- Intergenic_region_transcription_RCO_bed %>%
  mutate(CO_n = (end - str)/(end_CO-str_CO)) %>%
  group_by(transcription, group_label) %>%
  summarise(sum_CO = sum(CO_n))

#based on two data above, calculate CO rate of IRs and divide IRs into groups based on every 500 bps
Intergenic_region_transcription_bed_RCO_group_500bp <- Intergenic_region_transcription_bed %>%
  left_join(RCO_sum_transcription) %>%
  select(Chr, str, end, feature, transcription, group_label, sum_CO) %>%
  replace_na(list(sum_CO = 0)) %>%
  mutate(size = end - str) %>%
  mutate(group_size = if_else(size %% 500 == 0, size %/% 500, size %/% 500 + 1)) %>%
  mutate(group_size = if_else(group_size > 20, 21, group_size)) %>%
  arrange(group_size) %>%
  split(.$group_size)
  
#create the function to generate error bars of each IR group
jackknife_function_IR <- function(data1){
  length_group_size = length(data1$group_size)
  values_repeats = rep("NA", length_group_size)
  for(i in 1:length_group_size) {
    # remove i'th item
    y <- data1[-i,]
    values_repeats[i] <- sum(y$sum_CO)/sum(y$size)/2182/2*10^8
  }
  mean_all_Rec = sum(data1$sum_CO)/sum(data1$size)/2182/2*10^8
  values_repeats <- as.double(values_repeats)
  values_repeats_f <- length_group_size*mean_all_Rec - (length_group_size - 1)*values_repeats
  stderr <- sqrt(var(values_repeats_f)/length_group_size)
  tibble(
    lower_conf = round(mean_all_Rec - 1.96*stderr, 2),
    mean = round(mean_all_Rec, 2), 
    upper_conf = round(mean_all_Rec + 1.96*stderr, 2),
    group_size = data1$group_size[[1]]
  )
}

#produce confidence intervals of different IR groups
conf_Recrate_IR_500bp <- Intergenic_region_transcription_bed_RCO_group_500bp %>%
  map(~ jackknife_function_IR(data = .)) %>%
  bind_rows()

#produce confidence intervals of different IR groups based on 3 transcription types
Intergenic_region_transcription_bed_RCO_group_500bp_transcription <- Intergenic_region_transcription_bed_RCO_group_500bp %>%
  map(. %>% split(.$transcription))

conf_Recrate_IR_500bp_transcription_f <- Intergenic_region_transcription_bed_RCO_group_500bp_transcription %>%
  map(. %>% map (. %>% jackknife_function_IR())) %>%
  map(. %>% bind_rows()) %>%
  bind_rows() %>%
  mutate(transcription = rep(c("convergent", "divergent", "parallel"), 21))




#####the calculation of predicted CO rate using 10 states with or without modulation of IR size#####

###create the table of predicted CO rate (which are blue curves later) with genome-wide CO rate of 10 states only###
#based on the fraction of 10 states, we can calculate the theoretical recrate of each group of IR by using observed genome-wide recrate of 10 state
#calculate genome-wide recrate of 10 states
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

#sum of CO in 9 states
state9_CO_sum <- chromatin_state_total_R_CO_raw_noSV %>%
  mutate(CO_n = (end-str)/(end_CO-str_CO)) %>%
  group_by(state) %>%
  summarise(sum_CO = sum(CO_n)) 
sum(state9_CO_sum$sum_CO)
#import CO in SVs
SV_raw_R_CO <- read_delim("./data/Fig2/SV_raw_R_CO", delim = "\t", col_names = c("Chr", "str", "end", "state", "Chr_CO","str_CO", "end_CO", "Sel_420", "CO_l"))

#sum of CO in SV
SV_CO_sum <- SV_raw_R_CO %>%
  mutate(CO_n = (end-str)/(end_CO-str_CO)) %>%
  group_by(state) %>%
  summarise(sum_CO = sum(CO_n)) 

#the table of the genome-wide recrate of 10 states
experimental_rec_state <- state_bp_sum %>%
  left_join(state9_CO_sum) %>%
  replace_na(list(sum_CO = 0)) %>%
  mutate(sum_CO = if_else(state == "SV", 17077 - sum(sum_CO), sum_CO)) %>%
  mutate(experimental_rec = sum_CO/sum_bp/2182/2*10^8) %>%
  select(state, experimental_rec)

#load the IR information with state
Intergenic_region_transcription_bed_state <- read_delim("./data/Fig3/TAIR10_protein_coding_genes_IR_bed_trimmed_state", delim = "\t", col_names = c("Chr", "str", "end", "feature", "group_label", "transcription", "state"))

#create the list of each IR with 10 states which can connect with the sum of bp of states of each IR later
IR_size_group_state_list <- tibble(
  group_size = rep(c(1:21), each = 30),
  transcription = rep(c("convergent", "divergent", "parallel"), each = 10, times = 21),
  state = rep(c(str_c("state", c(1:9)), "SV"), 63)
) 

#calculate the sum of base pairs of states for each IR
IR_size_group_state_bp_sum <- Intergenic_region_transcription_bed %>%
  mutate(size = end - str) %>%
  mutate(group_size = if_else(size %% 500 == 0, size %/% 500, size %/% 500 + 1)) %>%
  mutate(group_size = if_else(group_size > 20, 21, group_size)) %>%
  select(group_label, group_size) %>%
  left_join(Intergenic_region_transcription_bed_state) %>%
  group_by(transcription, group_size, state) %>%
  summarise(sum_state_bp = sum(end-str))

#join the table of summed state data in the backbone list table
IR_size_group_state_prediction_raw <- IR_size_group_state_list %>%
  left_join(IR_size_group_state_bp_sum) %>%
  replace_na(list(sum_state_bp = 0)) %>%
  group_by(group_size, transcription) %>%
  mutate(sum_bp_bin = sum(sum_state_bp)) %>%
  mutate(frac_state_bin = sum_state_bp/sum_bp_bin) %>%
  left_join(experimental_rec_state)
 
##produce the table with predicted CO rate for each group of IR (all event & 3 transcription types)##

#3 types
IR_size_group_state_prediction_df_types <- IR_size_group_state_prediction_raw %>%
  mutate(pre_CO_rate_part = experimental_rec*frac_state_bin) %>%
  group_by(transcription, group_size) %>%
  summarise(pre_CO_rate = sum(pre_CO_rate_part))

#all types merges
IR_size_group_state_prediction_all <- IR_size_group_state_prediction_raw %>%
  group_by(group_size, state, experimental_rec) %>%
  #filter(group_size == 1 & state == "state1")
  summarise(sum_state_bp = sum(sum_state_bp)) %>%
  group_by(group_size) %>%
  mutate(sum_bp_bin = sum(sum_state_bp)) %>%
  mutate(frac_state_bin = sum_state_bp/sum_bp_bin) %>%
  mutate(pre_CO_rate_part = experimental_rec*frac_state_bin) %>%
  summarise(pre_CO_rate = sum(pre_CO_rate_part))


###create the table of predicted CO rate (which are red curves later) with genome-wide CO rate of 10 states and the modulation of IR size###

#produce a table with IR size, predicted CO rate before the modulation of IR size, and the experimental CO rate, and the sum of experimental COs

#produce the table with summed information of states for each IR
IR_state_sum_event <- Intergenic_region_transcription_bed_state %>%
  group_by(group_label, state) %>%
  summarise(sum_state_bp = sum(end-str))

#all the IR label for the following table production
IR_group_label <- Intergenic_region_transcription_bed$group_label

Intergenic_region_transcription_bed %>%
  group_by(transcription) %>%
  summarise(count = n())

#create the table with predicted and experimental CO rate and IR size for the following optimization
IR_pre_CO_rate_event <- tibble(
  group_label = rep(IR_group_label, each = 10),
  state = rep(c(str_c("state", c(1:9)), "SV"), times = length(IR_group_label))
) %>%
  left_join(IR_state_sum_event) %>%
  replace_na(list(sum_state_bp = 0)) %>%
  left_join(experimental_rec_state) %>%
  mutate(pre_CO_rate_part = experimental_rec*sum_state_bp) %>%
  group_by(group_label) %>%
  summarise(size = sum(sum_state_bp), pre_CO_rate = sum(pre_CO_rate_part)/sum(sum_state_bp))  %>%
  left_join(RCO_sum_transcription) %>%
  select(-transcription) %>%
  replace_na(list(sum_CO = 0)) %>%
  mutate(CO_rate = sum_CO/size/2182/2/10^-8)


#create the function of optimization based on IR size
IR_effect_all_events <-  function(a, data) {
  x = data %>%
    mutate(IR_para1 = a[1], IR_para2 =a[2], IR_para3 = a[3]) %>%
    mutate(value = 1/(IR_para1+IR_para2*exp(-IR_para3*size/1000))) %>%
    mutate(pre_CO_rate_modu = pre_CO_rate*value) 
  
  diff = x$pre_CO_rate_modu - x$CO_rate
  sum(diff^2)
}

IR_effect_all_events(rep(1,3), IR_pre_CO_rate_event)

IR_effect_all_event_optim_sigmoid <- optim(
  rep(1,3),
  IR_effect_all_events,
  data = IR_pre_CO_rate_event,
  method = c("L-BFGS-B"),
  lower = rep(10^-8, 3)
)

#After optimization, we add the optimized parameters to create table with the predicted CO rate modified by the modulation of IR size

#create IR label with transcription 
IR_label_transcription <- Intergenic_region_transcription_bed %>%
  select(group_label, transcription)

#produce the table before calculating the predicted CO rate with the modulation of IR size for all events merged and 3 types of IRs 
IR_pre_CO_rate_event_raw <- IR_pre_CO_rate_event %>%
  mutate(IR_para1 = IR_effect_all_event_optim_sigmoid$par[[1]], IR_para2 =IR_effect_all_event_optim_sigmoid$par[[2]], IR_para3 = IR_effect_all_event_optim_sigmoid$par[[3]]) %>%
  mutate(value = 1/(IR_para1+IR_para2*exp(-IR_para3*size/1000))) %>%
  mutate(pre_CO_rate_modu = pre_CO_rate*value) %>%
  left_join(IR_label_transcription)  %>%
  mutate(group_size = if_else(size %% 500 == 0, size %/% 500, size %/% 500 + 1)) %>%
  mutate(group_size = if_else(group_size > 20, 21, group_size))

#create the table with predicted CO rate with the IR-size modulation (3 transcription types separated)
IR_pre_CO_rate_event_modu_df_types <- IR_pre_CO_rate_event_raw %>%
  group_by(transcription, group_size) %>%
  summarise(pre_CO_rate_modu = sum(pre_CO_rate_modu*size)/sum(size))

#create the table with predicted CO rate with the IR-size modulation (3 transcription types merged)  
IR_pre_CO_rate_event_modu_all <- IR_pre_CO_rate_event_raw %>%
  group_by(group_size) %>%
  summarise(pre_CO_rate_modu = sum(pre_CO_rate_modu*size)/sum(size))


###draw plots with experimental (histograms) and predicted CO rate (red/blue curves) ####
Fig3_plot_function <- function(data){
  data %>%
    ggplot() +
    geom_col(aes(x = bp/1000, Recrate), fill = "grey50", colour = "black", size=0.7, width = 0.4) +
    geom_line(aes(x = bp/1000, pre_CO_rate), color = "blue", size = 0.8) +
    geom_line(aes(x = bp/1000, pre_CO_rate_modu), color = "red", size = 0.8) +
    geom_errorbar(aes(x = bp/1000, ymin=mean, ymax=upper_conf), width=.2, position = position_dodge(width = 0.9)) +
    labs(y = "recombination rate (cM/Mb)") +
    facet_wrap(~transcription) +
    theme_bw() +
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
          legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), 
          axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18)) +
    scale_x_continuous(name = "kb", breaks = c(seq(0, 10, 0.5)))
}


Fig3_top_f <- bind_rows(Intergenic_region_transcription_bed_RCO_group_500bp) %>%
  group_by(group_size) %>%
  summarise(sum_CO = sum(sum_CO), sum_size = sum(size)) %>%
  mutate(Recrate = sum_CO/sum_size/2/2182*10^8) %>%
  mutate(bp_str = (group_size-1)*500, bp_end = (group_size)*500, bp = (bp_str+bp_end)/2) %>%
  left_join(conf_Recrate_IR_500bp, by = c("group_size")) %>%
  mutate(group_size_l = c(str_c(seq(0.5, 10, 0.5), " kb"), "> 10 kb")) %>%
  mutate(group_size_l = factor(group_size_l, levels = c(str_c(seq(0.5, 10, 0.5), " kb"), "> 10 kb")), transcription = "all events") %>%
  mutate(bp = if_else(group_size == 21, 10500, bp)) %>%
  left_join(IR_size_group_state_prediction_all) %>%
  left_join(IR_pre_CO_rate_event_modu_all) %>%
  Fig3_plot_function() + theme(axis.text.x = element_text(size = 16, angle = 15, hjust = 0.5, vjust = 1, face = "bold"))


Fig3_bottom_f <- bind_rows(Intergenic_region_transcription_bed_RCO_group_500bp) %>%
  group_by(transcription, group_size) %>%
  summarise(sum_CO = sum(sum_CO), sum_size = sum(size)) %>%
  mutate(Recrate = sum_CO/sum_size/2/2182*10^8) %>%
  mutate(bp_str = (group_size-1)*500, bp_end = (group_size)*500, bp = (bp_str+bp_end)/2) %>%
  left_join(conf_Recrate_IR_500bp_transcription_f, by = c("transcription", "group_size")) %>%
  mutate(group_size_l = c(str_c(seq(0.5, 10, 0.5), " kb"), "> 10 kb")) %>%
  mutate(group_size_l = factor(group_size_l, levels = c(str_c(seq(0.5, 10, 0.5), " kb"), "> 10 kb"))) %>%
  mutate(bp = if_else(group_size == 21, 10500, bp)) %>%
  left_join(IR_size_group_state_prediction_df_types) %>%
  left_join(IR_pre_CO_rate_event_modu_df_types) %>%
  Fig3_plot_function() + theme(axis.text.x = element_text(size = 8, angle = 15, hjust = 0.5, vjust = 1, face = "bold"))


Fig3_final_f <- ggarrange(Fig3_top_f, Fig3_bottom_f, nrow = 2)

ggsave("./analysis/Fig3_final_f.jpeg", Fig3_final_f, width = 330, height = 288, units = c("mm"), dpi = 320)
550*0.6
480*0.6


#produce different stages of IR fig 3 top to show how we improve the model
Fig3_top_data_plot <- bind_rows(Intergenic_region_transcription_bed_RCO_group_500bp) %>%
  group_by(group_size) %>%
  summarise(sum_CO = sum(sum_CO), sum_size = sum(size)) %>%
  mutate(Recrate = sum_CO/sum_size/2/2182*10^8) %>%
  mutate(bp_str = (group_size-1)*500, bp_end = (group_size)*500, bp = (bp_str+bp_end)/2) %>%
  left_join(conf_Recrate_IR_500bp, by = c("group_size")) %>%
  mutate(group_size_l = c(str_c(seq(0.5, 10, 0.5), " kb"), "> 10 kb")) %>%
  mutate(group_size_l = factor(group_size_l, levels = c(str_c(seq(0.5, 10, 0.5), " kb"), "> 10 kb")), transcription = "all events") %>%
  mutate(bp = if_else(group_size == 21, 10500, bp)) %>%
  left_join(IR_size_group_state_prediction_all) %>%
  left_join(IR_pre_CO_rate_event_modu_all) 

Fig3_top_f_raw <- Fig3_top_data_plot %>%
  ggplot() +
  geom_col(aes(x = bp/1000, Recrate), fill = "grey50", colour = "black", size=0.7, width = 0.4) +
  #geom_line(aes(x = bp/1000, pre_CO_rate), color = "blue", size = 0.8) +
  #geom_line(aes(x = bp/1000, pre_CO_rate_modu), color = "red", size = 0.8) +
  geom_errorbar(aes(x = bp/1000, ymin=mean, ymax=upper_conf), width=.2, position = position_dodge(width = 0.9)) +
  labs(y = "recombination rate (cM/Mb)") +
  facet_wrap(~transcription) +
  theme_bw() +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), 
        axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18)) +
  scale_x_continuous(name = "kb", breaks = c(seq(0, 10, 0.5))) + theme(axis.text.x = element_text(size = 16, angle = 15, hjust = 0.5, vjust = 1, face = "bold"))

Fig3_top_f_10states_fit <- Fig3_top_data_plot %>%
  ggplot() +
  geom_col(aes(x = bp/1000, Recrate), fill = "grey50", colour = "black", size=0.7, width = 0.4) +
  geom_line(aes(x = bp/1000, pre_CO_rate), color = "blue", size = 0.8) +
  #geom_line(aes(x = bp/1000, pre_CO_rate_modu), color = "red", size = 0.8) +
  geom_errorbar(aes(x = bp/1000, ymin=mean, ymax=upper_conf), width=.2, position = position_dodge(width = 0.9)) +
  labs(y = "recombination rate (cM/Mb)") +
  facet_wrap(~transcription) +
  theme_bw() +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), 
        axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18)) +
  scale_x_continuous(name = "kb", breaks = c(seq(0, 10, 0.5))) + theme(axis.text.x = element_text(size = 16, angle = 15, hjust = 0.5, vjust = 1, face = "bold"))

Fig3_top_f_IR_fit <- Fig3_top_data_plot %>%
  ggplot() +
  geom_col(aes(x = bp/1000, Recrate), fill = "grey50", colour = "black", size=0.7, width = 0.4) +
  geom_line(aes(x = bp/1000, pre_CO_rate), color = "blue", size = 0.8) +
  geom_line(aes(x = bp/1000, pre_CO_rate_modu), color = "red", size = 0.8) +
  geom_errorbar(aes(x = bp/1000, ymin=mean, ymax=upper_conf), width=.2, position = position_dodge(width = 0.9)) +
  labs(y = "recombination rate (cM/Mb)") +
  facet_wrap(~transcription) +
  theme_bw() +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), 
        axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18)) +
  scale_x_continuous(name = "kb", breaks = c(seq(0, 10, 0.5))) + theme(axis.text.x = element_text(size = 16, angle = 15, hjust = 0.5, vjust = 1, face = "bold"))


ggsave("./analysis/Fig3_top_f_raw.jpeg", Fig3_top_f_raw, width = 300, height = 180, units = c("mm"), dpi = 320)
ggsave("./analysis/Fig3_top_f_10states_fit.jpeg", Fig3_top_f_10states_fit, width = 300, height = 180, units = c("mm"), dpi = 320)
ggsave("./analysis/Fig3_top_f_IR_fit.jpeg", Fig3_top_f_IR_fit, width = 300, height = 180, units = c("mm"), dpi = 320)
ggsave("./analysis/Fig3_bottom_f_presentation.jpeg", Fig3_bottom_f, width = 300, height = 180, units = c("mm"), dpi = 320)


Fig3_bottom_f

#######################old plot without any predicted CO rate ##############################
Fig3_top <- bind_rows(Intergenic_region_transcription_bed_RCO_group_500bp) %>%
  group_by(group_size) %>%
  summarise(sum_CO = sum(sum_CO), sum_size = sum(size)) %>%
  mutate(Recrate = sum_CO/sum_size/2/2182*10^8) %>%
  mutate(bp_str = (group_size-1)*500, bp_end = (group_size)*500, bp = (bp_str+bp_end)/2) %>%
  left_join(conf_Recrate_IR_500bp, by = c("group_size")) %>%
  mutate(group_size_l = c(str_c(seq(0.5, 10, 0.5), " kb"), "> 10 kb")) %>%
  mutate(group_size_l = factor(group_size_l, levels = c(str_c(seq(0.5, 10, 0.5), " kb"), "> 10 kb")), transcription = "all events") %>%
  mutate(bp = if_else(group_size == 21, 10500, bp)) %>%
  ggplot() +
  geom_col(aes(x = bp/1000, Recrate), fill = "grey50", colour = "black", size=0.7) +
  geom_errorbar(aes(x = bp/1000, ymin=mean, ymax=upper_conf), width=.2, position = position_dodge(width = 0.9)) +
  labs(y = "recombination rate (cM/Mb)") +
  facet_wrap(~transcription) +
  theme_bw() +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), 
        axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 16, angle = 15, hjust = 0.5, vjust = 1, face = "bold")) +
  scale_x_continuous(name = "kb", breaks = c(seq(0, 10, 0.5)))


Fig3_bottom <- bind_rows(Intergenic_region_transcription_bed_RCO_group_500bp) %>%
  group_by(transcription, group_size) %>%
  summarise(sum_CO = sum(sum_CO), sum_size = sum(size)) %>%
  mutate(Recrate = sum_CO/sum_size/2/2182*10^8) %>%
  mutate(bp_str = (group_size-1)*500, bp_end = (group_size)*500, bp = (bp_str+bp_end)/2) %>%
  left_join(conf_Recrate_IR_500bp_transcription_f, by = c("transcription", "group_size")) %>%
  mutate(group_size_l = c(str_c(seq(0.5, 10, 0.5), " kb"), "> 10 kb")) %>%
  mutate(group_size_l = factor(group_size_l, levels = c(str_c(seq(0.5, 10, 0.5), " kb"), "> 10 kb"))) %>%
  mutate(bp = if_else(group_size == 21, 10500, bp)) %>%
  ggplot() +
  geom_col(aes(x = bp/1000, Recrate), fill = "grey50", colour = "black", size=0.7) +
  geom_errorbar(aes(x = bp/1000, ymin=mean, ymax=upper_conf), width=.2, position = position_dodge(width = 0.9)) +
  labs(y = "recombination rate (cM/Mb)") +
  facet_wrap(~transcription, nrow = 1) + 
  theme_bw() +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), panel.spacing.x = unit(0,"line"), 
        axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 12, angle = 15, hjust = 0.5, vjust = 1, face = "bold")) +
  scale_x_continuous(name = "kb", breaks = c(seq(0, 10, 0.5)))


Fig3_final <- ggarrange(Fig3_top, Fig3_bottom, nrow = 2)

ggsave("./analysis/Fig3_top_final.jpeg", Fig3_top, width = 300, height = 200, units = c("mm"), dpi = 320)
ggsave("./analysis/Fig3_bottom_final.jpeg", Fig3_bottom, width = 450, height = 300, units = c("mm"), dpi = 320)
ggsave("./analysis/Fig3_final.jpeg", Fig3_final, width = 550, height = 480, units = c("mm"), dpi = 320)

#######################old plot without any predicted CO rate ##############################

Intergenic_region_transcription_bed_syn_f_state_10 <- read_delim("./data/Fig3/Intergenic_region_transcription_bed_syn_f_state_10", delim = "\t", col_names = c("Chr", "str", "end", "type", "transcription", "group_label", "Chr_state", "str_state", "end_state", "state", "state_label")) %>%
  select(-Chr_state, -str_state, -end_state, -state_label)

Intergenic_region_transcription_bed_syn_f_state_10_sum <- read_delim("./data/Fig3/Intergenic_region_transcription_bed_syn_f_state_10_RCO", delim = "\t", col_names = c("Chr", "str", "end", "type", "transcription", "group_label", "Chr_state", "str_state", "end_state", "state", "state_label", "Chr_CO", "str_CO", "end_CO", "sel", "CO_label")) %>%
  mutate(CO_n = (end-str)/(end_CO-str_CO)) %>%
  group_by(transcription, group_label, state) %>%
  summarise(sum_CO = sum(CO_n))


sum(Intergenic_region_transcription_bed_syn_f_state_10_sum$sum_CO)

Intergenic_region_transcription_bed_RCO_group_500bp_list <- bind_rows(Intergenic_region_transcription_bed_RCO_group_500bp) %>%
  select(transcription, group_label, group_size)


library(ggrepel)
frac_state_IR <- Intergenic_region_transcription_bed_syn_f_state_10 %>%
  group_by(state) %>%
  summarise(sum_size = sum(end-str)) %>%
  mutate(total = sum(sum_size)) %>%
  mutate(fraction = sum_size/total*100) %>%
  mutate(percent_fraction = percent(fraction/100, accuracy = 0.2)) %>%
  mutate(cul_frac = rev(cumsum(rev(fraction)))) %>%
  mutate(half_frac = fraction/2) %>%
  mutate(posi_y = cul_frac-fraction+half_frac) %>%
  mutate(type = "The_fraction_in_intergenic_regions") %>%
  select(-cul_frac, -half_frac) %>%
  mutate(label_pie = str_c(state, ": ", percent_fraction)) %>%
  ggplot(aes(x="", y=fraction, fill=state))+
  geom_bar(width = 1, stat = "identity", show.legend = FALSE) + 
  coord_polar("y", start=0) +
  theme(axis.text.x=element_blank()) +
  geom_label_repel(aes(label = label_pie, y = posi_y), fontface = "bold", size=6, show.legend = F, min.segment.length = 0.1, seed = 42, box.padding = 0.5, nudge_x = 0, nudge_y = 0) +
  scale_fill_manual(values=c(pal_npg("nrc", alpha = 0.7)(10))) +
  facet_wrap(~type, nrow = 1) +
  labs(x = NULL, y = NULL, fill = NULL) + theme_bw() + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), strip.text.x = element_text(colour = "black", face = "bold", size = 20), legend.title = element_blank(), panel.spacing.x = unit(0,"line"), plot.margin = margin(0,0,0,0, "cm"))


CO_rate_state <- Intergenic_region_transcription_bed_syn_f_state_10 %>%
  group_by(transcription, group_label, state) %>%
  mutate(size = end - str) %>%
  summarise(sum_size = sum(size)) %>%
  left_join(Intergenic_region_transcription_bed_syn_f_state_10_sum) %>%
  replace_na(list(sum_CO = 0)) %>%
  #group_by(group_size) %>%
  #group_by(group_size, state) %>%
  #group_by(group_size) %>%
  group_by(transcription, group_label, state, sum_size) %>%
  summarise(sum_CO = sum(sum_CO)) %>%
  ungroup() %>%
  left_join(Intergenic_region_transcription_bed_RCO_group_500bp_list) %>%
  mutate(state_n = if_else(state == "state8" | state == "state9", state, if_else(state == "state2" | state == "state4"| state == "state5", "state245", "state1367"))) %>%
  group_by(group_size, state_n) %>%
  summarise(sum_size = sum(sum_size), sum_CO = sum(sum_CO)) %>%
  mutate(Recrate = sum_CO/sum_size/2182/2*10^8) %>%
  mutate(x = if_else(group_size == 21, group_size/2, group_size/2-0.25)) %>%
  #View()
  ggplot() +
  geom_col(aes(x, Recrate)) +
  #geom_col(aes(group_size, sum_size)) +
  labs(y = "recombination rate (cM/Mb)") +
  facet_wrap(~state_n, scales = "free_y") + 
  theme_bw() +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), panel.spacing.x = unit(0,"line"), 
        axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 12, hjust = 0.5, vjust = 1, face = "bold")) +
  scale_x_continuous(name = "kb", breaks = c(seq(0, 10, 0.5)))
  #facet_wrap(scales = "free_y")

IR_state_analysis <- ggarrange(frac_state_IR, CO_rate_state, widths = c(1,1.5), heights = c(2,0.8)) 

IR_state_analysis_f <- IR_state_analysis + bgcolor("white")

ggsave("./analysis/IR_state_analysis_test.jpeg", IR_state_analysis_f, width = 540, height = 360, units = c("mm"), dpi = 320)
?ggarrange
