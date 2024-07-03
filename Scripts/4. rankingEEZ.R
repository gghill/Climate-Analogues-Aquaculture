## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
##      
##                  Ranking the best and worst EEZs 
##
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------

# best/worst EEZ ranking by proportion for: 
# a) combined all  
# b) per group 

closeAllConnections()
rm()
gc(reset=TRUE)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))    # sets directory to script path
getwd()

#### SETUP ####

# load packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(gt)
library(stringr)
library(sf)
library(knitr)
library(kableExtra)
library(gridExtra)
library(ggrepel)
library(patchwork)
library(purrr)


path_out <- "../Results/EEZs/"

dpi = 300

# cleaning up MegaDF
MegaDF_EEZall <- MegaDF %>%
  select(-x, -y) %>%
  drop_na(EEZ)

#### EEZs with No Dissimilarity ####

NoDissEEZs <- data.frame(EEZ = character(0), Average = numeric(0),
                         Sp_count = numeric(0), Species = character(0),
                         Scenario = character(0))

scenarios <- unique(MegaDF_EEZall$Scenario)

for (scenario in scenarios) {
  
  # average dissim per eez
  eez_av <- MegaDF_EEZall %>%
    filter(Scenario == scenario) %>%
    select(-Scenario) %>%
    group_by(EEZ) %>%
    summarise(Average = round(mean(sigmaNovelty), 2),
              Sp_count = n_distinct(Species),
              Species = paste(unique(Species), collapse = ","),
              .groups = "drop") 
  
  # eezs with no dissim
  eez_no <- eez_av %>%
    filter(Average == 0) %>%
    mutate(Scenario = scenario) %>%
    arrange(desc(Scenario), desc(Sp_count))
  
  # add to NoDissEEZs data frame
  NoDissEEZs <- rbind(NoDissEEZs, eez_no)

}

NoDissEEZs$Species <- str_replace_all(NoDissEEZs$Species, 
                                             ",", ", ")
a <- dplyr::slice(NoDissEEZs, 1:10)
b <- dplyr::slice(NoDissEEZs, 105:111)
c <- bind_rows(a, b)
write.csv(NoDissEEZs, "../Results/EEZs/NoDissEEZs.csv")


#### Ranking ALL by EEZ loop ####

scenarios <- unique(MegaDF_EEZall$Scenario)
scenario_data_frames <- list()

for (scenario in scenarios) {
  
  # mean dissim, count & spp per eez for all sigma novelty
  eez_scenario <- MegaDF_EEZall %>%
    filter(Scenario == scenario) %>%
    select(-Scenario, -sigmaNovelty) %>%
    group_by(EEZ) %>%
    summarise(
      Total_count = n_distinct(Species),
      .groups = "drop"
    )
  
  # count & spp per eez for dissim >2
  eez2 <- MegaDF_EEZall %>%
    filter(Scenario == scenario) %>%
    filter(sigmaNovelty >= 2) %>%
    group_by(EEZ, Species) %>%
    summarise(
      Average = round(mean(sigmaNovelty), 2),
      .groups = "drop"
    ) %>%
    group_by(EEZ) %>%
    summarise(
      Count = n_distinct(Species),
      Average_Scenario = round(mean(Average), 2),
      .groups = "drop"
    )
  
  # combine all and filtered data frames
  combined <- full_join(eez_scenario, eez2, by = "EEZ") %>%
    mutate(Proportion = round(if_else(Total_count == 0, 0, Count / Total_count), 2)) %>%
    select(EEZ, Average_Scenario, Proportion)
  
  # create scenario-specific column names
  combined <- combined %>%
    rename(
      !!paste0("Average_", scenario) := Average_Scenario,
      !!paste0("Proportion_", scenario) := Proportion
    )
  
  scenario_data_frames[[scenario]] <- combined

}

# Combine all scenario data frames
combined_all_scenarios <- reduce(scenario_data_frames, full_join, by = "EEZ") 

# Optionally, merge with GDP data and handle NA values
if (exists("EEZ_value")) {
  combined_all_scenarios <- full_join(combined_all_scenarios, EEZ_value, by = "EEZ")  %>%
  relocate(Value, .after = EEZ)
}

# Export the combined data frame
write.csv(combined_all_scenarios, paste0(path_out, "EEZrankingALL_combined.csv"))

# Only EEZs with value
EEZswValue <- combined_all_scenarios %>%
  filter(!is.na(Value))

write.csv(EEZswValue, paste0(path_out, "EEZswValue_rank.csv"))

#### Ranking Groups by EEZ loop ####

groups <- unique(MegaDF_EEZall$Group)

MegaRank_group <- data.frame(EEZ = character(0), Average = numeric(0), Proportion = numeric(0), Group = character (0), Scenario = character(0))

# Outer loop
for (group in groups) {
  
  group_rank <- MegaDF_EEZall %>%
    filter(Group == group) 

# Inner loop only for all groups combined
scenarios <- unique(group_rank$Scenario)

for (scenario in scenarios) {
  
  # mean dissim, count & spp per eez for all sigma novelty
  eez_scenario <- group_rank %>%
    filter(Scenario == scenario) %>%
    group_by(EEZ) %>%
    select(-Scenario, -Group, -sigmaNovelty) %>%
    #group_by(EEZ, Species) %>% arrange(EEZ) %>%
    #summarise(Average = round(mean(sigmaNovelty), 3), .groups = "drop") %>%
    group_by(EEZ) %>% 
    summarise(Total_count = n_distinct(Species), 
             .groups = "drop") 
  
  # count & spp per eez for dissim >2
  eez2 <- group_rank %>%
    filter(Scenario == scenario) %>%
    select(-Scenario, -Group) %>%
    filter(sigmaNovelty >= 2) %>%
    group_by(EEZ, Species) %>% arrange(EEZ) %>%
    summarise(Average = round(mean(sigmaNovelty), 2), .groups = "drop") %>%
    group_by(EEZ) %>%
    summarise(Count = n_distinct(Species), Average = round(mean(Average), 2),
              .groups = "drop") 
  
  # combine all and filtered data frames
  combined <- full_join(eez_scenario, eez2, by = "EEZ")
  
  # get the proportions (global prop - divide count by 327, regional prop - divide by Total_count)
  combined <- combined %>% 
    mutate(Proportion = Count/Total_count) %>%
    select(-Count, -Total_count)
  
  combined$Proportion <- round(combined$Proportion, 2)
  
  combined <- combined %>%
    mutate_if(is.numeric, funs(ifelse(is.na(.), 0, .))) %>%
    filter(Average > 0)
  
  combined <- combined[!duplicated(combined), ]
    
  combined <- combined %>%
    arrange(EEZ) %>%
    mutate(Group = group, Scenario = scenario)
  
  MegaRank_group <- rbind(MegaRank_group, combined)

}

}

# save, name based on proportion
write.csv(MegaRank_group, "../Results/MegaRank_group.csv")

#### NEW EEZ ranking and value figure ####

# loop to rank EEZ by mean dissim, prop per scenario
scenarios <- unique(MegaRank_group$Scenario)

for(scenario in scenarios) {
  
  # rank by dissim and prop per scenario
  scen_rank <- MegaRank_group %>%
    filter(Scenario == scenario) %>%
    arrange(desc(Average), desc(Proportion)) 
  
  # add GDP data per EEZ
  scen_rank <- left_join(scen_rank, EEZ_value, by = "EEZ")
  
  # add 1 to each value column so no zeros exist, only NA values
  scen_rank <- scen_rank %>%
    mutate(PreLog = Value+1) %>%
    drop_na(Scenario) %>%
    na.omit
  
  # log the value column and replace NA values with - number
  scen_rank <- scen_rank %>% mutate(logPC = log10(PreLog))
  
  df_name <- paste0("Rank_", scenario)
  assign(df_name, scen_rank)
  
  write.csv(scen_rank, paste0(path_out, scenario, "_rank.csv"))
  
  # slice top 10 rows and remove from plotting df
  rank_slice <- scen_rank %>%
    #mutate(Weight = (logPC^2)+(Average^2)) %>% ORIGINAL
    mutate(Weight = scale(logPC^2) + scale(Average)) %>% #normalized weighting
    arrange(desc(Weight)) %>%
    slice(1:20)
  
  df_name1 <- paste0("Rank_slice", scenario)
  assign(df_name1, rank_slice)
  
  unique_EEZs <- rank_slice %>% distinct(EEZ, .keep_all = TRUE)
  
  rank_plot <- anti_join(scen_rank, 
                          rank_slice, by = c("EEZ", "Average", "Proportion", "Group", "Scenario"))
  
  # plot 
  plot <- ggplot() +
    geom_point(data = rank_plot, aes(x = Average, y = logPC), colour = "grey", size = 1, shape = 16) +
    geom_point(data = rank_slice, aes(x = Average, y = logPC, 
                                          colour = Group, size = Proportion),
                                      alpha = 0.9, shape = 21) +
    #scale_colour_manual(values = group_colors) +
    geom_text_repel(data = unique_EEZs, aes(x = Average, y = logPC, label = EEZ),
                    max.overlaps = 30, size = 3) +
    labs(x = "Mean Dissimilarity per EEZ", y = "Aquaculture Value 2021 (USD) (log scale)") +
    scale_fill_brewer(palette = "Dark2") +
    scale_y_continuous(limits = c(0, 11)) +  
    scale_x_continuous(limits = c(2, 8.3)) + 
    scale_size_continuous(breaks=c(0, 0.25, 0.50, 0.75, 1), labels=c(0, 0.25, 0.50, 0.75, 1), 
                          limits = c(0, 1), range = c(2, 8)) +
    theme_light() +
    theme(legend.position="bottom") +
    ggtitle(paste("Scenario:", scenario)) +
    coord_fixed(ratio = 0.8) +
    geom_vline(xintercept = 4, 
               linetype = "dashed", color = "darkgrey", size = 0.5) +
    geom_hline(yintercept = 4, 
               linetype = "dashed", color = "darkgrey", size = 0.5) +
    guides(colour = guide_legend(override.aes = list(shape = 16))) 
  
  plot
  #plotslist[[as.character(scenario)]] <- plot
  
  plot_name <- paste0(scenario, "Value_diss_20.svg")
  ggsave(path = "../Results", plot_name, plot, 
         height=8,width=8,dpi=dpi, device = 'svg')
  
}

# points per quadrant
#119
q1Dev119 <- Rank_ssp119 %>%
  filter(Average<4, logPC<6) 
n_distinct(q1Dev119$EEZ) #3

q2Opt119 <- Rank_ssp119 %>%
  filter(Average<4, logPC>6) 
n_distinct(q2Opt119$EEZ) #19

# 370 
q1Dev370 <- Rank_ssp370 %>%
  filter(Average<4, logPC<6) 
n_distinct(q1Dev370$EEZ) #43

q2Opt370 <- Rank_ssp370 %>%
  filter(Average<4, logPC>6) 
n_distinct(q2Opt370$EEZ) #75

q3Mit370 <- Rank_ssp370 %>%
  filter(Average>4, logPC>6) 
n_distinct(q3Mit370$EEZ) #28

q4Ad370 <- Rank_ssp370 %>%
  filter(Average>4, logPC<6) 
n_distinct(q4Ad370$EEZ) #16

#585

q1Dev585 <- Rank_ssp585 %>%
  filter(Average<4, logPC<6) 
n_distinct(q1Dev585$EEZ) #39

q2Opt585 <- Rank_ssp585 %>%
  filter(Average<4, logPC>6) 
n_distinct(q2Opt585$EEZ) #73

q3Mit585 <- Rank_ssp585 %>%
  filter(Average>4, logPC>6) 
n_distinct(q3Mit585$EEZ) #51

q4Ad585 <- Rank_ssp585 %>%
  filter(Average>4, logPC<6) 
n_distinct(q4Ad585$EEZ) #34




#### EEZ TABLE APLPHABETICAL ####

# table based on new quadrant figure, including all EEZs and aquaculture value
ValuexDissim <- bind_rows(ssp119_rank, ssp370_rank)
ValuexDissim <- bind_rows(ValuexDissim, ssp585_rank)

ValuexDissim <- ValuexDissim %>%
  select(-Proportion, -PreLog) %>%
  mutate(Metric = (logPC^3)+(Average^3)) %>%
  unite("Group", Group:Scenario, remove = TRUE) %>%
  mutate(across(Metric, round, 0))

ValuexDissim_table <- tidyr::pivot_wider(ValuexDissim, id_cols = c("EEZ"),
                                         names_from = c("Group"), 
                                         values_from = "Metric")



ValuexDissim_table <- dplyr::full_join(ValuexDissim_table, EEZ_GDP_data, by = "EEZ")

write.csv(ValuexDissim_table, paste0(path_out, "Full_ValuexDissim_metrictable.csv"))

# OLD CODE
# ranking by sum of proportion
rankTable <- tidyr::pivot_wider(MegaRank_all, id_cols = c("EEZ"),
                                     names_from = c("Group", "Scenario"), 
                                     values_from = "Proportion")


sum_row <- rowSums(rankTable[, -1], na.rm = TRUE)
rankTable$pSum <- sum_row

rankTable <- rankTable %>%
  select(EEZ, pSum)


MegaRank_table <- MegaRank_all
MegaRank_all$Proportion <- round(MegaRank_all$Proportion, 2)
MegaRank_all$Average <- round(MegaRank_all$Average, 2)
MegaRank_table$Proportion <- paste0(MegaRank_table$Proportion, " (", MegaRank_table$Average, ")")
MegaRank_table <- MegaRank_table %>% select(-Average)

MegaRank_table <- tidyr::pivot_wider(MegaRank_table, id_cols = c("EEZ"),
                                  names_from = c("Group", "Scenario"), 
                                  values_from = "Proportion")

MegaRank_table <- full_join(MegaRank_table, rankTable, by = "EEZ")

arrange(MegaRank_table, desc(pSum))
unique(MegaRank_table$EEZ)
  

write.csv(MegaRank_table, paste0(path_out, "MegaRank_tableNEW.csv"))

