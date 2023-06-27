## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
##      
##                  Creating the Mega DF and Mapping
##
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------

# MegaDF which combines all climate dissimilarity DFs (results from Jorge) for all species
# MegaDF includes 333 species analysed under 3 climate change scenarios
# MegaDF is then modified to only show points > 2

closeAllConnections()
rm(list=(ls()[ls()]))
gc(reset=TRUE)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))    # sets directory to script path
getwd()

#### SETUP ####

# load packages
library(raster)
library(sf)
library(ggplot2)
library(geodata)
library(tidyr)
library(dplyr)
library(rgdal)
library(dggridR)
library(cowplot)

# plot function for EEZs
plot_data_factor <- function (column, df) {
  if (startsWith(column, 'P')) {
    ggplot(eez_sf) +
      geom_map(
        data = world, map = world,
        aes(map_id = region),
        color = "grey", fill = "lightgray", linewidth = 0.01
      ) +
      geom_sf(aes(fill=eez_sf[[column]])) +
      labs(fill = column) + xlab("Longitude") +
      scale_fill_gradientn(colours = c("blue", "green", "yellow", "red"),
                           limits = c(0, 1), na.value = "white") +
      scale_x_continuous(labels = c("-120", "-60", "0", "60", "120")) +
      theme(plot.margin=unit(c(0, 0, 0, 0), "pt"))
  } else {
    ggplot(eez_sf) +
      geom_map(
        data = world, map = world,
        aes(map_id = region),
        color = "grey", fill = "lightgray", linewidth = 0.01
      ) +
      geom_sf(aes(fill=eez_sf[[column]])) +
      labs(fill = column) + xlab("Longitude") +
      scale_fill_gradientn(colours = c("blue", "green", "yellow", "red"),
                           limits = c(2, 8.3), na.value = "white") +
      scale_x_continuous(labels = c("-120", "-60", "0", "60", "120")) +
      theme(plot.margin=unit(c(0, 0, 0, 0), "pt"))
  }
}

# to plot all eezs 
eez_all_sf <- st_as_sf(eez)

# hex grid size
hex_grid <- dgconstruct(res = 6)

#### ALL DISSIMILARITY FOR ALL SPECIES - MEGADF ####
# run this section on the server

# setup for looping through scenarios
resultsScenario <- list.files("../Results/finalDissim/", recursive=FALSE, full.names=TRUE)
resultsScenario <- gsub("../Results/finalDissim/","",resultsScenario)
resultsScenario

#empty data frame to populate
MegaDF <- data.frame()

Sys.time()
for( h in resultsScenario) {
  
  scenario <- h
  
  #resultsFolderMain <- paste0("../Results/climateDissimilarity/", scenario)
  resultsFolderMain <- paste0("../Results/finalDissim/", scenario)
  resultsScope <- "sigmaNovelty" # sigmaNovelty sigmaDisappearance 
  resultsFolders <- list.files(resultsFolderMain, recursive=TRUE, pattern="climateDissimilarity.RData", full.names=TRUE)
  resultsNames <- list.files(resultsFolderMain, recursive=TRUE, pattern="climateDissimilarity.RData", full.names=FALSE)
  resultsNames <- gsub("/climateDissimilarity.RData","",resultsNames)
  
  
  # reads in climateDissimilarity rdata files for each species, extracts relevant info and adds to MegaDF
  for(i in 1:length(resultsFolders)) {
    
    active_species <- resultsNames[i]
    
    cat("\n")
    cat("# ---------------------------------\n")
    cat(i,"out of",length(resultsFolders),"\n")       # progress bar
    cat("# ---------------------------------\n")
    cat("\n")
    
    DissimilarityData <- load(resultsFolders[i])
    DissimilarityData <- dataStructureResult
    
    df <- dataStructureResult %>% 
      #filter(sigmaNovelty > 2) %>%
      select(-cell, -analogNovelty, -analogDisappearance, -sigmaDisappearance, -analogNovelty.x, 
             -analogNovelty.y, -analogDisappearance.x, -analogDisappearance.y, -analogNoveltyDist,
             -analogDisappearanceDist, -refugiaWithinFocalPoll, -refugia, -messDissimilarity) %>%
      mutate(Species = active_species, Scenario = scenario)
    
    MegaDF <-rbind(MegaDF, df)
    
  }
}
Sys.time()

# add in EEZ col based on coords
MegaDF_sp <- SpatialPointsDataFrame(coords = MegaDF[, c("x", "y")], data = MegaDF)
proj4string(MegaDF_sp) <- CRS("+proj=longlat +datum=WGS84")

# convert eez to same CRS as dataframe
eez <- spTransform(eez, CRS(proj4string(MegaDF_sp)))     
proj4string(MegaDF_sp)
proj4string(eez)

# this is the most time-intensive part, run on the server:
eez_index <- over(MegaDF_sp[, c('x', 'y')], eez)
MegaDF <- cbind(MegaDF, eez_index[2])
colnames(MegaDF)[ncol(MegaDF)] <- "EEZ" 

# assign hex cells to each coordinate point
MegaDF <- MegaDF[(MegaDF$x<177) & (MegaDF$x>-177),]
MegaDF$cell <- dgGEO_to_SEQNUM(hex_grid, MegaDF$x, MegaDF$y)$seqnum

# add in the taxonomic group of each species
MegaDF <- left_join(MegaDF, AquacultureSpeciesList, by = "Species") %>%
  relocate(Group, .after = Species)

# save MegaDF
write.csv(MegaDF, file="../Results/MegaDF.csv", col.names = FALSE, row.names = FALSE)


### MEAN DISSIM PER EEZ AND SCENARIO ---------- 

# clean MegaDF
MegaDF <- MegaDF %>%
  na.omit() %>%
  select(-x, -y) %>%
  mutate(Scenario = str_replace(Scenario, "/", ""))

# split into scenarios
df1 <- MegaDF %>%
  filter(Scenario == "ssp119") %>%
  select(-Scenario, -Group) %>%
  group_by(EEZ, Species) %>% arrange(EEZ) %>%
  summarise(Average = round(mean(sigmaNovelty), 3), .groups = "drop") %>%
  group_by(EEZ) %>% 
  summarise(Total_count = n_distinct(Species), Average = round(mean(Average), 3), 
            Species = paste(unique(Species), collapse = ",")
            , .groups = "drop")
count_all1 <- df1 %>%
  select(EEZ, Total_count)

df3 <- MegaDF %>%
  filter(Scenario == "ssp370") %>%
  select(-Scenario, -Group) %>%
  group_by(EEZ, Species) %>% arrange(EEZ) %>%
  summarise(Average = round(mean(sigmaNovelty), 3), .groups = "drop") %>%
  group_by(EEZ) %>% 
  summarise(Total_count = n_distinct(Species), Average = round(mean(Average), 3), 
            Species = paste(unique(Species), collapse = ",")
            , .groups = "drop")
count_all3 <- df3 %>%
  select(EEZ, Total_count)

df5 <- MegaDF %>%
  filter(Scenario == "ssp585") %>%
  select(-Scenario, -Group) %>%
  group_by(EEZ, Species) %>% arrange(EEZ) %>%
  summarise(Average = round(mean(sigmaNovelty), 3), .groups = "drop") %>%
  group_by(EEZ) %>% 
  summarise(Total_count = n_distinct(Species), Average = round(mean(Average), 3), 
            Species = paste(unique(Species), collapse = ",")
            , .groups = "drop")
count_all5 <- df5 %>%
  select(EEZ, Total_count)

# total species count per group for proportion per group (>2) maps

### TOTAL COUNT PER GROUP PER SCENARIO --------

groups <- unique(MegaDF$Group)


for (group in groups) {
  
  # Filter and clean group data from MegaDF
  groupdf <- MegaDF %>%
    filter(Group == group) %>%
    select(-Group) 
  
  scenarios <- unique(groupdf$Scenario)
  # Create a list to store the plots for each scenario
  scenario_plots <- list()
  
  # Split according to scenario
  for (scenario in scenarios) {
    # Subset the data frame based on the scenario using filter
    scenario_df <- groupdf %>%
      filter(grepl(scenario, Scenario)) %>%
      group_by(EEZ) %>%
      summarise(
        Total_count = n_distinct(Species),
      )
    
    colnames(scenario_df)[2] = paste0(group, scenario)
    
    # Set the name of the data frame as group + scenario
    df_name <- paste0(group, scenario)
    assign(df_name, scenario_df)
    
  }
}


### MAPPING DISSIM > 2 COMBINED ALL MAPS --------------

# rewrite to run as a loop

# dissimilarity >= 2
MegaDF2 <- MegaDF %>%
  filter(sigmaNovelty >= 2)

# ssp119 
DF2eez1 <- MegaDF2 %>%
  filter(Scenario == "ssp119") %>%
  select(-Group, -Scenario) %>%
  group_by(EEZ, Species) %>% arrange(EEZ) %>%
  summarise(Average = round(mean(sigmaNovelty), 3), .groups = "drop") %>%
  group_by(EEZ) %>% 
  summarise(Count = n_distinct(Species), Average = round(mean(Average), 3), 
            Species = paste(unique(Species), collapse = ",")
            , .groups = "drop")

DF2eez1 <- left_join(DF2eez1, count_all1, by = "EEZ") 
DF2eez1 <- DF2eez1 %>% mutate(Prop = Count/Total_count, Scenario = "ssp119")


# ssp370 
DF2eez3 <- MegaDF2 %>%
  filter(Scenario == "ssp370") %>%
  select(-Group, -Scenario) %>%
  group_by(EEZ, Species) %>% arrange(EEZ) %>%
  summarise(Average = round(mean(sigmaNovelty), 3), .groups = "drop") %>%
  group_by(EEZ) %>% 
  summarise(Count = n_distinct(Species), Average = round(mean(Average), 3), 
            Species = paste(unique(Species), collapse = ",")
            , .groups = "drop")

DF2eez3 <- left_join(DF2eez3, count_all3, by = "EEZ") 
DF2eez3 <- DF2eez3 %>% mutate(Prop = Count/Total_count, Scenario = "ssp370")

# ssp585 
DF2eez5 <- MegaDF2 %>%
  filter(Scenario == "ssp585") %>%
  select(-Group, -Scenario) %>%
  group_by(EEZ, Species) %>% arrange(EEZ) %>%
  summarise(Average = round(mean(sigmaNovelty), 3), .groups = "drop") %>%
  group_by(EEZ) %>% 
  summarise(Count = n_distinct(Species), Average = round(mean(Average), 3), 
            Species = paste(unique(Species), collapse = ",")
            , .groups = "drop")

DF2eez5 <- left_join(DF2eez5, count_all5, by = "EEZ") 
DF2eez5 <- DF2eez5 %>% mutate(Prop = Count/Total_count, Scenario = "ssp585")

DF2eezAll <- rbind(DF2eez1, DF2eez3)
DF2eezAll <- rbind(DF2eezAll, DF2eez5)


# setting up the spatial df per scenario

active_eezs119 <- eez[eez$EEZ %in% unique(DF2eez1$EEZ),]
active_eezs119@data <- cbind(active_eezs119@data, DF2eez1[,-1])
active_eezs_vis1 <- st_as_sf(active_eezs119)


active_eezs370 <- eez[eez$EEZ %in% unique(DF2eez3$EEZ),]
active_eezs370@data <- cbind(active_eezs370@data, DF2eez3[,-1])
active_eezs_vis3 <- st_as_sf(active_eezs370)


active_eezs585 <- eez[eez$EEZ %in% unique(DF2eez5$EEZ),]
active_eezs585@data <- cbind(active_eezs585@data, DF2eez5[,-1])
active_eezs_vis5 <- st_as_sf(active_eezs585)

# Assign the columns we want to map
factors <- c("Average", "Prop")

# function take column + df as input, ggplot(df) throughout

plot_data_factor <- function (column, df) {
  if (startsWith(column, 'P')) {
    ggplot(active_eezs_vis5) +
      geom_map(
        data = world, map = world,
        aes(map_id = region),
        color = "grey", fill = "lightgray", linewidth = 0.01
      ) +
      geom_sf(aes(fill=active_eezs_vis5[[column]])) +
      labs(fill = column) + xlab("Longitude") +
      scale_fill_gradientn(colours = c("blue", "green", "yellow", "red"),
                           limits = c(0, 1), na.value = NULL) +
      scale_x_continuous(labels = c("-120", "-60", "0", "60", "120")) +
      theme(plot.margin=unit(c(0, 0, 0, 0), "pt"))
  } else {
    ggplot(active_eezs_vis5) +
      geom_map(
        data = world, map = world,
        aes(map_id = region),
        color = "grey", fill = "lightgray", linewidth = 0.01
      ) +
      geom_sf(aes(fill=active_eezs_vis5[[column]])) +
      labs(fill = column) + xlab("Longitude") +
      scale_fill_gradientn(colours = c("blue", "green", "yellow", "red"),
                           limits = c(2, 8.3), na.value = NULL) +
      scale_x_continuous(labels = c("-120", "-60", "0", "60", "120")) +
      theme(plot.margin=unit(c(0, 0, 0, 0), "pt"))
  }
}

# loop start
myplots <- lapply(factors, plot_data_factor, df=loop_eez)
title <- ggdraw() + 
  draw_label(
    paste("Mean Dissimilarity >2 and proportion of species >2 per EEZ - SSP5-8.5"),
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(15, 0, 10, 7)
  )

grid <- plot_grid(plotlist = myplots)
save_grid <- plot_grid(title,grid, ncol = 1, axis='b',rel_heights = c(0.03,1))
# save_grid

# write to file
dpi=300
save_path <- "../Results/Figures/All_585.jpg"
ggsave(filename = save_path, plot = save_grid,
       height = 1200/dpi, width = 3000/dpi, dpi = dpi, device = "jpeg")

# loop end


### MAPPING PER GROUP Per SCENARIO ####

groups <- unique(MegaDF$Group)

for (group in groups) {
  
  # Filter and clean group data from MegaDF
  groupdf <- MegaDF %>%
    filter(Group == group) %>%
    select(-x, -y, -Group, -cell) %>%
    na.omit()
  
  scenarios <- unique(groupdf$Scenario)
  scenarios <- gsub("/", "", scenarios)
  
  # Create a list to store the plots for each scenario
  scenario_plots <- list()
  
  # Split according to scenario
  for (scenario in scenarios) {
    
    # total species count from all dissimilarity
    t_count <- groupdf %>%
      filter(grepl(scenario, Scenario)) %>%
      group_by(EEZ) %>%
      select(-sigmaNovelty, -Scenario) %>%
      summarise(Total_count = n_distinct(Species))
    
    
    # Subset the data frame based on the scenario using filter
    scenario_df <- groupdf %>%
      filter(grepl(scenario, Scenario), sigmaNovelty >= 2) %>%
      group_by(EEZ) %>%
      summarise(
        Species_count = n_distinct(Species),
        Average = round(mean(sigmaNovelty), 2),
        Species = paste(unique(Species), collapse = ","))
    
    scenario_df <- full_join(scenario_df, t_count, by = "EEZ")
    
    scenario_df <- scenario_df %>%
      group_by(EEZ) %>%
      mutate(Proportion = Species_count/Total_count) %>%
      na.omit()
    
    # Set the name of the data frame as group + scenario
    df_name <- paste0(group, scenario)
    assign(df_name, scenario_df)
    
    # mapping the data
    
    # Retrieve the data frame
    group_scenario_df <- get(df_name)
    
    # Add the code for plotting and saving the map here
    active_eez <- eez[eez$EEZ %in% unique(group_scenario_df$EEZ),]
    active_eez@data <- cbind(active_eez@data, group_scenario_df[,-1])
    eez_sf <- st_as_sf(active_eez)
    
    # assign the columns we want to map
    factors <- c("Average", "Proportion")
    
    EEZmap <- lapply(factors, plot_data_factor)
    title <- ggdraw() + 
      draw_label(
        paste0(group, scenario),
        fontface = 'bold',
        x = 0,
        hjust = 0
      ) +
      theme(
        plot.margin = margin(15, 0, 0, 7)
      )
    
    EEZgrid <- plot_grid(plotlist = EEZmap)
    save_grid <- plot_grid(title, EEZgrid, ncol = 1, axis='b', rel_heights = c(0.05, 0.95))
    # Save the plot to a file
    plot_filename <- paste0(group, scenario, "_plot.jpg")
    dpi = 300
    ggsave(path = "../Results/Figures/", plot_filename, save_grid, 
           height=1200/dpi,width=3000/dpi,dpi=dpi, device = 'jpeg') 
    
  }
}



### HEX PLOTS ---------------

hex = TRUE
hex_grid <- dgconstruct(res = 6)

# all eezs active from previous EEZ code - to plot EEZ layer filled by hex cells

for (group in groups) {
  
  # Filter and clean group data from MegaDF
  groupcell <- MegaDF %>%
    filter(Group == group) %>%
    select(-Group, -EEZ) %>%
    na.omit()
  
  scenarios <- unique(groupdf$Scenario)
  scenarios <- gsub("/", "", scenarios)
  # Create a list to store the plots for each scenario
  scenario_plots <- list()
  
  # Split according to scenario
  for (scenario in scenarios) {
    
    # Subset the data frame based on the scenario using filter
    scenariocell <- groupcell %>%
      filter(grepl(scenario, Scenario), sigmaNovelty >= 2) %>%
      group_by(cell) %>%
      summarise(
        Average = round(mean(sigmaNovelty), 2),) %>%
      na.omit()
    
    # Set the name of the data frame as group + scenario
    celldf_name <- paste0(group, scenario, "cell")
    assign(celldf_name, scenariocell)
    
    # mapping the data
    
    # Retrieve the data frame
    cell_df <- get(celldf_name)
    
    hex_active_plot <- dgcellstogrid(hex_grid, cell_df$cell)
    names(hex_active_plot) <- c('cell', 'geometry') # labeling with highest overlap EEZ would be nice here
    hex_active_plot <- merge(hex_active_plot, cell_df, by.x='cell', by.y='cell')
    
    # assign the columns we want to map
    factor <- c("Average")
    
    # add proportion to plot function
    hex_plot_funtion = function (column, hex=TRUE) {
      if (hex==TRUE) {
        ggplot(hex_active_plot) + 
          geom_map(
            data = world, map = world,
            aes(map_id = region),
            color = "grey", fill = "lightgray", size = 0.01
          ) + 
          geom_sf(data= eez_all_sf, aes(), fill = NA) +
          geom_sf(aes(fill=hex_active_plot[[column]]), alpha=0.8, linewidth = 0.01, color='white') + theme_minimal() +
          labs(fill = column) +
          scale_x_continuous(labels = c("-120", "-60", "0", "60", "120")) +
          scale_fill_gradientn(colours = c("blue", "green", "yellow", "red"),
                               limits = c(2, 8.3), na.value = 'grey') +
          theme(plot.margin=unit(c(0, 0, 0, 0), "pt"))
      } else {
        ggplot(eez_sf) +
          geom_map(
            data = world, map = world,
            aes(map_id = region),
            color = "grey", fill = "lightgray", linewidth = 0.01
          ) +
          geom_sf(aes(fill=eez_sf[[column]])) +
          labs(fill = column) + xlab("Longitude") +
          scale_fill_gradientn(colours = c("blue", "green", "yellow", "red"),
                               limits = c(2, scale_max), na.value = 'grey') +
          scale_x_continuous(labels = c("-120", "-60", "0", "60", "120")) +
          theme(plot.margin=unit(c(0, 0, 0, 0), "pt"))
      }
    }
    
    EEZmap <- lapply(factor, hex_plot_funtion)
    title <- ggdraw() + 
      draw_label(
        paste0(group, " ", scenario),
        fontface = 'bold',
        x = 0,
        hjust = 0
      ) +
      theme(
        plot.margin = margin(15, 0, 0, 7)
      )
    
    EEZgrid <- plot_grid(plotlist = EEZmap)
    save_grid <- plot_grid(title, EEZgrid, ncol = 1, axis='b', rel_heights = c(0.05, 0.95))
    # Save the plot to a file
    plot_filename <- paste0(group, scenario, "_hexplot.jpg")
    dpi = 300
    ggsave(path = "../Results/Figures/", plot_filename, save_grid, 
           height=1200/dpi,width=3000/dpi,dpi=dpi, device = 'jpeg') 
    
  }
}

