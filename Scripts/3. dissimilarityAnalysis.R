## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
##      
##                  Creating the MegaDF and Mapping
##
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------

# MegaDF which combines all climate dissimilarity DFs (results from Jorge Assis) for all species and scenarios
# MegaDF includes 330 species analysed under 3 climate change scenarios (3 sturgeons removed in code below)
# MegaDF is then modified to only show points > 2

closeAllConnections()
rm(all_plot)
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
library(stringr)
library(sp)
library(maps)
library(grid)
library(rgdal)


install.packages("rgdal")

world <- shapefile("../Data/spatialInformation/GSHHS_i_L1.shp")        # Downloaded from NOAA
antarctica <- shapefile("../Data/spatialInformation/GSHHS_i_L5.shp")   # Downloaded from NOAA
eez <- shapefile("../Data/spatialInformation/eez.shp")
projection <- CRS("+proj=robin +over")

# generate bounding box and clip spatal data
bb <- sf::st_union(sf::st_make_grid(
  st_bbox(c(xmin = -180,
            xmax = 180,
            ymax = 90,
            ymin = -90), crs = st_crs(4326)), n = 100))
bb <- st_transform(bb, projection)
bb_sf <- st_as_sf(bb)

world_sf <- st_as_sf(world)
world_sf <- st_transform(world_sf, projection)
world_sf <- st_buffer(world_sf, dist =0.001)
world_sf <- st_intersection(world_sf, bb_sf)
ant_sf <- st_as_sf(antarctica)
ant_sf <- st_transform(ant_sf, projection)
ant_sf <- st_buffer(ant_sf, dist =0.001)
ant_sf <- st_intersection(ant_sf, bb_sf)

# specifications for plotting
hex_grid <- dgconstruct(res = 9)
dpi =  300

my_colours <- c("#6FBBE8","#A1ECD8","#F6F9AB","#FCB46D","#B21414","#D278E4","#9914B3")
my_theme <- theme(text = element_text(family = "Times New Roman", color = "#22211d"),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.background = element_rect(fill = "#FFFFFF", color = NA),
                  panel.background = element_rect(fill = "#FFFFFF", color = NA),
                  panel.border = element_blank(),
                  legend.background = element_rect(fill = "#FFFFFF", color = NA),
                  legend.position = "bottom",
                  legend.box = "horizontal",
                  legend.key.height = unit(0.25, 'cm'),
                  legend.key.width = unit(0.75, 'cm'),
                  legend.margin = ggplot2::margin(t = -12, unit = "pt"),
                  plot.margin = unit(c(0, 0, 0, 0), "pt"))

# remove  species
MegaDF <- MegaDF %>%
  filter(!Species %in% c("Acipenser gueldenstaedtii", "Acipenser nudiventris", "Acipenser transmontanus", "Naso vlamingii", "Melicertus kerathurus", "Coryphaena hippurus"))

#### ALL DISSIMILARITY FOR ALL SPECIES -> MEGADF ####
# run this section on the server

# setup for looping through scenarios
resultsScenario <- list.files("../Data/SpeciesDissim/", 
                              recursive=FALSE, full.names=TRUE)
resultsScenario <- gsub("/DataSpeciesDissim//","",resultsScenario) # include full file path
resultsScenario

#empty data frame to populate
MegaDF <- data.frame()

Sys.time()
for( h in resultsScenario) {
  
  scenario <- h
  
  #resultsFolderMain <- paste0("/data/amm/Chapter1/SpeciesDissim/", scenario)
  resultsFolderMain <- paste0("../Data/SpeciesDissim/", scenario)
  resultsScope <- "sigmaNovelty" # sigmaNovelty sigmaDisappearance 
  resultsFolders <- list.files(resultsFolderMain, recursive=TRUE, pattern="climateDissimilarity.RData", full.names=TRUE)
  resultsNames <- list.files(resultsFolderMain, recursive=TRUE, pattern="climateDissimilarity.RData", full.names=FALSE)
  resultsNames <- gsub("/climateDissimilarity.RData","",resultsNames)
  
  
  # reads in climateDissimilarity rdata files for each species, extracts info and adds to MegaDF
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

# assign eezs to coordinates
# this is the most time-intensive part, run on the server
eez_index <- over(MegaDF_sp, eez)
MegaDF <- cbind(MegaDF, eez_index[2])
colnames(MegaDF)[ncol(MegaDF)] <- "EEZ" 

# create sf object of eezs for mapping bounadries
active_eez <- eez[eez$EEZ %in% unique(MegaDF$EEZ),]
active_eez@data <- cbind(active_eez@data, MegaDF[,-1])
eez_sf <- st_as_sf(active_eez)
eez_sf_rob <- st_transform(eez_sf, projection)

# create SpatialPointsDataFrame
MegaDF_sp <- SpatialPointsDataFrame(coords = MegaDF[, c("x", "y")], data = MegaDF)
# set the initial CRS to WGS84
proj4string(MegaDF_sp) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
# transform to Robinson Projection
MegaDF_sp <- spTransform(MegaDF_sp, CRSobj = CRS("+proj=robin +over"))
MegaDF_sf <- st_as_sf(MegaDF_sp)

# convert eez to same CRS as dataframe
eez <- spTransform(eez, CRS(proj4string(MegaDF_sp)))     
proj4string(MegaDF_sp)
proj4string(eez)

# assign hex cells to coordinates
MegaDF$cell <- dgGEO_to_SEQNUM(hex_grid, MegaDF$x, MegaDF$y)$seqnum
# to wrap cells mapped on the date line
wrapped_grid = st_wrap_dateline(MegaDF_sf, options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"), quiet = TRUE)
wrapped_grid <- st_transform(wrapped_grid, crs = "+proj=robin +over")

# add in the taxonomic group of each species
MegaDF <- left_join(MegaDF, AquacultureSpeciesList, by = "Species") %>%
  relocate(Group, .after = Species)

# save MegaDF
write.csv(MegaDF, file="../Results/MegaDF.csv", col.names = TRUE, row.names = FALSE)

### SPECIES RICHNESS MAPS ####

DF <- MegaDF[(MegaDF$x<179) & (MegaDF$x>-179),]

# filter to one scenario (makes no difference which) to remove scenario duplicates
ssp119 <- DF %>%
  filter(Scenario == "ssp119") %>%
  select(-sigmaNovelty, -x, -y) 

# total species richness per EEZ
EEZrichness <- ssp119 %>% 
  select(-Scenario, -cell) %>%
  group_by(EEZ) %>%
  summarise(Count = n_distinct(Species)) %>%
  na.omit()

# total species richness per hex cell
HEXrichness <- ssp119 %>% 
  select(-Scenario, -EEZ) %>%
  group_by(cell) %>%
  summarise(Count = n_distinct(Species)) %>%
  na.omit()


# sf object of hex cells for mapping
hex_active_plot <- dgcellstogrid(hex_grid, HEXrichness$cell)
names(hex_active_plot) <- c('cell', 'geometry') 
#hex_active_plot <- st_as_sf(hex_active_plot)
#HEXrichness <- as.data.frame(HEXrichness)
#hex_active_plot <- hex_active_plot %>%
  #left_join(HEXrichness, by = "cell")
hex_active_plot <- merge(hex_active_plot, HEXrichness, by.x='cell', by.y='cell')
#hex_active_plot <- merge(hex_active_plot, HEXrichness, by = "cell")
hex_active_plot <- st_transform(hex_active_plot, projection)
hex_active_plot <- st_buffer(hex_active_plot, dist =0.001)
hex_active_plot <- st_intersection(hex_active_plot, bb_sf)

hex_active_plot <- st_make_valid(hex_active_plot)
eez_sf <- st_make_valid(eez_sf)
hex_masked <- st_intersection(hex_active_plot, eez_sf)
hex_masked <- hex_masked[!st_is_empty(hex_masked), ]
hex_masked <- st_transform(hex_masked, crs = 4326) #transform to appropriate crs for wrapping hex cells

#combined_geometry <- do.call(rbind, lapply(list(world_sf, ant_sf, hex_masked, eez_sf), st_geometry))
#bbox <- st_bbox(combined_geometry)

hex_wrapped = st_wrap_dateline(hex_masked, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE)
hex_wrapped <- st_transform(hex_wrapped, "+proj=robin +over")


# mapping total richness per cell
Hexrich_map <- ggplot(hex_wrapped) +
  geom_sf(data = world_sf, fill = "lightgrey", alpha = 0.9) +
  geom_sf(data = ant_sf, fill = "lightgrey", alpha = 0.9) +
  geom_sf(aes(fill=Count, colours = my_colours), alpha = 0.8, linewidth = 0.03) +
  scale_fill_gradientn(guide = guide_colorbar(title="Species Richness",
                                              direction = "horizontal", 
                                              title.position = "top", 
                                              title.hjust = 0.5),
                       colours = my_colours,
                       limits = c(1, 150),
                       breaks = c(1, 50, 100, 150)) +
  labs(fill = "Species Richness") +
  my_theme

Hexrich_map

ggsave("Hexrich_map.png", plot = Hexrich_map, height=1800/dpi,width=3000/dpi,dpi=dpi, device = 'png', bg = "white")

# species richness per group
# separate "Miscellaneous" group into actual taxons
ssp119 <- ssp119 %>%
  mutate(Group = case_when(
    Species == "Pyura stolonifera" ~ "Tunicate",
    Species == "Apostichopus japonicus" ~ "Echinoderm",
    Species == "Loxechinus albus" ~ "Echinoderm",
    Species == "Paracentrotus lividus" ~ "Echinoderm",
    TRUE ~ Group
  ))

# for loop for mapping species richness by group
groups <- unique(ssp119$Group)

for (group in groups) {
  
  Group_rich <- ssp119 %>%
    filter(Group == group) %>%
    select(-Scenario, -EEZ) %>%
    group_by(cell) %>%
    summarise(Count = n_distinct(Species)) %>%
    na.omit()
  
  hex_active_plot <- dgcellstogrid(hex_grid, Group_rich$cell)
  names(hex_active_plot) <- c('cell', 'geometry') 
  hex_active_plot <- merge(hex_active_plot, Group_rich, by.x='cell', by.y='cell')
  #hex_active_plot <- merge(hex_active_plot, HEXrichness, by = "cell")
  hex_active_plot <- st_transform(hex_active_plot, projection)
  
  Hexrich_map <- ggplot(hex_active_plot) +
  geom_sf(data = world1_sf, fill = "grey", alpha = 0.9) +
  geom_sf(data = ant_sf, fill = "grey", alpha = 0.9) +
  geom_sf(aes(fill=Count), alpha = 0.8, linewidth = 0.08) +
  scale_fill_gradientn(guide = guide_colorbar(title="Species Richness",
                                            direction = "horizontal", 
                                            title.position = "top", 
                                            title.hjust = 0.5), 
                       colours = my_colours,
                       limits = c(1, 150),
                       breaks = c(1, 50, 100, 150)) +
  labs(fill = "Species Richness") +
  scale_x_continuous(breaks = c(-120, -60, 0, 60, 120),
                     labels = c("-120", "-60", "0", "60", "120")) +
  my_theme

  Grouprich_name <- paste0(group, "GroupRich.png")
  ggsave(path = "../Results/", Grouprich_name, Hexrich_map, height=1200/dpi,width=3000/dpi,dpi=dpi, device = 'png', bg = "white")
}

### MAPPING DISSIM > 2 COMBINED ALL MAPS --------------

MegaDF_allscen <- MegaDF %>%
  select(-x, -y) %>%
  na.omit() 
  
scenarios <- unique(MegaDF_allscen$Scenario)

for (scenario in scenarios) {
  
  # total species count from all dissimilarity
  t_count_hex <- MegaDF_allscen %>%
    filter(Scenario == scenario) %>%
    select(-Scenario, -sigmaNovelty) %>%
    group_by(cell) %>%
    summarise(Total_count = n_distinct(Species)) 
  
  # Subset the data frame based on the scenario using filter
  scenariocell <- MegaDF_allscen %>%
    filter(grepl(scenario, Scenario), sigmaNovelty >= 2) %>%
    group_by(cell) %>%
    summarise(
      Species_count = n_distinct(Species),
      Average = round(mean(sigmaNovelty), 2),
      Species = paste(unique(Species), collapse = ","))
  
  scenario_hexdf <- full_join(scenariocell, t_count_hex, by = "cell")
  
  scenario_hexdf <- scenario_hexdf %>%
    group_by(cell) %>%
    mutate(Proportion = Species_count/Total_count, PropT = Species_count/330) %>%
    select(-Species_count, -Total_count)
  
  scenario_hexdf$Proportion <- round(scenario_hexdf$Proportion, 3)
  scenario_hexdf$PropT <- round(scenario_hexdf$PropT, 3)
  
  # Set the name of the data frame as group + scenario
  celldf_name <- paste0(scenario, "cell")
  assign(celldf_name, scenario_hexdf)
  
  # mapping the data
  
  # Retrieve the data frame
  cell_df <- get(celldf_name)
  
  hex_active_plot <- dgcellstogrid(hex_grid, cell_df$cell)
  names(hex_active_plot) <- c('cell', 'geometry') 
  hex_active_plot <- merge(hex_active_plot, cell_df, by.x='cell', by.y='cell')
  hex_active_plot = st_wrap_dateline(hex_active_plot, 
                                  options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"), 
                                  quiet = TRUE)
  
  # Transform to Robinson Projection
  hex_active_plot <- st_transform(hex_active_plot, crs = "+proj=robin +over")
  
  Propplot <- ggplot(hex_active_plot) +
    geom_sf(data = world_sf, fill = "grey", alpha = 0.9) +
    geom_sf(data = ant_sf, fill = "grey", alpha = 0.9) +
    geom_sf(data= eez_sf_rob, aes(), fill = NA) +
    geom_sf(aes(fill= .data$Proportion), alpha = 0.8, linewidth = 0.06) +
    scale_fill_gradientn(guide = guide_colorbar(title="Proportion", 
                                                direction = "horizontal", 
                                                title.position = "top", 
                                                title.hjust = 0.5),
                         colours = my_colours,
                         limits = c(0, 1), na.value = "white") +
    my_theme
    
  Propplot_name <- paste0(scenario, "Prophex.png")
  ggsave(path = "../Results/", Propplot_name, Propplot, 
         height=1200/dpi,width=3000/dpi,dpi=dpi, device = 'png', bg = "white")

  # "global" proportion
  Propplot <- ggplot(hex_active_plot) +
    geom_sf(data = world_sf, fill = "grey", alpha = 0.9) +
    geom_sf(data = ant_sf, fill = "grey", alpha = 0.9) +
    geom_sf(data= eez_sf_rob, aes(), fill = NA) +
    geom_sf(aes(fill= .data$PropT), alpha = 0.8, linewidth = 0.06) +
    scale_fill_gradientn(guide = guide_colorbar(title="Proportion", 
                                                direction = "horizontal", 
                                                title.position = "top", 
                                                title.hjust = 0.5),
                         colours = my_colours,
                         limits = c(0, 1), na.value = "white") +
    my_theme
  
  Propplot_name <- paste0(scenario, "GlobalProphex.png")
  ggsave(path = "../Results/", Propplot_name, Propplot, 
         height=1200/dpi,width=3000/dpi,dpi=dpi, device = 'png', bg = "white")
  
  Avplot <- ggplot(hex_active_plot) +
    geom_sf(data = world_sf, fill = "grey", alpha = 0.9) +
    geom_sf(data = ant_sf, fill = "grey", alpha = 0.9) +
    geom_sf(data= eez_sf_rob, aes(), fill = NA) +
    geom_sf(aes(fill= .data$Average), alpha = 0.8, linewidth = 0.06) +
    scale_fill_gradientn(guide = guide_colorbar(title="Average", 
                                                direction = "horizontal", 
                                                title.position = "top", 
                                                title.hjust = 0.5),
                         colours = my_colours,
                         limits = c(2, 8.3), na.value = "white") +
    my_theme
  
  Avplot_name <- paste0(scenario, "Avhex.png")
  ggsave(path = "../Results/", Avplot_name, Avplot, 
         height=1200/dpi,width=3000/dpi,dpi=dpi, device = 'png', bg = "white")
  
}

# summary stats for results write-up

ssp370cell <- ssp370cell %>%
  na.omit(EEZ)

summary(ssp119cell)
1643/3358*100



### MAPPING PER GROUP Per SCENARIO ####

groups <- unique(MegaDF_allscen$Group)

for (group in groups) {
  
  # Filter and clean group data from MegaDF
  groupdf <- MegaDF %>%
    filter(Group == group)
  
  scenarios <- unique(groupdf$Scenario)

  # Split according to scenario
  for (scenario in scenarios) {
    
    # total species count from all dissimilarity
    t_count_hex <- groupdf %>%
      filter(Scenario == scenario) %>%
      select(-Scenario, -sigmaNovelty) %>%
      group_by(cell) %>%
      summarise(Total_count = n_distinct(Species)) 
    
    # Subset the data frame based on the scenario using filter
    scenariocell <- groupdf %>%
      filter(grepl(scenario, Scenario), sigmaNovelty >= 2) %>%
      group_by(cell) %>%
      summarise(
        Species_count = n_distinct(Species),
        Average = round(mean(sigmaNovelty), 2),
        Species = paste(unique(Species), collapse = ","))
    
    scenario_hexdf <- full_join(scenariocell, t_count_hex, by = "cell")
    
    scenario_hexdf <- scenario_hexdf %>%
      group_by(cell) %>%
      mutate(Proportion = Species_count/Total_count, PropT = Species_count/330) %>%
      select(-Species_count, -Total_count)
    
    scenario_hexdf$Proportion <- round(scenario_hexdf$Proportion, 3)
    scenario_hexdf$PropT <- round(scenario_hexdf$PropT, 3)
    
    # Set the name of the data frame as group + scenario
    celldf_name <- paste0(scenario, group)
    assign(celldf_name, scenario_hexdf)
    
    # mapping the data
    
    # Retrieve the data frame
    cell_df <- get(celldf_name)
    
    hex_active_plot <- dgcellstogrid(hex_grid, cell_df$cell)
    names(hex_active_plot) <- c('cell', 'geometry') 
    hex_active_plot <- merge(hex_active_plot, cell_df, by.x='cell', by.y='cell')
    hex_active_plot = st_wrap_dateline(hex_active_plot, 
                                       options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"), 
                                       quiet = TRUE)
    
    # Transform to Robinson Projection
    hex_active_plot <- st_transform(hex_active_plot, crs = "+proj=robin +over")
    
    Propplot <- ggplot(hex_active_plot) +
      geom_sf(data = world1_sf, fill = "grey", alpha = 0.9) +
      geom_sf(data = ant_sf, fill = "grey", alpha = 0.9) +
      geom_sf(data= eez_sf_rob, aes(), fill = NA) +
      geom_sf(aes(fill= .data$Proportion), alpha = 0.8, linewidth = 0.06) +
      scale_fill_gradientn(guide = guide_colorbar(title="Proportion", 
                                                  direction = "horizontal", 
                                                  title.position = "top", 
                                                  title.hjust = 0.5),
                           colours = my_colours,
                           limits = c(0, 1), na.value = "white") +
      my_theme
    
    Propplot_name <- paste0(scenario, group, "Prophex.png")
    ggsave(path = "../Results/", Propplot_name, Propplot, 
           height=1200/dpi,width=3000/dpi,dpi=dpi, device = 'png', bg = "white")
    
    
    Avplot <- ggplot() +
      geom_sf(data = world1_sf, fill = "grey", alpha = 0.9) +
      geom_sf(data = ant_sf, fill = "grey", alpha = 0.9) +
      geom_sf(data= eez_sf_rob, aes(), fill = NA) +
      geom_sf(data= hex_active_plot, aes(fill= .data$Average), alpha = 0.8, linewidth = 0.06) +
      scale_fill_gradientn(guide = guide_colorbar(title="Average", 
                                                  direction = "horizontal", 
                                                  title.position = "top", 
                                                  title.hjust = 0.5),
                           colours = my_colours,
                           limits = c(2, 8.3), na.value = "white") +
      my_theme
    
    Avplot_name <- paste0(scenario, group, "Avhex.png")
    ggsave(path = "../Results/", Avplot_name, Avplot, 
           height=1200/dpi,width=3000/dpi,dpi=dpi, device = 'png', bg = "white")
    
  }
}

