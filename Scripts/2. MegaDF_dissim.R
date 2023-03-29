## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
##      
##                  Creating the Mega DF and Mapping
##
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------

# MegaDF which combines all climate dissimilarity DFs (results from Jorge) for all species
# MegaDF selects analogNovelty >2 (shows dissimilarity) at a point level but can be modified 
# to include points <2

closeAllConnections()
rm(list=(ls()[ls()]))
gc(reset=TRUE)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))    # sets directory to script path
getwd()

#### SETUP ####

scenario <- "ssp585"                  # can incorporate into loop, to run through all scenarios
hex_grid <- dgconstruct(res = 6)      # for hex plots 

resultsFolderMain <- paste0("../Results/climateDissimilarity/", scenario)
resultsScope <- "sigmaNovelty" # sigmaNovelty sigmaDisappearance 
resultsFolders <- list.files(resultsFolderMain, recursive=TRUE, pattern=".RData", full.names=TRUE)
resultsNames <- list.files(resultsFolderMain, recursive=TRUE, pattern=".RData", full.names=FALSE)
resultsNames <- gsub("/climateDissimilarity.RData","",resultsNames)

resultsNames


### PIXEL DISSIM > 2 DF ALL SPECIES ------------

# set up empty data frame for loop to populate
BigDF <- data.frame()

# reads in climateDissimilarity rdata files for each species, extracts relevant info and adds to BigDF
for(i in 1:length(resultsFolders)) {
  
  active_species <- resultsNames[i]
  
  cat("\n")
  cat("# ---------------------------------\n")
  cat(i,"out of",length(resultsFolders),"\n")       # progress bar
  cat("# ---------------------------------\n")
  cat("\n")
  
  df <- dataStructureResult %>% 
  filter(sigmaNovelty > 2) %>%
  select(-cell, -analogNovelty, -analogDisappearance, -sigmaDisappearance, -analogNovelty.x, 
         -analogNovelty.y, -analogDisappearance.x, -analogDisappearance.y, -analogNoveltyDist,
         -analogDisappearanceDist, -refugiaWithinFocalPoll, -refugia) %>%
  mutate(Species = active_species)

BigDF <-rbind(BigDF, df)

}

## Run the EEZ and hex code below on the server

# assign EEZ to each point of coordinates
BigDF_sp <- SpatialPointsDataFrame(coords = BigDF[, c("x", "y")], data = BigDF)
proj4string(BigDF_sp) <- CRS("+proj=longlat +datum=WGS84")
eez_index <- over(BigDF_sp[, c('x', 'y')], eez)
BigDF <- cbind(BigDF, eez_index[2])
colnames(BigDF)[ncol(BigDF)] <- "EEZ" 

# assign hex cells to each point of coordinates
df <- df[(df$x<177) & (df$x>-177),]
df$cell <- dgGEO_to_SEQNUM(hex_grid, df$x, df$y)$seqnum
BigDF <- left_join(BigDF, Species_groups, by = "Species") %>%
  relocate(Group, .after = Species)


# Averaging the dissimilarity per EEZ 

BigDFxEEZ <- BigDF %>%
  select(-cell, -Group, -x, -y) %>%
  group_by(EEZ, Species) %>% arrange(EEZ) %>%
  summarise(Average = round(mean(sigmaNovelty), 3), .groups = "drop") %>%
  group_by(EEZ) %>% 
  summarise(Species_count = n_distinct(Species), Average = round(mean(Average), 3), 
            Species = paste(unique(Species), collapse = ",")
            , .groups = "drop")


write.csv(BigDFxEEZ, file="../Results/MeanDiss2EEZ585.csv", col.names = FALSE, row.names = FALSE)

### MAPPING DISSIM > 2 --------------

summary(active_DF)

# setting up the spatial df
active_DF <- na.omit(BigDFxEEZ)
active_eezs <- eez[eez$EEZ %in% unique(active_DF$EEZ),]
active_eezs@data <- cbind(active_eezs@data, active_DF[,-1])
active_eezs_vis <- st_as_sf(active_eezs)

# assign the columns we want to map
factors <- c("Average", "Species_count")

hex = FALSE
# hex resolution changed in setup

plot_data_factor = function (column, hex=FALSE) {
  if (startsWith(column,'S')) {
    scale_max = 100
  } else {
    scale_max = 8.3
  }
  if (hex==TRUE) {
    ggplot(hex_map_all) + 
      geom_map(
        data = world, map = world,
        aes(map_id = region),
        color = "grey", fill = "lightgray", size = 0.01
      ) +
      geom_sf(data = active_eezs_vis, aes(fill=NULL)) +
      geom_sf(aes(fill=hex_map_all[[column]]), alpha=0.8, linewidth = 0.01, color='white')    +
      labs(fill = column) + xlab("Longitude") +
      scale_x_continuous(labels = c("-120", "-60", "0", "60", "120")) +
      scale_fill_gradientn(colours = c("blue", "green", "yellow", "red"),
                                              limits = c(0, scale_max), na.value = 'grey') +
      theme(plot.margin=unit(c(0, 0, 0, 0), "pt"))
  } else {
    ggplot(active_eezs_vis) +
      geom_map(
        data = world, map = world,
        aes(map_id = region),
        color = "grey", fill = "lightgray", linewidth = 0.01
      ) +
      geom_sf(aes(fill=active_eezs_vis[[column]])) +
      labs(fill = column) + xlab("Longitude") +
      scale_fill_gradientn(colours = c("blue", "green", "yellow", "red"),
                                              limits = c(2, scale_max), na.value = 'grey') +
      scale_x_continuous(labels = c("-120", "-60", "0", "60", "120")) +
      theme(plot.margin=unit(c(0, 0, 0, 0), "pt"))
  }
}


EEZmap <- lapply(factors, plot_data_factor, hex=hex)
title <- ggdraw() + 
  draw_label(
    "Mean dissimilarity >2 and number of species with mean dissimilarity >2 per EEZ for SSP5-8.5",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(15, 0, 0, 7)
  )

# grid plots together including title as plot
EEZgrid <- plot_grid(plotlist = EEZmap)
save_grid <- plot_grid(title,EEZgrid, ncol = 1, axis='b',rel_heights = c(0.05,0.95))

save_grid

dpi = 300
ggsave(save_grid, path = "../Results/Figures/", filename ="Mean_diss_spp_count585.jpeg", height=1200/dpi,width=3000/dpi,dpi=dpi, device = 'jpeg')

