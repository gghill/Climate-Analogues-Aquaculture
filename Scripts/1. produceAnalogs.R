## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
##      
##                  Mapping climate analogs
##
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------

# First draft for mapping novel and disappearing climates AM, working from modified climate dissimilarity script
# First loop using test data set for aggregating dissimilarity

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

closeAllConnections()
rm(list=(ls()[ls()]))
gc(reset=TRUE)

library(raster)
library(sf)
library(sp)
library(ggplot2)
library(gridExtra)
library(geodata)
library(tidyverse)
library(cowplot)
require('dggridR')
library(tidyr)

## SETUP -----------------

# oceanBasins <- shapefile("../Data/spatialInformation/goas_v01.shp")
eez <- shapefile("../Data/spatialInformation/eez.shp")
shape <- raster("../Data/MaskRes025.tif")
world <- map_data("world")
# hex_grid <- dgconstruct(spacing=200, metric=TRUE, resround='down')
hex_grid <- dgconstruct(res = 6)

# Simplify and to make valid polygons
# sf::sf_use_s2(FALSE)
# oceanBasins <- st_as_sf(oceanBasins)
# eez <- st_as_sf(eez)

# oceanBasins <- st_simplify(oceanBasins, preserveTopology=TRUE, dTolerance = 0.001)
# eez <- st_simplify(eez, preserveTopology=TRUE, dTolerance = 0.001)
# oceanBasins <- st_make_valid(oceanBasins)
# eez <- st_make_valid(eez)

# oceanBasins <- as_Spatial(oceanBasins)
# eez <- as_Spatial(eez)

## -----------------------
## Subset eezs
# names(eez)
# myeez <- c("FRA","ESP","PRT","ITA")
# eez <- eez[ eez$cntryCd %in% myeez,]

# finalExtent <- c(-30,60,20,75)
# eez <- crop(eez,extent(finalExtent))
# plot(eez)

## -----------------------
# load in data 

load("../Results/climateDissimilarity/ssp585/Salmo salar/climateDissimilarity.RData")
sal <- dataStructureResult
  
# load("../Results/climateDissimilarity/ssp585/Salmo salar/climateDissimilarity.RData")
# s <- dataStructureResult

# plot function for summary plots
plot_data_factor = function (column, hex=FALSE) {
  if (startsWith(column,'p')) {
    scale_max = 1
  } else {
    scale_max = 8.3
  }
  if (hex==TRUE) {
    ggplot(hex_active_plot) + 
      geom_map(
        data = world, map = world,
        aes(map_id = region),
        color = "grey", fill = "lightgray", size = 0.01
      ) +
      geom_sf(data = active_eezs_vis, aes(fill=NULL)) +
      geom_sf(aes(fill=hex_active_plot[[column]]), alpha=0.8, linewidth = 0.01, color='white')    +
      labs(fill = column) + xlab("Longitude") +
      #scale_x_continuous(labels = c("-120", "-60", "0", "60", "120")) +
      scale_fill_gradient(limits = c(0, scale_max), low = '#5BA300', high = '#B51963', na.value = 'white') +
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
    scale_x_continuous(labels = c("-120", "-60", "0", "60", "120")) +
    scale_fill_gradient(limits = c(0, scale_max), low = '#5BA300', high = '#B51963', na.value = 'white') +
    theme(plot.margin=unit(c(0, 0, 0, 0), "pt"))
  }
}


# Begin dissimilarity visualization and aggregation ----

resultsScenario <- list.files("../Results/climateDissimilarity/", recursive=FALSE, full.names=TRUE)
resultsScenario <- gsub("../Results/climateDissimilarity/","",resultsScenario)
resultsScenario


Sys.time()
for( h in resultsScenario) {
  
  scenario <- h
  
## List available results
#scenario <- "ssp119"
hex <- TRUE
# can add outer loop here to iterate over multiple scenarios if needed

# species <- 'Salmo salar'
#resultsFolderMain <- paste0("../Results/climateDissimilarity/", scenario)
resultsFolderMain <- paste0("../Results/climateDissimilarity/", scenario)
resultsScope <- "sigmaNovelty" # sigmaNovelty sigmaDisappearance 
resultsFolders <- list.files(resultsFolderMain, recursive=TRUE, pattern=".RData", full.names=TRUE)
resultsNames <- list.files(resultsFolderMain, recursive=TRUE, pattern=".RData", full.names=FALSE)
resultsNames <- gsub("/climateDissimilarity.RData","",resultsNames)


#resultsFolders
#resultsNames

# test dataset
#test_main <- paste0(resultsFolderMain,'/','test_species')
#test_folders <- list.files(test_main, recursive=TRUE, pattern=".RData", full.names=TRUE)
#testNames <- list.files(test_main, recursive=TRUE, pattern=".RData", full.names=FALSE)
#testNames <- gsub("/climateDissimilarity.RData","",testNames)

#testNames
#resultsFolders <- test_folders
#resultsNames <- testNames

## BIG LOOP -----------------------

## Make a data.frame ready to be populated by species / eez / ocean
myResultsDF <- data.frame()

# determine factors of interest for plotting
# in this case p2 and p4 are the percent above 2 sigma and 4 sigma within that EEZ respectively
factors <- c('p2', 'p4', 'Average', 'Max')

# test run 060223:
# 1 hour run time for 4 species
# CSVs look good, plots named wrong (fixed)
# needs to be optimized

# test run 090223:
# hex = TRUE
# scenario not updated properly in folder path
# plots being labeled correctly, but dataframe and plot does not reflect the different scenario

for( i in 1:length(resultsFolders)) {
  active_species <- resultsNames[i]
  #cat("\014")
  cat("\n")
  cat("# ---------------------------------\n")
  cat(i,"out of",length(resultsFolders),"\n")
  cat("# ---------------------------------\n")
  cat("\n")
  
  DissimilarityData <- load(resultsFolders[i])
  #DissimilarityData <- load(paste0(resultsFolderMain,'/',species,'/climateDissimilarity.Rdata'))
  DissimilarityData <- dataStructureResult
  # head(DissimilarityData)
  
  # creates a raster where pixel value = no. of points in each cell
  ROIRaster <- rasterize(DissimilarityData[,c("x","y")], shape, field=DissimilarityData[,resultsScope], fun=mean) 
  # writeRaster(ROIRaster, gsub(".RData",".tif",resultsFolders[i]), overwrite=TRUE)
  
  # Crop to match only the extent of where the species is
  ROIPts <- xyFromCell(ROIRaster,Which(!is.na(ROIRaster), cells=TRUE)) 
  # oceanBasins.i <- crop(oceanBasins, extent(c(min(ROIPts[,1])-1,max(ROIPts[,1])+1,min(ROIPts[,2])-1,max(ROIPts[,2])+1)))
  eez.i <- crop(eez, extent(c(min(ROIPts[,1])-1,max(ROIPts[,1])+1,min(ROIPts[,2])-1,max(ROIPts[,2])+1)))
  
  # iterate between all oceans and all eez inside
  
  # Not using oceanBasins for now, large files, difficulty iterating over
  # for( j in 1:nrow(oceanBasins.i) ) {
  #   
  #   ROIRaster.j <- crop(ROIRaster,oceanBasins.i[j,])
  #   ROIRaster.j <- mask(ROIRaster.j,oceanBasins.i[j,])
  
  for( k in 1:nrow(eez.i) ) {
    
    ROIRaster.k <- crop(ROIRaster,eez.i[k,])
    ROIRaster.k <- crop(ROIRaster.k,eez.i[k,])
    
    valuesInsideOceanInsideEEZ <- getValues(ROIRaster.k)
    valuesInsideOceanInsideEEZ <- valuesInsideOceanInsideEEZ[!is.na(valuesInsideOceanInsideEEZ)]
    
    if( length(valuesInsideOceanInsideEEZ) > 0) {
      
      myResultsDF <- rbind(myResultsDF,
                           data.frame(Species=resultsNames[i],
                                      #Ocean=oceanBasins.i[j,]$name[1],
                                      EEZ=eez.i[k,]$EEZ[1],
                                      p2=(length(valuesInsideOceanInsideEEZ[valuesInsideOceanInsideEEZ>2])/length(valuesInsideOceanInsideEEZ)),
                                      p4=(length(valuesInsideOceanInsideEEZ[valuesInsideOceanInsideEEZ>4])/length(valuesInsideOceanInsideEEZ)),
                                      Average=mean(valuesInsideOceanInsideEEZ),
                                      Max=max(valuesInsideOceanInsideEEZ)
                           ))
    }
  }
  #}
  
# construct dataframe and sf object of active species
active_DF <- myResultsDF[myResultsDF$Species == active_species,]
active_eezs <- eez[eez$EEZ %in% unique(active_DF$EEZ),]
active_eezs@data <- cbind(active_eezs@data, active_DF[,-2])
active_eezs_vis <- st_as_sf(active_eezs)
# construct name and write df to file
write.csv(active_DF, file=paste0(resultsFolderMain,'/',active_species,'/',active_species,'_',scenario,"_DF.csv"), row.names = FALSE)

hex_active_DF <- DissimilarityData[,c('x','y','sigmaNovelty')]
# very hacky solution to wrap around plotting issue
hex_active_DF <- hex_active_DF[(hex_active_DF$x<177) & (hex_active_DF$x>-177),]

# assign pixels to cells and summarize dissimilarity within
hex_active_DF$cell <- dgGEO_to_SEQNUM(hex_grid, hex_active_DF$x, hex_active_DF$y)$seqnum
hex_active_cells <- dgSEQNUM_to_GEO(hex_grid,hex_active_DF$cell)
hex_test <- hex_active_DF %>% group_by(cell)
hex_active_dissim <- hex_active_DF %>% group_by(cell) %>% summarise(Average=mean(sigmaNovelty),
                                                                    Max=max(sigmaNovelty),
                                                                    p2=sum(sigmaNovelty > 2)/sum(sigmaNovelty),
                                                                    p4=sum(sigmaNovelty > 4)/sum(sigmaNovelty))
hex_active_plot <- dgcellstogrid(hex_grid, hex_active_dissim$cell)
names(hex_active_plot) <- c('cell', 'geometry') # labeling with highest overlap EEZ would be nice here
hex_active_plot <- merge(hex_active_plot, hex_active_dissim, by.x='cell', by.y='cell')

# add EEZ labels by row 
hex_active_DF_sp <- SpatialPointsDataFrame(coords = hex_active_DF[, c("x", "y")], data = hex_active_DF)
proj4string(hex_active_DF_sp) <- CRS("+proj=longlat +datum=WGS84")
eez_index <- over(hex_active_DF_sp[, c('x', 'y')], eez)
hex_active_DF <- cbind(hex_active_DF, eez_index[2])
colnames(hex_active_DF)[ncol(hex_active_DF)] <- "EEZ" 

write.csv(hex_active_DF, file=paste0(resultsFolderMain,'/',active_species,'/',active_species,'_',scenario,"_hex_DF.csv"), row.names = FALSE)


myplots <- lapply(factors, plot_data_factor, hex=hex)
title <- ggdraw() + 
  draw_label(
    paste(active_species, scenario),
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

# write to file
dpi=300
if (hex==TRUE) {
  ggsave(paste0(resultsFolderMain,'/',active_species,'/',active_species,'_',scenario,'_EEZ_hex_plots.jpeg'),height=1200/dpi,width=3000/dpi,dpi=dpi, save_grid, device = 'jpeg')
  
} else {
  ggsave(paste0(resultsFolderMain,'/',active_species,'/',active_species,'_',scenario,'_EEZ_plots.jpeg'),height=1200/dpi,width=3000/dpi,dpi=dpi, save_grid, device = 'jpeg')
}
}
}
Sys.time()

# ~2 hours per test run (7 species) with hex data

# EEZ Aggregation ----

# Dissim metrics per EEZ averaged by all species in that EEZ
eezDF <- myResultsDF %>% group_by(EEZ) %>% 
  summarise(Species = n_distinct(Species), p2 = sum(Average > 2)/sum(Average),
            p4 = sum(Average > 4)/sum(Average),
            Average = mean(Average), Max = mean(Max), .groups = "drop")

# Dissim metrics per EEZ by species
eez_sp_mean <- myResultsDF %>% group_by(EEZ) %>% arrange(EEZ) %>% 
  pivot_wider(names_from = Species, values_from = Average) %>% select(-p2, -p4, -Max) %>%
  summarise(across(everything(), ~first(na.omit(.))))

eez_sp_p2 <- myResultsDF %>% group_by(EEZ) %>% arrange(EEZ) %>% 
  pivot_wider(names_from = Species, values_from = p2) %>% select(-Average, -p4, -Max) %>%
  summarise(across(everything(), ~first(na.omit(.)))) %>%
  mutate(prop_sp = rowMeans(across(-EEZ) > 0, na.rm = TRUE))

eez_sp_p4 <- myResultsDF %>% group_by(EEZ) %>% arrange(EEZ) %>% 
  pivot_wider(names_from = Species, values_from = p4) %>% select(-Average, -p2, -Max) %>%
  summarise(across(everything(), ~first(na.omit(.))))%>%
  mutate(prop_sp = rowMeans(across(-EEZ) > 0, na.rm = TRUE))
              
# Climate Analog visualisation ----

sal <- s[1:1000,]
sal1<- sal[1:100,]
sal1<- sal[1:1000,]

names(sal)[names(sal1) == 'x'] <- 'long'
names(sal)[names(sal1) == 'y'] <- 'lat'

# Mapping Salmo salar subset SSP585

countries <- world(path = "../Data/countries")

download_window <- ext(c(-60, 60, 50, 80))
plot(countries, ext = download_window)
points(salar585_subset[ , c("x", "y")], pch = 20, col = "red")
points(salar585_subset[ , c("analogNovelty.x", "analogNovelty.y")], pch = 20, col = "blue")
points(salar585_subset[ , c("analogDisappearance.x", "analogDisappearance.y")], pch = 20, col = "green")


# Combine coords into one df (but can't connect points based on cells)
current <- data.frame(long = sal1$long, lat = sal1$lat)
current$climate <- "current"
novel <- data.frame(long = sal1$analogNovelty.x, lat = sal1$analogNovelty.y)
novel$climate <- "novel"
dissappearing <- data.frame(long = sal1$analogDisappearance.x, lat = sal1$analogDisappearance.y)
dissappearing$climate <- "dissappearing"

dat <- rbind(current, novel, dissappearing) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)


ggplot() + 
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) + theme_bw()  +
  xlim(-80, 80) + ylim(10, 80) +
  geom_sf(data = dat, aes(color = climate)) +
  labs(title = "Salmo salar SSP5-8.5", x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_brewer(palette = "Dark2", name = "Climate", 
                      labels = c("Current",  "Disappearing", "Novel")) 

ggsave(path = "../Results/climateDissimilarity/ssp585/Salmo salar/", filename ="S.salar_subset585.png", width = 6, height = 4)


# Adding lines to connect climates

# new df with cell, lat and long values  
df <- data.frame(cell = sal1$cell, c.long = sal1$long, c.lat = sal1$lat, n.long = sal1$analogNovelty.x, 
                 n.lat = sal1$analogNovelty.y, d.long = sal1$analogDisappearance.x, 
                 d.lat = sal1$analogDisappearance.y)

#  create new data frames for the connecting lines added between coord cols - geom_path extracting x and y length (300) must be the same as the df length (100)
line_data <- data.frame(x = c(df$c.long, df$n.long), y = c(df$c.lat, df$n.lat), cell = df$cell)
line_data2 <- data.frame(x = c(df$c.long, df$d.long), y = c(df$c.lat, df$d.lat), cell = df$cell)

ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) + theme_bw()  +
  geom_path(data = line_data, aes(x = x, y = y, group = cell), color = "#B2ABD2", alpha = 0.4) +
  geom_path(data = line_data2, aes(x = x, y = y, group = cell), color = "#F4A582", alpha = 0.4) +
  geom_point(data=df, aes(x = c.long, y = c.lat, color = "current")) +
  geom_point(data=df, aes(x = n.long, y = n.lat, color = "novel")) +
  geom_point(data=df, aes(x = d.long, y = d.lat, color = "disappearing")) + 
  labs(title = "Salmo salar SSP5-8.5", x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(-80, 80) + ylim(10, 80) +
  scale_colour_brewer(palette = "Dark2", name = "Climate", 
                      labels = c("Current",  "Disappearing", "Novel")) 

ggsave(path = "../Results/climateDissimilarity/ssp585/Salmo salar/", filename ="S.salar_subset_lines585.png", width = 6, height = 4)

# Mapping Salmo salar full dataset

full_dat <- data.frame(cell = sal$cell, c.long = sal$long, c.lat = sal$lat, 
                       n.long = sal$analogNovelty.x, 
                       n.lat = sal$analogNovelty.y, d.long = sal$analogDisappearance.x, 
                       d.lat = sal$analogDisappearance.y)

ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) + theme_bw()  +
  geom_point(data=full_dat, aes(x = c.long, y = c.lat, color = "current")) +
  geom_point(data=full_dat, aes(x = n.long, y = n.lat, color = "novel")) +
  geom_point(data=full_dat, aes(x = d.long, y = d.lat, color = "disappearing")) + 
  labs(title = "Salmo salar SSP5-8.5", x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(-60, 60) + ylim(10, 80) +
  scale_colour_brewer(palette = "Dark2", name = "Climate", 
                      labels = c("Current",  "Disappearing", "Novel")) 



