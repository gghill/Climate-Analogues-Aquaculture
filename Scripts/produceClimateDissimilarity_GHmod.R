## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
##      
##  biodiversityDS [ biodiversitydatascience.com ]
##
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------

# Modified by GH to create example analyses for submission with Oceanography Supplement
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

closeAllConnections()
rm(list=(ls()[ls()]))
gc(reset=TRUE)

library(raster)
library(sf)

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## -----------------

oceanBasins <- shapefile("../Data/spatialInformation/goas_v01.shp")
eez <- shapefile("../Data/spatialInformation/eez.shp")
shape <- raster("../Data/MaskRes025.tif")

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
## List available results
scenario <- 'ssp585'
resultsFolderMain <- paste0("../Results/climateDissimilarity/", scenario)
resultsScope <- "sigmaNovelty" # sigmaNovelty sigmaDisappearance 
resultsFolders <- list.files(resultsFolderMain, recursive=TRUE, pattern=".RData", full.names=TRUE)
resultsNames <- list.files(resultsFolderMain, recursive=TRUE, pattern=".RData", full.names=FALSE)
resultsNames <- gsub("/climateDissimilarity.RData","",resultsNames)

resultsFolders
resultsNames

## -----------------------
## Make a data.frame ready to be populated by species / eez / ocean
myResultsDF <- data.frame()

for( i in 1:length(resultsFolders)) {
  
  cat("\014")
  cat("\n")
  cat("# ---------------------------------\n")
  cat(i,"out of",length(resultsFolders),"\n")
  cat("# ---------------------------------\n")
  cat("\n")
  
  DissimilarityData <- loadRData(resultsFolders[i])
  # head(DissimilarityData)
  
  ROIRaster <- rasterize(DissimilarityData[,c("x","y")], shape, field=DissimilarityData[,resultsScope], fun=mean) 
  writeRaster(ROIRaster, gsub(".RData",".tif",resultsFolders[i]), overwrite=TRUE)
  
  # Crop to match only the extent of where the species is
  ROIPts <- xyFromCell(ROIRaster,Which(!is.na(ROIRaster), cells=TRUE)) 
  oceanBasins.i <- crop(oceanBasins, extent(c(min(ROIPts[,1])-1,max(ROIPts[,1])+1,min(ROIPts[,2])-1,max(ROIPts[,2])+1)))
  eez.i <- crop(eez, extent(c(min(ROIPts[,1])-1,max(ROIPts[,1])+1,min(ROIPts[,2])-1,max(ROIPts[,2])+1)))

  # iterate between all oceans and all eez inside
  
  for( j in 1:nrow(oceanBasins.i) ) {
  
    ROIRaster.j <- crop(ROIRaster,oceanBasins.i[j,])
    ROIRaster.j <- mask(ROIRaster.j,oceanBasins.i[j,])
    
    for( k in 1:nrow(eez.i) ) {
      
      ROIRaster.k <- crop(ROIRaster,eez.i[k,])
      ROIRaster.k <- crop(ROIRaster.k,eez.i[k,])
      
      valuesInsideOceanInsideEEZ <- getValues(ROIRaster.k)
      valuesInsideOceanInsideEEZ <- valuesInsideOceanInsideEEZ[!is.na(valuesInsideOceanInsideEEZ)]
            
              if( length(valuesInsideOceanInsideEEZ) > 0) {
                
                myResultsDF <- rbind(myResultsDF,
                                     data.frame(Species=resultsNames[i],
                                                Ocean=oceanBasins.i[j,]$name[1],
                                                EEZ=eez.i[k,]$EEZ[1],
                                                Average=mean(valuesInsideOceanInsideEEZ),
                                                Max=max(valuesInsideOceanInsideEEZ)
                                     ))
              }
      }
  }
}

write.csv(myResultsDF , file=paste0("../Results/",resultsScope,scenario,"_DF.csv"), row.names = FALSE)


## EEZ level-----------------------

# next steps: color EEZs by values, identify areas of interest (max and/or min sigma novelty)
# myResultsDF aggregated to EEZ, only average and max available per EEZ
myResultsDF <- read.csv(paste0("../Results/",resultsScope,scenario,"_DF.csv"))
countries <- world(path = "../Data/countries")
cols <- c('EEZ',	'seaBasn',	'country', 'Average',	'Max'
)


salmo <- myResultsDF[myResultsDF$Species=='Salmo salar',]
salmo_eezs <- eez[eez$EEZ %in% unique(salmo$EEZ),]
salmo_agg <- aggregate.data.frame(x=salmo, by=list(salmo$EEZ), FUN = 'mean')
salmo_agg <- salmo_agg[,c('Group.1', 'Average', 'Max')]
salmo_eezs@data <- cbind(salmo_eezs@data, salmo_agg)
write_rds(salmo_eezs, paste0('../Results/salmon_eezs_',scenario,'.rds'))
write.csv(as.data.frame(salmo_eezs[,cols]), paste0('../Results/salmon_eezs_',scenario,'.csv'))

chanos <- myResultsDF[myResultsDF$Species=='Chanos chanos',]
chanos_eezs <- eez[eez$EEZ %in% unique(chanos$EEZ),]
chanos_agg <- aggregate.data.frame(x=chanos, by=list(chanos$EEZ), FUN = 'mean')
chanos_agg <- chanos_agg[,c('Group.1', 'Average', 'Max')]
chanos_eezs@data <- cbind(chanos_eezs@data, chanos_agg)
write_rds(chanos_eezs, paste0('../Results/chanos_eezs_',scenario,'.rds'))
write.csv(as.data.frame(chanos_eezs[,cols]), paste0('../Results/chanos_eezs_',scenario,'.csv'))

trout <- myResultsDF[myResultsDF$Species=='Oncorhynchus mykiss',]
trout_eezs <- eez[eez$EEZ %in% unique(trout$EEZ),]
trout_agg <- aggregate.data.frame(x=trout, by=list(trout$EEZ), FUN = 'mean')
trout_agg <- trout_agg[,c('Group.1', 'Average', 'Max')] # duplicated EEZ column
trout_eezs@data <- cbind(trout_eezs@data, trout_agg)
write_rds(trout_eezs, paste0('../Results/trout_eezs_',scenario,'.rds'))
write.csv(as.data.frame(trout_eezs[,cols]), paste0('../Results/trout_eezs_',scenario,'.csv'))

# 
# pal <- colorRamp(c("green", "yellow", "red"))    # 1) choose colors
# col <- rgb(pal((salmo$Average - min(salmo$Average)) / diff(range(salmo$Average))), max=255)
# AVERAGE
plot(countries)
plot(salmo_eezs, col=col, add=T)

plot(salmo_eezs[salmo_eezs$EEZ %in% unique(salmo[salmo$Average>2,]$EEZ),], col = 'yellow', add=T)
plot(salmo_eezs[salmo_eezs$EEZ %in% unique(salmo[salmo$Average>4,]$EEZ),], col = 'red', add=T)

# Plot at EEZ level using ggplot
library(ggplot2)
require(gridExtra)
library(geodata)
library(tidyverse)
library(sf)

# Salmo ----
df <- st_as_sf(salmo_eezs)
df_hot_salmon <- readRDS('../Results/salmon_eezs_ssp585.rds')
df_hot_salmon <- st_as_sf(df_hot_salmon)
world <- map_data("world")

avg_salmo <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data=df, aes(fill=Average)) +
  # labs(title = 'Average sigma dissimilarity by EEZ (2100, SSP 119') +
  labs(title = 'Atlantic Salmon (Salmo salar) SSP 119 (Paris Agreement)') +
  scale_fill_gradient(limits = c(0, max(df_hot_chanos$Max)), low = '#5BA300', high = '#B51963')

avg_hot_salmo <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data=df_hot_salmon, aes(fill=Average)) +
  # labs(title = 'Average sigma dissimilarity by EEZ (2100, SSP 119') +
  labs(title = 'Atlantic Salmon (Salmo salar) SSP 585') +
  scale_fill_gradient(limits = c(0, max(df_hot_chanos$Max)), low = '#5BA300', high = '#B51963')

max_salmo <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data=df, aes(fill=Max)) +
  labs(title = 'Atlantic Salmon (Salmo salar) SSP 119 (Paris Agreement)') +
  #labs(title = 'Max sigma dissimilarity by EEZ (2100, SSP 119)') +
  scale_fill_gradient(limits = c(0, max(df_hot_chanos$Max)), low = '#5BA300', high = '#B51963')

max_hot_salmo <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data=df_hot_salmon, aes(fill=Max)) +
  labs(title = 'Atlantic Salmon (Salmo salar) SSP 585') +
  #labs(title = 'Max sigma dissimilarity by EEZ (2100, SSP 119)') +
  scale_fill_gradient(limits = c(0, max(df_hot_chanos$Max)), low = '#5BA300', high = '#B51963')

grid.arrange(avg_salmo, max_salmo,  nrow=2, ncol=1, top = 'Salmo Salar')

# Chanos ----
df_chanos <- st_as_sf(chanos_eezs)
df_hot_chanos <- readRDS('../Results/chanos_eezs_ssp585.rds')
df_hot_chanos <- st_as_sf(df_hot_chanos)

avg_chanos <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data=df_chanos, aes(fill=Average)) +
  # labs(title = 'Average sigma dissimilarity by EEZ (2100, SSP 119') +
  labs(title = 'Milkfish (Chanos chanos) SSP 119 (Paris Agreement)') +
  scale_fill_gradient(limits = c(0, max(df_hot_chanos$Max)), low = '#5BA300', high = '#B51963')

avg_hot_chanos <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data=df_hot_chanos, aes(fill=Average)) +
  # labs(title = 'Average sigma dissimilarity by EEZ (2100, SSP 119') +
  labs(title = 'Milkfish (Chanos chanos) SSP 585') +
  scale_fill_gradient(limits = c(0, max(df_hot_chanos$Max)), low = '#5BA300', high = '#B51963')

max_chanos <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data=df_chanos, aes(fill=Max)) +
  #labs(title = 'Max sigma dissimilarity by EEZ (2100, SSP 119)') +
  labs(title = 'Milkfish (Chanos chanos) SSP 119 (Paris Agreement)') +
  scale_fill_gradient(limits = c(0, max(df_hot_chanos$Max)), low = '#5BA300', high = '#B51963')

max_hot_chanos <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data=df_hot_chanos, aes(fill=Max)) +
  #labs(title = 'Max sigma dissimilarity by EEZ (2100, SSP 119)') +
  labs(title = 'Milkfish (Chanos chanos) SSP 585') +
  scale_fill_gradient(limits = c(0, max(df_hot_chanos$Max)), low = '#5BA300', high = '#B51963')

grid.arrange(avg_chanos, max_chanos,  nrow=2, ncol=1, top = 'Chanos chanos')

# trout ----

df_trout <- st_as_sf(trout_eezs)
df_hot_trout <- readRDS('../Results/trout_eezs_ssp585.rds')
df_hot_trout <- st_as_sf(df_hot_trout)

avg_trout <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data=df_trout, aes(fill=Average)) +
  labs(title = 'Rainbow trout (Oncorhynchus mykiss) SSP 119 (Paris Agreement)') +
  scale_fill_gradient(limits = c(0, max(df_hot_chanos$Max)), low = '#5BA300', high = '#B51963')

avg_hot_trout <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data=df_hot_trout, aes(fill=Average)) +
  # labs(title = 'Average sigma dissimilarity by EEZ (2100, SSP 585)') +
  labs(title = 'Rainbow trout (Oncorhynchus mykiss) SSP 585') +
  scale_fill_gradient(limits = c(0, max(df_hot_chanos$Max)), low = '#5BA300', high = '#B51963')

max_trout <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data=df_trout, aes(fill=Max)) +
  # labs(title = 'Max sigma dissimilarity by EEZ (2100, SSP 119)') +
  labs(title = 'Rainbow trout (Oncorhynchus mykiss) SSP 119 (Paris Agreement)') +
  scale_fill_gradient(limits = c(0, max(df_hot_chanos$Max)), low = '#5BA300', high = '#B51963')

max_hot_trout <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data=df_hot_trout, aes(fill=Max)) +
  # labs(title = 'Max sigma dissimilarity by EEZ (2100, SSP 119)') +
  labs(title = 'Rainbow trout (Oncorhynchus mykiss) SSP 585') +
  scale_fill_gradient(limits = c(0, max(df_hot_chanos$Max)), low = '#5BA300', high = '#B51963')

grid.arrange(avg_trout, max_trout,  nrow=2, ncol=1, top = 'Oncorhynchus mykiss')

grid.arrange(max_salmo, max_chanos, max_trout, nrow=3, ncol=1, top = 'Cultured finfish Max Sigma Dissimilarity by EEZ (2100, SSP 119)')

grid.arrange(max_hot_salmo, max_hot_chanos, max_hot_trout, nrow=3, ncol=1, top = 'Cultured finfish Max Sigma Dissimilarity by EEZ (2100, SSP 585)')


grid.arrange(max_salmo, max_hot_salmo, max_chanos, max_hot_chanos, max_trout, max_hot_trout, nrow=3, ncol=2, 
             top = 'Cultured finfish end of century Max Sigma Dissimilarity by EEZ under two climate change scenarios')

grid.arrange(avg_salmo, avg_hot_salmo, avg_chanos, avg_hot_chanos, avg_trout, avg_hot_trout, nrow=3, ncol=2, 
             top = 'Cultured finfish end of century Avg Sigma Dissimilarity by EEZ under two climate change scenarios')

# full dissim data set ----
# not aggregating immediately to EEZ level
require('dggridR')
hex_grid <- dgconstruct(spacing=200, metric=TRUE, resround='down')
data_folder <- "../Results/climateDissimilarity"

# Salmo full dataset ----
scenario <- 'ssp119'
load(paste0(data_folder, '/',scenario, '/Salmo salar/climateDissimilarity.RData'))
salmo_dis <- dataStructureResult
scenario <- 'ssp585'
load(paste0(data_folder, '/',scenario, '/Salmo salar/climateDissimilarity.RData'))
salmo_hot_dis  <- dataStructureResult
# salmo_dis <- salmo_dis[(salmo_dis$x >= -30) & (salmo_dis$x <= 50),] # restrict to a test area
# salmo_dis <- salmo_dis[(salmo_dis$y >= 45) & (salmo_dis$y <= 85),]
# plot(salmo_ras, xlim=c(-30,50), ylim=c(45,85))


salmo_df <- salmo_dis[,c('x','y','sigmaNovelty')]
salmo_hot_df <- salmo_hot_dis[,c('x','y','sigmaNovelty')]

# plot to grid

salmo_hot_df$cell <- dgGEO_to_SEQNUM(hex_grid, salmo_hot_df$x, salmo_hot_df$y)$seqnum
salmo_hot_cells   <- dgSEQNUM_to_GEO(hex_grid,salmo_hot_df$cell)
salmo_avg_hex_hot_dissim   <- salmo_hot_df %>% group_by(cell) %>% summarise(Dissim=mean(sigmaNovelty))
salmo_plot_hot_grid <- dgcellstogrid(hex_grid, salmo_avg_hex_hot_dissim$cell)
names(salmo_plot_hot_grid) <- c('cell', 'geometry')
salmo_plot_hot_grid <- merge(salmo_plot_hot_grid, salmo_avg_hex_hot_dissim, by.x='cell', by.y='cell')
salmo_plot_hot_grid <- as.data.frame(salmo_plot_hot_grid)

salmon_gg_hot_grid <- ggplot() + 
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data = st_as_sf(salmo_eezs), aes(fill=NULL)) +
  geom_sf(data=st_as_sf(salmo_plot_hot_grid), aes(fill=Dissim), alpha=0.8, color='white')    +
  #geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  xlim(-100,100) +
  ylim(30,80) +
  scale_fill_gradient(limits = c(0, max(chanos_hot_df$sigmaNovelty)), low = '#5BA300', high = '#B51963') +
  labs(title = 'Salmo salar 585',
       fill = 'Sigma Dissim.')

salmon_gg_hot_grid

salmo_df$cell <- dgGEO_to_SEQNUM(hex_grid, salmo_df$x, salmo_df$y)$seqnum
salmo_cells   <- dgSEQNUM_to_GEO(hex_grid,salmo_df$cell)
salmo_avg_hex_dissim   <- salmo_df %>% group_by(cell) %>% summarise(Dissim=mean(sigmaNovelty))
salmo_plot_grid <- dgcellstogrid(hex_grid, salmo_avg_hex_dissim$cell)
names(salmo_plot_grid) <- c('cell', 'geometry')
salmo_plot_grid <- merge(salmo_plot_grid, salmo_avg_hex_dissim, by.x='cell', by.y='cell')
salmo_plot_grid <- as.data.frame(salmo_plot_grid)

salmon_gg_grid <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data = st_as_sf(salmo_eezs), aes(fill=NULL)) +
  geom_sf(data=st_as_sf(salmo_plot_grid), aes(fill=Dissim), alpha=0.8, color='white')    +
  #geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  xlim(-100,100) +
  ylim(30,80) +
  scale_fill_gradient(limits = c(0, max(chanos_hot_df$sigmaNovelty)), low = '#5BA300', high = '#B51963') +
  labs(title = 'Salmo salar 119',
       fill = 'Sigma Dissim.')
salmon_gg_grid

# Pixel scale

salmon_plt <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_tile(data=salmo_df, aes(x=x,y=y,fill=sigmaNovelty)) +
  xlim(-180,180) +
  ylim(30,80) +
  scale_fill_gradient(limits = c(0, max(chanos_df$sigmaNovelty)), low = '#5BA300', high = '#B51963') +
  theme(axis.title.x=element_blank(), #remove x axis labels
        axis.title.y=element_blank(),  #remove y axis labels
  ) +
  labs(title = 'Salmo salar 119',
       fill = 'Sigma Dissim.')
salmon_plt

salmon_hot_plt <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_tile(data=salmo_hot_df, aes(x=x,y=y,fill=sigmaNovelty)) +
  xlim(-180,180) +
  ylim(30,80) +
  scale_fill_gradient(limits = c(0, max(chanos_hot_df$sigmaNovelty)), low = '#5BA300', high = '#B51963') +
  theme(axis.title.x=element_blank(), #remove x axis labels
        axis.title.y=element_blank(),  #remove y axis labels
  ) +
  labs(title = 'Salmo salar 585',
       fill = 'Sigma Dissim.')
salmon_plt

# Chanos full dataset ----
scenario <- 'ssp119'
load(paste0(data_folder, '/',scenario, '/Chanos chanos/climateDissimilarity.RData'))
chanos_dis <- dataStructureResult
scenario <- 'ssp585'
load(paste0(data_folder, '/',scenario, '/Chanos chanos/climateDissimilarity.RData'))
chanos_hot_dis <- dataStructureResult

chanos_df <- chanos_dis[,c('x','y','sigmaNovelty')]
chanos_df <- chanos_df[(chanos_df$x<175) & (chanos_df$x>-175),]
chanos_hot_df <- chanos_hot_dis[,c('x','y','sigmaNovelty')]
chanos_hot_df <- chanos_hot_df[(chanos_hot_df$x<175) & (chanos_hot_df$x>-175),]

# Chanos hex plots

chanos_hot_df$cell <- dgGEO_to_SEQNUM(hex_grid, chanos_hot_df$x, chanos_hot_df$y)$seqnum
chanos_hot_cells   <- dgSEQNUM_to_GEO(hex_grid,chanos_hot_df$cell)
chanos_avg_hex_hot_dissim   <- chanos_hot_df %>% group_by(cell) %>% summarise(Dissim=mean(sigmaNovelty))
chanos_plot_hot_grid <- dgcellstogrid(hex_grid, chanos_avg_hex_hot_dissim$cell)
names(chanos_plot_hot_grid) <- c('cell', 'geometry')
chanos_plot_hot_grid <- merge(chanos_plot_hot_grid, chanos_avg_hex_hot_dissim, by.x='cell', by.y='cell')
chanos_plot_hot_grid <- as.data.frame(chanos_plot_hot_grid)
# chanos_plot_hot_grid <- st_as_sf(chanos_plot_hot_grid)
# x1 <- st_crop(chanos_plot_hot_grid, ext(-180, 0, -90, 90))
# x2 <- st_crop(chanos_plot_hot_grid, ext(0, 180, -90, 90))   
# #extent(x1) <- c(-180, 180, -90, 90)
# chanos_plot_hot_grid <- merge(as.data.frame(x1),as.data.frame(x2))


chanos_gg_hot_grid <- ggplot() + 
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data = st_as_sf(chanos_eezs), aes(fill=NULL)) +
  geom_sf(data=st_as_sf(chanos_plot_hot_grid), aes(fill=Dissim), alpha=0.8, color='white')    +
  #geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  #xlim(0,180) +
  #ylim(30,80) +
  scale_fill_gradient(limits = c(0, max(chanos_hot_df$sigmaNovelty)), low = '#5BA300', high = '#B51963') +
  theme(legend.position = "none") +
  labs(title = 'Chanos chanos 585',
       fill = 'Sigma Dissim.')
chanos_gg_hot_grid


chanos_df$cell <- dgGEO_to_SEQNUM(hex_grid, chanos_df$x, chanos_df$y)$seqnum
chanos_cells   <- dgSEQNUM_to_GEO(hex_grid,chanos_df$cell)
chanos_avg_hex_dissim   <- chanos_df %>% group_by(cell) %>% summarise(Dissim=mean(sigmaNovelty))
chanos_plot_grid <- dgcellstogrid(hex_grid, chanos_avg_hex_dissim$cell)
names(chanos_plot_grid) <- c('cell', 'geometry')
chanos_plot_grid <- merge(chanos_plot_grid, chanos_avg_hex_dissim, by.x='cell', by.y='cell')
chanos_plot_grid <- as.data.frame(chanos_plot_grid)


chanos_gg_grid <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data = st_as_sf(chanos_eezs), aes(fill=NULL)) +
  geom_sf(data=st_as_sf(chanos_plot_grid),      aes(fill=Dissim), alpha=0.8, color='white')    +
  #geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  # xlim(-100,100) +
  # ylim(30,80) +
  scale_fill_gradient(limits = c(0, max(chanos_hot_df$sigmaNovelty)), low = '#5BA300', high = '#B51963') +
  theme(legend.position = "none") +
  labs(title = 'Chanos chanos 119',
       fill = 'Sigma Dissim.')

# chanos_gg_grid

# Pixel based

chanos_plt <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data = st_as_sf(chanos_eezs), aes(fill=NULL)) +
  geom_tile(data=chanos_df, aes(x=x,y=y,fill=sigmaNovelty)) +
  scale_fill_gradient(limits = c(0, max(chanos_hot_df$sigmaNovelty)), low = '#5BA300', high = '#B51963') +
  theme(axis.title.x=element_blank(), #remove x axis labels
        axis.title.y=element_blank(),  #remove y axis labels
  ) +
  labs(title = 'Chanos chanos 119',
       fill = 'Sigma Dissim.')
chanos_plt

chanos_hot_plt <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data = st_as_sf(chanos_eezs), aes(fill=NULL)) +
  geom_tile(data=chanos_hot_df, aes(x=x,y=y,fill=sigmaNovelty)) +
  scale_fill_gradient(limits = c(0, max(chanos_hot_df$sigmaNovelty)), low = '#5BA300', high = '#B51963') +
  theme(axis.title.x=element_blank(), #remove x axis labels
        axis.title.y=element_blank(),  #remove y axis labels
  ) +
  labs(title = 'Chanos chanos 585',
       fill = 'Sigma Dissim.')

# trout full dataset ----
scenario <- 'ssp119'
load(paste0(data_folder, '/',scenario, '/Oncorhynchus mykiss/climateDissimilarity.RData'))
trout_dis <- dataStructureResult
trout_df <- trout_dis[,c('x','y','sigmaNovelty')]
trout_df <- trout_df[(trout_df$x<175) & (trout_df$x>-175),]
scenario <- 'ssp585'
load(paste0(data_folder, '/',scenario, '/Oncorhynchus mykiss/climateDissimilarity.RData'))

trout_hot_dis <- dataStructureResult
trout_hot_df <- trout_dis[,c('x','y','sigmaNovelty')]
trout_hot_df <- trout_hot_df[(trout_hot_df$x<175) & (trout_hot_df$x>-175),]

# gridded

trout_hot_df$cell <- dgGEO_to_SEQNUM(hex_grid, trout_hot_df$x, trout_hot_df$y)$seqnum
trout_hot_cells   <- dgSEQNUM_to_GEO(hex_grid,trout_hot_df$cell)
trout_avg_hex_hot_dissim   <- trout_hot_df %>% group_by(cell) %>% summarise(Dissim=mean(sigmaNovelty))
trout_plot_hot_grid <- dgcellstogrid(hex_grid, trout_avg_hex_hot_dissim$cell)
names(trout_plot_hot_grid) <- c('cell', 'geometry')
trout_plot_hot_grid <- merge(trout_plot_hot_grid, trout_avg_hex_hot_dissim, by.x='cell', by.y='cell')
trout_plot_hot_grid <- as.data.frame(trout_plot_hot_grid)

trout_gg_hot_grid <- ggplot() + 
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data = st_as_sf(trout_eezs), aes(fill=NULL), lwd=.1) +
  geom_sf(data=st_as_sf(trout_plot_hot_grid), aes(fill=Dissim), alpha=0.8, color='white') +
  #geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  # xlim(-100,100) +
  # ylim(30,80) +
  scale_fill_gradient(limits = c(0, max(chanos_hot_df$sigmaNovelty)), low = '#5BA300', high = '#B51963') +
  theme(legend.position = "none") +
  labs(title = 'Oncorhynchus mykiss 585',
       fill = 'Sigma Dissim.')

# trout_gg_hot_grid

trout_df$cell <- dgGEO_to_SEQNUM(hex_grid, trout_df$x, trout_df$y)$seqnum
trout_cells   <- dgSEQNUM_to_GEO(hex_grid,trout_df$cell)
trout_avg_hex_dissim   <- trout_df %>% group_by(cell) %>% summarise(Dissim=mean(sigmaNovelty))
trout_plot_grid <- dgcellstogrid(hex_grid, trout_avg_hex_dissim$cell)
names(trout_plot_grid) <- c('cell', 'geometry')
trout_plot_grid <- merge(trout_plot_grid, trout_avg_hex_dissim, by.x='cell', by.y='cell')
trout_plot_grid <- as.data.frame(trout_plot_grid)

trout_gg_grid <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data = st_as_sf(trout_eezs), aes(fill=NULL)) +
  geom_sf(data=st_as_sf(trout_plot_grid),      aes(fill=Dissim), alpha=0.8, color='white')    +
  #geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  # xlim(-100,100) +
  # ylim(30,80) +
  scale_fill_gradient(limits = c(0, max(chanos_hot_df$sigmaNovelty)), low = '#5BA300', high = '#B51963') +
  theme(legend.position = "none") +
  labs(title = 'Oncorhynchus mykiss 119',
       fill = 'Sigma Dissim.')

# trout_gg_grid

# Pixel based

trout_plt <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data = st_as_sf(trout_eezs), aes(fill=NULL)) +
  geom_tile(data=trout_df, aes(x=x,y=y,fill=sigmaNovelty)) +
  scale_fill_gradient(limits = c(0, max(chanos_hot_df$sigmaNovelty)), low = '#5BA300', high = '#B51963') +
  xlim(-180,180) +
  ylim(30,80) +
  theme(axis.title.x=element_blank(), #remove x axis labels
        axis.title.y=element_blank(),  #remove y axis labels
  ) +
  labs(title = 'Trout 119',
       fill = 'Sigma Dissim.')
trout_plt

trout_hot_plt <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf(data = st_as_sf(trout_eezs), aes(fill=NULL)) +
  geom_tile(data=trout_hot_df, aes(x=x,y=y,fill=sigmaNovelty)) +
  scale_fill_gradient(limits = c(0, max(chanos_hot_df$sigmaNovelty)), low = '#5BA300', high = '#B51963') +
  xlim(-180,180) +
  ylim(30,80) +
  theme(axis.title.x=element_blank(), #remove x axis labels
        axis.title.y=element_blank(),  #remove y axis labels
  ) +
  labs(title = 'Trout 585',
       fill = 'Sigma Dissim.')

# full dissim grids ----
# pixel plots
full_119 <- grid.arrange(salmon_plt, chanos_plt, trout_plt, nrow=3, ncol=1, top='Sigma dissimilarity by 2100 under SSP 119 (Paris Agreement)')

full_585 <- grid.arrange(salmon_hot_plt, chanos_hot_plt, trout_hot_plt, nrow=3, ncol=1, top='Sigma dissimilarity by 2100 under SSP 585')

# pixel plot of plots
grid.arrange(full_119, full_585, nrow=1, ncol=2)

# all hex grids
grid.arrange(salmon_gg_grid, salmon_gg_hot_grid, chanos_gg_grid, chanos_gg_hot_grid, trout_gg_grid, trout_gg_hot_grid, 
             nrow=3, ncol=2, top='Sigma dissimilarity by 2100')

# hybrid
grid.arrange(salmon_gg_grid, salmon_gg_hot_grid, chanos_plt, chanos_hot_plt, trout_plt, trout_hot_plt,
             nrow=3, ncol=2,
             top='Sigma dissimilarity by 2100 under two different climate change scenarios')

# salmon hex
grid.arrange(salmon_gg_grid, salmon_gg_hot_grid,
             nrow=2, ncol=1,
             top='Sigma dissimilarity by 2100 under two different climate change scenarios (Salmo salar)')

# boxplots ----
box <- data.frame()
salmo_hot_box <- as.data.frame(salmo_plot_hot_grid[,'Dissim'])
salmo_hot_box$species <- 'Salmo salar'
salmo_hot_box$scenario <- 'SSP5-8.8'
names(salmo_hot_box) <- c('Dissim', 'Species', 'Scenario')

salmo_box <- as.data.frame(salmo_plot_grid[,'Dissim'])
salmo_box$species <- 'Salmo salar'
salmo_box$scenario <- 'SSP1-1.9'
names(salmo_box) <- c('Dissim', 'Species', 'Scenario')

trout_hot_box <- as.data.frame(trout_plot_hot_grid[,'Dissim'])
trout_hot_box$species <- 'Oncorhynchus mykiss'
trout_hot_box$scenario <- 'SSP5-8.8'
names(trout_hot_box) <- c('Dissim', 'Species', 'Scenario')

trout_box <- as.data.frame(trout_plot_grid[,'Dissim'])
trout_box$species <- 'Oncorhynchus mykiss'
trout_box$scenario <- 'SSP1-1.9'
names(trout_box) <- c('Dissim', 'Species', 'Scenario')

chanos_hot_box <- as.data.frame(chanos_plot_hot_grid[,'Dissim'])
chanos_hot_box$species <- 'Chanos chanos'
chanos_hot_box$scenario <- 'SSP5-8.8'
names(chanos_hot_box) <- c('Dissim', 'Species', 'Scenario')

chanos_box <- as.data.frame(chanos_plot_grid[,'Dissim'])
chanos_box$species <- 'Chanos chanos'
chanos_box$scenario <- 'SSP1-1.9'
names(chanos_box) <- c('Dissim', 'Species', 'Scenario')

box <- rbind(salmo_box, trout_box, salmo_hot_box, trout_hot_box, chanos_hot_box, chanos_box)

ggplot() +
  geom_boxplot(data = box, aes(x=Scenario, y=Dissim, color=Species))

ggplot(box, aes(x=Scenario, y=Dissim, fill=Species)) +
  geom_violin(trim = FALSE) +
  stat_boxplot(geom = "errorbar", position = position_dodge(.9), width = 0.2, color='black') +
  geom_boxplot(width=.1, position = position_dodge(.9), color='black', outlier.shape = NA) +
  labs(title='Dissimilarity by EEZ, scenario, and species')

NZ <- subset(world, region == 'New Zealand')

ggplot(data = st_as_sf(eez), aes(fill=NULL)) +
  geom_map(
    data = NZ, map = NZ,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_sf()
  # xlim(160, 190) +
  # ylim(-70,-29)









