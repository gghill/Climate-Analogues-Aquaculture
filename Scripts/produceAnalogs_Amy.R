## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
##      
##                  Mapping climate analogs
##
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------

# First draft for mapping novel and disappearing climates AM, working from modified climate dissimilarity script

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

closeAllConnections()
rm(list=(ls()[ls()]))
gc(reset=TRUE)

library(raster)
library(sf)
library(ggplot2)
theme_set(theme_bw())
require(gridExtra)
library(geodata)
library(tidyverse)

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
resultsFolderMain <- paste0("../Data/climateDissimilarity/", scenario)
resultsScope <- "sigmaNovelty" # sigmaNovelty sigmaDisappearance 
resultsFolders <- list.files(resultsFolderMain, recursive=TRUE, pattern=".RData", full.names=TRUE)
resultsNames <- list.files(resultsFolderMain, recursive=TRUE, pattern=".RData", full.names=FALSE)
resultsNames <- gsub("/climateDissimilarity.RData","",resultsNames)

resultsFolders
resultsNames

## -----------------------
# load in data 

load("../Results/climateDissimilarity/ssp585/Salmo salar/climateDissimilarity.RData")
sal <- dataStructureResult
sal1<- sal[1:100,]

names(sal1)[names(sal1) == 'x'] <- 'long'
names(sal1)[names(sal1) == 'y'] <- 'lat'



# Mapping Salmo salar subset ----

countries <- world(path = "../Data/countries")

download_window <- ext(c(-60, 60, 50, 80))
plot(countries, ext = download_window)
points(salar585_subset[ , c("x", "y")], pch = 20, col = "red")
points(salar585_subset[ , c("analogNovelty.x", "analogNovelty.y")], pch = 20, col = "blue")
points(salar585_subset[ , c("analogDisappearance.x", "analogDisappearance.y")], pch = 20, col = "green")


# ggplot
world <- map_data("world")  

sal <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_point(data=sal1, aes(x = long, y = lat)) +
  geom_point(data=sal1, aes(x = analogNovelty.x, y = analogNovelty.y)) +
  geom_point(data=sal1, aes(x = analogDisappearance.x, y = analogDisappearance.y)) +
  xlim(-80, 80) + ylim(10, 80) +
  theme(legend.position = bottom)  + 
  scale_fill_manual(values = c('darkred', 'blue', 'darkgreen'), 
                      labels = c('Current climate', 'Novel climate', 'Dissappearing climate'))

  
sal      #legend is not working



# Combine coords into one df (but can't connect points based on cells)
current <- data.frame(long = sal1$long, lat = sal1$lat)
current$climate <- "current"
novel <- data.frame(long = sal1$analogNovelty.x, lat = sal1$analogNovelty.y)
novel$climate <- "novel"
dissappearing <- data.frame(long = sal1$analogDisappearance.x, lat = sal1$analogDisappearance.y)
dissappearing$climate <- "dissappearing"

dat <- rbind(current, novel, dissappearing) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)

install.packages("wesanderson")
library(wesanderson)
names(wes_palettes)

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
  scale_colour_manual(values = wes_palette("GrandBudapest1", n = 3), name = "Climate", 
                      labels = c("Current",  "Disappearing", "Novel"))
ggsave(path = "../Results/climateDissimilarity/ssp585/Salmo salar/", filename ="S.salar_subset.png", width = 6, height = 4)


# Adding lines to connect climates ----

#have to change points to sf before you can denote by colour...

df_c <- data.frame(cell = sal1$cell, c.long = sal1$long, c.lat = sal1$lat)
df_c <-st_as_sf(df_c %>% select(cell, c.long, c.lat), coords = c("c.long", "c.lat"), crs = 4326) %>%
  rename(current = geometry) 

df_n <- data.frame(cell = sal1$cell, n.long = sal1$analogNovelty.x, n.lat = sal1$analogNovelty.y)
df_n <-st_as_sf(df_n %>% select(cell, n.long, n.lat), coords = c("n.long", "n.lat"), crs = 4326) %>%
  rename(novel = geometry)
                   
df_d <- data.frame(cell = sal1$cell, d.long = sal1$analogDisappearance.x, 
                   d.lat = sal1$analogDisappearance.y)
df_d <-st_as_sf(df_d %>% select(cell, d.long, d.lat), coords = c("d.long", "d.lat"), crs = 4326) %>%
  rename(disappearing = geometry)

df_sf <- cbind(df_c, df_n, df_d)
df_sf <- df_sf[ , -c(2, 3)]

df_sf

# 

sal <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_point(data=sal1, aes(x = long, y = lat)) +
  geom_point(data=sal1, aes(x = analogNovelty.x, y = analogNovelty.y)) +
  geom_point(data=sal1, aes(x = analogDisappearance.x, y = analogDisappearance.y)) +
  xlim(-80, 80) + ylim(10, 80) +
  scale_colour_manual(values = wes_palette("GrandBudapest1", n = 3))




sal

# convert coords to sf for each?



# Mapping Salmo salar full dataset ----

# full_dat <- data.frame(long = rbind(sal$))      #I'll come back to this

# map all data
# represent data better with shape data 
# link points through cells using geom_line() - subset
# group points (size of point based on scale)




