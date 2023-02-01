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


# please ignore, I'm working on this
current <- data.frame(long = sal1$long, lat = sal1$lat)
current$climate <- "current"
novel <- data.frame(long = sal1$analogNovelty.x, lat = sal1$analogNovelty.y)
novel$climate <- "novel"
dissappearing <- data.frame(long = sal1$analogDisappearance.x, lat = sal1$analogDisappearance.y)
dissappearing$climate <- "dissappearing"

dat <- merge(current, novel, dissappearing, by.y = "long", all.x = TRUE, all.y = TRUE)

merge


countries <- world(path = "../Data/countries")

world <- map_data("world")



# Salmo maps ----

download_window <- ext(c(-60, 60, 50, 80))
plot(countries, ext = download_window)
points(salar585_subset[ , c("x", "y")], pch = 20, col = "red")
points(salar585_subset[ , c("analogNovelty.x", "analogNovelty.y")], pch = 20, col = "blue")
points(salar585_subset[ , c("analogDisappearance.x", "analogDisappearance.y")], pch = 20, col = "green")

  
sal <- ggplot() +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    color = "grey", fill = "lightgray", size = 0.01
  ) +
  geom_point(data=sal1, aes(x = long, y = lat), color = "darkred") +
  geom_point(data=sal1, aes(x = analogNovelty.x, y = analogNovelty.y), color = "blue") +
  geom_point(data=sal1, aes(x = analogDisappearance.x, y = analogDisappearance.y), color = "darkgreen") +
  xlim(-60, 60) + ylim(50, 80) +
  theme(legend.position = "right")

  
sal










