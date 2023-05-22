## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
##      
##                  Visualizing EEZ connections
##
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------

# take metrics from MegaDF to visualize climate connections between EEZs
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

closeAllConnections()
rm(list=(ls()[ls()]))
gc(reset=TRUE)

library(sf)
library(rgeos)
library(ggplot2)

eez_v11 <- shapefile("../Data/spatialInformation/eez_v11.shp")

# extract centroids from EEZs for plotting
# sp_cent <- gCentroid(eez, byid=TRUE, id=eez@data[["EEZ"]])
# 
# sp_cent <- st_as_sf(sp_cent)
# sp_cent <- cbind(sp_cent, eez@data[["EEZ"]])
# names(sp_cent) <- c('EEZ', 'geometry')

# use EEZs that already have centroids
sf_cent <- eez_v11@data[,c('TERRITORY1','X_1','Y_1')]
names(sf_cent) <- c('EEZ', 'x', 'y')

ggplot() +
  geom_map(
  data = world, map = world,
  aes(map_id = region),
  color = "grey", fill = "lightgray", linewidth = 0.01
) +
  geom_sf(data = st_as_sf(eez), aes(fill=NULL)) +
  geom_point(data = sf_cent, aes(x=x,y=y), color='red') +
  geom_text(data = sf_cent, aes(x=x,y=y, label=EEZ), check_overlap = TRUE )
  