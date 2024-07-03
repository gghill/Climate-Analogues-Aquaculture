
library(raster)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))    # sets directory to script path
getwd()

myFolder <- "../Data/Shapefiles/Group 5 Shapefiles"
myFiles <- list.files(myFolder, pattern="shp", full.names = TRUE)
myFilesNames <- list.files(myFolder, pattern="shp", full.names = FALSE)
myFilesNames <- gsub(".shp","",myFilesNames)

dumpFolder <- "../Data/rangeMaps/"

shape <- raster("../Data/spatialInformation/coastLine025.tif")

for( f in 1:length(myFiles)) {
  
  myShp <- shapefile(myFiles[f])
  myShpRaster <- rasterize(myShp,shape)
  myShpRaster[!is.na(myShpRaster)] <- 1
  myShpRaster <- mask(myShpRaster,shape)
  # plot(myShpRaster, col="black")
  
  if( ! dir.exists(paste0(dumpFolder,myFilesNames[f],"/"))) { dir.create(paste0(dumpFolder,myFilesNames[f],"/")) }
  if( ! dir.exists(paste0(dumpFolder,myFilesNames[f],"/_ source/"))) { dir.create(paste0(dumpFolder,myFilesNames[f],"/_ source/") ) }
  
  writeRaster(myShpRaster,filename=paste0(dumpFolder,myFilesNames[f],"/expertRangeMap.tif"),format="GTiff",overwrite=T)
  shapefile(myShp, file=paste0(dumpFolder,myFilesNames[f],"/_ source/manualRangeMap.shp"))
  
}
