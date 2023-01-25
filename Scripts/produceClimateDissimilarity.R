## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
##      
##  biodiversityDS [ biodiversitydatascience.com ]
##
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------

setwd("/Volumes/StingRay/Dropbox/Manuscripts/Exposure of global aquaculture to projected climate change/sharedFolder/Scripts")

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
resultsFolderMain <- "../Data/climateDissimilarity/ssp119/"
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
  writeRaster(ROIRaster, gsub(".RData",".tif",resultsFolders[i]))
  
  # Crop to match only the extent od where the species is
  ROIPts <- xyFromCell(ROIRaster,Which(!is.na(ROIRaster), cells=TRUE)) 
  oceanBasins.i <- crop(oceanBasins, extent(c(min(ROIPts[,1])-1,max(ROIPts[,1])+1,min(ROIPts[,2])-1,max(ROIPts[,2])+1)))
  eez.i <- crop(eez, extent(c(min(ROIPts[,1])-1,max(ROIPts[,1])+1,min(ROIPts[,2])-1,max(ROIPts[,2])+1)))

  # iteract between all oceans and all eez inside
  
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

write.csv(myResultsDF , file=paste0("../Results/",resultsScope,"DF.csv"), row.names = FALSE)

## -----------------------
## -----------------------