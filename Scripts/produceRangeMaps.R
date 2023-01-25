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

## -----------------

myTaxaDB <- read.csv("../Data/AquacultureSpeciesList.csv", sep=";")
speciesToGetRecords <- unique(myTaxaDB$Species)

baseMapFile <- "../Data/MaskRes025.tif"
baseMapFile <- raster(baseMapFile)

fileAquamaps <- "../Data/Biodiversity database Aquamaps.db"

directoryMF <- "../Data/marineForestsDB/"

# -----------------------
# -----------------------
# List species available in aquamaps

sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver, dbname = fileAquamaps)
taxaDB <- dbReadTable(db,"taxa")
taxaDB$scientificName <- taxa <- paste0(taxaDB$Genus," ",taxaDB$Species)

nativemaps <- dbReadTable(db,"nativemaps")
cells <- dbReadTable(db,"hcaf")
occurrence <- dbReadTable(db,"occ")

# -----------------------
# List species available in marineForestsDB

filesMF <- list.files(directoryMF)
filesMF <- gsub(".RData","",filesMF)

# -----------------------

speciesToGetRecords <- sort(speciesToGetRecords[speciesToGetRecords %in% taxaDB$scientificName | speciesToGetRecords %in% filesMF ])

# -----------------------
# -----------------------

rgdal::setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE")

for( sp in speciesToGetRecords ) {
  
  cat("\014")  
  cat("\n")  
  cat("# --------------")
  cat("\n")  
  cat("Process",which(speciesToGetRecords == sp),"out of",length(speciesToGetRecords),"\n")
  cat("\n")  
  cat("# --------------")
  cat("\n")  

  # -----------------------------
  # Aquamaps 

  sp.id <- taxaDB[which(taxaDB$scientificName == sp),"SPECIESID"]
  
  if( length(sp.id) != 0 ) {  
        
        sp.records <- occurrence[occurrence$SpeciesID == sp.id & occurrence$GoodCell == 1,"CsquareCode"]
        sp.records <- cells[cells$CsquareCode  %in% sp.records,c("CenterLong","CenterLat")]
        
        if(nrow(sp.records) > 0) {
          databaseAquamaps.occ <- data.frame( decimalLongitude=sp.records[,1],decimalLatitude=sp.records[,2] , speciesName=sp, sourceType="Occurrence Data")
        }
        
        # Native maps
        
        sp.records <- nativemaps[nativemaps$SpeciesID == sp.id & nativemaps$FAOAreaYN == 1 & nativemaps$probability >= 0.5 ,c("CsquareCode")]
        sp.records <- cells[cells$CsquareCode %in% sp.records,c("CenterLong","CenterLat")]
        
        if(nrow(sp.records) > 0) {
          databaseAquamaps.range <- data.frame( decimalLongitude=sp.records[,1],decimalLatitude=sp.records[,2] , speciesName=sp, sourceType="Native Maps")
        }
      
        # ------------
        
        RangeMap <- rasterize(databaseAquamaps.range[,1:2],baseMapFile, field=1)

  }
        
  # -----------------------------
  # Marine Forests
  
  sp.id <- which( filesMF == sp)
  
  if( length(sp.id) != 0) {  
       
        load(paste0(directoryMF,"/",filesMF[sp.id],".RData"))
        RangeMap <- raster::resample(ensembleReclassReachGlobal,baseMapFile,fun=max)

  }
  
  # -----------------------------
  # -----------------------------

  RangeMap[RangeMap > 0] <- 1
  RangeMap[RangeMap != 1] <- 0
  RangeMap[is.na(RangeMap)] <- 0
  RangeMap <- mask(RangeMap,baseMapFile)
  # plot(RangeMap)
  
  resultsDirectory.i <- paste0("../Data/rangeMapsPerSpecies/")
  if( ! dir.exists(resultsDirectory.i) ) { dir.create(resultsDirectory.i, recursive = T) }
  raster::writeRaster(RangeMap,filename=paste0(resultsDirectory.i,"/",sp,".tif"),format="GTiff",overwrite=T)
  
}

# --------------------------------------------------
# --------------------------------------------------