## ANALOGS MEGA DF ##

# 17.04.23
# Produces MegaDF, but still struggles with % A -> B
# which.max() is functional, but does not appear to handle ties or non-values very well, so additional filtering is required at each step
# Summarizing or aggregating top origins to EEZ level is difficult, as most will have different values for each species, how to weight fairly? Can use % contribution, but will still result in a lot of lost information if presented visually via directional links on a map

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

closeAllConnections()
rm(list=(ls()[ls()]))
gc(reset=TRUE)

require(raster)
require(sf)
require(sp)
require(ggplot2)
require(geodata)
require(tidyverse)
require(ggplot2)
require(rgeos)

eez <- shapefile("../Data/spatialInformation/eez.shp")
checkDF <- data.frame()


BigDF_Analogs <- function(test = FALSE, filtered = TRUE) {
  scenarios <- list.files(paste0("../Results/climateDissimilarity/"))
  bigDF <- data.frame()
  for (scenario in scenarios) {
    resultsFolderMain <- paste0("../Results/climateDissimilarity/", scenario)
    scenario_species <- list.files(resultsFolderMain)
    if(test==T) {
      test_main <- paste0(resultsFolderMain,'/','test_species')
      scenario_species <- list.files(test_main)
      print(scenario_species)
    }
    
    for (species in scenario_species) {
      print(species)
      active_species <- species
      DissimData <- load(paste0(test_main,'/', active_species,'/','climateDissimilarity.RData'))
      DissimData <- dataStructureResult
      
      test_dissim <- DissimData
      if(filtered==T) test_dissim <- DissimData[DissimData$sigmaNovelty>2,]
      print(nrow(test_dissim))
      if (nrow(test_dissim) > 5) {
        
        test_dissim <- test_dissim[,c('x', 'y', 'analogNovelty.x', 'analogNovelty.y')]
        points_dissim <- as.data.frame(cbind(test_dissim$x, test_dissim$y))
        names(points_dissim) <- c('x','y')
        points_dissim <- SpatialPointsDataFrame(coords = points_dissim[,c('x','y')], data = points_dissim)
        proj4string(points_dissim) <- proj4string(eez)
        
        points_analogs <- as.data.frame(cbind(test_dissim$analogNovelty.x, test_dissim$analogNovelty.y))
        names(points_analogs) <- c('x','y')
        points_analogs <- SpatialPointsDataFrame(coords = points_analogs[,c('x','y')], data = points_analogs)
        proj4string(points_analogs) <- proj4string(eez)
        
        test_dissim[,"originEEZ"] <- over(points_dissim,eez)[2]
        test_dissim[,'destEEZ'] <- over(points_analogs,eez)[2]
        
        # assigning points outside of EEZs as "International Waters"
        for (i in 1:nrow(test_dissim)) {
          if (is.na(test_dissim$destEEZ[[i]])) {
            test_dissim[i,'destEEZ'] <- 'International Waters'
          }
          if (is.na(test_dissim$originEEZ[[i]])) {
            test_dissim[i,'originEEZ'] <- 'International Waters'
          }
        }
        ret <- data.frame(matrix(nrow=length(unique(test_dissim$originEEZ)), ncol = 2))
        names(ret) <- c('EEZ', 'Retention %')
        row.names(ret) <- unique(test_dissim$originEEZ)
        # simple data frame of origin EEZ and % retention of analogs (already            filtered to > 2)
        for (e in unique(test_dissim$originEEZ)) {
          e_df <- test_dissim[test_dissim$originEEZ==e,]
          ret[e,1] <- e
          ret[e,2] <- nrow(e_df[e_df$destEEZ==e,])/nrow(e_df)
        }
        # construct intermediate DF to be populated with EEZ rows and later              rbound to bigDF
        origins <- unique(test_dissim$originEEZ)
        dests <- unique(test_dissim$destEEZ)
        only_dests <- unique(test_dissim$destEEZ)[!(unique(test_dissim$destEEZ) %in% unique(test_dissim$originEEZ))]
        all_lands <- c(origins, only_dests)
        inter_cols <- c('EEZ', 'dest_nrow', 'origin_nrow', 'ret_prop', 'species', 'scenario', 'top_origin', 'top_origin_prop', 'top_dest', 'top_dest_prop')
        interDF <- data.frame(matrix(ncol = length(inter_cols), nrow = length(all_lands)))
        # still need columns to capture % flow A -> B
        # can use names(which.max(table(test_dissim[test_dissim$destEEZ==land,]$originEEZ)))
        # or something similar to get the most prevalent element in a column, but capturing and then summarizing/aggregating these values to EEZ level will be difficult to do meaningfully
        colnames(interDF) <- inter_cols
        rownames(interDF) <- all_lands
        for (land in origins) {
          nrow_dest <- nrow(test_dissim[test_dissim$destEEZ==land,])
          df <- test_dissim[test_dissim$originEEZ==land,]
          nrow_orig <- nrow(df)
          prop_ret <- ret[land, 2]
          print(paste('Current land is', land))
          top_origin <- names(which.max(table(test_dissim[test_dissim$destEEZ==land & test_dissim$originEEZ!=land,]$originEEZ)))
          print(paste('The top origin for', land, 'is', top_origin))
          if (!is.null(top_origin)) {
            top_origin_prop <- (nrow(test_dissim[test_dissim$destEEZ==land & test_dissim$originEEZ==top_origin,])/nrow(test_dissim[test_dissim$originEEZ==top_origin,]))
            print(paste('The proportion of analogs coming from the top origin is', top_origin_prop))# % analogs headed for this EEZ (number bound from top_origin divided by total from top_origin)
          } else {
            top_origin_prop <- NA
            top_origin <- NA
          }
          top_dest <- names(which.max(table(test_dissim[test_dissim$destEEZ==land & test_dissim$originEEZ!=land,]$destEEZ)))
          if (!is.null(top_dest)) {
            top_dest_prop <- (nrow(test_dissim[test_dissim$originEEZ==land,])/nrow(test_dissim[test_dissim$destEEZ==top_dest,])) # % analogs exported to top destination (including retained analogs)
          } else {
            top_dest_prop <- NA
            top_dest <- NA
          }
          
          interDF[land, 'EEZ'] <- land
          interDF[land, 'dest_nrow'] <- nrow_dest
          interDF[land, 'origin_nrow'] <- nrow_orig
          interDF[land, 'ret_prop'] <- prop_ret
          interDF[land, 'species'] <- species
          interDF[land, 'scenario'] <- scenario
          interDF[land, 'top_origin'] <- top_origin
          interDF[land, 'top_origin_prop'] <- top_origin_prop
          interDF[land, 'top_dest'] <- top_dest
          interDF[land, 'top_dest_prop'] <- top_dest_prop
        }
        for (dest in only_dests) {
          nrow_dest <- nrow(test_dissim[test_dissim$destEEZ==dest,])
          df <- test_dissim[test_dissim$originEEZ==land,]
          nrow_orig <- nrow(df)
          prop_ret <- ret[land, 2]
          interDF[dest, 'EEZ'] <- dest
          interDF[dest, 'dest_nrow'] <- nrow_dest
          interDF[dest, 'origin_nrow'] <- 0
          interDF[dest, 'ret_prop'] <- -1
          interDF[dest, 'species'] <- species
          interDF[dest, 'scenario'] <- scenario
          interDF[land, 'top_origin'] <- top_origin
          interDF[land, 'top_origin_prop'] <- top_origin_prop
        }
        # rbind to bigDF
        bigDF <- rbind(bigDF, interDF)
        checkDF <- interDF
      }
    }
  }
  return(bigDF)
}

final <- BigDF_Analogs(test = T)
# still need refinement for a which.max tiebreaker

BigDFxEEZ <- final %>%
  group_by(EEZ, scenario) %>% arrange(EEZ) %>%
  summarise(Average_retention = round(mean(ret_prop[ret_prop >= 0]), 3), 
            n_gained = length(dest_nrow[dest_nrow > 0 & origin_nrow == 0]), 
            n_lost = length(dest_nrow[dest_nrow == 0 & origin_nrow > 0]), 
            n_retained = length(ret_prop[ret_prop > 0]),
            Species_count = n_distinct(species),
            Species = paste(unique(species), collapse = ","),
            Top_origins = paste(names(which.max(table(unique(top_origin)))), collapse = ","),
            Top_origin_prop = max(top_origin_prop)) #%>%
# for an EEZ, number of species with destination there, but not origin (n species gained)
# for an EEZ, number of species with origin and a destination there aka retention > 0 (n species retained)
# for an EEZ, number of species with origin there, but not destination (n species lost)
group_by(EEZ) %>% 
  summarise(Species_count = n_distinct(species), Average_retention = round(mean(na.omit(Average_retention)), 3), 
            Species = paste(unique(species), collapse = ",")
            , .groups = "drop")

