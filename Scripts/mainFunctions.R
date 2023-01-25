## ---------------------------------------------------------------------
## ---------------------------------------------------------------------

library(stringr)
library(RSQLite)
library(gdata)
library(raster)
library(RSQLite)
library(DBI)
library(stringi)
library(RCurl)
library(rvest)
library(dplyr)
library(scales)

# -----------
# -----------

firstToupper <- function(string) { return(paste(toupper(substr(string, 1, 1)), substr(tolower(string), 2, nchar(string)), sep="")) }

# -----------
# -----------

countWords <- function(x) { str_count(x, '\\w+') }

# -----------
# -----------

# getSpeciesTraitsAquamapsDepthRange("Pomadasys stridens", filename="/Volumes/Jellyfish/Dropbox/Data/Biodiversity Data/Occurrence records/Biodiversity database Aquamaps.db")

getSpeciesTraitsAquamaps <- function(speciesList, filename) {
  
  sqlite.driver <- RSQLite::SQLite() # dbDriver("SQLite")
  db <- dbConnect(sqlite.driver, dbname = filename)
  
  # dbListTables(db)
  
  taxa <- dbReadTable(db,"taxa")
  taxaTraits <- dbReadTable(db,"hspen")
  
  resultDF <- data.frame()
  
  for(species in speciesList) {
    
    taxa.sp <- taxa[taxa$Species == strsplit(species, " ")[[1]][2] & taxa$Genus == strsplit(species, " ")[[1]][1] ,]
    
    if( nrow(taxa.sp) > 0) {
      
      resultDF <- rbind(resultDF,data.frame( taxa.sp[,c("Kingdom","Phylum","Class","Order","Family","Genus")],
                                             Species = apply( taxa.sp[,c("Genus","Species")] , 1 , function(x) { paste( x[1],x[2]) } ),
                                             minDepth = sapply( taxa.sp[,c("SPECIESID")] , function(x) { taxaTraits[taxaTraits$SpeciesID == x , "DepthMin"] } ),
                                             maxDepth = sapply( taxa.sp[,c("SPECIESID")] , function(x) { taxaTraits[taxaTraits$SpeciesID == x , "DepthMax"] } ),
                                             distributionFAO = gsub(" ","",sapply( taxa.sp[,c("SPECIESID")] , function(x) { taxaTraits[taxaTraits$SpeciesID == x , "FAOAreas"] } )),
                                             distributionFAOReference = gsub(" ","",sapply( taxa.sp[,c("SPECIESID")] , function(x) { taxaTraits[taxaTraits$SpeciesID == x , "FAOAreasRef"] } )),
                                             isPleagic = sapply( taxa.sp[,c("SPECIESID")] , function(x) { taxaTraits[taxaTraits$SpeciesID == x , "Pelagic"] } ) == 1))
      
    }
  }
  
  return(resultDF)
  
}

# -----------
# -----------

# webScrapDepthRange("Diplodus sargus",fishBaseURL)
# webScrapDepthRange("Paramuricea clavata",seaLifeURL)

webScrapDepthRange <- function(species,baseURL) {
  
  defaultW <- getOption("warn")
  options(warn = -1)
  
  resultDF <- data.frame()

  pass <- TRUE
  tryCatch(
    expr = {
      web.html <- read_html(paste0(baseURL,"/summary/",gsub(" ","-",species),".html"))
    },
    error = function(e){
      pass <<- FALSE
    }
  )
  
  if( pass ) {
    
    depthRange <- html_elements(web.html,"span")
    
    if( length(grep("depth range ", depthRange)) > 0 ) { 
      
      depthRange <- depthRange[grep("depth range ", depthRange)[1]] 
      
      depthRangeValues <- substr(depthRange,gregexpr("depth range ",depthRange)[[1]][1]+12,gregexpr("m \\(Ref",depthRange)[[1]][1]-2)
      depthRangeValues <- strsplit(depthRangeValues, "-")[[1]]
      depthRangeValues <- as.numeric(stri_trim(depthRangeValues))
      
      depthRangeReference <- paste0(baseURL,"/references/FBRefSummary.php?ID=",substr(depthRange,regexpr("?ID",depthRange)[1]+3, ifelse(regexpr(" target=",depthRange)[1]-2 > 0,regexpr(" target=",depthRange)[1]-2,regexpr("</a>",depthRange)[1]-7)   ))
      depthRangeReference <- gsub('"',"",depthRangeReference, fixed=T)
      depthRangeReference <- gsub('>',"",depthRangeReference, fixed=T)
      depthRangeReference <- read_html(depthRangeReference)
      depthRangeReference <- html_elements(depthRangeReference,"td")
      
      if(length(depthRangeReference) != 0) {
        
        depthRangeReference <- depthRangeReference[grep("Citation", depthRangeReference)[1]+1] 
        depthRangeReference <- substr(depthRangeReference,as.vector(gregexpr("\t<b>",depthRangeReference)[[1]]) + 4, 
                                      as.vector(gregexpr("\t",depthRangeReference)[[1]])[as.vector(gregexpr("\t",depthRangeReference)[[1]]) > as.vector(gregexpr("\t<b>",depthRangeReference)[[1]])][1] - 2 )
        
        depthRangeReference <- gsub("<i>","",depthRangeReference)
        depthRangeReference <- gsub("</i>","",depthRangeReference)
        depthRangeReference <- gsub("</b>","",depthRangeReference)
        depthRangeReference <- gsub("&amp;"," ",depthRangeReference)
        
      }
      
      if(length(depthRangeReference) == 0) { depthRangeReference <- NA } 
      
      options(warn = defaultW)
      
      resultDF <- data.frame(species=species,
                             minDepth=depthRangeValues[1],
                             maxDepth=depthRangeValues[2],
                             reference=depthRangeReference)
      
    }

  }
  
  
  
  return(resultDF)
  
}

# -----------
# -----------

# webScrapFAORegions("Diplodus sargus",fishBaseURL)
# webScrapFAORegions("Paramuricea clavata",seaLifeURL)

webScrapFAORegions <- function(species,baseURL) {
  
  defaultW <- getOption("warn")
  options(warn = -1)
  
  pass <- TRUE
  tryCatch(
    expr = {
      web.html <- read_html(paste0(baseURL,"/summary/",gsub(" ","-",species),".html"));
      FAOAreas <- html_elements(web.html,"a")
    },
    error = function(e){
      pass <<- FALSE
    }
  )
  
  if(! pass ) { FAOAreas <- 0 }
  
  if( pass & length(FAOAreas) > 3 & length(grep("FAO areas", FAOAreas)) > 0 ) {
    
    FAOAreas <- FAOAreas[grep("FAO areas", FAOAreas)[1]] 
    FAOAreas <- substr(FAOAreas,as.vector(gregexpr("a href=",FAOAreas)[[1]]) + 11, as.vector(gregexpr("alt",FAOAreas)[[1]]) - 3 )
    FAOAreas <- paste0(baseURL,"/",FAOAreas)
    FAOAreas <- read_html(FAOAreas)
    FAOAreas <- html_elements(FAOAreas,"td")
    FAOAreas <- FAOAreas[grep("href=", FAOAreas)] 
    
    FAOAreasResult <- character(0)
    
    for( f in 1:length(FAOAreas)) {
      
      FAOAreasResult <- c(FAOAreasResult,
                          substr(FAOAreas[f],as.vector(gregexpr("AreaCode=",FAOAreas[f])[[1]]) + 9,
                                 min(as.vector(gregexpr("&amp;",FAOAreas[f])[[1]]))[as.vector(gregexpr("&amp;",FAOAreas[f])[[1]]) > as.vector(gregexpr("AreaCode=",FAOAreas[f])[[1]]) ] - 1))
      
    }
    
    FAOAreasResult <- paste(FAOAreasResult, collapse = ',')
    
    # ---------
    
    FAOAreasReference <- html_elements(web.html,"a")
    FAOAreasReference <- FAOAreasReference[(grep("/references/FBRefSummary.php\\?ID", FAOAreasReference))[grep("/references/FBRefSummary.php\\?ID", FAOAreasReference) > ifelse( length(grep("Upload your references", FAOAreasReference)) > 0 , grep("Upload your references", FAOAreasReference) , grep("Collaborators", FAOAreasReference) ) ][1]]
    
    FAOAreasReference <- paste0(baseURL,"/references/FBRefSummary.php?ID=",substr(FAOAreasReference,regexpr("?ID",FAOAreasReference)[1]+3, ifelse(regexpr(" target=",FAOAreasReference)[1]-2 > 0,regexpr(" target=",FAOAreasReference)[1]-2,regexpr("</a>",FAOAreasReference)[1]-7) ))
    FAOAreasReference <- gsub('"',"",FAOAreasReference, fixed=T)
    FAOAreasReference <- read_html(FAOAreasReference)
    FAOAreasReference <- html_elements(FAOAreasReference,"td")
    
    if( length(FAOAreasReference) != 0 ) {
      
      FAOAreasReference <- FAOAreasReference[grep("Citation", FAOAreasReference)[1]+1] 
      FAOAreasReference <- substr(FAOAreasReference,as.vector(gregexpr("\t<b>",FAOAreasReference)[[1]]) + 4, 
                                  as.vector(gregexpr("\t",FAOAreasReference)[[1]])[as.vector(gregexpr("\t",FAOAreasReference)[[1]]) > as.vector(gregexpr("\t<b>",FAOAreasReference)[[1]])][1] - 2 )
      
      FAOAreasReference <- gsub("</b>","",FAOAreasReference)
      FAOAreasReference <- gsub("&amp;"," ",FAOAreasReference)
      FAOAreasReference <- gsub("<i>","",FAOAreasReference)
      FAOAreasReference <- gsub("</i>","",FAOAreasReference)
      
    }
    
    if( length(FAOAreasReference) == 0 ) { FAOAreasReference <- NA }
      
    options(warn = defaultW)
    
    return(data.frame(Species=species,
                      FAOAreas=FAOAreasResult,
                      Reference=FAOAreasReference))
    
    
  }
  
  if( ! pass | length(FAOAreas) <= 3 | length(grep("FAO areas", FAOAreas)) == 0) { return( data.frame() ) } 
  
}

# -----------
# -----------

# webScrapGuild("Diplodus sargus",fishBaseURL)
# webScrapGuild("Paracentrotus lividus",seaLifeURL)

webScrapGuild <- function(species,baseURL) {
  
  GuildResultL1 <- NA
  GuildResultL2 <- NA
  GuildResult <- NA
  GuildReferenceL <- NA
  GuildReferences <- NA
  
  defaultW <- getOption("warn")
  options(warn = -1)
  
  pass <- TRUE
  tryCatch(
    expr = {
      web.html <- read_html(paste0(baseURL,"/summary/",gsub(" ","-",species),".html"));
      Guild <- html_elements(web.html,"a")
    },
    error = function(e){
      pass <<- FALSE
    }
  )
  
  if( pass & length(Guild) > 3 & length(grep("Diet", Guild)) > 0 ) {
    
    Guild <- html_elements(web.html,"a")
    Guild <- Guild[grep("Diet", Guild)[1]] 
    Guild <- substr(Guild,as.vector(gregexpr("a href=",Guild)[[1]]) + 9, as.vector(gregexpr("alt",Guild)[[1]]) - 3 )
    Guild <- gsub("&amp;","&",Guild)
    Guild <- paste0(baseURL,"/",Guild)
    Guild <- read_html(Guild)
    Guild <- html_elements(Guild,"td")
    Guild <- Guild[grep("href=", Guild)] 
    
    GuildResult <- character(0)
    
    for( f in 1:length(Guild)) {
      
      if( grepl("/references/",Guild[f])) { next }
      
      GuildResult <- c(GuildResult,
                       
                       substr(Guild[f],as.vector(gregexpr('">',Guild[f])[[1]]) + 2,
                              min(as.vector(gregexpr("</a>",Guild[f])[[1]]))[as.vector(gregexpr("</a>",Guild[f])[[1]]) > as.vector(gregexpr("dietcode=",Guild[f])[[1]]) ] - 1)   
                       
                       )
      
    }
    
    GuildResult <- paste( unique(GuildResult), collapse = ',')
    
    # ---------
    
    GuildReferences <- character(0)
    
    for( f in 1:length(Guild)) {
      
      if( grepl("/trophiceco/",Guild[f])) { next }
      
      GuildReference <- substr(Guild[f],as.vector(gregexpr("FBRefSummary.php\\?id=",Guild[f])[[1]]) + 0,as.vector(gregexpr('target="_blank"',Guild[f])[[1]]) - 2)   
      GuildReference <- gsub('"',"",GuildReference, fixed=T)
      GuildReference <- paste0(baseURL,"/references/",GuildReference)
      
      GuildReference <- read_html(GuildReference)
      GuildReference <- html_elements(GuildReference,"td")
      GuildReference <- GuildReference[grep("Citation", GuildReference)[1]+1] 
      GuildReference <- substr(GuildReference,as.vector(gregexpr("\t<b>",GuildReference)[[1]]) + 4, 
                                  as.vector(gregexpr("\t",GuildReference)[[1]])[as.vector(gregexpr("\t",GuildReference)[[1]]) > as.vector(gregexpr("\t<b>",GuildReference)[[1]])][1] - 2 )
      
      GuildReference <- gsub("</b>","",GuildReference)
      GuildReference <- gsub("&amp;"," ",GuildReference)
      
      GuildReferences <- c(GuildReferences, GuildReference)
     
    }
    
    GuildReferences <- paste( unique(GuildReferences), collapse = ',')
    
    options(warn = defaultW)

  }
  
  if( pass & length(Guild) > 3 & length(grep("Diet", Guild)) == 0 & length(grep("Food items", Guild)) > 0 ) { 
    
    Guild <- html_elements(web.html,"a")
    Guild <- Guild[grep("Food items", Guild)[1]] 
    Guild <- substr(Guild,as.vector(gregexpr("a href=",Guild)[[1]]) + 11, as.vector(gregexpr("alt",Guild)[[1]]) - 3 )
    Guild <- gsub("&amp;","&",Guild)
    Guild <- paste0(baseURL,"/",Guild)
    Guild <- read_html(Guild)
    
    GuildStructure <- html_elements(Guild,"tr")
    GuildStructure <- length(unlist(gregexpr("th width",GuildStructure[1])))
    
    Guild <- html_elements(Guild,"td")

    GuildResultL1 <- character(0)
    
    for( f in seq(1,length(Guild),by=GuildStructure) ) {
    
      GuildResult.f <- substr(Guild[f],as.vector(gregexpr('">',Guild[f])[[1]]) + 2,as.vector(gregexpr('</td>',Guild[f])[[1]]) - 2)
      if( grepl("/references/",Guild[f])) { next }
      
      GuildResultL1 <- c(GuildResultL1,GuildResult.f)   
                      
    }
    
    GuildResultL1 <- paste( unique(GuildResultL1), collapse = ',')

    # ----------------
    
    GuildResultL2 <- character(0)
    
    for( f in seq(2,length(Guild),by=GuildStructure) ) {
      
      GuildResult.f <- substr(Guild[f],as.vector(gregexpr('">',Guild[f])[[1]]) + 2,as.vector(gregexpr('</td>',Guild[f])[[1]]) - 2)
      if( grepl("/references/",Guild[f])) { next }
      
      GuildResultL2 <- c(GuildResultL2,GuildResult.f)   
      
    }
    
    GuildResultL2 <- paste( unique(GuildResultL2), collapse = ',')
    
    # ----------------
    
    
    GuildReferenceL <- html_elements(web.html,"a")
    GuildReferenceL <- GuildReferenceL[(grep("/references/FBRefSummary.php\\?ID", GuildReferenceL))[grep("/references/FBRefSummary.php\\?ID", GuildReferenceL) > ifelse( length(grep("Upload your references", GuildReferenceL)) > 0 , grep("Upload your references", GuildReferenceL) , grep("Collaborators", GuildReferenceL) ) ][1]]
    
    GuildReferenceL <- paste0(baseURL,"/references/FBRefSummary.php?ID=",substr(GuildReferenceL,regexpr("?ID",GuildReferenceL)[1]+3, ifelse(regexpr(" target=",GuildReferenceL)[1]-2 > 0,regexpr(" target=",GuildReferenceL)[1]-2,regexpr("</a>",GuildReferenceL)[1]-7) ))
    GuildReferenceL <- gsub('"',"",GuildReferenceL, fixed=T)
    GuildReferenceL <- read_html(GuildReferenceL)
    GuildReferenceL <- html_elements(GuildReferenceL,"td")
    
    if( length(GuildReferenceL) != 0) {
      
      GuildReferenceL <- GuildReferenceL[grep("Citation", GuildReferenceL)[1]+1] 
      GuildReferenceL <- substr(GuildReferenceL,as.vector(gregexpr("\t<b>",GuildReferenceL)[[1]]) + 4, 
                                as.vector(gregexpr("\t",GuildReferenceL)[[1]])[as.vector(gregexpr("\t",GuildReferenceL)[[1]]) > as.vector(gregexpr("\t<b>",GuildReferenceL)[[1]])][1] - 2 )
      
      GuildReferenceL <- gsub("</b>","",GuildReferenceL)
      GuildReferenceL <- gsub("&amp;"," ",GuildReferenceL)
      GuildReferenceL <- gsub("<i>","",GuildReferenceL)
      GuildReferenceL <- gsub("</i>","",GuildReferenceL)
      
    }
    
    if( length(GuildReferenceL) == 0) { GuildReferenceL <- NA } 
    
    
    
  }
  
  resultRequest <- data.frame(Species=species,
                              Diet=GuildResult,
                              DietReference=GuildReferences,
                              FoodItemsResultL1=GuildResultL1,
                              FoodItemsResultL2=GuildResultL2,
                              FoodItemsReference=GuildReferenceL )
  
  return(resultRequest)
  
}

# -----------
# -----------

# getTaxonomyWorms("Laminaria", FALSE)

getTaxonomyWorms <- function(name,verbose) {
  
  packages.to.use <- "worms"
  
  options(warn=-1)
  
  for(package in packages.to.use) {
    if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
    if( ! package %in% rownames(installed.packages()) ) { stop("Error on package instalation") }
    library(package, character.only = TRUE)
  }
  
  tryCatch( worms.query <- wormsbynames( name , verbose = FALSE , marine_only = "true" )
            , error=function(e) { } )
  
  options(warn=0)
  
  if( exists("worms.query") ) { 
    
    if( verbose ) { 
      
      cat('\n')
      cat('\n')
      cat('----------------------------------------------------------------------')
      cat('\n')
      cat('Processing taxonomy for', name )
      cat('\n')
      cat('\n')
      cat('\n')
      
    }
    
    return(data.frame(
      aphiaID = worms.query$AphiaID,
      name = worms.query$scientificname,
      authority = worms.query$authority,
      status = worms.query$status,
      taxKingdom = worms.query$kingdom,
      taxPhylum = worms.query$phylum,
      taxClass = worms.query$class,
      taxOrder = worms.query$order,
      taxFamily = worms.query$family,
      taxGenus = worms.query$genus,
      revisionByWormsDate = Sys.Date(),
      acceptedAphiaID = worms.query$valid_AphiaID,
      acceptedName = worms.query$valid_name
    ))
    
  } else {         
    
    if( verbose ) { 
      
      cat('\n')
      cat('\n')
      cat('----------------------------------------------------------------------')
      cat('\n')
      cat('Taxonomy for',name,'not found')
      cat('\n')
      cat('\n')
      cat('\n')
      
    }
    
  }
}

## -------------
## -------------

# getChildrenWorms(141433)

getChildrenWorms <- function(id) {
  
  packages.to.use <- "worrms"
  
  options(warn=-1)
  
  for(package in packages.to.use) {
    if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
    if( ! package %in% rownames(installed.packages()) ) { stop("Error on package instalation") }
    library(package, character.only = TRUE)
  }
  
  tryCatch( worms.query <- wormsbynames( name , verbose = FALSE , marine_only = "true" )
            , error=function(e) { } )
  
  options(warn=0)
  
  pass <- TRUE
  
  tryCatch( worms.query <- wm_children(id, marine_only = TRUE)
            , error=function(e) { pass <<- FALSE } )
  
  if(pass) { return(as.data.frame(worms.query)) }
  if(!pass) { return(data.frame()) }
  
  
  
}




