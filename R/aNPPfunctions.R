prepPSPaNPP <- function(studyAreaANPP, PSPgis, PSPmeasure, PSPplot,
                        useHeight, biomassModel, PSPperiod) {
  
  #Crop points to studyArea
  if (!is.null(studyAreaANPP)) {
    tempSA <- spTransform(x = studyAreaANPP, CRSobj = crs(PSPgis)) %>%
      st_as_sf(.)
    message(yellow("Filtering PSPs for ANPP to study Area..."))
    PSP_sa <- PSPgis[tempSA,] %>% 
      setkey(., OrigPlotID1)
    message(yellow(paste0("There are "), nrow(PSP_sa), " PSPs in your study area"))
    
    
    #Filter other PSP datasets to those in study Area
    PSPmeasure <- PSPmeasure[OrigPlotID1 %in% PSP_sa$OrigPlotID1,]
    PSPplot <- PSPplot[OrigPlotID1 %in% PSP_sa$OrigPlotID1,]
  }
  
  #length(PSPclimData)/length(PSP_sa) should always yield a whole number.
  #Filter data by study period
  message(yellow("Filtering PSPs for ANPP by study period..."))
  PSPmeasure <- PSPmeasure[MeasureYear > min(PSPperiod) &
                             MeasureYear < max(PSPperiod),]
  PSPplot <- PSPplot[MeasureYear > min(PSPperiod) &
                       MeasureYear < max(PSPperiod),]
  
  #Join data (should be small enough by now)
  PSPmeasure <- PSPmeasure[PSPplot, on = c('MeasureID', 'OrigPlotID1', 'MeasureYear')]
  PSPmeasure[, c('Longitude', 'Latitude', 'Easting', 'Northing', 'Zone'):= NULL]
  
  #Filter by > 30 trees at first measurement (P) to ensure forest.
  forestPlots <- PSPmeasure[, .(measures = .N), OrigPlotID1] %>%
    .[measures >= 30,]
  PSPmeasure <- PSPmeasure[OrigPlotID1 %in% forestPlots$OrigPlotID1,]
  PSPplot <- PSPplot[OrigPlotID1 %in% PSPmeasure$OrigPlotID1,]
  
  #Restrict to trees > 10 DBH (P) This gets rid of some big trees. Some 15 metres tall. Necessary because they are inconsistently recorded
  PSPmeasure <- PSPmeasure[DBH >= 10,]
  
  #Calculate biomass
  tempOut <- biomassCalculation(species = PSPmeasure$newSpeciesName,
                                DBH = PSPmeasure$DBH,
                                height = PSPmeasure$Height,
                                includeHeight = useHeight,
                                equationSource = biomassModel)
  message(yellow("No PSP biomass estimate possible for these species: "))
  print(tempOut$missedSpecies)
  PSPmeasure$biomass <- tempOut$biomass
  PSPmeasure <- PSPmeasure[biomass != 0]
  
  PSPmeasure$newSpeciesName <- as.factor(PSPmeasure$newSpeciesName)
  PSPmeasure[, standAge := baseSA + MeasureYear - baseYear]
  PSPmeasure <- PSPmeasure[standAge > 0]
  PSPmeasure <- PSPmeasure[!is.na(biomass)]
  
  return(PSPmeasure)
}

buildGrowthCurves <- function(PSPdata, speciesTable){
  
}