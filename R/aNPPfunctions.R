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
  # PSPmeasure <- PSPmeasure[DBH >= 10,]
  #decide what to do about above line and stem/density. PES approach is to just fix the data...
  
  #Calculate biomass
  tempOut <- biomassCalculation(species = PSPmeasure$newSpeciesName,
                                DBH = PSPmeasure$DBH,
                                height = PSPmeasure$Height,
                                includeHeight = useHeight,
                                equationSource = biomassModel)
  message(yellow("No PSP biomass estimate possible for these species: "))
  print(tempOut$missedSpecies)
  PSPmeasure$biomass <- tempOut$biomass
  
  #clean up - added a catch for incorrect plotSize affecting stem density
  densities <- PSPmeasure[, .(.N, PlotSize = mean(PlotSize)), MeasureID]
  densities[, density := N/PlotSize]
  junkPlots <- densities[density > 4000]$MeasureID
  PSPmeasure <- PSPmeasure[!MeasureID %in% junkPlots]
  
  #bad biomass estimate
  PSPmeasure <- PSPmeasure[biomass != 0]
  PSPmeasure$newSpeciesName <- as.factor(PSPmeasure$newSpeciesName)
  
  #add stand age estimate
  PSPmeasure[, standAge := baseSA + MeasureYear - baseYear]
  PSPmeasure <- PSPmeasure[standAge > 0]
  PSPmeasure <- PSPmeasure[!is.na(biomass)]
  
  #divide biomass by 10 to get kg/ha into g/m2, the LandR unit
  PSPmeasure[, biomass := biomass/10]
  
  return(PSPmeasure)
}

buildGrowthCurves <- function(PSPdata, speciesCol, sppEquiv, quantileAgeSubset =  95){
  
  #Must filter PSPdata by all sppEquiv$PSP with same sppEquivCol
  gcSpecies <- unique(sppEquiv[[speciesCol]])
  message("building growth curves from PSP data for these species: ")
  print(gcSpecies)
  
  #summed psp - need year and origPlotID as random effects
  spsp <- copy(PSPdata)
  spsp[, 'areaAdjustedB' := biomass/PlotSize]
  spsp  <- spsp[, 'plotBiomass' := sum(areaAdjustedB), .(MeasureID)]
  spsp[, 'spPlotBiomass' := sum(areaAdjustedB), .(MeasureID, newSpeciesName)]
  spsp[, "spDom" := spPlotBiomass/plotBiomass, .(MeasureID)]
  
  outputGCs <- lapply(gcSpecies, FUN = function(species, psp = spsp, speciesEquiv = sppEquiv, sppCol = speciesCol) {

    matchingSpecies <- speciesEquiv[speciesEquiv[[speciesCol]] == species, .(PSP),]
    #Need to calculate stand dominance first - remove all stands < 50% dominance, and of wrong species 
    spDominant <- spsp[newSpeciesName %in% unique(matchingSpecies) & spDom > 0.5,]
    spDominant <- spDominant[, .(PlotSize = mean(PlotSize), standAge = mean(standAge), 
                                 biomass = mean(plotBiomass), spDom = mean(spDom)), 
                             .(MeasureYear, OrigPlotID1)]
    
    #Try removing the 95th percentile of age - these points are usually too scattered to produce reliable estimates
    spDominant <- spDominant[standAge < quantile(spDominant$standAge, probs = quantileAgeSubset/100),]
    simulatedData <- simulateYoungStands(cohortData = spDominant, N = 50)
    simData <- rbind(spDominant, simulatedData)
    
    #This weights the real data by spDominance, without distorting the mean of fake data. Weights should center around 1 for gamm I think
    Realweights <- spDominant$spDom/mean(spDominant$spDom)
    Fakeweights <- rep(1, times = nrow(simulatedData))
    simData$Weights <- c(Realweights, Fakeweights)
    
    #What the hell is this non-contrasts warning - ask Ceres or Eliot
    speciesGamm <- suppressWarnings(try(expr = gamm(data = simData, formula = biomass ~ s(standAge, k = 4, pc = 0), 
                                   random = list(MeasureYear = ~1, OrigPlotID1 = ~1), weights = varFunc(~Weights), verbosePQL = FALSE, 
                                   niterPQL = 10),
                       silent = TRUE))
    
    #Append the true data to speciesGamm, so we don't have the 0s involved when we subset by age quantile
    if (!length(speciesGamm) == 1){ #this means it was a try-error as converged gamm is length 2
      speciesGamm <- append(speciesGamm, list(originalData = spDominant))
    }
    
    return(speciesGamm)
    
  })
  
  names(outputGCs) <- gcSpecies
  return(outputGCs)
}

modifySpeciesTable <- function(gamms, speciesTable, factorialTraits, factorialBiomass, sppEquiv, 
                               sppEquivCol, mortConstraints, growthConstraints) {
  #Join these two tables to make the full speciesTraits table
  factorialBiomass[, species := speciesCode]
  setkey(factorialTraits, species)
  setkey(factorialBiomass, species)
  
  #for each species (ie Gamm), find best fit
  outputTraits <- lapply(names(gamms), function(name, gamm = gamms, traits = speciesTable, fT = factorialTraits,
                                                fB = factorialBiomass, speciesEquiv = sppEquiv, sppCol = sppEquivCol) {
    
    Gamm <- gamm[[name]]
    
    #Subset traits to PSP species, return unchanged if no Gamm present
    traits <- traits[species == name]
    
    if (class(Gamm) == 'try-error') {
      return(traits)
    }
    
    predData <- data.table(standAge = min(Gamm$originalData$standAge):max(Gamm$originalData$standAge))
    output <- predict(Gamm$gam, predData, se.fit = TRUE)
    predData <- data.table("age" = predData$standAge, "predBiomass" = output$fit, 'predSE' = output$se.fit)
    
    
    closestLongevity <- abs(fT$longevity - traits$longevity) == min(abs(fT$longevity - traits$longevity))
    #subset traits by closest longevity 
    CandidateTraits <- fT[closestLongevity]
    
    #Constrain growth curve - this is because the effect is conflated with maxANPP
    CandidateTraits <- CandidateTraits[growthcurve >= min(growthConstraints)
                                       & growthcurve <= max(growthConstraints),]
    
    #Constrain mortality shape - limited available information on mortalitity in PSPs, too low adds computation strain
    CandidateTraits <- CandidateTraits[mortalityshape >= min(mortConstraints) 
                                       & mortalityshape <= max(mortConstraints)]
    
    #subset the simulation values by potential species
    CandidateValues <- fB[CandidateTraits]
    
    #This has to happen before age is subset, or it will underestimate
    CandidateValues[, inflationFactor := 5000/max(B), .(speciesCode)]
    #subset the simulation values by age for which curve is represented
    CandidateValues <- CandidateValues[age >= min(predData$age) & age <= max(predData$age),]
    
    setkey(predData, age)
    setkey(CandidateValues, age)
    CandidateValues <- na.exclude(CandidateValues[predData])
    scaleFactors <- CandidateValues[, .(scaleFactor = mean(predBiomass/B),
                                        sMaxB = max(B),
                                        inflationFactor = mean(inflationFactor)), 'speciesCode']

    #scale factor is the achieved maxB in the simulation / PSP maxB. We use this to scale simulation values to PSP
    #inflationFactor is the simulation's achieved maxB / the LANDIS speciesTrait maxB that was used (always 5000)
    #scale factor is NOT returned. inflation factor is returned to 'inflate' Boreal_LBMRDataPrep estimates
    CandidateValues <- CandidateValues[scaleFactors, on = 'speciesCode']
    
    #Find best possible candidate species
    CandidateValues[, se := sd(B * scaleFactor - predBiomass), by = 'speciesCode']
    Candidates <- CandidateValues[, .(LogLikelihood = sum(dnorm(x = B * scaleFactor, mean = predBiomass,
                                                                 sd = se, log = TRUE)),
                                      inflationFactor = mean(inflationFactor)), .(speciesCode)]
    
    bestCandidate <- CandidateTraits[Candidates, on = c("species" = 'speciesCode')] %>%
      .[LogLikelihood == max(LogLikelihood)]

    # What to do witht he LogLikelihood? Report? 
    bestTraits <- bestCandidate[, .(mortalityshape, growthcurve, mANPPproportion, inflationFactor)]
    #if there are multiple rows due to undifferentiated curves, then mortality hasn't kicked in. Take the max mortalityshape
    bestTraits <- bestTraits[mortalityshape == max(mortalityshape)]
    traits[, c('mortalityshape', 'growthcurve', 'mANPPproportion', 'inflationFactor') :=  bestTraits]
    return(traits)
  })

  #Collapse new traits and replace old traits
  newTraits <- rbindlist(outputTraits, fill = TRUE)
  newTraits[is.na(mANPPproportion), c('mANPPproportion', 'inflationFactor') := .(3.33, 1)] #default mANPP
  
  return(newTraits)
}

modifySpeciesEcoregionTable <- function(speciesEcoregion, speciesTable) {
  message("modifying speciesEcoregion table based on traits derived from PSP Gamms")
  
  #modify things by species
  newSpeciesEcoregion <- speciesEcoregion[speciesTable, on = c('speciesCode' = 'species')]
  newSpeciesEcoregion[, maxB := maxB * inflationFactor]
  newSpeciesEcoregion[, maxANPP := maxB * mANPPproportion/100]
  cols <- names(speciesEcoregion)
  newSpeciesEcoregion <- newSpeciesEcoregion[, .SD, .SDcols = cols]
  
  return(newSpeciesEcoregion)
}