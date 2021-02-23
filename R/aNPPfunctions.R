prepPSPaNPP <- function(studyAreaANPP, PSPgis, PSPmeasure, PSPplot,
                        useHeight, biomassModel, PSPperiod, minDBH) {

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

  #Restrict to trees >= ## DBH. Maybe necessary to get rid of small trees when they are inconsistently recorded
  PSPmeasure <- PSPmeasure[DBH >= minDBH,]
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

buildGrowthCurves <- function(PSPdata, speciesCol, sppEquiv, quantileAgeSubset,
                              minimumSampleSize, NoOfIterations, knots){
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
  outputGCs <- lapply(gcSpecies, FUN = makeGAMMdata, psp = spsp, speciesEquiv = sppEquiv,
                      sppCol = speciesCol, NoOfIters = NoOfIterations,
                      K = knots, minSize = minimumSampleSize, q = quantileAgeSubset)

  names(outputGCs) <- gcSpecies
  return(outputGCs)
}

modifySpeciesTable <- function(gamms, speciesTable, factorialTraits, factorialBiomass, sppEquiv,
                               sppEquivCol, mortConstraints, growthConstraints, mANPPconstraints) {

  #Join these two tables to make the full speciesTraits table
  factorialBiomass[, species := speciesCode]
  setkey(factorialTraits, species)
  setkey(factorialBiomass, species)

  #for each species (ie Gamm), find best fit
  outputTraits <- lapply(names(gamms), FUN = editSpeciesTraits, gamm = gamms,
                         traits = speciesTable, fT = factorialTraits, fB = factorialBiomass,
                         speciesEquiv = sppEquiv, sppCol = sppEquivCol, mortConstraints = mortConstraints,
                         growthConstraints = growthConstraints, mANPPconstraints = mANPPconstraints)

  #Collapse new traits and replace old traits
  newTraits <- rbindlist(outputTraits, fill = TRUE)

  if (!is.null(newTraits$mANPPproportion)) {
    newTraits[is.na(mANPPproportion), c('mANPPproportion', 'inflationFactor') := .(3.33, 1)] #default mANPP
  } else {
    message("no GAMMs converged :( Consider revising your parameters")
  }
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
  newSpeciesEcoregion[, speciesCode := as.factor(speciesCode)]
  newSpeciesEcoregion[, maxB := asInteger(maxB)]
  return(newSpeciesEcoregion)
}

makeGAMMdata <- function(species, psp, speciesEquiv,
                         sppCol, NoOfIters, K, minSize, q) {

  matchingSpecies <- unique(speciesEquiv[speciesEquiv[[sppCol]] == species, .(PSP),])

  #subset the parameters that may be lists
  if (class(NoOfIters) == "list"){
    NoOfIters <- NoOfIters[[species]]
  }
  if (class(K) == "list"){
    K <- K[[species]]
  }

  if (class(q) == "list") {
    q <- q[[species]]
  }

  #Need to calculate stand dominance first - remove all stands < 50% dominance, and of wrong species
  spDominant <- matchingSpecies[psp, nomatch = 0, on = "PSP==newSpeciesName"]  ## filter species first
  spDominant <- spDominant[spDom > 0.5]
  spDominant <- spDominant[, .(PlotSize = mean(PlotSize), standAge = mean(standAge),
                               biomass = mean(plotBiomass), spDom = mean(spDom)),
                           by = .(MeasureYear, OrigPlotID1)]

  #test if there are sufficient plots to estimate traits
  if (nrow(spDominant) < minSize) {
    speciesGamm <- "insufficient data"
    names(speciesGamm) <- species
    return(speciesGamm)
  }
  #By default removing the 95th percentile of age - these points are usually too scattered to produce reliable estimates
  spDominant <- spDominant[standAge < quantile(spDominant$standAge, probs = q/100),]
  simulatedData <- simulateYoungStands(cohortData = spDominant, N = 50)
  simData <- rbind(spDominant, simulatedData)

  #This weights the real data by spDominance, without distorting the mean of fake data. Weights should center around 1 for gamm I think
  Realweights <- spDominant$spDom/mean(spDominant$spDom)
  Fakeweights <- rep(1, times = nrow(simulatedData))
  simData$Weights <- c(Realweights, Fakeweights)

  localEnv <- environment()

  speciesGamm <- suppressWarnings(try(expr = mgcv::gamm(data = simData, formula = eval(gammFormula, enclos = localEnv),
                                                  random = list(MeasureYear = eval(randomFormula, envir = baseenv()),
                                                                OrigPlotID1 = eval(randomFormula, envir = baseenv())
                                                                ),
                                                  weights = nlme::varFunc(~Weights), verbosePQL = FALSE,
                                                  niterPQL = NoOfIters),
                                      silent = TRUE))

  if (FALSE) {
    #this exists for manually debugging whether the encompassing environment will be cached
    object_size(speciesGamm)
    # 838 kB  # Same as above
    tf <- tempfile(); system.time(saveRDS(speciesGamm, file = tf))
    #check time elapsed
    e1 <- new.env(parent = emptyenv())
    e1$manualVers <- readRDS(tf)
    object_size(e1$manualVers)
    # This should be comparable to above (+/- 50%)
    rm(e1, tf)
  }

  #Append the true data to speciesGamm, so we don't have the 0s involved when we subset by age quantile
  if (!length(speciesGamm) == 1){ #this means it was a try-error as converged gamm is length 2
    speciesGamm <- append(speciesGamm, list(originalData = spDominant))
  }

  return(speciesGamm)

}

editSpeciesTraits <- function(name, gamm, traits, fT, fB, speciesEquiv, sppCol,
                              mortConstraints, growthConstraints, mANPPconstraints) {

  Gamm <- gamm[[name]]

  #Subset traits to PSP species, return unchanged if no Gamm present
  traits <- traits[species == name]

  if (class(Gamm) == 'try-error' | class(Gamm) == 'character') {
    message(paste("not estimating growth/mortality traits for", name))
    return(traits)
  }

  predData <- data.table(standAge = min(Gamm$originalData$standAge):max(Gamm$originalData$standAge))
  library('mgcv')
  output <- predict(Gamm$gam, predData, se.fit = TRUE)
  predData <- data.table("age" = predData$standAge, "predBiomass" = output$fit, 'predSE' = output$se.fit)


  closestLongevity <- abs(fT$longevity - traits$longevity) == min(abs(fT$longevity - traits$longevity))
  #subset traits by closest longevity
  CandidateTraits <- fT[closestLongevity]

  #Constrain growth curve - this is because the effect is conflated with maxANPP
  if (class(growthConstraints) == 'list') {
    growthConstraint <- growthConstraints[[name]]
  } else {
    growthConstraint <- growthConstraints
  }

  #growth constraint may be NULL if user passed constraints on partial list of species
  if (!is.null(growthConstraint)){
    CandidateTraits <- CandidateTraits[growthcurve %>=% min(growthConstraint)
                                       & growthcurve %<=% max(growthConstraint),]
  }
  #Constrain mortality shape - limited available information on mortalitity in PSPs, too low adds computation strain
  if (class(mortConstraints) == 'list') {
    mortConstraint <- mortConstraints[[name]]
  } else {
    mortConstraint <- mortConstraints
  }

  if (!is.null(mortConstraints)) {
    CandidateTraits <- CandidateTraits[mortalityshape >= min(mortConstraint)
                                       & mortalityshape <= max(mortConstraint)]
  }

  #Constrain growth curve - this is because the effect is conflated with maxANPP
  if (class(mANPPconstraints) == 'list') {
    mANPPconstraint <- mANPPconstraints[[name]]
  } else {
    mANPPconstraint <- mANPPconstraints
  }

  #growth constraint may be NULL if user passed constraints on partial list of species
  if (!is.null(mANPPconstraint)){
    CandidateTraits <- CandidateTraits[mANPPproportion %>=% min(mANPPconstraint)
                                       & mANPPproportion %<=% max(mANPPconstraint),]
  }

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
  #inflationFactor is the the LANDIS speciesTrait maxB that was used (always 5000) / simulation's achieved maxB
  #scale factor is NOT returned. inflation factor is returned to 'inflate' Biomass_borealDataPrep estimates
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
  bestTraits[, mortalityshape := asInteger(mortalityshape)]
  traits[, c('mortalityshape', 'growthcurve', 'mANPPproportion', 'inflationFactor') :=  bestTraits]
  return(traits)
}

makePSPgamms <- function(studyAreaANPP, PSPperiod, PSPgis, PSPmeasure,
                         PSPplot, useHeight, biomassModel, speciesCol,
                         sppEquiv, NoOfIterations, knots, minimumSampleSize,
                         quantileAgeSubset, minDBH) {

  #this function is just a wrapper around these functions, for caching purposess
  psp <- prepPSPaNPP(studyAreaANPP = studyAreaANPP, PSPperiod = PSPperiod,
                     PSPgis = PSPgis, PSPmeasure = PSPmeasure, PSPplot = PSPplot,
                     useHeight = useHeight, biomassModel = biomassModel, minDBH = minDBH)

  #Wrapper used to avoid caching psp object - too large
  speciesGAMMs <- buildGrowthCurves(PSPdata = psp, speciesCol = speciesCol,
                                    sppEquiv = sppEquiv, NoOfIterations = NoOfIterations,
                                    knots = knots, minimumSampleSize = minimumSampleSize,
                                    quantileAgeSubset = quantileAgeSubset)
  return(speciesGAMMs)
}

gammFormula <- quote(biomass ~ s(standAge, k = K, pc = 0))
randomFormula <- quote(~1)
