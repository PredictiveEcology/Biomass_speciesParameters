prepPSPaNPP <- function(studyAreaANPP, PSPgis, PSPmeasure, PSPplot,
                        useHeight, biomassModel, PSPperiod, minDBH) {
  #Crop points to studyArea
  if (!is.null(studyAreaANPP)) {
    studyAreaANPP <- st_as_sf(studyAreaANPP) # in case SPDF
    studyAreaANPP <- st_transform(x = studyAreaANPP, crs = st_crs(PSPgis))
    message(yellow("Filtering PSPs for ANPP to study Area..."))
    PSP_sa <- PSPgis[studyAreaANPP,] %>%
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
  PSPmeasure <- PSPmeasure[PSPplot, on = c("MeasureID", "OrigPlotID1", "MeasureYear")]
  PSPmeasure[, c("Longitude", "Latitude", "Easting", "Northing", "Zone") := NULL]

  #Filter by > 30 trees at first measurement (P) to ensure forest.
  forestPlots <- PSPmeasure[, .(measures = .N), OrigPlotID1] %>%
    .[measures >= 30,]

  PSPmeasure <- PSPmeasure[OrigPlotID1 %in% forestPlots$OrigPlotID1,]
  PSPplot <- PSPplot[OrigPlotID1 %in% PSPmeasure$OrigPlotID1,]

  #Restrict to trees >= ## DBH. Maybe necessary to get rid of small trees when they are inconsistently recorded
  PSPmeasure <- PSPmeasure[DBH >= minDBH,]
  #decide what to do about above line and stem/density. PES approach is to just fix the data...

  #Calculate biomass
  #Height must be calculated separately if there are NA heights -
  if (useHeight) {
    PSPmeasureNoHeight <- PSPmeasure[is.na(Height)]
    PSPmeasureHeight <- PSPmeasure[!is.na(Height)]
    tempOut <- biomassCalculation(species = PSPmeasureHeight$newSpeciesName,
                                  DBH = PSPmeasureHeight$DBH,
                                  height = PSPmeasureHeight$Height,
                                  includeHeight = TRUE,
                                  equationSource = biomassModel)
    #check if height is missing, join if so -- function fails if data.table is empty
    if (nrow(PSPmeasureNoHeight) > 0) {
      tempOutNoHeight <- biomassCalculation(species = PSPmeasureNoHeight$newSpeciesName,
                                            DBH = PSPmeasureNoHeight$DBH,
                                            height = PSPmeasureNoHeight$Height,
                                            includeHeight = FALSE,
                                            equationSource = biomassModel)
      tempOut$biomass <- c(tempOut$biomass, tempOutNoHeight$biomass)
      tempOut$missedSpecies <- c(tempOut$missedSpecies, tempOutNoHeight$missedSpecies)
      PSPmeasure <- rbind(PSPmeasureHeight, PSPmeasureNoHeight)
    }
    PSPmeasure[, biomass := tempOut$biomass]
    setkey(PSPmeasure, MeasureID, OrigPlotID1, TreeNumber)
    #do stuff
  } else {
    tempOut <- biomassCalculation(species = PSPmeasure$newSpeciesName,
                                  DBH = PSPmeasure$DBH,
                                  height = PSPmeasure$Height,
                                  includeHeight = useHeight,
                                  equationSource = biomassModel)
    PSPmeasure$biomass <- tempOut$biomass
  }
  message(yellow("No PSP biomass estimate possible for these species: "))
  message(crayon::yellow(paste(unique(tempOut$missedSpecies), collapse = ", ")))


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
                              minimumSampleSize, NoOfIterations, knots,
                              speciesFittingApproach = "focal"){
  #Must filter PSPdata by all sppEquiv$PSP with same sppEquivCol
  if (!isTRUE(speciesCol %in% colnames(sppEquiv))) {
    stop("sppEquivCol not in sppEquiv")
  }
  gcSpecies <- speciesFittingApproach
  if (identical(speciesFittingApproach, "single")) {
    gcSpecies <- unique(sppEquiv[[speciesCol]])
  }
  gcSpecies <- setNames(nm = gcSpecies)

  spsp <- copy(PSPdata)
  spsp[, "areaAdjustedB" := biomass/PlotSize]
  spsp  <- spsp[, "plotBiomass" := sum(areaAdjustedB), .(MeasureID)]
  spsp[, "spPlotBiomass" := sum(areaAdjustedB), .(MeasureID, newSpeciesName)]
  spsp[, "spDom" := spPlotBiomass/plotBiomass, .(MeasureID)]
  spsp[, speciesTemp := equivalentName(value = newSpeciesName, df = sppEquiv,
                                       column = speciesCol, searchColumn = "PSP")]
  whNA <- is.na(spsp$speciesTemp)
  message(crayon::yellow("Removing ", paste(unique(spsp$newSpeciesName[whNA]), collapse = ", ")))
  message(crayon::yellow("   ... because they are not the sppEquiv"))
  spsp <- spsp[!whNA]
  spsp <- spsp[newSpeciesName %in% sppEquiv[["PSP"]]]
  freq <- spsp[, .(N = .N, spDom = spDom[1]), .(speciesTemp, MeasureID)]
  if (isTRUE(gcSpecies == "pairwise") || isTRUE(gcSpecies == "focal"))
    freq <- freq[spDom > 0.2] # minimum 20% dominance for pairwise or focal
  else
    spsp <- spsp[spDom > 0.5] # minimum 50% dominance
  # setorderv(freq, "speciesTemp")
  if (isTRUE(gcSpecies == "pairwise") || isTRUE(gcSpecies == "focal")) {
    speciesComp <- freq[, .(numSp = .N, spComp = paste(speciesTemp, collapse = "__")), by = "MeasureID"]
    # Pick only 2 spcies plots when using "pairwise"
    if (isTRUE(gcSpecies == "pairwise"))
      speciesComp <- speciesComp[numSp == 2]
    speciesCompN <- speciesComp[, .N, by = "spComp"]
    # speciesPairs <- speciesCompN[N > 100]
    # speciesComp <- speciesComp[spComp %in% speciesPairs$spComp]
    gcSpecies2 <- unique(sppEquiv[[speciesCol]]) %>% setNames(nm = .)
    gcSpecies3 <- unique(speciesComp$spComp) %>% setNames(nm = .)
    speciesForSplit <- if (isTRUE(gcSpecies == "pairwise")) gcSpecies3 else gcSpecies2
    speciesCompListAll <- split(speciesComp, speciesComp$spComp)
    speciesCompList <- lapply(speciesForSplit, function(spName) {
      rbindlist(speciesCompListAll[grepl(spName, names(speciesCompListAll))])
    })
    spspList <- lapply(speciesCompList, function(speciesCompInner) {
      spsp[MeasureID %in% speciesCompInner$MeasureID][spDom > 0.2]
    })
    if (isTRUE(gcSpecies == "focal")) {
      spspList <- Map(psp = spspList, spName = names(spspList), function(psp, spName) {
        psp[speciesTemp != spName, speciesTemp := "Other"]
      })
    }
  } else {
    spspList <- split(spsp, spsp$speciesTemp)
  }

  speciesForCurves <- names(spspList) %>% setNames(nm = .)
  message(crayon::yellow("-----------------------------------------------"))
  message(crayon::yellow("building growth curves from PSP data: "))
  outputGCs <- Map(species = speciesForCurves, makeGAMMdata, psp = spspList,
                   MoreArgs = list(speciesEquiv = sppEquiv,
                                   sppCol = speciesCol, NoOfIters = NoOfIterations,
                                   K = knots, minSize = minimumSampleSize, q = quantileAgeSubset))

  if (isTRUE("all" == speciesFittingApproach)) {
    outputGCs <- lapply(gcSpecies2, function(x) {
      matchingSpecies <- equivalentName(x, sppEquiv, column = "PSP")
      ogc <- Copy(outputGCs[[1]])
      #ogc$originalData <- ogc$originalData[species %in% matchingSpecies]
      #ogc$simData <- ogc$simData[species %in% matchingSpecies | is.na(species)]
      ogc
    })

  }
  return(outputGCs)
}

modifySpeciesTable <- function(gamms, speciesTable, factorialTraits, factorialBiomass, sppEquiv,
                               sppEquivCol, inflationFactorKey, standAgesForFitting,
                               approach, maxBInFactorial) {

  #set(factorialBiomass, NULL, "speciesCode", factorialBiomass$speciesCode)
  #setkey(factorialTraits, species)
  #setkey(factorialBiomass, species)

  #for each species (ie Gamm), find best fit
  #first subset by speices that had data
  gamms <- gamms[sapply(gamms, function(x) class(x) == "list")]
  species <- names(gamms)
  names(species) <- species
  outputTraits <- list()

  # This next try is because the factorial traits may have recovered from a memoised state and will fail this:
  try(setnames(factorialTraits, old = "species", new = "speciesCode"), silent = TRUE)
  # Make column with Sp which is Sp1 and Sp2 in both datasets
  set(factorialBiomass, NULL, "Sp", gsub(".+_(Sp.)", "\\1", factorialBiomass$species, perl = TRUE))
  set(setDT(factorialTraits), NULL, "Sp", gsub(".+_(Sp.)", "\\1", factorialTraits$species, perl = TRUE))

  factorialTraitsThatVary <- sapply(factorialTraits, function(x) length(unique(x)) > 1)
  factorialTraitsThatVary <- names(factorialTraitsThatVary)[factorialTraitsThatVary]
  factorialTraitsVarying <- factorialTraits[, ..factorialTraitsThatVary]

  gammsT <- purrr::transpose(gamms)
  message("starting digest")
  # digFB <- CacheDigest(list(factorialBiomass))
  # digFT <- CacheDigest(list(factorialTraitsVarying))
  dig <- CacheDigest(list(factorialBiomass, factorialTraitsVarying,
                          gammsT$speciesGamm, gammsT$NonLinearModel, speciesTable))
  factorialBiomass <- factorialBiomass[startsWith(factorialBiomass$Sp, "Sp")]

  #join with inflationFactorKey - it's possible this data.table::copy is unnecessary
  set(inflationFactorKey, NULL, "species", NULL)
  tempTraits <- copy(factorialTraitsVarying)
  tempTraits <- inflationFactorKey[tempTraits, on = c("growthcurve", "mortalityshape", "longevity", "mANPPproportion")]
  tempTraits <- tempTraits[, .(speciesCode, inflationFactor)]
  factorialBiomass <- tempTraits[factorialBiomass, on = "speciesCode"]

  # Take only 2-cohort pixels -- they will start with Sp
  factorialTraits <- factorialTraits[startsWith(factorialTraits$Sp, "Sp")]
  setnames(factorialBiomass, "age", "standAge")

  message("Estimate species parameters; minimizing diff between statistical fit and Biomass_core experiment")

  for (name in species) {
    # use for loop to allow for Cache on each species
    outputTraits[[name]] <- Cache(editSpeciesTraits, name = name, gamm = gamms[[name]],
                                  traits = speciesTable, fT = factorialTraitsVarying, fB = factorialBiomass,
                                  speciesEquiv = sppEquiv, sppCol = sppEquivCol,
                                  standAgesForFitting = standAgesForFitting,
                                  approach = approach,
                                  maxBInFactorial = maxBInFactorial,
                                  userTags = name, #debugCache = "quick",
                                  .cacheExtra = dig$outputHash, omitArgs = c("fT", "fB", "gamm", "traits"))
  }
  # outputTraits <- lapply(species, FUN = editSpeciesTraits, gamm = gamms,
  #                        traits = speciesTable, fT = factorialTraits, fB = factorialBiomass,
  #                        speciesEquiv = sppEquiv, sppCol = sppEquivCol, mortConstraints = mortConstraints,
  #                        growthConstraints = growthConstraints, mANPPconstraints = mANPPconstraints)

  outputTraitsT <- purrr::transpose(outputTraits)
  gammsList <- rbindlist(gammsT$originalData, idcol = "Pair")
  fullDataAll <- rbindlist(outputTraitsT$fullData, idcol = "Pair")
  newTraits <- rbindlist(outputTraitsT$bestTraits, idcol = "Pair")
  llAll <- rbindlist(outputTraitsT$ll, idcol = "Pair")

  # limit best traits to only those that are nearest to longevity provided in SpeciesTable
  # Take next higher longevity (the -Inf in the rolling joing, on the "last" join column i.e., longevity)
  # set(setDT(speciesTable), NULL, "longevityOrig", speciesTable$longevity)
  set(setDT(newTraits), NULL, "longevityOrigFac", newTraits$longevity)

  bt2 <- newTraits[speciesTable[, c("species", "longevity")],
                   #                 bt2 <- newTraits[speciesTable[, c("species", "longevity", "longevityOrig")],
                   on = c("species", "longevity"), roll = -Inf]
  bt2[, longevity := longevityOrigFac]
  bt2 <- unique(bt2[, c("species", "longevity")], by = c("species"))
  bt1 <- newTraits[speciesTable[, c("species", "longevity")],
                   #bt1 <- newTraits[speciesTable[, c("species", "longevity", "longevityOrig")],
                   on = c("species", "longevity"), roll = Inf]
  bt1[, longevity := longevityOrigFac]
  bt1 <- unique(bt1[, c("species", "longevity")], by = c("species"))
  speciesTableNew <- na.omit(rbindlist(list(bt1, bt2))  )
  newTraits <- newTraits[speciesTableNew, on = c("species", "longevity")]

  # Take next lower longevity (the Inf in the rolling joing, on the "last" join column i.e., longevity)
  ll2 <- fullDataAll[, list(BscaledNonLinear = mean(BscaledNonLinear),
                            predNonLinear = mean(predNonLinear)#,
                            #Pair = Pair[1]
                            #mortalityshape = mean(mortalityshape),
                            #mANPPproportion = mean(mANPPproportion)
  ), c("species", "standAge", "Pair")]
  gg <- ggplot(ll2, aes(standAge, BscaledNonLinear, colour = species)) +
    geom_line(size = 2) +
    geom_point(data = gammsList, aes(standAge, biomass, colour = speciesTemp), size = 0.25, alpha = 0.3) +
    geom_line(size = 2, aes(standAge, predNonLinear, col = species), lty = "dashed") +
    facet_wrap(~ Pair, nrow = ceiling(sqrt(length(outputTraits))), scales = "fixed") +
    xlim(c(0, 150)) + # ggplot2::scale_y_log() +
    ylim(c(0, 30000)) +
    ggtitle("Comparing best LandR curves (solid) with best Non-Linear fit (dashed)") +
    theme_bw()

  #Collapse new traits and replace old traits
  # newTraits <- rbindlist(outputTraitsT$bestTraits, fill = TRUE, idcol = "Pair")
  # newTraits <- na.omit(newTraits)
  newTraits[, AICWeights := exp( -0.5 * llNonLinDelta)]
  newTraits[, AICWeightsStd := AICWeights/sum(AICWeights), by = "species"]

  # Showing species-level averages -- this is not assigned to object
  bestWeighted <- newTraits[, .(growthcurve = round(sum(AICWeightsStd * growthcurve), 2),
                                longevity = round(sum(AICWeightsStd * longevity), 0),
                                mortalityshape = round(sum(AICWeightsStd * mortalityshape), 0),
                                mANPPproportion = round(sum(AICWeightsStd * mANPPproportion), 3),
                                inflationFactor = round(sum(AICWeightsStd * inflationFactor), 3)),
                            by = "species"]
  # numSp <- length(unique(newTraits$species))
  # best <- newTraits[, .(growthcurve = round(sum(AICWeightsStd * growthcurve)/numSp, 2)
  #                       #, mortalityshape = round(sum(AICWeightsStd * mortalityshape)/numSp, 0)
  #                       #, mANPPproportion = round(sum(AICWeightsStd * mANPPproportion)/numSp, 3)
  #                       #, inflationFactor = round(sum(AICWeightsStd * inflationFactor)/numSp, 3)
  # )]
  # pick best
  # newTraits <- newTraits[, whmin := .I[which.min(llNonLinDelta)], by = "speciesCode"][unique(whmin)]

  # update growthcurve so all are same
  # newTraits[ , c(names(best)) := best]

  # Find associated "species" in factorial --> remove existing inflationFactor and
  #  replace with the equivalent inflationFactor from the factorial
  # newTraits <- updateInflationFactor(best, factorialTraits, factorialBiomass)

  # Remove the loglikelihood and AIC columns
  # set(newTraits, NULL, grep("^ll|^AIC", colnames(newTraits), value = TRUE), NULL)
  #
  #
  #
  # if (!is.null(newTraits$mANPPproportion)) {
  #   newTraits[is.na(mANPPproportion), c("mANPPproportion", "inflationFactor") := .(3.33, 1)] #default mANPP
  # } else {
  # }
  # Update speciesTable with bestWeighted
  speciesTable <- copy(speciesTable) # This is needed to allow the Cache on editSpeciesTraits to work because of pass-by-reference
  bestWeighted <- speciesTable[match(bestWeighted$species, species),
                               c(names(bestWeighted)) := bestWeighted]

  return(list(best = bestWeighted, gg = gg))
}

modifySpeciesEcoregionTable <- function(speciesEcoregion, speciesTable) {

  if (nrow(speciesTable[is.na(inflationFactor),]) > 0) {
    missing <- speciesTable[is.na(inflationFactor)]$species
    message("averaging traits for these species: ", missing)
    #
    averageOfEstimated <- speciesTable[!is.na(inflationFactor),
                                       .(growthcurve = mean(growthcurve),
                                         mortalityshape = asInteger(mean(mortalityshape)),
                                         mANPPproportion = mean(mANPPproportion))]
    #I don't think inflationFactor should be averaged if longevity isn't..
    speciesTable[is.na(inflationFactor), `:=`(
      growthcurve = averageOfEstimated$growthcurve,
      mortalityshape = averageOfEstimated$mortalityshape,
      mANPPproportion = averageOfEstimated$mANPPproportion
    )]
  }

  message("modifying speciesEcoregion table based on newly estimated traits")
  #modify things by species

  newSpeciesEcoregion <- speciesEcoregion[speciesTable, on = c("speciesCode" = "species")]
  newSpeciesEcoregion[!is.na(inflationFactor), maxB := asInteger(maxB * inflationFactor)]

  newSpeciesEcoregion[, maxANPP := asInteger(maxB * mANPPproportion/100)]
  cols <- names(speciesEcoregion)
  newSpeciesEcoregion <- newSpeciesEcoregion[, .SD, .SDcols = cols]
  newSpeciesEcoregion[, speciesCode := as.factor(speciesCode)]
  newSpeciesEcoregion[, maxB := asInteger(maxB)]
  return(newSpeciesEcoregion)
}

makeGAMMdata <- function(species, psp, speciesEquiv,
                         sppCol, NoOfIters, K, minSize, q) {
  # psp already contains all the correct species -- next chunk is now obsolete & commented
  if (identical(species, "all")) {
    #matchingSpecies <- unique(speciesEquiv[, .(PSP),])
    K <- mean(unlist(K))
    q <- mean(unlist(q))
  } else if (grepl(".+__.+", species)) { # this is pairwise
    #species2 <- strsplit(species, split = "__")[[1]]
    # matchingSpecies <- unique(speciesEquiv[speciesEquiv[[sppCol]] %in% species2, .(PSP),])
  } else {
    # matchingSpecies <- unique(speciesEquiv[speciesEquiv[[sppCol]] == species, .(PSP),])
  }

  if (class(K) == "list") {
    K <- K[[species]]
  }

  if (class(q) == "list") {
    q <- q[[species]]
  }

  #subset the parameters that may be lists
  if (class(NoOfIters) == "list") {
    NoOfIters <- NoOfIters[[species]]
  }

  # #Need to calculate stand dominance first - remove all stands < 50% dominance, and of wrong species
  # if (NROW(matchingSpecies) != 2) {
  #   browser()
  #   spDominant <- matchingSpecies[psp, nomatch = 0, on = "PSP==newSpeciesName"]  ## filter species first
  #   spDominant[, length(unique(PSP)), by = "MeasureID"]
  #   spDominant <- spDominant[spDom > 0.5]
  #   standData <- spDominant[, .(PlotSize = mean(PlotSize), standAge = mean(standAge),
  #                               biomass = mean(plotBiomass), spDom = mean(spDom),
  #                               species = unique(PSP)),
  #                           by = .(MeasureYear, OrigPlotID1)]
  # } else {
  standData <- psp[, .(PlotSize = PlotSize[1], standAge = standAge[1],
                       biomass = spPlotBiomass[1], spDom = spDom[1]),
                   by = .(speciesTemp, MeasureYear, OrigPlotID1)]
  # }

  #test if there are sufficient plots to estimate traits
  if (nrow(standData) < minSize) {
    speciesGamm <- "insufficient data"
    names(speciesGamm) <- species
    return(speciesGamm)
  }
  #By default removing the 95th percentile of age - these points are usually too scattered to produce reliable estimates
  standData <- standData[standAge < quantile(standData$standAge, probs = q/100),]
  simulatedData <- simulateYoungStands(cohortData = standData, N = 50)
  simData <- rbindlist(list(standData, simulatedData), fill = TRUE)

  #This weights the real data by spDominance, without distorting the mean of fake data. Weights should center around 1 for gamm I think
  Realweights <- standData$spDom/mean(standData$spDom)
  Fakeweights <- rep(1, times = nrow(simulatedData))
  simData$Weights <- c(Realweights, Fakeweights)

  species <- na.omit(unique(simData$speciesTemp))

  # localEnv <- environment()
  localEnv <- new.env(parent = emptyenv())

  # Spread out the 0s so they are in the first 10 years.
  simData2 <- copy(simData)
  len <- sum(simData2$standAge <= 10)
  simData2[standAge <= 10, standAge := sample(1:10, size = len, replace = TRUE)]
  len <- sum(simData2$biomass <= 100)
  simData2[biomass <= 100, biomass := sample(1:100, size = len, replace = TRUE)]

  models <- list()
  models$CR <- list(
    eqn = "biomass ~ A * (1 - exp(-k * standAge))^p",
    plim = c(min = 1, max = 80),
    Alim = c(max(simData2$biomass) * c(min = 0.3, max = 0.9)),
    klim = c(min = 0.0001, max = 0.13))

  # models$Weibull <- list(
  #   eqn = "biomass ~ A * (1 - exp(-k * standAge^p))",
  #   plim = c(min = 1, max = 80),
  #   Alim = c(max(simData2$biomass) * c(min = 0.3, max = 0.9)),
  #   klim = c(min = 0.0001, max = 0.13))
  # models$Korf <- list(
  #   eqn = "biomass ~ A * (exp(-k * standAge^(-p)))",
  #   plim = c(min = 0.01, max = 1),
  #   Alim = c(max(simData2$biomass) * c(min = 0.3, max = 0.9)),
  #   klim = c(min = 2, max = 20))

  models$Gompertz <- list(
    eqn = "biomass ~ A * (exp(-k * exp(-p * standAge)))",
    plim = c(min = 0.01, max = 1),
    Alim = c(max(simData2$biomass) * c(min = 0.3, max = 0.9)),
    klim = c(min = 2, max = 20))

  models$Logistic <- list(
    eqn = "biomass ~ A / (1 + k * exp(-p * standAge))",
    plim = c(min = 0.001, max = 1),
    Alim = c(max(simData2$biomass) * c(min = 0.3, max = 0.9)) )
  models$Logistic$klim <- c(min = 10, max = max(models$Logistic$Alim))

  # models$Levakovic <- list(
  #   eqn = "biomass ~ A * (standAge ^ 2 / (k + standAge ^ 2)) ^ p",
  #   plim = c(min = 1, max = 6),
  #   Alim = c(max(simData2$biomass) * c(min = 0.3, max = 0.9)) )
  # models$Levakovic$klim <- c(min = 150, max = max(models$Levakovic$Alim) * 0.3)

  nlsouts <- list()
  if (FALSE) {
    plot(simData2$standAge, simData2$biomass, pch = 19, cex = 0.5)
    model <- "Korf"
    standAge <- 1:200;
    A <- 23500
    k <- 4
    p <- 0.3
    points(col = "green", standAge, eval(parse(text = models[[model]])[[1]][[3]]), pch = 19, cex = 0.5)
  }
  species <- as.character(species) %>% setNames(nm = .)
  speciesForFits <- setdiff(species, "Other") %>% setNames(nm = .)
  speciesForFitsMessage <- paste(speciesForFits, collapse = ", ")
  message(crayon::yellow(speciesForFitsMessage, ": fitting Non-linear equations (Chapman-Richards, Logistic, Gomopertz)"))
  nlsouts <- lapply(speciesForFits, function(sp) {
    datForFit <- dataForFit(simData2, sp)
    nlsoutInner <- list()
    for (model in names(models)) {
      eqnChar <- models[[model]]$eqn
      plim <- models[[model]]$plim
      Alim <- models[[model]]$Alim
      klim <- models[[model]]$klim
      for (ii in 1:700) {
        # Chapman Richards
        # https://www.srs.fs.usda.gov/pubs/gtr/gtr_srs092/gtr_srs092-068-coble.pdf
        nlsoutInner[[model]] <- try({
          robustbase::nlrob(as.formula(eqnChar, env = .GlobalEnv),
                            data = datForFit, #maxit = 200,
                            weights = Weights,
                            maxit = 200,
                            start = list(A = runif(1, Alim[["min"]], Alim[["max"]]),
                                         k = runif(1, klim[["min"]], klim[["max"]]),
                                         p = runif(1, plim[["min"]], plim[["max"]])
                            ),
                            trace = FALSE)
        }, silent = TRUE
        )
        #})

        if (!is(nlsoutInner[[model]], "try-error")) {
          break
        }
      }
    }
    nlsoutInner
  })

  nlsout <- lapply(nlsouts, function(nlsIn) {
    nlsIn[[which.min(sapply(nlsIn, function(x) if (is(x, "try-error")) 1e12 else AIC(x)))]]
  })

  #  if (is(nlsout, "try-error")) {
  message(crayon::yellow(paste(collapse = "", rep(" ", nchar(speciesForFitsMessage))), ": fitting gamm"))
  envOnGlobal <- new.env(parent = .GlobalEnv)
  envOnGlobal$K <- K
  speciesGamm <- lapply(speciesForFits, function(sp)  {
    datForFit <- dataForFit(simData2, sp)
    suppressWarnings(try(expr = mgcv::gamm(data = datForFit,
                                           formula = as.formula(gammFormula, env = envOnGlobal),
                                           random = list(MeasureYear = eval(randomFormula, envir = baseenv()),
                                                         OrigPlotID1 = eval(randomFormula, envir = baseenv())),
                                           weights = nlme::varFunc(~Weights), verbosePQL = FALSE,
                                           niterPQL = NoOfIters),
                         silent = TRUE))
  }
  )

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
  if (!any(lengths(speciesGamm) == 1)) { #this means it was a try-error as converged gamm is length 2
    speciesGamm <- list(speciesGamm = speciesGamm,
                        originalData = standData,
                        simData = simData,
                        NonLinearModel = nlsout)
  }


  return(speciesGamm)
}

editSpeciesTraits <- function(name, gamm, traits, fT, fB, speciesEquiv, sppCol, maxBInFactorial,
                              standAgesForFitting = c(0, 150), approach, inflationFactorKey) {

  # Gamm <- gamm[[name]]
  nameOrig <- name
  if (grepl("__", name)) {
    name <- strsplit(name, "__")[[1]]
    names(name) <- name
  }
  message("Estimating fit for ", nameOrig)


  #Subset traits to PSP species, return unchanged if no gamm present
  traits <- traits[species %in% name]
  #with two species - the gamm might converge for one only
  #this structure is to catch try-errors in both pairwise and single
  if (class(gamm) == "try-error" | class(gamm) == "character") {
    message("not estimating traits for ", name)
    # return(list(bestTraits = traits, fullData = NULL, ll = NULL))
    return(NULL)
  } else {
    #catch when not all models converged
    classesNonLinear <- unlist(lapply(gamm$NonLinearModel, class))
    classesGAMM <- unlist(lapply(gamm$speciesGamm, class))
    #decided to allow non-converged gamms, if non-linear converged
    if (any("try-error" %in% c(classesNonLinear, classesGAMM))) {
      message("not estimating traits for ", name)
      return(NULL)
    }

    # if (any("try-error" %in% classesGAMM)){
    #   browser()
    # }
  }

  maxBiomass <- gamm$originalData[, .(maxBiomass = max(biomass)), "speciesTemp"]
  setorderv(maxBiomass, "maxBiomass", order = -1L)
  set(maxBiomass, NULL, "Sp", paste0("Sp", 1:nrow(maxBiomass)))
  SpNames <- maxBiomass$Sp[match(names(gamm$NonLinearModel), maxBiomass$speciesTemp)]
  SpMapping <- data.table(Sp = SpNames, species = name)
  names(gamm$NonLinearModel) <- SpNames
  names(gamm$speciesGamm) <- SpNames
  standAge <- unique(fB$standAge)
  standAge <- standAge[standAge <= max(standAgesForFitting) & standAge >= min(standAgesForFitting)]
  predGrid <- as.data.table(expand.grid(Sp = SpNames, standAge = standAge))

  # Predict from statistical fits to data
  if (all(names(gamm$speciesGamm[[1]]) == c("lme", "gam"))) {# means from gamm (ie. mixed effect, not gam)
    predGrid[, `:=`(
      predNonLinear = predict(gamm$NonLinearModel[[unlist(.BY)]], .SD),
      predGamm = as.vector(predict(gamm$speciesGamm[[unlist(.BY)]][["gam"]], .SD, se.fit = FALSE))
    ), by = "Sp", .SDcols = "standAge"]
  } else {
    predGrid[, `:=`(
      predNonLinear = predict(gamm$NonLinearModel[[unlist(.BY)]], .SD),
      predGamm = as.vector(predict(gamm$speciesGamm[[unlist(.BY)]], .SD, se.fit = FALSE))
    ), by = "Sp", .SDcols = "standAge"]
  }
  # misses points past longevity -- have to expand explicitly the standAges for each "speciesCode"
  dt <- as.data.table(expand.grid(speciesCode = unique(fB$speciesCode), standAge = unique(predGrid$standAge)))
  set(dt, NULL, "Sp", gsub(".+_", "", dt$speciesCode))
  dt <- dt[predGrid, on = c("standAge", "Sp")]
  candFB <- fB[dt, on = c("standAge", "speciesCode", "Sp")] # misses points past longevity
  whNA <- which(is.na(candFB$B))
  set(candFB, whNA, "B", 0L)
  candFB[, `:=`(
    pixelGroup = pixelGroup[1],
    inflationFactor = inflationFactor[1]), by = "speciesCode"]

  scaleFactorGamm <- max(predGrid$predGamm)/maxBInFactorial
  scaleFactorNonLinear <- max(predGrid$predNonLinear)/maxBInFactorial
  stdevGamm <- sd(predGrid$predGamm)
  stdevNonLinear <- sd(predGrid$predNonLinear)

  set(candFB, NULL, "BscaledGamm", candFB$B * scaleFactorGamm * candFB$inflationFactor)
  set(candFB, NULL, "BscaledNonLinear", candFB$B * scaleFactorNonLinear * candFB$inflationFactor)

  # Curve matching
  ll <- candFB[, .(llGamm = sum(dnorm(x = BscaledGamm, mean = predGamm,
                                      sd = stdevGamm, log = TRUE)),
                   llNonLinear = sum(dnorm(x = BscaledNonLinear, mean = predNonLinear,
                                           sd = stdevNonLinear, log = TRUE))), .(pixelGroup)]
  # Sort them so that best is on top
  setorderv(ll, "llNonLinear", order = -1L)
  ll[, c("llNonLinDelta", "llGammDelta") := list(
    abs(llNonLinear - max(llNonLinear)),
    abs(llGamm - max(llGamm)))]
  set(ll, NULL, c("llGamm", "llNonLinear"), NULL)

  deltaDiff <- 2 # because we
  ll <- ll[#llDirectDelta < deltaDiff |
    llGammDelta < deltaDiff |
      llNonLinDelta < deltaDiff]
  candFB <- candFB[ll, on = "pixelGroup"]
  varsToInteger <- c("BscaledGamm", "BscaledNonLinear", "predNonLinear", "predGamm")
  set(candFB, NULL, varsToInteger, lapply(varsToInteger, function(v) asInteger(candFB[[v]])))

  rr <- candFB[standAge == 1][, c("speciesCode", "llNonLinDelta", "inflationFactor")]

  # Take the average of the best
  best <- fT[rr, on = c("speciesCode" = "speciesCode")]
  bestTraits <- SpMapping[best, on = "Sp"]
  candFB <- SpMapping[candFB, on = "Sp"]


  return(list(bestTraits = bestTraits, fullData = candFB, ll = ll))
}

makePSPgamms <- function(studyAreaANPP, PSPperiod, PSPgis, PSPmeasure,
                         PSPplot, useHeight, biomassModel, speciesCol,
                         sppEquiv, NoOfIterations, knots, minimumSampleSize,
                         quantileAgeSubset, minDBH, speciesFittingApproach) {

  #this function is just a wrapper around these functions, for caching purposess
  psp <- prepPSPaNPP(studyAreaANPP = studyAreaANPP, PSPperiod = PSPperiod,
                     PSPgis = PSPgis, PSPmeasure = PSPmeasure, PSPplot = PSPplot,
                     useHeight = useHeight, biomassModel = biomassModel, minDBH = minDBH)

  #Wrapper used to avoid caching psp object - too large
  speciesGAMMs <- buildGrowthCurves(PSPdata = psp, speciesCol = speciesCol,
                                    sppEquiv = sppEquiv, NoOfIterations = NoOfIterations,
                                    knots = knots, minimumSampleSize = minimumSampleSize,
                                    quantileAgeSubset = quantileAgeSubset,
                                    speciesFittingApproach = speciesFittingApproach)
  return(speciesGAMMs)
}

gammFormula <- quote(biomass ~ s(standAge, k = K, pc = 0))
randomFormula <- quote(~1)

dataForFit <- function(simDat2, sp) {
  simDat2[is.na(simDat2$speciesTemp) | simDat2$speciesTemp %in% sp]
}
