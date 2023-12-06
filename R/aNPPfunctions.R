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

  #Filter data by study period
  message(yellow("Filtering PSPs for ANPP by study period..."))
  PSPmeasure <- PSPmeasure[MeasureYear > min(PSPperiod) &
                             MeasureYear < max(PSPperiod),]
  PSPplot <- PSPplot[MeasureYear > min(PSPperiod) &
                       MeasureYear < max(PSPperiod),]

  PSPmeasure <- PSPmeasure[PSPplot, on = c("MeasureID", "OrigPlotID1", "MeasureYear")]

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
                              speciesFittingApproach = "focal") {
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
  message(crayon::yellow("   ... because they are not in sppEquiv"))
  spsp <- spsp[!whNA]
  spsp <- spsp[newSpeciesName %in% sppEquiv[["PSP"]]]
  freq <- spsp[, .(N = .N, spDom = spDom[1]), .(speciesTemp, MeasureID)]
  
     # setorderv(freq, "speciesTemp")
  if (isTRUE(speciesFittingApproach == "pairwise") || isTRUE(speciesFittingApproach == "focal")) {
    #subset to species-of-interest using relative biomass (dominance)
    freq <- freq[spDom > 0.2] # minimum 20% dominance for pairwise or focal
    
    speciesComp <- freq[, .(numSp = .N, spComp = paste(speciesTemp, collapse = "__")), by = "MeasureID"]
    # Pick only 2 spcies plots when using "pairwise"
    if (isTRUE(speciesFittingApproach == "pairwise")) {
      speciesComp <- speciesComp[numSp == 2]
    }
    speciesCompN <- speciesComp[, .N, by = "spComp"]

    gcSpecies1 <- unique(sppEquiv[[speciesCol]]) %>% setNames(nm = .)
    gcSpecies2 <- unique(speciesComp$spComp) %>% setNames(nm = .)
    speciesForSplit <- if (isTRUE(speciesFittingApproach == "pairwise")) gcSpecies2 else gcSpecies1
    speciesCompListAll <- split(speciesComp, speciesComp$spComp)
    speciesCompList <- lapply(speciesForSplit, function(spName) {
      rbindlist(speciesCompListAll[grepl(spName, names(speciesCompListAll))])
    })
    spspList <- lapply(speciesCompList, function(speciesCompInner) {
      spsp[MeasureID %in% speciesCompInner$MeasureID][spDom > 0.2]
    })
    if (isTRUE(speciesFittingApproach == "focal")) {
      spspList <- Map(psp = spspList, spName = names(spspList), function(psp, spName) {
        psp[speciesTemp != spName, speciesTemp := "Other"]
      })
    }
  } else {
    #subset to species-of-interest using relative biomass (dominance)
    spsp <- spsp[spDom > 0.5] # minimum 50% dominance for single 
    spspList <- split(spsp, spsp$speciesTemp)
  }
  
  speciesForCurves <- names(spspList) %>% setNames(nm = .)
  message(crayon::yellow("-----------------------------------------------"))
  message(crayon::yellow("building growth curves from PSP data: "))
  outputGCs <- Map(species = speciesForCurves, buildModels, psp = spspList,
                   MoreArgs = list(speciesEquiv = sppEquiv,
                                   sppCol = speciesCol, NoOfIters = NoOfIterations,
                                   K = knots, minSize = minimumSampleSize, q = quantileAgeSubset))
  
  if (isTRUE("all" == speciesFittingApproach)) {
    outputGCs <- lapply(gcSpecies1, function(x) {
      matchingSpecies <- equivalentName(x, sppEquiv, column = "PSP")
      ogc <- Copy(outputGCs[[1]])
      #ogc$originalData <- ogc$originalData[species %in% matchingSpecies]
      #ogc$simData <- ogc$simData[species %in% matchingSpecies | is.na(species)]
      ogc
    })
  }
  return(outputGCs)
}

modifySpeciesTable <- function(gcs, speciesTable, factorialTraits, factorialBiomass, sppEquiv,
                               sppEquivCol, inflationFactorKey, standAgesForFitting,
                               approach, maxBInFactorial) {

  #set(factorialBiomass, NULL, "speciesCode", factorialBiomass$speciesCode)
  #setkey(factorialTraits, species)
  #setkey(factorialBiomass, species)

  #for each species (ie Gamm), find best fit
  #first subset by speices that had data
  gcs <- gcs[sapply(gcs, function(x) class(x) == "list")]
  species <- names(gcs)
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

  gammsT <- purrr::transpose(gcs)
  message("starting digest")
  # digFB <- CacheDigest(list(factorialBiomass))
  # digFT <- CacheDigest(list(factorialTraitsVarying))
  dig <- CacheDigest(list(factorialBiomass, factorialTraitsVarying,
                          gammsT$speciesGamm, gammsT$NonLinearModel, speciesTable))
  factorialBiomass <- factorialBiomass[startsWith(factorialBiomass$Sp, "Sp")]
  gc()
  #join with inflationFactorKey - it's possible this data.table::copy is unnecessary
  suppressWarnings(set(inflationFactorKey, NULL, "species", NULL))
  tempTraits <- copy(factorialTraitsVarying)
  tempTraits <- inflationFactorKey[tempTraits, on = c("growthcurve", "mortalityshape", "longevity", "mANPPproportion")]
  tempTraits <- tempTraits[, .(speciesCode, inflationFactor)]
  factorialBiomass <- tempTraits[factorialBiomass, on = "speciesCode"]

  # Take only 2-cohort pixels -- they will start with Sp
  factorialTraits <- factorialTraits[startsWith(factorialTraits$Sp, "Sp")]
  setnames(factorialBiomass, "age", "standAge")

  message("Estimate species parameters; minimizing diff between statistical fit and Biomass_core experiment")

  # use for loop to allow for Cache on each species
  outputTraits <- Map(name = species, gamm = gcs, f = editSpeciesTraits, 
                      MoreArgs = list(traits = speciesTable, fT = factorialTraitsVarying, fB = factorialBiomass,
                                       speciesEquiv = sppEquiv, sppCol = sppEquivCol,
                                       standAgesForFitting = standAgesForFitting,
                                       approach = approach, maxBInFactorial = maxBInFactorial))
  gc()

  outputTraitsT <- purrr::transpose(outputTraits)
  gammsList <- rbindlist(gammsT$originalData, idcol = "Pair")
  fullDataAll <- rbindlist(outputTraitsT$fullData, idcol = "Pair")
  newTraits <- rbindlist(outputTraitsT$bestTraits, idcol = "Pair")
  llAll <- rbindlist(outputTraitsT$ll, idcol = "Pair")
  rm(gammsT, factorialBiomass, factorialTraits, factorialTraitsVarying)
  gc()
  # limit best traits to only those that are nearest to longevity provided in SpeciesTable
  # Take next higher longevity (the -Inf in the rolling joing, on the "last" join column i.e., longevity)
  # set(setDT(speciesTable), NULL, "longevityOrig", speciesTable$longevity)
  suppressWarnings(set(setDT(newTraits), NULL, "longevityOrigFac", newTraits$longevity))

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
                            predNonLinear = mean(predNonLinear)),
                     c("species", "standAge", "Pair")]
  rm(fullDataAll)
  ymaxes <- max(ll2$BscaledNonLinear, ll2$predNonLinear)
  gg <- ggplot(ll2, aes(standAge, BscaledNonLinear, colour = species)) +
    geom_line(size = 2) +
    geom_point(data = gammsList, aes(standAge, biomass, colour = speciesTemp), size = 0.25, alpha = 0.3) +
    geom_line(size = 2, aes(standAge, predNonLinear, col = species), lty = "dashed") +
    facet_wrap(~ Pair, nrow = ceiling(sqrt(length(outputTraits))), scales = "fixed") +
    xlim(c(0, max(ll2$standAge))) + # ggplot2::scale_y_log() +
    ylim(c(0, ymaxes)) +
    ylab(label = "biomass") +
    xlab(label = "stand age") +
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


  speciesTable <- copy(speciesTable) # This is needed to allow the Cache on editSpeciesTraits to work because of pass-by-reference
  bestWeighted <- speciesTable[match(bestWeighted$species, species),
                               c(names(bestWeighted)) := bestWeighted]

  return(list(best = bestWeighted, gg = gg))
}

buildModels <- function(species, psp, speciesEquiv,
                       sppCol, NoOfIters, K, minSize, q) {

  if (identical(species, "all")) {
    K <- mean(unlist(K))
    q <- mean(unlist(q))
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
  
  standData <- psp[, .(PlotSize = PlotSize[1], standAge = standAge[1],
                       biomass = spPlotBiomass[1], spDom = spDom[1]),
                   by = .(speciesTemp, MeasureYear, OrigPlotID1)]
  
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
 
  #   plot(simData2$standAge, simData2$biomass, pch = 19, cex = 0.5)
  #   model <- "Korf"
  #   standAge <- 1:200;
  #   A <- 23500
  #   k <- 4
  #   p <- 0.3
  #   points(col = "green", standAge, eval(parse(text = models[[model]])[[1]][[3]]), pch = 19, cex = 0.5)
  species <- as.character(species) %>% setNames(nm = .)
  speciesForFits <- setdiff(species, "Other") %>% setNames(nm = .)
  speciesForFitsMessage <- paste(speciesForFits, collapse = ", ")
  message(crayon::yellow(
    speciesForFitsMessage, ": fitting Non-linear equations (Chapman-Richards, Logistic, Gompertz)"
  ))
  nlsouts <- lapply(speciesForFits, function(sp, spFitData = simData2) {
    
    datForFit <- spFitData[is.na(spFitData$speciesTemp) | spFitData$speciesTemp %in% sp]
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
  speciesGamm <- lapply(speciesForFits, tryGAMM, gammData = simData2, envir = envOnGlobal,
                        randomFormula = randomFormula, gammFormula = gammFormula
                        )
  
  #Append the true data to speciesGamm, so we don't have the 0s involved when we subset by age quantile
  sppGrowthCurves <- list(speciesGamm = speciesGamm,
                      originalData = standData,
                      simData = simData,
                      NonLinearModel = nlsout)
  
  return(sppGrowthCurves)
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
    message("not estimating traits for ", name, " as model was not fit. Output of fitting attempt:")
    message(paste(gamm))
    # return(list(bestTraits = traits, fullData = NULL, ll = NULL))
    return(NULL)
  } else {
    #catch when not all models converged
    classesNonLinear <- unlist(lapply(gamm$NonLinearModel, class))
    classesGAMM <- unlist(lapply(gamm$speciesGamm, class))
    #decided to allow non-converged gamms, if non-linear converged
    if (any("try-error" %in% c(classesNonLinear, classesGAMM))) {
      message("not estimating traits for ", name, " as model was not fit. Output of fitting attempt:")
      message(paste(gamm$NonLinearModel))
      message(paste(gamm$speciesGamm))
      return(NULL)
    }
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
  rm(rr,fT, SpMapping, deltaDiff, gamm, fB)
  gc()
  return(list(bestTraits = bestTraits, fullData = candFB, ll = ll)) 
}

buildGrowthCurves_Wrapper <- function(studyAreaANPP, PSPperiod, PSPgis, PSPmeasure,
                         PSPplot, useHeight, biomassModel, speciesCol,
                         sppEquiv, NoOfIterations, knots, minimumSampleSize,
                         quantileAgeSubset, minDBH, speciesFittingApproach) {

  ## this function is just a wrapper around these functions, for caching purposes
  psp <- prepPSPaNPP(studyAreaANPP = studyAreaANPP, PSPperiod = PSPperiod,
                     PSPgis = PSPgis, PSPmeasure = PSPmeasure, PSPplot = PSPplot,
                     useHeight = useHeight, biomassModel = biomassModel, minDBH = minDBH)

  ## wrapper used to avoid caching PSP object - too large
  speciesGAMMs <- buildGrowthCurves(PSPdata = psp, speciesCol = speciesCol,
                                    sppEquiv = sppEquiv, NoOfIterations = NoOfIterations,
                                    knots = knots, minimumSampleSize = minimumSampleSize,
                                    quantileAgeSubset = quantileAgeSubset,
                                    speciesFittingApproach = speciesFittingApproach)
  return(speciesGAMMs)
}

gammFormula <- quote(biomass ~ s(standAge, k = K, pc = 0))
randomFormula <- quote(~1)


tryGAMM <- function(sp, gammFormula, randomFormula, 
                    gammData, iters = noOfIters, envir) {
  datForFit <- gammData[is.na(gammData$speciesTemp) | gammData$speciesTemp %in% sp]
  suppressWarnings(try(expr = mgcv::gamm(data = datForFit,
                                         formula = as.formula(gammFormula, env = envir),
                                         random = list(MeasureYear = eval(randomFormula, envir = baseenv()),
                                                       OrigPlotID1 = eval(randomFormula, envir = baseenv())),
                                         weights = nlme::varFunc(~Weights), verbosePQL = FALSE,
                                         niterPQL = iters),
                       silent = TRUE))
}
