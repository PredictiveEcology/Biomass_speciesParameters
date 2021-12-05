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
  PSPmeasure <- PSPmeasure[PSPplot, on = c("MeasureID", "OrigPlotID1", "MeasureYear")]
  PSPmeasure[, c("Longitude", "Latitude", "Easting", "Northing", "Zone"):= NULL]
  
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
      tempOutNoHeight <- biomassCalculation(species =PSPmeasureNoHeight$newSpeciesName,
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
  print(unique(tempOut$missedSpecies))
  
  
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
                              allSpInOne = FALSE){
  #Must filter PSPdata by all sppEquiv$PSP with same sppEquivCol
  if (!isTRUE(speciesCol %in% colnames(sppEquiv))) {
    stop("sppEquivCol not in sppEquiv")
  }
  gcSpecies <- allSpInOne
  if (identical(allSpInOne, "single")) {
    gcSpecies <- unique(sppEquiv[[speciesCol]])
  } 
  
  
  #summed psp - need year and origPlotID as random effects
  spsp <- copy(PSPdata)
  spsp[, "areaAdjustedB" := biomass/PlotSize]
  spsp  <- spsp[, "plotBiomass" := sum(areaAdjustedB), .(MeasureID)]
  spsp[, "spPlotBiomass" := sum(areaAdjustedB), .(MeasureID, newSpeciesName)]
  spsp[, "spDom" := spPlotBiomass/plotBiomass, .(MeasureID)]
  spsp[, speciesTemp := equivalentName(value = newSpeciesName, df = sppEquiv, column = speciesCol, searchColumn = "PSP")]
  if (identical(gcSpecies, "pairwise")) {
    spsp <- spsp[newSpeciesName %in% sppEquiv[["PSP"]]]
    freq <- spsp[, .(N = .N, spDom = spDom[1]), .(speciesTemp, MeasureID)]
    freq <- freq[spDom > 0.2] # minimum 10 trees in a plot
    setorderv(freq, "speciesTemp")
    speciesComp <- freq[, .(numSp = .N, spComp = paste(speciesTemp, collapse = "__")), by = "MeasureID"]
    speciesComp <- speciesComp[numSp == 2]
    speciesCompN <- speciesComp[, .N, by = "spComp"]
    speciesPairs <- speciesCompN[N > 100]
    speciesComp <- speciesComp[spComp %in% speciesPairs$spComp]
    speciesCompList <- split(speciesComp, speciesComp$spComp)
    spspList <- lapply(speciesCompList, function(speciesCompInner) {
      spsp[MeasureID %in% speciesCompInner$MeasureID][spDom > 0.2]
    })
    gcSpecies <- names(spspList)
    message("building growth curves from PSP data for these species combinations: ")
    outputGCs <- Map(species = gcSpecies, makeGAMMdata, psp = spspList, 
                     MoreArgs = list(speciesEquiv = sppEquiv,
                                     sppCol = speciesCol, NoOfIters = NoOfIterations,
                                     K = knots, minSize = minimumSampleSize, q = quantileAgeSubset))
  } else {
    message("building growth curves from PSP data for these species: ")
    print(gcSpecies)
    browser()
    outputGCs <- lapply(gcSpecies, FUN = makeGAMMdata, psp = spsp, speciesEquiv = sppEquiv,
                        sppCol = speciesCol, NoOfIters = NoOfIterations,
                        K = knots, minSize = minimumSampleSize, q = quantileAgeSubset)
  }
  names(outputGCs) <- gcSpecies
  if (isTRUE(allSpInOne)) {
    gcSpecies2 <- unique(sppEquiv[[speciesCol]])
    names(gcSpecies2) <- gcSpecies2
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
                               sppEquivCol, mortConstraints, growthConstraints, mANPPconstraints) {
  
  #Join these two tables to make the full speciesTraits table
  factorialBiomass[, species := speciesCode]
  setkey(factorialTraits, species)
  setkey(factorialBiomass, species)
  
  #for each species (ie Gamm), find best fit
  species <- names(gamms)
  names(species) <- species
  outputTraits <- lapply(species, FUN = editSpeciesTraits, gamm = gamms,
                         traits = speciesTable, fT = factorialTraits, fB = factorialBiomass,
                         speciesEquiv = sppEquiv, sppCol = sppEquivCol, mortConstraints = mortConstraints,
                         growthConstraints = growthConstraints, mANPPconstraints = mANPPconstraints)
  outputTraitsT <- purrr::transpose(outputTraits)
  gammsT <- purrr::transpose(gamms)
  gammsList <- rbindlist(gammsT$originalData, idcol = "Pair")
  llAll <- rbindlist(outputTraitsT$fullData, idcol = "Pair")
  ll2 <- llAll[, list(BscaledNonLinear = mean(BscaledNonLinear),
                      predNonLinear = mean(predNonLinear)#,
                      #Pair = Pair[1]
                      #mortalityshape = mean(mortalityshape),
                      #mANPPproportion = mean(mANPPproportion)
  ),
  c("species", "standAge", "Pair")]
  gg <- ggplot(ll2, aes(standAge, BscaledNonLinear, colour = species)) + 
    geom_line(size = 2) + 
    geom_point(data = gammsList, aes(standAge, biomass, colour = speciesTemp), size = 0.5) +
    geom_line(size = 2, aes(standAge, predNonLinear, col = species), lty = "dashed") + 
    facet_wrap(~ Pair, nrow = ceiling(sqrt(length(outputTraits))), scales = "fixed") +
    xlim(c(0, 150)) + # ggplot2::scale_y_log() +
    ylim(c(0, 20000)) +
    ggtitle("Comparing best LandR curves (solid) with best Non-Linear fit (dashed)") +
    theme_bw()
  
  #Collapse new traits and replace old traits
  newTraits <- rbindlist(outputTraitsT$bestTraits, fill = TRUE, idcol = "Pair")
  newTraits <- na.omit(newTraits)
  newTraits[, AICWeights := exp( -0.5 * llNonLinDelta)]
  newTraits[, AICWeightsStd := AICWeights/sum(AICWeights), by = "species"]
  
  # Showing species-level averages -- this is not assigned to object
  best <- newTraits[, .(growthcurve = round(sum(AICWeightsStd * growthcurve), 2),
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
  #   message("no GAMMs converged :( Consider revising your parameters")
  # }
  return(list(best = best, gg = gg))
}

modifySpeciesEcoregionTable <- function(speciesEcoregion, speciesTable) {
  message("modifying speciesEcoregion table based on traits derived from PSP Gamms")
  #modify things by species
  newSpeciesEcoregion <- speciesEcoregion[speciesTable, on = c("speciesCode" = "species")]
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
  
  if (identical(species, "all")) {
    matchingSpecies <- unique(speciesEquiv[, .(PSP),])  
    K <- mean(unlist(K))
    q <- mean(unlist(q))
  } else if (grepl(".+_.+", species)) { # this is pairwise
    species2 <- strsplit(species, split = "__")[[1]]
    matchingSpecies <- unique(speciesEquiv[speciesEquiv[[sppCol]] %in% species2, .(PSP),])  
  } else {
    matchingSpecies <- unique(speciesEquiv[speciesEquiv[[sppCol]] == species, .(PSP),])  
  }
  
  if (class(K) == "list"){
    K <- K[[species]]
  }
  
  if (class(q) == "list") {
    q <- q[[species]]
  }
  
  #subset the parameters that may be lists
  if (class(NoOfIters) == "list"){
    NoOfIters <- NoOfIters[[species]]
  }
  
  #Need to calculate stand dominance first - remove all stands < 50% dominance, and of wrong species
  if (NROW(matchingSpecies) != 2) {
    spDominant <- matchingSpecies[psp, nomatch = 0, on = "PSP==newSpeciesName"]  ## filter species first
    spDominant[, length(unique(PSP)), by = "MeasureID"]
    spDominant <- spDominant[spDom > 0.5]
    standData <- spDominant[, .(PlotSize = mean(PlotSize), standAge = mean(standAge),
                                biomass = mean(plotBiomass), spDom = mean(spDom),
                                species = unique(PSP)),
                            by = .(MeasureYear, OrigPlotID1)]
  } else {
    standData <- psp[, .(PlotSize = PlotSize[1], standAge = standAge[1],
                         biomass = spPlotBiomass[1], spDom = spDom[1]),
                     by = .(speciesTemp, MeasureYear, OrigPlotID1)]
  }
  
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
  
  species <- unique(simData$speciesTemp)
  
  localEnv <- environment()
  
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
  species <- as.character(species)
  names(species) <- species
  nlsouts <- lapply(species, function(sp) {
    
    nlsoutInner <- list()
    for (model in names(models)) {
      eqnChar <- models[[model]]$eqn
      plim <- models[[model]]$plim
      Alim <- models[[model]]$Alim
      klim <- models[[model]]$klim
      
      for (ii in 1:700) {
        # Chapman Richards
        # https://www.srs.fs.usda.gov/pubs/gtr/gtr_srs092/gtr_srs092-068-coble.pdf
        nlsoutInner[[model]] <- try(
          robustbase::nlrob(as.formula(paste(eqnChar)),
                            data = simData2[speciesTemp == sp], #maxit = 200,
                            weights = Weights,
                            maxit = 200,
                            start = list(A = runif(1, Alim[["min"]], Alim[["max"]]),
                                         k = runif(1, klim[["min"]], klim[["max"]]),
                                         p = runif(1, plim[["min"]], plim[["max"]])
                            ),
                            trace = FALSE)
          #  , silent = TRUE
        )
        if (!is(nlsoutInner[[model]], "try-error")) {
          break
        }
      }  
    }
    nlsoutInner
  })
  
  
  if (FALSE) {
    for (sp in as.character(species)) {
      plot(simData2$standAge, simData2$biomass, pch = 19, cex = 0.5)
      points(simData2$standAge, fitted(nlsouts[[sp]]$CR), pch = 19, col = "green", cex = 0.5)
      points(simData2$standAge, fitted(nlsouts[[sp]]$Gompertz), pch = 19, col = "blue", cex = 0.5)
      points(simData2$standAge, fitted(nlsouts[[sp]]$Logistic), pch = 19, col = "red", cex = 0.5)
      legend(col = c("green", "blue", "red"),
             legend = paste(names(nlsouts[[sp]]), "AIC: ", sapply(nlsouts[[sp]], function(x) round(AIC(x), 1 ))),
             "topleft", pch = 19)
    }
  }
  nlsout <- lapply(nlsouts, function(nlsIn) {
    nlsIn[[which.min(sapply(nlsIn, function(x) if (is(x, "try-error")) 1e12 else AIC(x)))]]
  })
  # nlsout <- nlsouts[[sp]][[which.min(sapply(nlsouts[[sp]], function(x) if (is(x, "try-error")) 1e12 else AIC(x)))]]
  
  #  if (is(nlsout, "try-error")) {
  speciesGamm <- lapply(species, function(sp) 
    suppressWarnings(try(expr = mgcv::gamm(data = simData2[speciesTemp == sp], formula = eval(gammFormula, enclos = localEnv),
                                           random = list(MeasureYear = eval(randomFormula, envir = baseenv()),
                                                         OrigPlotID1 = eval(randomFormula, envir = baseenv())),
                                           weights = nlme::varFunc(~Weights), verbosePQL = FALSE,
                                           niterPQL = NoOfIters),
                         silent = TRUE))
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

editSpeciesTraits <- function(name, gamm, traits, fT, fB, speciesEquiv, sppCol,
                              mortConstraints, growthConstraints, mANPPconstraints) {
  
  Gamm <- gamm[[name]]
  nameOrig <- name
  if (grepl("__", name)) {
    name <- strsplit(name, "__")[[1]]
    names(name) <- name
  }
  message("Estimating fit for ", nameOrig)
  
  
  #Subset traits to PSP species, return unchanged if no Gamm present
  traits <- traits[species %in% name]
  
  if (class(Gamm) == "try-error" | class(Gamm) == "character") {
    message(paste("not estimating growth/mortality traits for", name))
    return(traits)
  }
  
  predData <- data.table(standAge = min(Gamm$originalData$standAge):max(Gamm$originalData$standAge))
  
  # output <- lapply(name, function(sp) {
  #   out <- try(predict(Gamm$speciesGamm[[sp]], predData, se.fit = TRUE))
  #   if (is(out, "try-error")) browser()
  #   
  # })
  # 
  # nlsoutput <- lapply(name, function(sp) {
  #   tryCatch(predict(Gamm$NonLinearModel[[sp]], predData, se.fit = TRUE), error = function(x) 0)
  # })
  # 
  # if (!is(output, "list")) {
  #   o <- output
  #   output <- list()
  #   output$fit <- o
  #   rm(o)
  # }
  set(fB, NULL, "Sp", gsub(".+_(Sp.)", "\\1", fB$species))
  set(setDT(fT), NULL, "Sp", gsub(".+_(Sp.)", "\\1", fT$species))
  
  # Just 2-species combinations
  candFB <- fB[startsWith(fB$Sp, "Sp")]
  maxBInSims <- 5000
  candFB[, inflationFactor := maxBInSims/max(B), .(pixelGroup)] 
  # candFB[, inflationFactor := maxBInSims/max(B), .(speciesCode)] 
  candFT <- fT[startsWith(fT$Sp, "Sp")]
  maxBiomass <- Gamm$originalData[, .(maxBiomass = max(biomass)), "speciesTemp"]
  setorderv(maxBiomass, "maxBiomass", order = -1L)
  set(maxBiomass, NULL, "Sp", paste0("Sp", 1:2))
  SpNames <- maxBiomass$Sp[match(names(Gamm$NonLinearModel), maxBiomass$speciesTemp)]
  SpMapping <- data.table(Sp = SpNames, species = name)
  names(Gamm$NonLinearModel) <- SpNames
  names(Gamm$speciesGamm) <- SpNames
  setnames(candFB, "age", "standAge")
  constraints <- list(growthcurve = growthConstraints,
                      mortalityshape = mortConstraints,
                      #longevityConstraints,
                      mANPPproportion = mANPPconstraints)
  standAge <- unique(candFB$standAge)
  standAge <- standAge[standAge < 150]
  predGrid <- as.data.table(expand.grid(Sp = SpNames, standAge = standAge))
  
  # Predict from statistical fits to data
  predGrid[, `:=`(
    predNonLinear = predict(Gamm$NonLinearModel[[unlist(.BY)]], .SD),
    predGamm = as.vector(predict(Gamm$speciesGamm[[unlist(.BY)]], .SD, , se.fit = FALSE))
  ), "Sp", .SDcols = "standAge"]
  # misses points past longevity -- have to expand explicitly the standAges for each "speciesCode"
  dt <- as.data.table(expand.grid(speciesCode = unique(candFB$speciesCode), standAge = unique(predGrid$standAge)))
  set(dt, NULL, "Sp", gsub(".+_", "", dt$speciesCode))
  dt <- dt[predGrid, on = c("standAge", "Sp")] 
  candFB <- candFB[dt, on = c("standAge", "speciesCode", "Sp")] # misses points past longevity
  whNA <- which(is.na(candFB$B))
  set(candFB, whNA, "B", 0L)
  candFB[, `:=`(
    pixelGroup = pixelGroup[1],
    inflationFactor = inflationFactor[1]), by = "speciesCode"]
  
  scaleFactorGamm <- max(predGrid$predGamm)/maxBInSims
  scaleFactorNonLinear <- max(predGrid$predNonLinear)/maxBInSims
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
  ll[, c("llNonLinDelta", "llGammDelta") := list(#abs(llDirect - max(llDirect)), 
    abs(llNonLinear - max(llNonLinear)),
    abs(llGamm - max(llGamm)))]
  deltaDiff <- 2 # because we 
  ll <- ll[#llDirectDelta < deltaDiff | 
    llGammDelta < deltaDiff |
      llNonLinDelta < deltaDiff]
  ll <- candFB[ll, on = "pixelGroup"]
  candFTShort <- candFT[unique(ll[, c("species", "pixelGroup")], by = c("species")), on = "species"]
  
  bestGrowthData <- ll[pixelGroup == ll$pixelGroup[1]][, -"species"]
  bestGrowthData <- bestGrowthData[SpMapping, on = "Sp"]
  gg <- ggplot(bestGrowthData, aes(standAge, BscaledNonLinear, colour = species)) + 
    geom_line() + 
    geom_line(aes(standAge, predNonLinear, col = species), lty = "dashed") + 
    geom_point(data = Gamm$originalData, aes(standAge, biomass, colour = speciesTemp)) +
    ggtitle("Comparing best LandR curves (solid) with best Non-Linear fit (dashed)") +
    theme_bw()
  
  # bestPGs <- head(unique(ll$pixelGroup), 6)
  # rr <- ll[pixelGroup %in% bestPGs & standAge == 1][, "speciesCode"]
  rr <- ll[standAge == 1][, c("speciesCode", "llNonLinDelta", "inflationFactor")]
  
  # Take the average of the best
  best <- fT[rr, on = c("species"="speciesCode")]
  bestTraits <- SpMapping[best, on = "Sp"]
  ll <- SpMapping[ll[, -"species"], on = "Sp"]
  
  # bestTraits <- lapply(name, function(sp) {
  #   predData <- data.table("age" = predData$standAge, "predBiomass" = output[[sp]]$fit, 
  #              "predSE" = output[[sp]]$se.fit, "predBiomassNonLin" = nlsoutput[[sp]])
  #   closestLongevity <- abs(fT$longevity - traits$longevity) == min(abs(fT$longevity - traits$longevity)) |
  #     (fT$longevity > max(Gamm$simData$standAge) & (fT$longevity < traits$longevity))
  #   #subset traits by closest longevity
  #   CandidateTraits <- fT[closestLongevity]
  #   
  #   #Constrain growth curve - this is because the effect is conflated with maxANPP
  #   if (class(growthConstraints) == "list") {
  #     growthConstraint <- growthConstraints[[name]]
  #   } else {
  #     growthConstraint <- growthConstraints
  #   }
  #   
  #   #growth constraint may be NULL if user passed constraints on partial list of species
  #   if (!is.null(growthConstraint)){
  #     CandidateTraits <- CandidateTraits[growthcurve %>=% min(growthConstraint)
  #                                        & growthcurve %<=% max(growthConstraint),]
  #   }
  #   #Constrain mortality shape - limited available information on mortalitity in PSPs, too low adds computation strain
  #   if (class(mortConstraints) == "list") {
  #     mortConstraint <- mortConstraints[[name]]
  #   } else {
  #     mortConstraint <- mortConstraints
  #   }
  #   
  #   if (!is.null(mortConstraints)) {
  #     CandidateTraits <- CandidateTraits[mortalityshape >= min(mortConstraint)
  #                                        & mortalityshape <= max(mortConstraint)]
  #   }
  #   
  #   #Constrain growth curve - this is because the effect is conflated with maxANPP
  #   if (class(mANPPconstraints) == "list") {
  #     mANPPconstraint <- mANPPconstraints[[name]]
  #   } else {
  #     mANPPconstraint <- mANPPconstraints
  #   }
  #   
  #   #growth constraint may be NULL if user passed constraints on partial list of species
  #   if (!is.null(mANPPconstraint)){
  #     CandidateTraits <- CandidateTraits[mANPPproportion %>=% min(mANPPconstraint)
  #                                        & mANPPproportion %<=% max(mANPPconstraint),]
  #   }
  #   
  #   
  #   
  #   #subset the simulation values by potential species
  #   CandidateValues <- fB[CandidateTraits]
  #   
  #   #This has to happen before age is subset, or it will underestimate
  #   #  5000 is hard coded because that was the maxB used in the simulations
  #   maxBInSims <- 5000
  #   CandidateValues[, inflationFactor := maxBInSims/max(B), .(speciesCode)] 
  #   #subset the simulation values by age for which curve is represented
  #   CandidateValues <- CandidateValues[age >= min(predData$age) & age <= max(predData$age),]
  #   
  #   setkey(predData, age)
  #   setkey(CandidateValues, age)
  #   CandidateValues <- na.exclude(CandidateValues[predData])
  #   scaleFactor <- max(CandidateValues$predBiomass)/maxBInSims
  #   scaleFactorNonLin <- max(CandidateValues$predBiomassNonLin)/maxBInSims
  #   
  #   #scale factor is the achieved maxB in the simulation / PSP maxB. We use this to scale simulation values to PSP
  #   #inflationFactor is the the LANDIS speciesTrait maxB that was used (always 5000) / simulation's achieved maxB
  #   #scale factor is NOT returned. inflation factor is returned to 'inflate' Biomass_borealDataPrep estimates
  #   # CandidateValues <- CandidateValues[scaleFactors, on = c("speciesCode", "inflationFactor")]
  #   
  #   #Find best possible candidate species
  #   stdev <- sd(predData$predBiomass)
  #   if (FALSE) {
  #     c1 <- CandidateValues[growthcurve %>=% 0.9 & growthcurve %<=% 0.9 &
  #                             (mANPPproportion >= 6 & mANPPproportion <= 6)]
  #     inflFact <- sort(unique(c1$inflationFactor))
  #     c1 <- c1[inflationFactor == inflFact[1]]
  #     par(mfrow = c(2,1))
  #     plot(c1$age, c1$predBiomassNonLin / (max(c1$predBiomassNonLin)/maxBInSims) / c1$inflationFactor, col = "green", pch = 19, cex = .3)
  #     points(c1$age, c1$B, col = "blue", pch = 19, cex = .3)
  #     plot(c1$age, c1$predBiomass / (max(c1$predBiomass)/maxBInSims) / c1$inflationFactor, col = "green", pch = 19, cex = .3)
  #     points(c1$age, c1$B, col = "blue", pch = 19, cex = .3)
  #   }
  #   # c1$predBiomassNonLin / (max(c1$predBiomassNonLin)/maxBInSims) / c1$inflationFactor
  #   #CandidateValues[, `:=`(#LogLikelihood = sum(dnorm(x = B * scaleFactor, mean = predBiomass,
  #   #                                  #                          sd = stdev, log = TRUE)),
  #   #                                  inflationFactor = mean(inflationFactor)), .(speciesCode)]
  #   browser()
  #   Candidates <- CandidateValues[, .(llGamm = sum(dnorm(x = B * scaleFactor * inflationFactor, mean = predBiomass,
  #                                                        sd = stdev, log = TRUE)),
  #                                     llNonLin = sum(dnorm(x = B * scaleFactor * inflationFactor, mean = predBiomassNonLin,
  #                                                          sd = stdev, log = TRUE))), .(speciesCode)]
  #   CandidatesAll <- CandidateTraits[Candidates, on = c("species" = "speciesCode")] 
  #   
  #   CandidateValuesDirectData <- Gamm$simData[CandidateValues, on = c("standAge"="age"), allow.cartesian = TRUE, nomatch = 0]
  #   CandidateDirectData <- CandidateValuesDirectData[, .(llDirect = sum(dnorm(x = biomass, mean = B * scaleFactor * inflationFactor,
  #                                                                             sd = stdev, log = TRUE))
  #   ), .(speciesCode)]
  #   CandidatesAll <- CandidatesAll[CandidateDirectData, on = c("species" = "speciesCode")]
  #   setorderv(CandidatesAll, "llNonLin", order = -1L)
  #   CandidatesAll[, llNonLinRank := rev(frank(llNonLin))]
  #   setorderv(CandidatesAll, "llGamm", order = -1L)
  #   CandidatesAll[, llGammRank := rev(frank(llGamm))]
  #   setorderv(CandidatesAll, "llDirect", order = -1L)
  #   CandidatesAll[, llDirectRank := rev(frank(llDirect))]
  #   
  #   bestCandidates <- CandidatesAll[, c("llDirectDelta", "llNonLinDelta", "llGammDelta") := list(abs(llDirect - max(llDirect)), 
  #                                                                                                abs(llNonLin - max(llNonLin)),
  #                                                                                                abs(llGamm - max(llGamm)))]
  #   deltaDiff <- 5
  #   bestCandidates <- bestCandidates[llDirectDelta < deltaDiff | 
  #                                      llGammDelta < deltaDiff |
  #                                      llNonLinDelta < deltaDiff]
  #   # bestCandidates <- bestCandidates[ llNonLinDelta %==% 0 | llGammDelta %==% 0 | llDirectDelta %==% 0]
  #   ccv <- CandidateValues[bestCandidates[, c("species", "llGammDelta", "llNonLinDelta", "llDirectDelta")], 
  #                          on = c("speciesCode" = "species")]
  #   ccvNonLin <- ccv[llNonLinDelta %==% 0]
  #   ccvGamm <- ccv[llGammDelta %==% 0]
  #   ccvDirect <- ccv[llDirectDelta %==% 0]
  #   
  #   # bestCandidate <- bestCandidates[ llNonLinDelta %==% 0]
  #   bestCandidate <- bestCandidates
  #   CandidateValues <- CandidateValues[speciesCode %in% bestCandidate$species]
  #   CandidateValues <- unique(CandidateValues[, -c("age", "B", "predBiomass", "predBiomassNonLin", "predSE")], by = c("growthcurve", "mortalityshape", "inflationFactor"))
  #   bestCandidate <- CandidateValues[bestCandidate[, c("species", "llGammDelta", "llNonLinDelta", "llDirectDelta")], 
  #                                    on = c("speciesCode" = "species")]
  #   
  #   if (FALSE) { # This is about plotting (Eliot Nov 18, 2021) -- couldn't fit into the module in a modular way 
  #     ma <- match(name, names(gamm))
  #     col <- factor(names(gamm))[ma]
  #     if (FALSE) { # This is about plotting (Eliot Nov 18, 2021) -- couldn't fit into the module in a modular way 
  #       if (ma == 1)
  #         plot(Gamm$simData$standAge, Gamm$simData$biomass, pch = 19, cex = 0.5, col = ma, xlim = c(0, 200), ylim = c(0, 40000))
  #       else 
  #         points(Gamm$simData$standAge, Gamm$simData$biomass, pch = 19, cex = 0.5, col = ma)
  #       
  #       # THese are the LandR equivalents
  #       if (!isTRUE(nlsoutput == 0))
  #         # direct from data
  #         lines(ccvDirect$age, ccvDirect$B * scaleFactor * ccvDirect$inflationFactor, lty = "solid", lwd = 4, col = col )
  #       
  #       lines(ccvGamm$age, ccvGamm$B * scaleFactor * ccvGamm$inflationFactor, lty = "dotted", lwd = 4, col = col )
  #       lines(ccvNonLin$age, ccvNonLin$B * scaleFactor * ccvNonLin$inflationFactor, lty = "dashed", lwd = 2, col = col )
  #       
  #       # These are the Gamm or NonLin
  #       lines(ccvGamm$age, ccvGamm$predBiomass, col = col, lty = "dotted", lwd = 1)
  #       lines(ccvNonLin$age, ccvNonLin$predBiomassNonLin, col = col, lty = "dashed", lwd = 1)
  #       
  #       if (ma == length(gamm)) {
  #         legend("topleft",
  #                lwd = c(4, 4, 2, 1, 1),
  #                lty = c("solid", "dotted", "dashed", "dotted", "dashed"),
  #                legend = c("LandR curve directly from points",
  #                           "LandR curve from gamm (current method)", 
  #                           "LandR curve from Chapman Richards",
  #                           "gamm",
  #                           "Chapman Richards"))
  #         # legend("bottomright",
  #         #        pch = 19,
  #         #        col = c(4, 5, 7),
  #         #        legend = c("W Spruce", "B Spruce", "T Aspen"))
  #       }
  #     }
  #     
  #   }
  #   
  #   
  #   
  #   # What to do with the LogLikelihood? Report?
  #   bestTraits <- bestCandidate[, .(mortalityshape, growthcurve, mANPPproportion, inflationFactor, llNonLinDelta, llGammDelta, llDirectDelta)]
  #   #if there are multiple rows due to undifferentiated curves, then mortality hasn't kicked in. Take the max mortalityshape
  #   # bestTraits <- bestTraits[mortalityshape == max(mortalityshape)]
  #   bestTraits[, mortalityshape := asInteger(mortalityshape)]
  # })
  # if (NROW(bestTraits) > 1) {
  #   traits <- rbindlist(lapply(seq(NROW(bestTraits$Abie_bal)), function(x) traits))
  # } 
  # 
  # traits[, colnames(bestTraits) :=  bestTraits] 
  return(list(bestTraits = bestTraits, fullData = ll))
}

makePSPgamms <- function(studyAreaANPP, PSPperiod, PSPgis, PSPmeasure,
                         PSPplot, useHeight, biomassModel, speciesCol,
                         sppEquiv, NoOfIterations, knots, minimumSampleSize,
                         quantileAgeSubset, minDBH, allSpInOne) {
  
  #this function is just a wrapper around these functions, for caching purposess
  psp <- prepPSPaNPP(studyAreaANPP = studyAreaANPP, PSPperiod = PSPperiod,
                     PSPgis = PSPgis, PSPmeasure = PSPmeasure, PSPplot = PSPplot,
                     useHeight = useHeight, biomassModel = biomassModel, minDBH = minDBH)
  
  #Wrapper used to avoid caching psp object - too large
  speciesGAMMs <- buildGrowthCurves(PSPdata = psp, speciesCol = speciesCol,
                                    sppEquiv = sppEquiv, NoOfIterations = NoOfIterations,
                                    knots = knots, minimumSampleSize = minimumSampleSize,
                                    quantileAgeSubset = quantileAgeSubset,
                                    allSpInOne = allSpInOne)
  return(speciesGAMMs)
}

gammFormula <- quote(biomass ~ s(standAge, k = K, pc = 0))
randomFormula <- quote(~1)


updateInflationFactor <- function(traits, factorialTraits, factorialBiomass) {
  longs <- traits$longevity
  names(longs) <- traits$species
  closestLongevity <- sapply(longs, function(long) {
    unique(factorialTraits$longevity[abs(factorialTraits$longevity - long) == 
                                       min(abs(factorialTraits$longevity - long))])
  })
  newTraitsClosestLongevity <- copy(traits)
  newTraitsClosestLongevity[, longevity := closestLongevity]
  
  maxBInSims <- 5000
  fT <- copy(factorialTraits)
  fB <- copy(factorialBiomass)
  setnames(fT, old = "species", new = "speciesDummy")
  setnames(fB, old = "speciesCode", new = "speciesDummy")
  set(fB, NULL, "species", NULL)
  factorialSp <- fT[newTraitsClosestLongevity, on = c("growthcurve", "mortalityshape", "mANPPproportion", "longevity")]
  factorialSp <- fB[speciesDummy %in% factorialSp$speciesDummy, 
                    list(inflationFactor = maxBInSims/max(B)), by = "speciesDummy"][
                      factorialSp[, -"inflationFactor"], on = "speciesDummy"
                    ]
  
  # Update traits with new inflationFactor
  traits <- factorialSp[, c("species", "inflationFactor")][traits[, -"inflationFactor"], on = "species"]
  traits[]
}