prepPSPaNPP <- function(studyAreaANPP, PSPgis, PSPmeasure, PSPplot, sppEquivLong,
                        useHeight, biomassModel, PSPperiod, minDBH, sppEquiv = NULL,
                        mergeHybridSprucePSP = NA_character_, excludeBECzones = character(0),
                        BECzones = NULL) {
  ## crop points to studyArea
  if (!is.null(studyAreaANPP)) {
    studyAreaANPP <- st_as_sf(studyAreaANPP) # in case SPDF
    studyAreaANPP <- st_transform(x = studyAreaANPP, crs = st_crs(PSPgis))
    message(cli::col_yellow("Filtering PSPs for ANPP to study Area..."))
    PSP_sa <- PSPgis[studyAreaANPP,] |>
      as.data.table() |> # PSPgis may be a plain sf data.frame, which setkey() rejects
      setkey(OrigPlotID1)
    message(cli::col_yellow(paste0("There are "), nrow(PSP_sa), " PSPs in your study area"))
    if (nrow(PSP_sa) == 0) {
      stop("subsetting PSPs to those within studyAreaANPP appears to ",
      "have removed all of them...please review")
    }
    ## filter other PSP datasets to those in study Area
    PSPmeasure <- PSPmeasure[OrigPlotID1 %in% PSP_sa$OrigPlotID1, ]
    PSPplot <- PSPplot[OrigPlotID1 %in% PSP_sa$OrigPlotID1, ]
  }

  ## filter data by study period
  message(cli::col_yellow("Filtering PSPs for ANPP by study period..."))
  PSPmeasure <- PSPmeasure[MeasureYear > min(PSPperiod) & MeasureYear < max(PSPperiod), ]
  PSPplot <- PSPplot[MeasureYear > min(PSPperiod) & MeasureYear < max(PSPperiod), ]
  PSPmeasure <- PSPmeasure[PSPplot, on = c("MeasureID", "OrigPlotID1", "MeasureYear", "source")]

  ## TODO: this should be parameterized - besides its tree density you want not the raw count
  ## filter by > 30 trees at first measurement (P) to ensure forest.
  forestPlots <- PSPmeasure[MeasureYear == baseYear, .(measures = .N), OrigPlotID1] %>%
    .[measures >= 30, ]

  PSPmeasure <- PSPmeasure[OrigPlotID1 %in% forestPlots$OrigPlotID1, ]
  PSPplot <- PSPplot[OrigPlotID1 %in% PSPmeasure$OrigPlotID1, ]

  ## Restrict to trees >= ## DBH. Maybe necessary to get rid of small trees when they are inconsistently recorded
  PSPmeasure <- PSPmeasure[DBH >= minDBH,]
  ## decide what to do about above line and stem/density. PES approach is to just fix the data...
  
  #get column for estimating biomass
  PSPmeasure <- sppEquivLong[PSPmeasure, on = c("Latin_full" = "Species")]
  PSPmeasure[is.na(SpBiomassEq), SpBiomassEq := ""] #biomassCalculation errors with NA

  #Calculate biomass
  #Height must be calculated separately if there are NA heights -
  if (useHeight & !is.null(PSPmeasure$Height)) {
    PSPmeasureNoHeight <- PSPmeasure[is.na(Height)]
    PSPmeasureHeight <- PSPmeasure[!is.na(Height)]
    tempOut <- biomassCalculation(species = PSPmeasureHeight$SpBiomassEq,
                                  DBH = PSPmeasureHeight$DBH,
                                  height = PSPmeasureHeight$Height,
                                  includeHeight = TRUE,
                                  equationSource = biomassModel)
    ## check if height is missing, join if so -- function fails if data.table is empty
    if (nrow(PSPmeasureNoHeight) > 0) {
      tempOutNoHeight <- biomassCalculation(species = PSPmeasureNoHeight$SpBiomassEq,
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
  } else {
    tempOut <- biomassCalculation(species = PSPmeasure$SpBiomassEq,
                                  DBH = PSPmeasure$DBH,
                                  height = PSPmeasure$Height,
                                  includeHeight = useHeight,
                                  equationSource = biomassModel)
    PSPmeasure$biomass <- tempOut$biomass
  }
  # message(cli::col_yellow("No PSP biomass estimate possible for these species: "))
  # message(cli::col_yellow(paste(unique(tempOut$missedSpecies), collapse = ", ")))

  ## BC hybrid-zone "Picea glauca" -> the hybrid that sppEquiv merged into Pice_eng. After the
  ## biomass calculation: biomass keeps the species' own equation, only the species grouping follows
  if (!is.null(sppEquiv)) {
    PSPmeasure <- relabelHybridSprucePSP(PSPmeasure, PSPgis, sppEquiv, mergeHybridSprucePSP,
                                         excludeBECzones, BECzones, speciesCol = "Latin_full")
  }

  #bad biomass estimate
  PSPmeasure <- PSPmeasure[biomass != 0]
  #keep objects as small as possible
  PSPmeasure[, c("status", "first_tree_year", "last_tree_year", "diff_dbh") := NULL]

  ## add stand age estimate
  PSPmeasure[, standAge := baseSA + MeasureYear - baseYear]
  PSPmeasure <- PSPmeasure[standAge > 0]
  PSPmeasure <- PSPmeasure[!is.na(biomass)]

  ## divide biomass by 10 to convert kg/ha to g/m2, the LandR unit
  PSPmeasure[, biomass := biomass/10]

  return(PSPmeasure)
}

buildGrowthCurves <- function(PSPdata, speciesCol, sppEquiv, quantileAgeSubset,
                              minimumSampleSize, speciesFittingApproach = "single") {
  if (identical(speciesFittingApproach, "all")) {
    stop("speciesFittingApproach = 'all' is no longer supported: it has failed since December 2021, ",
         "and the current fitting has no pooled all-species fit. ",
         "Use 'focal' (default), 'pairwise' or 'single'.")
  }
  if (!isTRUE(speciesFittingApproach %in% c("focal", "pairwise", "single"))) {
    stop("speciesFittingApproach must be one of 'focal', 'pairwise' or 'single'.")
  }
  ## must filter PSPdata by all sppEquiv$PSP with same sppEquivCol
  if (!isTRUE(speciesCol %in% colnames(sppEquiv))) {
    stop("sppEquivCol not in sppEquiv")
  }

  gcSpecies <- speciesFittingApproach
  if (identical(speciesFittingApproach, "single")) {
    gcSpecies <- unique(sppEquiv[[speciesCol]])
  }
  gcSpecies <- setNames(nm = gcSpecies)

  SpPSP <- copy(PSPdata)
  #SpPSP has SpBiomassEq, Latin_full and PSP
  sppEquivShort <- unique(sppEquiv[, .SD, .SDcols = c("Latin_full", speciesCol)])
  
  SpPSP[, speciesTemp := LandR::equivalentName(value = Latin_full, df = sppEquiv,
                                               column = speciesCol)]
  #probably do setnames, Species, and make things NULL
  
  
  
  
  SpPSP[, "areaAdjustedB" := biomass/PlotSize]
  SpPSP  <- SpPSP[, "plotBiomass" := sum(areaAdjustedB), .(MeasureID)]
  SpPSP[, "spPlotBiomass" := sum(areaAdjustedB), .(MeasureID, speciesTemp)]
  SpPSP[, "spDom" := spPlotBiomass/plotBiomass, .(MeasureID)]

  #remove these AFTER calculating stand biomass and dominance
  whNA <- is.na(SpPSP$speciesTemp)
  if (any(whNA)) {
    message(cli::col_yellow("Removing ", paste(unique(SpPSP$Latin_full[whNA]), collapse = ", ")))
    message(cli::col_yellow("   ... because they are not in sppEquiv"))
    SpPSP <- SpPSP[!whNA]
  }
  SpPSP <- SpPSP[Latin_full %in% sppEquiv[["Latin_full"]]]
  freq <- SpPSP[, .(N = .N, spDom = spDom[1]), .(speciesTemp, MeasureID)]
  SpPSP[, c("SpBiomassEq", "source", "PSP", "Latin_full", "Elevation") := NULL]
  SpPSP[, speciesTemp := as.factor(speciesTemp)]
  if (speciesFittingApproach %in% c("focal", "pairwise")) {
    ## subset to species-of-interest using relative biomass (dominance)
    freq <- freq[spDom > 0.2] # minimum 20% dominance for pairwise or focal

    speciesComp <- freq[, .(numSp = .N, spComp = paste(speciesTemp, collapse = "__")), by = "MeasureID"]
    ## Pick only 2 species plots when using "pairwise"
    if (isTRUE(speciesFittingApproach == "pairwise")) {
      speciesComp <- speciesComp[numSp == 2]
    }
    speciesCompN <- speciesComp[, .N, by = "spComp"]

    gcSpecies1 <- unique(sppEquiv[[speciesCol]]) |> setNames(nm = _)
    gcSpecies2 <- unique(speciesComp$spComp) |> setNames(nm = _)
    speciesForSplit <- if (isTRUE(speciesFittingApproach == "pairwise")) gcSpecies2 else gcSpecies1
    speciesCompListAll <- split(speciesComp, speciesComp$spComp)
    speciesCompList <- lapply(speciesForSplit, function(spName) {
      rbindlist(speciesCompListAll[grepl(spName, names(speciesCompListAll))])
    })
    SpPSPList <- lapply(speciesCompList, function(speciesCompInner) {
      SpPSP[MeasureID %in% speciesCompInner$MeasureID][spDom > 0.2]
    })
    if (isTRUE(speciesFittingApproach == "focal")) {
      SpPSPList <- Map(psp = SpPSPList, spName = names(SpPSPList), function(psp, spName) {
        psp[speciesTemp != spName, speciesTemp := "Other"]
        ## pool the non-focal species: "Other" biomass and dominance are their sum in each plot measurement
        psp[speciesTemp == "Other", spPlotBiomass := sum(areaAdjustedB), by = "MeasureID"]
        psp[speciesTemp == "Other", spDom := spPlotBiomass / plotBiomass]
      })
    }
  } else {
    ## subset to species-of-interest using relative biomass (dominance)
    SpPSP <- SpPSP[spDom > 0.5] # minimum 50% dominance for single
    SpPSPList <- split(SpPSP, SpPSP$speciesTemp)
  }
  
  speciesForCurves <- names(SpPSPList) |> setNames(nm = _)
  message(cli::col_yellow("-----------------------------------------------"))
  message(cli::col_yellow("building growth curves from PSP data: "))
  outputGCs <- Map(
    f = buildModels,
    species = speciesForCurves,
    psp = SpPSPList,
    MoreArgs = list(
      speciesEquiv = sppEquiv,
      sppCol = speciesCol,
      minSize = minimumSampleSize,
      q = quantileAgeSubset
    )
  )

  return(outputGCs)
}

buildModels <- function(species, psp, speciesEquiv, sppCol, minSize, q) {
  if (class(q) == "list") {
    q <- q[[species]]
  }

  standData <- psp[, .(PlotSize = PlotSize[1], standAge = standAge[1],
                       biomass = spPlotBiomass[1], spDom = spDom[1]),
                   by = .(speciesTemp, MeasureYear, OrigPlotID1)]

  ## test if there are sufficient plots to estimate traits
  if (nrow(standData) < minSize) {
    GC <- "insufficient data"
    names(GC) <- species

    warning("insufficient PSP data to estimate traits:\n",
            "(have ", nrow(standData), " rows; requested minimum ", minSize, " rows)")

    return(GC)
  }

  ## By default removing the 95th percentile of age - these points are usually too scattered to produce reliable estimates
  standData <- standData[standAge < quantile(standData$standAge, probs = q/100),]
  simulatedData <- simulateYoungStands(cohortData = standData, N = 50)
  simData <- rbindlist(list(standData, simulatedData), fill = TRUE)

  ## This weights the real data by spDominance, without distorting the mean of fake data.
  Realweights <- standData$spDom/mean(standData$spDom)
  Fakeweights <- rep(1, times = nrow(simulatedData))
  simData$Weights <- c(Realweights, Fakeweights)

  species <- na.omit(unique(simData$speciesTemp))

  ## Spread out the 0s so they are in the first 10 years.
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
  species <- as.character(species) |> setNames(nm = _)
  speciesForFits <- setdiff(species, "Other") |> setNames(nm = _)
  speciesForFitsMessage <- paste(speciesForFits, collapse = ", ")
  message(cli::col_yellow(
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

  ## Append the true data so we don't have the 0s involved when we subset by age quantile
  sppGrowthCurves <- list(originalData = standData,
                          simData = simData,
                          NonLinearModel = nlsout)

  return(sppGrowthCurves)
}
