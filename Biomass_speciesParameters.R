defineModule(sim, list(
  name = "Biomass_speciesParameters",
  description = paste(
    "Estimates the `growthcurve`, `mANPPproportion` and `inflationFactor` species traits of `Biomass_core`",
    "from Permanent Sample Plot (PSP) and Temporary Sample Plot (TSP) data. A non-linear biomass-by-age model is",
    "fitted to each species' PSP data and matched to a single-species `Biomass_speciesFactorial` (which gives",
    "`mANPPproportion` and `inflationFactor` for every `growthcurve`); `growthcurve` is then calibrated with",
    "`Biomass_core` simulations of remeasured, mixed-species PSPs."
  ),
  keywords = NA, # TODO: add keywords
  authors = c(
    person(c("Ian"), "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person(c("Eliot"), "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut")),
    person(c("Ceres"), "Barros", email = "ceres.barros@ubc.ca", role = c("ctb"))
  ),
  childModules = character(0),
  version = list(Biomass_speciesParameters = "4.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "Biomass_speciesParameters.Rmd"),
  loadOrder = list(after = c("Biomass_speciesFactorial", "Biomass_borealDataPrep"),
                   before = c("Biomass_core")),
  reqdPkgs = list(
    "arrow", "bcdata", "cli", "data.table", "dplyr", "fpCompare", "fs", "ggplot2", "gridExtra",
    "googledrive", "httr2", "mgcv", "nlme", "purrr", "qs2", "robustbase", "sf", "terra", "viridis",
    "reproducible (>= 2.1.0)",
    "SpaDES.core (>= 2.1.4)",
    "PredictiveEcology/LandR@development (>= 1.1.0.9077)",
    "PredictiveEcology/pemisc@development (>= 0.0.3.9002)",
    "ianmseddy/PSPclean@development(>= 1.0.0.9001)"
  ),
  parameters = rbind(
    defineParameter("biomassModel", "character", "Lambert2005", NA, NA,
                    desc =  paste("The model used to calculate biomass from DBH. Can be either 'Lambert2005' or 'Ung2008'.")),
    defineParameter("calibrateGrowthcurve", "logical", TRUE, NA, NA,
                    desc = paste("If `TRUE`, `growthcurve` is calibrated with `Biomass_core` simulations of remeasured PSPs",
                                 "(see the module manual). This needs `ecoregionMap` and `speciesEcoregion` from",
                                 "`Biomass_borealDataPrep`. If `FALSE`, or if there are too few remeasured plots, the species'",
                                 "`growthcurve` in `species` is kept and only `mANPPproportion` and `inflationFactor` are",
                                 "estimated, for that `growthcurve`.")),
    defineParameter("calibrationIterations", "integer", 2L, 1L, NA,
                    desc = paste("Number of passes of the `growthcurve` calibration. Competitors are simulated",
                                 "with the `growthcurve` found by the previous pass.")),
    defineParameter("calibrationMaxInterval", "numeric", 30, NA, NA,
                    desc = "Longest interval (years) between PSP measurements that is simulated in the calibration."),
    defineParameter("calibrationMinFocalProportion", "numeric", 0.2, 0, 1,
                    desc = paste("A species is the focal species of a remeasured PSP, i.e., its `growthcurve` is scored",
                                 "on it, if it is at least this proportion of the plot's biomass at the first measurement.")),
    defineParameter("calibrationMinInterval", "numeric", 3, NA, NA,
                    desc = "Shortest interval (years) between PSP measurements that is simulated in the calibration."),
    defineParameter("landis", "logical", FALSE, NA, NA,
                    desc = paste("If `TRUE`, run in 'LANDIS mode': expose the fitted LANDIS-version growth curves",
                                 "(the scaled non-linear biomass-over-age curves shown in the LandR-vs-non-linear plot)",
                                 "per species as the output object `speciesGrowthCurvesLandis`, for use as inputs to",
                                 "LANDIS-II Biomass Succession. The per-species growth-curve parameters",
                                 "(`growthcurve`, `mortalityshape`, etc.) are written to `species` regardless.",
                                 "Default `FALSE` preserves the standard behaviour.")),
    defineParameter("maxBInFactorial", "integer", 5000L, NA, NA,
                    desc = paste("The arbitrary maximum biomass for the factorial simulations.",
                                 "This is a per-species maximum within a pixel")),
    defineParameter("maxCalibrationPlots", "integer", 500L, 1L, NA,
                    desc = paste("Maximum number of pairs of consecutive PSP measurements simulated in each pass of the",
                                 "`growthcurve` calibration. Each is simulated once per focal species and candidate",
                                 "`growthcurve`. A random subset is used if there are more.")),
    defineParameter("excludeBECzonesHybridSpruce", "character", c("BWBS", "SWB"), NA, NA,
                    desc = paste("BC BEC zones whose PSPs are NOT relabelled by `mergeHybridSprucePSP`:",
                                 "the boreal zones, where 'Picea glauca' is true white spruce.")),
    defineParameter("mergeHybridSprucePSP", "character", getOption("LandR.mergeHybridSpruce", "engelmann"), NA, NA,
                    desc = paste("Mirrors `LandR.mergeHybridSpruce`, and follows it by default. With 'engelmann',",
                                 "and the hybrid spruce merged into Pice_eng in `sppEquiv`, BC PSP records of",
                                 "'Picea glauca' (how BC's interior hybrid-zone spruce is recorded) are relabelled",
                                 "'Picea engelmannii x glauca' so they count as Pice_eng, except in the BEC zones",
                                 "of `excludeBECzonesHybridSpruce`. 'white' and NA do not relabel.")),
    defineParameter("minimumPlots", "numeric", 50, 10, NA,
                    desc = paste("Minimum number of PSP plots per species")),
    defineParameter("minDBH", "integer", 0L, 0L, NA,
                    desc = paste("Minimum diameter at breast height (DBH) in cm used to filter PSP data.",
                                 "Defaults to 0 cm, i.e. all tree measurements are used.")),
    defineParameter("PSPdataTypes", "character", "all", NA, NA,
                    desc = paste("Which PSP datasets to source, defaulting to all. Other available options include",
                                 "'BC', 'AB', 'SK', 'ON', 'NB', 'NFI', and 'dummy'.",
                                 "'dummy' should be used for unauthorized users.")),
    defineParameter("PSPperiod", "numeric", c(1920, 2019), NA, NA,
                    desc = paste("The years by which to subset sample plot data, if desired. Must be a vector of length 2")),
    defineParameter("quantileAgeSubset", "numeric", 99, 1, 100,
                    desc = paste("Quantile by which to subset PSP data. As older stands are sparsely represented",
                                 "the oldest measurements become vastly more influential. This parameter accepts",
                                 "both a single value and a list of vectors, named according to `sppEquivCol`.")),
    defineParameter("sppEquivCol", "character", "LandR", NA, NA,
                    paste("The column in `sim$sppEquiv` data.table that defines individual species.",
                          "The names should match those in the species table.")),
    defineParameter("standAgesForFitting", "integer", c(21L, 91L), NA, NA,
                    desc = paste("The minimum and maximum ages of the biomass-by-age curves used to compare the fitted",
                                 "models with the factorial. The factorial's ages (every 10 years from 1) within this",
                                 "range are used. It is generally recommended to keep this param under 200, given the low data",
                                 "availability of stands aged 200+, with some exceptions.")),
    defineParameter("useHeight", "logical", TRUE, NA, NA,
                    desc = paste("Should height be used to calculate biomass (in addition to DBH).",
                                 "DBH is used by itself when height is missing.")),
    defineParameter(".plots", "character", "screen", NA, NA,
                    desc = "Used by Plots function, which can be optionally used here"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    desc = "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time interval between save events"),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    desc = paste("Human-readable name for the growth curve filename.",
                                 "If `NA`, a hash of sppEquiv[[sppEquivCol]] will be used.")),
    defineParameter(".useCache", "character", c(".inputObjects", "init"), NA, NA,
                    desc = paste("Should this entire module be run with caching activated?",
                                 "This is generally intended for data-type modules,",
                                 "where stochasticity and time are not relevant")),
    defineParameter(".useParallel", "integer", 2L, NA, NA,
                    desc = paste("maximum number of threads/workers to use for data.table operations;",
                                 "passed to `data.table::setDTthreads` and should be <= 4."))
  ),
  inputObjects = bindrows(
    expectsInput("cohortDataFactorial_path", "fs_path",
                 desc = paste(
                   "Path where the `cohortDataFactorial` object is saved as an `arrow` dataset.",
                   "A `cohortData` table (**sensu** `Biomass_core`) with columns `age`, `B`,",
                   "and `speciesCode` that joins with `speciesTableFactorial`. It must come from",
                   "`Biomass_speciesFactorial` with `factorialDesign = 'perSpecies'`: one species per pixelGroup, no competition.",
                   "See `PredictiveEcology/Biomass_speciesFactorial` for further information."
                  ),
                 sourceURL = NA),
    expectsInput("ecoregionMap", "SpatRaster",
                 desc = paste("Categorical raster of `ecoregionGroup`, from `Biomass_borealDataPrep`. Used to place PSPs in",
                              "ecoregions when calibrating `growthcurve`. Optional: without it `growthcurve` is not calibrated.")),
    expectsInput("PSPmeasure_sppParams", "data.table",
                 desc = paste("Merged PSP and TSP individual tree measurements. Must include the following columns:",
                              "`MeasureID`, `OrigPlotID1`, `MeasureYear`, `TreeNumber`, `Species`, `DBH` and `PSP`,",
                              "where `Species` corresponds to species names in `LandR::sppEquivalencies_CA$Latin_full`.",
                              "Defaults to randomized PSP data stripped of real `plotID`s"),
                 sourceURL = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing"),
    expectsInput("PSPplot_sppParams", "data.table",
                 desc = paste("Merged PSP and TSP plot data. Defaults to randomized PSP data stripped of real `plotID`s.",
                              "Must contain columns `MeasureID`, `MeasureYear`, `OrigPlotID1`, and `baseSA`,",
                              "the latter being stand age at year of first measurement"),
                 sourceURL = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing"),
    expectsInput("BECzonesBC", "sf",
                 desc = paste("Optional. BC BEC zone polygons, with a `ZONE` column, used to exclude the zones of",
                              "`excludeBECzonesHybridSpruce` from the hybrid-spruce relabel of PSPs.",
                              "If not supplied, they are fetched from the BC Data Catalogue",
                              "(WHSE_FOREST_VEGETATION.BEC_BIOGEOCLIMATIC_POLY); if that fails there is no relabel."),
                 sourceURL = NA),
    expectsInput("PSPgis_sppParams", "sf",
                 desc = paste("Plot location `sf` object. Defaults to PSP data stripped of real `plotID`s/location.",
                              "Must include field `OrigPlotID1` for joining to `PSPplot` object"),
                 sourceURL = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing"),
    expectsInput("species", "data.table",
                 desc = paste("A table of invariant species traits with the following trait colums:",
                              "'species', 'Area', 'longevity', 'sexualmature', 'shadetolerance',",
                              "'firetolerance', 'seeddistance_eff', 'seeddistance_max', 'resproutprob',",
                              "'mortalityshape', 'growthcurve', 'resproutage_min', 'resproutage_max',",
                              "'postfireregen', 'wooddecayrate', 'leaflongevity' 'leafLignin', and 'hardsoft'.",
                              "'longevity' and 'mortalityshape' are taken as known; 'growthcurve' is the starting value for",
                              "the calibration and is replaced, with 'mANPPproportion' and 'inflationFactor', for fitted species.",
                              "Default is from Dominic Cyr and Yan Boulanger's applications of LANDIS-II"),
                 sourceURL = "https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv"),
    expectsInput("speciesEcoregion", "data.table",
                 desc = paste("Table of spatially-varying species traits (`maxB`, `maxANPP`,",
                              "`establishprob`), defined by species and `ecoregionGroup`).",
                              "Defaults to a dummy table based on dummy data of biomass, age,",
                              "ecoregion and land cover class.")),
    expectsInput("speciesTableFactorial_path", "fs_path",
                 desc = paste(
                   "Path where the `speciesTableFactorial` object is saved as an `arrow` dataset.",
                   "A species table (**sensu** `Biomass_core`) with the traits of `cohortDataFactorial`: `species`",
                   "(or `speciesCode`), `pixelGroup`, `targetSpecies` (the species the pixelGroup stands in for),",
                   "`longevity`, `growthcurve`, `mortalityshape` and `mANPPproportion`.",
                   "See `PredictiveEcology/Biomass_speciesFactorial` for further information."
                 ),
                 sourceURL = NA),
    expectsInput("sppEquiv", "data.table",
                 desc = paste("Table of species equivalencies - see `?LandR::sppEquivalencies_CA`.",
                              "Traits will be estimated for each unique entry in the `sppEquivCol` column.")),
    expectsInput("sppEquivLong", "data.table",
                 desc = paste("The full table of species equivalencies - see `?LandR::sppEquivalencies_CA`.",
                              "Biomass will be estimated for each species based on the `sp_Biomass_eq' column,",
                              "which uses `pemisc::biomassCalculation` to derive AGB from DBH and height",
                              "(based on the equations from https://doi.org/10.1139/x05-112). The full table",
                              "is used to improve stand biomass estimates even if some species are not of interest.")),
    expectsInput("studyAreaANPP", "sf",                                         # nolint: in_no_default
                 desc = paste("Optional study area used to crop PSP data before building growth curves.",
                              "If supplied, an ecoregion-scale object is recommended, at a minimum."))
  ),
  outputObjects = bindrows(
    createsOutput("species", "data.table",
                  desc = paste("The updated invariant species traits table (see above), with `growthcurve`,",
                               "`mANPPproportion` and `inflationFactor` of the fitted species.")),
    createsOutput("speciesEcoregion", "data.table",
                  desc = paste("The updated spatially-varying species traits table",
                               "(see description for this object in inputs)")),
    createsOutput("speciesGrowthCurves", "list",
                  desc = paste("list containing each species' non-linear model,",
                               "model data, and the unfiltered PSP data")),
    createsOutput("speciesGrowthcurveScore", "data.table",
                  desc = paste("An empty `data.table` unless `growthcurve` was calibrated, when it has the scores of the last",
                               "pass for every species and candidate `growthcurve`: the number of simulated plots `n`, `mse`",
                               "(mean squared error of the focal species' biomass at the second measurement, relative to its",
                               "`maxB`), `ll` (the Gaussian profile log-likelihood) and `dLL` (`ll` minus the best `ll` of the species).")),
    createsOutput("speciesGrowthCurvesLandis", "data.table",
                  desc = paste("An empty `data.table` unless `P(sim)$landis` is `TRUE`, when it holds the fitted",
                               "LANDIS-version growth curves (`BscaledNonLinear`) by `species` and `standAge` (the",
                               "scaled non-linear curves shown in the `LandR_VS_NLM_growthCurves` plot), for use as",
                               "LANDIS-II Biomass Succession inputs.")),
    createsOutput("speciesGrowthCurvesPSP", "data.table",
                  desc = paste("An empty `data.table` unless `P(sim)$landis` is `TRUE`, when it holds the PSP",
                               "observations used to fit the growth curves (`biomass` by `standAge` and species, with",
                               "`OrigPlotID1` for joining to plot locations / ecoregion), a diagnostic to plot against",
                               "`speciesGrowthCurvesLandis`.")),
    createsOutput("speciesTraitMatch", "data.table",
                  desc = paste("For each fitted species and each candidate `growthcurve` of the factorial, the",
                               "`mANPPproportion` and `inflationFactor` for which the factorial's net peak biomass increment",
                               "(relative to maxB, `r`) matches that of the fitted non-linear model (`rTarget`), `rSim`",
                               "(the `r` that was reached) and `clamped` (`TRUE` if `rTarget` is outside the range of the",
                               "factorial for that `growthcurve`), and `nlmMaxB`, the fitted model's maximum biomass in",
                               "the fitting window."))
  )
))

## event types
#   - type `init` is required for initialization

doEvent.Biomass_speciesParameters = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      ## No tree species in this study area (sppEquiv has no rows, established by fireSense_ELFs):
      ## there are no species parameters to estimate. Biomass_borealDataPrep already supplied the
      ## empty species / speciesEcoregion tables, so leave them and schedule nothing.
      if (is.data.frame(sim$sppEquiv) && nrow(sim$sppEquiv) == 0L) {
        message("Biomass_speciesParameters: no tree species in this study area; nothing to estimate")
        return(invisible(sim))
      }

      ## estimate the traits
      sim <- Init(sim)

      ## update tables
      sim <- updateSpeciesTables(sim)
    },
    warning(
      paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
            "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = "")
    )
  )
  return(invisible(sim))
}

## event functions
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
Init <- function(sim) {
  origDTthreads <- data.table::getDTthreads()
  if (getDTthreads() > P(sim)$.useParallel) {
    data.table::setDTthreads(P(sim)$.useParallel)
  }
  on.exit(data.table::setDTthreads(origDTthreads))

  ## outputs are always created so the `createsOutput` contract holds; populated below when applicable
  sim$speciesGrowthCurvesLandis <- data.table::data.table()
  sim$speciesGrowthCurvesPSP <- data.table::data.table()
  sim$speciesGrowthcurveScore <- data.table::data.table()
  sim$speciesTraitMatch <- data.table::data.table()

  ## if no PSP data supplied, simList returned unchanged
  if (all(P(sim)$PSPdataTypes == "none")) {
    message("P(sim)$PSPdataTypes is 'none' -- bypassing species traits estimation from PSP data.")
    return(sim)
  }
  if (is.na(P(sim)$sppEquivCol)) {
    stop("Please supply 'sppEquivCol' in parameters of Biomass_speciesParameters.")
  }
  paramCheckOtherMods(sim, "maxBInFactorial")
  paramCheckOtherMods(sim, paramToCheck = "sppEquivCol", ifSetButDifferent = "error")

  ## 1. the single-species factorial ---------------------------------------------------------------
  ## A local, not a `mod$` object, so it is freed when `Init` returns.
  factorial <- loadSingleSpeciesFactorial(sim)
  fcs <- factorialCurveStats(factorial, P(sim)$standAgesForFitting, P(sim)$maxBInFactorial)

  ## 2. non-linear models of biomass by age, one per species ---------------------------------------
  biomassKey <- unique(sim$sppEquivLong[, .(PSP, Latin_full)])
  setnames(biomassKey, old = "PSP", new = "SpBiomassEq")

  message("preparing PSPs for growth curves")
  psp <- prepPSPaNPP(studyAreaANPP = sim$studyAreaANPP, PSPperiod = P(sim)$PSPperiod,
                     PSPgis =  sim$PSPgis_sppParams, PSPmeasure = sim$PSPmeasure_sppParams,
                     PSPplot = sim$PSPplot_sppParams, useHeight = P(sim)$useHeight,
                     biomassModel = P(sim)$biomassModel, minDBH = P(sim)$minDBH,
                     sppEquivLong = biomassKey, sppEquiv = sim$sppEquiv,
                     mergeHybridSprucePSP = P(sim)$mergeHybridSprucePSP,
                     excludeBECzones = P(sim)$excludeBECzonesHybridSpruce,
                     BECzones = sim$BECzonesBC) |>
    Cache(userTags = c(currentModule(sim), "prepPSPaNPP"))

  message("building growth curves") # this cache call takes several minutes to process..
  sim$speciesGrowthCurves <- buildGrowthCurves(PSPdata = psp,
                                               speciesCol = P(sim)$sppEquivCol,
                                               sppEquiv = sim$sppEquiv,
                                               minimumSampleSize = P(sim)$minimumPlots,
                                               quantileAgeSubset = P(sim)$quantileAgeSubset,
                                               speciesFittingApproach = "single") |>
    Cache(userTags = c(currentModule(sim), "buildGrowthCurves"))

  noData <- !vapply(sim$speciesGrowthCurves, is.list, logical(1))
  if (any(noData)) {
    message(cli::col_yellow("Insufficient data to estimate species parameters for ",
                            paste(names(noData)[noData], collapse = ", "),
                            " - will keep original user-supplied parameters"))
  }
  fitted <- sim$speciesGrowthCurves[!noData]
  nlms <- lapply(fitted, function(x) x$NonLinearModel[[1]])

  ## 3. mANPPproportion and inflationFactor that make the factorial grow like the fitted models -----
  ## (for every growthcurve of the factorial; see `matchManppProportion()`)
  facAges <- sort(unique(factorial$cohortData$age))
  matched <- rbindlist(lapply(names(nlms), function(sp) {
    f <- fcs[targetSpecies == sp]
    if (!nrow(f)) {
      message(cli::col_yellow(sp, " is not in the factorial - will keep original user-supplied parameters"))
      return(NULL)
    }
    ages <- comparisonAges(facAges, P(sim)$standAgesForFitting, f$longevity[1])
    ns <- nlmCurveStats(nlms[[sp]], ages)
    if (is.na(ns$r) || ns$r <= 0) {
      message(cli::col_yellow(sp, ": the fitted curve has no biomass increase in the fitting window",
                              " - will keep original user-supplied parameters"))
      return(NULL)
    }
    m <- matchManppProportion(f, ns$r)
    if (!nrow(m)) {
      return(NULL)
    }
    m[, `:=`(species = sp, nlmMaxB = ns$maxB, rTarget = ns$r)]
  }))
  if (!NROW(matched)) {
    message("No species could be matched to the factorial; species traits are not changed.")
    return(sim)
  }
  setcolorder(matched, "species")
  sim$speciesTraitMatch <- matched
  fittedSpecies <- unique(matched$species)

  ## 4. growthcurve ---------------------------------------------------------------------------------
  startGC <- setNames(sim$species$growthcurve[match(fittedSpecies, sim$species$species)], fittedSpecies)
  startGC[is.na(startGC)] <- median(matched$growthcurve)
  finalGC <- startGC

  canCalibrate <- isTRUE(P(sim)$calibrateGrowthcurve)
  if (canCalibrate && is.null(sim$ecoregionMap)) {
    message("No `ecoregionMap`: `growthcurve` is not calibrated.")
    canCalibrate <- FALSE
  }
  if (canCalibrate) {
    cal <- calibrateGrowthcurves(
      psp = psp, sppEquiv = sim$sppEquiv, sppEquivCol = P(sim)$sppEquivCol, PSPgis = sim$PSPgis_sppParams,
      ecoregionMap = sim$ecoregionMap, speciesEcoregion = sim$speciesEcoregion,
      speciesTable = sim$species, matched = matched[, .(species, growthcurve, mANPPproportion, inflationFactor)],
      initialGrowthcurve = startGC, iterations = P(sim)$calibrationIterations,
      minInterval = P(sim)$calibrationMinInterval, maxInterval = P(sim)$calibrationMaxInterval,
      maxUnits = P(sim)$maxCalibrationPlots, focalMinProportion = P(sim)$calibrationMinFocalProportion,
      paths = paths(sim), corePath = findCoreModule(modulePath(sim))
    )
    if (!is.null(cal)) {
      finalGC[names(cal$growthcurve)] <- cal$growthcurve
      sim$speciesGrowthcurveScore <- cal$score
      gcPlot <- ggplot(cal$score, aes(growthcurve, dLL)) +
        geom_line() + geom_point() +
        facet_wrap(~ focal, scales = "free_y") +
        geom_vline(data = data.table(focal = names(finalGC), gc = finalGC), aes(xintercept = gc), linetype = "dashed") +
        labs(x = "growthcurve", y = "log-likelihood relative to the best", title = "Calibration of growthcurve") +
        theme_bw()
      suppressWarnings(Plots(gcPlot, usePlot = FALSE, fn = print, ggsaveArgs = list(width = 10, height = 7),
                             filename = "growthcurve_calibration"))
    }
  }

  ## 5. species table -------------------------------------------------------------------------------
  final <- rbindlist(lapply(fittedSpecies, function(sp) {
    cbind(species = sp, traitsAtGrowthcurve(matched[species == sp], finalGC[[sp]]))
  }))
  message("Species traits:")
  print(final)
  species <- copy(sim$species)
  if (!"inflationFactor" %in% names(species)) {
    species[, inflationFactor := NA_real_]
  }
  idx <- match(final$species, species$species)
  set(species, idx, "growthcurve", round(final$growthcurve, 3))
  set(species, idx, "mANPPproportion", round(final$mANPPproportion, 3))
  set(species, idx, "inflationFactor", round(final$inflationFactor, 3))
  sim$species <- species

  ## 6. the fitted and the LandR curves ----------------------------------------------------------------
  curves <- landrVsNlmCurves(factorial, final, nlms, P(sim)$standAgesForFitting)
  originalData <- rbindlist(lapply(fitted, function(x) x$originalData), idcol = "species")
  gg <- ggplot(curves, aes(x = standAge, colour = species)) +
    geom_point(data = originalData, aes(standAge, biomass), size = 0.25, alpha = 0.3) +
    geom_line(aes(y = BscaledNonLinear, linetype = "LandR"), linewidth = 1.5) +
    geom_line(aes(y = predNonLinear, linetype = "Non-Linear"), linewidth = 1.5) +
    scale_linetype_manual(name = "Curve", values = c("LandR" = "solid", "Non-Linear" = "twodash")) +
    facet_wrap(~ species, scales = "free_y") +
    coord_cartesian(xlim = c(0, max(P(sim)$standAgesForFitting) * 1.5)) +
    labs(y = "Biomass", x = "Stand Age", colour = "Species") +
    ggtitle("Comparing best LandR and Non-linear curves") +
    theme_bw()
  suppressWarnings(Plots(gg, usePlot = FALSE, fn = print, ggsaveArgs = list(width = 10, height = 7),
                         filename = "LandR_VS_NLM_growthCurves"))

  if (isTRUE(P(sim)$landis)) {
    ## LANDIS mode: expose the fitted LANDIS-version growth curves (BscaledNonLinear by species and
    ## standAge) and the PSP observations behind them.
    sim$speciesGrowthCurvesLandis <- curves
    sim$speciesGrowthCurvesPSP <- originalData
  }
  return(sim)
}

#' Single-species factorial tables from the paths in the simList
loadSingleSpeciesFactorial <- function(sim) {
  fmt <- "feather" ## faster for small-med data compared to parquet
  openFactorial <- function(path, nm) {
    if (tools::file_ext(path) == "rds") {
      ## using rds default from .inputObjects
      prepInputs(targetFile = basename(path), destinationPath = inputPath(sim), fun = "readRDS",
                 useCache = FALSE, url = extractURL(nm, sim))
    } else {
      arrow::open_dataset(path, format = fmt)
    }
  }
  readSingleSpeciesFactorial(openFactorial(sim$cohortDataFactorial_path, "cohortDataFactorial_path"),
                             openFactorial(sim$speciesTableFactorial_path, "speciesTableFactorial_path"))
}

updateSpeciesTables <- function(sim) {
  modifiedTables <- modifySpeciesAndSpeciesEcoregionTable(
    speciesEcoregion = sim$speciesEcoregion,
    speciesTable = sim$species
  )
  sim$speciesEcoregion <- modifiedTables$newSpeciesEcoregion
  sim$species <- modifiedTables$newSpeciesTable
  return(sim)
}

### template for save events
Save <- function(sim) {
  sim <- saveFiles(sim)
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  origDTthreads <- data.table::getDTthreads()
  if (origDTthreads > 4) {
    data.table::setDTthreads(4)
  }
  on.exit(data.table::setDTthreads(origDTthreads))

  cacheTags <- c(currentModule(sim), "OtherFunction:.inputObjects")
  dPath <- asPath(inputPath(sim), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  if (!suppliedElsewhere("cohortDataFactorial_path", sim)) {
    sim$cohortDataFactorial_path <- file.path(dPath, "cohortDataFactorial_medium.rds")
  }

  if (!suppliedElsewhere("speciesTableFactorial_path", sim)) {
    sim$speciesTableFactorial_path <- file.path(dPath, "speciesTableFactorial_medium.rds")
  }

  if (!suppliedElsewhere("sppEquiv", sim)) {
    ## pass a default sppEquivalencies_CA for common species in western Canada
    sppEquiv <- LandR::sppEquivalencies_CA
    ## sppEquiv should be Boreal to match default P(sim)$sppEquivCol
    sim$sppEquiv <- sppEquiv[Boreal %in% c("Pice_Mar", "Pice_Gla", "Pinu_Con",
                                           "Pinu_Ban", "Popu_Tre", "Lari_Lar",
                                           "Betu_Pap", "Abie_Bal"), ]
  }

  if (!suppliedElsewhere("sppEquivLong", sim)) {
    ## pass a default sppEquivalencies_CA for common species in western Canada
    sim$sppEquivLong <- LandR::sppEquivalencies_CA
  }

  if (!suppliedElsewhere("speciesEcoregion", sim)) {
    warning("generating dummy speciesEcoregion data - run Biomass_borealDataPrep for table with real speciesEcoregion attributes")
    sim$speciesEcoregion <- data.table(
      speciesCode = unique(sim$sppEquiv[[P(sim)$sppEquivCol]]),
      ecoregionGroup = "x",
      establishprob = 0.5,
      maxB = P(sim)$maxBInFactorial,
      maxANPP = P(sim)$maxBInFactorial / 30,
      year = 0
    )
  }

  ## check parameter consistency across modules
  paramCheckOtherMods(sim, "sppEquivCol", ifSetButDifferent = "error")

  if (!suppliedElsewhere("species", sim)) {
    message("generating dummy species data - run Biomass_borealDataPrep for table with real species attributes")
    speciesTable <- getSpeciesTable()
    sim$species <- prepSpeciesTable(
      speciesTable,
      sppEquiv = sim$sppEquiv,
      sppEquivCol = P(sim)$sppEquivCol
    )
  }

  if (!suppliedElsewhere("PSPmeasure_sppParams", sim) ||
      !suppliedElsewhere("PSPplot_sppParams", sim) ||
      !suppliedElsewhere("PSPgis_sppParams", sim)) {
    message("one or more PSP objects not supplied. Generating PSP data...")

    PSPdata <- getPSP(
      PSPdataTypes = P(sim)$PSPdataTypes,
      destinationPath = dPath,
      forGMCS = FALSE
    ) |>
      Cache(userTags = c(cacheTags, P(sim)$PSPdataTypes))

    sim$PSPmeasure_sppParams <- PSPdata$PSPmeasure
    sim$PSPplot_sppParams <- PSPdata$PSPplot
    sim$PSPgis_sppParams <- PSPdata$PSPgis
  }

  return(invisible(sim))
}
