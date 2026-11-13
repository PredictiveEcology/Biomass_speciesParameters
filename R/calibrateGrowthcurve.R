## Calibrating `growthcurve` against remeasured PSPs.
##
## With `maxANPP` fixed by the fitted biomass-by-age curve (see `matchManppProportion()`), `growthcurve` is
## the one trait left to estimate. Single-species curves cannot separate it from the rest, but competition
## can: each remeasured plot is simulated in Biomass_core from its first measurement, once per candidate
## `growthcurve` of the focal species, and the `growthcurve` whose simulated biomass is closest to the second
## measurement wins.

#' Biomass of each species in each PSP measurement
#'
#' @param psp PSP tree data from `prepPSPaNPP()`.
#' @return `data.table` (`OrigPlotID1`, `MeasureYear`, `species`, `standAge`, `B`), `B` in g/m2, for species in `sppEquiv`
#' @noRd
speciesPlotBiomass <- function(psp, sppEquiv, sppEquivCol) {
  d <- copy(psp)
  d[, species := LandR::equivalentName(value = Latin_full, df = sppEquiv, column = sppEquivCol)]
  d <- d[!is.na(species)]
  d[, areaAdjustedB := biomass / PlotSize]
  d[, .(B = sum(areaAdjustedB), standAge = standAge[1]), by = c("OrigPlotID1", "MeasureYear", "species")]
}

#' Pairs of consecutive PSP measurements to simulate
#'
#' @param spBiomass Output of `speciesPlotBiomass()`.
#' @param minInterval,maxInterval Range of years between measurements to keep.
#' @param maxBiomassDrop Largest fraction of total plot biomass lost between measurements; larger
#'   losses are treated as disturbance, which Biomass_core is not asked to reproduce here.
#' @param maxUnits Random subsample size.
#' @return `data.table` (`unit`, `OrigPlotID1`, `year1`, `year2`, `interval`, `age1`)
#' @noRd
remeasurementUnits <- function(spBiomass, minInterval = 3, maxInterval = 30, maxBiomassDrop = 0.5,
                               maxUnits = 2000L) {
  plots <- spBiomass[, .(B = sum(B), standAge = standAge[1]), by = c("OrigPlotID1", "MeasureYear")]
  setorderv(plots, c("OrigPlotID1", "MeasureYear"))
  units <- plots[, .(year1 = MeasureYear[-.N], year2 = MeasureYear[-1], age1 = standAge[-.N],
                     B1 = B[-.N], B2 = B[-1]), by = "OrigPlotID1"]
  units <- units[year2 - year1 >= minInterval & year2 - year1 <= maxInterval &
                   B2 >= (1 - maxBiomassDrop) * B1]
  if (nrow(units) > maxUnits) {
    units <- units[sample.int(nrow(units), maxUnits)]
  }
  units[, `:=`(unit = .I, interval = year2 - year1)]
  units[, c("B1", "B2") := NULL]
  units[]
}

#' Ecoregion group of each PSP plot
#'
#' @param PSPgis `sf` points with `OrigPlotID1`.
#' @param ecoregionMap Categorical `SpatRaster` with an `ecoregionGroup` level (Biomass_borealDataPrep).
#' @return `data.table` (`OrigPlotID1`, `ecoregionGroup`); plots outside the map get `NA`.
#' @noRd
plotEcoregions <- function(PSPgis, ecoregionMap) {
  pts <- sf::st_transform(PSPgis, terra::crs(ecoregionMap))
  vals <- terra::extract(ecoregionMap, terra::vect(pts), ID = FALSE)[[1]]
  if (is.numeric(vals)) {
    lev <- as.data.table(terra::cats(ecoregionMap)[[1]])
    vals <- lev$ecoregionGroup[match(vals, lev$ID)]
  }
  data.table(OrigPlotID1 = pts$OrigPlotID1, ecoregionGroup = as.character(vals))
}

#' Per-pixelGroup simulation input for the calibration
#'
#' One pixelGroup for each remeasurement unit, each focal species in it, and each candidate
#' `growthcurve`. The focal species has the candidate `growthcurve` (and the `mANPPproportion` and
#' `inflationFactor` that go with it); every other species has its current traits.
#'
#' @param units Output of `remeasurementUnits()`.
#' @param spBiomass Output of `speciesPlotBiomass()`.
#' @param plotEco Output of `plotEcoregions()`.
#' @param speciesEcoregion `maxB` by `speciesCode` and `ecoregionGroup`, as made by Biomass_borealDataPrep.
#' @param speciesTable `species`, `longevity`, `mortalityshape`.
#' @param matched `data.table` (`species`, `growthcurve`, `mANPPproportion`, `inflationFactor`) of the candidates,
#'   several `growthcurve` per species.
#' @param current `data.table` (`species`, `growthcurve`, `mANPPproportion`, `inflationFactor`): traits of the
#'   species that compete with the focal species.
#' @param focalMinProportion The smallest share of plot biomass a species needs to be a focal species of the plot.
#' @param minB Cohorts smaller than this (g/m2) at the first measurement are dropped.
#'
#' @return A list: `pixelGroups` (`pixelGroup`, `unit`, `focal`, `growthcurve`, `ecoregionGroup`, `interval`,
#'   `maxB`, `observedB`), `cohortData`, `species` (the Biomass_core species table, one row per variant) and
#'   `speciesEcoregion`.
#' @noRd
calibrationInput <- function(units, spBiomass, plotEco, speciesEcoregion, speciesTable, matched, current,
                             focalMinProportion = 0.2, minB = 10) {
  se <- as.data.table(speciesEcoregion)
  se[, speciesCode := as.character(speciesCode)]
  se[, ecoregionGroup := as.character(ecoregionGroup)]
  if ("year" %in% names(se)) se <- se[year == min(year)]
  maxBnominal <- unique(se[, .(species = speciesCode, ecoregionGroup, maxB)], by = c("species", "ecoregionGroup"))

  ## species of each unit at both measurements; a species absent at the second one died
  m1 <- merge(units, spBiomass, by.x = c("OrigPlotID1", "year1"), by.y = c("OrigPlotID1", "MeasureYear"))
  m1 <- merge(m1, plotEco, by = "OrigPlotID1")
  m1 <- m1[!is.na(ecoregionGroup)]
  m1 <- merge(m1, speciesTable[, .(species, longevity)], by = "species")
  m1 <- merge(m1, maxBnominal, by = c("species", "ecoregionGroup"))
  m1 <- m1[B >= minB & age1 < longevity]
  m1[, share := B / sum(B), by = "unit"]
  m2 <- spBiomass[, .(OrigPlotID1, year2 = MeasureYear, species, observedB = B)]
  m1 <- merge(m1, m2, by = c("OrigPlotID1", "year2", "species"), all.x = TRUE)
  m1[is.na(observedB), observedB := 0]

  focal <- m1[share >= focalMinProportion & species %in% matched$species,
              .(unit, focal = species, ecoregionGroup, interval, maxB, observedB)]
  cands <- matched[, .(focal = species, growthcurve)]
  pg <- merge(focal, cands, by = "focal", allow.cartesian = TRUE)
  setorderv(pg, c("unit", "focal", "growthcurve"))
  pg[, pixelGroup := .I]

  ## variants: the focal species at each candidate growthcurve, and every species at its current traits
  traits <- rbind(
    matched[, .(variant = paste0(species, "__g", growthcurve), species, growthcurve, mANPPproportion, inflationFactor)],
    current[, .(variant = paste0(species, "__cur"), species, growthcurve, mANPPproportion, inflationFactor)]
  )
  sp <- merge(traits, speciesTable[, .(species, longevity, mortalityshape)], by = "species")

  ## cohorts: all species of the unit, the focal one as its candidate variant
  cohorts <- merge(pg[, .(pixelGroup, unit, focal, growthcurve, ecoregionGroup)],
                   m1[, .(unit, species, age = age1, B)], by = "unit", allow.cartesian = TRUE)
  cohorts[, variant := fifelse(species == focal, paste0(species, "__g", growthcurve), paste0(species, "__cur"))]
  cohorts <- cohorts[variant %in% sp$variant]
  cohortData <- cohorts[, .(pixelGroup, ecoregionGroup, speciesCode = variant, age = as.integer(age),
                            B = as.integer(round(B)))]

  ## maxB, maxANPP of each variant in each ecoregion that is used
  seOut <- CJ(variant = unique(cohortData$speciesCode), ecoregionGroup = unique(pg$ecoregionGroup))
  seOut <- merge(seOut, sp[, .(variant, species, mANPPproportion, inflationFactor)], by = "variant")
  seOut <- merge(seOut, maxBnominal, by = c("species", "ecoregionGroup"))
  seOut[, maxB := as.integer(round(maxB * inflationFactor))]
  seOut[, maxANPP := as.integer(round(maxB * mANPPproportion / 100))]
  seOut <- seOut[, .(speciesCode = variant, ecoregionGroup, establishprob = 0.5, maxB, maxANPP, year = 0)]

  list(pixelGroups = pg[, .(pixelGroup, unit, focal, growthcurve, ecoregionGroup, interval, maxB, observedB)],
       cohortData = cohortData,
       species = sp[variant %in% cohortData$speciesCode, .(species = variant, longevity, growthcurve,
                                                           mortalityshape, mANPPproportion)],
       speciesEcoregion = seOut)
}

#' Run Biomass_core on the calibration input
#'
#' The mechanics follow `Biomass_speciesFactorial::RunExperiment()`: pixelGroups are laid out as cells of a
#' raster, Biomass_core runs annually with no seeding, and `cohortData` is saved at the remeasurement intervals.
#'
#' @param input Output of `calibrationInput()`.
#' @param paths As `SpaDES.core::paths()`; `outputPath` receives the saved `cohortData`.
#' @param corePath Directory containing the `Biomass_core` module.
#' @param minCohortBiomass,initialB Passed to Biomass_core.
#'
#' @return `data.table` (`pixelGroup`, `speciesCode`, `age`, `B`, `saveTime`).
#' @noRd
runCoreExperiment <- function(input, paths, corePath, minCohortBiomass = 9, initialB = 10) {
  pg <- input$pixelGroups
  egLevels <- sort(unique(pg$ecoregionGroup))
  egFactor <- function(x) factor(x, levels = egLevels)

  nPG <- max(pg$pixelGroup)
  nr <- ceiling(sqrt(nPG))
  nc <- ceiling(nPG / nr)
  pixelGroupMap <- terra::rast(nrows = nr, ncols = nc, xmin = 0, xmax = nc, ymin = 0, ymax = nr)
  terra::values(pixelGroupMap) <- c(seq_len(nPG), rep(NA_integer_, nr * nc - nPG))
  ecoregionMap <- terra::rast(pixelGroupMap)
  terra::values(ecoregionMap) <- c(as.integer(egFactor(pg$ecoregionGroup[order(pg$pixelGroup)])),
                                   rep(NA_integer_, nr * nc - nPG))
  levels(ecoregionMap) <- data.frame(ID = seq_along(egLevels), ecoregion = seq_along(egLevels),
                                     ecoregionGroup = egFactor(egLevels))
  studyArea <- terra::as.polygons(terra::ext(pixelGroupMap), crs = "")
  terra::crs(studyArea) <- terra::crs(pixelGroupMap)

  cohortData <- copy(input$cohortData)
  cohortData[, `:=`(speciesCode = factor(speciesCode, levels = input$species$species),
                    ecoregionGroup = egFactor(ecoregionGroup))]
  setcolorder(cohortData, c("speciesCode", "pixelGroup", "ecoregionGroup", "age", "B"))
  speciesEcoregion <- copy(input$speciesEcoregion)
  speciesEcoregion[, `:=`(speciesCode = factor(speciesCode, levels = input$species$species),
                          ecoregionGroup = egFactor(ecoregionGroup))]

  species <- copy(input$species)
  species[, c("sexualmature", "seeddistance_eff", "seeddistance_max", "resproutprob", "resproutage_min",
              "resproutage_max", "postfireregen", "leaflongevity", "wooddecayrate", "leafLignin", "hardsoft",
              "Area", "firetolerance", "shadetolerance") :=
            list(30L, 0L, 0L, 0.5, 0L, 0L, factor("none"), 3L, 0.07, 0.1, factor("soft"), factor("BP"), 3L, 1L)]
  species[, mANPPproportion := NULL]

  ecoregion <- data.table(ecoregionGroup = egFactor(egLevels), active = "yes")
  minRelativeB <- data.table(ecoregionGroup = egFactor(egLevels),
                             LandR::minRelativeBDefaults()[rep(1L, length(egLevels))])

  intervals <- sort(unique(pg$interval))
  outputs <- data.frame(objectName = "cohortData", saveTime = intervals, eventPriority = 1,
                        fun = "qs2::qs_save", stringsAsFactors = FALSE)

  objects <- list(
    cohortData = cohortData, ecoregion = ecoregion, ecoregionMap = ecoregionMap,
    minRelativeB = minRelativeB, pixelGroupMap = pixelGroupMap, species = species,
    speciesEcoregion = speciesEcoregion, speciesLayers = "species",
    sppColorVect = setNames(viridis::viridis(nrow(species)), species$species),
    sppEquiv = data.table(B_factorial = species$species), sppNameVector = species$species,
    studyArea = studyArea, rasterToMatch = pixelGroupMap
  )
  parameters <- list(Biomass_core = list(
    .maxMemory = 1e9, .plots = NULL, .saveInitialTime = NA, .saveInterval = NA, .useParallel = 1,
    .useCache = NULL, calcSummaryBGM = NULL, initialBiomassSource = "cohortData",
    seedingAlgorithm = "noSeeding", sppEquivCol = "B_factorial", successionTimestep = 10,
    vegLeadingProportion = 0, minCohortBiomass = minCohortBiomass, initialB = initialB
  ))

  opts <- options(LandR.assertions = FALSE, LandR.verbose = 0, spades.moduleCodeChecks = FALSE,
                  spades.recoveryMode = FALSE)
  on.exit(options(opts), add = TRUE)

  coreModule <- basename(corePath)
  paths$modulePath <- dirname(corePath)
  out <- SpaDES.core::simInitAndSpades(
    times = list(start = 0, end = max(intervals)), params = parameters, modules = coreModule,
    paths = paths, objects = objects, outputs = outputs, debug = 1, outputObjects = "pixelGroupMap"
  )
  of <- as.data.table(SpaDES.core::outputs(out))[objectName == "cohortData"]
  rbindlist(lapply(seq_len(nrow(of)), function(i) {
    cd <- as.data.table(qs2::qs_read(of$file[i]))[, .(pixelGroup, speciesCode, age, B)]
    cd[, saveTime := of$saveTime[i]]
  }))
}

#' Score each candidate `growthcurve` of each focal species
#'
#' @param pixelGroups `input$pixelGroups` from `calibrationInput()`.
#' @param simulated Output of `runCoreExperiment()`.
#' @return `data.table` (`focal`, `growthcurve`, `n`, `mse`, `ll`, `dLL`): `mse` is the mean squared error of the
#'   focal species' biomass at the second measurement, relative to its `maxB`; `ll` the Gaussian profile
#'   log-likelihood `-n/2 * log(mse)`; `dLL` is `ll` minus the species' best `ll`.
#' @noRd
scoreCalibration <- function(pixelGroups, simulated) {
  sim <- simulated[, .(pixelGroup, speciesCode = as.character(speciesCode), saveTime, B)]
  pg <- copy(pixelGroups)
  pg[, variant := paste0(focal, "__g", growthcurve)]
  pg <- merge(pg, sim, by.x = c("pixelGroup", "variant", "interval"), by.y = c("pixelGroup", "speciesCode", "saveTime"),
              all.x = TRUE)
  pg[is.na(B), B := 0]
  score <- pg[, .(n = .N, mse = mean(((B - observedB) / maxB)^2)), by = c("focal", "growthcurve")]
  score[, ll := -n / 2 * log(pmax(mse, 1e-12))]
  score[, dLL := ll - max(ll), by = "focal"]
  setorderv(score, c("focal", "growthcurve"))
  score[]
}

#' The `growthcurve` that best fits, with a parabola through the best score and its neighbours
#'
#' @param score `score` for a single species from `scoreCalibration()`.
#' @return The estimated `growthcurve`; the grid value with the best score if it is at the edge of the grid
#'   or the neighbouring scores do not bend upwards.
#' @noRd
bestGrowthcurve <- function(score) {
  score <- score[order(growthcurve)]
  i <- which.max(score$ll)
  if (i == 1L || i == nrow(score)) {
    return(score$growthcurve[i])
  }
  x <- score$growthcurve[(i - 1):(i + 1)]
  y <- score$mse[(i - 1):(i + 1)]
  fit <- coef(lm(y ~ x + I(x^2)))
  if (is.na(fit[3]) || fit[3] <= 0) {
    return(score$growthcurve[i])
  }
  min(max(-fit[2] / (2 * fit[3]), x[1]), x[3])
}

#' Find the Biomass_core module to run the calibration with
#' @param modulePath The module path of the simulation.
#' @return Path to a directory called `Biomass_core`.
#' @noRd
findCoreModule <- function(modulePath) {
  candidates <- c(file.path(modulePath, "Biomass_core"),
                  file.path(modulePath, "Biomass_speciesFactorial", "submodules", "Biomass_core"))
  found <- candidates[file.exists(file.path(candidates, "Biomass_core.R"))]
  if (!length(found)) {
    stop("Cannot find the Biomass_core module in ", paste(dirname(candidates), collapse = " or "),
         ". Biomass_speciesParameters needs it to calibrate `growthcurve`.")
  }
  found[1]
}

#' Calibrate `growthcurve` of each species with remeasured PSPs
#'
#' Iterates `calibrationInput()`, `runCoreExperiment()`, `scoreCalibration()`: competitors are simulated with
#' the traits found by the previous pass, so that each pass improves on the preceding one.
#'
#' @param psp Output of `prepPSPaNPP()`.
#' @param matched `data.table` (`species`, `growthcurve`, `mANPPproportion`, `inflationFactor`) from
#'   `matchManppProportion()` for all species.
#' @param initialGrowthcurve Named vector, the starting `growthcurve` of each species in `matched`.
#' @param iterations Number of passes.
#' @inheritParams calibrationInput
#' @inheritParams remeasurementUnits
#' @inheritParams plotEcoregions
#' @inheritParams runCoreExperiment
#'
#' @return A list: `growthcurve` (named vector), `score` (the last pass' `scoreCalibration()`, with `pass`), or
#'   `NULL` if there were no usable remeasurements.
#' @noRd
calibrateGrowthcurves <- function(psp, sppEquiv, sppEquivCol, PSPgis, ecoregionMap, speciesEcoregion, speciesTable,
                                  matched, initialGrowthcurve, iterations = 2L, minInterval = 3, maxInterval = 30,
                                  maxUnits = 500L, focalMinProportion = 0.2, paths, corePath) {
  spBiomass <- speciesPlotBiomass(psp, sppEquiv, sppEquivCol)
  plotEco <- plotEcoregions(PSPgis, ecoregionMap)
  units <- remeasurementUnits(spBiomass, minInterval = minInterval, maxInterval = maxInterval,
                              maxUnits = maxUnits)
  gcNow <- initialGrowthcurve
  scores <- NULL
  for (pass in seq_len(iterations)) {
    current <- rbindlist(lapply(names(gcNow), function(sp) {
      cbind(species = sp, traitsAtGrowthcurve(matched[species == sp], gcNow[[sp]]))
    }))
    input <- calibrationInput(units, spBiomass, plotEco, speciesEcoregion, speciesTable, matched, current,
                              focalMinProportion = focalMinProportion)
    if (!nrow(input$pixelGroups)) {
      message("No remeasured PSPs could be matched to the species and ecoregions; `growthcurve` not calibrated.")
      return(NULL)
    }
    message("Calibrating growthcurve, pass ", pass, " of ", iterations, ": ", max(input$pixelGroups$pixelGroup),
            " simulated pixelGroups")
    simulated <- runCoreExperiment(input, paths = paths, corePath = corePath)
    scores <- scoreCalibration(input$pixelGroups, simulated)
    scores[, pass := pass]
    best <- scores[, .(growthcurve = bestGrowthcurve(.SD)), by = "focal"]
    gcNow[best$focal] <- best$growthcurve
  }
  list(growthcurve = gcNow, score = scores)
}
