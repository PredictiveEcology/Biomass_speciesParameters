## When no species of a study area gets a growth curve (every species below `minimumPlots`, or
## every fit failed) there are no traits to estimate. The run used to stop in modifySpeciesTable(),
## whose limitToSpeciesLongevity() joined an empty trait table on `species` ("non-existing
## column(s): cols[1]='species'"). It now warns and ends as with PSPdataTypes = "none": the species
## traits are those given as input.

## Run the module on tiny inputs: the helper's factorial tables, and PSP records of 3 Alberta plots
## where Pice_mar and Popu_tre are both co-dominant, i.e. 3 plots each, below minimumPlots 10.
runOnTinyInputs <- function(PSPdataTypes, tmp) {
  inp <- makeModifySpeciesTableInputs()
  cdfPath <- file.path(tmp, "cohortDataFactorial.feather")
  stfPath <- file.path(tmp, "speciesTableFactorial.feather")
  arrow::write_feather(inp$factorialBiomass, cdfPath)
  arrow::write_feather(inp$factorialTraits, stfPath)

  plots <- data.table::data.table(MeasureID = 1:3, OrigPlotID1 = paste0("p", 1:3),
                                  MeasureYear = 2000L, source = "AB", baseYear = 2000L,
                                  baseSA = c(40, 60, 80), PlotSize = 0.04)
  ## 25 black spruce and 15 aspen per plot (the forest filter needs >= 30 trees)
  trees <- plots[, .(TreeNumber = 1:40,
                     Species = rep(c("Picea mariana", "Populus tremuloides"), c(25, 15)),
                     DBH = 20, Height = 15, status = "A", first_tree_year = 2000L,
                     last_tree_year = 2000L, diff_dbh = 0),
                 by = c("MeasureID", "OrigPlotID1", "MeasureYear", "source")]
  gis <- sf::st_as_sf(data.frame(OrigPlotID1 = plots$OrigPlotID1, lon = -115, lat = 55),
                      coords = c("lon", "lat"), crs = 4326)
  sppEquiv <- data.table::data.table(LandR = c("Pice_mar", "Popu_tre"),
                                     Latin_full = c("Picea mariana", "Populus tremuloides"),
                                     PSP = c("black spruce", "trembling aspen"))
  species <- data.table::data.table(species = c("Pice_mar", "Popu_tre"),
                                    longevity = c(250L, 150L), growthcurve = c(0.5, 0.6),
                                    mortalityshape = c(15L, 20L), hardsoft = c("soft", "hard"))
  speciesEcoregion <- data.table::data.table(speciesCode = c("Pice_mar", "Popu_tre"),
                                             ecoregionGroup = "x", establishprob = 0.5,
                                             maxB = 5000L, maxANPP = 166L, year = 0)

  ## modulePath must contain a directory literally named Biomass_speciesParameters
  modRoot <- normalizePath(file.path("..", ".."))
  modPath <- dirname(modRoot)
  if (!identical(basename(modRoot), "Biomass_speciesParameters")) {
    modPath <- file.path(tmp, "modules")
    dir.create(modPath, showWarnings = FALSE)
    file.symlink(modRoot, file.path(modPath, "Biomass_speciesParameters"))
  }

  SpaDES.core::simInitAndSpades(
    times = list(start = 0, end = 1),
    modules = "Biomass_speciesParameters",
    params = list(Biomass_speciesParameters = list(PSPdataTypes = PSPdataTypes, minimumPlots = 10,
                                                   .plots = NA, .studyAreaName = "tinyArea",
                                                   .useCache = FALSE)),
    objects = list(cohortDataFactorial_path = cdfPath, speciesTableFactorial_path = stfPath,
                   sppEquiv = sppEquiv, sppEquivLong = data.table::copy(sppEquiv),
                   species = species, speciesEcoregion = speciesEcoregion,
                   PSPmeasure_sppParams = trees, PSPplot_sppParams = plots,
                   PSPgis_sppParams = gis),
    paths = list(modulePath = modPath, outputPath = file.path(tmp, "outputs"),
                 cachePath = file.path(tmp, "cache"), inputPath = tmp),
    debug = FALSE
  )
}

test_that("regression: a run where no species gets a growth curve ends as with PSPdataTypes none", {
  skip_on_cran()
  skip_if_not_installed("arrow")
  attachModuleDeps()
  withr::local_package("SpaDES.core")
  withr::local_options(spades.useRequire = FALSE, spades.moduleCodeChecks = FALSE,
                       reproducible.useCache = FALSE, reproducible.useMemoise = FALSE)

  w <- capture_warnings(simPSP <- suppressMessages(runOnTinyInputs("all", withr::local_tempdir())))
  simNone <- suppressWarnings(suppressMessages(runOnTinyInputs("none", withr::local_tempdir())))

  expect_match(w, "No growth curve could be fitted from the PSP data in tinyArea", fixed = TRUE,
               all = FALSE)
  expect_match(w, "fewer than minimumPlots (10) PSP plots: Pice_mar, Popu_tre", fixed = TRUE,
               all = FALSE)
  expect_identical(simPSP$species, simNone$species)
  expect_identical(simPSP$speciesEcoregion, simNone$speciesEcoregion)
  ## the input traits, not ones estimated from PSP data
  expect_identical(simPSP$species$growthcurve, c(0.5, 0.6))
  expect_identical(simPSP$species$mortalityshape, c(15L, 20L))
})

test_that("isFittedGC() is TRUE only for a growth curve fitted for every species of the fit", {
  fit <- structure(list(), class = "nlrob")
  failed <- structure("Error : failed\n", class = "try-error")
  expect_true(isFittedGC(list(NonLinearModel = list(Pice_mar = fit))))
  expect_false(isFittedGC(c(Pice_mar = "insufficient data")))
  expect_false(isFittedGC(failed))
  expect_false(isFittedGC(list(NonLinearModel = list(Pice_mar = failed))))
  expect_false(isFittedGC(list(NonLinearModel = list())))
  ## "pairwise": both species of the pair
  expect_false(isFittedGC(list(NonLinearModel = list(Abie_las = fit, Pice_eng = failed))))
})

test_that("the warning names each species and why it has no growth curve", {
  failed <- structure("Error : failed\n", class = "try-error")
  GCs <- list(Pice_mar = c(Pice_mar = "insufficient data"),
              Pinu_ban = list(NonLinearModel = list(Pinu_ban = failed)))
  msg <- noFittedCurvesMessage(GCs, minimumPlots = 25, studyAreaName = NA_character_)
  expect_match(msg, "fitted from the PSP data;", fixed = TRUE)
  expect_match(msg, "fewer than minimumPlots (25) PSP plots: Pice_mar;", fixed = TRUE)
  expect_match(msg, "every fit failed: Pinu_ban.", fixed = TRUE)
  expect_match(msg, "as with PSPdataTypes = 'none'", fixed = TRUE)
})
