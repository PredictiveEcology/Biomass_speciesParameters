## The two factorial tables (~8 GB + ~1 GB at full size) are read inside `Init` and used nowhere
## else. Until 3.0.2 they were stored as `mod$` objects, which kept them in the simList's `.modObjs`
## (and in every `Copy(sim)`) for the whole simulation. This runs the module on two tiny feather
## tables and checks that nothing is left behind in `.modObjs`.
test_that("factorial tables do not stay in the module's .modObjs after init", {
  skip_on_cran()
  skip_if_not_installed("arrow")
  library(SpaDES.core)
  library(data.table)

  opts <- options(spades.useRequire = FALSE, spades.moduleCodeChecks = FALSE,
                  reproducible.useCache = FALSE, reproducible.useMemoise = FALSE)
  on.exit(options(opts), add = TRUE)

  tmp <- file.path(tempdir(), "bsp-modObjs")
  dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  ## tiny stand-ins for the factorial tables, same columns `Init` touches
  cohortDataFactorial <- data.table(pixelGroup = c(1L, 1L, 2L), speciesCode = c("A", "A", "B"),
                                    age = c(1L, 2L, 1L), B = c(100, 200, 150))
  speciesTableFactorial <- data.table(pixelGroup = c(1L, 2L), speciesCode = c("A", "B"),
                                      longevity = c(100L, 150L), growthcurve = c(0.5, 0.6),
                                      mortalityshape = c(15L, 20L), mANPPproportion = c(3.3, 3.5))
  cdfPath <- file.path(tmp, "cohortDataFactorial.feather")
  stfPath <- file.path(tmp, "speciesTableFactorial.feather")
  arrow::write_feather(cohortDataFactorial, cdfPath)
  arrow::write_feather(speciesTableFactorial, stfPath)

  ## enough of the other inputs that `.inputObjects` downloads nothing and `updateSpeciesTables` runs
  sppEquiv <- data.table(LandR = c("Pice_Mar", "Popu_Tre"), Boreal = c("Pice_Mar", "Popu_Tre"),
                         PSP = c("Pice_Mar", "Popu_Tre"), Latin_full = c("Picea mariana", "Populus tremuloides"))
  species <- data.table(species = c("Pice_Mar", "Popu_Tre"), longevity = c(250L, 150L),
                        growthcurve = c(0.5, 0.5), mortalityshape = c(15L, 15L),
                        hardsoft = c("soft", "hard"), inflationFactor = c(1, 1), mANPPproportion = c(3.3, 3.3))
  speciesEcoregion <- data.table(speciesCode = c("Pice_Mar", "Popu_Tre"), ecoregionGroup = "x",
                                 establishprob = 0.5, maxB = 5000L, maxANPP = 166L, year = 0)

  ## modulePath must contain a directory literally named Biomass_speciesParameters; a checkout
  ## under another name (a git worktree) gets a symlink
  modRoot <- normalizePath(file.path("..", ".."))
  if (identical(basename(modRoot), "Biomass_speciesParameters")) {
    modPath <- dirname(modRoot)
  } else {
    modPath <- file.path(tmp, "modules")
    dir.create(modPath, showWarnings = FALSE)
    file.symlink(modRoot, file.path(modPath, "Biomass_speciesParameters"))
  }

  sim <- simInitAndSpades(
    times = list(start = 0, end = 1),
    modules = "Biomass_speciesParameters",
    params = list(Biomass_speciesParameters = list(PSPdataTypes = "none", .plots = NA,
                                                   .useCache = FALSE)),
    objects = list(cohortDataFactorial_path = cdfPath, speciesTableFactorial_path = stfPath,
                   sppEquiv = sppEquiv, sppEquivLong = sppEquiv, species = species,
                   speciesEcoregion = speciesEcoregion,
                   PSPmeasure_sppParams = data.table(), PSPplot_sppParams = data.table(),
                   PSPgis_sppParams = data.table()),
    paths = list(modulePath = modPath, outputPath = file.path(tmp, "outputs"),
                 cachePath = file.path(tmp, "cache"), inputPath = tmp),
    debug = FALSE
  )

  expect_s4_class(sim, "simList")
  modObjs <- sim@.xData$.modObjs$Biomass_speciesParameters
  expect_null(modObjs$cohortDataFactorial)
  expect_null(modObjs$speciesTableFactorial)
})
