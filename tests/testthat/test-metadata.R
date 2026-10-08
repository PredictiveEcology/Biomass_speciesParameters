## The module's metadata is its public contract: a project using this module binds
## to these object names and classes. Renaming or retyping one breaks every caller,
## which is exactly the class of change the raster -> terra migration makes, so it is
## worth asserting here rather than discovering downstream.
##
## When a change is deliberate, update this file in the same commit and bump the
## module version to match: removed, renamed or retyped is a MAJOR bump.

test_that("module metadata parses", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_type(md, "list")
  expect_identical(md$name, moduleName)
})

test_that("inputs are the expected names and classes", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  inputs <- stats::setNames(md$inputObjects$objectClass, md$inputObjects$objectName)
  expect_identical(
    inputs[order(names(inputs))],
    c(BECzonesBC                 = "sf",
      cohortDataFactorial_path   = "fs_path",
      PSPgis_sppParams           = "sf",
      PSPmeasure_sppParams       = "data.table",
      PSPplot_sppParams          = "data.table",
      species                    = "data.table",
      speciesEcoregion           = "data.table",
      speciesTableFactorial_path = "fs_path",
      sppEquiv                   = "data.table",
      sppEquivLong               = "data.table",
      studyAreaANPP              = "sf")
  )
})

test_that("outputs are the expected names and classes", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  outputs <- stats::setNames(md$outputObjects$objectClass, md$outputObjects$objectName)
  expect_identical(
    outputs[order(names(outputs))],
    c(species                   = "data.table",
      speciesBalanceCheck       = "data.table",
      speciesEcoregion          = "data.table",
      speciesGrowthCurves       = "list",
      speciesGrowthCurvesLandis = "data.table",
      speciesGrowthCurvesPSP    = "data.table")
  )
})

test_that("parameters are the expected names", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_identical(
    sort(md$parameters$paramName),
    sort(c(".plotInitialTime", ".plotInterval", ".plots", ".saveInitialTime",
           ".saveInterval", ".studyAreaName", ".useCache", ".useParallel",
           "balanceGrowth", "biomassModel", "excludeBECzonesHybridSpruce", "landis",
           "maxBInFactorial", "mergeHybridSprucePSP", "minDBH", "minimumPlots",
           "PSPdataTypes", "PSPperiod", "quantileAgeSubset", "sharedGrowthcurve",
           "speciesFittingApproach", "sppEquivCol", "standAgesForFitting", "targetK",
           "useHeight"))
  )
})

test_that("balanceGrowth is off by default, with sharedGrowthcurve and targetK unset", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  defaults <- stats::setNames(md$parameters$default, md$parameters$paramName)
  expect_false(defaults[["balanceGrowth"]])
  expect_true(is.na(defaults[["sharedGrowthcurve"]]))
  expect_true(is.na(defaults[["targetK"]]))
})
