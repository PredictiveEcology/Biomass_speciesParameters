## speciesFittingApproach = "all" failed with "object 'gcSpecies1' not found" from December 2021
## (c387671) on. It was removed; it must now fail at once, with a message that says so.

test_that("speciesFittingApproach = 'all' gives a clear error", {
  expect_error(
    buildGrowthCurves(PSPdata = NULL, speciesCol = "LandR", sppEquiv = NULL, quantileAgeSubset = 99,
                      minimumSampleSize = 50, speciesFittingApproach = "all"),
    "'all' is no longer supported"
  )
})

test_that("an unknown speciesFittingApproach gives a clear error", {
  expect_error(
    buildGrowthCurves(PSPdata = NULL, speciesCol = "LandR", sppEquiv = NULL, quantileAgeSubset = 99,
                      minimumSampleSize = 50, speciesFittingApproach = "pooled"),
    "must be one of 'focal', 'pairwise' or 'single'"
  )
})
