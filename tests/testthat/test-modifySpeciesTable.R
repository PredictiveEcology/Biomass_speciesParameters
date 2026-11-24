## `longevity` is an input to the fit, not a result of it: it picks the factorial curves a species'
## traits are fitted from, those at the grid longevities either side of its own
## (limitToSpeciesLongevity()). modifySpeciesTable() used to write the AIC-weighted mean
## of those grid longevities over the species' own, so a fitted species ended up with a value between
## its two grid neighbours, or with the grid's largest longevity if it lived longer than that. The
## helper's grid has longevities 150 and 300: the old code gave Abie_las 225 (not 250) and Pice_eng 300
## (not 450).
test_that("regression: modifySpeciesTable keeps each species' own longevity", {
  skip_if_not_installed("LandR")
  skip_if_not_installed("purrr")
  skip_if_not_installed("ggplot2")
  attachModuleDeps()
  inputs <- makeModifySpeciesTableInputs()
  res <- suppressMessages(suppressWarnings(do.call(modifySpeciesTable, inputs)))
  expect_equal(res$best$species, inputs$speciesTable$species)
  ## both species were fitted (growthcurve moved off the input 0.9), so the check below is not vacuous
  expect_true(all(res$best$growthcurve[1:2] != inputs$speciesTable$growthcurve[1:2]))
  expect_equal(res$best$longevity, inputs$speciesTable$longevity)
})
