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

test_that("modifySpeciesTable warns for a species whose longevity is outside the factorial's range", {
  skip_if_not_installed("LandR")
  skip_if_not_installed("purrr")
  skip_if_not_installed("ggplot2")
  attachModuleDeps()
  inputs <- makeModifySpeciesTableInputs()
  ## the helper's grid has longevities 150 and 300: Abie_las (250) is inside it, Pice_eng (450) above
  w <- capture_warnings(suppressMessages(do.call(modifySpeciesTable, inputs)))
  expect_match(w, "outside the factorial's range of longevities: Pice_eng (450; traits fitted at 300)",
               fixed = TRUE, all = FALSE)
  expect_false(any(grepl("Abie_las (", w, fixed = TRUE)))
})

test_that("longevityOutsideGrid() names the species above or below the grid, not those inside or on it", {
  withr::local_package("data.table")
  ## the grid longevities each species' traits were fitted at, as limitToSpeciesLongevity() leaves them:
  ## a between 150 and 300, b above the grid (300 only), c below it (150 only), d on a grid value
  newTraits <- data.table(species = c("a", "a", "b", "c", "d"), longevity = c(150, 300, 300, 150, 300))
  ## e has no fitted traits
  speciesTable <- data.table(species = c("a", "b", "c", "d", "e"), longevity = c(250, 450, 100, 300, 999))
  expect_equal(as.data.frame(longevityOutsideGrid(newTraits, speciesTable)),
               data.frame(species = c("b", "c"), longevity = c(450, 100), gridLongevity = c(300, 150)))
})
