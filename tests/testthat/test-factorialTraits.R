## A synthetic single-species factorial: biomass is a logistic curve whose speed rises with
## `mANPPproportion` up to 6 and falls beyond it (as it does in Biomass_core, where mortality
## scales with maxANPP), and whose scale falls with `growthcurve`.
synthFactorial <- function() {
  grid <- data.table::CJ(growthcurve = c(0.4, 0.8), mANPPproportion = 1:10)
  grid[, pixelGroup := .I]
  grid[, `:=`(species = paste0("A", pixelGroup), targetSpecies = "Pice_mar", longevity = 200L,
              mortalityshape = 20L)]
  ages <- seq(1L, 191L, by = 10L)
  cd <- grid[, {
    speed <- 0.02 * (1 + 0.5 * min(mANPPproportion, 12 - mANPPproportion)) * (1 + growthcurve / 4)
    data.table::data.table(speciesCode = species, age = ages, B = 5000 / (1 + exp(-speed * (ages - 60))))
  }, by = "pixelGroup"]
  list(cohortData = cd[, .(speciesCode, age, B, pixelGroup)],
       speciesTable = grid[, .(species, targetSpecies, longevity, growthcurve, mortalityshape,
                               mANPPproportion, pixelGroup)])
}

test_that("curveStats gives the net peak increment relative to the maximum", {
  cs <- curveStats(c(10, 20, 30), c(100, 300, 400))
  expect_equal(cs$maxB, 400)
  expect_equal(cs$peakInc, 20)
  expect_equal(cs$r, 0.05)
  expect_true(is.na(curveStats(10, 100)$r))
})

test_that("comparisonAges keeps the fitting window and ages up to longevity", {
  expect_equal(comparisonAges(seq(1, 301, 10), c(21, 91), 60), c(21, 31, 41, 51))
})

test_that("factorialCurveStats gives r over the window and the inflation factor", {
  fs <- factorialCurveStats(synthFactorial(), c(21, 91), 5000)
  expect_equal(nrow(fs), 20L)
  expect_false(anyNA(fs$r))
  ## the factorial's maxB is 5000 and the curve reaches (almost) all of it
  expect_true(all(fs$inflationFactor >= 1 & fs$inflationFactor < 1.1))
})

test_that("matchManppProportion recovers mANPPproportion on the rising branch only", {
  fs <- factorialCurveStats(synthFactorial(), c(21, 91), 5000)
  one <- fs[growthcurve == 0.4 & mANPPproportion == 4]
  m <- matchManppProportion(fs, one$r)
  expect_equal(m[growthcurve == 0.4, mANPPproportion], 4, tolerance = 1e-6)
  expect_false(m[growthcurve == 0.4, clamped])
  ## beyond the largest r of the factorial the traits at the largest r are returned and flagged
  big <- matchManppProportion(fs, max(fs$r) * 2)
  expect_true(all(big$clamped))
  expect_true(all(big$mANPPproportion <= 6))
})

test_that("traitsAtGrowthcurve interpolates and clamps", {
  m <- data.table::data.table(growthcurve = c(0.4, 0.8), mANPPproportion = c(2, 4), inflationFactor = c(1.2, 1))
  expect_equal(traitsAtGrowthcurve(m, 0.6)$mANPPproportion, 3)
  expect_equal(traitsAtGrowthcurve(m, 0.9)$inflationFactor, 1)
})

test_that("readSingleSpeciesFactorial needs a factorial made with factorialDesign = 'perSpecies'", {
  f <- synthFactorial()
  expect_error(readSingleSpeciesFactorial(f$cohortData, f$speciesTable[, !"targetSpecies"]), "perSpecies")
  expect_named(readSingleSpeciesFactorial(f$cohortData, f$speciesTable), c("cohortData", "speciesTable"))
})

test_that("remeasurementUnits drops disturbed plots and bad intervals", {
  sb <- data.table::data.table(
    OrigPlotID1 = rep(c("a", "b", "c"), each = 3), MeasureYear = rep(c(2000, 2010, 2011), 3),
    species = "Pice_mar", standAge = 50,
    B = c(1000, 1200, 1250, 1000, 300, 320, 1000, 1100, 1150))
  u <- remeasurementUnits(sb, minInterval = 3, maxInterval = 30)
  ## plot b lost 70% of its biomass; the 1-year interval is too short
  expect_setequal(u$OrigPlotID1, c("a", "c"))
  expect_true(all(u$interval == 10))
})

test_that("bestGrowthcurve refines between grid values and stays on the grid at its edge", {
  gcs <- seq(0.3, 0.7, 0.1)
  mse <- (gcs - 0.52)^2 + 0.01
  score <- data.table::data.table(growthcurve = gcs, mse = mse, ll = -10 * log(mse))
  expect_equal(bestGrowthcurve(score), 0.52, tolerance = 1e-6)
  score[, `:=`(mse = (growthcurve - 0.1)^2 + 0.01)][, ll := -10 * log(mse)]
  expect_equal(bestGrowthcurve(score), 0.3)
})
