## `balanceGrowth`: share growthcurve across species and set mANPPproportion so that
## K = mANPPproportion * maxB^(1 - growthcurve) is equal (see R/balanceGrowth.R).

## likelihood table of one species: a bowl around (gc0, m0), `slope` log-likelihood units per unit of
## mANPPproportion and twice that per unit of growthcurve; kinks on grid points, so linear interpolation is exact
bowl <- function(species, gc0, m0, slope = 1) {
  g <- data.table::CJ(growthcurve = seq(0.6, 0.8, by = 0.02), mANPPproportion = seq(2, 8, by = 0.5))
  g[, `:=`(species = species, llNonLinDelta = slope * abs(mANPPproportion - m0) + 2 * slope * abs(growthcurve - gc0))]
}

fitted <- data.table::data.table(species = c("A", "B", "C"), growthcurve = c(0.70, 0.72, 0.74),
                     mANPPproportion = c(5, 4, 3))
maxB <- c(A = 1000, B = 4000, C = 9000, D = 500)
ll <- data.table::rbindlist(list(bowl("A", 0.70, 5), bowl("B", 0.72, 4), bowl("C", 0.74, 3)))

test_that("balanceGrowthK shares growthcurve and equalises K", {
  chk <- suppressMessages(balanceGrowthK(fitted, maxB, ll))
  expect_equal(chk$growthcurveNew, rep(0.72, 3))             # median of the fitted values
  expect_equal(chk$growthcurveFitted, fitted$growthcurve)
  kFitted <- fitted$mANPPproportion * maxB[fitted$species]^(1 - 0.72)
  expect_equal(chk$K, rep(median(kFitted), 3), tolerance = 1e-3)       # equal, up to the 3 dp rounding
  expect_equal(chk$mANPPproportionNew,
               unname(round(median(kFitted) / maxB[fitted$species]^(1 - 0.72), 3)))
  expect_equal(chk$maxB, unname(maxB[c("A", "B", "C")]))
})

test_that("balanceGrowthK rounds to 3 decimals, and honours sharedGrowthcurve and targetK", {
  ## (outside the grid, so deltaLL is NA and there is a warning)
  chk <- suppressWarnings(suppressMessages(balanceGrowthK(fitted, maxB, ll, sharedGrowthcurve = 0.7234, targetK = 80)))
  expect_equal(chk$growthcurveNew, rep(0.723, 3))
  expect_equal(chk$mANPPproportionNew, round(chk$mANPPproportionNew, 3))
  expect_equal(chk$mANPPproportionNew, round(80 / unname(maxB[c("A", "B", "C")])^(1 - 0.723), 3))
  expect_equal(chk$K, rep(80, 3), tolerance = 1e-3)
})

test_that("applyBalancedTraits leaves unfitted species untouched and reports them", {
  chk <- suppressMessages(balanceGrowthK(fitted, maxB, ll))
  sp <- data.table::data.table(species = c("A", "D", "B", "C"), growthcurve = c(1, 2, 3, 4),
                   mANPPproportion = c(9, 8, 7, 6), longevity = c(100L, 200L, 300L, 400L))
  expect_message(out <- applyBalancedTraits(sp, chk), "not fitted.*D")
  expect_equal(out[species == "D"], sp[species == "D"])
  expect_equal(out$longevity, sp$longevity)
  expect_equal(out$growthcurve[out$species != "D"], rep(0.72, 3))
  expect_equal(out[match(chk$species, species)]$mANPPproportion, chk$mANPPproportionNew)
})

test_that("regression: a fitted species missing from speciesEcoregion is left out, not an error", {
  ## stopped TSA04, TSA08, ... with "no maxB in speciesEcoregion for fitted species Betu_pap"
  noC <- maxB[c("A", "B", "D")]
  expect_message(chk <- suppressWarnings(balanceGrowthK(fitted, noC, ll)), "not in speciesEcoregion.*C")
  expect_equal(chk$species, c("A", "B"))
  expect_equal(chk$growthcurveNew, rep(0.71, 2))             # median over A and B only
  sp <- data.table::data.table(species = c("A", "B", "C"), growthcurve = c(1, 2, 3), mANPPproportion = c(9, 8, 7))
  out <- suppressMessages(applyBalancedTraits(sp, chk))
  expect_equal(out[species == "C"], sp[species == "C"])
})

test_that("balanceGrowthK says what it did for each species", {
  expect_message(balanceGrowthK(fitted, maxB, ll), "balanceGrowth B: growthcurve 0.720 -> 0.720, mANPPproportion 4.000 -> ")
})

test_that("deltaLL is the loss in log-likelihood against the species' best, and flags species far from it", {
  ## C is sharply peaked at its own fit, so moving it to the shared growthcurve costs more than 2 units
  llSharp <- data.table::rbindlist(list(bowl("A", 0.70, 5), bowl("B", 0.72, 4), bowl("C", 0.74, 3, slope = 100)))
  w <- tryCatch(suppressMessages(balanceGrowthK(fitted, maxB, llSharp)), warning = function(w) conditionMessage(w))
  expect_match(w, "C \\(")                                   # names the species ...
  expect_no_match(w, "A \\(|B \\(")                        # ... and only that one
  chk <- suppressWarnings(suppressMessages(balanceGrowthK(fitted, maxB, llSharp)))
  expect_true(chk[species == "C", deltaLL] > 2)
  expect_true(all(chk[species != "C", deltaLL] < 2))
  ## the loss is the bowl evaluated at the new point (interpolation is exact for this table)
  expect_equal(chk[species == "A", deltaLL],
               abs(chk[species == "A", mANPPproportionNew] - 5) + 2 * abs(0.72 - 0.70))
  ## no warning when everything is within 2 units
  expect_no_warning(suppressMessages(balanceGrowthK(fitted, maxB, ll)))
})

test_that("a point outside the factorial grid has deltaLL NA and a warning", {
  expect_warning(chk <- suppressMessages(balanceGrowthK(fitted, maxB, ll, targetK = 5000)), "outside the factorial grid")
  expect_true(all(is.na(chk$deltaLL)))
})

test_that("medianMaxB takes the median over each species' speciesEcoregion rows", {
  se <- data.table::data.table(speciesCode = rep(c("A", "B"), c(3, 2)), ecoregionGroup = c(1:3, 1:2),
                   maxB = c(100, 300, 200, 50, 70))
  expect_equal(medianMaxB(se), c(A = 200, B = 60))
})

test_that("regression: medianMaxB works on integer maxB with odd and even row counts per species", {
  ## median() of an integer is an integer for an odd count and a double for an even one, which
  ## data.table refused across groups in box D
  se <- data.table::data.table(speciesCode = rep(c("A", "B"), c(3, 2)), maxB = c(100L, 300L, 200L, 50L, 71L))
  expect_equal(medianMaxB(se), c(A = 200, B = 60.5))
})

test_that("modifySpeciesTable with default settings gives the unchanged species table", {
  skip_if_not_installed("LandR")
  skip_if_not_installed("purrr")
  skip_if_not_installed("ggplot2")
  attachModuleDeps()
  inputs <- makeModifySpeciesTableInputs()
  res <- suppressMessages(suppressWarnings(do.call(modifySpeciesTable, inputs)))
  ## captured from the module before `balanceGrowth` existed (development at ba70f74); Betu_pap has no fit
  expect_equal(as.data.frame(res$best),
               data.frame(species = c("Abie_las", "Pice_eng", "Betu_pap"), longevity = c(225L, 300L, 150L),
                          growthcurve = c(0.7, 0.7, 0.9), mortalityshape = c(20L, 20L, 15L),
                          mANPPproportion = c(5.398, 3.808, 9), inflationFactor = c(1.044, 1.038, 1)))
  ## the unrounded fits behind it are what `balanceGrowth` starts from
  expect_equal(round(res$bestUnrounded$growthcurve, 2), res$best$growthcurve[1:2])
  expect_equal(round(res$bestUnrounded$mANPPproportion, 3), res$best$mANPPproportion[1:2])
})

test_that("balancing the fits of modifySpeciesTable equalises K and checks it against the likelihood table", {
  skip_if_not_installed("LandR")
  skip_if_not_installed("purrr")
  skip_if_not_installed("ggplot2")
  attachModuleDeps()
  res <- suppressMessages(suppressWarnings(do.call(modifySpeciesTable, makeModifySpeciesTableInputs())))
  expect_setequal(res$llTable$species, c("Abie_las", "Pice_eng"))
  expect_true(nrow(res$llTable) > 0)
  chk <- suppressMessages(suppressWarnings(
    balanceGrowthK(res$bestUnrounded, c(Abie_las = 18000, Pice_eng = 12000), res$llTable)))
  expect_equal(chk$K[1], chk$K[2], tolerance = 1e-3)
  expect_equal(chk$growthcurveNew, rep(0.7, 2))
  expect_false(anyNA(chk$deltaLL))
})
