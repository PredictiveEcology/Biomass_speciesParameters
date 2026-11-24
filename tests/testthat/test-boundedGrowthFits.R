## buildModels() fits the growth curves within bounds and rejects a fit whose optimum sits on a bound.
## When a species' youngest real stands are already at the biomass plateau, the curve's rate is not
## identified: unbounded, the fit failed or converged by luck of the random starts and jitter, so
## whether the species got its own traits depended on the RNG.

## focal PSP table as buildGrowthCurves() hands it to buildModels(): the focal species' plots only
synthFocalPSP <- function(age, biomass, species = "Pice_gla") {
  data.table::data.table(speciesTemp = species, MeasureYear = 2000L,
                         OrigPlotID1 = paste0("p", seq_along(age)), PlotSize = 0.04,
                         standAge = age, spPlotBiomass = biomass, spDom = 0.8)
}
fitFor <- function(psp, seed) {
  ## a simulation attaches the module's reqdPkgs; simulateYoungStands() calls data.table() unqualified
  withr::local_package("data.table")
  set.seed(seed)
  gc <- suppressMessages(buildModels("Pice_gla", psp, speciesEquiv = NULL, sppCol = "LandR",
                                     minSize = 10, q = 99))
  gc$NonLinearModel[["Pice_gla"]]
}

test_that("a gradual growth curve is fitted, with every parameter inside its bounds", {
  set.seed(42)
  age <- rep(seq(10, 200, by = 5), each = 2)
  biomass <- 20000 * (1 - exp(-0.02 * age))^2 * exp(rnorm(length(age), 0, 0.1))
  psp <- synthFocalPSP(age, biomass)
  for (s in 1:5) {
    fit <- fitFor(psp, s)
    expect_false(inherits(fit, "try-error"))
    expect_s3_class(fit, "nlrob")
  }
})

test_that("a curve whose plateau is reached by the youngest stands is not fitted, for any seed", {
  set.seed(42)
  age <- c(31, rep(seq(45, 210, by = 5), each = 2))
  biomass <- 5000 * exp(rnorm(length(age), 0, 0.3))
  psp <- synthFocalPSP(age, biomass)
  for (s in 1:5) {
    fit <- fitFor(psp, s)
    expect_s3_class(fit, "try-error")
    expect_match(as.character(fit), "not identifiable")
  }
})

test_that("the fit does not depend on the RNG seed", {
  set.seed(42)
  age <- c(rep(c(0, 5, 8), 3), rep(seq(10, 200, by = 5), each = 2))
  biomass <- 20000 * (1 - exp(-0.02 * age))^2 * exp(rnorm(length(age), 0, 0.1))
  psp <- synthFocalPSP(age, biomass)
  cf <- lapply(c(1, 2, 3), function(s) coef(fitFor(psp, s)))
  expect_equal(cf[[2]], cf[[1]], tolerance = 1e-4)
  expect_equal(cf[[3]], cf[[1]], tolerance = 1e-4)
})

test_that("buildModels() leaves the caller's RNG state as it found it", {
  withr::local_package("data.table")
  set.seed(42)
  age <- rep(seq(10, 200, by = 5), each = 2)
  psp <- synthFocalPSP(age, 20000 * (1 - exp(-0.02 * age))^2 * exp(rnorm(length(age), 0, 0.1)))
  set.seed(7)
  expected <- runif(3)
  set.seed(7)
  suppressMessages(buildModels("Pice_gla", psp, speciesEquiv = NULL, sppCol = "LandR", minSize = 10, q = 99))
  expect_identical(runif(3), expected)
})

test_that("paramsOnBound() names the parameters within tolerance of a bound", {
  lower <- c(A = 0, k = 0.0001, p = 1)
  upper <- c(A = 10, k = 0.13, p = 80)
  expect_identical(paramsOnBound(c(A = 5, k = 0.05, p = 3), lower, upper), character(0))
  expect_identical(paramsOnBound(c(A = 5, k = 0.13, p = 3), lower, upper), "k")
  expect_identical(paramsOnBound(c(A = 0, k = 0.05, p = 80), lower, upper), c("A", "p"))
})
