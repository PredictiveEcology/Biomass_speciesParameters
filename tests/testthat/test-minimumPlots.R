## buildModels() fits a growth curve only if it has at least `minimumPlots` plot-years (plot x
## measurement year) to fit it to. The check used to count rows of species x plot-year instead. Under
## "focal", a plot-year where another species was also co-dominant (> 20% of plot biomass) counted
## twice: once for the focal species and once for the other species, pooled as "Other", whose rows
## are not fitted. Under "pairwise", every plot-year counted once per species of the pair.

## biomass along a growth curve the fits pin down at once
growth <- function(age) 20000 * (1 - exp(-0.02 * age))^2 * exp(stats::rnorm(length(age), 0, 0.1))

## Pice_mar on 17 plots, 8 of them measured twice: 25 plot-years. Popu_tre is co-dominant (a third of
## plot biomass) in 17 of them, those on the first 9 plots. The old count gave Pice_mar 25 + 17 = 42
## rows (as for black spruce in one study area: 42 rows, 25 plot-years on 17 plots).
focalPSP <- function() {
  withr::local_seed(1)
  plots <- sprintf("p%02d", 1:17)
  plot <- c(plots, plots[1:8])
  year <- rep(c(2000L, 2010L), c(17L, 8L))
  age <- seq(15, 175, by = 10)[match(plot, plots)] + (year - 2000L)
  biomass <- growth(age)
  coDom <- plot %in% plots[1:9]
  rbind(pspRecords("Pice_mar", plot, year, age, biomass),
        pspRecords("Popu_tre", plot[coDom], year[coDom], age[coDom], biomass[coDom] / 2))
}

fitFocal <- function(psp, minimumPlots) {
  suppressMessages(buildGrowthCurves(psp, speciesCol = "LandR",
                                     sppEquiv = equivFor(c("Pice_mar", "Popu_tre")),
                                     quantileAgeSubset = 99, minimumSampleSize = minimumPlots,
                                     speciesFittingApproach = "focal"))
}

test_that("minimumPlots counts the focal species' plot-years, not the rows of co-dominant 'Other'", {
  ## regression: the old count gave Pice_mar 42 >= 40 rows, so its curve was fitted
  withr::local_package("data.table")
  w <- capture_warnings(gcs <- fitFocal(focalPSP(), minimumPlots = 40))
  expect_identical(unname(gcs$Pice_mar), "insufficient data")
  expect_identical(unname(gcs$Popu_tre), "insufficient data")
  expect_match(w, "Pice_mar: 25 plot-years", fixed = TRUE, all = FALSE)
  expect_match(w, "Popu_tre: 17 plot-years", fixed = TRUE, all = FALSE)
})

test_that("a species with exactly minimumPlots plot-years is fitted, with one fewer it is not", {
  withr::local_package("data.table")
  psp <- focalPSP()
  gcs <- suppressWarnings(fitFocal(psp, minimumPlots = 25))
  expect_named(gcs$Pice_mar, c("originalData", "simData", "NonLinearModel"))
  expect_s3_class(gcs$Pice_mar$NonLinearModel$Pice_mar, "nlrob")
  ## regression: 17 plot-years, but 34 rows with the co-dominant Pice_mar as "Other"
  expect_identical(unname(gcs$Popu_tre), "insufficient data")
  gcs <- suppressWarnings(fitFocal(psp, minimumPlots = 26))
  expect_identical(unname(gcs$Pice_mar), "insufficient data")
})

test_that("a species of sppEquiv without plot-years is not fitted", {
  withr::local_package("data.table")
  w <- capture_warnings(gcs <- suppressMessages(buildGrowthCurves(
    focalPSP(), speciesCol = "LandR", sppEquiv = equivFor(c("Pice_mar", "Popu_tre", "Abie_las")),
    quantileAgeSubset = 99, minimumSampleSize = 40, speciesFittingApproach = "focal")))
  expect_identical(unname(gcs$Abie_las), "insufficient data")
  expect_match(w, "Abie_las: 0 plot-years", fixed = TRUE, all = FALSE)
})

test_that("under pairwise, minimumPlots counts the pair's plot-years once, not once per species", {
  ## regression: the old count gave the pair 2 x 17 = 34 rows, so its curves were fitted
  withr::local_package("data.table")
  withr::local_seed(1)
  age <- seq(15, 175, by = 10)
  biomass <- growth(age)
  plots <- sprintf("p%02d", seq_along(age))
  psp <- rbind(pspRecords("Abie_las", plots, 2000L, age, biomass / 2),
               pspRecords("Pice_eng", plots, 2000L, age, biomass / 2))
  w <- capture_warnings(gcs <- suppressMessages(buildGrowthCurves(
    psp, speciesCol = "LandR", sppEquiv = equivFor(c("Abie_las", "Pice_eng")),
    quantileAgeSubset = 99, minimumSampleSize = 18, speciesFittingApproach = "pairwise")))
  expect_identical(unname(gcs[["Abie_las__Pice_eng"]]), "insufficient data")
  expect_match(w, "Abie_las__Pice_eng: 17 plot-years", fixed = TRUE, all = FALSE)
})

