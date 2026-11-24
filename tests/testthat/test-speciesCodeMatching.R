## buildGrowthCurves() picks the PSP plots each growth curve is fitted to. Species codes used to be
## matched as substrings (grepl), so a code that is part of another (Pice_eng in Pice_eng_gla) drew
## that species' plots, and a pair whose name is part of another pair's drew that pair's plots.

test_that("focal: a species code that is part of another does not draw that species' plots", {
  ## regression: matched as a substring, Pice_eng also drew the Pice_eng_gla plots, as "Other"
  withr::local_package("data.table")
  ## the PSP records each species' fit is given
  local_mocked_bindings(buildModels = function(species, psp, ...) psp)
  age <- seq(20, 100, by = 20)
  psp <- rbind(pspRecords("Pice_eng", paste0("eng", 1:5), 2000L, age, 100 * age),
               pspRecords("Pice_eng_gla", paste0("gla", 1:5), 2000L, age, 100 * age),
               pspRecords(c("Pice_eng", "Pice_eng_gla"), "mixed", 2000L, 60, c(4000, 3000)))
  out <- suppressMessages(buildGrowthCurves(
    psp, speciesCol = "LandR", sppEquiv = equivFor(c("Pice_eng", "Pice_eng_gla")),
    quantileAgeSubset = 99, minimumSampleSize = 1, speciesFittingApproach = "focal"))
  expect_setequal(unique(out$Pice_eng$OrigPlotID1), c(paste0("eng", 1:5), "mixed"))
  expect_setequal(unique(out$Pice_eng_gla$OrigPlotID1), c(paste0("gla", 1:5), "mixed"))
  ## the other species is "Other" only on the plot it shares with the focal species
  expect_identical(unique(out$Pice_eng$OrigPlotID1[out$Pice_eng$speciesTemp == "Other"]), "mixed")
})

test_that("pairwise: a pair whose name is part of another pair's does not draw that pair's plots", {
  ## regression: matched as a substring, "Abie_las__Pice_eng" also drew the plots of
  ## "Abie_las__Pice_eng_gla"
  withr::local_package("data.table")
  local_mocked_bindings(buildModels = function(species, psp, ...) psp)
  ## Abie_las is recorded first on every plot, so the pairs are named "Abie_las__<other>"
  stands <- function(other, plots) {
    rbind(pspRecords("Abie_las", plots, 2000L, 50, 5000), pspRecords(other, plots, 2000L, 50, 5000))
  }
  psp <- rbind(stands("Pice_eng", paste0("eng", 1:3)), stands("Pice_eng_gla", paste0("gla", 1:3)))
  out <- suppressMessages(buildGrowthCurves(
    psp, speciesCol = "LandR", sppEquiv = equivFor(c("Abie_las", "Pice_eng", "Pice_eng_gla")),
    quantileAgeSubset = 99, minimumSampleSize = 1, speciesFittingApproach = "pairwise"))
  expect_setequal(names(out), c("Abie_las__Pice_eng", "Abie_las__Pice_eng_gla"))
  expect_setequal(unique(out[["Abie_las__Pice_eng"]]$OrigPlotID1), paste0("eng", 1:3))
  expect_setequal(as.character(unique(out[["Abie_las__Pice_eng"]]$speciesTemp)),
                  c("Abie_las", "Pice_eng"))
  expect_setequal(unique(out[["Abie_las__Pice_eng_gla"]]$OrigPlotID1), paste0("gla", 1:3))
})
