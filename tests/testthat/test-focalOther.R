## `buildGrowthCurves`: under speciesFittingApproach = "focal", the species other than the focal one
## that are above the 20% dominance threshold in a plot are pooled into "Other".

## synthetic PSP plots, one measurement each: Pice_mar 40%, Pinu_ban 30% and Popu_tre 30% of the plot biomass,
## two trees per species; biomass follows a Chapman-Richards curve of stand age, so the fits converge
synthPSPfocal <- function(nPlots = 40) {
  ages <- seq(20, 160, length.out = nPlots)
  plotB <- 15000 * (1 - exp(-0.03 * ages))^2
  data.table::rbindlist(lapply(seq_len(nPlots), function(i) {
    data.table::data.table(
      MeasureID = i, OrigPlotID1 = paste0("p", i), MeasureYear = 2000L, standAge = ages[i], PlotSize = 1,
      Latin_full = rep(c("Picea mariana", "Pinus banksiana", "Populus tremuloides"), each = 2),
      biomass = plotB[i] * rep(c(0.4, 0.3, 0.3), each = 2) / 2,
      SpBiomassEq = NA, source = "x", PSP = NA, Elevation = NA)
  }))
}

sppEquivFocal <- data.table::data.table(
  Latin_full = c("Picea mariana", "Pinus banksiana", "Populus tremuloides"),
  LandR = c("Pice_mar", "Pinu_ban", "Popu_tre"))

growthCurves <- function(psp, approach) {
  ## buildModels() calls data.table functions unqualified, as they are attached when the module runs
  withr::local_package("data.table")
  withr::local_seed(2718)
  suppressWarnings(suppressMessages(
    buildGrowthCurves(psp, "LandR", sppEquivFocal, quantileAgeSubset = 99, minimumSampleSize = 1,
                      speciesFittingApproach = approach)))
}

test_that("regression: focal 'Other' is the summed biomass of the non-focal species, one row per plot", {
  psp <- synthPSPfocal()
  gcs <- growthCurves(psp, "focal")
  plotB <- psp[, .(plotB = sum(biomass)), by = "OrigPlotID1"]
  od <- gcs$Pice_mar$originalData
  other <- od[speciesTemp == "Other"]
  focal <- od[speciesTemp == "Pice_mar"]
  ## one "Other" row per plot measurement (it counts once towards minimumPlots and the weights)
  expect_equal(nrow(other), nrow(focal))
  expect_false(anyDuplicated(other[, c("OrigPlotID1", "MeasureYear")]) > 0)
  ## Pinu_ban + Popu_tre = 60% of the plot, not the 30% of one of them
  other <- other[plotB, on = "OrigPlotID1", nomatch = NULL]
  expect_equal(other$biomass, 0.6 * other$plotB)
  expect_equal(other$spDom, rep(0.6, nrow(other)))
  ## the focal species is unchanged
  focal <- focal[plotB, on = "OrigPlotID1", nomatch = NULL]
  expect_equal(focal$biomass, 0.4 * focal$plotB)
  ## so "Other" has the higher maximum biomass, and the focal species is matched to the second cohort
  maxB <- od[, .(maxB = max(biomass)), by = "speciesTemp"]
  expect_gt(maxB[speciesTemp == "Other"]$maxB, maxB[speciesTemp == "Pice_mar"]$maxB)
})

test_that("pairwise and single keep each species' own biomass", {
  psp <- synthPSPfocal()
  ## pairwise uses plots with exactly two species above 20%: drop Popu_tre from half the plots
  pspPair <- psp[!(Latin_full == "Populus tremuloides" & MeasureID %% 2 == 0)]
  pair <- growthCurves(pspPair, "pairwise")
  od <- pair$Pice_mar__Pinu_ban$originalData
  plotB <- pspPair[, .(plotB = sum(biomass)), by = "OrigPlotID1"]
  od <- od[plotB, on = "OrigPlotID1", nomatch = NULL]
  expect_setequal(as.character(unique(od$speciesTemp)), c("Pice_mar", "Pinu_ban"))
  expect_equal(od[speciesTemp == "Pice_mar"]$biomass, 0.4 / 0.7 * od[speciesTemp == "Pice_mar"]$plotB)
  expect_equal(od[speciesTemp == "Pinu_ban"]$biomass, 0.3 / 0.7 * od[speciesTemp == "Pinu_ban"]$plotB)

  ## single uses plots where one species is above 50%: make Pice_mar 80% of the plot
  pspSingle <- data.table::copy(psp)[Latin_full == "Picea mariana", biomass := biomass * 6]
  single <- growthCurves(pspSingle, "single")
  od <- single$Pice_mar$originalData
  plotB <- pspSingle[, .(plotB = sum(biomass)), by = "OrigPlotID1"]
  od <- od[plotB, on = "OrigPlotID1", nomatch = NULL]
  expect_identical(as.character(unique(od$speciesTemp)), "Pice_mar")
  expect_equal(od$biomass, 0.8 * od$plotB)
})
