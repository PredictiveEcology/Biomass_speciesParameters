## BC interior (hybrid-zone) spruce plots are recorded as "Picea glauca". When the study area's
## sppEquiv has the hybrid merged into Engelmann (LandR.mergeHybridSpruce = "engelmann"), they match
## no species and prepPSPaNPP drops them. relabelHybridSprucePSP() relabels them as the hybrid,
## except in the BEC zones where white spruce is the true species (BWBS, SWB).

hybridLatin <- "Picea engelmannii x glauca"

## sppEquiv as LandR::speciesInStudyArea() returns it with the hybrid merged into Engelmann
mergedEquiv <- function() {
  data.table::data.table(
    Latin_full = c("Picea engelmannii", hybridLatin, "Abies lasiocarpa"),
    LandR = c("Pice_eng", "Pice_eng", "Abie_las"))
}

## four BC plots (SBS, ICH, BWBS, SWB), one Alberta plot and one plot outside every BEC polygon
synthPSP <- function() {
  plots <- data.table::data.table(
    OrigPlotID1 = c("sbs", "ich", "bwbs", "swb", "ab", "sea"),
    source = c("BC", "BC", "BC", "BC", "AB", "BC"),
    lon = c(-122, -118, -121, -126, -112, -135), lat = c(54, 50, 57, 58, 54, 45))
  gis <- sf::st_as_sf(plots, coords = c("lon", "lat"), crs = 4326)
  measure <- plots[, .(Species = c("Picea glauca", "Abies lasiocarpa"), DBH = c(20, 10)),
                   by = c("OrigPlotID1", "source")]
  list(measure = measure, gis = gis[, "OrigPlotID1"])
}

synthBEC <- function() {
  box <- function(x, y) sf::st_polygon(list(rbind(c(x[1], y[1]), c(x[2], y[1]), c(x[2], y[2]),
                                                  c(x[1], y[2]), c(x[1], y[1]))))
  sf::st_sf(ZONE = c("SBS", "ICH", "BWBS", "SWB", "AB"),
            geometry = sf::st_sfc(box(c(-123, -121.5), c(53, 55)), box(c(-119, -117), c(49, 51)),
                                  box(c(-122, -120), c(56, 58)), box(c(-127, -125), c(57, 59)),
                                  box(c(-113, -111), c(53, 55)), crs = 4326))
}

## the module parameters' defaults, so they are written down only in the module metadata
paramDefault <- function(name) {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  md$parameters$default[[match(name, md$parameters$paramName)]]
}

relabel <- function(sppEquiv = mergedEquiv(), mergeHybridSprucePSP = paramDefault("mergeHybridSprucePSP"),
                    excludeBECzones = paramDefault("excludeBECzonesHybridSpruce"), BECzones = synthBEC()) {
  p <- synthPSP()
  relabelHybridSprucePSP(p$measure, p$gis, sppEquiv, mergeHybridSprucePSP, excludeBECzones, BECzones)
}
species <- function(m, plot) m[OrigPlotID1 == plot & DBH == 20]$Species

test_that("BC Picea glauca outside BWBS/SWB becomes the hybrid, and then Pice_eng", {
  out <- suppressMessages(relabel())
  expect_identical(species(out, "sbs"), hybridLatin)
  expect_identical(species(out, "ich"), hybridLatin)
  ## the existing Latin_full match in prepPSPaNPP/buildGrowthCurves then gives Pice_eng
  expect_identical(LandR::equivalentName(species(out, "sbs"), mergedEquiv(), "LandR"), "Pice_eng")
  ## other species are untouched
  expect_true(all(out[DBH == 10]$Species == "Abies lasiocarpa"))
})

test_that("regression: factor plot IDs, as PSPclean returns them, are matched by name", {
  ## zoneOf[<factor>] indexed by the factor codes, so in box D no plot was relabelled
  p <- synthPSP()
  p$measure[, OrigPlotID1 := factor(OrigPlotID1)]
  p$gis$OrigPlotID1 <- factor(p$gis$OrigPlotID1)
  out <- suppressMessages(relabelHybridSprucePSP(p$measure, p$gis, mergedEquiv(),
                                                 paramDefault("mergeHybridSprucePSP"),
                                                 paramDefault("excludeBECzonesHybridSpruce"), synthBEC()))
  expect_identical(species(out, "sbs"), hybridLatin)
  expect_identical(species(out, "ich"), hybridLatin)
  expect_identical(species(out, "bwbs"), "Picea glauca")
})

test_that("the defaults are the documented ones", {
  expect_identical(paramDefault("excludeBECzonesHybridSpruce"), c("BWBS", "SWB"))
  expect_identical(paramDefault("mergeHybridSprucePSP"), getOption("LandR.mergeHybridSpruce", "engelmann"))
})

test_that("BWBS and SWB plots keep Picea glauca", {
  out <- suppressMessages(relabel())
  expect_identical(species(out, "bwbs"), "Picea glauca")
  expect_identical(species(out, "swb"), "Picea glauca")
  ## the excluded zones are a parameter
  out2 <- suppressMessages(relabel(excludeBECzones = "BWBS"))
  expect_identical(species(out2, "swb"), hybridLatin)
  expect_identical(species(out2, "bwbs"), "Picea glauca")
})

test_that("non-BC Picea glauca and BC plots outside any BEC polygon are not relabelled", {
  out <- suppressMessages(relabel())
  expect_identical(species(out, "ab"), "Picea glauca")
  expect_identical(species(out, "sea"), "Picea glauca")
})

test_that("a message reports the relabelled plots and records by BEC zone", {
  expect_message(relabel(), "SBS")
  expect_message(relabel(), "ICH")
})

test_that("mergeHybridSprucePSP = NA, or 'white', changes nothing", {
  p <- synthPSP()
  expect_identical(suppressMessages(relabel(mergeHybridSprucePSP = NA_character_)), p$measure)
  expect_identical(suppressMessages(relabel(mergeHybridSprucePSP = "white")), p$measure)
})

test_that("sppEquiv without the hybrid merged into Pice_eng changes nothing", {
  p <- synthPSP()
  noHybrid <- mergedEquiv()[Latin_full != hybridLatin]
  unmerged <- data.table::copy(mergedEquiv())[Latin_full == hybridLatin, LandR := "Pice_eng_gla"]
  expect_identical(suppressMessages(relabel(noHybrid)), p$measure)
  expect_identical(suppressMessages(relabel(unmerged)), p$measure)
})

test_that("the input table is not modified by reference", {
  p <- synthPSP()
  before <- data.table::copy(p$measure)
  suppressMessages(relabelHybridSprucePSP(p$measure, p$gis, mergedEquiv(), "engelmann", "BWBS",
                                          synthBEC()))
  expect_identical(p$measure, before)
})

test_that("the default follows LandR.mergeHybridSpruce", {
  withr::local_options(LandR.mergeHybridSpruce = NA_character_)
  expect_identical(suppressMessages(relabel()), synthPSP()$measure)
  withr::local_options(LandR.mergeHybridSpruce = "engelmann")
  expect_identical(species(suppressMessages(relabel()), "sbs"), hybridLatin)
})

test_that("with no BEC polygons available, nothing is relabelled and a warning says so", {
  p <- synthPSP()
  local_mocked_bindings(getBECzonesBC = function(...) NULL)
  expect_warning(out <- relabelHybridSprucePSP(p$measure, p$gis, mergedEquiv(), "engelmann", "BWBS"),
                 "BEC")
  expect_identical(out, p$measure)
})

test_that("regression: a study area without BC plots is unchanged by default", {
  p <- synthPSP()
  nonBC <- p$measure[source != "BC"]
  expect_identical(suppressMessages(relabelHybridSprucePSP(nonBC, p$gis, mergedEquiv(), "engelmann",
                                                           "BWBS", synthBEC())), nonBC)
})

test_that("a relabelled tree keeps the biomass of its own species, only its species grouping changes", {
  ## prepPSPaNPP calls biomassCalculation(), setkey() ... unqualified, as it does attached in a simList
  withr::local_package("pemisc")
  withr::local_package("data.table")
  trees <- function(sp, ids) data.table::data.table(
    MeasureID = 1L, OrigPlotID1 = ids, MeasureYear = 2000L, source = "BC", TreeNumber = seq_along(ids),
    DBH = 20, Height = 15, Species = sp, status = "A", first_tree_year = 2000L, last_tree_year = 2010L,
    diff_dbh = 0)
  ## 40 trees per plot (the forest filter needs >= 30 at the first measurement)
  ids <- rep("sbs", 40)
  plot <- data.table::data.table(MeasureID = 1L, OrigPlotID1 = "sbs", MeasureYear = 2000L, source = "BC",
                                 baseYear = 2000L, baseSA = 50, PlotSize = 0.1)
  gis <- sf::st_as_sf(data.table::data.table(OrigPlotID1 = "sbs", lon = -122, lat = 54),
                      coords = c("lon", "lat"), crs = 4326)
  equivLong <- data.table::data.table(
    Latin_full = c("Picea glauca", hybridLatin), SpBiomassEq = c("white spruce", "spruce"))
  prep <- function(mergeHybridSprucePSP) suppressMessages(prepPSPaNPP(
    studyAreaANPP = NULL, PSPgis = gis, PSPmeasure = trees("Picea glauca", ids), PSPplot = plot,
    sppEquivLong = equivLong, useHeight = TRUE, biomassModel = "Lambert2005",
    PSPperiod = c(1990, 2020), minDBH = 0, sppEquiv = mergedEquiv(),
    mergeHybridSprucePSP = mergeHybridSprucePSP, excludeBECzones = "BWBS", BECzones = synthBEC()))
  off <- prep(NA_character_)
  on <- prep("engelmann")
  expect_identical(unique(off$Latin_full), "Picea glauca")
  expect_identical(unique(on$Latin_full), hybridLatin)
  expect_equal(on$biomass, off$biomass)
  ## and that equation is the white spruce one, not the hybrid's generic spruce one
  asSpruce <- pemisc::biomassCalculation(species = "spruce", DBH = 20, height = 15, includeHeight = TRUE,
                                        equationSource = "Lambert2005")$biomass / 10
  expect_false(isTRUE(all.equal(unique(on$biomass), asSpruce)))
})
