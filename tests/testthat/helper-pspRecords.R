## PSP records for the buildGrowthCurves() tests, in the LandR species codes of `latinOf`
latinOf <- c(Abie_las = "Abies lasiocarpa", Pice_eng = "Picea engelmannii",
             Pice_eng_gla = "Picea engelmannii x glauca", Pice_mar = "Picea mariana",
             Popu_tre = "Populus tremuloides")

## sppEquiv for the LandR codes `spp`
equivFor <- function(spp) data.table::data.table(Latin_full = unname(latinOf[spp]), LandR = spp)

## PSP records as prepPSPaNPP() hands them to buildGrowthCurves(), with the columns it reads or drops:
## here one record per species x plot x measurement year, with that species' biomass on a 1-ha plot
pspRecords <- function(species, plot, year, age, biomass) {
  data.table::data.table(
    MeasureID = paste(plot, year, sep = "_"), OrigPlotID1 = plot, MeasureYear = year,
    Latin_full = unname(latinOf[species]), SpBiomassEq = "x", source = "x", PSP = "x",
    Elevation = 0, PlotSize = 1, standAge = age, biomass = biomass)
}
