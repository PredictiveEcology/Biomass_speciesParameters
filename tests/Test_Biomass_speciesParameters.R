library(SpaDES)
library(magrittr)

modulePath <- "../"

inputDir <- file.path(dirname(tempdir()), "Biomass_speciesParameters", "inputs") %>% checkPath(create = TRUE)
outputDir <- file.path(dirname(tempdir()), "Biomass_speciesParameters", "outputs")
times <- list(start = 0, end = 10)
parameters <- list()
modules <- list("Biomass_speciesParameters")
paths <- list(
  cachePath = file.path(outputDir, "cache"),
  modulePath = modulePath,
  inputPath = inputDir,
  outputPath = outputDir
)


mySim <- simInit(times = times, params = parameters, modules = modules, paths = paths)

mySimOut <- spades(mySim, debug = TRUE)
