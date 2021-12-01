#This script is for running the species traits factorial design. 
# it requires a version of Biomass_core that is gitIgnored. The versino had many changes to accomodate the 'no regeneration' scenario
# Many of the other changes have been subsequently incorporated into Biomass_core, so it may work with a newer version
rootDir <- file.path("~/Yield")
if (isTRUE(identical(Sys.info()[["user"]], "elmci1"))) {
  options(repos = c(CRAN = "https://cloud.r-project.org/"),
          "Require.install" = FALSE)
  rootDir <- "~/projects/def-elmci1-ab/elmci1/Yield"
}
options(reproducible.showSimilar = TRUE, reproducible.inputPaths = file.path(dirname(rootDir),"/Eliot/data"),
        reproducible.useMemoise = TRUE, reproducible.cacheSaveFormat = "qs")
if (!require("Require")) {install.packages("Require"); require("Require")}
Require("PredictiveEcology/SpaDES.install")
installSpaDES()
doExperiment <- TRUE

Require(c("PredictiveEcology/SpaDES.core@development", 
          "PredictiveEcology/LandR@development (>= 1.0.6.9005)",
          "data.table", "raster", "viridis"), upgrade = FALSE)

source("~/Biomass_speciesParameters/R/factorialGenerators.R")
species1And2 <- Cache(factorialSpeciesTable, cohortsPerPixel = 1:2,
                                     growthcurve = seq(0.65, 0.85, 0.02),
                                     mortalityshape = seq(18, 25, 1),
                                     longevity = seq(125, 300, 25),
                                     mANPPproportion = seq(3.5, 6, 0.25))
speciesEcoregion <- factorialSpeciesEcoregion(species1And2)
cohortData <- factorialCohortData(species1And2, speciesEcoregion)
speciesTable <- factorialSpeciesTableFillOut(species1And2)



pixelGroupMap <- raster(res = c(1,1))
nrow(pixelGroupMap) <- round(sqrt(max(cohortData$pixelGroup)), 0)
ncol(pixelGroupMap) <- round(sqrt(max(cohortData$pixelGroup)), 0) + 1
vals <- c(1:max(cohortData$pixelGroup), rep(NA, times = ncell(pixelGroupMap) - max(cohortData$pixelGroup)))
pixelGroupMap[] <- vals

minRelativeB <- data.table("ecoregionGroup" = 1, X1 = 0.2, X2 = 0.4, X3 = 0.5, X4 = 0.7, X5 = 0.9)
ecoregion <- data.table("ecoregionGroup" = as.factor(1), 'active' = 'yes')

##Make ecoregion Map
ecoregionMap <- pixelGroupMap
levels(ecoregionMap) <- data.frame(ID = 1:max(cohortData$pixelGroup, na.rm = TRUE), ecoregion = 1, ecoregionGroup = 1, stringsAsFactors = TRUE)

#Make sppColors
sppColors <- viridis::viridis(n = nrow(speciesTable))
names(sppColors) <-  speciesTable$species

####RUN IT#####

rasterOptions(tmpdir = "temp")
times <- list(start = 0, end = 700)
#Do this so one cohort is alive at time == 700. This cohort is unlikely to ever match with anything real,
# speciesTable[longevity == 700 & growthcurve == 0 & mortalityshape == 25 & maxANPP == 250, longevity := 701]

studyArea <- as(extent(pixelGroupMap), 'SpatialPolygons')
crs(studyArea) <- crs(pixelGroupMap)
rasterToMatch <- pixelGroupMap


parameters <- list(
  Biomass_core = list(.plotInitialTime = NA,
              .saveInitialTime = NA,
              .saveInterval = NA,
              .useParallel = 1,
              seedingAlgorithm = "noSeeding",
              calcSummaryBGM = NULL,
              .plots = NULL,
              useCache = TRUE,
              successionTimestep = 10,
              initialBiomassSource = "cohortData",
              vegLeadingProportion = 0
  ))

## Paths are not workign with multiple module paths yet
setPaths(rasterPath = "temp",
         cachePath =  file.path("temp/Cache"),
         modulePath = file.path("modules"),
         inputPath = file.path(getwd(), "inputs"),
         outputPath = file.path(getwd(),"outputs"))


moduleNameAndBranch <- c("Biomass_core@development") #, "Biomass_speciesParameters@EliotTweaks")
lapply(moduleNameAndBranch, function(modName) {
  Cache(getModule, file.path("PredictiveEcology", modName), #modulePath = getPaths()$modulePath, 
        overwrite = TRUE)
})
modules <- gsub("@.+", "", moduleNameAndBranch)


#Tree species that are important to us
speciesLayers <- "species"

#sppEquiv needed or module stops, but object unused, likewise with speciesLayers
objects <- list(
  "studyArea" = studyArea,
  "rasterToMatch" = rasterToMatch,
  cohortData = cohortData,
  species = speciesTable,
  speciesEcoregion = speciesEcoregion,
  pixelGroupMap = pixelGroupMap,
  speciesLayers = speciesLayers,
  minRelativeB = minRelativeB,
  ecoregion = ecoregion,
  ecoregionMap = ecoregionMap,
  sppEquiv = data.table(),
  sppColorVect = sppColors
)


opts <- options(
  "future.globals.maxSize" = 1000*1024^2,
  "LandR.assertions" = TRUE,
  "LandR.verbose" = 1,
  "reproducible.futurePlan" = FALSE,
  "reproducible.overwrite" = TRUE,
  "reproducible.useMemoise" = TRUE, # Brings cached stuff to memory during the second run
  "reproducible.useNewDigestAlgorithm" = TRUE,  # use the new less strict hashing algo
  "reproducible.useCache" = TRUE,
  "spades.moduleCodeChecks" = FALSE, # Turn off all module's code checking
  "spades.useRequire" = TRUE, # assuming all pkgs installed correctly
  "pemisc.useParallel" = FALSE
)

outputs <- data.frame(expand.grid(objectName = "cohortData",
                                  saveTime = unique(seq(times$start, times$end, by = 10)),
                                  eventPriority = 1, fun = "qs::qsave",
                                  stringsAsFactors = FALSE))


set.seed(161616)

#EDIT ALGO 2 IN Biomass_core/HELPERS TO ALGO 1. #Also made succesionTimeStep 1 so calculateSumB wouldnt' return NA
#Also had to completely remove lines in NoDiserpsal, to shut off all regeneration.
#removed all LandR.CS references Nov 28th
####NOTE: This will fail at year 700, because every cohort is dead. Not sure why that fails yet...
#files are still output so it isn't a problem

mySimOut <- simInitAndSpades(times = times, params = parameters, modules = modules, 
                 objects = objects, outputs = outputs, debug = 1)

# mySimOut <- spades(mySim, debug = 1)

#####Pull in the files#####

cdFiles <- list.files("outputs/", full.names = TRUE) 
cds <- lapply(cdFiles, FUN = 'readRDS') 
cds <- rbindlist(cds)


#This is dumb, the merge shouldn't be necessary on rerun
# temp <- species1[, .(species, maxANPPpct)]
# cds <- cds[temp, on = c("speciesCode" = 'species')]
# saveRDS(cds, file = "factorialCohortData.Rdat")
reducedFactorial <- cds[, .(speciesCode, age, B)]
factorialSpeciesTable <- species1[,.(species, longevity, growthcurve, mortalityshape, mANPPproportion)]

saveRDS(reducedFactorial, file = "data/reducedFactorialCD.Rdat")
saveRDS(factorialSpeciesTable, file = "data/factorialSpeciesTable.Rdat")
