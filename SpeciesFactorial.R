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
        reproducible.useMemoise = TRUE, reproducible.cacheSaveFormat = "qs",
        spades.moduleCodeChecks = FALSE, LandR.assertions = FALSE)
if (!require("Require")) {install.packages("Require"); require("Require")}
Require("PredictiveEcology/SpaDES.install")
installSpaDES()
doExperiment <- TRUE

Require(c("PredictiveEcology/SpaDES.core@development", 
          "PredictiveEcology/LandR@development (>= 1.0.6.9004)",
          "data.table", "raster", "viridis"), upgrade = FALSE)

# library(crayon)
#set out the combinations for factorial
growthcurves <- seq(0, 1, 0.1)
mortalityshapes <- seq(5, 25, 1)
longevity <- seq(150, 700, 25)
mANPPproportion <- seq(0.25,10, 0.25)

shadetolerance <- c(1) #turns out shade tolerance is unimportant for this particular scenario
#as it only controls probability of establishment, not ANPP (previously factorial had understories)

Attributes <- c('longevity', 'growthcurve', 'mortalityshape', "shadetolerance", 'mANPPproportion')

species1 <- expand.grid(longevity, growthcurves, mortalityshapes, shadetolerance, mANPPproportion)
names(species1) <- Attributes
species1$species <- paste0("A", 1:nrow(species1))
species1 <- data.table(species1)
#age, B, totalB, speciesProportion - added later

#Get number of pixel groups
cohortData <- data.table('speciesCode' = species1$species)
cohortData$pixelGroup <- 1:nrow(cohortData)

#Now set the extra species column as null, set up the biomass and age
cohortData$age <- 1
species1$maxB <- 5000
species1$maxANPP <- asInteger(species1$maxB * species1$mANPPproportion/100)
cohortData$B <- species1$maxANPP

cohortData$speciesProportion <- 100
cohortData$sumB <- cohortData$B


#####Make LANDR Inputs####
cohortData$ecoregionGroup <- 1
cohortData <- setcolorder(cohortData, c('speciesCode', 'pixelGroup', 'ecoregionGroup', 'age', "B", 'sumB', 'speciesProportion'))

speciesEcoregion <- copy(species1)
speciesEcoregion[, c("ecoregionGroup", "establishprob", "maxB", "maxANPP", "year") := .(1, 0.5, species1$maxB, species1$maxANPP, 0)]
speciesEcoregion[, c("mANPPproportion", "growthcurve", "mortalityshape", "longevity", "species") := NULL]
speciesEcoregion[, "speciesCode" := species1$species]

SpeciesTable <- copy(species1)
SpeciesTable[, c("sexualmature", 'SeedEffDist', 'SeedMaxDist', 'VegProb', 'MinAgeVeg', 'MaxAgeVeg', 'PostFireRegen',
                 'Leaflongevity', 'WoodDecayRate', 'LeafLignin', 'HardSoft') :=
               list(30, 0, 0, 0.5, 0, 0, 'none', 3, 0.07, 0.1, 'soft'), 'species']
SpeciesTable[, shade := shadetolerance]


pixelGroupMap <- raster(res = c(1,1))
nrow(pixelGroupMap) <- round(sqrt(max(cohortData$pixelGroup)), 0)
ncol(pixelGroupMap) <- round(sqrt(max(cohortData$pixelGroup)), 0) + 1
vals <- c(1:max(cohortData$pixelGroup), rep(NA, times = ncell(pixelGroupMap) - max(cohortData$pixelGroup)))
pixelGroupMap[] <- vals

minRelativeB <- data.table("ecoregionGroup" = 1, X1 = 0.2, X2 = 0.4, X3 = 0.5, X4 = 0.7, X5 = 0.9)
ecoregion <- data.table("ecoregionGroup" = 1, 'active' = 'yes')

##Make ecoregion Map
ecoregionMap <- pixelGroupMap
levels(ecoregionMap) <- data.frame(ID = 1:max(cohortData$pixelGroup, na.rm = TRUE), ecoregion = 1, ecoregionGroup = 1, stringsAsFactors = TRUE)

#Make sppColors
sppColors <- viridis::viridis(n = nrow(SpeciesTable))
names(sppColors) <-  SpeciesTable$species

####RUN IT#####

rasterOptions(tmpdir = "temp")
times <- list(start = 0, end = 700)
#Do this so one cohort is alive at time == 700. This cohort is unlikely to ever match with anything real,
SpeciesTable[longevity == 700 & growthcurve == 0 & mortalityshape == 25 & maxANPP == 250, longevity := 701]

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

options("spades.moduleCodeChecks" = FALSE)

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
  species = SpeciesTable,
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
  "LandR.assertions" = FALSE,
  "LandR.verbose" = 1,
  "reproducible.futurePlan" = FALSE,
  "reproducible.inputPaths" = NULL,
  "reproducible.quick" = FALSE,
  "reproducible.overwrite" = TRUE,
  "reproducible.useMemoise" = TRUE, # Brings cached stuff to memory during the second run
  "reproducible.useNewDigestAlgorithm" = TRUE,  # use the new less strict hashing algo
  "reproducible.useCache" = TRUE,
  # "reproducible.cachePath" = paths$cachePath,
  "reproducible.useCloud" = FALSE,
  "spades.moduleCodeChecks" = FALSE, # Turn off all module's code checking
  "spades.useRequire" = TRUE, # assuming all pkgs installed correctly
  "pemisc.useParallel" = FALSE
)

set.seed(161616)

#EDIT ALGO 2 IN Biomass_core/HELPERS TO ALGO 1. #Also made succesionTimeStep 1 so calculateSumB wouldnt' return NA
#Also had to completely remove lines in NoDiserpsal, to shut off all regeneration.
#removed all LandR.CS references Nov 28th
####NOTE: This will fail at year 700, because every cohort is dead. Not sure why that fails yet...
#files are still output so it isn't a problem
mySim <- simInit(times = times, params = parameters, modules = modules, 
                 objects = objects
                 #paths = paths, loadOrder = unlist(modules)
                 )

mySimOut <- spades(mySim, debug = 1)

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
