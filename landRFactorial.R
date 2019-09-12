library(data.table)
library(raster)
#Set up the overstory 'A' species
growthcurves <- seq(0, 1, 0.1)
MortCurves <- c(5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 20, 22, 25)
longevity <- seq(150, 700, 25)
mANPPproportion <- seq(1,5, 0.25)


shadetolerance <- c(1)
Attributes <- c('longevity', 'growthcurve', 'mortalityshape', "shadetolerance", 'mANPPproportion')

#Pure stands 1:175 There will be double this amount of pure stands because I will run for different ages
species1 <- expand.grid(longevity, growthcurves, MortCurves, shadetolerance, mANPPproportion)
names(species1) <- Attributes
species1$species <- paste0("A", 1:nrow(species1))
species1 <- data.table(species1)
#age, B, totalB, speciesProportion - added later

#Get number of pixel groups
cohortData <- data.table('species' = species1$species)
cohortData$pixelGroup <- 1:nrow(cohortData)

#Now set the extra species column as null, set up the biomass and age
cohortData$age <- 1
species1$maxB <- 5000
species1$maxANPP <- species1$maxB * species1$mANPPproportion/100
cohortData$B <- species1$maxANPP 

cohortData$speciesProportion <- 100
cohortData$totalB <- cohortData$B #Fix 
#Species B will always be the understory, and has identical traits to species A


#####Make LANDR Inputs####
cohortData$ecoregionGroup <- 1
cohortData <- setcolorder(cohortData, c('species', 'pixelGroup', 'ecoregionGroup', 'age', "B", 'totalB', 'speciesProportion'))

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
library(SpaDES)
library(raster)
library(LandR)
library(sp)

rasterOptions(tmpdir = "temp")
spadesModulesDirectory <-  c(file.path("../Land-R/modules/")) # where modules are 
modules <- list("LBMR")
times <- list(start = 0, end = 701)
#Do this so one cohort is alive at time == 700. This cohort is unlikely to ever match with anything real, 
SpeciesTable[longevity == 700 & growthcurve == 0 & mortalityshape == 25 & maxANPP == 250, longevity := 701]

# "PSP_Clean", "gmcsDataPrep",


studyArea <- as(extent(pixelGroupMap), 'SpatialPolygons') 
crs(studyArea) <- crs(pixelGroupMap)
rasterToMatch <- pixelGroupMap


parameters <- list(
  LBMR = list(.plotInitialTime = NA,
              .saveInitialTime = 0,
              .saveInterval = 10,
              seedingAlgorithm = "noDispersal",
              useCache = TRUE,
              successionTimestep = 10,
              initialBiomassSource = "cohortData",
              vegLeadingProportion = 0,
              growthAndMortalityDrivers = "LandR")
)

## Paths are not workign with multiple module paths yet
setPaths(cachePath =  file.path("temp/Cache"),
         modulePath = file.path("../Land-R/modules/"),
         inputPath = file.path(getwd(), "inputs"),
         outputPath = file.path(getwd(),"outputs/sensitivityCohortData/"))
paths <- SpaDES.core::getPaths()

options("spades.moduleCodeChecks" = FALSE)

#Tree species that are important to us
speciesLayers <- "species"

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
  "reproducible.cachePath" = paths$cachePath,
  "reproducible.useCloud" = FALSE,
  "spades.moduleCodeChecks" = FALSE, # Turn off all module's code checking
  "spades.useRequire" = TRUE, # assuming all pkgs installed correctly
  "pemisc.useParallel" = FALSE
)

set.seed(161616)

mySim <- simInit(times = times, params = parameters, modules = modules, objects = objects,
                 paths = paths, loadOrder = unlist(modules))

mySimOut <- spades(mySim, debug = 2)

