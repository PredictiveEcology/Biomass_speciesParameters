library(data.table)
library(raster)
#Set up the overstory 'A' species
growthcurves <- c(0, 0.25, 0.5, 0.75, 1)
MortCurves <- c(5, 10, 15, 20, 25)
longevity <- c(150, 200, 250, 300)
shadetolerance <- c(1,3,5)
Attributes <- c('longevity', 'growthcurve', 'mortalityshape', "shadetolerance")

#Pure stands 1:175 There will be double this amount of pure stands because I will run for different ages
species1 <- expand.grid(longevity, growthcurves, MortCurves, shadetolerance)
names(species1) <- Attributes
species1$species <- paste0("A", 1:nrow(species1))
species1 <- data.table(species1)
#age, B, totalB, speciesProportion - added later

#Make mixed stands
species2 <- expand.grid(longevity, growthcurves, MortCurves, shadetolerance)
names(species2) <- c("longevity", "growthcurve","mortalityshape", "shadetolerance")
species2$species <- paste0("B", 1:nrow(species1))
species2 <- data.table(species2)


#Get number of pixel groups
cohortData <- data.table(expand.grid(species1$species, species2$species))
names(cohortData) <- c("Species1", 'Species2')
cohortData$pixelGroup <- 1:nrow(cohortData)
setkey(cohortData, Species1, Species2)

firstSpecies <- species1[cohortData, on = c(species = "Species1")]
secondSpecies <- species2[cohortData, on = c(species = "Species2")]

#Now set the extra species column as null, set up the biomass and age
firstSpecies[, Species2 := NULL]
secondSpecies[, Species1 := NULL]
firstSpecies$B <- 500
firstSpecies$age <- 50
secondSpecies$B <- 50
secondSpecies$age <- 10
firstSpecies$speciesProportion <- 500/550
secondSpecies$speciesProportion <- 50/550
#Species B will always be the understory, and has identical traits to species A

cohortData <- rbind(firstSpecies, secondSpecies)
cohortData$totalB <- 550
#Adjust pixelGroups to account for 1100 more species
cohortData$pixelGroup <- cohortData$pixelGroup + nrow(species1)*2


#Finally, add back in the pure stands (PureStands)
species1$pixelGroup <- 1:nrow(species1)
species2$pixelGroup <- species1$pixelGroup + nrow(species1)
species2$B <- 50
species2$age <- 10
species1$B <- 500
species1$age <- 50

PureStands <- rbind(species1, species2)

PureStands$speciesProportion <- 1
PureStands$totalB <- PureStands$B

#####Make LANDR Inputs####

cohortData <- rbind(PureStands, cohortData)
cohortData[, speciesCode := species] #
cohortData$ecoregionGroup <- 1
cohortData <- setcolorder(cohortData, c('species', 'pixelGroup', 'ecoregionGroup', 'age', "B", 'totalB', 'speciesProportion'))

speciesEcoregion <- copy(PureStands)
speciesEcoregion[, c("ecoregionGroup", "establishprob", "maxB", "maxANPP", "year") := .(1, 0.5, 5000, 5000/30, 0)]
speciesEcoregion[, c("speciesProportion", "totalB", "growthcurve", "mortalityshape", "age", "B", "longevity", "pixelGroup") := NULL]
speciesEcoregion[, "speciesCode" := species]

SpeciesTable <- copy(PureStands)
SpeciesTable[, c('totalB', 'speciesProportion', 'age', 'B', 'pixelGroup') := NULL]
SpeciesTable[, c("sexualmature", 'SeedEffDist', 'SeedMaxDist', 'VegProb', 'MinAgeVeg', 'MaxAgeVeg', 'PostFireRegen',
                 'Leaflongevity', 'WoodDecayRate', 'LeafLignin', 'HardSoft') :=
               list(30, 0, 0, 0.5, 0, 0, 'none', 3, 0.07, 0.1, 'soft'), 'species']
SpeciesTable[, shade := shadetolerance]


pixelGroupMap <- raster(res = c(1,1))
nrow(pixelGroupMap) <- round(sqrt(max(cohortData$pixelGroup)), 0)
ncol(pixelGroupMap) <- round(sqrt(max(cohortData$pixelGroup)), 0)
vals <- c(1:max(cohortData$pixelGroup), rep(NA, ncell(pixelGroupMap) - max(cohortData$pixelGroup)))
pixelGroupMap[] <- vals

minRelativeB <- data.table("ecoregionGroup" = 1, X1 = 0.2, X2 = 0.4, X3 = 0.5, X4 = 0.7, X5 = 0.9)
ecoregion <- data.table("ecoregionGroup" = 1, 'active' = 'yes')

##Make ecoregion Map
ecoregionMap <- pixelGroupMap
levels(ecoregionMap) <- data.frame(ID = 1:max(cohortData$pixelGroup, na.rm = TRUE), ecoregion = 1, ecoregionGroup = 1, stringsAsFactors = TRUE)

#Make sppColors
sppColors <- viridis::viridis(n = nrow(SpeciesTable))
names(sppColors) <-  SpeciesTable$species

rm(species1, species2, PureStands, firstSpecies, secondSpecies)



####RUN IT#####
library(SpaDES)
library(raster)
library(LandR)
library(sp)

rasterOptions(tmpdir = "temp/")
# library(sf)
spadesModulesDirectory <-  c(file.path("modules")) # where modules are 
modules <- list("LBMR")
times <- list(start = 0, end = 1)
# "PSP_Clean", "gmcsDataPrep",


studyArea <- as(extent(pixelGroupMap), 'SpatialPolygons') 
crs(studyArea) <- crs(pixelGroupMap)
rasterToMatch <- pixelGroupMap


parameters <- list(
  LBMR = list(.plotInitialTime = NA,
              .saveInitialTime = 0,
              .saveInterval = 20,
              seedingAlgorithm = "noDispersal",
              useCache = TRUE,
              successionTimestep = 10,
              initialBiomassSource = "cohortData",
              vegLeadingProportion = 0,
              growthAndMortalityDrivers = "LandR")
)

## Paths are not workign with multiple module paths yet
setPaths(cachePath =  file.path("temp/Cache"),
         modulePath = file.path(getwd(), "modules"),
         inputPath = file.path(getwd(), "inputs"),
         outputPath = file.path(getwd(),"outputs/sensitivityCohortData/"))
paths <- SpaDES.core::getPaths()

options("spades.moduleCodeChecks" = FALSE)

#Tree species that are important to us
speciesLayers <- "species"

speciesEcoregion[shadetolerance == 1, maxB := 3000,]
speciesEcoregion[shadetolerance == 5, maxB := 7000,]
speciesEcoregion$maxANPP <- speciesEcoregion$maxB/30

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

