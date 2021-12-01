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

# library(crayon)
#set out the combinations for factorial
if (FALSE) {
  growthcurves <- seq(0, 1, 0.1)
  mortalityshapes <- seq(5, 25, 1)
  longevity <- seq(150, 700, 25)
  mANPPproportion <- seq(0.25,10, 0.25)
}
growthcurves <- seq(0.65, 0.85, 0.02)
mortalityshapes <- seq(18, 25, 1)
longevity <- seq(125, 300, 25)
mANPPproportion <- seq(3.5, 6, 0.25)

shadetolerance <- c(1) #turns out shade tolerance is unimportant for this particular scenario
#as it only controls probability of establishment, not ANPP (previously factorial had understories)

Attributes <- c('longevity', 'growthcurve', 'mortalityshape', 'mANPPproportion')
species1 <- expand.grid(longevity, growthcurves, mortalityshapes, mANPPproportion)
names(species1) <- Attributes
species1$species <- paste0("A", 1:nrow(species1))
species1 <- data.table(species1)

cohortData <- data.table('speciesCode' = species1$species)
cohortData$pixelGroup <- 1:nrow(cohortData)

species1 <- speciesTable2CohortGenerate()
# speciesTable2CohortGenerate <- function(Attributes = list(growthcurve = seq(0.65, 0.85, 0.02),
#                                                           mortalityshape = seq(18, 25, 1),
#                                                           longevity = seq(125, 300, 25),
#                                                           mANPPproportion = seq(3.5, 6, 0.25)) ) {
#   Atts <- names(Attributes)
#   Attributes1 <- paste0(Atts, "1")
#   Attributes2 <- paste0(Atts, "2")
#   species2 <- expand.grid(growthcurves, growthcurves, 
#                           longevity, longevity, 
#                           mANPPproportion, mANPPproportion,
#                           mortalityshapes, mortalityshapes)
#   species2 <- as.data.table(species2)
#   colnames2sp <- sort(c(Attributes1, Attributes2))
#   setnames(species2, new = colnames2sp)
#   species2a <- species2[species2$longevity1 > species2$longevity2 &
#                           species2$growthcurve1 > species2$growthcurve2 &
#                           species2$mortalityshape1 > species2$mortalityshape2 &
#                           species2$mANPPproportion1 > species2$mANPPproportion2]
#   
#   set(species2a, NULL, "pixelGroup", seq(NROW(species2a)))
#   set(species2a, NULL, "species", paste0("A", species2a$pixelGroup))
#   
#   species2b <- melt(species2a, id.vars = c("species", "pixelGroup")) 
#   set(species2b, NULL, "Sp", gsub(".+(1$|2$)", "Sp\\1", species2b$variable))
#   set(species2b, NULL, "variable", gsub("1$|2$", "", species2b$variable))
#   species2b <- dcast(species2b, pixelGroup + species + Sp ~ variable )
#   set(species2b, NULL, "speciesCode", paste0(species2b$species, "_", species2b$Sp))
#   set(species2b, NULL, c("species", "Sp"), NULL)
#   cohortData2 <- copy(species2b[, c("speciesCode", "pixelGroup")])
#   set(cohortData2, NULL, "pixelGroup", cohortData2$pixelGroup + max(cohortData$pixelGroup))
#   cohortData <- rbindlist(list(cohortData, cohortData2), use.names = TRUE)
#   setnames(species2b, old = "speciesCode", new = "species")
#   colsToKeep <- c(Atts, "species")
#   species2 <- species2b[, ..colsToKeep]
#   species1 <- rbindlist(list(species1, species2), use.names = TRUE)
# }


#age, B, totalB, speciesProportion - added later

#Get number of pixel groups

#Now set the extra species column as null, set up the biomass and age
set(species1, NULL, "mortalityshape", asInteger(species1$mortalityshape))
set(species1, NULL, "longevity", asInteger(species1$longevity))
set(cohortData, NULL, "age", 1)
set(species1, NULL, "maxB", 5000)
set(species1, NULL, "maxANPP", asInteger(species1$maxB * species1$mANPPproportion/100))
set(cohortData, NULL, "B", species1$maxANPP)
# cohortData$B <- species1$maxANPP

# cohortData$speciesProportion <- 100
# cohortData$sumB <- cohortData$B


#####Make LANDR Inputs####
set(cohortData, NULL, "ecoregionGroup", factor(1))
# cohortData$ecoregionGroup <- 1
cohortData <- setcolorder(cohortData, c('speciesCode', 'pixelGroup', 'ecoregionGroup', 'age', "B"))#, 'sumB', 'speciesProportion'))

speciesEcoregion <- copy(species1[, c("maxB", "maxANPP", "species")])
speciesEcoregion[, c("ecoregionGroup", "establishprob", "maxB", "maxANPP", "year") := .(1, 0.5, species1$maxB, species1$maxANPP, 0)]
setnames(speciesEcoregion, old = "species", new = "speciesCode")

# Change classes
set(speciesEcoregion, NULL, c("maxB", "maxANPP", "speciesCode"), 
    list(
      asInteger(speciesEcoregion$maxB), 
      as.numeric(speciesEcoregion$maxANPP),
      as.factor(speciesEcoregion$speciesCode)
      ))
# speciesEcoregion[, c("mANPPproportion", "growthcurve", "mortalityshape", "longevity", "species") := NULL]
# speciesEcoregion[, "speciesCode" := species1$species]

SpeciesTable <- copy(species1[, c("species", "longevity", "growthcurve", "mortalityshape")])
SpeciesTable[, c("sexualmature", 'seeddistance_eff', 'seeddistance_max', 'resproutprob', 
                 'resproutage_min', 'resproutage_max', 'postfireregen',
                 'leaflongevity', 'wooddecayrate', 'leafLignin', 'hardsoft', "Area", "firetolerance") :=
               list(30L, 0L, 0L, 0.5, 0L, 0L, factor('none'), 3L, 0.07, 0.1, factor('soft'), factor("BP"), 3L)]
set(SpeciesTable, NULL, "shadetolerance", asInteger(shadetolerance))


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
sppColors <- viridis::viridis(n = nrow(SpeciesTable))
names(sppColors) <-  SpeciesTable$species

####RUN IT#####

rasterOptions(tmpdir = "temp")
times <- list(start = 0, end = 700)
#Do this so one cohort is alive at time == 700. This cohort is unlikely to ever match with anything real,
# SpeciesTable[longevity == 700 & growthcurve == 0 & mortalityshape == 25 & maxANPP == 250, longevity := 701]

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
