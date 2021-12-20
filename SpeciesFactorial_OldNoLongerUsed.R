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
        reproducible.showSimilarDepth = 5)
if (!require("Require")) {install.packages("Require"); require("Require")}
Require("PredictiveEcology/SpaDES.install")
installSpaDES()
doExperiment <- TRUE

Require(c("PredictiveEcology/SpaDES.core@development", 
          "PredictiveEcology/LandR@minRelativeB (>= 1.0.6.9007)",
          "data.table", "raster", "viridis"), upgrade = FALSE)

source("~/Biomass_speciesParameters/R/factorialGenerators.R")
argsForFactorial <- list(cohortsPerPixel = 1:2,
             growthcurve = seq(0.65, 0.85, 0.02),
             mortalityshape = seq(20, 25, 1),
             longevity = seq(125, 300, 25),
             mANPPproportion = seq(3.5, 6, 0.25))

# Next sequence is all dependent on argsForFactorial, so do digest once
dig <- fastdigest::fastdigest(args)
species1And2 <- Cache(do.call, factorialSpeciesTable, argsForFactorial, 
                      .cacheExtra = dig, omitArgs = c("args"))
speciesEcoregion <- Cache(factorialSpeciesEcoregion, species1And2, .cacheExtra = dig, omitArgs = c("speciesTable"))
cohortData <- Cache(factorialCohortData, species1And2, speciesEcoregion, .cacheExtra = dig, 
                    omitArgs = c("speciesTable", "speciesEcoregion"))
speciesTable <- Cache(factorialSpeciesTableFillOut, species1And2, .cacheExtra = dig,
                      omitArgs = "speciesTable")

# Maps
pixelGroupMap <- factorialPixelGroupMap(cohortData)
studyArea <- as(extent(pixelGroupMap), 'SpatialPolygons')
crs(studyArea) <- crs(pixelGroupMap)
rasterToMatch <- pixelGroupMap
ecoregionMap <- pixelGroupMap
levels(ecoregionMap) <- data.frame(ID = 1:max(cohortData$pixelGroup, na.rm = TRUE), 
                                   ecoregion = 1, ecoregionGroup = 1, stringsAsFactors = TRUE)

# Simple Tables
minRelativeB <- data.table("ecoregionGroup" = factor(1), minRelativeBDefaults())
ecoregion <- data.table("ecoregionGroup" = as.factor(1), 'active' = 'yes')

#Make sppColors
sppColors <- viridis::viridis(n = nrow(speciesTable))
names(sppColors) <-  speciesTable$species

####RUN IT#####
times <- list(start = 0, end = 300)

parameters <- list(
  Biomass_core = list(.plotInitialTime = NA,
              .saveInitialTime = NA,
              .saveInterval = NA,
              .useParallel = 1,
              seedingAlgorithm = "noSeeding",
              calcSummaryBGM = NULL,
              .plots = NULL,
              .maxMemory = 1e9,
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


# Get modules
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
  "LandR.assertions" = FALSE,
  "LandR.verbose" = 1,
  "reproducible.overwrite" = TRUE,
  "reproducible.useMemoise" = TRUE, # Brings cached stuff to memory during the second run
  "reproducible.useCache" = TRUE,
  "spades.moduleCodeChecks" = FALSE, # Turn off all module's code checking
  "spades.useRequire" = TRUE
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

#profvis::profvis(

mySimOut <- Cache(simInitAndSpades, times = times, params = parameters, modules = modules, 
                  objects = objects, outputs = outputs, debug = 1, .cacheExtra = dig, 
                  omitArgs = c("objects"))
# )

# mySimOut <- spades(mySim, debug = 1)

#####Pull in the files#####

cdFiles <- as.data.table(outputs(mySimOut))[saved == TRUE]
fEs <- .fileExtensions()
cds <- by(cdFiles, cdFiles[, "saveTime"], function(x) {
  fE <- reproducible:::fileExt(x$file)
  wh <- fEs[fEs$exts %in% fE,]
  getFromNamespace(wh$fun, ns = asNamespace(wh$package))(x$file)
})  
cds <- rbindlist(cds, use.names = TRUE)


#This is dumb, the merge shouldn't be necessary on rerun
# temp <- species1[, .(species, maxANPPpct)]
# cds <- cds[temp, on = c("speciesCode" = 'species')]
# saveRDS(cds, file = "factorialCohortData.Rdat")
reducedFactorial <- cds[, .(speciesCode, age, B)]
# set(cds, NULL, c("ecoregionGroup", "mortalitly", "aANPPAct"), NULL)
factorialSpeciesTable <- speciesTable[,.(species, longevity, growthcurve, mortalityshape, mANPPproportion)]

Cache(saveRDS, reducedFactorial, file = "data/reducedFactorialCD.rds", quick = "file")
Cache(saveRDS, factorialSpeciesTable, file = "data/factorialSpeciesTable.rds", quick = "file")


uniq <- unique(cds$pixelGroup)
cds[, Sp := gsub(".+_", "", speciesCode)]; 
cds[, maxB := max(B), by = "pixelGroup"]
setkeyv(cds, c("pixelGroup", "Sp", "age"))
cds[, diffB := diff(B), by = c("age", "pixelGroup")]
cds[is.na(diffB), diffB := 0]
cds[, maxDiffB := max(diffB, na.rm = TRUE), by = c("pixelGroup")]
setorderv(cds, "maxDiffB")

# maxDiffBs <- head(unique(cds$maxDiffB), 40)

sam <- Cache(sample, uniq, 64)
ff <- cds[pixelGroup %in% sam]; 
# ff <- cds[maxDiffB %in% maxDiffBs]; 
ff <- factorialSpeciesTable[ff, on = c("species" = "speciesCode")]

# ff <- ff[pixelGroup %in% sample(unique(ff$pixelGroup), 49)]; 
ff[, Title := paste0(maxDiffB, "_", pixelGroup)]
ff[, 
   params := paste0(unique(Sp),"(l=",unique(longevity),";g=",unique(growthcurve), ";m=",unique(mortalityshape),";p=", unique(mANPPproportion ),")"), 
   by = c("Sp", "pixelGroup")]
ff[, Title := paste0(unique(params), collapse = "\n"), by = "pixelGroup"]
title <- paste0("Factorial Experiment: ", length(sam), "random")
gg1 <- ggplot(ff, aes(x = age, y = B, colour = Sp)) + 
  geom_line() + 
  facet_wrap(~ Title, nrow = ceiling(sqrt(length(sam))), scales = "fixed") +
  ggtitle(title) +
  theme(strip.text.x = element_text(size = 5)) #+


filename <- paste0("factorial2Sp_", Sys.time())
ggsave(paste0(filename, ".pdf"), gg1, width = 12, height = 7)
