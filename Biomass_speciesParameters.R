
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "Biomass_speciesParameters",
  description = "For estimating LANDIS-II species traits from PSP-derived growth curves",
  keywords = NA, # c("insert key words here"),
  authors = c(person(c("Ian"), "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = c("aut", "cre")),
              person(c("Eliot"), "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut")),
              person(c("Ceres"), "Barros", email = "ceres.barros@ubc.ca", role = c("ctb"))),
  childModules = character(0),
  version = list(Biomass_speciesParameters = "0.0.2"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_speciesParameters.Rmd"),
  reqdPkgs = list("mgcv", "nlme", "fpCompare", "crayon",
                  "PredictiveEcology/LandR@development (>= 1.0.5)",
                  "PredictiveEcology/pemisc@development",
                  "PredictiveEcology/reproducible@development (>= 1.2.10.9001)",
                  "PredictiveEcology/SpaDES.core@development (>= 0.2.6)"),
  parameters = rbind(
    defineParameter("biomassModel", "character", "Lambert2005", NA, NA,
                    desc =  paste("The model used to calculate biomass from DBH. Can be either 'Lambert2005' or 'Ung2008'")),
    defineParameter("constrainGrowthCurve", "numeric", c(0.5, 0.5), 0, 1,
                    desc = paste("Upper and lower bounds on range of potential growth curves when fitting traits. This module accepts a",
                                 "list of vectors, with names equal to `P(sim)$sppEquivCol`, so that traits are customizable")),
    defineParameter("constrainMortalityShape", 'numeric', c(15, 25), 5, 25,
                    desc = paste("Upper and lower bounds on mortality shape when fitting traits. Low mortality curve values lead to numerous",
                                 "cohorts with very little biomass as longevity is approached, adding computation strain.",
                                 "Alternatively accepts a list of vectors, with names equal to those in `sim$sppEquiv[, P(sim)$sppEquivCol]`")),
    defineParameter("constrainMaxANPP", "numeric", c(2.0, 5.0), 1, 10,
                    desc = paste("Upper and lower bounds on `maxANPP` when fitting traits. When cohorts are initiated with `B = maxANPP`",
                                 "(see *Biomass_core* parameter `P(sim)$initialB`), this can lead to unreasonably initial biomass",
                                 "if `maxANPP` is also high. Both `maxANPP` and `growthcurve` parameters control when `maxB`",
                                 "is reached. High `maxANPP` results in earlier peaks. Alternatively accepts a list of vectors,",
                                 "with names equal to those in `sim$sppEquiv[, P(sim)$sppEquivCol]`")),
    defineParameter("GAMMiterations", "numeric", 8, 1, NA,
                    desc = paste("Number of iterations for GAMMs. Accepts a list of vectors, with names equal to those in",
                                 "`sim$sppEquiv[, P(sim)$sppEquivCol]`, so that GAMMS are customizable per species")),
    defineParameter("GAMMknots", "numeric", 3, NA, NA,
                    desc = paste("The number of knots to use in the GAMM. Either 3 or 4 is recommended. Accepts a",
                                 "list of vectors, with names equal to those in `sim$sppEquiv[, P(sim)$sppEquivCol]`, so that GAMMS",
                                 "are customizable per species")),
    defineParameter("minDBH", "integer", 0L, 0L, NA,
                    desc = paste("Minimum diameter at breast height (DBH) in cm used to filter PSP data.",
                                 "Defaults to 0cm, i.e. all tree measurements are used.")),
    defineParameter("minimumPlotsPerGamm", "numeric", 50, 10, NA,
                    desc = paste("Minimum number of PSP plots before building GAMM")),
    defineParameter("PSPperiod", "numeric", c(1920, 2019), NA, NA,
                    desc = paste("The years by which to subset sample plot data, if desired. Must be a vector of length 2")),
    defineParameter("quantileAgeSubset", "numeric", 95, 1, 100,
                    desc = paste("Quantile by which to subset PSP data. As older stands are sparsely represented, the oldest measurements",
                                 "become vastly more influential. This parameter accepts both a single value and a list of vectors",
                                 "named by `sim$sppEquiv[, P(sim)$sppEquivCol]`. The PSP stand ages are found in",
                                 "`sim$speciesGAMMs$SPECIES$originalData`, where SPECIES is the species ID")),
    defineParameter("sppEquivCol", "character", "default", NA, NA,
                    paste("The column in `sim$sppEquiv` data.table to group species by. This parameter should share the same",
                          "name as in *Biomass_borealDataPrep*. PSPs are aggregated by names in the PSP column and traits estimated",
                          "for species with corresponding names in the `sim$sppEquiv[, P(sim)$sppEquivCol]`")),
    defineParameter("useHeight", "logical", FALSE, NA, NA,
                    desc = paste("Should height be used to calculate biomass (in addition to DBH)?",
                                 "We advise against including height unless you are certain it is present in every PSP")),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    desc = paste("Should this entire module be run with caching activated?",
                                 "This is generally intended for data-type modules, where stochasticity and time are not relevant"))
  ),
  inputObjects = bindrows(
    expectsInput(objectName  = "factorialSpeciesTable", objectClass = "data.table",
                 desc = paste("Table with species traits to be matched with `sim$reducedFactorialCohortData`"),
                 sourceURL = "https://drive.google.com/open?id=1q0ou0CBzD9GqGSparpHqf318IWK6ycty"),
    expectsInput(objectName = "reducedFactorialCohortData", objectClass = "data.table",
                 desc = paste("Results of factorial species trait simulation. This can be found by running",
                              "`SpeciesFactorial.R` but requires a specific commit of *Biomass_core*"),
                 sourceURL = "https://drive.google.com/open?id=1h8StXE0vm8xyDycRomCkwIaL7wfh5Irj"),
    expectsInput(objectName = "PSPmeasure", objectClass = "data.table",
                 desc = paste("Merged PSP and TSP individual tree measurements. Must include the following columns:",
                              "'MeasureID', 'OrigPlotID1', 'MeasureYear', 'TreeNumber', 'Species', 'DBH' and 'newSpeciesName'",
                              "the latter corresponding to species names in `LandR::sppEquivalencies_CA$PSP`.",
                              "Defaults to randomized PSP data stripped of real plotIDs"),
                 sourceURL = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing"),
    expectsInput(objectName = "PSPplot", objectClass = "data.table",
                 desc = paste("Merged PSP and TSP plot data. Defaults to randomized PSP data stripped of real plotIDs.",
                              "Must contain columns 'MeasureID', 'MeasureYear', 'OrigPlotID1', and 'baseSA',",
                              "the latter being stand age at year of first measurement"),
                 sourceURL = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing"),
    expectsInput(objectName = "PSPgis", objectClass = "sf",
                 desc = paste("Plot location `sf` object. Defaults to PSP data stripped of real plotIDs/location.",
                              "Must include field 'OrigPlotID1' for joining to PSPplot object"),
                 sourceURL = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing"),
    expectsInput(objectName = "species", objectClass = "data.table",
                 desc = paste("A table of invariant species traits with the following trait colums:",
                              "'species', 'Area', 'longevity', 'sexualmature', 'shadetolerance',",
                              "'firetolerance', 'seeddistance_eff', 'seeddistance_max', 'resproutprob',",
                              "'mortalityshape', 'growthcurve', 'resproutage_min', 'resproutage_max',",
                              "'postfireregen', 'wooddecayrate', 'leaflongevity' 'leafLignin',",
                              "'hardsoft'. Only 'growthcurve' and 'mortalityshape' are used in this module.",
                              "Default is from Dominic Cyr and Yan Boulanger's applications of LANDIS-II"),
                 sourceURL = "https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv"),
    expectsInput(objectName = "speciesEcoregion", objectClass = "data.table",
                 desc = paste("Table of spatially-varying species traits ('maxB', 'maxANPP',",
                              "'establishprob'), defined by species and 'ecoregionGroup')",
                              "Defaults to a dummy table based on dummy data os biomass, age, ecoregion and land cover class")),
    expectsInput(objectName = "sppEquiv", objectClass = "data.table",
                 desc = "Table of species equivalencies. See `?LandR::sppEquivalencies_CA`."),
    expectsInput(objectName = "studyAreaANPP", objectClass = "SpatialPolygonsDataFrame",
                 desc = "Study area used to crop PSP data before building growth curves")
  ),
  outputObjects = bindrows(
    createsOutput("species", "data.table",
                  desc = "The updated invariant species traits table (see description for this object in inputs)"),
    createsOutput(objectName = "speciesEcoregion", "data.table",
                  desc = paste("The updated spatially-varying species traits table",
                               "(see description for this object in inputs)")),
    createsOutput(objectName = "speciesGAMMs", objectClass = "list",
                  desc = paste("A list of mixed-effect general additive models (GAMMs) for each tree species",
                               "modeling biomass as a function of age"))
  )
))

## event types
#   - type `init` is required for initialization

doEvent.Biomass_speciesParameters = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "Biomass_speciesParameters", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "Biomass_speciesParameters", "save")
    },
    plot = {

      #sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "Biomass_speciesParameters", "plot")

      # ! ----- STOP EDITING ----- ! #
    },
    save = {


      # sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "Biomass_speciesParameters", "save")

      # ! ----- STOP EDITING ----- ! #
    },

    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

## event functions
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
Init <- function(sim) {

  if (is.na(P(sim)$sppEquivCol)) {
    stop("Please supply sppEquivCol in parameters of Biomass_speciesParameters")
  }

  #prepare PSPdata
  speciesGAMMs <- Cache(makePSPgamms, studyAreaANPP = sim$studyAreaANPP,
                        PSPperiod = P(sim)$PSPperiod,
                        PSPgis = sim$PSPgis,
                        PSPmeasure = sim$PSPmeasure,
                        PSPplot = sim$PSPplot,
                        useHeight = P(sim)$useHeight,
                        biomassModel = P(sim)$biomassModel,
                        speciesCol = P(sim)$sppEquivCol,
                        sppEquiv = sim$sppEquiv,
                        NoOfIterations = P(sim)$GAMMiterations,
                        knots = P(sim)$GAMMknots,
                        minimumSampleSize = P(sim)$minimumPlotsPerGamm,
                        quantileAgeSubset = P(sim)$quantileAgeSubset,
                        minDBH = P(sim)$minDBH,
                        userTags = c(currentModule(sim), "makePSPgamms"))
  sim$speciesGAMMs <- speciesGAMMs

  classes <- lapply(sim$speciesGAMMs, FUN = "class")
  badModels <- classes[classes == "try-error"]

  if (!length(badModels) == 0) {
    message("convergence failures for these PSP growth curve models: ")
    print(names(badModels))
  }

  noData <- vapply(sim$speciesGAMMs[classes == "character"], FUN = function(x) {
    x == "insufficient data"
  }, FUN.VALUE = logical(1))

  if (any(noData)) {
    message("The following species did not have sufficient data for model estimation: ")
    print(names(noData))
  }

  modifiedSpeciesTables <- modifySpeciesTable(gamms = sim$speciesGAMMs,
                                              speciesTable = sim$species,
                                              factorialTraits = sim$factorialSpeciesTable,
                                              factorialBiomass = sim$reducedFactorialCohortData,
                                              sppEquiv = sim$sppEquiv,
                                              sppEquivCol = P(sim)$sppEquivCol,
                                              mortConstraints = P(sim)$constrainMortalityShape,
                                              growthConstraints = P(sim)$constrainGrowthCurve,
                                              mANPPconstraints = P(sim)$constrainMaxANPP)

  sim$species <- modifiedSpeciesTables

  modifiedSpeciesEcoregion <- modifySpeciesEcoregionTable(speciesEcoregion = sim$speciesEcoregion,
                                                          speciesTable = sim$species)
  sim$speciesEcoregion <- modifiedSpeciesEcoregion

  return(sim)
}

### template for save events
Save <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for plot events
plotFun <- function(sim) {

  # not sure we need to plot anything
  return(invisible(sim))
}


.inputObjects <- function(sim) {

  cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  if (!suppliedElsewhere("reducedFactorialCohortData", sim)) {
    sim$reducedFactorialCohortData <- prepInputs(targetFile = "reducedFactorialCD.Rdat",
                                                 destinationPath = dPath,
                                                 fun = "readRDS", overwrite = TRUE,
                                                 url = extractURL("reducedFactorialCohortData", sim),
                                                 useCache = TRUE, userTags = c(cacheTags, "reducedFactorial"))
  }

  if (!suppliedElsewhere("factorialSpeciesTable", sim)) {
    sim$factorialSpeciesTable <- prepInputs(targetFile = "factorialSpeciesTable.Rdat",
                                            destinationPath = dPath,
                                            url = extractURL("factorialSpeciesTable", sim),
                                            fun = "readRDS", overwrite = TRUE,
                                            useCache = TRUE, userTags = c(cacheTags, "factorialSpecies"))
  }

  if (!suppliedElsewhere("speciesEcoregion", sim)) {
    warning("generating dummy speciesEcoregion data - run Biomass_borealDataPrep for table with real speciesEcoregion attributes")
    sim$speciesEcoregion <- data.table(ecoregionGroup = "x",
                                       speciesCode = c("Abie_las", "Abie_bal", "Betu_pap", "Lari_lar", "Pice_eng",
                                                       "Pice_gla", "Pice_mar", "Pinu_ban",
                                                       "Pinu_con", "Pseu_men", "Popu_tre"),
                                       establishprob = 0.5, maxB = 5000, maxANPP = 5000/30, year = 0)
  }

  if (!suppliedElsewhere("species", sim)) {
    warning("generating dummy species data - run Biomass_borealDataPrep for table with real species attributes")
    sim$species <- data.table(species = c("Abie_las", "Abie_bal", "Betu_pap", "Lari_lar", "Pice_eng",
                                          "Pice_gla", "Pice_mar", "Pinu_ban",
                                          "Pinu_con", "Pseu_men", "Popu_tre"),
                              longevity = c(300, 300, 150, 140, 450, 400, 250, 150, 325, 600, 200),
                              mortalityshape = 15, growthcurve = 0)
  }

  if (!suppliedElsewhere("sppEquiv", sim)) {
    #pass a default sppEquivalencies_CA for common species in western Canada
    sppEquivalencies_CA <-  LandR::sppEquivalencies_CA
    sppEquivalencies_CA[, default := c(Pice_mar = "Pice_mar", Pice_gla = "Pice_gla",
                                       Pinu_con = "Pinu_con", Popu_tre = "Popu_tre",
                                       Betu_pap = "Betu_pap", Pice_eng = "Pice_eng",
                                       Pseu_men = "Pseu_men", Abie_bal = "Abie_bal",
                                       Pinu_ban = "Pinu_ban", Lari_lar = "Lari_lar")[LandR]]
    sppEquivalencies_CA[LANDIS_traits == "ABIE.LAS"]$default <- "Abie_las"
    sppEquivalencies_CA <- sppEquivalencies_CA[!LANDIS_traits == "PINU.CON.CON"]
    sppEquivalencies_CA <- sppEquivalencies_CA[!is.na(default)]
    sppEquivalencies_CA[LANDIS_traits == "ABIE.LAS", LandR := "Abie_las"]
    sim$sppEquiv <- sppEquivalencies_CA
  }

  if (!suppliedElsewhere("PSPmeasure", sim)) {
    sim$PSPmeasure <- Cache(prepInputs,
                            targetFile = "randomizedPSPmeasure.rds",
                            archive = "randomized_LandR_speciesParameters_Inputs.zip",
                            url =  extractURL("PSPmeasure", sim),
                            destinationPath = dPath,
                            fun = "readRDS")
  }
  if (!suppliedElsewhere("PSPplot", sim)) {
    sim$PSPplot <- Cache(prepInputs,
                         targetFile = "randomizedPSPplot.rds",
                         archive = "randomized_LandR_speciesParameters_Inputs.zip",
                         url = extractURL("PSPplot", sim),
                         destinationPath = dPath,
                         fun = "readRDS")

  }
  if (!suppliedElsewhere("PSPgis", sim)) {
    sim$PSPgis <- Cache(prepInputs,
                        targetFile = "randomizedPSPgis.rds",
                        archive = "randomized_LandR_speciesParameters_Inputs.zip",
                        url = extractURL("PSPgis", sim),
                        overwrite = TRUE,
                        destinationPath = dPath,
                        fun = "readRDS")
  }

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
### add additional events as needed by copy/pasting from above
