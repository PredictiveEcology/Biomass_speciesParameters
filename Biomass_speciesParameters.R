
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "Biomass_speciesParameters",
  description = "For estimating LANDIS-II species traits from PSP-derived growth curves",
  keywords = NA, # c("insert key words here"),
  authors = c(person(c("Ian"), "Eddy", email = "ian.eddy@example.com", role = c("aut", "cre"))),
  childModules = character(0),
  version = list(SpaDES.core = "0.2.6", Biomass_speciesParameters = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_speciesParameters.Rmd"),
  reqdPkgs = list("PredictiveEcology/pemisc@development", 'mgcv', 'nlme', 'fpCompare',
                  'PredictiveEcology/LandR@development', 'crayon'),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    desc = paste("Should this entire module be run with caching activated?",
                                 "This is generally intended for data-type modules, where stochasticity and time are not relevant")),
    defineParameter("biomassModel", "character", "Lambert2005", NA, NA,
                    desc =  paste("The model used to calculate biomass from DBH. Can be either 'Lambert2005' or 'Ung2008'")),
    defineParameter("constrainGrowthCurve", "numeric", c(0.5, 0.5), 0, 1,
                    desc = paste("constraints on range of potential growth curves when fitting traits. This module accepts a",
                                 "list of vectors, with names equal to sppEquivCol, so that traits are customizable")),
    defineParameter("constrainMortalityShape", 'numeric', c(15, 25), 5, 25,
                    desc = paste('constraints on mortality shape when fitting traits. low mortality curve needs to excessive',
                                 'cohorts with very little biomass as longevity is approached, adding computation strain.',
                                 'This module accepts a list of vectors, with names equal to sppEquivCol',
                                 'so that GAMMS are customizable')),
    defineParameter("GAMMiterations", "numeric", 8, 1, NA,
                    desc = paste("number of iterations for GAMMs. This module accepts a",
                                 "list of vectors, with names equal to sppEquivCol, so that GAMMS are customizable")),
    defineParameter("GAMMknots", "numeric", 3, NA, NA,
                    desc = paste("the number of knots to use in the GAMM. Either 3 or 4 is recommended. This module accepts a",
                                 "list of vectors, with names equal to sppEquivCol, so that GAMMS are customizable")),
    defineParameter("minimumPlotsPerGamm", "numeric", 50, 10, NA, desc = paste("minimum number of PSP plots before building GAMM")),
    defineParameter("PSPperiod", "numeric", c(1920, 2019), NA, NA,
                    desc = paste("The years by which to subset sample plot data, if desired. Must be a vector of length 2")),
    defineParameter("quantileAgeSubset", "numeric", 95, 1, 100,
                    desc = paste("quantile by which to subset PSP data. As older stands are sparsely represented, the oldest measurements",
                                 "become vastly more influential. This parameter accepts both a single value and a list of vectors",
                                 "named by sppEquivCol. The PSP stand ages are found in sim$speciesGAMMs$<species>$originalData")),
    defineParameter("sppEquivCol", "character", 'default', NA, NA,
                    paste("The column in sim$sppEquiv data.table to group species by. This parameter should share the same",
                          "name as in Biomass_borealDataPrep. PSPs are aggregated by names in the PSP column and traits estimated",
                          "for the corresponding names in the sppEquivCol")),
    defineParameter("useHeight", "logical", FALSE, NA, NA,
                    desc = paste("Should height be used to calculate biomass (in addition to DBH).
                    Advise against including height unless you are certain it is present in every PSP"))
  ),
  inputObjects = bindrows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput(objectName = "factorialSpeciesTable", objectClass = "data.table",
                 desc = paste("table with species traits for matching to factorialCohortData"),
                 sourceURL = "https://drive.google.com/open?id=1q0ou0CBzD9GqGSparpHqf318IWK6ycty"),
    expectsInput(objectName = "reducedFactorialCohortData", objectClass = "data.table",
                 desc = paste("results of factorial species trait simulation. This can be found by running",
                              "SpeciesFactorial.R but requires a specific commit of Boreal_Biomass"),
                 sourceURL = "https://drive.google.com/open?id=1h8StXE0vm8xyDycRomCkwIaL7wfh5Irj"),
    expectsInput(objectName = "PSPmeasure", objectClass = "data.table", desc = "merged PSP and TSP individual measurements"),
    expectsInput(objectName = "PSPplot", objectClass = "data.table", desc = "merged PSP and TSP plot data"),
    expectsInput(objectName = "PSPgis", objectClass = "sf", desc = "Plot location sf object. Contains duplicates"),
    expectsInput(objectName = "species", "data.table",
                 desc = paste("a table that has species traits such as longevity, shade tolerance, etc.",
                              "Default is partially based on Dominic Cir and Yan's project"),
                 sourceURL = "https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv"),
    expectsInput(objectName = "speciesEcoregion", "data.table",
                 desc = paste("table defining the maxANPP, maxB and SEP, which can change with both ecoregion and simulation time.",
                              "Defaults to a dummy table based on dummy data os biomass, age, ecoregion and land cover class")),
    expectsInput(objectName = "sppEquiv", "data.table",
                 desc = "table of species equivalencies. See LandR::sppEquivalencies_CA."),
    expectsInput(objectName = "studyAreaANPP", "SpatialPolygonsDataFrame",
                 desc = "study area used to crop PSP data before building growth curves")
  ),
  outputObjects = bindrows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = "speciesEcoregion", "data.table",
                  desc = paste("table defining the maxANPP, maxB and SEP, which can change with both ecoregion and simulation time.",
                               "Defaults to a dummy table based on dummy data os biomass, age, ecoregion and land cover class")),
    createsOutput(objectName = 'speciesGAMMs', objectClass = 'list',
                  desc = paste('a list of mixed-effect general additive models (gamm) for each tree species',
                               'modeling biomass as a function of age')),
    createsOutput("species", "data.table",
                  desc = "a table that has species traits such as longevity..."),
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
                        userTags = c(currentModule(sim), "makePSPgamms"))
  sim$speciesGAMMs <- speciesGAMMs

  classes <- lapply(sim$speciesGAMMs, FUN = 'class')
  badModels <- classes[classes == 'try-error']

  if (!length(badModels) == 0) {
    message("convergence failures for these PSP growth curve models: ")
    print(names(badModels))
  }

  modifiedSpeciesTables <- modifySpeciesTable(gamms = sim$speciesGAMMs,
                                              speciesTable = sim$species,
                                              factorialTraits = sim$factorialSpeciesTable,
                                              factorialBiomass = sim$reducedFactorialCohortData,
                                              sppEquiv = sim$sppEquiv,
                                              sppEquivCol = P(sim)$sppEquivCol,
                                              mortConstraints = P(sim)$constrainMortalityShape,
                                              growthConstraints = P(sim)$constrainGrowthCurve)

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
                                                 url = extractURL('reducedFactorialCohortData', sim),
                                                 useCache = TRUE, userTags = c(cacheTags, "reducedFactorial"))
  }

  if (!suppliedElsewhere("factorialSpeciesTable", sim)) {
    sim$factorialSpeciesTable <- prepInputs(targetFile = "factorialSpeciesTable.Rdat",
                                            destinationPath = dPath,
                                            url = extractURL('factorialSpeciesTable', sim),
                                            fun = "readRDS", overwrite = TRUE,
                                            useCache = TRUE, userTags = c(cacheTags, "factorialSpecies"))
  }

  if (!suppliedElsewhere("speciesEcoregion", sim)) {
    sim$speciesEcoregion <- data.table(ecoregionGroup = "x",
                                       speciesCode = c("Abie_las", 'Abie_bal', 'Betu_pap', 'Lari_lar', 'Pice_eng',
                                                       'Pice_gla', 'Pice_mar', 'Pinu_ban',
                                                       'Pinu_con', 'Pseu_men', "Popu_tre"),
                                       establishprob = 0.5, maxB = 5000, maxANPP = 5000/30, year = 0)
  }

  if (!suppliedElsewhere("species", sim)) {
    sim$species <- data.table(species = c("Abie_las", 'Abie_bal', 'Betu_pap', 'Lari_lar', 'Pice_eng',
                                          'Pice_gla', 'Pice_mar', 'Pinu_ban',
                                          'Pinu_con', 'Pseu_men', "Popu_tre"),
                              longevity = c(300, 300, 150, 150, 450, 400, 250, 150, 325, 600, 200),
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

  if (!suppliedElsewhere("PSPmeasure", sim) |
      !suppliedElsewhere("PSPplot", sim) |
      !suppliedElsewhere("PSPgis", sim)) {
    stop("This module requires PSP objects. Try including the module PSP_Clean or contact ian.eddy@canada.ca")
  }
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
### add additional events as needed by copy/pasting from above
