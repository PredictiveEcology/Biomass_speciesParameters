
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "LandR_speciesParameters",
  description = NA, #"insert module description here",
  keywords = NA, # c("insert key words here"),
  authors = c(person(c("Ian"), "Eddy", email = "ian.eddy@example.com", role = c("aut", "cre"))),
  childModules = character(0),
  version = list(SpaDES.core = "0.2.6", LandR_speciesParameters = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "LandR_speciesParameters.Rmd"),
  reqdPkgs = list("PredictiveEcology/pemisc@development", 'mgcv', 'fpCompare'),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, 
                    desc = paste("Should this entire module be run with caching activated?",
                          "This is generally intended for data-type modules, where stochasticity and time are not relevant")),
    defineParameter("PSPperiod", "numeric", c(1958, 2011), NA, NA, 
                    desc = paste("The years by which to compute climate normals and subset sampling plot data. Must be a vector of at least length 2")),
    defineParameter("sppEquivCol", "character", NA, NA, NA,
                    paste("The column in sim$specieEquivalency data.table to group species by. This parameter should share the same",
                          "name as in Boreal_LBMRDataPrep. PSPs are aggregated by names in the PSP column and traits estimated",
                          "for the corresponding names in the sppEquivCol")),
    defineParameter("useHeight", "logical", FALSE, NA, NA, 
                    desc = paste("Should height be used to calculate biomass (in addition to DBH).
                    Advise against including height unless you work ONLY in BC")),
    defineParameter("biomassModel", "character", "Lambert2005", NA, NA, 
                    desc =  paste("The model used to calculate biomass from DBH. Can be either 'Lambert2005' or 'Ung2008'")),
    defineParameter("constrainGrowthCurve", "numeric", c(0.5, 0.5), 0, 1, 
                    desc = paste("constraints on range of potential growth curves when fitting traits. This parameter exists",
                                 "because growth curve is confounded by aNPPproportion when estimating traits")),
    defineParameter("constrainMortalityShape", 'numeric', c(15, 25), 5, 25, 
                    desc = paste("constraints on mortality shape when fitting traits. low mortality curve needs to excessive",
                                 "cohorts with very little biomass as longevity is approached, adding computation strain"))
  ),
  inputObjects = bind_rows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput(objectName = "factorialSpeciesTable", objectClass = "data.table", 
                 desc = paste("table with species traits for matching to factorialCohortData"),
                 sourceURL = "https://drive.google.com/open?id=1u0aXpJKbFYisEPjczTHYfXu6xZjz9O_8"),
    expectsInput(objectName = "reducedFactorialCohortData", objectClass = "data.table", 
                 desc = paste("results of factorial species trait simulation. This can be found by running",
                              "the script in landRFactorial.R but you need a special LBMR"), 
                 sourceURL = "https://drive.google.com/open?id=1rv_ROwBhCmu4K6ZiicsphjKgsHthrmBn"),
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
                 desc = "table of species equivalencies. See LandR::sppEquivalencies_CA.",
                 sourceURL = ""),
    expectsInput(objectName = "studyAreaANPP", "SpatialPolygonsDataFrame", desc = "study area used to crop PSP data before building growth curves")
  ),
  outputObjects = bind_rows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = 'speciesGAMMs', objectClass = 'list', 
                  desc = paste('a list of mixed-effect general additive models (gamm) for each tree species',
                               'modeling biomass as a function of age')),
    createsOutput("species", "data.table",
                  desc = "a table that has species traits such as longevity..."),
  )
))

## event types
#   - type `init` is required for initialization

doEvent.LandR_speciesParameters = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "LandR_speciesParameters", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "LandR_speciesParameters", "save")
    },
    plot = {
     
      #sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "LandR_speciesParameters", "plot")

      # ! ----- STOP EDITING ----- ! #
    },
    save = {


    # sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "LandR_speciesParameters", "save")

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

  #prepare PSPdata
  psp <- prepPSPaNPP(studyAreaANPP = sim$studyAreaANPP, PSPperiod = P(sim)$PSPperiod,
                     PSPgis = sim$PSPgis, PSPmeasure = sim$PSPmeasure, PSPplot = sim$PSPplot,
                     useHeight = P(sim)$useHeight, biomassModel = P(sim)$biomassModel)

  sim$speciesGAMMs <- buildGrowthCurves(PSPdata = psp, 
                                        speciesCol = P(sim)$sppEquivCol,
                                        sppEquiv = sim$sppEquiv)
  classes <- lapply(sim$speciesGAMMs, FUN = 'class')
  badModels <- classes[classes == 'try-error']
  if (!is.null(badModels)) {
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
  browser()
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
                                                 destinationPath = dataPath(sim),
                                                 fun = "readRDS",
                                                 url = extractURL('reducedFactorialCohortData', sim),
                                                 useCache = TRUE, userTags = c(cacheTags, "reducedFactorial"))
  }
  
  if (!suppliedElsewhere("factorialSpeciesTable", sim)) {
    sim$factorialSpeciesTable <- prepInputs(targetFile = "factorialSpeciesTable.Rdat",
                                            destinationPath = dataPath(sim),
                                            url = extractURL('factorialSpeciesTable', sim), 
                                            fun = "readRDS",
                                            useCache = TRUE, userTags = c(cacheTags, "factorialSpecies"))
  }
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
### add additional events as needed by copy/pasting from above
