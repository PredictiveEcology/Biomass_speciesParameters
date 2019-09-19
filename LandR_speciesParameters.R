
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
  reqdPkgs = list(),
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
    defineParameter("useHeight", "logical", FALSE, NA, NA, 
                    desc = paste("Should height be used to calculate biomass (in addition to DBH).
                    Don't use if studyAreaPSP includes Alberta")),
    defineParameter("biomassModel", "character", "Lambert2005", NA, NA, 
                    desc =  paste("The model used to calculate biomass from DBH. Can be either 'Lambert2005' or 'Ung2008'"))
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
    expectsInput("species", "data.table",
                 desc = paste("a table that has species traits such as longevity, shade tolerance, etc.",
                              "Default is partially based on Dominic Cir and Yan's project"),
                 sourceURL = "https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv"),
    expectsInput("speciesEcoregion", "data.table",
                 desc = paste("table defining the maxANPP, maxB and SEP, which can change with both ecoregion and simulation time.",
                              "Defaults to a dummy table based on dummy data os biomass, age, ecoregion and land cover class")),
    expectsInput("studyAreaANPP", "SpatialPolygonsDataFrame", desc = "study area used to crop PSP data before building growth curves")
  ),
  outputObjects = bind_rows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = NA, objectClass = NA, desc = NA)
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
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      #plotFun(sim) # uncomment this, replace with object to plot
      # schedule future event(s)

      # e.g.,
      #sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "LandR_speciesParameters", "plot")

      # ! ----- STOP EDITING ----- ! #
    },
    save = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "LandR_speciesParameters", "save")

      # ! ----- STOP EDITING ----- ! #
    },
    event1 = {
      # ! ----- EDIT BELOW ----- ! #

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
  PSPanPPP <- prepPSPaNPP(studyAreaANPP = sim$studyAreaANPP, PSPperiod = P(sim)$PSPperiod,
                          PSPgis = sim$PSPgis, PSPmeasure = sim$PSPmeasure, PSPplot = sim$PSPplot,
                          useHeight = P(sim)$useHeight, biomassModel = P(sim)$biomasssModel)
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
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  #Plot(sim$object)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

prepPSPaNPP <- function(studyAreaANPP, PSPgis, PSPmeasure, PSPplot,
                          useHeight, biomassModel, PSPperiod) {
  
  #Crop points to studyArea
  tempSA <- spTransform(x = studyAreaANPP, CRSobj = crs(PSPgis)) %>%
    st_as_sf(.)
  message(yellow("Filtering PSPs for ANPP to study Area..."))
  PSP_sa <- PSPgis[tempSA,] %>% #Find how to cache this. '[' did not work
    setkey(., OrigPlotID1)
  message(yellow(paste0("There are "), nrow(PSP_sa), " PSPs in your study area"))
  
  #Filter other PSP datasets to those in study Area
  PSPmeasure <- PSPmeasure[OrigPlotID1 %in% PSP_sa$OrigPlotID1,]
  PSPplot <- PSPplot[OrigPlotID1 %in% PSP_sa$OrigPlotID1,]
  
  #length(PSPclimData)/length(PSP_sa) should always yield a whole number.
  #Filter data by study period
  message(yellow("Filtering PSPs for ANPP by study period..."))
  PSPmeasure <- PSPmeasure[MeasureYear > min(PSPperiod) &
                             MeasureYear < max(PSPperiod),]
  PSPplot <- PSPplot[MeasureYear > min(PSPperiod) &
                       MeasureYear < max(PSPperiod),]
  
  #Join data (should be small enough by now)
  PSPmeasure <- PSPmeasure[PSPplot, on = c('MeasureID', 'OrigPlotID1', 'MeasureYear')]
  PSPmeasure[, c('Longitude', 'Latitude', 'Easting', 'Northing', 'Zone'):= NULL]
  
  #Filter by > 30 trees at first measurement (P) to ensure forest.
  forestPlots <- PSPmeasure[, .(measures = .N), OrigPlotID1] %>%
    .[measures >= 30,]
  PSPmeasure <- PSPmeasure[OrigPlotID1 %in% forestPlots$OrigPlotID1,]
  PSPplot <- PSPplot[OrigPlotID1 %in% PSPmeasure$OrigPlotID1,]
  
  #Restrict to trees > 10 DBH (P) This gets rid of some big trees. Some 15 metres tall
  PSPmeasure <- PSPmeasure[DBH >= 10,]
  

  #Calculate biomass
  tempOut <- biomassCalculation(species = PSPmeasure$newSpeciesName,
                                DBH = PSPmeasure$DBH,
                                height = PSPmeasure$Height,
                                includeHeight = useHeight,
                                equationSource = biomassModel)
  message(yellow("No PSP biomass estimate possible for these species: "))
  print(tempOut$missedSpecies)
  PSPmeasure$biomass <- tempOut$biomass
  PSPmeasure <- PSPmeasure[biomass != 0]
 
  #Need to join the psps with PSPgis to get the actual coords with same crs
  browser() #you need to join plots and measures, to get age
  PSPmeasure <- PSPgis[psp, on = c("OrigPlotID1" = 'OrigPlotID1')]
  #remove some redundant columns from before the merge + useless GIS
  psp[, c('geometry', 'Latitude', 'Longitude', 'Easting', 'Northing','i.MeasureYear', 'i.OrigPlotID1') := NULL]
  
  #standardize biomass by plotsize
  #Note: there are still obvious errors in PSP data. e.g. plotSize 0.0029 ha, but 90 trees? yeah right buddy
  psp$Biomass <- psp$Biomass/psp$PlotSize/10 #puts in LandR units of g/m2 and scales by plot size
  psp$newSpeciesName <- as.factor(psp$newSpeciesName)
  psp <- psp[standAge > 0]
  psp <- psp[!is.na(Biomass)]
  
  return(PSPmodelData)
}


### template for your event1
Event1 <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event1Test1 <- " this is test for event 1. " # for dummy unit test
  # sim$event1Test2 <- 999 # for dummy unit test


  # ! ----- STOP EDITING ----- ! #
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
