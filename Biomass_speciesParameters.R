defineModule(sim, list(
  name = "Biomass_speciesParameters",
  description = paste("For estimating LANDIS-II species traits based on growth curves derived",
                      "from Permanent Sample Plot (PSP) and Temporary Sample Plot (TSP) data"),
  keywords = NA, # c("insert key words here"),
  authors = c(person(c("Ian"), "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = c("aut", "cre")),
              person(c("Eliot"), "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut")),
              person(c("Ceres"), "Barros", email = "ceres.barros@ubc.ca", role = c("ctb"))),
  childModules = character(0),
  version = list(Biomass_speciesParameters = "2.0.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_speciesParameters.Rmd"),
  loadOrder = list(after = c("Biomass_speciesFactorial", "Biomass_borealDataPrep"),
                   before = c("Biomass_core")),
  reqdPkgs = list("crayon", "data.table", "disk.frame", "fpCompare", "ggplot2", "gridExtra",
                  "magrittr", "mgcv", "nlme", "purrr", "robustbase", "sf",
                  "PredictiveEcology/LandR@development (>= 1.1.0.9054)",
                  "PredictiveEcology/pemisc@development (>= 0.0.3.9002)",
                  "PredictiveEcology/reproducible@development (>= 1.2.10.9001)",
                  "PredictiveEcology/SpaDES.core@development (>= 2.0.2.9004)",
                  "ianmseddy/PSPclean@development (>= 0.1.3.9003)"),
  parameters = rbind(
    defineParameter("biomassModel", "character", "Lambert2005", NA, NA,
                    desc =  paste("The model used to calculate biomass from DBH. Can be either 'Lambert2005' or 'Ung2008'.")),
    defineParameter("maxBInFactorial", "integer", 5000L, NA, NA,
                    desc = paste("The arbitrary maximum biomass for the factorial simulations.",
                                 "This is a per-species maximum within a pixel")),
    defineParameter("minimumPlots", "numeric", 50, 10, NA,
                    desc = paste("Minimum number of PSP plots per species")),
    defineParameter("minDBH", "integer", 0L, 0L, NA,
                    desc = paste("Minimum diameter at breast height (DBH) in cm used to filter PSP data.",
                                 "Defaults to 0 cm, i.e. all tree measurements are used.")),
    defineParameter("PSPdataTypes", "character", "all", NA, NA,
                    desc = paste("Which PSP datasets to source, defaulting to all. Other available options include",
                                 "'BC', 'AB', 'SK', 'ON', 'NFI', and 'dummy'. 'dummy' should be used for unauthorized users.")),
    defineParameter("PSPperiod", "numeric", c(1920, 2019), NA, NA,
                    desc = paste("The years by which to subset sample plot data, if desired. Must be a vector of length 2")),
    defineParameter("speciesFittingApproach", "character", "focal", NA, NA,
                    desc =  paste("Either 'all', 'pairwise', 'focal' or 'single', indicating whether to pool ",
                                  "all species into one fit, do pairwise species (for multiple cohort situations), do",
                                  "pairwise species, but using a focal species approach where all other species are ",
                                  "pooled into 'other' or do one species at a time. If 'all', all species will have",
                                  "identical species-level traits")),
    defineParameter("sppEquivCol", "character", "Boreal", NA, NA,
                    paste("The column in `sim$sppEquiv` data.table that defines individual species.",
                          "The names should match those in the species table.")),
    defineParameter("standAgesForFitting", "integer", c(0L, 150L), NA, NA,
                    desc = paste("The minimum and maximum ages of the biomass-by-age curves used in fitting.",
                                 "It is generally recommended to keep this param under 200, given the low data",
                                 "availability of stands aged 200+, with some exceptions.")),
    defineParameter("useHeight", "logical", TRUE, NA, NA,
                    desc = paste("Should height be used to calculate biomass (in addition to DBH).",
                                 "DBH is used by itself when height is missing.")),
    defineParameter(".plots", "character", "screen", NA, NA,
                    desc = "Used by Plots function, which can be optionally used here"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    desc = "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time interval between save events"),
    defineParameter(".studyAreaName", "character", NA, NA, NA,  
                    desc = paste("Human-readable name for the growth curve filename.",
                                 "If `NA`, a hash of sppEquiv[[sppEquivCol]] will be used.")),
    defineParameter(".useCache", "character", c(".inputObjects", "init"), NA, NA,
                    desc = paste("Should this entire module be run with caching activated?",
                                 "This is generally intended for data-type modules, where stochasticity and time are not relevant"))
  ),
  inputObjects = bindrows(
    expectsInput(objectName = "speciesTableFactorial", objectClass = "data.table",
                 desc = paste("A large species table (sensu Biomass_core) with all columns used by",
                              "Biomass_core, e.g., longevity, growthcurve, mortalityshape, etc., when",
                              "it was used to generate `cohortDataFactorial`.",
                              "See PredictiveEcology/Biomass_factorial for futher information.",
                              "It will be written to `disk.frame` following sim completion, to preserve RAM."),
                 sourceURL = "https://drive.google.com/file/d/1NH7OpAnWtLyO8JVnhwdMJakOyapBnuBH/"),
    expectsInput(objectName = "cohortDataFactorial", objectClass = "data.table",
                 desc = paste("A large `cohortData` table (sensu Biomass_core) with columns age, B, and speciesCode",
                              "that joins with `speciesTableFactorial`. See PredictiveEcology/Biomass_factorial",
                              "for further information. It will be written to `disk.frame` following sim",
                              "completion, to preserve RAM."),
                 sourceURL = "https://drive.google.com/file/d/1NH7OpAnWtLyO8JVnhwdMJakOyapBnuBH/"),
    expectsInput(objectName = "PSPmeasure_sppParams", objectClass = "data.table",
                 desc = paste("Merged PSP and TSP individual tree measurements. Must include the following columns:",
                              "'MeasureID', 'OrigPlotID1', 'MeasureYear', 'TreeNumber', 'Species', 'DBH' and 'newSpeciesName'",
                              "the latter corresponding to species names in `LandR::sppEquivalencies_CA$PSP`.",
                              "Defaults to randomized PSP data stripped of real plotIDs"),
                 sourceURL = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing"),
    expectsInput(objectName = "PSPplot_sppParams", objectClass = "data.table",
                 desc = paste("Merged PSP and TSP plot data. Defaults to randomized PSP data stripped of real plotIDs.",
                              "Must contain columns 'MeasureID', 'MeasureYear', 'OrigPlotID1', and 'baseSA',",
                              "the latter being stand age at year of first measurement"),
                 sourceURL = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing"),
    expectsInput(objectName = "PSPgis_sppParams", objectClass = "sf",
                 desc = paste("Plot location `sf` object. Defaults to PSP data stripped of real plotIDs/location.",
                              "Must include field 'OrigPlotID1' for joining to PSPplot object"),
                 sourceURL = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing"),
    expectsInput(objectName = "species", objectClass = "data.table",
                 desc = paste("A table of invariant species traits with the following trait colums:",
                              "'species', 'Area', 'longevity', 'sexualmature', 'shadetolerance',",
                              "'firetolerance', 'seeddistance_eff', 'seeddistance_max', 'resproutprob',",
                              "'mortalityshape', 'growthcurve', 'resproutage_min', 'resproutage_max',",
                              "'postfireregen', 'wooddecayrate', 'leaflongevity' 'leafLignin', and 'hardsoft'.",
                              "Only 'growthcurve', 'hardsoft',  and 'mortalityshape' are used in this module.",
                              "Default is from Dominic Cyr and Yan Boulanger's applications of LANDIS-II"),
                 sourceURL = "https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv"),
    expectsInput(objectName = "speciesEcoregion", objectClass = "data.table",
                 desc = paste("Table of spatially-varying species traits ('maxB', 'maxANPP',",
                              "'establishprob'), defined by species and 'ecoregionGroup')",
                              "Defaults to a dummy table based on dummy data os biomass, age, ecoregion and land cover class")),
    expectsInput(objectName = "sppEquiv", objectClass = "data.table",
                 desc = "Table of species equivalencies. See `?LandR::sppEquivalencies_CA`."),
    expectsInput(objectName = "studyAreaANPP", objectClass = "sf",
                 desc = paste("Optional study area used to crop PSP data before building growth curves.",
                              "If supplied, an ecoregion-scale object is recommended, at a minimum."))
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "cohortDataFactorial", objectClass = "disk.frame",
                  desc = "This object is converted to a `disk.frame` to save memory. Read using `as.data.table()`."),
    createsOutput("species", "data.table",
                  desc = "The updated invariant species traits table (see description for this object in inputs)"),
    createsOutput(objectName = "speciesEcoregion", "data.table",
                  desc = paste("The updated spatially-varying species traits table",
                               "(see description for this object in inputs)")),
    createsOutput(objectName = "speciesTableFactorial", objectClass = "disk.frame",
                  desc = "This object is converted to a `disk.frame` to save memory. Read using `as.data.table()`.")
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

      # build growth curves if applicable
      sim <- Init(sim)
      #update tables

      # schedule future event(s)
      # sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "Biomass_speciesParameters", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "Biomass_speciesParameters", "save")
      sim <- scheduleEvent(sim, start(sim), "Biomass_speciesParameters",
                           "updateSpeciesTables", eventPriority = 1)
      sim <- scheduleEvent(sim, start(sim), "Biomass_speciesParameters",
                           "writeFactorialToDisk", eventPriority = 2)
    },

    updateSpeciesTables = {
      sim <- updateSpeciesTables(sim)
    },

    writeFactorialToDisk = {
      sim <- useDiskFrame(sim)
    },

    plot = {
      ## plotting happens in Init - it could be moved if relevant objects are assigned to mod
    },
    save = {
      # sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "Biomass_speciesParameters", "save")

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
  origDTthreads <- data.table::getDTthreads()
  if (getDTthreads() > 4) {
    data.table::setDTthreads(4)
  }
  on.exit(data.table::setDTthreads(origDTthreads))

  ## if no PSP data supplied, simList returned unchanged

  if (all(P(sim)$PSPdataTypes != "none")) {
    if (is.na(P(sim)$sppEquivCol)) {
      stop("Please supply 'sppEquivCol' in parameters of Biomass_speciesParameters.")
    }

    paramCheckOtherMods(sim, "maxBInFactorial")
    paramCheckOtherMods(sim, paramToCheck = "sppEquivCol", ifSetButDifferent = "error")

    #find the max biomass achieved by each species when growing with no competition
    tempMaxB <- sim$cohortDataFactorial[age == 1, .N, .(pixelGroup)]
    #take the pixelGroups with only 1 species at start of factorial
    tempMaxB <- tempMaxB[N == 1,]
    tempMaxB <- sim$cohortDataFactorial[pixelGroup %in% tempMaxB$pixelGroup,
                                        .(inflationFactor = P(sim)$maxBInFactorial/max(B)),
                                        , .(pixelGroup, speciesCode)]
    # in some cases the speciesTableFactorial doesn't have "species" column; just "speciesCode"
    if (!("species" %in% colnames(sim$speciesTableFactorial))) { # TODO this is a work around -- the speciesTableFactorial should be stable
      setnames(sim$speciesTableFactorial, old = "speciesCode", new = "species")
    }
    tempMaxB <- sim$speciesTableFactorial[tempMaxB, on = c("species" = "speciesCode", "pixelGroup")]
    #pair-wise species will be matched with traits, as the species code won't match
    tempMaxB <- tempMaxB[, .(species, longevity, growthcurve, mortalityshape, mANPPproportion, inflationFactor)]
    gc()
    #prepare PSPdata
    speciesGrowthCurves <- buildGrowthCurves_Wrapper(
      studyAreaANPP = sim$studyAreaANPP,
      PSPperiod = P(sim)$PSPperiod,
      PSPgis = sim$PSPgis_sppParams,
      PSPmeasure = sim$PSPmeasure_sppParams,
      PSPplot = sim$PSPplot_sppParams,
      useHeight = P(sim)$useHeight,
      biomassModel = P(sim)$biomassModel,
      speciesCol = P(sim)$sppEquivCol,
      sppEquiv = sim$sppEquiv,
      minimumSampleSize = P(sim)$minimumPlots,
      minDBH = P(sim)$minDBH,
      speciesFittingApproach = P(sim)$speciesFittingApproach)
     # userTags = c(currentModule(sim), "buildGrowthCurves_Wrapper"))
    
    if (is.na(P(sim)$.studyAreaName)) {
      studyAreaName <- reproducible::studyAreaName(sim$sppEquiv[[P(sim)$sppEquivCol]])
    } else {
      studyAreaName <- P(sim)$.studyAreaName
    }

    saveRDS(speciesGrowthCurves, file.path(outputPath(sim), paste0("speciesGrowthCurves_", studyAreaName, ".rds")))
    gc()
    classes <- lapply(speciesGrowthCurves, FUN = "class")

    noData <- vapply(speciesGrowthCurves[classes == "character"], FUN = function(x) {
      x == "insufficient data"
    }, FUN.VALUE = logical(1))

    if (any(noData)) {
      message("The following species did not have sufficient data for model estimation: ")
      print(names(noData))
    }
    speciesWithNewlyEstimated <- unique(unlist(strsplit(names(speciesGrowthCurves), "__")))
    speciesWithoutNewlyEstimated <- setdiff(sim$sppEquiv[[Par$sppEquivCol]], speciesWithNewlyEstimated)
    if (length(speciesWithoutNewlyEstimated)) {
      message(crayon::yellow(paste(speciesWithoutNewlyEstimated, collapse = ", "),
                             "have insufficient data to estimate species parameters; using original user supplied"))
    }

    modifiedSpeciesTables <- modifySpeciesTable(
      GCs = speciesGrowthCurves,
      speciesTable = sim$species,
      factorialTraits = setDT(sim$speciesTableFactorial),
      # setDT to deal with reload from Cache (no effect otherwise)
      factorialBiomass = setDT(sim$cohortDataFactorial),
      # setDT to deal with reload from Cache (no effect otherwise)
      sppEquiv = sim$sppEquiv,
      approach = P(sim)$speciesFittingApproach,
      sppEquivCol = P(sim)$sppEquivCol,
      maxBInFactorial = P(sim)$maxBInFactorial,
      inflationFactorKey = tempMaxB,
      standAgesForFitting = P(sim)$standAgesForFitting)
    
    gg <- modifiedSpeciesTables$gg
    Plots(gg, usePlot = FALSE, fn = print, ggsaveArgs = list(width = 10, height = 7),
          filename = paste("LandR_VS_NLM_growthCurves"))
    sim$species <- modifiedSpeciesTables$best
  } else {
    message("P(sim)$PSPdataTypes is 'none' -- bypassing species traits estimation from PSP data.")
  }
  return(sim)
}

updateSpeciesTables <- function(sim) {
  modifiedTables <- modifySpeciesAndSpeciesEcoregionTable(speciesEcoregion = sim$speciesEcoregion,
                                                          speciesTable = sim$species)
  sim$speciesEcoregion <- modifiedTables$newSpeciesEcoregion
  sim$species <- modifiedTables$newSpeciesTable
  return(sim)
}

useDiskFrame <- function(sim) {
  cdRows <- nrow(sim$cohortDataFactorial)
  # the rows of a factorial object will determine whether it is unique in 99.9% of cases
  sim$cohortDataFactorial <- as.disk.frame(sim$cohortDataFactorial, overwrite = TRUE,
                                           outdir = file.path(dataPath(sim),
                                                              paste0("cohortDataFactorial", cdRows)))
  stRows <- nrow(sim$speciesTableFactorial)
  sim$speciesTableFactorial <- as.disk.frame(sim$speciesTableFactorial, overwrite = TRUE,
                                             outdir = file.path(dataPath(sim),
                                                                paste0("speciesTableFactorial", stRows)))
  ## NOTE: disk.frame objects can be converted to data.table with as.data.table
  gc()
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

.inputObjects <- function(sim) {
  origDTthreads <- data.table::getDTthreads()
  if (getDTthreads() > 4) {
    data.table::setDTthreads(4)
  }
  on.exit(data.table::setDTthreads(origDTthreads))

  cacheTags <- c(currentModule(sim), "OtherFunction:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  if (!suppliedElsewhere("cohortDataFactorial", sim)) {
    sim$cohortDataFactorial <- prepInputs(targetFile = "cohortDataFactorial_medium.rds",
                                          destinationPath = dPath,
                                          fun = "readRDS", overwrite = TRUE,
                                          url = extractURL("cohortDataFactorial", sim),
                                          useCache = TRUE, userTags = c(cacheTags, "factorialCohort"))
  }

  if (!suppliedElsewhere("speciesTableFactorial", sim)) {
    sim$speciesTableFactorial <- prepInputs(targetFile = "speciesTableFactorial_medium.rds",
                                            destinationPath = dPath,
                                            url = extractURL("speciesTableFactorial", sim),
                                            fun = "readRDS", overwrite = TRUE,
                                            useCache = TRUE, userTags = c(cacheTags, "factorialSpecies"))
  }

  if (!suppliedElsewhere("sppEquiv", sim)) {
    #pass a default sppEquivalencies_CA for common species in western Canada
    sppEquiv <- LandR::sppEquivalencies_CA
    sim$sppEquiv <- sppEquiv[LandR %in% c(Pice_mar = "Pice_mar", Pice_gla = "Pice_gla",
                                          Pinu_con = "Pinu_con", Popu_tre = "Popu_tre",
                                          Betu_pap = "Betu_pap", Pice_eng = "Pice_eng",
                                          Abie_bal = "Abie_bal",
                                          Pinu_ban = "Pinu_ban", Lari_lar = "Lari_lar"), ]
  }

  if (!suppliedElsewhere("speciesEcoregion", sim)) {
    warning("generating dummy speciesEcoregion data - run Biomass_borealDataPrep for table with real speciesEcoregion attributes")
    sim$speciesEcoregion <- data.table(
      speciesCode = unique(sim$sppEquiv[[P(sim)$sppEquivCol]]),
      ecoregionGroup = "x",
      establishprob = 0.5,
      maxB = P(sim)$maxBInFactorial,
      maxANPP = P(sim)$maxBInFactorial/30,
      year = 0
    )
  }

  ## check parameter consistency across modules
  paramCheckOtherMods(sim, "sppEquivCol", ifSetButDifferent = "error")

  if (!suppliedElsewhere("species", sim)) {
    message("generating dummy species data - run Biomass_borealDataPrep for table with real species attributes")
    speciesTable <- getSpeciesTable()
    sim$species <- prepSpeciesTable(speciesTable,
                                    sppEquiv = sim$sppEquiv,
                                    sppEquivCol = P(sim)$sppEquivCol)
  }

  if (!suppliedElsewhere("PSPmeasure_sppParams", sim) |
      !suppliedElsewhere("PSPplot_sppParams", sim) |
      !suppliedElsewhere("PSPgis_sppParams", sim)) {
    message("one or more PSP objects not supplied. Generating PSP data...")

    if ("dummy" %in% P(sim)$PSPdataTypes) {
      message("generating randomized PSP data")
      sim$PSPmeasure_sppParams <- Cache(prepInputs,
                                        targetFile = "randomizedPSPmeasure_sppParams.rds",
                                        archive = "randomized_LandR_speciesParameters_Inputs.zip",
                                        url =  extractURL("PSPmeasure_sppParams", sim),
                                        destinationPath = dPath,
                                        fun = "readRDS")

      sim$PSPplot_sppParams <- Cache(prepInputs,
                                     targetFile = "randomizedPSPplot_sppParams.rds",
                                     archive = "randomized_LandR_speciesParameters_Inputs.zip",
                                     url = extractURL("PSPplot_sppParams", sim),
                                     destinationPath = dPath,
                                     fun = "readRDS")

      sim$PSPgis_sppParams <- Cache(prepInputs,
                                    targetFile = "randomizedPSPgis_sppParams.rds",
                                    archive = "randomized_LandR_speciesParameters_Inputs.zip",
                                    url = extractURL("PSPgis_sppParams", sim),
                                    overwrite = TRUE,
                                    destinationPath = dPath,
                                    fun = "readRDS")
    } else if (!any(P(sim)$PSPdataTypes %in% "none")) {
      if (!any(c("BC", "AB", "SK", "NFI", "ON", "all") %in% P(sim)$PSPdataTypes)) {
        stop("Please review P(sim)$dataTypes - incorrect value specified")
      }

      PSPmeasure_sppParams <- list()
      PSPplot_sppParams <- list()

      if ("BC" %in% P(sim)$PSPdataTypes | "all" %in% P(sim)$PSPdataTypes) {
        PSPbc <- prepInputsBCPSP(dPath = dPath)
        PSPbc <- dataPurification_BCPSP(treeDataRaw = PSPbc$treeDataRaw,
                                        plotHeaderDataRaw = PSPbc$plotHeaderDataRaw,
                                        damageAgentCodes = PSPbc$pspBCdamageAgentCodes,
                                        codesToExclude = NULL)
        PSPmeasure_sppParams[["BC"]] <- PSPbc$treeData
        PSPplot_sppParams[["BC"]] <- PSPbc$plotHeaderData
      }

      if ("AB" %in% P(sim)$PSPdataTypes | "all" %in% P(sim)$PSPdataTypes) {
        PSPab <- prepInputsAlbertaPSP(dPath = dPath)
        PSPab <- dataPurification_ABPSP(treeMeasure = PSPab$pspABtreeMeasure,
                                        plotMeasure = PSPab$pspABplotMeasure,
                                        tree = PSPab$pspABtree,
                                        plot = PSPab$pspABplot,
                                        codesToExclude = NULL)
        #TODO: confirm if they really didn't record species on 11K trees
        PSPmeasure_sppParams[["AB"]] <- PSPab$treeData
        PSPplot_sppParams[["AB"]] <- PSPab$plotHeaderData
      }

      if ("SK" %in% P(sim)$PSPdataTypes | "all" %in% P(sim)$PSPdataTypes) {
        PSPsk <- prepInputsSaskatchwanPSP(dPath = dPath)
        PSPsk <- dataPurification_SKPSP(SADataRaw = PSPsk$SADataRaw,
                                        plotHeaderRaw = PSPsk$plotHeaderRaw,
                                        measureHeaderRaw = PSPsk$measureHeaderRaw,
                                        treeDataRaw = PSPsk$treeDataRaw)
        PSPmeasure_sppParams[["SK"]] <- PSPsk$treeData
        PSPplot_sppParams[["SK"]] <- PSPsk$plotHeaderData

        TSPsk <- prepInputsSaskatchwanTSP(dPath = dPath)
        TSPsk <- dataPurification_SKTSP_Mistik(compiledPlotData = TSPsk$compiledPlotData,
                                               compiledTreeData = TSPsk$compiledTreeData)
        PSPmeasure_sppParams[["SKtsp"]] <- TSPsk$treeData
        PSPplot_sppParams[["SKtsp"]] <- TSPsk$plotHeaderData
      }

      if ("ON" %in% P(sim)$PSPdataTypes | "all" %in% P(sim)$PSPdataTypes) {
        PSPon <- prepInputsOntarioPSP(dPath = dPath)
        #sppEquiv should not be subset to species of interest the way LandR requires
        #the latin is used to translate species into common names for the biomass equations
        sppEquivForON <- LandR::sppEquivalencies_CA #make sure this is fresh from the package
        PSPon <- dataPurification_ONPSP(PSPon, sppEquiv = sppEquivForON)
        PSPmeasure_sppParams[["ON"]] <- PSPon$treeData
        PSPplot_sppParams[["ON"]] <- PSPon$plotHeaderData
      }

      if ("NFI" %in% P(sim)$PSPdataTypes | "all" %in% P(sim)$PSPdataTypes) {

        PSPnfi <- prepInputsNFIPSP(dPath = dPath)
        PSPnfi <- dataPurification_NFIPSP(lgptreeRaw = PSPnfi$pspTreeMeasure,
                                          lgpHeaderRaw = PSPnfi$pspHeader,
                                          approxLocation = PSPnfi$pspLocation,
                                          treeDamage = PSPnfi$pspTreeDamage,
                                          codesToExclude = NULL)
        PSPmeasure_sppParams[["NFI"]] <- PSPnfi$treeData
        PSPplot_sppParams[["NFI"]] <- PSPnfi$plotHeaderData
      }

      PSPmeasure_sppParams <- rbindlist(PSPmeasure_sppParams, fill = TRUE)
      PSPplot_sppParams <- rbindlist(PSPplot_sppParams, fill = TRUE)

      PSPgis_sppParams <- geoCleanPSP(Locations = PSPplot_sppParams)

      #clean up
      toRemove <- c("Zone", "Datum", "Easting", "Northing", "Latitude", "Longitude")
      toRemove <- toRemove[toRemove %in% colnames(PSPplot_sppParams)]
      set(PSPplot_sppParams, NULL, toRemove, NULL)

      #keep only plots with valid coordinates
      PSPmeasure_sppParams <- PSPmeasure_sppParams[OrigPlotID1 %in% PSPgis_sppParams$OrigPlotID1,]
      PSPplot_sppParams <- PSPplot_sppParams[OrigPlotID1 %in% PSPgis_sppParams$OrigPlotID1,]
      sim$PSPmeasure_sppParams <- PSPmeasure_sppParams
      sim$PSPplot_sppParams <- PSPplot_sppParams
      sim$PSPgis_sppParams <- PSPgis_sppParams
    }
  }

  return(invisible(sim))
}
