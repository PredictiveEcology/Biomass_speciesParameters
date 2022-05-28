# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "Biomass_speciesParameters",
  description = "For estimating LANDIS-II species traits from PSP-derived growth curves",
  keywords = NA, # c("insert key words here"),
  authors = c(person(c("Ian"), "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = c("aut", "cre")),
              person(c("Eliot"), "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut"))),
  childModules = character(0),
  version = list(Biomass_speciesParameters = "1.0.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_speciesParameters.Rmd"),
  reqdPkgs = list("mgcv", "nlme", "fpCompare", "crayon", "data.table", "sf", "magrittr",
                  "PredictiveEcology/LandR@development (>= 1.0.5)",
                  "PredictiveEcology/pemisc@development (>= 0.0.3.9002)",
                  "PredictiveEcology/SpaDES.core@development (>= 1.0.9.9004)",
                  "ianmseddy/PSPclean", "robustbase", "gridExtra", "ggplot2", "purrr"),
  parameters = rbind(
    defineParameter(".plots", "character", "screen", NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "character", c(".inputObjects", "init"), NA, NA,
                    desc = paste("Should this entire module be run with caching activated?",
                                 "This is generally intended for data-type modules, where stochasticity and time are not relevant")),
    defineParameter("speciesFittingApproach", "character", "pairwise", NA, NA,
                    desc =  paste("Either 'all', 'pairwise', 'focal' or 'single', indicating whether to pool ",
                                  "all species into one fit, do pairwise species (for multiple cohort situations), do",
                                  "pairwise species, but using a focal species approach where all other species are ",
                                  "pooled into 'other'",
                                  " or do one species at a time. If 'all', all species will have identical species-level traits")),
    defineParameter("biomassModel", "character", "Lambert2005", NA, NA,
                    desc =  paste("The model used to calculate biomass from DBH. Can be either 'Lambert2005' or 'Ung2008'")),
    defineParameter("constrainGrowthCurve", "numeric", c(0, 1), 0, 1,
                    desc = paste("upper and lower bounds on range of potential growth curves when fitting traits. This module accepts a",
                                 "list of vectors, with names equal to `sppEquivCol`, so that traits are customizable")),
    defineParameter("constrainMortalityShape", 'numeric', c(15, 25), 5, 25,
                    desc = paste('upper and lower bounds on mortality shape when fitting traits. low mortality curve needs to excessive',
                                 'cohorts with very little biomass as longevity is approached, adding computation strain.',
                                 'alternatively accepts a list of vectors, with names equal to `sppEquivCol`')),
    defineParameter("constrainMaxANPP", 'numeric', c(2.0, 5.0), 1, 10,
                    desc = paste("upper and lower bounds on `maxANPP` when fitting traits. cohorts are initiated with `B = maxANPP`",
                                 "which may be unreasonably high if `maxANPP` is also high. Both `maxANPP` and growthcurve params",
                                 "control when `maxB` is reached. High `maxANPP` results in earlier peaks",
                                 'alternatively accepts a list of vectors, with names equal to `sppEquivCol`')),
    defineParameter("GAMMiterations", "numeric", 8, 1, NA,
                    desc = paste("number of iterations for GAMMs. This module accepts a",
                                 "list of vectors, with names equal to `sppEquivCol`, so that GAMMS are customizable")),
    defineParameter("GAMMknots", "numeric", 3, NA, NA,
                    desc = paste("the number of knots to use in the GAMM. Either 3 or 4 is recommended. This module accepts a",
                                 "list of vectors, with names equal to `sppEquivCol`, so that GAMMS are customizable")),
    defineParameter("minimumPlotsPerGamm", "numeric", 50, 10, NA, desc = paste("minimum number of PSP plots before building GAMM")),
    defineParameter("minDBH", "integer", 0L, 0L, NA,
                    desc = paste("minimum diameter at breast height (DBH) in cm used to filter PSP data.",
                                 "Defaults to 0 cm, i.e. all tree measurements are used.")),
    defineParameter("PSPdataTypes", "character", "all", NA, NA,
                    desc = paste("Which PSP datasets to source, defaulting to all. Other available options include",
                                 "'BC', 'AB', 'SK', 'NFI', and 'dummy'. 'dummy' should be used for unauthorized users.")),
    defineParameter("PSPperiod", "numeric", c(1920, 2019), NA, NA,
                    desc = paste("The years by which to subset sample plot data, if desired. Must be a vector of length 2")),
    defineParameter("quantileAgeSubset", "numeric", 95, 1, 100,
                    desc = paste("quantile by which to subset PSP data. As older stands are sparsely represented, the oldest measurements",
                                 "become vastly more influential. This parameter accepts both a single value and a list of vectors",
                                 "named by `sppEquivCol.` The PSP stand ages are found in `sim$speciesGAMMs$SPECIES$originalData`, where",
                                 "SPECIES is the species ID")),
    defineParameter("sppEquivCol", "character", 'default', NA, NA,
                    paste("The column in `sim$sppEquiv` data.table to group species by. This parameter should share the same",
                          "name as in *Biomass_borealDataPrep*. PSPs are aggregated by names in the PSP column and traits estimated",
                          "for the corresponding names in the `sppEquivCol`")),
    defineParameter("useHeight", "logical", FALSE, NA, NA,
                    desc = paste("Should height be used to calculate biomass (in addition to DBH).
                    Advise against including height unless you are certain it is present in every PSP"))
  ),
  inputObjects = bindrows(
    expectsInput(objectName  = "factorialSpeciesTable", objectClass = "data.table",
                 desc = paste("table with species traits for matching to `reducedFactorialCohortData`"),
                 sourceURL = "https://drive.google.com/open?id=1q0ou0CBzD9GqGSparpHqf318IWK6ycty"),
    expectsInput(objectName = "reducedFactorialCohortData", objectClass = "data.table",
                 desc = paste("results of factorial species trait simulation. This can be found by running",
                              "`SpeciesFactorial.R` but requires a specific commit of *Biomass_core*"),
                 sourceURL = "https://drive.google.com/open?id=1h8StXE0vm8xyDycRomCkwIaL7wfh5Irj"),
    expectsInput(objectName = "PSPmeasure", objectClass = "data.table",
                 desc = paste("merged PSP and TSP individual tree measurements. Must include the following columns:",
                              "MeasureID, OrigPlotID1, MeasureYear, TreeNumber, Species, DBH and newSpeciesName",
                              "the latter corresponding to species names in `LandR::sppEquivalencies_CA$PSP`.",
                              "Defaults to randomized PSP data stripped of real plotIDs"),
                 sourceURL = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing"),
    expectsInput(objectName = "PSPplot_sppParams", objectClass = "data.table",
                 desc = paste("merged PSP and TSP plot data. Defaults to randomized PSP data stripped of real plotIDs.",
                              "Must contain fields 'MeasureID', 'MeasureYear', 'OrigPlotID1', and 'baseSA'",
                              "the latter being stand age at year of first measurement"),
                 sourceURL = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing"),
    expectsInput(objectName = "PSPgis", objectClass = "sf",
                 desc = paste("Plot location sf object. Defaults to PSP data stripped of real plotIDs/location.",
                              "Must include field 'OrigPlotID1' for joining to PSPplot object"),
                 sourceURL = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing"),
    expectsInput(objectName = "species", objectClass = "data.table",
                 desc = paste("a table of invariant species traits with the following trait colums:",
                              "'species', 'Area', 'longevity', 'sexualmature', 'shadetolerance',",
                              "'firetolerance', 'seeddistance_eff', 'seeddistance_max', 'resproutprob',",
                              "'mortalityshape', 'growthcurve', 'resproutage_min', 'resproutage_max',",
                              "'postfireregen', 'wooddecayrate', 'leaflongevity' 'leafLignin',",
                              "'hardsoft'. Only 'growthcurve' and 'mortalityshape' are used in this module.",
                              "Default is from Dominic Cyr and Yan Boulanger's project"),
                 sourceURL = "https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv"),
    expectsInput(objectName = "speciesEcoregion", objectClass = "data.table",
                 desc = paste("table of spatially-varying species traits (`maxB`, `maxANPP`,",
                              "`establishprob`), defined by species and `ecoregionGroup`)",
                              "Defaults to a dummy table based on dummy data os biomass, age, ecoregion and land cover class")),
    expectsInput(objectName = "sppEquiv", objectClass = "data.table",
                 desc = "table of species equivalencies. See `?LandR::sppEquivalencies_CA`."),
    expectsInput(objectName = "studyAreaANPP", objectClass = "SpatialPolygonsDataFrame",
                 desc = "study area used to crop PSP data before building growth curves")
  ),
  outputObjects = bindrows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput("species", "data.table",
                  desc = "a table that has species traits such as longevity..."),
    createsOutput(objectName = "speciesEcoregion", "data.table",
                  desc = paste("table of spatially-varying species traits (`maxB`, `maxANPP`,",
                              "`establishprob`), defined by species and `ecoregionGroup`)",
                               "Defaults to a dummy table based on dummy data os biomass, age, ecoregion and land cover class")),
    createsOutput(objectName = 'speciesGAMMs', objectClass = 'list',
                  desc = paste('a list of mixed-effect general additive models (gamm) for each tree species',
                               'modeling biomass as a function of age'))
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
      # sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "Biomass_speciesParameters", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "Biomass_speciesParameters", "save")
    },
    plot = {
      # lapply(names(sim$speciesGAMMs), FUN = function(spp, GAMMs = sim$speciesGAMMs){
      #   if (!class(GAMMs[[spp]]) == "character"){
      #     gam <- GAMMs[[spp]][["gam"]]
      #     jpeg(filename = file.path(outputPath(sim), paste0(spp, "_gamm.jpg")))
      #     plot(gam, xlab = "stand age", ylab = "biomass", main = spp)
      #     dev.off()
      #   }
      # })
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
    stop("Please supply 'sppEquivCol' in parameters of Biomass_speciesParameters.")
  }

  paramCheckOtherMods(sim, "maxBInFactorial")
  paramCheckOtherMods(sim, paramToCheck = "sppEquivCol",
                      ifSetButDifferent = "error")

  #find the max biomass achieved by each species when growing with no competition
  tempMaxB <- sim$cohortDataFactorial[age == 1, .N, .(pixelGroup)]
  #take the pixelGroups with only 1 species at start of factorial
  tempMaxB <- tempMaxB[N == 1,]
  tempMaxB <- sim$cohortDataFactorial[pixelGroup %in% tempMaxB$pixelGroup,
                                      .(inflationFactor = P(sim)$maxBInFactorial/max(B)),
                                      , .(pixelGroup, speciesCode)]

  # in some cases the speciesTableFactorial doesn't have "species" column; just "speciesCode"
  if (!("species" %in% colnames(sim$speciesTableFactorial))) {
    setnames(sim$speciesTableFactorial, old = "speciesCode", new = "species")
  }
  tempMaxB <- sim$speciesTableFactorial[tempMaxB, on = c("species" = "speciesCode", "pixelGroup")]
  #pair-wise species will be matched with traits, as the species code won't match
  tempMaxB <- tempMaxB[, .(species, longevity, growthcurve, mortalityshape, mANPPproportion, inflationFactor)]

  #prepare PSPdata
  speciesGAMMs <- Cache(makePSPgamms,
                        studyAreaANPP = sim$studyAreaANPP,
                        PSPperiod = P(sim)$PSPperiod,
                        PSPgis = sim$PSPgis_sppParams,
                        PSPmeasure = sim$PSPmeasure_sppParams,
                        PSPplot = sim$PSPplot_sppParams,
                        useHeight = P(sim)$useHeight,
                        biomassModel = P(sim)$biomassModel,
                        speciesCol = P(sim)$sppEquivCol,
                        sppEquiv = sim$sppEquiv,
                        NoOfIterations = P(sim)$GAMMiterations,
                        knots = P(sim)$GAMMknots,
                        minimumSampleSize = P(sim)$minimumPlotsPerGamm,
                        quantileAgeSubset = P(sim)$quantileAgeSubset,
                        minDBH = P(sim)$minDBH,
                        speciesFittingApproach = P(sim)$speciesFittingApproach,
                        userTags = c(currentModule(sim), "makePSPgamms"))
  sim$speciesGAMMs <- speciesGAMMs

  classes <- lapply(sim$speciesGAMMs, FUN = 'class')

  noData <- vapply(sim$speciesGAMMs[classes == "character"], FUN = function(x) {
    x == "insufficient data"
  }, FUN.VALUE = logical(1))

  if (any(noData)) {
    message("The following species did not have sufficient data for model estimation: ")
    print(names(noData))
  }
  speciesWithNewlyEstimated <- unique(unlist(strsplit(names(sim$speciesGAMMs), "__")))
  speciesWithoutNewlyEstimated <- setdiff(sim$sppEquiv[[Par$sppEquivCol]], speciesWithNewlyEstimated)
  if (length(speciesWithoutNewlyEstimated))
    message(crayon::yellow(paste(speciesWithoutNewlyEstimated,
                                 collapse = ", "),
                           "have insufficient data to estimate species parameters; using original user supplied"))
  modifiedSpeciesTables <- modifySpeciesTable(gamms = sim$speciesGAMMs,
                                              speciesTable = sim$species,
                                              factorialTraits = setDT(sim$speciesTableFactorial), # setDT to deal with reload from Cache (no effect otherwise)
                                              factorialBiomass = setDT(sim$cohortDataFactorial), # setDT to deal with reload from Cache (no effect otherwise)
                                              sppEquiv = sim$sppEquiv,
                                              approach = P(sim)$speciesFittingApproach,
                                              sppEquivCol = P(sim)$sppEquivCol,
                                              maxBInFactorial = P(sim)$maxBInFactorial,
                                              inflationFactorKey = tempMaxB,
                                              standAgesForFitting = P(sim)$standAgesForFitting,
                                              mortConstraints = P(sim)$constrainMortalityShape,
                                              growthConstraints = P(sim)$constrainGrowthCurve,
                                              mANPPconstraints = P(sim)$constrainMaxANPP)
  gg <- modifiedSpeciesTables$gg
  Plots(gg, usePlot = FALSE, fn = print, ggsaveArgs = list(width = 10, height = 7),
        filename = paste("Pairwise species fits ", gsub(":", "_", sim$._startClockTime)))
  sim$species <- modifiedSpeciesTables$best

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

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  if (!suppliedElsewhere("cohortDataFactorial", sim)) {
    sim$cohortDataFactorial <- prepInputs(targetFile = "reducedFactorialCD.Rdat",
                                          destinationPath = dPath,
                                          fun = "readRDS", overwrite = TRUE,
                                          url = extractURL('cohortDataFactorial', sim),
                                          useCache = TRUE, userTags = c(cacheTags, "reducedFactorial"))
  }

  if (!suppliedElsewhere("speciesTableFactorial", sim)) {
    sim$speciesTableFactorial <- prepInputs(targetFile = "speciesTableFactorial.Rdat",
                                            destinationPath = dPath,
                                            url = extractURL('speciesTableFactorial', sim),
                                            fun = "readRDS", overwrite = TRUE,
                                            useCache = TRUE, userTags = c(cacheTags, "factorialSpecies"))
  }

  if (!suppliedElsewhere("speciesEcoregion", sim)) {
    warning("generating dummy speciesEcoregion data - run Biomass_borealDataPrep for table with real speciesEcoregion attributes")
    sim$speciesEcoregion <- data.table(
      ecoregionGroup = "x",
      speciesCode = c("Abie_las", "Abie_bal", "Betu_pap", "Lari_lar", "Pice_eng",
                      "Pice_gla", "Pice_mar", "Pinu_ban",
                      "Pinu_con", "Pseu_men", "Popu_tre"),
      establishprob = 0.5,
      maxB = P(sim)$maxBInFactorial,
      maxANPP = P(sim)$maxBInFactorial/30,
      year = 0
    )
  }

  if (!suppliedElsewhere("species", sim)) {
    warning("generating dummy species data - run Biomass_borealDataPrep for table with real species attributes")
    sim$species <- data.table(
      species = c("Abie_las", "Abie_bal", "Betu_pap", "Lari_lar", "Pice_eng",
                  "Pice_gla", "Pice_mar", "Pinu_ban",
                  "Pinu_con", "Pseu_men", "Popu_tre"),
      longevity = c(300, 300, 150, 140, 450, 400, 250, 150, 325, 600, 200),
      mortalityshape = 15,
      growthcurve = 0
    )
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

  if (!suppliedElsewhere("PSPmeasure_sppParams", sim) |
      !suppliedElsewhere("PSPplot_sppParams", sim) |
      !suppliedElsewhere("PSPgis_sppParams", sim)) {
    message("one or more PSP objects not suppplied. Generating PSP data...")

    if ("dummy" %in% P(sim)$PSPdataTypes) {
      message("generating randomized PSP data")
      sim$PSPmeasure_sppParams <- Cache(prepInputs,
                                        targetFile = "randomizedPSPmeasure_sppParams.rds",
                                        archive = "randomized_LandR_speciesParameters_Inputs.zip",
                                        url =  extractURL('PSPmeasure_sppParams', sim),
                                        destinationPath = dPath,
                                        fun = "readRDS")

      sim$PSPplot_sppParams <- Cache(prepInputs,
                                     targetFile = "randomizedPSPplot_sppParams.rds",
                                     archive = "randomized_LandR_speciesParameters_Inputs.zip",
                                     url = extractURL('PSPplot_sppParams', sim),
                                     destinationPath = dPath,
                                     fun = "readRDS")

      sim$PSPgis_sppParams <- Cache(prepInputs,
                                    targetFile = "randomizedPSPgis_sppParams.rds",
                                    archive = "randomized_LandR_speciesParameters_Inputs.zip",
                                    url = extractURL('PSPgis_sppParams', sim),
                                    overwrite = TRUE,
                                    destinationPath = dPath,
                                    fun = "readRDS")
    } else {
      if (!any(c("BC", "AB", "SK", "NFI", "all") %in% P(sim)$PSPdataTypes)) {
        stop("Please review P(sim)$dataTypes - incorrect value specified")
      }

      PSPmeasure_sppParams <- list()
      PSPplot_sppParams <- list()

      if ("BC" %in% P(sim)$PSPdataTypes | "all" %in% P(sim)$PSPdataTypes) {
        PSPbc <- Cache(prepInputsBCPSP, dPath = dPath, userTags = c(cacheTags, "BCPSP"))
        PSPbc <- PSPclean::dataPurification_BCPSP(treeDataRaw = PSPbc$treeDataRaw,
                                                  plotHeaderDataRaw = PSPbc$plotHeaderDataRaw,
                                                  damageAgentCodes = PSPbc$pspBCdamageAgentCodes,
                                                  codesToExclude = NULL)
        PSPmeasure_sppParams[["BC"]] <- PSPbc$treeData
        PSPplot_sppParams[["BC"]] <- PSPbc$plotHeaderData
      }

      if ("AB" %in% P(sim)$PSPdataTypes | "all" %in% P(sim)$PSPdataTypes) {
        PSPab <- Cache(prepInputsAlbertaPSP, dPath = dPath, userTags = c(cacheTags, "ABPSP"))
        PSPab <- PSPclean::dataPurification_ABPSP(treeMeasure = PSPab$pspABtreeMeasure,
                                                  plotMeasure = PSPab$pspABplotMeasure,
                                                  tree = PSPab$pspABtree,
                                                  plot = PSPab$pspABplot,
                                                  codesToExclude = NULL)
        PSPmeasure_sppParams[["AB"]] <- PSPab$treeData
        PSPplot_sppParams[["AB"]] <- PSPab$plotHeaderData
      }

      if ("SK" %in% P(sim)$PSPdataTypes | "all" %in% P(sim)$PSPdataTypes) {
        PSPsk <- Cache(prepInputsSaskatchwanPSP, dPath = dPath, userTags = c(cacheTags, "SKPSP"))
        PSPsk <- PSPclean::dataPurification_SKPSP(SADataRaw = PSPsk$SADataRaw,
                                                  plotHeaderRaw = PSPsk$plotHeaderRaw,
                                                  measureHeaderRaw = PSPsk$measureHeaderRaw,
                                                  treeDataRaw = PSPsk$treeDataRaw)
        PSPmeasure_sppParams[["SK"]] <- PSPsk$treeData
        PSPplot_sppParams[["SK"]] <- PSPsk$plotHeaderData

        TSPsk <- Cache(prepInputsSaskatchwanTSP, dPath = dPath, userTags = c(cacheTags, "SKTSP"))
        TSPsk <- PSPclean::dataPurification_SKTSP_Mistik(compiledPlotData = TSPsk$compiledPlotData,
                                                         compiledTreeData = TSPsk$compiledTreeData)
        PSPmeasure_sppParams[["SKtsp"]] <- TSPsk$treeData
        PSPplot_sppParams[["SKtsp"]] <- TSPsk$plotHeaderData
      }

      if ("NFI" %in% P(sim)$PSPdataTypes | "all" %in% P(sim)$PSPdataTypes) {

        PSPnfi <- Cache(prepInputsNFIPSP, dPath = dPath, userTags = c(cacheTags, "NFIPSP"))
        PSPnfi <- PSPclean::dataPurification_NFIPSP(lgptreeRaw = PSPnfi$pspTreeMeasure,
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
### add additional events as needed by copy/pasting from above
