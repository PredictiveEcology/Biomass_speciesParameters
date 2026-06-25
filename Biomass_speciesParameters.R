defineModule(sim, list(
  name = "Biomass_speciesParameters",
  description = paste(
    "For estimating LANDIS-II species traits based on growth curves derived",
    "from Permanent Sample Plot (PSP) and Temporary Sample Plot (TSP) data"
  ),
  keywords = NA, # TODO: add keywords
  authors = c(
    person(c("Ian"), "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person(c("Eliot"), "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut")),
    person(c("Ceres"), "Barros", email = "ceres.barros@ubc.ca", role = c("ctb"))
  ),
  childModules = character(0),
  version = list(Biomass_speciesParameters = "3.0.1"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "Biomass_speciesParameters.Rmd"),
  loadOrder = list(after = c("Biomass_speciesFactorial", "Biomass_borealDataPrep"),
                   before = c("Biomass_core")),
  reqdPkgs = list(
    "arrow", "cli", "data.table", "fpCompare", "fs", "ggplot2", "gridExtra",
    "mgcv", "nlme", "purrr", "robustbase", "sf",
    "reproducible (>= 2.1.0)",
    "SpaDES.core (>= 2.1.4)",
    "PredictiveEcology/LandR (>= 1.1.0.9077)",
    "PredictiveEcology/pemisc@development (>= 0.0.3.9002)",
    "ianmseddy/PSPclean@development(>= 1.0.0.9001)"
  ),
  parameters = rbind(
    defineParameter("biomassModel", "character", "Lambert2005", NA, NA,
                    desc =  paste("The model used to calculate biomass from DBH. Can be either 'Lambert2005' or 'Ung2008'.")),
    defineParameter("landis", "logical", FALSE, NA, NA,
                    desc = paste("If `TRUE`, run in 'LANDIS mode': expose the fitted LANDIS-version growth curves",
                                 "(the scaled non-linear biomass-over-age curves shown in the LandR-vs-non-linear plot)",
                                 "per species as the output object `speciesGrowthCurvesLandis`, for use as inputs to",
                                 "LANDIS-II Biomass Succession. The per-species growth-curve parameters",
                                 "(`growthcurve`, `mortalityshape`, etc.) are written to `species` regardless.",
                                 "Default `FALSE` preserves the standard behaviour.")),
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
                                 "'BC', 'AB', 'SK', 'ON', 'NB', 'NFI', and 'dummy'.",
                                 "'dummy' should be used for unauthorized users.")),
    defineParameter("PSPperiod", "numeric", c(1920, 2019), NA, NA,
                    desc = paste("The years by which to subset sample plot data, if desired. Must be a vector of length 2")),
    defineParameter("quantileAgeSubset", "numeric", 99, 1, 100,
                    desc = paste("Quantile by which to subset PSP data. As older stands are sparsely represented",
                                 "the oldest measurements become vastly more influential. This parameter accepts",
                                 "both a single value and a list of vectors, named according to `sppEquivCol`.")),
    defineParameter("speciesFittingApproach", "character", "focal", NA, NA,
                    desc =  paste(
                      "Either 'all', 'pairwise', 'focal' or 'single', indicating whether to pool",
                      "all species into one fit, do pairwise species (for multiple cohort situations)",
                      "do pairwise species, but using a focal species approach where all other species",
                      "are pooled into 'other' or do one species at a time.",
                      "If 'all', all species will have identical species-level traits."
                    )),
    defineParameter("sppEquivCol", "character", "LandR", NA, NA,
                    paste("The column in `sim$sppEquiv` data.table that defines individual species.",
                          "The names should match those in the species table.")),
    defineParameter("standAgesForFitting", "integer", c(21L, 91L), NA, NA,
                    desc = paste("The minimum and maximum ages of the biomass-by-age curves used in fitting.",
                                 "It is generally recommended to keep this param under 200, given the low data",
                                 "availability of stands aged 200+, with some exceptions.",
                                 "For a closed interval, end with a 1, e.g. `c(31, 101)`.")),
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
    defineParameter(".useCache", "character", c(".inputObjects"), NA, NA,
                    desc = paste("Should this entire module be run with caching activated?",
                                 "This is generally intended for data-type modules,",
                                 "where stochasticity and time are not relevant")),
    defineParameter(".useParallel", "integer", 2L, NA, NA,
                    desc = paste("maximum number of threads/workers to use for data.table operations;",
                                 "passed to `data.table::setDTthreads` and should be <= 4."))
  ),
  inputObjects = bindrows(
    expectsInput("cohortDataFactorial_path", "fs_path",
                 desc = paste(
                   "Path where the `cohortDataFactorial` object is saved as an `arrow` dataset.",
                   "A large `cohortData` table (**sensu** `Biomass_core`) with columns `age`, `B`,",
                   "and `speciesCode` that joins with `speciesTableFactorial`.",
                   "See `PredictiveEcology/Biomass_factorial` for further information."
                  ),
                 sourceURL = "https://drive.google.com/file/d/1NH7OpAnWtLyO8JVnhwdMJakOyapBnuBH/"),
    expectsInput("PSPmeasure_sppParams", "data.table",
                 desc = paste("Merged PSP and TSP individual tree measurements. Must include the following columns:",
                              "`MeasureID`, `OrigPlotID1`, `MeasureYear`, `TreeNumber`, `Species`, `DBH` and `PSP`,",
                              "where `Species` corresponds to species names in `LandR::sppEquivalencies_CA$Latin_full`.",
                              "Defaults to randomized PSP data stripped of real `plotID`s"),
                 sourceURL = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing"),
    expectsInput("PSPplot_sppParams", "data.table",
                 desc = paste("Merged PSP and TSP plot data. Defaults to randomized PSP data stripped of real `plotID`s.",
                              "Must contain columns `MeasureID`, `MeasureYear`, `OrigPlotID1`, and `baseSA`,",
                              "the latter being stand age at year of first measurement"),
                 sourceURL = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing"),
    expectsInput("PSPgis_sppParams", "sf",
                 desc = paste("Plot location `sf` object. Defaults to PSP data stripped of real `plotID`s/location.",
                              "Must include field `OrigPlotID1` for joining to `PSPplot` object"),
                 sourceURL = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing"),
    expectsInput("species", "data.table",
                 desc = paste("A table of invariant species traits with the following trait colums:",
                              "'species', 'Area', 'longevity', 'sexualmature', 'shadetolerance',",
                              "'firetolerance', 'seeddistance_eff', 'seeddistance_max', 'resproutprob',",
                              "'mortalityshape', 'growthcurve', 'resproutage_min', 'resproutage_max',",
                              "'postfireregen', 'wooddecayrate', 'leaflongevity' 'leafLignin', and 'hardsoft'.",
                              "Only 'growthcurve', 'hardsoft',  and 'mortalityshape' are used in this module.",
                              "Default is from Dominic Cyr and Yan Boulanger's applications of LANDIS-II"),
                 sourceURL = "https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv"),
    expectsInput("speciesEcoregion", "data.table",
                 desc = paste("Table of spatially-varying species traits (`maxB`, `maxANPP`,",
                              "`establishprob`), defined by species and `ecoregionGroup`).",
                              "Defaults to a dummy table based on dummy data of biomass, age,",
                              "ecoregion and land cover class.")),
    expectsInput("speciesTableFactorial_path", "fs_path",
                 desc = paste(
                   "Path where the `speciesTableFactorial` object is saved as an `arrow` dataset.",
                   "A large species table (**sensu** `Biomass_core`) with all columns used by",
                   "Biomass_core, e.g., `longevity`, `growthcurve`, `mortalityshape`, etc., when",
                   "it was used to generate `cohortDataFactorial`.",
                   "See `PredictiveEcology/Biomass_factorial` for futher information."
                 ),
                 sourceURL = "https://drive.google.com/file/d/1NH7OpAnWtLyO8JVnhwdMJakOyapBnuBH/"),
    expectsInput("sppEquiv", "data.table",
                 desc = paste("Table of species equivalencies - see `?LandR::sppEquivalencies_CA`.", 
                              "Traits will be estimated for each unique entry in the `sppEquivCol` column.")),
    expectsInput("sppEquivLong", "data.table", 
                 desc = paste("The full table of species equivalencies - see `?LandR::sppEquivalencies_CA`.",
                              "Biomass will be estimated for each species based on the `sp_Biomass_eq' column,",
                              "which uses `pemisc::biomassCalculation` to derive AGB from DBH and height",
                              "(based on the equations from https://doi.org/10.1139/x05-112). The full table",
                              "is used to improve stand biomass estimates even if some species are not of interest.")),
    expectsInput("studyAreaANPP", "sf",                                         # nolint: in_no_default
                 desc = paste("Optional study area used to crop PSP data before building growth curves.",
                              "If supplied, an ecoregion-scale object is recommended, at a minimum."))
  ),
  outputObjects = bindrows(
    createsOutput("species", "data.table",
                  desc = paste("The updated invariant species traits table (see above).")),
    createsOutput("speciesEcoregion", "data.table",
                  desc = paste("The updated spatially-varying species traits table",
                               "(see description for this object in inputs)")),
    createsOutput("speciesGrowthCurves", "list",
                  desc = paste("list containing each species' non-linear model,",
                               "model data, and the unfiltered PSP data")),
    createsOutput("speciesGrowthCurvesLandis", "data.table",
                  desc = paste("Only when `P(sim)$landis` is `TRUE`. The fitted LANDIS-version growth curves",
                               "(`BscaledNonLinear`) by `species` and `standAge` (the scaled non-linear curves shown",
                               "in the `LandR_VS_NLM_growthCurves` plot), for use as LANDIS-II Biomass Succession inputs."))
  )
))

## event types
#   - type `init` is required for initialization

doEvent.Biomass_speciesParameters = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      ## build growth curves if applicable
      sim <- Init(sim)

      ## update tables
      sim <- updateSpeciesTables(sim)
    },
    warning(
      paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
            "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = "")
    )
  )
  return(invisible(sim))
}

## event functions
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
Init <- function(sim) {
  origDTthreads <- data.table::getDTthreads()
  if (getDTthreads() > P(sim)$.useParallel) {
    data.table::setDTthreads(P(sim)$.useParallel)
  }
  on.exit(data.table::setDTthreads(origDTthreads))
  ## load factorial tables -------------------------------------------------------------------------
  fmt <- "feather" ## faster for small-med data compared to parquet

  if (tools::file_ext(sim$cohortDataFactorial_path) == "rds") {
    ## using rds default from .inputObjects
    mod$cohortDataFactorial <- prepInputs(
      targetFile = basename(sim$cohortDataFactorial_path),
      destinationPath = inputPath(sim),
      fun = "readRDS",
      overwrite = TRUE, useCache = FALSE, # don't internal cache as it is a waste of time
      url = extractURL("cohortDataFactorial_path", sim)
    ) |> Cache(.functionName = "prepInputs_cohortDataFactorial")
  } else {
    ## connect to arrow dataset
    ## TODO: consider adding try-catch and update Biomass_speciesFactorial if fails occur
    mod$cohortDataFactorial <- arrow::open_dataset(sim$cohortDataFactorial_path, format = fmt)
  }

  if (tools::file_ext(sim$speciesTableFactorial_path) == "rds") {
    ## using rds default from .inputObjects
    mod$speciesTableFactorial <- prepInputs(
      targetFile = basename(sim$speciesTableFactorial_path),
      destinationPath = inputPath(sim),
      url = extractURL("speciesTableFactorial_path", sim),
      fun = "readRDS", useCache = FALSE, # don't internal cache as it is a waste of time
      overwrite = TRUE
    ) |> Cache(.functionName = "prepInputs_speciesTableFactorial")
  } else {
    ## connect to arrow dataset
    ## TODO: consider adding try-catch and update Biomass_speciesFactorial if fails occur
    mod$speciesTableFactorial <- arrow::open_dataset(sim$speciesTableFactorial_path, format = fmt)
  }

  ## if no PSP data supplied, simList returned unchanged

  if (all(P(sim)$PSPdataTypes != "none")) {
    if (is.na(P(sim)$sppEquivCol)) {
      stop("Please supply 'sppEquivCol' in parameters of Biomass_speciesParameters.")
    }

    paramCheckOtherMods(sim, "maxBInFactorial")
    paramCheckOtherMods(sim, paramToCheck = "sppEquivCol", ifSetButDifferent = "error")

    ## 2025-10: removed deprecated package disk.frame, which no longer supports data.table syntax;
    ##          switched to using arrow, which uses dplyr syntax;
    ##          (retained data.table versions for comparison only).

    ## find the max biomass achieved by each species when growing with no competition
    # tempMaxB <- mod$cohortDataFactorial[age == 1, .N, .(pixelGroup)]
    tempMaxB <- mod$cohortDataFactorial |>
      dplyr::filter(age == 1) |>
      dplyr::group_by(pixelGroup) |>
      dplyr::summarise(N = n()) |>
      dplyr::collect()

    ## take the pixelGroups with only 1 species at start of factorial.
    # tempMaxB <- tempMaxB[N == 1, ]
    tempMaxB <- tempMaxB |> dplyr::filter(N == 1)

    # tempMaxB <- mod$cohortDataFactorial[pixelGroup %in% tempMaxB$pixelGroup,
    #                                     .(inflationFactor = P(sim)$maxBInFactorial/max(B)),
    #                                     , .(pixelGroup, speciesCode)]
    tempMaxB <- mod$cohortDataFactorial |>
      dplyr::filter(pixelGroup %in% tempMaxB$pixelGroup) |>
      dplyr::group_by(pixelGroup, speciesCode) |>
      dplyr::summarise(inflationFactor = P(sim)$maxBInFactorial / max(B)) |>
      dplyr::collect()

    ## speciesTableFactorial sometimes doesn't have 'species' column (only 'speciesCode');
    ## TODO this is a work around -- the speciesTableFactorial should be stable
    if (!("species" %in% names(mod$speciesTableFactorial))) {
      # setnames(speciesTableFactorial, old = "speciesCode", new = "species")
      mod$speciesTableFactorial <- mod$speciesTableFactorial |>
        dplyr::rename(species = speciesCode) |>
        dplyr::collect()
    }

    ## NOTE: disk.frame doesn't support right_join, so need to rework as left_join below
    # tempMaxB <- mod$speciesTableFactorial[tempMaxB, on = c("species" = "speciesCode", "pixelGroup")]
    tempMaxB <- arrow_table(tempMaxB) |>
      dplyr::rename(species = speciesCode) |>
      dplyr::left_join(mod$speciesTableFactorial, by = c("species", "pixelGroup")) |>
      dplyr::collect()

    ## pair-wise species will be matched with traits, as the species code won't match
    # tempMaxB <- tempMaxB[, .(species, longevity, growthcurve, mortalityshape, mANPPproportion, inflationFactor)]
    tempMaxB <- tempMaxB |>
      dplyr::select(species, longevity, growthcurve, mortalityshape, mANPPproportion, inflationFactor) |>
      setDT()

    ## bring tables to RAM for use below
    mod$cohortDataFactorial <- dplyr::collect(mod$cohortDataFactorial) |> setDT()
    mod$speciesTableFactorial <- dplyr::collect(mod$speciesTableFactorial) |> setDT()

    gc()

    #TODO:  #change PSP to be a dedicated biomass equation column in sppEquiv
    biomassKey <- unique(sim$sppEquivLong[, .(PSP, Latin_full)])
    setnames(biomassKey, old = "PSP", new = "SpBiomassEq")

    message("preparing PSPs for growth curves")
    psp <- prepPSPaNPP(studyAreaANPP = sim$studyAreaANPP, PSPperiod = P(sim)$PSPperiod,
                       PSPgis =  sim$PSPgis_sppParams, PSPmeasure = sim$PSPmeasure_sppParams, 
                       PSPplot = sim$PSPplot_sppParams, useHeight = P(sim)$useHeight, 
                       biomassModel = P(sim)$biomassModel, minDBH = P(sim)$minDBH, 
                       sppEquivLong = biomassKey) |>
      Cache(userTags = c(currentModule(sim), "prepPSPaNPP"))

    message("building growth curves") # this cache call takes several minutes to process..
    sim$speciesGrowthCurves <-  buildGrowthCurves(PSPdata = psp, 
                                                  speciesCol =  P(sim)$sppEquivCol, 
                                                  sppEquiv = sim$sppEquiv,
                                                  minimumSampleSize = P(sim)$minimumPlots,
                                                  quantileAgeSubset = P(sim)$quantileAgeSubset,
                                                  speciesFittingApproach = P(sim)$speciesFittingApproach) |>
    Cache(userTags = c(currentModule(sim), "buildGrowthCurves"))

    classes <- lapply(sim$speciesGrowthCurves, FUN = "class")
    
    noDataSpp <- vapply(sim$speciesGrowthCurves[classes == "character"], FUN = function(x) {
      x == "insufficient data"
    }, FUN.VALUE = logical(1))
    
    if (any(noDataSpp)) {
      message(cli::col_yellow(
        "Insufficient data to estimate species parameters for ",
        paste(names(noDataSpp), collapse = ", "),
        " - will keep original user-supplied parameters"
      ))
    }

    message("Digesting growth curve library and factorial...")
    # Digesting of the growth curves has a stats object which changes every time;
    #   just take the summary of the stats object
    speciesGrowthCurves <-
      sim$speciesGrowthCurves[!names(sim$speciesGrowthCurves) %in% names(noDataSpp)]
    spgForDigest <- Map(sgc = speciesGrowthCurves, function(sgc) {
      sgc$NonLinearModel <- summary(sgc$NonLinearModel[[1]]) # has no "function"
      sgc
    }
      )
    toDigest <- list(
      speciesGrowthCurves = spgForDigest, # this has functions in it, so needs to be dealt with
      speciesTableFactorial = setDT(mod$speciesTableFactorial),
      cohortDataFactorial = setDT(mod$cohortDataFactorial)
    )
    cacheExtra <- .robustDigest(toDigest)
    message("Done!")

    modifiedSpeciesTables <- modifySpeciesTable(
      GCs = speciesGrowthCurves,
      speciesTable = sim$species,
      factorialTraits = setDT(mod$speciesTableFactorial),
      ## setDT to deal with reload from Cache (no effect otherwise)
      factorialBiomass = setDT(mod$cohortDataFactorial),
      ## setDT to deal with reload from Cache (no effect otherwise)
      sppEquiv = sim$sppEquiv,
      sppEquivCol = P(sim)$sppEquivCol,
      inflationFactorKey = tempMaxB,
      standAgesForFitting = P(sim)$standAgesForFitting,
      approach = P(sim)$speciesFittingApproach,
      maxBInFactorial = P(sim)$maxBInFactorial
    ) |> Cache(
      omitArgs = c("GCs", "factorialTraits", "factorialBiomass"),
      .cacheExtra = cacheExtra,
      userTags = c(currentModule(sim), "modifiedSpeciesTables")
    )

    gg <- modifiedSpeciesTables$gg
    suppressWarnings(Plots(gg, usePlot = FALSE, fn = print, ggsaveArgs = list(width = 10, height = 7),
          filename = paste("LandR_VS_NLM_growthCurves")))
    sim$species <- modifiedSpeciesTables$best
    if (isTRUE(P(sim)$landis)) {
      ## LANDIS mode: expose the fitted LANDIS-version growth curves (BscaledNonLinear by species and
      ## standAge) for use as LANDIS-II Biomass Succession inputs.
      sim$speciesGrowthCurvesLandis <- modifiedSpeciesTables$landisCurves
    }
  } else {
    message("P(sim)$PSPdataTypes is 'none' -- bypassing species traits estimation from PSP data.")
  }
  return(sim)
}

updateSpeciesTables <- function(sim) {
  modifiedTables <- modifySpeciesAndSpeciesEcoregionTable(
    speciesEcoregion = sim$speciesEcoregion,
    speciesTable = sim$species
  )
  sim$speciesEcoregion <- modifiedTables$newSpeciesEcoregion
  sim$species <- modifiedTables$newSpeciesTable
  return(sim)
}

### template for save events
Save <- function(sim) {
  sim <- saveFiles(sim)
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  origDTthreads <- data.table::getDTthreads()
  if (origDTthreads > 4) {
    data.table::setDTthreads(4)
  }
  on.exit(data.table::setDTthreads(origDTthreads))

  cacheTags <- c(currentModule(sim), "OtherFunction:.inputObjects")
  dPath <- asPath(inputPath(sim), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  if (!suppliedElsewhere("cohortDataFactorial_path", sim)) {
    sim$cohortDataFactorial_path <- file.path(dPath, "cohortDataFactorial_medium.rds")
  }

  if (!suppliedElsewhere("speciesTableFactorial_path", sim)) {
    sim$speciesTableFactorial_path <- file.path(dPath, "speciesTableFactorial_medium.rds")
  }

  if (!suppliedElsewhere("sppEquiv", sim)) {
    ## pass a default sppEquivalencies_CA for common species in western Canada
    sppEquiv <- LandR::sppEquivalencies_CA
    ## sppEquiv should be Boreal to match default P(sim)$sppEquivCol
    sim$sppEquiv <- sppEquiv[Boreal %in% c("Pice_Mar", "Pice_Gla", "Pinu_Con",
                                           "Pinu_Ban", "Popu_Tre", "Lari_Lar",
                                           "Betu_Pap", "Abie_Bal"), ]
  }
  
  if (!suppliedElsewhere("sppEquivLong", sim)) {
    ## pass a default sppEquivalencies_CA for common species in western Canada
    sim$sppEquivLong <- LandR::sppEquivalencies_CA
  }

  if (!suppliedElsewhere("speciesEcoregion", sim)) {
    warning("generating dummy speciesEcoregion data - run Biomass_borealDataPrep for table with real speciesEcoregion attributes")
    sim$speciesEcoregion <- data.table(
      speciesCode = unique(sim$sppEquiv[[P(sim)$sppEquivCol]]),
      ecoregionGroup = "x",
      establishprob = 0.5,
      maxB = P(sim)$maxBInFactorial,
      maxANPP = P(sim)$maxBInFactorial / 30,
      year = 0
    )
  }

  ## check parameter consistency across modules
  paramCheckOtherMods(sim, "sppEquivCol", ifSetButDifferent = "error")

  if (!suppliedElsewhere("species", sim)) {
    message("generating dummy species data - run Biomass_borealDataPrep for table with real species attributes")
    speciesTable <- getSpeciesTable()
    sim$species <- prepSpeciesTable(
      speciesTable,
      sppEquiv = sim$sppEquiv,
      sppEquivCol = P(sim)$sppEquivCol
    )
  }

  if (!suppliedElsewhere("PSPmeasure_sppParams", sim) ||
      !suppliedElsewhere("PSPplot_sppParams", sim) ||
      !suppliedElsewhere("PSPgis_sppParams", sim)) {
    message("one or more PSP objects not supplied. Generating PSP data...")

    PSPdata <- getPSP(
      PSPdataTypes = P(sim)$PSPdataTypes,
      destinationPath = dPath,
      forGMCS = FALSE
    ) |>
      Cache(userTags = c(cacheTags, P(sim)$PSPdataTypes))

    sim$PSPmeasure_sppParams <- PSPdata$PSPmeasure
    sim$PSPplot_sppParams <- PSPdata$PSPplot
    sim$PSPgis_sppParams <- PSPdata$PSPgis
  }

  return(invisible(sim))
}
