#' @export
#' @importFrom data.table as.data.table
factorialSpeciesTable <- function(growthcurve = seq(0.65, 0.85, 0.02),
                                 mortalityshape = seq(18, 25, 1),
                                 longevity = seq(125, 300, 25),
                                 mANPPproportion = seq(3.5, 6, 0.25),
                                 ...,
                                 cohortsPerPixel = 1:2) {

  if (FALSE) { # These are the original ones
    growthcurves <- seq(0, 1, 0.1)
    mortalityshapes <- seq(5, 25, 1)
    longevity <- seq(150, 700, 25)
    mANPPproportion <- seq(0.25,10, 0.25)
  }
  
  forms <- formals()
  formNames <- setdiff(names(forms), c("...", "cohortsPerPixel"))
  params <- append(forms[formNames], list(...))
  defaultAttributes <- lapply(params, eval)
  Attributes <- modifyList2(defaultAttributes, append(mget(formNames), list(...)))
  
  Atts <- names(Attributes)
  
  if (any(cohortsPerPixel %in% 1)) {
    species1 <- expand.grid(Attributes)
    species1 <- as.data.table(species1)
    setnames(species1, new = Atts)
  }
  
  if (any(cohortsPerPixel %in% 2)) {
    Attributes1 <- paste0(Atts, "1")
    Attributes2 <- paste0(Atts, "2")
    A <- Attributes # for ease of writing
    species2 <- expand.grid(append(A, A))
    species2 <- as.data.table(species2)
    colnames2sp <- c(Attributes1, Attributes2)
    setnames(species2, new = colnames2sp)
    
    # Take top right of matrix
    Map(A1 = Attributes1, A2 = Attributes2, function(A1, A2) {
      species2 <<- species2[species2[[A1]] > species2[[A2]]]
    })
  }
  
  prevMaxPG <- 0
  if (any(cohortsPerPixel %in% 1)) {
    set(species1, NULL, "pixelGroup", seq(NROW(species1)))
    set(species1, NULL, "species", paste0("A", species1$pixelGroup))
    prevMaxPG <- max(species1$pixelGroup)
  }
  
  if (any(cohortsPerPixel %in% 2)) {
    
    set(species2, NULL, "pixelGroup", seq(NROW(species2)))
    set(species2, NULL, "species", paste0("B", species2$pixelGroup))
    species2b <- melt(species2, id.vars = c("species", "pixelGroup")) 
    set(species2b, NULL, "Sp", gsub(".+(1$|2$)", "Sp\\1", species2b$variable))
    set(species2b, NULL, "variable", gsub("1$|2$", "", species2b$variable))
    species2b <- dcast(species2b, pixelGroup + species + Sp ~ variable )
    set(species2b, NULL, "speciesCode", paste0(species2b$species, "_", species2b$Sp))
    set(species2b, NULL, c("species", "Sp"), NULL)
    # cohortData2 <- copy(species2b[, c("speciesCode", "pixelGroup")])
    # set(cohortData2, NULL, "pixelGroup", cohortData2$pixelGroup + max(cohortData$pixelGroup))
    # cohortData <- rbindlist(list(cohortData, cohortData2), use.names = TRUE)
    setnames(species2b, old = "speciesCode", new = "species")
    set(species2b, NULL, "pixelGroup", species2b$pixelGroup + prevMaxPG)
  }
  
  speciesOut <- data.table()
  if (any(cohortsPerPixel %in% 1)) {
    speciesOut <- species1
  }
  if (any(cohortsPerPixel %in% 2)) {
    speciesOut <- rbindlist(list(speciesOut, species2b), use.names = TRUE)
  }
  
  set(speciesOut, NULL, "mortalityshape", asInteger(speciesOut$mortalityshape))
  set(speciesOut, NULL, "longevity", asInteger(speciesOut$longevity))
  speciesOut
}

#' 
#' Will create a cohortData table from a `speciesTable`, and combine it with another 
#' `cohortData` table (e.g., from singles)
#' Must have `species` and `pixelGroup`
factorialCohortData <- function(speciesTable, speciesEcoregion) {
  if (!identical(speciesTable$species, as.character(speciesEcoregion$speciesCode)))
    stop("speciesTable and speciesEcoregion must have identical species and as.character(speciesCode)")
  cohortData2 <- speciesTable[, c("species", "pixelGroup")]#, "maxANPP")]
  set(cohortData2, NULL, "speciesCode", as.factor(cohortData2$species))
  set(cohortData2, NULL, "species", NULL)
  
  set(cohortData2, NULL, "age", 1L)
  
  set(cohortData2, NULL, "ecoregionGroup", factor(1))
  set(cohortData2, NULL, "B", speciesEcoregion$maxANPP)
  setcolorder(cohortData2, c('speciesCode', 'pixelGroup', 'ecoregionGroup', 'age', "B"))
}


factorialSpeciesEcoregion <- function(speciesTable) {
  speciesEcoregion <- speciesTable[, c("species", "mANPPproportion")]
  
  set(speciesEcoregion, NULL, "maxB", 5000L)
  set(speciesEcoregion, NULL, "maxANPP", speciesEcoregion$maxB * speciesEcoregion$mANPPproportion/100)
  speciesEcoregion[, c("ecoregionGroup", "establishprob", "year") := .(factor(1), 0.5, 0)]
  setnames(speciesEcoregion, old = "species", new = "speciesCode")
  
  # Change classes
  set(speciesEcoregion, NULL, "speciesCode", as.factor(speciesEcoregion$speciesCode))
  set(speciesEcoregion, NULL, "mANPPproportion", NULL)
  speciesEcoregion[]
}


factorialSpeciesTableFillOut <- function(speciesTable) {
  speciesTableInner <- copy(speciesTable[, c("species", "longevity", "growthcurve", "mortalityshape")])
  speciesTableInner[, c("sexualmature", 'seeddistance_eff', 'seeddistance_max', 'resproutprob',
                   'resproutage_min', 'resproutage_max', 'postfireregen',
                   'leaflongevity', 'wooddecayrate', 'leafLignin', 'hardsoft', "Area", "firetolerance",
                   "shadetolerance") :=
                 list(30L, 0L, 0L, 0.5, 0L, 0L, factor('none'), 3L, 0.07, 0.1, factor('soft'), factor("BP"), 3L, 1L)]
  speciesTableInner[]
}

factorialPixelGroupMap <- function(cohortData) {
  pixelGroupMap <- raster(res = c(1,1))
  nrow(pixelGroupMap) <- round(sqrt(max(cohortData$pixelGroup)), 0)
  ncol(pixelGroupMap) <- round(sqrt(max(cohortData$pixelGroup)), 0) + 1
  vals <- c(1:max(cohortData$pixelGroup), rep(NA, times = ncell(pixelGroupMap) - max(cohortData$pixelGroup)))
  pixelGroupMap[] <- vals
  pixelGroupMap
}