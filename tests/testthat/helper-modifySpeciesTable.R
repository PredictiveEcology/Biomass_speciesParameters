## A small, deterministic stand-in for the inputs of `modifySpeciesTable()` (focal fitting, two species):
## growth curves fitted to synthetic PSP biomass, and a factorial grid of traits with a simulated biomass-by-age
## curve per trait combination.
## `modifySpeciesTable()` and `editSpeciesTraits()` call these unqualified, as they are attached when the module runs
attachModuleDeps <- function(env = parent.frame()) {
  for (p in c("data.table", "ggplot2", "LandR")) withr::local_package(p, .local_envir = env)
}

makeModifySpeciesTableInputs <- function() {
  set.seed(1)
  ages <- seq(1, 200, by = 10)
  spp <- c(Abie_las = 0.035, Pice_eng = 0.018)   # true rate behind the synthetic PSP data

  ## synthetic PSP biomass (what the nls fits see) and the nls fits
  curve <- function(age, A, k, p) A * (1 - exp(-k * age))^p
  gcs <- lapply(names(spp), function(s) {
    age <- sample(25:120, 60, replace = TRUE)
    biomass <- curve(age, 20000 + 3000 * (s == "Pice_eng"), spp[[s]], 2) * exp(rnorm(60, 0, 0.05))
    d <- data.table::data.table(speciesTemp = factor(s, levels = c(s, "Other")), standAge = age, biomass = biomass,
                    OrigPlotID1 = "p")
    other <- data.table::copy(d)[, `:=`(speciesTemp = factor("Other", levels = c(s, "Other")), biomass = biomass / 2)]
    nlm <- nls(biomass ~ curve(standAge, A, k, 2), data = d, start = list(A = 15000, k = 0.03))
    list(originalData = rbind(d, other), NonLinearModel = setNames(list(nlm), s))
  })
  names(gcs) <- names(spp)

  ## factorial grid: one pixelGroup per trait combination, each with a "Sp1" and a "Sp2" row
  grid <- data.table::CJ(growthcurve = c(0.6, 0.65, 0.7, 0.75, 0.8), mANPPproportion = c(3, 4, 5, 6),
             longevity = c(150L, 300L), mortalityshape = c(18L, 22L))
  grid[, pixelGroup := .I]
  fT <- data.table::rbindlist(lapply(c("Sp1", "Sp2"), function(sp) {
    x <- data.table::copy(grid); x[, Sp := sp]; x[, speciesCode := paste0("B", pixelGroup, "_", sp)]
  }))
  fT[, species := speciesCode]
  fB <- fT[, .(age = ages[ages <= longevity],
               B = 5000 * (1 - exp(-0.003 * mANPPproportion * (1 + 2 * growthcurve) * ages[ages <= longevity]))^2 *
                 (1 - 0.1 * (Sp == "Sp2"))),
           by = .(speciesCode, pixelGroup, Sp, species)]
  fB[, species := speciesCode][, Sp := NULL]
  ## as `tempMaxB` in Init: the traits plus the inflation factor
  inflation <- data.table::copy(grid)[, `:=`(species = "x", inflationFactor = 1 + pixelGroup / 1000)]

  speciesTable <- data.table::data.table(species = c("Abie_las", "Pice_eng", "Betu_pap"), longevity = c(250L, 450L, 150L),
                             growthcurve = c(0.9, 0.9, 0.9), mortalityshape = c(15L, 15L, 15L),
                             mANPPproportion = c(9, 9, 9), inflationFactor = 1)
  list(GCs = gcs, speciesTable = speciesTable, factorialTraits = fT[, !"speciesCode"],  # as read from the feather: `species` only, renamed to `speciesCode` inside
       factorialBiomass = fB, sppEquiv = data.table::data.table(LandR = names(spp)), sppEquivCol = "LandR",
       inflationFactorKey = inflation, standAgesForFitting = c(21L, 91L), approach = "focal",
       maxBInFactorial = 5000L)
}
