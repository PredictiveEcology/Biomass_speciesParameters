## Growth statistics shared by the non-linear models fitted to PSP data and by the single-species
## Biomass_core factorial, so that the two are compared on identical ages.

#' Growth statistics of a biomass-by-age curve
#'
#' @param age Stand ages, strictly increasing.
#' @param B Biomass at `age`.
#'
#' @return A one-row `data.table`: `maxB` (the largest `B`), `peakInc` (the largest net biomass
#'   increment per year, `diff(B) / diff(age)`), and `r = peakInc / maxB`. `r` is unit-free, which
#'   is what makes it comparable between PSP biomass and the (arbitrarily scaled) factorial.
#' @noRd
curveStats <- function(age, B) {
  if (length(age) < 2L || anyNA(B)) {
    return(data.table(maxB = NA_real_, peakInc = NA_real_, r = NA_real_))
  }
  inc <- diff(B) / diff(age)
  data.table(maxB = max(B), peakInc = max(inc), r = max(inc) / max(B))
}

#' Ages at which fitted and simulated curves are compared
#'
#' Ages in `ages` that lie in the closed interval `standAgesForFitting` and do not exceed `longevity`.
#' @noRd
comparisonAges <- function(ages, standAgesForFitting, longevity) {
  ages <- sort(unique(ages))
  ages[ages >= min(standAgesForFitting) & ages <= max(standAgesForFitting) & ages <= longevity]
}

#' `curveStats()` of a fitted non-linear model
#' @param nlm A fitted model with `standAge` as its predictor.
#' @noRd
nlmCurveStats <- function(nlm, ages) {
  B <- as.numeric(predict(nlm, newdata = data.frame(standAge = ages)))
  curveStats(ages, B)
}

#' Read the single-species factorial into memory
#'
#' @param cohortDataFactorial,speciesTableFactorial `arrow` datasets, or data.frames.
#'
#' @return A list with `cohortData` (`speciesCode`, `age`, `B`, `pixelGroup`) and `speciesTable`
#'   (`species`, `targetSpecies`, `longevity`, `growthcurve`, `mortalityshape`, `mANPPproportion`,
#'   `pixelGroup`), both `data.table`s, restricted to the single-species pixelGroups.
#' @noRd
readSingleSpeciesFactorial <- function(cohortDataFactorial, speciesTableFactorial) {
  speciesTable <- as.data.table(dplyr::collect(speciesTableFactorial))
  if ("speciesCode" %in% names(speciesTable) && !"species" %in% names(speciesTable)) {
    setnames(speciesTable, "speciesCode", "species")
  }
  if (!"targetSpecies" %in% names(speciesTable)) {
    stop("speciesTableFactorial has no `targetSpecies` column. Biomass_speciesParameters needs the ",
         "factorial from Biomass_speciesFactorial with `factorialDesign = 'perSpecies'`.")
  }
  speciesTable <- speciesTable[, c("species", "targetSpecies", "longevity", "growthcurve",
                                   "mortalityshape", "mANPPproportion", "pixelGroup")]
  cohortData <- dplyr::collect(cohortDataFactorial) |> as.data.table()
  cohortData <- cohortData[, c("speciesCode", "age", "B", "pixelGroup")]
  cohortData[, speciesCode := as.character(speciesCode)]
  cohortData <- cohortData[speciesCode %in% speciesTable$species]
  list(cohortData = cohortData, speciesTable = speciesTable)
}

#' Summarise each single-species factorial pixelGroup
#'
#' For every combination of species, `growthcurve` and `mANPPproportion`: `r` (see `curveStats()`)
#' over the ages that are used for the PSP fit, and `inflationFactor`, `maxBInFactorial` over the
#' largest biomass the cohort reached at any age.
#'
#' @param factorial Output of `readSingleSpeciesFactorial()`.
#' @inheritParams comparisonAges
#' @param maxBInFactorial The `maxB` of the factorial simulations.
#'
#' @return A `data.table` with one row per factorial pixelGroup. `r` is `NA` when the cohort
#'   did not live through the comparison ages.
#' @noRd
factorialCurveStats <- function(factorial, standAgesForFitting, maxBInFactorial) {
  cd <- merge(factorial$cohortData, factorial$speciesTable, by.x = c("speciesCode", "pixelGroup"),
              by.y = c("species", "pixelGroup"))
  setorderv(cd, c("pixelGroup", "age"))

  ## a species' comparison ages are those of its factorial cohorts that fall in the fitting window
  refAges <- cd[, .(ages = list(comparisonAges(age, standAgesForFitting, longevity[1]))),
                by = "targetSpecies"]
  refAges <- setNames(refAges$ages, refAges$targetSpecies)

  stats <- cd[, {
    a <- refAges[[targetSpecies[1]]]
    inWindow <- age %in% a
    ## a cohort that died before the end of the window has no complete curve to compare
    cs <- if (sum(inWindow) == length(a)) curveStats(age[inWindow], B[inWindow]) else curveStats(0, NA_real_)
    cs[, lifetimeMaxB := max(B)]
  }, by = c("targetSpecies", "pixelGroup", "growthcurve", "mANPPproportion", "longevity", "mortalityshape")]
  stats[, inflationFactor := maxBInFactorial / lifetimeMaxB]
  stats[, lifetimeMaxB := NULL]
  stats[]
}

#' Find the `mANPPproportion` that reproduces a target `r`, for each growthcurve
#'
#' LandR's `maxANPP` is the gross ANPP, whereas the slope of a biomass-by-age curve fitted to PSP
#' data is net of mortality, so a fitted curve does not give `maxANPP` directly. How much smaller the
#' net peak is depends on `growthcurve`, `longevity` and `mortalityshape`. The factorial gives
#' `r` over a grid of `mANPPproportion` for each `growthcurve`; this inverts it.
#'
#' `r` is not monotone in `mANPPproportion`: past some value mortality, which scales with `maxANPP`,
#' removes more than the extra growth adds, and `r` falls again. Only the rising branch (up to the
#' `mANPPproportion` with the largest `r`) is inverted, as the falling branch has mortality close to
#' ANPP and cohorts that reach `maxB` within a few decades.
#'
#' @param curveStatsFactorial Output of `factorialCurveStats()`, for one species.
#' @param rTarget The target `r` (from the fitted PSP curve).
#'
#' @return A `data.table` with one row per `growthcurve`: `mANPPproportion` and `inflationFactor`
#'   (linearly interpolated along `mANPPproportion`), `rSim` (the `r` actually reached) and
#'   `clamped`, `TRUE` where `rTarget` is above the largest `r` of the factorial for that `growthcurve`
#'   (the traits at the largest `r` are returned) or below the smallest.
#' @noRd
matchManppProportion <- function(curveStatsFactorial, rTarget) {
  cs <- curveStatsFactorial[!is.na(r)]
  out <- cs[, {
    o <- order(mANPPproportion)
    prop <- mANPPproportion[o]
    rr <- r[o]
    infl <- inflationFactor[o]
    rising <- seq_len(which.max(rr))
    if (length(rising) < 2L) {
      list(mANPPproportion = NA_real_, inflationFactor = NA_real_, rSim = NA_real_, clamped = NA)
    } else {
      rng <- range(rr[rising])
      list(mANPPproportion = approx(rr[rising], prop[rising], xout = rTarget, rule = 2, ties = mean)$y,
           inflationFactor = approx(rr[rising], infl[rising], xout = rTarget, rule = 2, ties = mean)$y,
           rSim = min(max(rTarget, rng[1]), rng[2]),
           clamped = rTarget < rng[1] | rTarget > rng[2])
    }
  }, by = "growthcurve"]
  out[!is.na(mANPPproportion)][order(growthcurve)]
}

#' Interpolate matched traits at a given `growthcurve`
#' @param matched Output of `matchManppProportion()`.
#' @param growthcurve A single value; clamped to the range of `matched$growthcurve`.
#' @noRd
traitsAtGrowthcurve <- function(matched, growthcurve) {
  gc <- min(max(growthcurve, min(matched$growthcurve)), max(matched$growthcurve))
  data.table(growthcurve = growthcurve,
             mANPPproportion = approx(matched$growthcurve, matched$mANPPproportion, xout = gc)$y,
             inflationFactor = approx(matched$growthcurve, matched$inflationFactor, xout = gc)$y)
}

#' The LandR curve that corresponds to the chosen traits, scaled to the fitted curve
#'
#' For each species, the factorial pixelGroup closest in `growthcurve` and `mANPPproportion` to the chosen
#' traits, with its biomass scaled so that its largest value in the fitting window equals that of the fitted
#' non-linear model; and the model's own prediction at the same ages.
#'
#' @param factorial Output of `readSingleSpeciesFactorial()`.
#' @param chosen `data.table` (`species`, `growthcurve`, `mANPPproportion`).
#' @param nlms Named list of fitted models (one per species).
#' @inheritParams comparisonAges
#' @return `data.table` (`species`, `standAge`, `BscaledNonLinear`, `predNonLinear`).
#' @noRd
landrVsNlmCurves <- function(factorial, chosen, nlms, standAgesForFitting) {
  rbindlist(lapply(chosen$species, function(sp) {
    ch <- chosen[species == sp]
    st <- factorial$speciesTable[targetSpecies == sp]
    if (!nrow(st) || is.null(nlms[[sp]])) return(NULL)
    d <- (st$growthcurve - ch$growthcurve)^2 / var(st$growthcurve) +
      (st$mANPPproportion - ch$mANPPproportion)^2 / var(st$mANPPproportion)
    pg <- st$pixelGroup[which.min(d)]
    cd <- factorial$cohortData[pixelGroup == pg][order(age)]
    ages <- comparisonAges(cd$age, standAgesForFitting, st$longevity[1])
    nlmB <- as.numeric(predict(nlms[[sp]], newdata = data.frame(standAge = cd$age)))
    winB <- cd$B[cd$age %in% ages]
    nlmWin <- nlmB[cd$age %in% ages]
    data.table(species = sp, standAge = cd$age, BscaledNonLinear = cd$B * max(nlmWin) / max(winB),
               predNonLinear = nlmB)
  }))
}
