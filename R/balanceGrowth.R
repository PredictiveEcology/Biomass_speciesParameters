## Shared growthcurve + K-balanced mANPPproportion (`P(sim)$balanceGrowth`).
##
## Under LandR competition a small cohort grows with K = mANPPproportion * maxB^(1 - growthcurve), and
## the competition makes any difference in K between species winner-take-all. Fitting each species alone
## ("focal") gives growthcurves that the data cannot tell apart, yet the small differences exclude
## species that coexist in plots. Here growthcurve is shared and mANPPproportion is set so K is equal.

## `[.data.table` only does its thing for callers it considers aware; the package rendition of the module
## (tests, R CMD check) does not import data.table, so say so
.datatable.aware <- TRUE

## digits `growthcurve` and `mANPPproportion` are rounded to (the default path rounds `growthcurve` to 2,
## see `modifySpeciesTable`)
balanceDigits <- 3L

#' Median maxB per species across `speciesEcoregion` rows
#'
#' @param speciesEcoregion data.table with `speciesCode` and `maxB`.
#' @return a named numeric vector, names are species codes.
medianMaxB <- function(speciesEcoregion) {
  mb <- speciesEcoregion[, list(maxB = stats::median(as.numeric(maxB), na.rm = TRUE)), by = "speciesCode"]
  stats::setNames(mb$maxB, mb$speciesCode)
}

#' Log-likelihood cost of a (growthcurve, mANPPproportion) pair, from the fit's likelihood table
#'
#' The factorial grid is coarse, so the cost is linearly interpolated: along `mANPPproportion` within the
#' two grid `growthcurve` values that bracket `growthcurve`, then linearly between those two. At each grid
#' point the best (smallest `llNonLinDelta`) over all other traits (partner traits, mortalityshape,
#' longevity) is used. `NA` if the pair is outside the grid.
#'
#' @param ll data.table of one species with `growthcurve`, `mANPPproportion`, `llNonLinDelta`.
#' @param gc,mANPP the `growthcurve` and `mANPPproportion` to evaluate.
#' @return `llNonLinDelta` at the point, relative to the best row of `ll`.
llAtPoint <- function(ll, gc, mANPP) {
  ll <- data.table::copy(ll)[, growthcurve := round(growthcurve, 6)]
  ll[, llNonLinDelta := llNonLinDelta - min(llNonLinDelta)]
  grid <- unique(ll$growthcurve)
  atGrid <- function(g) {
    p <- ll[growthcurve == g, list(d = min(llNonLinDelta)), by = "mANPPproportion"][order(mANPPproportion)]
    if (nrow(p) < 2L) return(NA_real_)
    stats::approx(p$mANPPproportion, p$d, mANPP)$y
  }
  lo <- suppressWarnings(max(grid[grid <= gc]))
  hi <- suppressWarnings(min(grid[grid >= gc]))
  if (!is.finite(lo) || !is.finite(hi)) return(NA_real_)
  dLo <- atGrid(lo)
  if (lo == hi) return(dLo)
  dLo + (gc - lo) / (hi - lo) * (atGrid(hi) - dLo)
}

#' Share growthcurve across species and set mANPPproportion so K is equal
#'
#' @param fitted data.table of the weighted, *unrounded* focal fits: `species`, `growthcurve`,
#'   `mANPPproportion`.
#' @param maxB named numeric, maxB per species (see `medianMaxB`).
#' @param ll data.table of the fit's likelihood table: `species`, `growthcurve`, `mANPPproportion`,
#'   `llNonLinDelta`.
#' @param sharedGrowthcurve the shared growthcurve; `NA` for the median of the fitted values.
#' @param targetK the K every species gets; `NA` for the median over species of K at the shared
#'   growthcurve and each species' fitted mANPPproportion.
#' @param warnDeltaLL warn for species whose new traits are worse than their best fit by more than this.
#' @return a data.table with one row per fitted species: fitted and new `growthcurve` and
#'   `mANPPproportion`, `maxB` used, `K`, `deltaLL` (log-likelihood units worse than that species' best).
balanceGrowthK <- function(fitted, maxB, ll, sharedGrowthcurve = NA_real_, targetK = NA_real_,
                           warnDeltaLL = 2) {
  chk <- data.table::data.table(species = fitted$species,
                                growthcurveFitted = fitted$growthcurve,
                                mANPPproportionFitted = fitted$mANPPproportion)
  chk[, maxB := unname(maxB[species])]
  ## a fitted species with no speciesEcoregion row (e.g. below the ecoregion support threshold
  ## everywhere) cannot establish, so it is left out of the balance and keeps its fitted traits
  if (anyNA(chk$maxB)) {
    message("balanceGrowth: not in speciesEcoregion, so not balanced, traits unchanged: ",
            paste(chk$species[is.na(chk$maxB)], collapse = ", "))
    chk <- chk[!is.na(maxB)]
  }
  if (is.na(sharedGrowthcurve)) sharedGrowthcurve <- stats::median(chk$growthcurveFitted)
  if (is.na(targetK)) {
    targetK <- stats::median(chk$mANPPproportionFitted * chk$maxB^(1 - sharedGrowthcurve))
  }
  chk[, growthcurve := round(sharedGrowthcurve, balanceDigits)]
  chk[, mANPPproportion := round(targetK / maxB^(1 - growthcurve), balanceDigits)]
  chk[, K := mANPPproportion * maxB^(1 - growthcurve)]
  chk[, deltaLL := mapply(function(s, g, m) llAtPoint(ll[ll$species == s], g, m),
                          species, growthcurve, mANPPproportion)]
  data.table::setcolorder(chk, c("species", "growthcurveFitted", "growthcurve",
                                 "mANPPproportionFitted", "mANPPproportion", "maxB", "K", "deltaLL"))
  data.table::setnames(chk, c("growthcurve", "mANPPproportion"), c("growthcurveNew", "mANPPproportionNew"))

  for (i in seq_len(nrow(chk))) {
    with(chk[i], message(sprintf(
      "balanceGrowth %s: growthcurve %.3f -> %.3f, mANPPproportion %.3f -> %.3f, K = %.1f (deltaLL %s)",
      species, growthcurveFitted, growthcurveNew, mANPPproportionFitted, mANPPproportionNew, K,
      if (is.na(deltaLL)) "NA" else sprintf("%.2f", deltaLL))))
  }
  worse <- chk[!is.na(deltaLL) & deltaLL > warnDeltaLL]
  if (nrow(worse)) {
    warning("balanceGrowth: the shared growthcurve with K-balanced mANPPproportion fits the PSP data worse ",
            "than the species' own best fit by more than ", warnDeltaLL, " log-likelihood units for: ",
            paste0(worse$species, " (", round(worse$deltaLL, 1), ")", collapse = ", "), call. = FALSE)
  }
  outside <- chk[is.na(deltaLL)]
  if (nrow(outside)) {
    warning("balanceGrowth: the new traits are outside the factorial grid, so the fit could not be ",
            "checked for: ", paste(outside$species, collapse = ", "), call. = FALSE)
  }
  chk
}

#' Apply the balanced traits to the species table
#'
#' @param speciesTable the `species` table.
#' @param check output of `balanceGrowthK`.
#' @return `speciesTable` with `growthcurve` and `mANPPproportion` replaced for the fitted species;
#'   other species are untouched, and reported.
applyBalancedTraits <- function(speciesTable, check) {
  speciesTable <- data.table::copy(speciesTable)
  idx <- match(check$species, speciesTable$species)
  data.table::set(speciesTable, idx, "growthcurve", check$growthcurveNew)
  data.table::set(speciesTable, idx, "mANPPproportion", check$mANPPproportionNew)
  unfitted <- setdiff(speciesTable$species, check$species)
  if (length(unfitted)) {
    message("balanceGrowth: not fitted, traits unchanged: ", paste(unfitted, collapse = ", "))
  }
  speciesTable
}
