## prepPSPaNPP() keeps a tree if its DBH reaches the larger of `minDBH` and its plot's own recorded
## minimum DBH, PSPplot's `minDBH` column (as PSPclean returns it). The filter used to run after the
## join that adds that column, which inside `[` masked the argument: each tree was compared with its
## plot's minimum only, so `minDBH` had no effect (BC plots, which record trees from 2 cm, kept them at
## `minDBH = 9`). Comparing each tree with `minDBH` alone would instead make the default `minDBH = 0`
## keep the trees under their plot's minimum (those of Alberta plots under 9.1 cm, for one).

## plots measured once, each with 10 trees of each DBH class (the forest filter needs >= 30 trees at
## the first measurement): one plot per recorded minimum DBH in `plotMinDBH`, or, if it is NULL, one
## plot and a plot table without the `minDBH` column
dbhClasses <- c(3, 7, 8.9, 9, 9.1, 10, 12, 15)
smallTreesPSP <- function(plotMinDBH) {
  plots <- if (is.null(plotMinDBH)) "noMin" else paste0("min", plotMinDBH)
  dbh <- rep(dbhClasses, each = 10)
  n <- length(dbh)
  plot <- list(MeasureID = seq_along(plots), OrigPlotID1 = plots, MeasureYear = 2000L, source = "BC",
               baseYear = 2000L, baseSA = 50, PlotSize = 0.1)
  plot$minDBH <- plotMinDBH ## NULL leaves the column out
  list(
    measure = data.table::data.table(
      MeasureID = rep(seq_along(plots), each = n), OrigPlotID1 = rep(plots, each = n), MeasureYear = 2000L,
      source = "BC", TreeNumber = rep(seq_len(n), length(plots)), DBH = rep(dbh, length(plots)),
      Height = 10, Species = "Picea glauca", status = "A", first_tree_year = 2000L,
      last_tree_year = 2010L, diff_dbh = 0),
    plot = data.table::as.data.table(plot)
  )
}

## prepPSPaNPP() on those plots
prepSmallTrees <- function(minDBH, plotMinDBH) {
  p <- smallTreesPSP(plotMinDBH)
  suppressMessages(prepPSPaNPP(
    studyAreaANPP = NULL, PSPgis = NULL, PSPmeasure = p$measure, PSPplot = p$plot,
    sppEquivLong = data.table::data.table(Latin_full = "Picea glauca", SpBiomassEq = "white spruce"),
    useHeight = TRUE, biomassModel = "Lambert2005", PSPperiod = c(1990, 2020), minDBH = minDBH))
}

## the DBH classes of the trees prepPSPaNPP() keeps, by plot
keptDBH <- function(minDBH, plotMinDBH) {
  out <- prepSmallTrees(minDBH, plotMinDBH)
  plots <- smallTreesPSP(plotMinDBH)$plot$OrigPlotID1
  lapply(split(out$DBH, factor(out$OrigPlotID1, levels = plots)), function(x) sort(unique(x)))
}

## the recorded minimum DBH of plots that record trees from 2 cm (as BC's do), from 7.1 cm (as
## Saskatchewan's since 1977), from 9.1 cm (as Alberta's), from 12 cm, and of a plot without one
plotMin <- c(2, 7.1, 9.1, 12, NA)

test_that("regression: the default minDBH = 0 keeps the trees it kept before the fix", {
  ## prepPSPaNPP calls biomassCalculation(), setkey() ... unqualified, as it does attached in a simList
  withr::local_package("pemisc")
  withr::local_package("data.table")
  ## before the fix each tree was compared with its plot's recorded minimum DBH only, so a plot
  ## without one kept no trees
  p <- smallTreesPSP(plotMin)
  preFix <- merge(p$measure, p$plot, by = c("MeasureID", "OrigPlotID1", "MeasureYear", "source"))
  preFix <- preFix[which(preFix$DBH >= preFix$minDBH), ]
  out <- prepSmallTrees(minDBH = 0, plotMin)
  expect_identical(sort(paste(out$OrigPlotID1, out$TreeNumber)),
                   sort(paste(preFix$OrigPlotID1, preFix$TreeNumber)))
  expect_identical(keptDBH(minDBH = 0, plotMin),
                   list(min2 = dbhClasses, min7.1 = c(8.9, 9, 9.1, 10, 12, 15), min9.1 = c(9.1, 10, 12, 15),
                        min12 = c(12, 15), minNA = numeric(0)))
})

test_that("regression: minDBH = 9 removes the smaller trees of the plots that record them, and no others", {
  withr::local_package("pemisc")
  withr::local_package("data.table")
  ## the plots recording trees from 2 and 7.1 cm lose those under 9 cm; the plots whose recorded
  ## minimum is above 9 cm keep their trees from it, as at minDBH = 0
  expect_identical(keptDBH(minDBH = 9, plotMin),
                   list(min2 = c(9, 9.1, 10, 12, 15), min7.1 = c(9, 9.1, 10, 12, 15),
                        min9.1 = c(9.1, 10, 12, 15), min12 = c(12, 15), minNA = numeric(0)))
})

test_that("without the plots' recorded minimum DBH, minDBH alone decides", {
  withr::local_package("pemisc")
  withr::local_package("data.table")
  expect_identical(keptDBH(minDBH = 0, plotMinDBH = NULL), list(noMin = dbhClasses))
  expect_identical(keptDBH(minDBH = 9, plotMinDBH = NULL), list(noMin = c(9, 9.1, 10, 12, 15)))
})
