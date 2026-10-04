## the module's package rendition (convertToPackage) does not import data.table; this makes the
## data.table syntax below, and in the tests, work there
.datatable.aware <- TRUE

## Latin names, defined once: the species BC's interior PSPs record, and the hybrid row of sppEquiv
## that LandR::speciesInStudyArea() merges into Pice_eng (see LandR.mergeHybridSpruce).
whiteSpruceLatin <- "Picea glauca"
hybridSpruceLatin <- "Picea engelmannii x glauca"

#' Relabel BC hybrid-zone white spruce PSP records as the hybrid spruce
#'
#' LandR merges the hybrid white x Engelmann spruce (`Pice_eng_gla`) into the species named by
#' `getOption("LandR.mergeHybridSpruce")`. BC's interior PSPs record the hybrid-zone spruce as
#' "Picea glauca", which is not in a `sppEquiv` where the hybrid is merged into Engelmann, so
#' those records were dropped. With the merge on `"engelmann"`, they are relabelled
#' "Picea engelmannii x glauca" here, which the existing `Latin_full` match gives `Pice_eng`.
#' BC plots in the BEC zones in `excludeBECzones` (boreal white spruce) are left alone, as are
#' plots from other sources and plots outside every BEC polygon. Only `"engelmann"` relabels:
#' with `"white"` the hybrid records already map to `Pice_gla`.
#'
#' @param PSPmeasure data.table of tree measurements, with the `speciesCol` column, `OrigPlotID1` and `source`.
#' @param PSPgis sf of plot locations, with `OrigPlotID1`.
#' @param sppEquiv the study area's `sppEquiv`; the relabel happens only if it has a row with
#'   `Latin_full` "Picea engelmannii x glauca" and `LandR` "Pice_eng" (the merged hybrid).
#' @param mergeHybridSprucePSP `"engelmann"` to relabel, `"white"` or `NA` for no relabel. The
#'   module parameter of that name defaults to `getOption("LandR.mergeHybridSpruce", "engelmann")`.
#' @param excludeBECzones BEC zones (`ZONE`) whose plots are not relabelled; the module parameter
#'   `excludeBECzonesHybridSpruce` holds the default.
#' @param BECzones optional sf of BEC polygons with a `ZONE` column; if `NULL`, they are fetched
#'   with [getBECzonesBC()].
#' @param speciesCol the column of `PSPmeasure` holding the Latin name (`Species`, or `Latin_full` after
#'   the `sppEquivLong` join in `prepPSPaNPP`).
#' @return `PSPmeasure`, as a copy, with `speciesCol` relabelled where it applies.
#' @keywords internal
relabelHybridSprucePSP <- function(PSPmeasure, PSPgis, sppEquiv,
                                   mergeHybridSprucePSP, excludeBECzones, BECzones = NULL,
                                   speciesCol = "Species") {
  hybridMerged <- nrow(sppEquiv[Latin_full == hybridSpruceLatin & LandR == "Pice_eng"]) > 0
  if (!isTRUE(mergeHybridSprucePSP == "engelmann") || !hybridMerged) return(PSPmeasure)
  candidate <- PSPmeasure$source %in% "BC" & PSPmeasure[[speciesCol]] %in% whiteSpruceLatin
  if (!any(candidate)) return(PSPmeasure)

  plots <- unique(PSPmeasure$OrigPlotID1[candidate])
  gis <- PSPgis[PSPgis$OrigPlotID1 %in% plots, "OrigPlotID1"]
  if (is.null(BECzones)) BECzones <- getBECzonesBC(gis)
  if (is.null(BECzones)) {
    warning("BC BEC zones are not available, so BC 'Picea glauca' PSP records were NOT relabelled ",
            "as the hybrid spruce (BWBS and SWB cannot be excluded). Supply `BECzonesBC`, or retry ",
            "when the BC Data Catalogue (WHSE_FOREST_VEGETATION.BEC_BIOGEOCLIMATIC_POLY) is up.",
            call. = FALSE)
    return(PSPmeasure)
  }
  zone <- sf::st_join(sf::st_transform(gis, sf::st_crs(BECzones)), BECzones["ZONE"])
  zone <- data.table::as.data.table(sf::st_drop_geometry(zone))[!duplicated(OrigPlotID1)]
  zoneOf <- setNames(as.character(zone$ZONE), zone$OrigPlotID1)

  PSPmeasure <- data.table::copy(PSPmeasure)
  plotZone <- zoneOf[as.character(PSPmeasure$OrigPlotID1)] # PSPclean gives a factor; its codes are not names
  relabel <- candidate & !is.na(plotZone) & !plotZone %in% excludeBECzones
  data.table::set(PSPmeasure, which(relabel), speciesCol, hybridSpruceLatin)

  counts <- PSPmeasure[relabel, .(plots = data.table::uniqueN(OrigPlotID1), records = .N),
                       by = .(zone = plotZone[relabel])][order(zone)]
  message(cli::col_yellow("Relabelled BC '", whiteSpruceLatin, "' PSP records as '", hybridSpruceLatin,
                          "' (merged into Pice_eng): ", sum(counts$records), " records in ",
                          sum(counts$plots), " plots (",
                          paste0(counts$zone, ": ", counts$plots, " plots", collapse = "; "),
                          "). Not relabelled: ", sum(candidate & !relabel), " records in BEC zones ",
                          paste(excludeBECzones, collapse = "/"), " or outside any BEC polygon."))
  PSPmeasure
}

#' BEC zone polygons of BC around some plots
#'
#' Queries the BC Data Catalogue layer `WHSE_FOREST_VEGETATION.BEC_BIOGEOCLIMATIC_POLY` for the
#' polygons in the bounding box of `pts`.
#'
#' @param pts sf of plot locations.
#' @return an sf with a `ZONE` column, or `NULL` (with a message) if the web service fails, as it
#'   did on 2026-10-03 ("HikariPool idwprod1 ... Connection is not available").
#' @keywords internal
getBECzonesBC <- function(pts) {
  ## the box is computed before the query: bcdata cannot translate a geometry built inside filter()
  box <- sf::st_as_sfc(sf::st_bbox(sf::st_transform(pts, 3005)))
  tryCatch(
    bcdata::bcdc_query_geodata("WHSE_FOREST_VEGETATION.BEC_BIOGEOCLIMATIC_POLY") |>
      bcdata::filter(bcdata::INTERSECTS(box)) |>
      bcdata::select("ZONE") |>
      bcdata::collect(),
    error = function(e) {
      message(cli::col_red("BEC polygons could not be fetched from the BC Data Catalogue: ",
                           conditionMessage(e)))
      NULL
    })
}
