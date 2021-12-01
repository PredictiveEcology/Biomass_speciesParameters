#' @importFrom data.table as.data.table
speciesTable2CohortGenerate <- function(Attributes = list(growthcurve = seq(0.65, 0.85, 0.02),
                                                          mortalityshape = seq(18, 25, 1),
                                                          longevity = seq(125, 300, 25),
                                                          mANPPproportion = seq(3.5, 6, 0.25)) ) {
  Atts <- names(Attributes)
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
  # species2a <- species2[species2$longevity1 > species2$longevity2 &
  #                         species2$growthcurve1 > species2$growthcurve2 &
  #                         species2$mortalityshape1 > species2$mortalityshape2 &
  #                         species2$mANPPproportion1 > species2$mANPPproportion2]
  
  set(species2a, NULL, "pixelGroup", seq(NROW(species2a)))
  set(species2a, NULL, "species", paste0("A", species2a$pixelGroup))
  browser()
  species2b <- melt(species2a, id.vars = c("species", "pixelGroup")) 
  set(species2b, NULL, "Sp", gsub(".+(1$|2$)", "Sp\\1", species2b$variable))
  set(species2b, NULL, "variable", gsub("1$|2$", "", species2b$variable))
  species2b <- dcast(species2b, pixelGroup + species + Sp ~ variable )
  set(species2b, NULL, "speciesCode", paste0(species2b$species, "_", species2b$Sp))
  set(species2b, NULL, c("species", "Sp"), NULL)
  cohortData2 <- copy(species2b[, c("speciesCode", "pixelGroup")])
  set(cohortData2, NULL, "pixelGroup", cohortData2$pixelGroup + max(cohortData$pixelGroup))
  cohortData <- rbindlist(list(cohortData, cohortData2), use.names = TRUE)
  setnames(species2b, old = "speciesCode", new = "species")
  colsToKeep <- c(Atts, "species")
  species2 <- species2b[, ..colsToKeep]
  species1 <- rbindlist(list(species1, species2), use.names = TRUE)
}
