# simulateYoungStands <- function(cohortData, N, x = 0.0001) {
#   #minAge was added as an argument so we can change it on a species basis
#   minAge = min(cohortData$standAge)
#   minB <-  min(cohortData[standAge == minAge, biomass])
#   simuAge <- round(runif(0, minAge, n = N), digits = 0)
#   
#   b = (minB/x) ^ (1/minAge)
#   
#   simuB <- x*b^simuAge
#   #add unique random effects 
#   potentialYears <- 1:2050
#   years <- sample(potentialYears[!potentialYears %in% cohortData$MeasureYear], size = N, replace = FALSE)
#   randomPlotIDs <- paste0('random', 1:N)
#   
#   simulatedData <- data.table(MeasureYear = years, OrigPlotID1 = randomPlotIDs, PlotSize = 1, 
#                               standAge = simuAge, biomass = simuB, spDom = 1)
#   
#   return(simulatedData)
# }

#Changed the function to simulate linear biomass based on minimum 2 points
simulateYoungStands <- function(cohortData, N) {
  
  #add unique random effects
  potentialYears <- 1:2050
  years <- sample(potentialYears[!potentialYears %in% cohortData$MeasureYear], size = N, replace = TRUE)
  randomPlotIDs <- paste0('random', 1:N)
  
  simulatedData <- data.table(MeasureYear = years, OrigPlotID1 = randomPlotIDs,
                              PlotSize = 1, standAge = 0, biomass = 0, spDom = 1)
  
  return(simulatedData)
}
