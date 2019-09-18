library(ggplot2)
library(magrittr)
library(data.table)
library(qgam)
library(reproducible)

#Exploring PSP data, see if I can make growth curves.
#To make dataset below, had to put browser in gmcsDataPrep and output a dataset midway through data cleaning
#Also saved PSPgis after getting coordinates using st_coordinates function
#There is still no data for Quebec, Ontario, Maritimes, and non-Boreal BC - Aug 01, 2019


PSPgis <- as.data.table(prepInputs(url = 'https://drive.google.com/open?id=1pJwNzzJXZCROtVug00cTfxpyaUyXmwWe', 
                                   destinationPath = tempdir(), fun = 'readRDS',
                                   targetFile = 'PSPgis.Rdat'))
#stored on 254 as C:/Ian/Data/PSP/PSPgis.Rdat

psp <- prepInputs(url = 'https://drive.google.com/open?id=1VwV3NUBiPuC3N7-WajqWUig3T0DJCfJZ',
                  destinationPath = tempdir(), fun = 'readRDS', targetFile = 'PSPs_withAgeandBiomass.RDat')
#stored on 254 as C:/Ian/Data/PSP/PSPs_withAgeandBiomass.Rdat

#Need to join the psps with PSPgis to get the actual coords with same crs
psp <- PSPgis[psp, on = c("OrigPlotID1" = 'OrigPlotID1')]
#remove some redundant columns from before the merge + useless GIS
psp[, c('geometry', 'Latitude', 'Longitude', 'Easting', 'Northing','i.MeasureYear', 'i.OrigPlotID1') := NULL]

#standardize biomass by plotsize
#Note: there are still obvious errors in PSP data. e.g. plotSize 0.0029 ha, but 90 trees? yeah right buddy
psp$Biomass <- psp$Biomass/psp$PlotSize/10 #puts in LandR units of g/m2 and scales by plot size
psp$newSpeciesName <- as.factor(psp$newSpeciesName)
psp <- psp[standAge > 0]
psp <- psp[!is.na(Biomass)]

#Need to take pure cohorts, but also sum of biomass in each stand, not individual trees
cpsp <- data.table::copy(psp)
cpsp[, "sumB" := sum(Biomass), .(MeasureID)]
cpsp[, "spSumB" := sum(Biomass), .(MeasureID, newSpeciesName)]
cpsp[, "spDom" := spSumB/sumB, .(MeasureID)]
cpsp$newSpeciesName <- as.factor(cpsp$newSpeciesName)
pureStands <- cpsp[spDom > 0.60]

X <- 10^-6 #this is a constant used to simulate B, idea developed by Ceres
#Create functions to simulate young stands (underrepresented in data)
simulateB <- function(cohortData, n, x = X, minAge = min(cohortData$standAge)) {
  #minAge was added as an argument so we can change it on a species basis
  minB <-  min(cohortData[standAge == minAge, sumB])
  simuAge <- round(runif(0, minAge, n = 1500), digits = 0)
  
  b = (minB/x) ^ (1/minAge)
  
  simuB <- x*b^simuAge
  
  fakeData <- data.table(spSumB = 0, sumB = simuB, MeasureID = 1, spDom = 1,
                         standAge = simuAge, lat = 55, lon = -9999)
  return(fakeData)
}


####trembling aspen qgam#####
PopulusTremuloides <- grep(pureStands$newSpeciesName, pattern = "trembling aspen") %>%
  pureStands[.] %>%
  .[, .(spSumB, sumB, MeasureID, spDom, standAge, lat, lon)] %>%
  .[!duplicated(.)]

fakeAspenData <- simulateB(cohortData = PopulusTremuloides, n = 500, x = X) #changed this to 500 due to some errors...
PopulusTrem <- rbind(fakeAspenData, PopulusTremuloides)

aspengam <- qgam(sumB ~ s(standAge), data = PopulusTrem, qu = 0.90) #I occasionally hit errors when adding the simulated data...
predAspen <- data.frame("standAge" = seq(from = 0, to = 150, by = .1))
predAspen$Biomass <- predict(aspengam, newdata = predAspen)
populusTaNPP <- round(max(diff(predAspen$Biomass, lag = 1))/0.1/max(predAspen$Biomass), digits = 4)
ageAtMaxANPP <- predAspen$standAge[diff(predAspen$Biomass, lag = 1) == max(diff(predAspen$Biomass, lag = 1))]

g1 <- ggplot(data = predAspen, aes(y = Biomass, x =  standAge)) +
  geom_point(data = PopulusTremuloides, aes(y = sumB, x = standAge, col = lat)) +
  scale_color_gradient(low = "green", high = "blue") +
  geom_point(data = fakeAspenData[1:5,], aes(y = sumB, x = standAge)) +
  geom_line(size = 2) +
  theme_bw() +
  labs(title = "90th quantile of Trembling aspen biomass by stand age",
       subtitle = paste0("max ANPP = ", populusTaNPP," of maxB at age ", ageAtMaxANPP),
       color = "latitude", ylab = "biomass (g/m2)")
g1
#check that it is reasonably close to 0
predAspen$Biomass[predAspen$standAge == 0]/max(predAspen$Biomass) 
#1.4% that's not great but I can't seem to lower it



#### lodgepole Pine qgam ####
PinusContorta <- grep(pureStands$newSpeciesName, pattern = "lodgepole pine") %>%
  pureStands[.] %>%
  .[, .(spSumB, sumB, MeasureID, spDom, standAge, lat, lon)] %>%
  .[!duplicated(.)]
#There is one obvious poor point (2 much older trees make up 90% of biomass in 20-year old stand of 'not pine')
PinusContorta <- PinusContorta[sumB != max(sumB)]

fakePineCData <- simulateB(cohortData = PinusContorta, n = 1500, x = 10^-2)

PinusCon <- rbind(PinusContorta, fakePineCData)

pinusCgam <- qgam(sumB ~ s(standAge), data = PinusCon, qu = 0.90)
predPinusC <- data.frame("standAge" = seq(from = 0, to = 150, by = .1))
predPinusC$Biomass <- predict(pinusCgam, newdata = predPinusC)
pinusCaNPP <- round(max(diff(predPinusC$Biomass, lag = 1))/0.1/max(predPinusC$Biomass), digits = 4)
ageAtMaxANPP <- predPinusC$standAge[diff(predPinusC$Biomass, lag = 1) == max(diff(predPinusC$Biomass, lag = 1))]


g1 <- ggplot(data = predPinusC, aes(y = Biomass, x =  standAge)) +
  geom_point(data = PinusContorta, aes(y = sumB, x = standAge, col = lat)) +
  geom_line(size = 2) +
    scale_color_gradient(low = "green", high = "blue") +
  theme_bw() +
  labs(title = "90th quantile of Lodgepole pine biomass by stand age",
       subtitle = paste0("max ANPP = ", pinusCaNPP," of maxB at age ", ageAtMaxANPP),
       color = "latitude", ylab = "biomass (g/m2)")
g1
#check that it is reasonably close to 0 - very close!
predPinusC$Biomass[predPinusC$standAge == 0]/max(predPinusC$Biomass) 


### Jack pine qgam####
PinusBanksiana <- grep(pureStands$newSpeciesName, pattern = "jack pine") %>%
  pureStands[.] %>%
  .[, .(spSumB, sumB, MeasureID, spDom, standAge, lat, lon)] %>%
  .[!duplicated(.)]

fakePineBData <- simulateB(cohortData = PinusBanksiana, n = 1500)

PinusBank <- rbind(PinusBanksiana, fakeData)

PinusBgam <- qgam(sumB ~ s(standAge), data = PinusBank, qu = 0.90)
predPinusB <- data.frame("standAge" = seq(from = 0, to = 100, by = .1))
predPinusB$Biomass <- predict(PinusBgam, newdata = predPinusB)
PinusBaNPP <- round(max(diff(predPinusB$Biomass, lag = 1))/0.1/max(predPinusB$Biomass), digits = 4)
ageAtMaxANPP <- predPinusB$standAge[diff(predPinusB$Biomass, lag = 1) == max(diff(predPinusB$Biomass, lag = 1))]

g1 <- ggplot(data = predPinusB, aes(y = Biomass, x =  standAge)) +
  scale_color_gradient(low = "green", high = "blue") +
  geom_point(data = PinusBanksiana, aes(y = sumB, x = standAge, col = lat)) +
  geom_line(size = 2) +
  theme_bw() +
  labs(title = "90th quantile of Jack pine biomass by stand age",
       subtitle = paste0("max ANPP = ", PinusBaNPP," of maxB at age ", ageAtMaxANPP),
       color = "latitude", ylab = "biomass (g/m2)")
g1

#check that it is reasonably close to 0
predPinusB$Biomass[predPinusB$standAge == 0]/max(predPinusB$Biomass) 

#### Black Spruce qgam ####
PiceaMariana <- grep(pureStands$newSpeciesName, pattern = "black spruce") %>%
  pureStands[.] %>%
  .[, .(spSumB, sumB, MeasureID, spDom, standAge, lat, lon)] %>%
  .[!duplicated(.)]

fakePiceaMData <- simulateB(PiceaMariana, n = 1500, x = 10^-8)
PiceaMar <- rbind(PiceaMariana, fakePiceaMData)

PiceaMgam <- qgam(sumB ~ s(standAge, k = 6), data = PiceaMar, qu = 0.90)
predPiceaM <- data.frame("standAge" = seq(from = 0, to = 150, by = .1))
predPiceaM$Biomass <- predict(PiceaMgam, newdata = predPiceaM)
PiceaMaNPP <- round(max(diff(predPiceaM$Biomass, lag = 1))/0.1/max(predPiceaM$Biomass), digits = 4)
ageAtMaxANPP <- predPiceaM$standAge[diff(predPiceaM$Biomass, lag = 1) == max(diff(predPiceaM$Biomass, lag = 1))]


g1 <- ggplot(data = predPiceaM, aes(y = Biomass, x =  standAge)) +
  geom_point(data = PiceaMariana, aes(y = sumB, x = standAge, col = lat)) +
  scale_color_gradient(low = "green", high = "blue") +
  geom_line(size = 2) +
  theme_bw() +
  labs(title = "90th quantile of Black spruce biomass by stand age",
       subtitle = paste0("max ANPP = ", PiceaMaNPP," of maxB at age ", ageAtMaxANPP),
       color = "latitude", ylab = "biomass (g/m2)")
g1
#check that it is reasonably close to 0
predPiceaM$Biomass[predPiceaM$standAge == 0]/max(predPiceaM$Biomass) 
# 2.9 percent  

#### White Spruce qgam####
PiceaGlauca <- grep(pureStands$newSpeciesName, pattern = "white spruce") %>%
  pureStands[.] %>%
  .[, .(spSumB, sumB, MeasureID, spDom, standAge, lat, lon)] %>%
  .[!duplicated(.)]

fakePiceaGData <- simulateB(PiceaGlauca, 1500)
PiceaGla <- rbind(PiceaGlauca, fakePiceaGData)

PiceaGgam <- qgam(sumB ~ s(standAge), data = PiceaGla, qu = 0.90)
predPiceaG <- data.frame("standAge" = seq(from = 0, to = 150, by = .1))
predPiceaG$Biomass <- predict(PiceaGgam, newdata = predPiceaG)
PiceaGaNPP <- round(max(diff(predPiceaG$Biomass, lag = 1))/0.1/max(predPiceaG$Biomass), digits = 5)
ageAtMaxANPP <- predPiceaG$standAge[diff(predPiceaG$Biomass, lag = 1) == max(diff(predPiceaG$Biomass, lag = 1))]


g1 <- ggplot(data = predPiceaG, aes(y = Biomass, x =  standAge)) +
  geom_point(data = PiceaGlauca, aes(y = sumB, x = standAge, col = lat)) +
  scale_color_gradient(low = "green", high = "blue") +
  geom_line(size =2 ) +
  theme_bw() +
  labs(title = "90th quantile of White spruce biomass by stand age",
       subtitle = paste0("max ANPP = ", PiceaGaNPP," of maxB at age ", ageAtMaxANPP),
       color = "latitude", ylab = "biomass (g/m2)")
g1
#check that it is reasonably close to 0
predPiceaG$Biomass[predPiceaG$standAge == 0]/max(predPiceaG$Biomass) 
#1.3%

#### sub-alpine fir qgam ####
AbiesLasiocarpa <- grep(pureStands$newSpeciesName, pattern = "alpine fir") %>%
  pureStands[.] %>%
  .[, .(spSumB, sumB, MeasureID, spDom, standAge, lat, lon)] %>%
  .[!duplicated(.)]

fakeAbiseLData <- simulateB(AbiesLasiocarpa, 1500)
AbiesLasio <- rbind(AbiesLasiocarpa, fakeAbiesLData)


AbiesLgam <- qgam(sumB ~ s(standAge), data = AbiesLasio, qu = 0.90)
predAbiesL <- data.frame("standAge" = seq(from = 0, to = 150, by = .1))
predAbiesL$Biomass <- predict(AbiesLgam, newdata = predAbiesL)
AbiesLaNPP <- round(max(diff(predAbiesL$Biomass, lag = 1))/0.1/max(predAbiesL$Biomass), digits = 5)
ageAtMaxANPP <- predAbiesL$standAge[diff(predAbiesL$Biomass, lag = 1) == max(diff(predAbiesL$Biomass, lag = 1))]


g1 <- ggplot(data = predAbiesL, aes(y = Biomass, x =  standAge)) +
  geom_point(data = AbiesLasiocarpa, aes(y = sumB, x = standAge, col = lat)) +
  geom_smooth(method = "loess", size = 2, colour = 'black') +
  scale_color_gradient(low = "green", high = "blue") +
  theme_bw() +
  labs(title = "90th quantile of Sub-alpine fir biomass by stand age",
       subtitle = paste0("max ANPP = ", AbiesLaNPP," of maxB at age ", ageAtMaxANPP),
       color = "latitude", ylab = "biomass (g/m2)")
g1

#check that it is reasonably close to 0
predAbiesL$Biomass[predAbiesL$standAge == 0]/max(predAbiesL$Biomass) 



##### balsam fir gqam####
AbiesBalsamea <- grep(pureStands$newSpeciesName, pattern = "balsam fir") %>%
  pureStands[.] %>%
  .[, .(spSumB, sumB, MeasureID, spDom, standAge, lat, lon)] %>%
  .[!duplicated(.)]
fakeAbiesBData <- simulateB(AbiesBalsamea, 1500)
AbiesBalsam <- rbind(AbiesBalsamea, fakeAbiesBData)

AbiesBgam <- qgam(sumB ~ s(standAge), data = AbiesBalsam, qu = 0.90)
predAbiesB <- data.frame("standAge" = seq(from = 0, to = 170, by = .1))
predAbiesB$Biomass <- predict(AbiesBgam, newdata = predAbiesB)
AbiesBaNPP <- round(max(diff(predAbiesB$Biomass, lag = 1))/0.1/max(predAbiesB$Biomass), digits = 5)
ageAtMaxANPP <- predAbiesB$standAge[diff(predAbiesB$Biomass, lag = 1) == max(diff(predAbiesB$Biomass, lag = 1))]


g1 <- ggplot(data = predAbiesB, aes(y = Biomass, x =  standAge)) +
  geom_point(data = AbiesBalsamea, aes(y = sumB, x = standAge, col = lat)) +
  scale_color_gradient(low = "green", high = "blue") +
  geom_smooth(method = "loess", size = 2, colour = 'black') +
  theme_bw() +
  labs(title = "90th quantile of Balsam fir biomass by stand age",
       subtitle = paste0("max ANPP = ", AbiesBaNPP," of maxB at age ", ageAtMaxANPP),
       color = "lat")
g1

#check that it is reasonably close to 0
predAbiesB$Biomass[predAbiesB$standAge == 0]/max(predAbiesB$Biomass) 


##### paper birch gqam####
BetulaPapyrifera <- grep(pureStands$newSpeciesName, pattern = "white birch") %>%
  pureStands[.] %>%
  .[, .(spSumB, sumB, MeasureID, spDom, standAge, lat, lon)] %>%
  .[!duplicated(.)]
fakeBetulaPData <- simulateB(BetulaPapyrifera, 1500)
BetulaPap <- rbind(BetulaPapyrifera, fakeBetulaPData)

BetulaPgam <- qgam(sumB ~ s(standAge), data = BetulaPap, qu = 0.90)
predBetulaP <- data.frame("standAge" = seq(from = 0, to = 100, by = .1))
predBetulaP$Biomass <- predict(BetulaPgam, newdata = predBetulaP)
BetulaPaNPP <- round(max(diff(predBetulaP$Biomass, lag = 1))/0.1/max(predBetulaP$Biomass), digits = 5)
ageAtMaxANPP <- predBetulaP$standAge[diff(predBetulaP$Biomass, lag = 1) == max(diff(predBetulaP$Biomass, lag = 1))]

g1 <- ggplot(data = predBetulaP, aes(y = Biomass, x =  standAge)) +
  geom_point(data = BetulaPapyrifera, aes(y = sumB, x = standAge, col = lat)) +
  scale_color_gradient(low = "green", high = "blue") +
  geom_smooth(method = "loess", size = 2, colour = 'black') +
  theme_bw() +
  labs(title = "90th quantile of White birch biomass by stand age",
       subtitle = paste0("max ANPP = ", BetulaPaNPP," of maxB at age ", ageAtMaxANPP),
       color = "lat")
g1

#check that it is reasonably close to 0
predBetulaP$Biomass[predBetulaP$standAge == 0]/max(predBetulaP$Biomass) 

#####old plots#########

cpsp <- data.table::copy(psp)
cpsp[, c('Easting', 'Northing', 'i.MeasureYear', 'i.OrigPlotID1', 'Longitude', 'Latitude') := NULL]
cpsp[, "sumB" := sum(Biomass), .(MeasureID)]
cpsp[, "spSumB" := sum(Biomass), .(MeasureID, newSpeciesName)]
cpsp[, "spDom" := spSumB/sumB, .(MeasureID)]
cpsp$newSpeciesName <- as.factor(cpsp$newSpeciesName)

pureStands <- cpsp[spDom > 0.6]

######Black Spruce Stand B#####
PiceaMariana <- grep(pureStands$newSpeciesName, pattern = "black spruce") %>%
  pureStands[.] %>%
  .[, .(spSumB, sumB, MeasureID, spDom, standAge)] %>%
  .[!duplicated(.)]

PiceaMariana$standAge <- round(PiceaMariana$standAge/10, digits = 0) * 10
PiceaMarianaQ <- PiceaMariana[, .("quantileB" = quantile(sumB, 0.90)), .(standAge)]


g1 <- ggplot(data = PiceaMarianaQ, aes(y = quantileB, x =  standAge)) +
  # geom_point() +
  geom_smooth(method = "loess") +
  geom_point(data = PiceaMariana, aes(y = sumB, x = standAge, col = spDom)) +
  theme_bw() +
  labs(title = "90th quantile of Black spruce biomass by stand age",
       color = "Prop. of stand")
g1


####### White Spruce Stand B#####
PiceaGlauca <- grep(pureStands$newSpeciesName, pattern = "white spruce") %>%
  pureStands[.] %>%
  .[, .(spSumB, sumB, MeasureID, spDom, standAge)] %>%
  .[!duplicated(.)]

PiceaGlauca$standAge <- round(PiceaGlauca$standAge/10, digits = 0) * 10
PiceaGlaucaQ <- PiceaGlauca[, .("quantileB" = quantile(sumB, 0.90)), .(standAge)]


g1 <- ggplot(data = PiceaGlaucaQ, aes(y = quantileB, x =  standAge)) +
  # geom_point() +
  geom_smooth(method = "loess") +
  geom_point(data = PiceaGlauca, aes(y = sumB, x = standAge, col = spDom)) +
  theme_bw() +
  labs(title = "90th quantile of White spruce biomass by stand age",
       color = "Prop. of stand")
g1

#####Lodgepole Pine Stand B ####
PinusContorta <- grep(pureStands$newSpeciesName, pattern = "lodgepole pine") %>%
  pureStands[.] %>%
  .[, .(spSumB, sumB, MeasureID, spDom, standAge)] %>%
  .[!duplicated(.)]

PinusContorta$standAge <- round(PinusContorta$standAge/10, digits = 0) * 10
PinusContortaQ <- PinusContorta[, .("quantileB" = quantile(sumB, 0.90)), .(standAge)]


g1 <- ggplot(data = PinusContortaQ, aes(y = quantileB, x =  standAge)) +
  # geom_point() +
  geom_smooth(method = "loess") +
  geom_point(data = PinusContorta, aes(y = sumB, x = standAge, col = spDom)) +
  theme_bw() +
  labs(title = "90th quantile of Lodgepole pine biomass by stand age",
       color = "Prop. of stand")
g1
PinusContortaQ[quantileB == max(quantileB)]

###Trembling Aspen Stand B####
PopulusTremuloides <- grep(pureStands$newSpeciesName, pattern = "trembling aspen") %>%
  pureStands[.] %>%
  .[, .(spSumB, sumB, MeasureID, spDom, standAge)] %>%
  .[!duplicated(.)]

PopulusTremuloides$standAge <- round(PopulusTremuloides$standAge/10, digits = 0) * 10
PopulusTremuloidesQ <- PopulusTremuloides[, .("quantileB" = quantile(sumB, 0.90)), .(standAge)]


g1 <- ggplot(data = PopulusTremuloidesQ, aes(y = quantileB, x =  standAge)) +
  # geom_point() +
  geom_smooth(method = "loess") +
  geom_point(data = PopulusTremuloides, aes(y = sumB, x = standAge, col = spDom)) +
  theme_bw() +
  labs(title = "90th quantile of Trembling aspen biomass by stand age",
       color = "Prop. of stand")
g1

##### Abies Lasoicarpa Stand B #####
AbiesLasiocarpa <- grep(pureStands$newSpeciesName, pattern = "alpine fir") %>%
  pureStands[.] %>%
  .[, .(spSumB, sumB, MeasureID, spDom, standAge)] %>%
  .[!duplicated(.)]

AbiesLasiocarpa$standAge <- round(AbiesLasiocarpa$standAge/10, digits = 0) * 10
AbiesLasiocarpaQ <- AbiesLasiocarpa[, .("quantileB" = quantile(sumB, 0.90)), .(standAge)]


g1 <- ggplot(data = AbiesLasiocarpaQ, aes(y = quantileB, x =  standAge)) +
  # geom_point() +
  geom_smooth(method = "loess") +
  geom_point(data = AbiesLasiocarpa, aes(y = sumB, x = standAge, col = spDom)) +
  theme_bw() +
  labs(title = "90th quantile of Subalpine fir biomass by stand age",
       color = "Prop. of stand")
g1


