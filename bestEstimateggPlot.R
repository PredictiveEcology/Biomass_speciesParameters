#big table is just the two tables joined after running simInit
library(gridExtra)
bigTable <- sim$factorialSpeciesTable[sim$reducedFactorialCohortData, on = c("species" = "speciesCode")]

Abies_las <- bigTable[growthcurve == 0.5 & mortalityshape == 23 & longevity == 250 & mANPPproportion == 2.50,]
Abies_las$species <- "Abies_las"

Betu_pap <-  bigTable[mortalityshape == 14 & longevity == 150 & mANPPproportion == 2.00 & growthcurve %==% 0.5,]
Betu_pap$species <- "Betu_pap"

Picea_engelmanii <- bigTable[growthcurve == 1 & mortalityshape == 15 & longevity == 450  & mANPPproportion == 3.25,] 
Picea_engelmanii$species <- "Pice_eng"

Pice_gla <- bigTable[growthcurve == 0.5 & mortalityshape == 8 & longevity == 400  & mANPPproportion == 2.50,] 
Pice_gla$species <- "Pice_gla"

Pice_mar <- bigTable[growthcurve %==% 0.5 & mortalityshape == 12 & longevity == 250 & mANPPproportion == 4.00,]
Pice_mar$species <- "Pice_mar"

Pine_con <- bigTable[growthcurve %==% 0.5 & mortalityshape %==% 7 & longevity == 325 & mANPPproportion %==% 1.25,]
Pine_con$species <- "Pine_con"

Popu_tre <- bigTable[growthcurve %==% 0.5 & mortalityshape %==% 15 & longevity == 200 & mANPPproportion %==% 1.50,]
Popu_tre$species <- "Popu_tre"

estimatedSpecies <- rbind(Abies_las, Betu_pap, Picea_engelmanii, Pice_gla, Pice_mar, Pine_con, Popu_tre)

Bmultiplier <- estimatedSpecies[, .('Bmultiplier' = 5000/max(B)), .(species)]
newEstimatedSpecies <- Bmultiplier[estimatedSpecies, on = "species"]
newEstimatedSpecies[, scaledB := B * Bmultiplier]

myPlot <- ggplot(data = newEstimatedSpecies, aes(y = scaledB, x = age, col = species)) + geom_line(size = 1) + theme_bw()
myPlot


####Some plots to show how growth curve, mortality shape, and anpp work

varyingGCs <- bigTable[longevity == 200 & mortalityshape == 15 & mANPPproportion %==% 3.0,]
varyingGCs$stat <- as.factor(varyingGCs$growthcurve)
varyingGCs$cat <- 'growthcurve'
varyingMS <- bigTable[longevity == 200 & growthcurve %==% 0.5 & mANPPproportion %==% 3.0,]
varyingMS$stat <- as.factor(varyingMS$mortalityshape)
varyingMS$cat <- 'mortalityshape'
varyingANPP <- bigTable[longevity == 200 & growthcurve %==% 0.5 & mortalityshape == 15,]
varyingANPP$stat <- as.factor(varyingANPP$mANPPproportion)
varyingANPP$cat <- 'anpp'

plotData <- rbind(varyingGCs, varyingMS, varyingANPP)
GCPlot <- ggplot(data = varyingGCs, aes(y = B, x = age)) + geom_line(size = 1, aes(col = stat)) + theme_bw() + 
  labs(name = "growth curve") + theme(legend.position = "bottom")
GCPlot

MSPlot <- ggplot(data = varyingMS, aes(y = B, x = age, col = stat)) + geom_line(size = 1) + theme_bw() + 
  labs(color = "mortality shape") + theme(legend.position = "bottom")
MSPlot

aNPPPlot <- ggplot(data = varyingANPP, aes(y = B, x = age, col = stat)) + geom_line(size = 1) + theme_bw() + 
  labs(color = "mANPP") + theme(legend.position = "bottom")
dev()
grid.arrange(GCPlot, MSPlot, aNPPPlot, ncol = 3)
