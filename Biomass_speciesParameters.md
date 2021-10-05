---
title: "Biomass_speciesParameters"
author: ""
date: "17 September 2019; updated Sept 30, 2021"
output:
  html_document:
    keep_md: yes
editor_options:
  chunk_output_type: console
---



[![Gitter](https://badges.gitter.im/PredictiveEcology/LandR_Biomass.svg)](https://gitter.im/PredictiveEcology/LandR_Biomass?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

# Overview

This module matches simulated growth curves of theoretical species with varying `maxANPP`, growth curve and mortality shape values (using LandR `Biomass_core`) against growth curves derived from Permanent Sample Plots (PSP data) across Western Canada to find the most likely combination of traits, and to adjust the `maxB` trait in respect to the observed maximum biomass on the landscape.

As of 2020-04-08, the PSP data needed for this module is not freely available, and data sharing agreements must be obtained from the governments of SK, AB, and BC. 

Each species grows as a single cohort with no understory (i.e., no dispersal, regeneration, or disturbance).
The full factorial included `mortalityshape` 5 to 25, in increments of 1; `growthcurve` 0 to 1, in increments of 0.1, 
`mANPPproportion` (the proportion of max ANPP to `maxB`) from .25 to 10 in increments of .25, and longevity from 150 to 700 in increments of 25. 
The results were saved in a combination of tables so that the module can be run without needing to simulate the factorial.

## Links with other modules

The module is intended to be used in combination with `Biomass_borealDataPrep.` 
However, users may want to experiment with the traits as well as examine the output GAMMs used to select traits. 
Therefore, it can be run as a stand-alone. The parameters are applied to all species (rather than each individually),
thus users may want to run several times, varying parameters, select the traits they feel best capture the true growth curve
(the PSPs are quite limited in their representation of stand dynamics for those stands older than 200 yrs), and modify the traits
directly using `Biomass_borealDataPrep.` 

# Usage


```r
library(data.table)
library(SpaDES.core)

setPaths(modulePath = file.path("../"))
getPaths() # shows where the 4 relevant paths are

times <- list(start = 0, end = 10)

modules <- list("Biomass_speciesParameters")

#the purpose of this table is experiment with modify longevity - longevity is not estimated by the module
#but it is used in trait estimation. 
inputSpecies <- data.table(species = c("Abie_bal", 'Abie_las', 'Betu_pap', 'Lari_lar',
                                        'Pice_eng', 'Pice_gla', 'Pice_mar', 'Pinu_ban',
                                       'Pinu_con', 'Pseu_men', "Popu_tre"),
                           longevity = c(300, 300, 170, 170, 330, 250, 250, 175, 300, 600, 200),
                           mortalityshape = 15, growthcurve = 0)
objects <- list(species = inputSpecies)

inputs <- list()
outputs <- list()

parameters <- list(
  Biomass_speciesParameters = 
    list(GAMMiterations = 2, 
         GAMMknots = list(
           "Abie_bal" = 3,
           "Abie_las" = 3,
           "Betu_pap" = 3,
           "Lari_lar" = 3,
           "Pice_eng" = 3,
           "Pice_gla" = 3,
           "Pice_mar" = 3,
           "Pinu_ban" = 3,
           "Pinu_con" = 3, 
           "Popu_tre" = 3,
           "Pseu_men" = 3),
         minimumPlotsPerGamm = 40, 
         constrainMortalityShape = c(8, 25), #also accepts a named list
         constrainGrowthCurve = c(0, 1), #also accepts a named list
         constrainMaxANPP = c(3.0, 4.0), #also accepts a named list
         quantileAgeSubset = list(
           "Abie_bal" = 98, 
           "Abie_las" = 98, 
           "Betu_pap" = 98, 
           "Lari_lar" = 98, 
           "Pice_eng" = 98, 
           "Pice_gla" = 98, 
           "Pice_mar" = 98, 
           "Pinu_ban" = 98, 
           "Pinu_con" = 98, 
           "Popu_tre" = 98, 
           "Pseu_men" = 95 #lower is better - very susceptible to outlying stand age
         )
    )
)

mySim <- simInit(times = times, params = parameters, modules = modules, objects = objects)
mySimOut <- spades(mySim)
```

# Events

Describe what happens for each event type.

## Plotting

Write what is plotted.

## Saving

Write what is saved.

# Downloads

During the `simInit` call, if the user does not provide alternatives for the expected inputs, the module will download 3 large `.tar` files (~2 GB each) and 1 `.zip` file (~45 MB) from the internet.

# Data dependencies

**NOTE:** all raster _inputs_ should be at the scale of `rasterToMatchLarge`/`studyAreaLarge` and all raster _outputs_ will be at the scale of `rasterToMatch`/`studyArea.`

## Module parameters


|paramName               |paramClass |default     |min |max |paramDesc                                                                                                                                                                                                                                                                                                                                        |
|:-----------------------|:----------|:-----------|:---|:---|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|.plotInitialTime        |numeric    |NA          |NA  |NA  |This describes the simulation time at which the first plot event should occur                                                                                                                                                                                                                                                                    |
|.plotInterval           |numeric    |NA          |NA  |NA  |This describes the simulation time interval between plot events                                                                                                                                                                                                                                                                                  |
|.saveInitialTime        |numeric    |NA          |NA  |NA  |This describes the simulation time at which the first save event should occur                                                                                                                                                                                                                                                                    |
|.saveInterval           |numeric    |NA          |NA  |NA  |This describes the simulation time interval between save events                                                                                                                                                                                                                                                                                  |
|.useCache               |logical    |FALSE       |NA  |NA  |Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant                                                                                                                                                                                         |
|biomassModel            |character  |Lambert2005 |NA  |NA  |The model used to calculate biomass from DBH. Can be either 'Lambert2005' or 'Ung2008'                                                                                                                                                                                                                                                           |
|constrainGrowthCurve    |numeric    |0, 1        |0   |1   |upper and lower bounds on range of potential growth curves when fitting traits. This module accepts a list of vectors, with names equal to `sppEquivCol`, so that traits are customizable                                                                                                                                                        |
|constrainMortalityShape |numeric    |12, 25      |5   |25  |upper and lower bounds on mortality shape when fitting traits. low mortality curve needs to excessive cohorts with very little biomass as longevity is approached, adding computation strain. alternatively accepts a list of vectors, with names equal to `sppEquivCol`.                                                                        |
|constrainMaxANPP        |numeric    |3, 4        |1   |10  |upper and lower bounds on maxANPP when fitting traits. Cohorts are initiated with `B = maxANPP`, which may be unreasonably high if mANPP is also high. Both `mANPP` and `growthcurve` params control when `maxB` is reached. High `mANPP` results in earlier peaks. Alternatively, accepts a list of vectors, with names equal to `sppEquivCol`. |
|GAMMiterations          |numeric    |8           |1   |NA  |number of iterations for GAMMs. This module accepts a list of vectors, with names equal to `sppEquivCol`, so that GAMMS are customizable                                                                                                                                                                                                         |
|GAMMknots               |numeric    |3           |NA  |NA  |the number of knots to use in the GAMM. Either 3 or 4 is recommended. This module accepts a list of vectors, with names equal to `sppEquivCol`, so that GAMMS are customizable                                                                                                                                                                   |
|minimumPlotsPerGamm     |numeric    |50          |10  |NA  |minimum number of PSP plots before building GAMM.                                                                                                                                                                                                                                                                                                |
|minDBH                  |integer    |0           |0   |NA  |minimum diameter at breast height (DBH) in cm used to filter PSP data. Defaults to 0 cm, i.e. all tree measurements are used.                                                                                                                                                                                                                    |
|PSPdataTypes            |character  |all         |NA  |NA  |Which PSP datasets to source, defaulting to all. Other available options include 'BC', 'AB', 'SK', 'NFI', and 'dummy'. 'dummy' should be used for unauthorized users.                                                                                                                                                                            |
|PSPperiod               |numeric    |1920, 2019  |NA  |NA  |The years by which to subset sample plot data, if desired. Must be a vector of length 2                                                                                                                                                                                                                                                          |
|quantileAgeSubset       |numeric    |98          |1   |100 |quantile by which to subset PSP data. As older stands are sparsely represented, the oldest measurements become vastly more influential. This parameter accepts both a single value and a list of vectors named by `sppEquivCol`. The PSP stand ages are found in `sim$speciesGAMMs$<species>$originalData`.                                      |
|sppEquivCol             |character  |default     |NA  |NA  |The column in `sim$sppEquiv` data.table to group species by. This parameter should share the same name as in Biomass_borealDataPrep. PSPs are aggregated by names in the PSP column and traits estimated for the corresponding names in the `sppEquivCol`                                                                                        |
|useHeight               |logical    |TRUE        |NA  |NA  |Should height be used to calculate biomass (in addition to DBH). DBH is used by itself when height is missing.                                                                                                                                                                                                                                   |

## Inputs

This module has several input requirements. 
One is a study area, which should be provided as a `SpatialPolygonsDataFrame`, and named `studyAreaLarge`.
This should be inside the boundaries of the boreal forest of Canada. 
When first running the code in this `.Rmd` file, you will be prompted to draw a polygon if none is provided as an input.


|objectName                 |objectClass              |desc                                                                                                                                                                                                                                                                                                          |sourceURL                                                                                    |
|:--------------------------|:------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------|
|factorialSpeciesTable      |data.table               |table with species traits for matching to factorialCohortData                                                                                                                                                                                                                                                 |https://drive.google.com/open?id=1q0ou0CBzD9GqGSparpHqf318IWK6ycty                           |
|reducedFactorialCohortData |data.table               |results of factorial species trait simulation. This can be found by running SpeciesFactorial.R but requires a specific commit of Boreal_Biomass                                                                                                                                                               |https://drive.google.com/open?id=1h8StXE0vm8xyDycRomCkwIaL7wfh5Irj                           |
|PSPmeasure_sppParams       |data.table               |merged PSP and TSP individual tree measurements. Must include the following columns: MeasureID, OrigPlotID1, MeasureYear, TreeNumber, Species, DBH and newSpeciesName the latter corresponding to species names in `LandR::sppEquivalencies_CA$PSP`. Defaults to randomized PSP data stripped of real plotIDs |https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing           |
|PSPplot_sppParams          |data.table               |merged PSP and TSP plot data. Defaults to randomized PSP data stripped of real plotIDs. Must contain fields MeasureID, MeasureYear, OrigPlotID1, and baseSA the latter being stand age at year of first measurement                                                                                           |https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing           |
|PSPgis_sppParams           |sf                       |Plot location sf object. Defaults to PSP data stripped of real plotIDs/location. Must include field OrigPlotID1 for joining to PSPplot_sppParams object                                                                                                                                                       |https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing           |
|species                    |data.table               |a table that has species traits such as longevity, shade tolerance, etc. Default is partially based on Dominic Cir and Yan's project                                                                                                                                                                          |https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv |
|speciesEcoregion           |data.table               |table defining the maxANPP, maxB and SEP, which can change with both ecoregion and simulation time. Defaults to a dummy table based on dummy data os biomass, age, ecoregion and land cover class                                                                                                             |NA                                                                                           |
|sppEquiv                   |data.table               |table of species equivalencies. See `LandR::sppEquivalencies_CA`.                                                                                                                                                                                                                                             |NA                                                                                           |
|studyAreaANPP              |SpatialPolygonsDataFrame |study area used to crop PSP data before building growth curves                                                                                                                                                                                                                                                |NA                                                                                           |

### Creates Inputs

Most of the inputs will be created automatically, if they are not provided by the user. 

## Outputs

This will show the outputs of this module, which can be used directly as the inputs for Biomass_core:


|objectName       |objectClass |desc                                                                                                                                                                                              |
|:----------------|:-----------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|speciesEcoregion |data.table  |table defining the maxANPP, maxB and SEP, which can change with both ecoregion and simulation time. Defaults to a dummy table based on dummy data os biomass, age, ecoregion and land cover class |
|speciesGAMMs     |list        |a list of mixed-effect general additive models (gamm) for each tree species modeling biomass as a function of age                                                                                 |
|species          |data.table  |a table that has species traits such as longevity...                                                                                                                                              |


```r
## species table
simOut$speciesTable
```


```r
Plot(simOut$biomassMap)
simOut$studyAreaLarge <- spTransform(simOut$studyAreaLarge, crs(simOut$biomassMap))
Plot(simOut$studyAreaLarge, addTo = "simOut$biomassMap")
```
