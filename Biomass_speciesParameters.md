---
title: "LandR _Biomass_speciesParameters_ Manual"
date: "Last updated: 2023-02-08"
output:
  bookdown::html_document2:
    toc: true
    toc_float: true
    toc_depth: 4
    theme: sandstone
    number_sections: false
    df_print: paged
    keep_md: yes
editor_options:
  chunk_output_type: console
  markdown: 
    wrap: 80
bibliography: citations/references_Biomass_speciesParameters.bib
citation-style: citations/ecology-letters.csl
link-citations: true
always_allow_html: true
---

<!-- the following are text references used in captions for LaTeX compatibility -->

(ref:Biomass-speciesParameters) *Biomass_speciesParameters*





[![module-version-Badge](/home/achubaty/Documents/GitHub/FOR-CAST/Ontario_AOU_ROF/modules/Biomass_speciesParameters/figures/moduleVersionBadge.png)](ssh://git@github.com/PredictiveEcology/Biomass_speciesParameters9362fa82d3f3b473ee441905dbee8ef2b68ef6a3)

[![Issues-badge](/home/achubaty/Documents/GitHub/FOR-CAST/Ontario_AOU_ROF/modules/Biomass_speciesParameters/figures/issuesBadge.png)](https://github.com/PredictiveEcology/Biomass_speciesParameters/issues)


<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:

Ian Eddy <ian.eddy@nrcan-rncan.gc.ca> [aut, cre], Eliot McIntire <eliot.mcintire@nrcan-rncan.gc.ca> [aut], Ceres Barros <ceres.barros@ubc.ca> [ctb]
<!-- ideally separate authors with new lines, '\n' not working -->

**This documentation is work in progress. Potential discrepancies and omissions
may exist for the time being. If you find any, contact us using the "Get help"
link above.**

## Module Overview

### Quick links

-   [General functioning](#bsppparam-general-functioning)

-   [List of input objects](#bsppparam-inputs-list)

-   [List of parameters](#bsppparam-params-list)

-   [List of outputs](#bsppparam-outputs-list)

-   [Simulation flow and module events](#bsppparam-sim-flow)

### Summary

LandR *Biomass_speciesParameters* (hereafter *Biomass_speciesParameters*)
calibrates species growth and mortality trait values used in *Biomass_core*, by
matching theoretical species' growth curves obtained with different trait values
(see [Simulated species data](#bsppparam-simdata)) against observed growth
curves derived from Permanent Sample Plots (PSP data) across Canada (see
[Permanent sample plot data](#bsppparam-PSPdata)), to find the combination of
trait values that allows a better match to the observed curves. In particular, 
it calibrates the `growthcurve`, `mortalityshape`, maximum biomass (`maxB`) and
maximum aboveground net primary productivity (`maxANPP`) traits (see [Parameter estimation/calibration](#bsppparam-calib)).

This module **will not** obtain other traits or parameters used in
*Biomass_core* and so it is meant to be used in conjunction with another data/calibration
module that does so (e.g., *Biomass_borealDataPrep*). However it can be used stand-alone 
in an initial developmental phase for easier inspection of the statistical 
calibration procedure employed.

As of February 08, 2023, the *raw* PSP data used in this module is not freely
available, and data sharing agreements must be obtained from the governments of
SK, AB, and BC to obtain it. However, the *processed and anonymized* PSP data is 
provided via a Google Drive folder accessed automatically by the module.

*Google Account is therefore necessary to access the data used for
calibration.**

### Links to other modules {#bsppparam-links-modules}

*Biomass_speciesParameters* is intended to be used with another data module,
like *Biomass_borealDataPrep*, that prepares all other traits and parameters
(including `maxB` and `maxANPP`) for *Biomass_core*. See
[here](https://rpubs.com/PredictiveEcology/LandR_Module_Ecosystem) for all
available modules in the LandR ecosystem and select *Biomass_speciesParameters*
from the drop-down menu to see potential linkages.

-   [*Biomass_borealDataPrep*](https://github.com/PredictiveEcology/Biomass_borealDataPrep):
prepares all parameters and inputs (including initial landscape conditions)
that *Biomass_core* needs to run a realistic simulation. Default
values/inputs produced are relevant for boreal forests of Western Canada.
Used upstream from *Biomass_speciesParameters*;

-   [*Biomass_core*](https://github.com/PredictiveEcology/Biomass_core): core
forest dynamics simulation module. Used downstream from
*Biomass_speciesParameters*;

-   [*Biomass_speciesFactorial*](https://github.com/PredictiveEcology/Biomass_core):
a module that generates theoretical species curves by running thousands of 
*Biomass_core* simulations on landscapes populated by one or more 
species, each simulation using a different set of species trait values.

-   [*Biomass_borealDataPrep*](https://github.com/PredictiveEcology/Biomass_borealDataPrep):
prepares all parameters and inputs (including initial landscape conditions)
that *Biomass_core* needs to run a realistic simulation. Default
values/inputs produced are relevant for boreal forests of Western Canada.
Used upstream from *Biomass_speciesParameters*;

-   [*Biomass_core*](https://github.com/PredictiveEcology/Biomass_core): core
forest dynamics simulation module. Used downstream from
*Biomass_speciesParameters*;

## Module manual

### General functioning {#bsppparam-general-functioning}

Tree cohort growth and mortality in *Biomass_core* are essentially determined by
five parameters: `growthcurve`, `mortalityshape`, maximum biomass (`maxB`), maximum 
aboveground net primary productivity (`maxANPP`) and `longevity`.

The `growthcurve` and `mortalityshape` parameters (called 'growth curve' and 
'mortality shape' in LANDIS-II Biomass Succession Extension v3.2, the base model
for *Biomass_core*) strongly modulate the shape of species growth curves and so 
it is important that they are calibrated to the study area in question.

Also, the growth and mortality equations used in *Biomass_core* are non-linear
and their resulting actual biomass accumulation curve is an emergent phenomenon
due to competition effects. This means that the ideal trait/parameter values
should not be estimated on pure single species growth conditions, as their
resulting dynamics will be different in a multi-species context.

*Biomass_speciesParameters* attempts to address these issues (at least partially)
using a "curve-matching" approach. It compares the best fit (according to their 
AIC) of three non-linear forms (Chapman-Richard's, Gompertz, and a logistic form) 
fitted to permanent sample plot (PSP) data to a large collection of theoretical 
(i.e. simulated) species curves, each representing a different set of the five 
key parameters that govern biomass increment in `Biomass_core`: `growthcurve`, 
`mortalityshape`, the ratio of `maxANPP` to `maxB`, and `longevity`. This library 
of curves is produced by the *Biomass_speciesFactorial* module.

*Biomass_speciesParameters* generally follows other LandR data modules, like 
*Biomass_boreaDataPrep*, which also attempts to calibrate previously estimated 
spatially varying species traits such as `maxB` and `maxANPP` from the input data
layers.

#### Permanent sample plot data {#bsppparam-PSPdata}

*Biomass_speciesParameters* can use all the PSP data available (note that it may
span several thousands of kilometres), or select the data based on a shapefile
(`studyAreaANPP`; see [List of input objects](#bsppparam-inputs-list)).

By default, the PSP data are obtained from the National Forest Inventory
(NFI), the Alberta Ministry of Agriculture, the Saskatchewan Ministry of the
Environment, and the British Columbia Ministry of Forests. These data were previously
treated for errors and standardized into a single dataset with the exact location and identifying
attributes anonymized.

The data include individual species, diameter at breast height (DBH), and
sometimes tree height measurements for each tree in a plot, as well as stand
age. As part of the standardization process, dead trees were removed from the
dataset. Tree biomass was then  per species using either a DBH-only model or a
DBH-height model from @LambertEtAl2005, in $g/m^2$.

Note that the model used to calculate biomass can also be changed to @UngEtAl2008
via the `P(sim)$biomassModel` module parameter (see [list of parameters](#bsppparam-params-list)).

#### Simulated species data {#bsppparam-simdata}

The *Biomass_speciesFactorial* module was used to create a library of
theoretical species curves (biomass accumulation curves, to be more precise) to 
which the best non-linear model form fit to the PSP-biomass will be matched for
each species and species combinations in the study area landscape. The library of curves was
created by running several *Biomass_core* simulations with no reproduction, competition,
disturbance, or dispersal effects, on the study area. Each simulation differed in
the combination of species trait values that influence growth and mortality
dynamics, namely: `growthcurve`, `mortalityshape`, `longevity`, `maxANPP` and
maximum biomass (`maxBiomass`, not to be confused with the data-driven `maxB`
which is later calibrated).

The values for `maxANPP` were explored via the `mANPPproportion`, the ratio of 
`maxANPP` to `maxBiomass` (the parameter used for theoretical curves), as it 
reflects their relationship.

`growthcurve` values varied from 0 to 1, in increments of 0.1; `mortalityshape`
varied from 5 to 25, in increments of 1; `longevity` varied from 150 to 700 in
increments of 25; `mANPPproportion` varied from 0.25 to 10 in increments of
0.25. `maxBiomass` was held constant at 5000.

This resulted in over 64,000,000 theoretical curves.

Results from these simulations were compiled into a table (`cohortDataFactorial`
; see [List of input objects](#bsppparam-inputs-list)) that is accessed by
*Biomass_speciesParameters*, so that the module can be run without needing to
re-simulate the theoretical curves.

#### Parameter estimation/calibration {#bsppparam-calib}

*Biomass_speciesParameters* calibrates `growthcurve`, `mortalityshape` and
`mANPPproportion` by matching the theoretical species curves produced by
*Biomass_speciesFactorial* (`cohortDataFactorial`) against observed
species growth curves from the PSP data.

Before calculating the *observed* species growth curves (i.e., the best of three
non-linear forms to match PSP data), the module subsets the PSP
data to stand ages below the 95th percent quantile for all species (this can be
changed via the `P(sim)$quantileAgeSubset` module parameter), as records for
larger age classes were limited and constituted statistical outliers. In some
species, changing the quantile value may improve results, however. Two examples
are *Pinus banksiana* and *Populus sp* (in western Canada), for which using the 
99th percent quantile improved the models, because these are short-lived species
for which data at advanced ages is scarce.

In addition, weights are added at the origin (age = 0 and biomass = 0) to force 
the intercept to be essentially at 0 age and 0 biomass.

The best fit of three non-linear forms, for each focal species, is then
calculated. Focal species are defined as either 50% of dominance in the plot, or
20% if we are looking to capture the multi-species dynamics (currently the
default). Three growth model forms are then fit to the observations for the
focal species: a Chapman-Richard's form [Equation \@ref(eq:Chapman); see, e.g.,
@CobleLee2006], a Gompertz form (Equation \@ref(eq:Gompertz)) and a Logistic
form [Equation \@ref(eq:Logistic); see @FekedulegnEtAl1999 for a complete
overview of these equations]. Multiple tries using the estimation methods from
the `robustbase::nlrob` function for each form are used, and the best model fit
is selected via Akaike Information Criterion (AIC).

```{=tex}
\begin{equation} 
  B \sim A \times (1 - e^{-k \times age})^{p}
  (\#eq:Chapman)
\end{equation}
```

```{=tex}
\begin{equation} 
  B \sim A \times e^{-k \times e^{-p \times age}}
  (\#eq:Gompertz)
\end{equation}
```
```{=tex}
\begin{equation} 
  B \sim \frac{A}{1 + k \times e^{-p \times age}} 
  (\#eq:Logistic)
\end{equation}
```

Species biomass ($B$) is estimated as a function of stand age ($age$), with the
best values of the $A$, $k$ and $p$ parameters to fit the PSP data.

It is possible that some selected species do not have enough data to allow for
model convergence. In this case, *Biomass_speciesParameters* skips trait
(re-)calibration, and values remain unchanged.

After each species best fit is selected (using AIC), *Biomass_speciesParameters*
compares it to the library of theoretical curves, and picks the best one based
on maximum likelihood. This best theoretical curve will be associated with a
given combination of `growthcurve`, `mortalityshape` and `maxANPPproportion`
values, which are then used directly as the calibrated values, in case of
`growthcurve` and `mortalityshape`, or to calibrate `maxANPP` in the case of
`maxANPPproportion` (see below).

Since simulated growth curves never achieve the maximum biomass parameter (the
`maxBiomass` parameter set to 5000 for all simulations of theoretical species
curves, or the `maxB` parameter in *Biomass_core* simulations), it acts as an
asymptotic limit that reflects the potential maximum biomass for a species in an
ecolocation (ecological zone and land cover combination).

*Biomass_speciesParameters* uses the ratio between the potential maximum biomass
(`maxBiomass`, always 5000) to the achieved maximum biomass in the theoretical
curves, to rescale `maxB`. This ratio is called the `inflationFactor` and it is
multiplied by `maxB` values previously estimated from data (e.g. by
*Biomass_borealDataPrep*). This way, species simulated in *Biomass_core* are
able to achieve the maximum observed biomasses used to *initially* estimate
`maxB`.

Finally, the module calibrates `maxANPP` using the `mANPPproportion` value from
the best matching theoretical growth curve as:

```{=tex}
\begin{equation}
maxB \times \frac{mANPPproportion}{100}
(\#eq:maxANPPcalib)
\end{equation}
```
where `maxB` is the already (re-)calibrated version.

In cases where there are not sufficient PSP data to fit the growth models and 
perform the calibration, `mANPPproportion` defaults to 3.33 (the value used in LANDIS-II
applications in Canada's boreal forests) and the `inflationFactor` defaults to 1.

### List of input objects {#bsppparam-inputs-list}

The full list of input objects required by the module is presented below (Table
\@ref(tab:moduleInputs2-Biomass-speciesParameters)). The only input that
**must** be provided is `studyAreaANPP` (the study area used extract the PSP
data from). All other input objects have internal defaults, but the user may
need to request access to their online files.

Of these inputs, the following are particularly important and deserve special
attention:

**Spatial layers**

-   `studyAreaANPP` -- shapefile. A `SpatialPolygonsDataFrame` with a single
polygon determining the where the PSP should be subset to simulation will
take place. This input object **must be supplied by the user or another module**.

**Tables**

-   `speciesTableFactorial` and `cohortDataFactorial` -- a tables of species
    trait combinations and the theoretical species grwoth curve data
    (respectively)
-   `PSPmeasure_sppParams`, `PSPplot_sppParams` and `PSPgis_sppParams` --
    tree measurement, biomass growth and geographical data of the PSP
    datasets used to buildi observed species growth curves.
-   `species` -- a table of invariant species traits that may have been
    produced by another module. It **must** contain the columns 'species',
    'growthcurve' and 'mortality shape', whose values will be calibrated.
-   `speciesEcoregion` -- table of spatially-varying species traits that may
    have been produced by another module. It **must** contain the columns
    'speciesCode', 'maxB' and 'maxANPP' and 'ecoregionGroup' (the
    ecolocation ID). 'maxB' and 'maxANPP' values will be calibrated by
    species.

\newpage
\blandscape

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs2-Biomass-speciesParameters)List of (ref:Biomass-speciesParameters) input objects and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
   <th style="text-align:left;"> sourceURL </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> cohortDataFactorial </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Results of factorial species trait simulation. </td>
   <td style="text-align:left;"> https://drive.google.com/file/d/1NH7OpAnWtLyO8JVnhwdMJakOyapBnuBH/ </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PSPmeasure_sppParams </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Merged PSP and TSP individual tree measurements. Must include the following columns: 'MeasureID', 'OrigPlotID1', 'MeasureYear', 'TreeNumber', 'Species', 'DBH' and 'newSpeciesName' the latter corresponding to species names in `LandR::sppEquivalencies_CA$PSP`. Defaults to randomized PSP data stripped of real plotIDs </td>
   <td style="text-align:left;"> https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PSPplot_sppParams </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Merged PSP and TSP plot data. Defaults to randomized PSP data stripped of real plotIDs. Must contain columns 'MeasureID', 'MeasureYear', 'OrigPlotID1', and 'baseSA', the latter being stand age at year of first measurement </td>
   <td style="text-align:left;"> https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PSPgis_sppParams </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Plot location `sf` object. Defaults to PSP data stripped of real plotIDs/location. Must include field 'OrigPlotID1' for joining to PSPplot object </td>
   <td style="text-align:left;"> https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing </td>
  </tr>
  <tr>
   <td style="text-align:left;"> species </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> A table of invariant species traits with the following trait colums: 'species', 'Area', 'longevity', 'sexualmature', 'shadetolerance', 'firetolerance', 'seeddistance_eff', 'seeddistance_max', 'resproutprob', 'mortalityshape', 'growthcurve', 'resproutage_min', 'resproutage_max', 'postfireregen', 'wooddecayrate', 'leaflongevity' 'leafLignin', 'hardsoft'. Only 'growthcurve' and 'mortalityshape' are used in this module. Default is from Dominic Cyr and Yan Boulanger's applications of LANDIS-II </td>
   <td style="text-align:left;"> https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv </td>
  </tr>
  <tr>
   <td style="text-align:left;"> speciesEcoregion </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Table of spatially-varying species traits ('maxB', 'maxANPP', 'establishprob'), defined by species and 'ecoregionGroup') Defaults to a dummy table based on dummy data os biomass, age, ecoregion and land cover class </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppEquiv </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Table of species equivalencies. See `?LandR::sppEquivalencies_CA`. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> speciesTableFactorial </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Table with species traits to be matched with `sim$cohortDataFactorial`. </td>
   <td style="text-align:left;"> https://drive.google.com/file/d/1NH7OpAnWtLyO8JVnhwdMJakOyapBnuBH/ </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyAreaANPP </td>
   <td style="text-align:left;"> SpatialPolygonsDataFrame </td>
   <td style="text-align:left;"> Study area used to crop PSP data before building growth curves </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

\elandscape

### List of parameters {#bsppparam-params-list}

The full list of parameters used by the module is presented below (Table
\@ref(tab:moduleParams2-Biomass-speciesParameters)), all of which have default
values specified in the module's metadata.

Of these parameters, the following are particularly important:

**Calibration parameters**

-   `biomassModel` -- the model used to calculate biomass from DBH

-   `speciesFittingApproach` -- should the calibration take into account species
growing in single- or multi-species context?

**Data processing**

-   `PSPperiod` -- PSP data period to use.

-   `quantileAgeSubset` -- upper quantile age value used to subset PSP data.

\newpage
\blandscape

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleParams2-Biomass-speciesParameters)List of (ref:Biomass-speciesParameters) parameters and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> paramName </th>
   <th style="text-align:left;"> paramClass </th>
   <th style="text-align:left;"> default </th>
   <th style="text-align:left;"> min </th>
   <th style="text-align:left;"> max </th>
   <th style="text-align:left;"> paramDesc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> biomassModel </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> Lambert2005 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> The model used to calculate biomass from DBH. Can be either 'Lambert2005' or 'Ung2008'. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GAMMiterations </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Number of iterations for GAMMs. Accepts a list of vectors, with names equal to those in `sim$sppEquiv[, P(sim)$sppEquivCol]`, so that GAMMS are customizable per species </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GAMMknots </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> The number of knots to use in the GAMM. Either 3 or 4 is recommended. Accepts a list of vectors, with names equal to those in `sim$sppEquiv[, P(sim)$sppEquivCol]`, so that GAMMS are customizable per species </td>
  </tr>
  <tr>
   <td style="text-align:left;"> maxBInFactorial </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 5000 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> The arbitrary maximum biomass for the factorial simulations. This is a per-species maximum within a pixel </td>
  </tr>
  <tr>
   <td style="text-align:left;"> minimumPlotsPerGamm </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 50 </td>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Minimum number of PSP plots before building GAMM </td>
  </tr>
  <tr>
   <td style="text-align:left;"> minDBH </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Minimum diameter at breast height (DBH) in cm used to filter PSP data. Defaults to 0 cm, i.e. all tree measurements are used. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PSPdataTypes </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> all </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Which PSP datasets to source, defaulting to all. Other available options include 'BC', 'AB', 'SK', 'NFI', and 'dummy'. 'dummy' should be used for unauthorized users. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PSPperiod </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1920, 2019 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> The years by which to subset sample plot data, if desired. Must be a vector of length 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> quantileAgeSubset </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 95 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 100 </td>
   <td style="text-align:left;"> Quantile by which to subset PSP data. As older stands are sparsely represented, the oldest measurements become vastly more influential. This parameter accepts both a single value and a list of vectors named by `sim$sppEquiv[, P(sim)$sppEquivCol]`. The PSP stand ages are found in `sim$speciesGAMMs$SPECIES$originalData`, where SPECIES is the species ID </td>
  </tr>
  <tr>
   <td style="text-align:left;"> speciesFittingApproach </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> focal </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Either 'all', 'pairwise', 'focal' or 'single', indicating whether to pool all species into one fit, do pairwise species (for multiple cohort situations), do pairwise species, but using a focal species approach where all other species are pooled into 'other' or do one species at a time. If 'all', all species will have identical species-level traits </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppEquivCol </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> default </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> The column in `sim$sppEquiv` data.table to group species by. This parameter should share the same name as in Biomass_borealDataPrep . PSPs are aggregated by names in the PSP column and traits estimated for species with corresponding names in the `sim$sppEquiv[, P(sim)$sppEquivCol]` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> standAgesForFitting </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 0, 150 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> The minimum and maximum ages to use while matching NonLinearFit (or GAMM) with LandR curves provided in the factorial. Since the majory of the data that went into fits for the NonLinearFit from PSPs is less than 200, it is likely wise to constrain the range to something smaller than 0 to 200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> useHeight </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Should height be used to calculate biomass (in addition to DBH). DBH is used by itself when height is missing. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plots </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> screen </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Used by Plots function, which can be optionally used here </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This describes the simulation time at which the first plot event should occur </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This describes the simulation time interval between plot events </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This describes the simulation time at which the first save event should occur </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This describes the simulation time interval between save events </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> .inputOb.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant </td>
  </tr>
</tbody>
</table>

\elandscape

### List of outputs {#bsppparam-outputs-list}

The module produces the following outputs (Table
\@ref(tab:moduleOutputs-Biomass-speciesParameters)). Note that `species` and
`speciesEcoregion` are modified versions of the input objects with the same
name.

**Tables**

-   `species` and `speciesEcoregion` -- tables with calibrated trait values.

-   `speciesGAMMs` -- the fitted GAMM model objects for each species.

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleOutputs-Biomass-speciesParameters)List of (ref:Biomass-speciesParameters) output objects and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> species </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> The updated invariant species traits table (see description for this object in inputs) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> speciesEcoregion </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> The updated spatially-varying species traits table (see description for this object in inputs) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> speciesGAMMs </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> A list of mixed-effect general additive models (GAMMs) for each tree species modeling biomass as a function of age </td>
  </tr>
</tbody>
</table>

### Simulation flow and module events {#bsppparam-sim-flow}

*Biomass_speciesParameters* initializes itself and prepares all inputs provided
there is an active internet connection and the user has access to the data (and
a Google Account to do so).

We advise future users to run *Biomass_speciesParameters* with defaults and
inspect what the objects are like before supplying their own data. The user does
not need to run *Biomass_speciesFactorial* to generate their own theoretical
curves (unless they wish to), as the module accesses pre-generated theoretical curves.

Note that this module only runs once (in one "time step") and only executes one
event (`init`). The general flow of *Biomass_speciesParameters* processes is:

1.  Preparation of all necessary data and input objects that do not require
parameter fitting (e.g., the theoretical species growth curve data);

2.  Sub-setting PSP data and calculating the observed species growth curves
    using non-linear growth models;

3.  Finding the theoretical species growth curve that best matches the observed
curve, for each species. Theoretical curves are subset to those with longevity
matching the species' longevity (in `species` table) and with
`growthcurve` and `mortalityshape` values;

4.  Calibrating `maxB` and `maxANPP`.

## Usage example {#bsppparam-example}

This module can be run stand-alone, but it won't do much more than calibrate
species trait values based on dummy input trait values. We provide an example of
this below, since it may be of value to run the module by itself to become
acquainted with the calibration process and explore the fitted non-linear
models. However, we remind that to run this example you will need a Google
Account, and to be granted access to the data.

A realistic usage example of this module and a few others can be found in [this
repository](https://github.com/CeresBarros/LandRBiomass_publication) and in
@BarrosEtAlinreview.

### Load `SpaDES` and other packages.

### Set up R libraries {#bsppparam-example-libs}


```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
tempDir <- tempdir()

pkgPath <- file.path(tempDir, "packages", version$platform,
                     paste0(version$major, ".", strsplit(version$minor, "[.]")[[1]][1]))
dir.create(pkgPath, recursive = TRUE)
.libPaths(pkgPath, include.site = FALSE)

if (!require(Require, lib.loc = pkgPath)) {
  remotes::install_github(
    paste0("PredictiveEcology/",
           "Require@5c44205bf407f613f53546be652a438ef1248147"),
    upgrade = FALSE, force = TRUE)
  library(Require, lib.loc = pkgPath)
}

setLinuxBinaryRepo()
```

### Get the module and module dependencies {#bsppparam-example-pkg-mods}


```r
Require(pasteo("PredictiveEcology/",
               "SpaDES.project@6d7de6ee12fc967c7c60de44f1aa3b04e6eeb5db"), 
        require = FALSE, upgrade = FALSE, standAlone = TRUE)

paths <- list(inputPath = normPath(file.path(tempDir, "inputs")), 
              cachePath = normPath(file.path(tempDir, "cache")), 
              modulePath = normPath(file.path(tempDir, "modules")), 
              outputPath = normPath(file.path(tempDir, "outputs")))

SpaDES.project::getModule(modulePath = paths$modulePath,
                          c("PredictiveEcology/Biomass_speciesParameters"),
                          overwrite = TRUE)

## make sure all necessary packages are installed:
outs <- SpaDES.project::packagesInModules(modulePath = paths$modulePath)
Require(c(unname(unlist(outs)), "SpaDES"),
        require = FALSE, standAlone = TRUE)

## load necessary packages
Require(c("SpaDES"), upgrade = FALSE, install = FALSE)
```

### Setup simulation


```r
times <- list(start = 0, end = 1)

modules <- list("Biomass_speciesParameters")

#the purpose of this table is experiment with modify longevity - longevity is not estimated by the module
#but it is used in trait estimation. 
inputSpecies <- data.table(species = c("Abie_bal", 'Abie_las', 'Betu_pap', 
                                       'Lari_lar', 'Pice_eng', 'Pice_gla', 
                                       'Pice_mar', 'Pinu_ban', 'Pinu_con', 
                                       'Pseu_men', "Popu_tre"),
                           longevity = c(300, 300, 170, 170, 330, 250, 
                                         250, 175, 300, 600, 200),
                           mortalityshape = 15, growthcurve = 0)
objects <- list(species = inputSpecies)

inputs <- list()
outputs <- list()

parameters <- list(Biomass_speciesParameters = 
                     list(quantileAgeSubset = list(
                            "Abie_bal" = 95, 
                            "Abie_las" = 95,
                            "Betu_pap" = 95,
                            "Lari_lar" = 95,
                            "Pice_eng" = 95,
                            "Pice_gla" = 95,
                            "Pice_mar" = 95,
                            "Pinu_ban" = 95,
                            "Pinu_con" = 99, 
                            "Popu_tre" = 99,
                            "Pseu_men" = 99
                          )
                     ))

mySim <- simInitAndSpades(times = times, 
                          params = parameters, 
                          modules = modules, 
                          paths = paths, 
                          objects = objects)

## to inspect the fitted GAMM models:
mySim$speciesGAMMs$Pice_mar
```

## References {#bsppparam-refs}
