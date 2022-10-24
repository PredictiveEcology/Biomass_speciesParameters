---
title: "LandR _Biomass_speciesParameters_ Manual"
date: "Last updated: 2022-10-24"
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





[![module-version-Badge](D:/GitHub/LandR-Manual/modules/Biomass_speciesParameters/figures/moduleVersionBadge.png)](https://github.com/CeresBarros/Biomass_speciesParameters/commit/6500980f1c38286d0c677ba608151b8521dc3e41)

[![Issues-badge](D:/GitHub/LandR-Manual/modules/Biomass_speciesParameters/figures/issuesBadge.png)](https://github.com/PredictiveEcology/Biomass_speciesParameters/issues)


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
trait values that allows a better match to the observed curves.

In particular, it directly calibrates the `growthcurve`, `mortalityshape`
invariant species traits and two new traits `inflationFactor` and
`mANPPproportion`, which are used to calibrate previously estimated species
maximum biomass (`maxB`) and maximum aboveground net primary productivity
(`maxANPP`) values (see [Parameter estimation/calibration](#bsppparam-calib)).

This module **will not** obtain other traits or parameters used in
*Biomass_core* and so must be used in conjunction with another data/calibration
module that does so (e.g., *Biomass_borealDataPrep*).

It can however be used stand-alone in an initial developmental phase for easier
inspection of the statistical calibration procedure employed.

As of October 24, 2022, the *raw* PSP data used in this module is not freely
available, and data sharing agreements must be obtained from the governments of
SK, AB, and BC to obtain it. However, the *processed and anonymized* PSP data is 
provided via a Google Drive folder accessed automatically by the module.

**A Google Account is therefore necessary to access the data used for
calibration.**

If you do not have a Google Account, or cannot access the data, please report an
issue by clicking on the "Get help" link above.

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

## Module manual

### General functioning {#bsppparam-general-functioning}

Tree cohort growth and mortality in *Biomass_core* are essentially determined by
five parameters: the invariant species traits 'growth curve' (`growthcurve`),
'mortality shape', (`mortalityshape`) and `longevity`, and the spatio-temporally
varying traits maximum biomass (`maxB`) and maximum aboveground net primary
productivity (`maxANPP`).

All five traits strongly modulate the shape of species growth curves and so it
is important that they are calibrated to the study area in question.

Also, the growth and mortality equations used in *Biomass_core* are non-linear
and their resulting actual biomass accumulation curve is an emergent phenomenon
due to competition effects. This means that the ideal trait/parameter values
should not be estimated on pure single species growth conditions, as their
resulting dynamics will be different in a multi-species context.

*Biomass_speciesParameters* attempts to address these issues (at
least partially) using a "curve-matching" approach. It compares a GAMM fitted to
permanent sample plot (PSP) data to a large collection of theoretical species
curves, each representing a different set of growth and mortality parameters.
This also provides a means to calibrate these traits using a dataset that is
independent from the one used to derive initial landscape conditions and initial
values of `maxB` and `maxANPP`.

While `longevity` is adjusted using published values (see
*Biomass_borealDataPrep* manual), the remaining four parameters are calibrated
using the PSP data. Hence, *Biomass_speciesParameters* generally follows other
data modules, like *Biomass_boreaDataPrep*, that prepare other traits such as
`longevity`, `maxB` and `maxANPP`.

#### Permanent sample plot data {#bsppparam-PSPdata}

*Biomass_speciesParameters* can use all the PSP data available (note that it may
span several thousands of kilometres), or select the data based on a polygon
(`studyAreaANPP`; see [List of input objects](#bsppparam-inputs-list)).

The default PSP data were initially obtained from the National Forest Inventory
(NFI), the Alberta Ministry of Agriculture, the Saskatchewan Ministry of the
Environment, and the British Columbia Ministry of Forests, treated for errors
and standardized into a single data set with the exact location and identifying
attributes anonymized. We only share the randomized and anonymized dataset, as data sharing
agreements must be met to access the raw data.

The data include individual species, diameter at breast height (DBH), and
sometimes tree height measurements for each tree in a plot, as well as stand
age. As part of the standardization process, dead trees were removed from the
dataset. Tree biomass was then estimated by tree species, in $g/m^2$, using either the DBH-only model or a
DBH-height model from either @LambertEtAl2005 or @UngEtAl2008 (see `P(sim)$biomassModel` module parameter
in [list of parameters](#bsppparam-params-list)).

#### Simulated species data {#bsppparam-simdata}

The *Biomass_speciesFactorial* module was used to create a library of
theoretical species curves (biomass accumulation curves, to be more precise) to
which the empirical species curves derived from PSP-biomass are matched for each
species trait combination in the study area. The library of curves was
created by running several *Biomass_core* simulations with no reproduction, competition,
disturbance, or dispersal effects, on the study area. Each simulation differed in
the combination of species trait values that influence growth and mortality
dynamics, namely: `growthcurve`, `mortalityshape`, `longevity`, `maxANPP` and
maximum biomass (`maxBiomass`, not to be confused with the data-driven `maxB`
which is later calibrated).

The values for `maxANPP` were explored via the `mANPPproportion`, the ratio of 
`maxANPP` to `maxBiomass` (the parameter used for theroetical curves), as it 
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
*Biomass_speciesFactorial* (`cohortDataFactorial` object) against observed
species growth curves from permanent sample plot (PSP) data.

Before fitting the *observed* species growth curves, the module subsets the PSP
data to stand ages below the 95th percent quantile for all species (this can be
changed via the `P(sim)$quantileAgeSubset` module parameter), as records for
larger age classes were limited and constituted statistical outliers. In some
species, changing the quantile value may improve results, however. Two examples
are *Pinus banksiana* and *Populus sp*, for which using the 99th percent
quantile improved the models, because these are short-lived species for which
data at advanced ages is scarce.

The module attempts to fit the models using stands where the focal species is 
dominant (but not monocultures), while balancing sample size (see [biomass
weighting](#bsppparam-calibbw) below). Hence, for a given species, it only includes plots where the
species' relative biomass is at least 50%. This is, when calibrating *Populus
tremuloides* traits, PSP daa plots are only included if 50% of the stand biomass
is composed of *P. tremuloides*.

In addition, 50 points are added at the origin (age = 0 and biomass = 0) to
force the intercept to be essentially 0 age and 0 biomass.

Observed growth curves for each species are then fit using generalized additive
mixed models (GAMMs) that relate species biomass ($B$) with stand age
($standAge$), accounting for the random effects of the measurement year
($measureYear$) and plot ($plotID$) on the intercept:

```{=tex}
\begin{equation}
B \sim f_{1}(standAge) + (\sim 1 | measureYear + plotID)
(\#eq:GAMM)
\end{equation}
```
where $f_{1}$ denotes the smoother function. To avoid overfitting, the module
constrains the smoother on stand age to a maximum smoothing degree of 3 (i.e. 3
knots and a polynomial degree of 2) and a default point constraint at 0 that attempts to
force the intercept to 0. The smoother degree constraint, however,
can be changed via the `P(sim)$GAMMknots` module parameter.

##### Biomass-weighting {#bsppparam-calibbw}
In addition, $B$ is weighted with respect to species dominance. This consisted
in 1) calculating the average biomass of each dominant species (i.e. relative
biomass in a plot \> 0.5; $domSpeciesB_{1}$), in each plot and measurement year,
and 2) dividing the species average biomass by the average biomass across all
*n* dominant species ($allDomSpeciesB$):

```{=tex}
\begin{equation}
\frac{\overline{\rm domSpeciesB_{1}}}{\overline{\rm allDomSpeciesB}}
(\#eq:Bweights)
\end{equation}
```
For the added 0 age and 0 biomass data the module uses weights equal to 1.

It is possible that some selected species do not have enough data to allow for
model convergence. In this case, *Biomass_speciesParameters* skips trait
(re-)calibration, and values remain unchanged.

After fitting each species GAMM, *Biomass_speciesParameters* compares it to the
theoretical curves obtained with a `longevity` value that matches the focal
species' longevity, and picks the best one based on maximum likelihood. This best
theoretical curve will be associated with a given combination of `growthcurve`,
`mortalityshape` and `mANPPproportion` values, which are then used directly as
the calibrated values, in case of `growthcurve` and `mortalityshape`, or to
calibrate `maxANPP` in the case of `mANPPproportion` (see below).

The user has the option to constrain the values of the `growthcurve` and
`mortalityshape` parameters. By default, `growthcurve` is forced to 0.5,
`mortalityshape` is allowed to vary between 15 and 25, and `mANPPproportion`
between 2.0 and 5.0 (see module parameters `P(sim)$constrainGrowthCurve`,
`P(sim)constrainMortalityShape` and `P(sim)constrainMaxANPP`). These boundary
values were based on preliminary runs and analyses using the default data and
may not apply to other data sets, or to different spatial subsets of the default
data.

If boundary values are used, *Biomass_speciesParameters* subsets the theoretical
species growth curves to those with trait values within the selected boundaries.

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
where `maxB` is the already (re-)calibrated version. As already stated above, the
final `maxANPP` value is then constrained between 2.0 and 5.0 by default.

In cases where insufficient PSP data prevent fitting the GAMMs and performing
the calibration, `mANPPproportion` defaults to 3.33 (the value used in LANDIS-II
applications in Canada's boreal forests) and the `inflationFactor` to 1.

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

-   `factorialSpeciesTable` and `reducedFactorialCohortData` -- a tables of
species trait combinations and the theoretical species grwoth curve data
(respectively).
-   `PSPmeasure`, `PSPplot` and `PSPgis` -- tree measurement, biomass growth and
geographical data of the PSP datasets used to build observed species growth
curves.
-   `species` -- a table of invariant species traits that may have been produced
by another module. It **must** contain the columns 'species', 'growthcurve'
and 'mortality shape', whose values will be calibrated.
-   `speciesEcoregion` -- table of spatially-varying species traits that may
have been produced by another module. It **must** contain the columns
'speciesCode', 'maxB' and 'maxANPP' and 'ecoregionGroup' (the ecolocation
ID). 'maxB' and 'maxANPP' values are (re-)calibrated by species.

\newpage
\blandscape

<table>
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
   <td style="text-align:left;"> factorialSpeciesTable </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Table with species traits to be matched with `sim$reducedFactorialCohortData` </td>
   <td style="text-align:left;"> https://drive.google.com/open?id=1q0ou0CBzD9GqGSparpHqf318IWK6ycty </td>
  </tr>
  <tr>
   <td style="text-align:left;"> reducedFactorialCohortData </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Results of factorial species trait simulation. This can be found by running `SpeciesFactorial.R` but requires a specific commit of Biomass_core </td>
   <td style="text-align:left;"> https://drive.google.com/open?id=1h8StXE0vm8xyDycRomCkwIaL7wfh5Irj </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PSPmeasure </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Merged PSP and TSP individual tree measurements. Must include the following columns: 'MeasureID', 'OrigPlotID1', 'MeasureYear', 'TreeNumber', 'Species', 'DBH' and 'newSpeciesName' the latter corresponding to species names in `LandR::sppEquivalencies_CA$PSP`. Defaults to randomized PSP data stripped of real plotIDs </td>
   <td style="text-align:left;"> https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PSPplot </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Merged PSP and TSP plot data. Defaults to randomized PSP data stripped of real plotIDs. Must contain columns 'MeasureID', 'MeasureYear', 'OrigPlotID1', and 'baseSA', the latter being stand age at year of first measurement </td>
   <td style="text-align:left;"> https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PSPgis </td>
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

-   `GAMMiterations` and `GAMMknots` -- control the number of iterations and
smoother degree used to fit the GAMMs, respectively.

-   `constrainGrowthCurve`, `constrainMortalityShape` and `constrainMaxANPP` --
determine the upper and lower boundaries of the calibrated values of
`growthcurve`, `mortalityshape` and `maxANPP`, respectively.

**Data processing**

-   `minimumPlotsPerGamm` -- define a minimum number of PSP plots needed to fit
the GAMMs.

-   `PSPperiod` -- PSP data period to use.

-   `quantileAgeSubset` -- upper quantile age value used to subset PSP data.

\newpage
\blandscape


<table>
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
   <td style="text-align:left;"> The model used to calculate biomass from DBH. Can be either 'Lambert2005' or 'Ung2008' </td>
  </tr>
  <tr>
   <td style="text-align:left;"> constrainGrowthCurve </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0.5, 0.5 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> Upper and lower bounds on range of potential growth curves when fitting traits. This module accepts a list of vectors, with names equal to `P(sim)$sppEquivCol`, so that traits are customizable </td>
  </tr>
  <tr>
   <td style="text-align:left;"> constrainMortalityShape </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 15, 25 </td>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> 25 </td>
   <td style="text-align:left;"> Upper and lower bounds on mortality shape when fitting traits. Low mortality curve values lead to numerous cohorts with very little biomass as longevity is approached, adding computation strain. Alternatively accepts a list of vectors, with names equal to those in `sim$sppEquiv[, P(sim)$sppEquivCol]` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> constrainMaxANPP </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 2, 5 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> Upper and lower bounds on `maxANPP` when fitting traits. When cohorts are initiated with `B = maxANPP` (see Biomass_core parameter `P(sim)$initialB`), this can lead to unreasonably initial biomass if `maxANPP` is also high. Both `maxANPP` and `growthcurve` parameters control when `maxB` is reached. High `maxANPP` results in earlier peaks. Alternatively accepts a list of vectors, with names equal to those in `sim$sppEquiv[, P(sim)$sppEquivCol]` </td>
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
   <td style="text-align:left;"> minDBH </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Minimum diameter at breast height (DBH) in cm used to filter PSP data. Defaults to 0cm, i.e. all tree measurements are used. </td>
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
   <td style="text-align:left;"> sppEquivCol </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> default </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> The column in `sim$sppEquiv` data.table to group species by. This parameter should share the same name as in Biomass_borealDataPrep . PSPs are aggregated by names in the PSP column and traits estimated for species with corresponding names in the `sim$sppEquiv[, P(sim)$sppEquivCol]` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> useHeight </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Should height be used to calculate biomass (in addition to DBH)? We advise against including height unless you are certain it is present in every PSP </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
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
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
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
`speciesEcoregion` are modified versions of the inputed objects with the same
name.

**Tables**

-   `species` and `speciesEcoregion` -- tables with calibrated trait values.

-   `speciesGAMMs` -- the fitted GAMM model objects for each species.

<table>
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

*Biomass_speciesParameters* initialies itself and prepares all inputs provided
there is an active internet connection and the user has access to the data (and
a Google Account to do so).

We advise future users to run *Biomass_speciesParameters* with defaults and
inspect what the objects are like before supplying their own data. The user does
not need to run *Biomass_speciesFactorial* to generate their own theoretical
curves (unless they wish to), as the module accesses a pre-generated on-line
library with these simulated data.

Note that this module only runs once (in one "time step") and only executes one
event (`init`). The general flow of *Biomass_speciesParameters* processes is:

1.  Preparation of all necessary data and input objects that do not require
parameter fitting (e.g., the theoretical species growth curve data);

2.  Sub-setting PSP data and calculating the observed species growth curves
using GAMMs;

3.  Finding the theoretical species growth curve that best matches the observed
curve, for each species. Theoretical curves are subset to those with longevity
matching the species' longevity (in `species` table) and with
`growthcurve` and `mortalityshape` values within the chosen boundaries
(`P(sim)$constrainGrowthCurve`, `P(sim)$constrainMortalityShape`);

4.  Calibrating `maxB` and `maxANPP`

5.  Adjusting `maxANPP` to match the chosen boundaries
(`P(sim)$constrainMaxANPP`)

## Usage example {#bsppparam-example}

This module can be run stand-alone, but it won't do much more than calibrate
species trait values based on dummy input trait values. We provide an example of
this below, since it may be of value to run the module by itself to become
acquainted with the calibration process and explore the fitted GAMMs. However,
we remind that to run this example you will need a Google Account, and access to
the data may need to be granted.

A realistic usage example of this module and a few others can be found in [this
repository](https://github.com/CeresBarros/LandRBiomass_publication) and in
@BarrosEtAlinreview.

### Load `SpaDES` and other packages.


```r
if (!require(Require)) {
  install.packages("Require")
  library(Require)
}

Require(c("PredictiveEcology/SpaDES.install",
          "SpaDES", "PredictiveEcology/SpaDES.core@development"), 
        install_githubArgs = list(dependencies = TRUE))
```

### Get module, necessary packages and set up folder directories


```r
tempDir <- tempdir()

paths <- list(inputPath = normPath(file.path(tempDir, "inputs")), 
              cachePath = normPath(file.path(tempDir, "cache")), 
              modulePath = normPath(file.path(tempDir, "modules")), 
              outputPath = normPath(file.path(tempDir, "outputs")))

getModule("PredictiveEcology/Biomass_speciesParameters@79896a4e3b785e34e5f509798ab6c2204bb334d7", modulePath = paths$modulePath, overwrite = TRUE)

## make sure all necessary packages are installed:
makeSureAllPackagesInstalled(paths$modulePath)
```

### Setup simulation


```r
library(SpaDES)

times <- list(start = 0, end = 1)

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

parameters <- list(Biomass_speciesParameters = 
                     list(GAMMiterations = 2, 
                          GAMMknots = list(
                            "Abie_bal" = 3,
                            "Abie_las" = 3,
                            "Betu_pap" = 3,
                            "Lari_lar" = 4,
                            "Pice_eng" = 4,
                            "Pice_gla" = 3,
                            "Pice_mar" = 4,
                            "Pinu_ban" = 3,
                            "Pinu_con" = 4, 
                            "Popu_tre" = 4,
                            "Pseu_men" = 3),
                          minimumPlotsPerGamm = 40,
                          constrainMortalityShape = list(
                            "Abie_bal" = c(15,25),
                            "Abie_las" = c(15,25),
                            "Betu_pap" = c(15,20),
                            "Lari_lar" = c(20,25),
                            "Pice_eng" = c(20,25),
                            "Pice_gla" = c(20,25),
                            "Pice_mar" = c(15,25),
                            "Pinu_ban" = c(15,25),
                            "Pinu_con" = c(15,25), 
                            "Popu_tre" = c(20,25),
                            "Pseu_men" = c(20,25)
                          ),
                          constrainGrowthCurve = list(
                            "Abie_bal" = c(0, 1),
                            "Abie_las" = c(0, 1),
                            "Betu_pap" = c(0, 1),
                            "Lari_lar" = c(0, 1),
                            "Pice_eng" = c(0, 1),
                            "Pice_gla" = c(0, 1),
                            "Pice_mar" = c(0, 1),
                            "Pinu_ban" = c(0, 1),
                            "Pinu_con" = c(0, 1), 
                            "Popu_tre" = c(0, 1),
                            "Pseu_men" = c(0, 1)
                          ),
                          quantileAgeSubset = list(
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
