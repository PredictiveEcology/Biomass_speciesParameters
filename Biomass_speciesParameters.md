---
title: "LandR _Biomass_speciesParameters_ Manual"
date: "Last updated: 2026-09-21"
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





[![module-version-Badge](/home/runner/work/Biomass_speciesParameters/Biomass_speciesParameters/figures/moduleVersionBadge.png)](https://github.com/PredictiveEcology/Biomass_speciesParameters6b1a39c099711e300318477797fc7b43ce049659)

[![Issues-badge](/home/runner/work/Biomass_speciesParameters/Biomass_speciesParameters/figures/issuesBadge.png)](https://github.com/PredictiveEcology/Biomass_speciesParameters/issues)

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
maximum aboveground net primary productivity (`maxANPP`) traits (see [Parameter
estimation/calibration](#bsppparam-calib)), the latter two in conjunction with
the module *Biomass_borealDataPrep*.

This module **will not** obtain other traits or parameters used in
*Biomass_core* and so it is meant to be used in conjunction with another
data/calibration module that does so (e.g., *Biomass_borealDataPrep*). However
it can be used stand-alone in an initial developmental phase for easier
inspection of the statistical calibration procedure employed.

As of September 21, 2026, the *raw* PSP data used in this
module is not freely available, and data sharing agreements must be obtained
from the governments of SK, AB, and BC to obtain it. However, the *processed and
anonymized* PSP data is provided via a Google Drive folder accessed
automatically by the module.

\*Google Account is therefore necessary to access the data used for
calibration.\*\*

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
    *Biomass_core* simulations on landscapes populated by one or more species,
    each simulation using a different set of species trait values.

## Module manual

### General functioning {#bsppparam-general-functioning}

Tree cohort growth and mortality in *Biomass_core* are essentially determined by
five parameters: `growthcurve`, `mortalityshape`, maximum biomass (`maxB`),
maximum aboveground net primary productivity (`maxANPP`) and `longevity`.

The `growthcurve` and `mortalityshape` parameters (called 'growth curve' and
'mortality shape' in LANDIS-II Biomass Succession Extension v3.2, the base model
for *Biomass_core*) strongly modulate the shape of species growth curves and so
it is important that they are calibrated to the study area in question.

Also, the growth and mortality equations used in *Biomass_core* are non-linear
and their resulting actual biomass accumulation curve is an emergent phenomenon
due to competition effects. This means that the ideal trait/parameter values
should not be estimated on pure single species growth conditions, as their
resulting dynamics will be different in a multi-species context.

*Biomass_speciesParameters* attempts to address these issues (at least
partially) using a "curve-matching" approach. It compares the best fit
(according to their AIC) of three non-linear forms (Chapman-Richard's, Gompertz,
and a logistic form) fitted to permanent sample plot (PSP) data to a large
collection of theoretical (i.e. simulated) species curves, each representing a
different set of the five key parameters that govern biomass increment in
`Biomass_core`: `growthcurve`, `mortalityshape`, the ratio of `maxANPP` to
`maxB`, and `longevity`. This library of curves is produced by the
*Biomass_speciesFactorial* module.

*Biomass_speciesParameters* generally follows other LandR data modules, like
*Biomass_boreaDataPrep*, which also attempts to calibrate previously estimated
spatially varying species traits such as `maxB` and `maxANPP` from the input
data layers.

#### Permanent sample plot data {#bsppparam-PSPdata}

*Biomass_speciesParameters* can use all the PSP data available (note that it may
span several thousands of kilometres), or select the data based on a shapefile
(`studyAreaANPP`; see [List of input objects](#bsppparam-inputs-list)).

By default, all available PSP is obtained, via the `PSPdataTypes` parameter. The
particular data sets will depend on the version of `ianmseddy/PSPclean` that is
installed (assuming requisite file-sharing permissions are available). This may
include data from the provinces of BC, AB, SK, ON, QC, and NB, as well as the
National Forest Inventory. If no file-sharing agreements are in place, the dummy
option can be used. This will obtain data from BC, AB, SK, and the NFI that were
previously treated for errors and standardized into a single data set with the
exact location and identifying attributes anonymized. However, it should be
noted that data for BC, Quebec, New Brunswick, and the NFI are freely available,
thus no file-sharing agreement is necessary for these jurisdictions.

The data include individual species, diameter at breast height (DBH), and
sometimes tree height measurements for each tree in a plot, as well as stand
age. As part of the standardization process, dead trees were removed from the
data set. Tree biomass was then per estimated species using either a DBH-only
model or a DBH-height model from @LambertEtAl2005, in $g/m^2$.

Note that the model used to calculate biomass can also be changed to
@UngEtAl2008 via the `P(sim)$biomassModel` module parameter (see [list of
parameters](#bsppparam-params-list)).

#### Simulated species data {#bsppparam-simdata}

The *Biomass_speciesFactorial* module was used to create a library of
theoretical species curves (biomass accumulation curves, to be more precise) to
which the best non-linear model form fit to the PSP-biomass will be matched for
each species and species combinations in the study area landscape. The library
of curves was created by running several *Biomass_core* simulations with no
reproduction, competition, disturbance, or dispersal effects. The species in the
simulations encompassed the full factorial of traits, growing alone and in
competition with one cohort of each other combination. Each simulation differed
in the combination of species trait values that influence growth and mortality
dynamics, namely: `growthcurve`, `mortalityshape`, `longevity`, `maxANPP` and
maximum biomass (`maxBiomass`, not to be confused with the data-driven `maxB`
which is later calibrated).

The values for `maxANPP` were explored via the `mANPPproportion`, the ratio of
`maxANPP` to `maxBiomass` (the parameter used for theoretical curves), as it
reflects their relationship. While any factorial combination is possible to run
using the module *Biomass_speciesFactorial*, the default objects used by
*Biomass_speciesParameters* utilized `growthcurve` values ranging from 0.65 to
0.85 in increments of 0.02, `mortalityshape` values of 20, 22, and 24,
`mANPPproportion` values from 3.5 to 6.5 in increments of 0.5, and `longevity`
values of 125 to 700 in increments of 25, all of which total over 6.5 million
unique combinations of traits.

Results from these simulations were compiled into a table (`cohortDataFactorial`
; see [List of input objects](#bsppparam-inputs-list)) that is accessed by
*Biomass_speciesParameters*, so that the module can be run without needing to
re-simulate the theoretical curves.

#### Parameter estimation/calibration {#bsppparam-calib}

*Biomass_speciesParameters* calibrates `growthcurve`, `mortalityshape` and
`mANPPproportion` by matching the theoretical species curves produced by
*Biomass_speciesFactorial* (`cohortDataFactorial`) against observed species
growth curves from the PSP data.

The parameter `P(sim)$speciesFittingApproach` determines which of four possible
fitting approaches to use: all, single, pair-wise, or focal. The `all` approach
combines all PSPs into a single species, thus parameterizing all species
identically. This is not intended to have an ecological application. The
`"single"` approach is the earliest method, where each species utilizes plots
where 50% or more of the biomass is composed solely of that species. Non-linear
models for each species are fit to these observations independently of each
other, with each species possessing its own subset of plots. However, this
approach is unable to accurately characterize the competitive effects that arise
when multiple species occupy the same plot. Therefore it is suitable only for
characterizing species that grow in pure stands, and has been retained largely
for backwards compatibility.

The `"pairwise"` approach retains all plots where exactly two species each
represent more than 20% of the total plot biomass, e.g. a plot with 40% Pinus
contorta and 39% Populus tremuloides is retained only if the remaining 21% of
biomass is composed of more than one species. Then, for each combination of 2
species-of-interest, separate non-linear models are fit for each species. #TODO:
I believe the models are ultimately combined but this might not happen until
traits are selected. Of the three approaches, the `pairwise` is best able to
account for competition, but is most sensitive to data availability and can be
confounded by combinations of species that seldom occur together.

The third approach `speciesFittingApproach`, `"focal"`, is the default. This
approach is similar to `"pairwise"`, but combines the species that are *not* the
species-of-interest (or 'focal' species). For example, in a plot with biomass
composition of 20% Pinus contorta, 40% Populus tremuloides, and 21% Picea
mariana, when fitting the P. contorta equations, the biomass of P. tremuloides
and P. mariana are combined, and when fitting P. tremuloides, the biomass of P.
contorta and P. mariana are combined. In comparison to the `pairwise` approach,
this approach sacrifices some detail but is less sensitive to data quality and
availability.

Before calculating the *observed* species growth curves (i.e., the best of three
non-linear forms to match PSP data), the module subsets the PSP data to stand
ages below the 95th percent quantile for all species (this can be changed via
the `P(sim)$quantileAgeSubset` module parameter), as records for larger age
classes were limited and constituted statistical outliers. In some species,
changing the quantile value may improve results, however. Two examples are
*Pinus banksiana* and *Populus sp* (in western Canada), for which using the 99th
percent quantile improved the models, because these are short-lived species for
which data at advanced ages is scarce. In addition, weights are added at the
origin (age = 0 and biomass = 0) to force the intercept to be essentially at 0
age and 0 biomass.

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
model convergence. In this case, *Biomass_speciesParameters* skips parameter
calibration. Consequently the module will interpolate the `mANPPproportion`,
`growthcurve`, and `mortalityshape` from the respective means of other species.
These species are tracked via the `source` column in `sim$species`, which will
be one of `"interpolated"` or `"estimated"`. This mechanism ensures these
species remain competitive, as the parameterized traits are often significant
departures from defaults used in many LANDIS-II applications to Canada's boreal
forests.

After each species best fit is selected (using AIC), *Biomass_speciesParameters*
compares it to the library of theoretical curves, and picks the best one based
on maximum likelihood. This best theoretical curve will be associated with a
given combination of `growthcurve`, `mortalityshape` and `maxANPPproportion`
values, which are then used directly as the calibrated values, in case of
`growthcurve` and `mortalityshape`, or to calibrate `maxANPP` in the case of
`maxANPPproportion` (see below).

Because simulated growth curves never achieve the maximum biomass parameter (the
`maxBiomass` parameter set to 5000 for all simulations of theoretical species
curves, or the `maxB` parameter in *Biomass_core* simulations), it acts as an
asymptotic limit that reflects the potential maximum biomass for a species in an
ecolocation (ecological zone and land cover combination).

*Biomass_speciesParameters* uses the ratio between the potential maximum biomass
(`maxBiomass`, always 5000) to the achieved maximum biomass in the theoretical
curves, to rescale `maxB`. This ratio is called the `inflationFactor` and it is
multiplied by `maxB` values previously estimated from data (e.g. by
*Biomass_borealDataPrep*). This way, species simulated in *Biomass_core* are
able to achieve the maximum observed biomass used to *initially* estimate
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

### List of input objects {#bsppparam-inputs-list}

The full list of input objects required by the module is presented below (Table
\@ref(tab:moduleInputs2-Biomass-speciesParameters)). The input `studyAreaANPP`
(the study area used extract the PSP data from) is optional, and therefore no
default is supplied. All other input objects have internal defaults, but the
user may need to request access to their online files.

Of these inputs, the following are particularly important and deserve special
attention:

**Spatial layers**

```         
-   `studyAreaANPP` -- shapefile. An `sf` or `SpatVector` object 
determining the geographic extent of the PSP data. 
```

-   **Tables**

    -   `speciesTableFactorial` and `cohortDataFactorial` -- a tables of species
        trait combinations and the theoretical species growth curve data
        (respectively)

    -   `PSPmeasure_sppParams`, `PSPplot_sppParams` and `PSPgis_sppParams` --
        tree measurement, biomass growth and geographical data of the PSP
        datasets used to build observed species growth curves.

    -   `species` -- a table of invariant species traits that may have been
        produced by another module. It **must** contain the columns 'species',
        'growthcurve' and 'mortality shape', whose values will be calibrated.

    -   `speciesEcoregion` -- table of spatially-varying species traits that may
        have been produced by another module. It **must** contain the columns
        'speciesCode', 'maxB' and 'maxANPP' and 'ecoregionGroup' (the
        ecolocation ID). 'maxB' and 'maxANPP' values will be calibrated by
        species.

```{=tex}
\newpage
\blandscape
```
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
   <td style="text-align:left;"> cohortDataFactorial_path </td>
   <td style="text-align:left;"> fs_path </td>
   <td style="text-align:left;"> Path where the `cohortDataFactorial` object is saved as an `arrow` dataset. A large `cohortData` table ( sensu `Biomass_core`) with columns `age`, `B`, and `speciesCode` that joins with `speciesTableFactorial`. See `PredictiveEcology/Biomass_factorial` for further information. </td>
   <td style="text-align:left;"> https://drive.google.com/file/d/1NH7OpAnWtLyO8JVnhwdMJakOyapBnuBH/ </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PSPmeasure_sppParams </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Merged PSP and TSP individual tree measurements. Must include the following columns: `MeasureID`, `OrigPlotID1`, `MeasureYear`, `TreeNumber`, `Species`, `DBH` and `PSP`, where `Species` corresponds to species names in `LandR::sppEquivalencies_CA$Latin_full`. Defaults to randomized PSP data stripped of real `plotID`s </td>
   <td style="text-align:left;"> https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PSPplot_sppParams </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Merged PSP and TSP plot data. Defaults to randomized PSP data stripped of real `plotID`s. Must contain columns `MeasureID`, `MeasureYear`, `OrigPlotID1`, and `baseSA`, the latter being stand age at year of first measurement </td>
   <td style="text-align:left;"> https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PSPgis_sppParams </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Plot location `sf` object. Defaults to PSP data stripped of real `plotID`s/location. Must include field `OrigPlotID1` for joining to `PSPplot` object </td>
   <td style="text-align:left;"> https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing </td>
  </tr>
  <tr>
   <td style="text-align:left;"> species </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> A table of invariant species traits with the following trait colums: 'species', 'Area', 'longevity', 'sexualmature', 'shadetolerance', 'firetolerance', 'seeddistance_eff', 'seeddistance_max', 'resproutprob', 'mortalityshape', 'growthcurve', 'resproutage_min', 'resproutage_max', 'postfireregen', 'wooddecayrate', 'leaflongevity' 'leafLignin', and 'hardsoft'. Only 'growthcurve', 'hardsoft', and 'mortalityshape' are used in this module. Default is from Dominic Cyr and Yan Boulanger's applications of LANDIS-II </td>
   <td style="text-align:left;"> https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv </td>
  </tr>
  <tr>
   <td style="text-align:left;"> speciesEcoregion </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Table of spatially-varying species traits (`maxB`, `maxANPP`, `establishprob`), defined by species and `ecoregionGroup`). Defaults to a dummy table based on dummy data of biomass, age, ecoregion and land cover class. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> speciesTableFactorial_path </td>
   <td style="text-align:left;"> fs_path </td>
   <td style="text-align:left;"> Path where the `speciesTableFactorial` object is saved as an `arrow` dataset. A large species table ( sensu `Biomass_core`) with all columns used by Biomass_core, e.g., `longevity`, `growthcurve`, `mortalityshape`, etc., when it was used to generate `cohortDataFactorial`. See `PredictiveEcology/Biomass_factorial` for futher information. </td>
   <td style="text-align:left;"> https://drive.google.com/file/d/1NH7OpAnWtLyO8JVnhwdMJakOyapBnuBH/ </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppEquiv </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Table of species equivalencies - see `?LandR::sppEquivalencies_CA`. Traits will be estimated for each unique entry in the `sppEquivCol` column. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppEquivLong </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> The full table of species equivalencies - see `?LandR::sppEquivalencies_CA`. Biomass will be estimated for each species based on the `sp_Biomass_eq' column, which uses `pemisc::biomassCalculation` to derive AGB from DBH and height (based on the equations from https://doi.org/10.1139/x05-112). The full table is used to improve stand biomass estimates even if some species are not of interest. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyAreaANPP </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Optional study area used to crop PSP data before building growth curves. If supplied, an ecoregion-scale object is recommended, at a minimum. </td>
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

```         
-   `standAgesForFitting` -- determines the range of ages for which the fit of the growth curves is evaluated against the non-linear model. It should be a subset of the full growth curve, because there is extremely limited plot data for stands aged 0-20, as well as towards the end of a species' longevity.

-   `speciesFittingApproach` -- should the calibration take into account 
species growing in single- or multi-species context?

-   `quantileAgeSubset` -- upper quantile age value used to subset PSP data.
It can be a vector named by species, or a single numeric. It is used to limit the outsized influence of old stands at the tail end of the stand age 
distribution. However, it will have varying effects by species and PSP data.
```

**Data processing**

```         
-   `PSPdataTypes` -- which jurisdictional plot data to use, important 
because the default will attempt to use PSP data that may be inaccessible

-   `sppEquivCol` -- the column name in `sim$sppEquiv` that ultimately determines which species are modeled, and whether any are combined. It should be complete for every row, and duplicates are merged together. For example, a table with `"Picea englemannii"` and `"Picea glauca"` in the "Latin" column, and "White spruce" in the "EN_generic_short" column would result in the two rows merged into one single growth curve if "EN_generic_short" were the sppEquivCol, but have two growth curves if `Latin` were the sppEquivCol. The estimate biomass estimation is always dependent upon the column "PSP", and missing entries will lead to the species' omission. 
```

```{=tex}
\newpage
\blandscape
```
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
   <td style="text-align:left;"> landis </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> If `TRUE`, run in 'LANDIS mode': expose the fitted LANDIS-version growth curves (the scaled non-linear biomass-over-age curves shown in the LandR-vs-non-linear plot) per species as the output object `speciesGrowthCurvesLandis`, for use as inputs to LANDIS-II Biomass Succession. The per-species growth-curve parameters (`growthcurve`, `mortalityshape`, etc.) are written to `species` regardless. Default `FALSE` preserves the standard behaviour. </td>
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
   <td style="text-align:left;"> minimumPlots </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 50 </td>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Minimum number of PSP plots per species </td>
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
   <td style="text-align:left;"> Which PSP datasets to source, defaulting to all. Other available options include 'BC', 'AB', 'SK', 'ON', 'NB', 'NFI', and 'dummy'. 'dummy' should be used for unauthorized users. </td>
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
   <td style="text-align:left;"> 99 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 100 </td>
   <td style="text-align:left;"> Quantile by which to subset PSP data. As older stands are sparsely represented the oldest measurements become vastly more influential. This parameter accepts both a single value and a list of vectors, named according to `sppEquivCol`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> speciesFittingApproach </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> focal </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Either 'all', 'pairwise', 'focal' or 'single', indicating whether to pool all species into one fit, do pairwise species (for multiple cohort situations) do pairwise species, but using a focal species approach where all other species are pooled into 'other' or do one species at a time. If 'all', all species will have identical species-level traits. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppEquivCol </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> LandR </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> The column in `sim$sppEquiv` data.table that defines individual species. The names should match those in the species table. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> standAgesForFitting </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 21, 91 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> The minimum and maximum ages of the biomass-by-age curves used in fitting. It is generally recommended to keep this param under 200, given the low data availability of stands aged 200+, with some exceptions. For a closed interval, end with a 1, e.g. `c(31, 101)`. </td>
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
   <td style="text-align:left;"> .studyAreaName </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Human-readable name for the growth curve filename. If `NA`, a hash of sppEquiv[[sppEquivCol]] will be used. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> .inputOb.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useParallel </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> maximum number of threads/workers to use for data.table operations; passed to `data.table::setDTthreads` and should be &lt;= 4. </td>
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

```         
-   `species` and `speciesEcoregion` -- tables with calibrated trait values.

-   `speciesGAMMs` -- the fitted GAMM model objects for each species.
```

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
   <td style="text-align:left;"> The updated invariant species traits table (see above). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> speciesEcoregion </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> The updated spatially-varying species traits table (see description for this object in inputs) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> speciesGrowthCurves </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> list containing each species' non-linear model, model data, and the unfiltered PSP data </td>
  </tr>
  <tr>
   <td style="text-align:left;"> speciesGrowthCurvesLandis </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> An empty `data.table` unless `P(sim)$landis` is `TRUE`, when it holds the fitted LANDIS-version growth curves (`BscaledNonLinear`) by `species` and `standAge` (the scaled non-linear curves shown in the `LandR_VS_NLM_growthCurves` plot), for use as LANDIS-II Biomass Succession inputs. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> speciesGrowthCurvesPSP </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> An empty `data.table` unless `P(sim)$landis` is `TRUE`, when it holds the PSP observations used to fit the growth curves (`biomass` by `standAge` and species, with `OrigPlotID1` for joining to plot locations / ecoregion), a diagnostic to plot against `speciesGrowthCurvesLandis`. </td>
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
@BarrosEtAl2023.

### Load `SpaDES` and other packages.

### Set up R libraries {#bsppparam-example-libs}


``` r
tempDir <- tempdir()

pkgPath <- file.path(tempDir, "packages", version$platform,
                     paste0(version$major, ".", strsplit(version$minor, "[.]")[[1]][1]))
dir.create(pkgPath, recursive = TRUE)
.libPaths(pkgPath, include.site = FALSE)
 
repos <- c("predictiveecology.r-universe.dev", getOption("repos"))
options(repos = repos)
install.packages("SpaDES.project")   ## gets Require too
```

### Get the module and module dependencies {#bsppparam-example-pkg-mods}


``` r
library(Require)

paths <- list(inputPath = normPath(file.path(tempDir, "inputs")), 
              cachePath = normPath(file.path(tempDir, "cache")), 
              modulePath = normPath(file.path(tempDir, "modules")), 
              outputPath = normPath(file.path(tempDir, "outputs")))

SpaDES.project::getModule(modulePath = paths$modulePath,
                          c("PredictiveEcology/Biomass_speciesParameters@main"),
                          overwrite = TRUE)
Require::Require("SpaDES.core (>= 2.1.4)")
Require("googledrive")
```

### Setup simulation


``` r
times <- list(start = 0, end = 1)

modules <- list("Biomass_speciesParameters")

objects <- list()
inputs <- list()
outputs <- list()
parameters <- list()
mySim <- SpaDES.core::simInitAndSpades(times = times, 
                          params = parameters, 
                          modules = modules, 
                          paths = paths, 
                          objects = objects)

## to inspect the fitted growth models:
mySim$speciesGrowthCurves$Pice_mar
```

## References {#bsppparam-refs}

