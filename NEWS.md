Known issues: <https://github.com/PredictiveEcology/Biomass_speciesParameters/issues>

# Biomass_speciesParameters 3.0.2.9003

* `balanceGrowth = TRUE` no longer stops when a fitted species has no `speciesEcoregion` rows (e.g. it is below the ecoregion support threshold everywhere in the study area, as Betu_pap in several BC TSAs). Such a species cannot establish, so it is left out of the balance, keeps its fitted traits, and a message names it.

# Biomass_speciesParameters 3.0.2.9002

* Fixed two bugs in the 3.0.2.9001 additions, found in the first FIIS run with them. `relabelHybridSprucePSP` matched plots by the factor codes of `OrigPlotID1` (PSPclean returns a factor) instead of by name, so almost no BC plot got its BEC zone and none was relabelled. `medianMaxB` failed with `balanceGrowth = TRUE` when `maxB` is integer, because `median()` returns an integer for an odd number of rows and a double for an even number.

# Biomass_speciesParameters 3.0.2.9001

* BC interior (hybrid-zone) spruce PSPs are recorded as "Picea glauca", so they matched nothing, and were dropped, when the study area's `sppEquiv` has the hybrid spruce merged into Engelmann (`LandR.mergeHybridSpruce = "engelmann"`). New parameter `mergeHybridSprucePSP` (default: follow that option; NA or "white" for no relabel) relabels those BC records as "Picea engelmannii x glauca", so they count as `Pice_eng`. BC plots in the BEC zones of the new parameter `excludeBECzonesHybridSpruce` (default `BWBS`, `SWB`: boreal white spruce) are not relabelled. The BEC zones come from the new optional input `BECzonesBC`, else from the BC Data Catalogue (`bcdata`, added to `reqdPkgs`); if neither is available there is a warning and no relabel.
* Missing packages in `reqdPkgs`: added `dplyr`.
* New opt-in `balanceGrowth` (default `FALSE`, so nothing changes unless set), with `sharedGrowthcurve` and `targetK`: all fitted species get one `growthcurve`, and `mANPPproportion` is set so K = `mANPPproportion` * maxB^(1 - `growthcurve`) is equal across species, since under LandR competition any difference in K is winner-take-all while the per-species fitted `growthcurve` values are statistically indistinguishable. Both traits are then rounded to 3 decimals. The fit of the new values to the PSP data is returned in the new output `speciesBalanceCheck`, with a warning for species more than 2 log-likelihood units worse than their own best.

# Biomass_speciesParameters 3.0.2

* `cohortDataFactorial` and `speciesTableFactorial` are now local variables of `Init` instead of `mod$` objects. They were only ever used inside `Init`, but as `mod$` objects the ~9 GB stayed in the simList's `.modObjs` (and in every `Copy(sim)`) for the rest of the simulation.

# Biomass_speciesParameters 3.0.1

* new "LANDIS mode": added a `landis` parameter (default `FALSE`) that exposes the fitted, scaled non-linear biomass-over-age growth curves per species as a new output `speciesGrowthCurvesLandis` for use as inputs to LANDIS-II Biomass Succession; also added a `speciesGrowthCurvesPSP` output holding the PSP points behind the fitted curves.
* migrated to the new PSPclean data format and approach; trim extraneous columns and bump the PSPclean dependency accordingly.
* robustness and cleanup: added error checks, fixed a `Cache()` call, partial-match fixes, and `#nolint` tags; merged contributor PRs #53 and #54.

# Biomass_speciesParameters 3.0.0 (2025-11-17)

* major backend and tooling modernization: replaced the disk.frame intermediate-data backend with an arrow-based backend, swapped crayon for cli messaging, and moved from magrittr to the native `|>` pipe (reqdPkgs: added `arrow`, `cli`, `fs`; removed `disk.frame`, `crayon`, `magrittr`).
* reworked intermediate data storage: disk.frame objects written to `outputPath` (fix #49) with a `.df` extension, ahead of the arrow migration.
* changed the default `sppEquivCol` to "LandR".
* fixes: guard against an empty PSP study-area polygon, a deltaDiff-by-quantile height fix, a forest-plot subset fix, and a join issue arising from the source column.
* extensive vignette/Rmd rewrite plus CI updates to the render-module-rmd workflow.

# Biomass_speciesParameters 2.0.2 (2024-06-06)

* cache the `modifiedSpeciesTables` step (PR #47).
* consolidate all module functions into `init` (PR #46).

# Biomass_speciesParameters 2.0.1 (2024-05-31)

* "deGAMM": replaced the GAMM-based fitting with growth curves (gamms to gcs), removed `speciesGAMMs` from the simList, and renamed the related objects and functions (large `aNPPfunctions.R` rewrite).
* added New Brunswick (NB) PSP data support and removed ~200 lines of inlined PSP input.
* performance and RAM management: sample `standAge` from the factorial instead of `unique()` over ~280 million rows, skip sampling when the age vector is below the 1e5 sample size, and add `gc()` calls.
* fixes for age-1 stands, improved bad/no-data model detection, plotting fixes, metadata and checksum cleanup, and added Biomass_speciesFactorial to `loadOrder`.

# Biomass_speciesParameters 2.0.0 (2023-10-03)

* non-linear growth models (PR #21) with a switch to `focal` as the default fitting approach, and separation of the maxB correction from trait estimation.
* switched PSP data handling to disk.frame and adopted PSPclean@development.
* use the new `speciesTableFactorial` / Biomass_speciesFactorial factorial object, and derive `speciesTable` via LandR when not supplied.
* removed the spatExtent dependency ("no more spatExtent"); default `sppEquivCol` "Boreal".
* added a GitHub Actions workflow to build the Rmd; fixes for integer coercion, trait averaging for non-estimated species, and Gompertz/`gamm` predict depth.
