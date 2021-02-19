
This module matches simulated growth curves of theoretical species with varying maxANPP, growth curve and mortality shape values (using LandR Biomass_core) against growth curves derived from Permanent Sample Plots (PSP data) across Western Canada to find the most likely combination of traits, and to adjust the maxB trait in respect to the observed maximum biomass on the landscape.

As of 2020-04-08, the PSP data needed for this module is not freely available, and data sharing agreements must be obtained from the governments of SK, AB, and BC. 

Each species grows as a single cohort with no understory (ie, no dispersal, regeneration, or disturbance).
The full factorial included mortalityshape 5 to 25, in increments of 1, growthcurve 0 to 1, in increments of 0.1, 
mANPPproportion (the proportion of max ANPP to maxB) from .25 to 10 in increments of .25, and longevity from 150 to 700 in increments of 25. 
The results were saved in a combination of tables so that the module can be run without needing to simulate the factorial.

The module is inteded to be used in combination with Biomass_borealDataPrep. 
However, users may want to experiment with the traits as well as examine the output GAMMs used to select traits. 
Therefore, it can be run as a stand-alone. The parameters are applied to all species (rather than each individually),
thus users may want to run several times, varying parameters, select the traits they feel best capture the true growth curve
(the PSPs are quite limited in their representation of stand dynamics for those stands older than 200 yrs), and modify the traits
directly using Biomass_borealDataPrep. 


