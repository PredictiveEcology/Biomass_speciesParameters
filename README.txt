
Any other details that a user may need to know, like where to get more information,
where to download data, etc.

This module uses LandR to simulate biomass for a factorial combination of species traits. 
Each species grows as a single cohort with no understory (ie, no dispersal, regeneration, or disturbance). 
It contains a 'frozen' version of Boreal_Biomass (called LBMR at the time), that has been modified
to allow for understory regeneration, as well as to avoid some bugs in the code (read, errors that are edge cases).
The full factorial included mortalityshape 5 to 25, in increments of 1, growthcurve 0 to 1, in increments of 0.1, 
mANPPproportion (the proportion of max ANPP to maxB) from 1 to 5 in increments of .25, and longevity from 150 to 700 in increments of 25. 
The results were saved in a combination of tables so that the module can be run without needing to simulate the factorial.


