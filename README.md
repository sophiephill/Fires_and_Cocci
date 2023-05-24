# Fires_and_Cocci

This repo contains the code associated with the paper, *Association between wildfires and coccidioidomycosis incidence in California, 2000–2018: a synthetic control analysis*.

Human case data are protected health information (PHI) with access restricted to authorized California Department of Public Health (CDPH) staff. Case data is thus not included in this repository. 

The code is divided into four folders:

1. **Analysis** contains the code to analyze the effect estimates obtained from the generalized synthetic control analysis and covariate data. Within this folder, the R scripts are as follows:
* **Results Data** contains csv files for the effect estimates obtained from our analyses

* **GsynthControlSensitivityAnalysis.R** provides the code for calculating the effect estimates when each of the controls for a fire is left out.

* **InteractionAnalysis.R** regresses the number of high wind days in a region on whether it occurs during the month of the fire and in a fire-exposed region, examining the significance of the interaction term. 

* **MesmaItsAnalysis.R** contains the code for conducting an Interrupted Time Series analysis to ascertain whether the levels of vegetative land cover (MESMA) were significantly different in any month following the fire. 

* PoolEffects.R

2. **Covariates**
* data 
* GetFireElevationData.R 
* prism.R

3. **Hexagon Boundaries**
* Shape files for control grid and hexagonal boundaries surrounding fires. The main analysis used boundaries with radius 20km
* getHexagonalBoundary
* RmPostFireOverlaps

4. **gsynth**
* getCaseandCovdata.R

* gsynth.R

* **MesmaPositiveControl.R** contains the script to run the generalized synthetic control method using percent vegetative land cover (MESMA) as outcome, and pool estimates using a fixed effects meta analysis.