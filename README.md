# Fires_and_Cocci

This repo contains the code associated with the paper *Association between wildfires and coccidioidomycosis incidence in California, 2000–2018: a synthetic control analysis*.

Human case data are protected health information (PHI) with access restricted to authorized California Department of Public Health (CDPH) staff. Case data is thus not included in this repository. 

The code is divided into four folders:

## Analysis 

* **Results Data** contains csv files for the effect estimates obtained from our gsynth analysis

* **GsynthControlSensitivityAnalysis.R** contains the code for the control sensitivity analysis in which each control is removed one at a time from the eligible control pool and calculating the effect of the wildfire with this control excluded.

* **InteractionAnalysis.R** regresses the number of high wind days in a region on whether it occurs in a fire-exposed region during the month of the fire. We reported on the significance of this interaction.

 * **MesmaItsAnalysis.R** contains the code for conducting an Interrupted Time Series analysis to ascertain whether the levels of vegetative land cover (MESMA) were significantly different in the months following the fire. 

* **PoolEffects.R** contains the script to pool effect estimates from multiple fires across different time periods using a fixed effects meta analysis from the R package meta.

## Covariates
* **data** contains various environmental variables extracted for fire and control regions, summarized monthly. Case data are PHI and thus not included.

* **GetFireElevationData.R** script to extract elevation data from USGS rasters to the centroid of the hexagonal boundary surrounding each fire.

* **prism.R** script to extract temperature and precipitation data from rasters obtained from Prism Climate Group. The data are summarized monthly within each fire and control region.

## Hexagon Boundaries
* Shape files for control grid and hexagonal boundaries surrounding fires. The main analysis used boundaries with radius 20km

* **getHexagonalBoundary.R** contains functions to draw a hexagonal boundary centered on the center of a burned area.

* **PostFireOverlaps.csv** identifies the fire-exposed regions that overlap with another fire during the three years pre- or post-ignition. Overlaps are reported for each of the hexagonal buffer sizes examined (15, 20, 25 km).

## gsynth
* **getCaseandCovdata.R** provides functions to format a data frame for use in gsynth. For each fire id, eligible control regions are identified, and the associated covariate and case data are formatted into a dataframe where the time variable is measured as the number of months before or after the fire ignition date.

* **gsynth.R** provides the script used to obtain effect and uncertainty estimates using the generalized synthetic control method implemented by the R package gsynth. Each fire effect is estimated separately as the differing ignition dates required each fire to have a separate pool of controls. 

* **MesmaPositiveControl.R** contains the script to run the generalized synthetic control method using percent vegetative land cover (MESMA) as outcome, and pool estimates using a fixed effects meta analysis.

