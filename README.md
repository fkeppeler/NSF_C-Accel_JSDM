# Are Joint Species Models better than Single Species Models to predict species- and community-level responses to climate change?

### Project goals
The goal of our project is to fit a joint species distribution model (JSDM) to the northeast USA bottom trawl data and compare its results with single-species models (SDM). A second objective is to use the created JSDM to predict species and community response to future climate scenarios in the northeast USA.

### Repository contents
This repository contains R scripts used for data extraction and cleaning  (script 1), model construction (script 2), model fitting (script 3), hindcast (script 4), model comparison (scripts 5), and model extrapolation for future scenarios. Below, I listed some information about the databases, modeling approach, the metrics used for model comparisons, and the future scenario used.    

### Trawl data
The used trawl data from compiled and pre-processed by James Morley for his paper published at Plos One in 2018 ([link here](https://doi.org/10.1371/journal.pone.0196127)). The data contain thousand of trawl samplings from twenty
long-term surveys encompassing continental shelf regions around the USA and Canada. For this project, we only used trawl data from the east coast from 1980 to 2014. 

### Climate data
Salinity, oxygen partial pressure (PO2), and temperature from the bottom and surface of the ocean were derived from a dynamically downscaled Regional Ocean Modeling System (ROMS) with integrated biogeochemistry (Zhang et al. [2018](https://doi.org/10.1002/2017JC013402), [2019](https://doi.org/10.1029/2018JC014308)). These climate variables were then used in our models as fixed effects. 

### Habitat data
Grain size was obtained for each trawl sampling from Morley et al. ([2018](https://doi.org/10.1371/journal.pone.0196127)) to distinguish seafloor sediment type. Sea bed forms and seafloor complexity were extracted from layers created by The Nature Conservancy’s Northwest Atlantic Ecoregional Assessment ([NAMERA](https://www.conservationgateway.org/ConservationByGeography/NorthAmerica/UnitedStates/edc/reportsdata/marine/namera/namera/Pages/default.aspx)). Sea bed forms reflect the shape of the seafloor (10 - Depression/valley, 15 - Valley peaks, 20 - Low flat , 30 - Mid flat, 40 - Upper flat/bank, 45 - Upper slope/peak) and are derived from broad- and fine-scale Benthic Position Indices (BPI). Seafloor complexity is derived only from fine-scale BPI and represents how diverse the seafloor is in a given site by calculating how many fine-scale features is found within a 10 km buffer. The three habitats metrics (grain size, sea bed form, and seafloor complexity) were used in our models as fixed effects.

### Model structure 
SDM and JSDM were hierarchical generalized linear mixed models fitted with Bayesian inference through the R package HMSC (Ovaskainen et al. [2017](https://doi.org/10.1111/2041-210X.12723); Tikhonov et al. [2017](https://doi.org/10.1111/ele.12757)). SDM and JSDM had the same general structure:

* Response : Presence-absence of species.

* Fixed effects :  PO2 monthly average - Bottom, Salinity monthly average - Bottom, Salinity monthly average - Surface, Temperature monthly average - Bottom, Temperature monthly average - Surface, Annual maximum temperature - Bottom, Annual maximum temperature - Surface, Annual minimum temperature - Bottom, Annual minimum temperature - Surface, Grain size, Sea bed forms, Seafloor complexity, and vessel type (Albatross / Bigelow).

* Random effects : Latitude and Longitude coordinates, and sampling year.

Each prediction contained linear and squared effects in accordance with niche theory. The fixed effect component included only variables that had Spearman correlations lower than 0.75 to reduce multicollinearity. The addition of latitude and longitude in the models increased exponentially the computer processing time and led to convergence problems. So, we decided to simplify the spatial structure of the model by clustering the sites into grids (1x1 decimal degrees) and assigning the samples to the appropriate grid centroids. This allowed us to explore large-scale spatial effects in a fast and efficient manner without loosing too much information.

JSDM models had either all species sampled (~260 species) or only the most abundant ones (~60 species, frequency of occurrence greater than 5%). For SDM, we conducted a model for each species individually and stacked them together for predictions about richness and species composition (see below). 

We configured HMSC models with default non-informative priors following Ovaskainen and Abrego ([2020](https://doi.org/10.1017/9781108591720)). The JSDM were built with a Probit distribution. The posterior distributions were sampled with 5 chains (Markov Chain Monte Carlo – MCMC). Each chain was composed of 150,000 iterations, from which 20% was burn-in and from the remaining only 250 were retained with a thin value of 500. 

### Model evaluation
We trained the models with ~ 90% of the data and tested it with the remaining portion of the data. The train - test division was done in two ways: 1) Interpolation, where the division was done randomly; and 2) Extrapolation, where train data was composed of all samples conducted from 1980 to 2011 and test data was composed of all samples conducted from 2012 to 2014. We evaluated our models using four metrics (accuracy, discrimination, calibration, and precision) at two different ecological scales (species-level and community-level [richness and community composition similarity]) following the approach proposed by Norberg et al. ([2019](https://doi.org/10.1002/ecm.1370)) 

### Future Scenarios
We predicted community responses to the SSP5-8.5 future socio-economic scenario using our JSDM. More specifically, salinity, temperature, and pO2 under the SSP5-8.5 scenario were simulated for the 2045-2055 period using ROMS and a dynamic delta forcing approach (Pozo Buil et al. [2021](https://doi.org/10.3389/fmars.2021.612874)) and a single CMIP6 model (GFDL) as a driver.  For this exercise, we assumed that habitat variables (grain size, sea bed forms, and seafloor complexity) did not change with time and that trawl samples were conducted with the Albatross vessel.



