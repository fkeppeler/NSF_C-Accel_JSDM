# Are Joint Species Models better than Single Species Models to predict species- and community-level responses to climate change?

## Project goals
The goal of our project is to fit a joint species distribution model (JSDM) to the northeast USA bottom trawl data and compare its results with single-species models (SDM). A secondary objective is to use the JSDM created to predict species and community response to future climate scenarios in the northeast USA.

## Repository details
This repository contains R scripts used for data extraction and cleaning  (script 1), model construction (script 2), model fitting (script 3), model comparison, and model extrapolation for future scenarios. Below, a listed some information about the databases, modeling approach, the metrics used for model comparisons, and the future scenario used.    

### Trawl data
The used trawl data from compiled and pre-processed by James Morley for his paper published at Plos One in 2018 ([link here](https://doi.org/10.1371/journal.pone.0196127)). The data contain thousand of trawl samplings from twenty
long-term surveys encompassing continental shelf regions around the USA and Canada. For this project, we only used trawl data from the east coast from 1980 to 2014. 

### Climate data
Salinity, oxygen partial pressure (PO2), and temperature from the bottom and surface of the ocean were derived from a dynamically downscaled Regional Ocean Modeling System (ROMS) with integrated biogeochemistry (Zhang et al. [2018](https://doi.org/10.1002/2017JC013402), [2019](https://doi.org/10.1029/2018JC014308)). These climate variables were then used in our models as fixed effects. 

### Habitat data
Grain size was obtained for each trawl sampling from Morley et al. ([2018](https://doi.org/10.1371/journal.pone.0196127)) to distinguish seafloor sediment type. Sea bed forms and seafloor complexity were extracted from layers created by The Nature Conservancy’s Northwest Atlantic Ecoregional Assessment ([NAMERA](https://www.conservationgateway.org/ConservationByGeography/NorthAmerica/UnitedStates/edc/reportsdata/marine/namera/namera/Pages/default.aspx)). Sea bed forms reflect the shape of the seafloor (10 - Depression/valley, 15 - Valley peaks, 20 - Low flat , 30 - Mid flat, 40 - Upper flat/bank, 45 - Upper slope/peak) and are derived from broad- and fine-scale Benthic Position Indices (BPI). Seafloor complexity is derived only from fine-scale BPI and represents how diverse the seafloor is in a given site by calculating how many fine-scale features is found within a 10 km buffer. The three habitats metrics (grain size, sea bed form, and seafloor complexity) were used in our models as fixed effects.


### Modeling approach 


### Model comparisons


### Future Scenarios






