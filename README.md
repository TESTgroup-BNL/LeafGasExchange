# LeafGasExchange
An R package for fitting and simulating leaf level gas exchange.

This package works using 4 main functions:
* The function f.make.param allows to build a list of photosynthetic parameters (Vcmax, Jmax, kinetic constants, etc..) which is used by the other functions of the package to simulate the photosynthesis. By default, the call of f.make.param() gives the list of parameters used in FATES model for tropical forests.

* The function f.Aci allows to simulate photosynthetic rate using as input variables the leaf temperature, the intercellular CO2, the light intensity.

* The function f.A couples the photosynthetic model with a conductance model (USO, USO_simpl, BWB, ...) and simulates the photosynthetic rate using the air temperature, leaf temperature, light intensity, leaf surface CO2 and humidity. A documentation for simulation and fitting stomatal conductance data is given in [Simulation_of_leaf_conductance](https://github.com/TESTgroup-BNL/LeafGasExchange/tree/master/vignettes/Simulation_of_leaf_conductance.md).

* Eventually, the function f.AT uses, in addition to the stomata model and the photosynthesis model, a leaf energy budget using package tealeaves. It allows to simulate the photosynthetic rates with environment variables as sole inputs (Air temperature, humidity, wind speed and light intensity).
A particular documentation for this function is given in [Simulation_of_leaf_photosynthesis](https://github.com/TESTgroup-BNL/LeafGasExchange/tree/master/vignettes/Simulation_of_leaf_photosynthesis.md). 

Other functions are also implemented to import gas exchange data from LiCor (f.import_licor6400 and f.import_licor6800) and to fit Aci curve and Aq curves (f.fitting). A particular documentation of this function is given in [Aci_fitting](https://github.com/TESTgroup-BNL/LeafGasExchange/tree/master/vignettes/Aci_fitting.md)

A documentation of the equations and parameters used for FATES like photosynthesis module for broadleaf tropical species is presented [here](https://github.com/TESTgroup-BNL/LeafGasExchange/tree/master/Description%20of%20the%20model.pdf). 

# Canopy scaling

The scaling of leaf level gas exchange to the canopy can also be done using the functions f.GPP and f.GPPT. The function f.GPP doesnt include the energy budget. The function f.GPPT includes an energy budget but is slower. Both functions use the Norman 1978 radiation model to descrive the light environment within the canopy. A particular documentation of this function is given in [Canopy_scaling](https://github.com/TESTgroup-BNL/LeafGasExchange/tree/master/vignettes/Canopy_scaling.md).


### Installation
`LeafGasExchange` is not currently on CRAN, but
you can install `LeafGasExchange` from [GitHub](https://github.com/) with:
``` r
devtools::install_github("TESTgroup-BNL/LeafGasExchange", dependencies=TRUE)
```
You will need to install the package devtools first.

### Source Code Citation/DOI
[![DOI](https://zenodo.org/badge/255713515.svg)](https://zenodo.org/badge/latestdoi/255713515)

### This package was used to analyse and simulate gas exchanges in several of our lab studies:

Davidson, K.J., Lamour, J., McPherran, A., Rogers, A. and Serbin, S.P. (2023), Seasonal trends in leaf-level photosynthetic capacity and water use efficiency in a North American Eastern deciduous forest and their impact on canopy-scale gas exchange. New Phytol, 240: 138-156. https://doi.org/10.1111/nph.1913

Davidson, K.J., Lamour, J., Rogers, A., Ely, K.S., Li, Q., McDowell, N.G., Pivovaroff, A.L., Wolfe, B.T., Wright, S.J., Zambrano, A. and Serbin, S.P., 2023. Short‐term variation in leaf‐level water use efficiency in a tropical forest. New Phytologist, 237(6), pp.2069-2087.

Lamour, J., Davidson, K.J., Ely, K.S., Le Moguédec, G., Leakey, A.D., Li, Q., Serbin, S.P. and Rogers, A., 2022. An improved representation of the relationship between photosynthesis and stomatal conductance leads to more stable estimation of conductance parameters and improves the goodness‐of‐fit across diverse data sets. Global change biology, 28(11), pp.3537-3556.

Lamour, J., Davidson, K.J., Ely, K.S., Le Moguédec, G., Anderson, J.A., Li, Q., Calderón, O., Koven, C.D., Wright, S.J., Walker, A.P. and Serbin, S.P., 2023. The effect of the vertical gradients of photosynthetic parameters on the CO2 assimilation and transpiration of a Panamanian tropical forest. New Phytologist.

