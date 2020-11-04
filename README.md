# LeafGasExchange
An R package for fitting and simulating leaf level gas exchange.

This package works using 4 main functions. 
The function f.make.param allows to build a list of photosynthetic parameters (Vcmax, Jmax, kinetic constants, etc..) which is used by the other functions of the package to simulate the photosynthesis. By default, the call of f.make.param() gives the list of parameters used in FATES model for tropical forests. Parameters from other TBMs are also implemented such as CLM4.5, ORCHIDEE and JULES. The choice of a TBM modifies the list of parameters but also the equations that are used to simulate the photosynthesis. For example, ORCHIDEE takes into account the mesophyll conductance while other TBM don't. 

The function f.Aci allows to simulate photosynthetic rate using as input variables the leaf temperature, the intracellular CO[2], the light intensity.

The function f.A couples the photosynthetic model with a conductance model (USO or USO_simpl) and simulates the photosynthetic rate using the air temperature, leaf temperature, light intensity, leaf surface CO2 and humidity.

Eventually, the function f.AT uses, in addition to the stomata model and the photosynthesis model, a leaf energy budget using package tealeaves. It allows to simulate the photosynthetic rates with environment variables as sole inputs (Air temperature, humidity, wind speed and light intensity).
A particular documentation for this function is given in [Simulation_of_leaf_photosynthesis](https://github.com/TESTgroup-BNL/LeafGasExchange/tree/master/vignettes/Simulation_of_leaf_photosynthesis.md). 

Other functions are also implemented to import gas exchange data from LiCor (f.import_licor6400 and f.import_licor6800) and to fit Aci curve and Aq curves (f.fitting). A particular documentation of this function is given in [Aci_fitting](https://github.com/TESTgroup-BNL/LeafGasExchange/tree/master/vignettes/Aci_fitting.md)

# Scaling the leaf gas exchange at the canopy scale

Functions to scale  the leaf gas exchange from the leaf level to the canopy scale are also present in the package. They use notably functions from the package BioCro. A presentation of those functions is given in [Canopy_scaling](https://github.com/TESTgroup-BNL/LeafGasExchange/tree/master/vignettes/Canopy_scaling.md)

### Installation
`LeafGasExchange` is not currently on CRAN, but
you can install `LeafGasExchange` from [GitHub](https://github.com/) with:
``` r
devtools::install_github("TESTgroup-BNL/LeafGasExchange", dependencies=TRUE)
```
