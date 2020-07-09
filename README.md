# LeafGasExchange
An R package for fitting and simulating leaf level gas exchange.

This package works using 4 main functions. 
The function f.make.param allows to build a list of photosynthetic parameters (Vcmax, Jmax, kinetic constants, etc..) which is used by the other functions of the package to simulate the photosynthesis. By default, the call of f.make.param() give the list of parameters used in FATES model for tropical forests. Parameters from other TBMs are also implemented such as CLM4.5, ORCHIDEE and JULES. The choice of a TBM modifies the list of parameters but also the equations that are used to simulate the photosynthesis. For example, ORCHIDEE takes into account the mesophyll conductance while other TBM don't. 

The function f.Aci allow to simulates photosynthetic rate using as input variables the leaf temperature, the intracellular CO2, the light intensity.

The function f.A couples the photosynthetic model with a conductance model (USO or USO_simpl) and simulate the photosynthetic rate using the air temperature, leaf temperature, light intensity, leaf surface CO2 and humidity.

Eventually, the function f.AT use a leaf energy budget (using package tealeaves) and allow to simulate the photosynthetic rates knowing only environment variables. (Air temperature, humidity, wind speed and light intensity)

Other functions are also implemented to import gas exchange data from LiCor (f.import_licor6400 and f.import_licor6800) and to fit Aci curve and Aq curves (f.fitting)

