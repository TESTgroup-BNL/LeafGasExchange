Aci\_fitting
================
Julien LAMOUR
10/16/2020

# Simulating an Aci curve

For this exemple we first simulate a photosynthesis curve, but it would
work the same if the data were not simulated but measured. The data
simulation is done using the function f.A. This function needs a list of
photosynthetic parameters which are produced using the function
f.make.param() and a list of input variables (CO2 at the surface of the
leaf, leaf temperature, incdent light, RH)

``` r
param=f.make.param(VcmaxRef = 50,JmaxRef=50*1.7,TpRef=50/10)
CO2=seq(50,1500,50)
Tleaf=30+273.16
Tair=27+273.16
PAR=1800
RH=80
simul=f.A(PFD = PAR,cs = CO2,Tleaf = Tleaf,Tair = Tair,RH = RH,param = param)
simul$A=simul$A+rnorm(n = length(CO2),mean = 0,sd = 0.4)
measures=data.frame(Tleaf=Tleaf,Ci=simul$ci,PARi=PAR,Photo=simul$A)
```

We display this simulated curve using the function f.plot

``` r
f.plot(measures = measures,type = 'Aci',list_legend = param[c('VcmaxRef','JmaxRef','TpRef','RdRef')],param = param)
```

![](Aci_fitting_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

\#Fitting an Aci curve

To fit an Aci curve, it is necessary to detail the parameter that we
want to estimate. All the parameters present in f.make.param can
potentially be fitted even if it would not always make sense. We do a
first fitting with only the parameters VcmaxRef, JmaxRef and RdRef.
Those parameters have to be given in the list Start, with initial
values. The method will look for different initial values around those
values so it is not necessary to give very good ones, just not too
stupid ones. The photosynthetic parameters have to be given in the list
param. This is used to determine what should be the parameters for the
temperature dependence, for the leaf absorbance, theta, etc.

``` r
fitting1=f.fitting(measures = measures,Start = list(JmaxRef = 30, VcmaxRef = 50, RdRef = 1),param=f.make.param(),modify.init=TRUE,do.plot=TRUE,type='Aci')
```

    ## $par
    ##   JmaxRef  VcmaxRef     RdRef 
    ## 79.685854 50.740471  1.380129 
    ## 
    ## $value
    ## [1] 13.46319
    ## 
    ## $counts
    ## function gradient 
    ##      150       NA 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## NULL
    ## 
    ## [1] "sd 0.66990512831346"

![](Aci_fitting_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

    ## Length  Class   Mode 
    ##      1   mle2     S4

![](Aci_fitting_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->
In a second example we now also fit TpRef

``` r
fitting2=f.fitting(measures = measures,Start = list(JmaxRef = 30, VcmaxRef = 50, RdRef = 1,TpRef=9),param=f.make.param(),modify.init=TRUE,do.plot=TRUE,type='Aci')
```

    ## $par
    ##   JmaxRef  VcmaxRef     RdRef     TpRef 
    ## 88.761768 52.172765  1.700766  5.013401 
    ## 
    ## $value
    ## [1] 3.940224
    ## 
    ## $counts
    ## function gradient 
    ##      343       NA 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## NULL
    ## 
    ## [1] "sd 0.36240969890763"

![](Aci_fitting_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

    ## Length  Class   Mode 
    ##      1   mle2     S4

![](Aci_fitting_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

The fitting corresponds to a list of 2 objects. The first object
corresponds to the fitting using a minimum square function whereas the
second object corresponds to a maximum likelihood which is made using
the package mle2. This latter method is usefull because it allowed us to
calculate the confidence interval of the parameters:

``` r
confint(fitting2[[2]])
```

    ##              2.5 %     97.5 %
    ## sigma     0.287082  0.4777669
    ## JmaxRef  84.935260 93.7495123
    ## VcmaxRef 49.149819 55.5009942
    ## TpRef     4.849933  5.1790073
    ## RdRef     1.254726  2.1588684
