---
title: "Guide to use the LeagGasExchange package to simulate leaf conductance"
author: "Julien LAMOUR"
date: "2/18/2021"
output: github_document
---





## Simulation of simple conductance values

The "BWB" model (Ball et al. 1987) is:
gs=g0+g1*A*RH/CO2s
With g0 the intercept for A == 0, and g1 the slope parameter


```r
A=-2:20
cs=400
RH=70
g0=0.02
g1_BWB=6
BWB_gs=f.gs(A = A,cs = cs,RH = RH,g0 = g0,g1 = g1_BWB,model='BWB')
```

The "USO" model (Medlyn et al. 2011) is: 

gs = g0+ 1.6 * (1 + g1 * A / (ds^0.5 * CO2s)) 

It can be simplified in the "USO_simpl" model :  

gs = g0+ 1.6 * g1 *A / (ds^0.5 * CO2s)


```r
ds=1000
USO_gs=f.gs(A = A,cs = cs,ds = ds,g0 = g0,g1 = g1_BWB,model='USO')
```


```r
USO_simpl_gs=f.gs(A = A,cs = cs,ds = ds,g0 = g0,g1 = g1_BWB,model='USO_simpl')
```

Finally, a "Nonlinear" version of the USO_simpl conductance model is implemented: 

gs = g0 + g1 * 1.6 * (A + Rd)^2 / (CO2s * ds^0.5)


```r
Rd=2
g1_Nonlinear=0.23
Nonlinear_gs=f.gs(A = A,cs = cs,ds = ds,g0 = g0,g1 = g1_Nonlinear, Rd=Rd,model='Nonlinear')
```

The parameters g0 and g1 can be estimated using linear regressions:


```r
regressor_BWB=A*RH/100/cs
regressor_USO_simpl=1.6*A/cs/sqrt(ds/1000)
regressor_Nonlinear=1.6*(A+Rd)^2/sqrt(ds/1000)/cs

lm(BWB_gs~regressor_BWB)
```

```
## 
## Call:
## lm(formula = BWB_gs ~ regressor_BWB)
## 
## Coefficients:
##   (Intercept)  regressor_BWB  
##          0.02           6.00
```

```r
lm(USO_simpl_gs~regressor_USO_simpl)
```

```
## 
## Call:
## lm(formula = USO_simpl_gs ~ regressor_USO_simpl)
## 
## Coefficients:
##         (Intercept)  regressor_USO_simpl  
##                0.02                 6.00
```

```r
lm(Nonlinear_gs~regressor_Nonlinear)
```

```
## 
## Call:
## lm(formula = Nonlinear_gs ~ regressor_Nonlinear)
## 
## Coefficients:
##         (Intercept)  regressor_Nonlinear  
##                0.02                 0.23
```
For the non simplified USO model, it is necessary to change change the variables to be able to perform a linear regression. Indeed: 

gs - 1.6 * A/ (ds^0.5 * CO2s)) = g0+ 1.6 * g1 * A / (ds^0.5 * CO2s)

corresponds to a linear model with Y = gs - 1.6 * A/ (ds^0.5 * CO2s)), and X = 1.6 * A / (ds^0.5 * CO2s)


```r
Y_USO= USO_gs - 1.6 * A/cs
regressor_USO=regressor_USO_simpl

lm(Y_USO~regressor_USO)
```

```
## 
## Call:
## lm(formula = Y_USO ~ regressor_USO)
## 
## Coefficients:
##   (Intercept)  regressor_USO  
##          0.02           6.00
```

## Simulation of conductance using a coupled photosynthesis - conductance model

In the previous examples, the conductance was simulated using the photosynthesis (A) as an input variable. However, in reality both the conductance and the photosynthesis are linked and influence each other. It is possible to simulate the photosynthesis and the conductance together using the f.A function which simulates the leaf gas exchange.



```r
PFD=0:2000
USO_simulation=f.A(PFD = PFD,cs = 400,Tleaf = 25+273.16,Tair = 25+273.16,RH = 70,param = f.make.param(g0=0.02,g1=3,model.gs = "USO"))
USO_simpl_simulation=f.A(PFD = PFD,cs = 400,Tleaf = 25+273.16,Tair = 25+273.16,RH = 70,param = f.make.param(g0=0.02,g1=2.67,model.gs = "USO_simpl"))
BWB_simulation=f.A(PFD = PFD,cs = 400,Tleaf = 25+273.16,Tair = 25+273.16,RH = 70,param = f.make.param(g0=0.02,g1=5,model.gs = "BWB"))
Nonlinear_simulation=f.A(PFD = PFD,cs = 400,Tleaf = 25+273.16,Tair = 25+273.16,RH = 70,param = f.make.param(g0=0.02,g1=0.23,model.gs = "Nonlinear",VcmaxRef=55,RdRef=0.015*55,JmaxRef=1.67*55,TpRef = 20,TBM = "CLM4.5"))
plot(x=PFD,y=USO_simulation$gs,type='l',ylab=expression(Conductance~ mol~m^-2~s^-1),col='lightblue')
lines(x=PFD,y=USO_simpl_simulation$gs,col='slateblue3')
lines(x=PFD,y=BWB_simulation$gs,col='orchid1')
lines(x=PFD,y=Nonlinear_simulation$gs,col='deeppink2')
legend('bottomright',col = c('lightblue','slateblue3','orchid1','deeppink2'),lty = c(1,1,1,1),legend=c("USO","USO_simpl","BWB","Nonlinear"))
```

![plot of chunk unnamed-chunk-7](Simulation_of_leaf_conductance_files/unnamed-chunk-7-1.png)
The same figure is now performed for the transpirated water:


```r
plot(x=PFD,y=USO_simulation$Transp,type='l',ylab=expression(Transpiration~ mL~m^-2~s^-1),col='lightblue')
lines(x=PFD,y=USO_simpl_simulation$Transp,col='slateblue3')
lines(x=PFD,y=BWB_simulation$Transp,col='orchid1')
lines(x=PFD,y=Nonlinear_simulation$Transp,col='deeppink2')
legend('bottomright',col = c('lightblue','slateblue3','orchid1','deeppink2'),lty = c(1,1,1,1),legend=c("USO","USO_simpl","BWB","Nonlinear"))
```

![plot of chunk unnamed-chunk-8](Simulation_of_leaf_conductance_files/unnamed-chunk-8-1.png)



## References
Ball, J. T., Woodrow, I. E., & Berry, J. A. (1987). A model predicting stomatal conductance and its contribution to the control of photosynthesis under different environmental conditions. In Progress in photosynthesis research (pp. 221-224). Springer, Dordrecht.

Medlyn, B.E., Duursma, R.A., Eamus, D., Ellsworth, D.S., Prentice, I.C., Barton, C.V.M., Crous, K.Y., De Angelis, P., Freeman, M. and Wingate, L. (2011), Reconciling the optimal and empirical approaches to modelling stomatal conductance. Global Change Biology, 17: 2134-2144. 

