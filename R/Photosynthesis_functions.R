#' @title Leaf to air water vapour pressure deficit calculation
#' @description This function calculates the leaf water pressure deficit (VPDleaf or Ds) using the temperature of the leaf, the temperature of the air and its relative humidity. The saturation vapor pressure of water is approximated usung the Tetens 1930 equation. 
#' @param Tleaf Temperature of the leaf in Kelvin.
#' @param Tair Temperature of the air in Kelvin.
#' @param RH Humidity of the air (0 to 100).
#' @references Tetens O. 1930. Uber einige meteorologische Begriffe. Z. geophys 6: 297-309.
#' @return Ds in Pascal
#' @export
#' @examples f.ds(Tleaf=273.16 + 30, Tair=273.16+28, RH=70)
f.ds<-function(Tleaf,Tair,RH){
  ds=(0.6108*exp(17.27*(Tleaf-273.16)/(Tleaf-273.16+237.3))-0.6108*exp(17.27*(Tair-273.16)/(Tair-273.16+237.3))*RH/100)*1000
  return(ds)
}


#' @title Conductance model for stomata conductance to water vapour
#' @description Semi-empirical model of the leaf conductance to water vapour
#' @param A Net assimilation in micromol.m-2.s-1, i-e, the assimilation in presence of respiration.
#' @param cs CO2 at the surface of the leaf in ppm.
#' @param ds Leaf surface to air vapour pressure deficit in Pa. Needed for the USO, USO_simpl and Nonlinear models.
#' @param Rd Dark respiration for the non linear model. Be carefull that the net assimilation is implemented in all the models except the non linear where the gross assimilation is used (A + Rd).
#' @param RH Humidity at the surface of the leaf (0 - 100). Only needed for the BWB model.
#' @param g0 Constant of the USO model, representing the conductance when A is 0, in mol.m-2.s-1
#' @param g1 Slope parameter, between 1.14 and 3.58 KPa^0.5 (Wu et al., 2019)
#' @param power Power of the VPDl in USO model. By default is is 0.5 as in Medlyn publication
#' @param model Stomatal model ("USO", "USO_simpl" or "BWB" or "Nonlinear")
#' @export
#' @details 
#' USO : gs=g0+1.6*(1+g1/(ds/1000)^power)*(A)/cs. Cf Medlyn et al. 2011
#' 
#' USO_simpl : gs=g0+1.6*(g1/(ds/1000)^power)*(A)/cs. Cf Medlyn et al. 2011
#' 
#' BWB : gs=g0+g1*(A*RH/100)/cs. Cf Ball et al. 1987
#' 
#' Nonlinear: gs=g0+1.6*(g1/(ds/1000)^power)*(A+Rd)^2/cs. Cf Lamour et al. 2022
#' 
#' @return This function returns the stomatal conductance to water vapour in mol.m-2.s-1
#' @references
#'  Medlyn, B.E., Duursma, R.A., Eamus, D., Ellsworth, D.S., Colin Prentice, I., Barton, C.V.M., Crous, K.Y., de Angelis, P., Freeman, M. and Wingate, L. (2012), Reconciling the optimal and empirical approaches to modelling stomatal conductance. Glob Change Biol, 18: 3476-3476. doi:10.1111/j.1365-2486.2012.02790.x.
#'  
#'  Wu, J, Serbin, SP, Ely, KS, et al. The response of stomatal conductance to seasonal drought in tropical forests. Glob Change Biol. 2020; 26: 823– 839. https://doi.org/10.1111/gcb.14820.
#'  
#'  Leuning, R., Kelliher, F. M., De Pury, D. G. G., & Schulze, E. D. (1995). Leaf nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to canopies. Plant, Cell & Environment, 18(10), 1183-1200.
#'  
#'  Ball, J. T., Woodrow, I. E., & Berry, J. A. (1987). A model predicting stomatal conductance and its contribution to the control of photosynthesis under different environmental conditions. In Progress in photosynthesis research (pp. 221-224). Springer, Dordrecht.
#'  
#'  Lamour, J., Davidson, K.J., Ely, K.S., Le Moguédec, G., Leakey, A.D., Li, Q., Serbin, S.P. and Rogers, A., 2022. An improved representation of the relationship between photosynthesis and stomatal conductance leads to more stable estimation of conductance parameters and improves the goodness‐of‐fit across diverse data sets. Global change biology, 28(11), pp.3537-3556.
#' @examples f.gs(A=30,cs=400,ds=1500,g0=0.01,g1=2,power=0.5)
f.gs<-function(A,cs,ds=NULL,RH=NULL,Rd=NULL,g0,g1,power=0.5,model="USO"){
  if(model=="USO"|model==0){
    gs=g0+1.6*(1+g1/(ds/1000)^power)*(A)/cs
  } else if(model=="USO_simpl"|model==1){
    gs=g0+1.6*(g1/(ds/1000)^power)*(A)/cs
  } else if(model=="BWB"|model==2){
    gs=g0+g1*(A*RH/100)/cs
  } else if(model=="Nonlinear"|model==3){
    gs=g0+1.6*(g1/(ds/1000)^power)*(A+Rd)^2/cs
  } 
  else{print(paste("Parameter model =",model,"is not in the list of implemented models"))}
  return(gs)
}

#' @title Calculation of the minimal conductance given by a particular coupled conductance and photosynthesis model
#' @description The minimal conductance of a model depends on the parameters of the model (ie g0 and g1) but also on the minimum A value, which corresponds to the dark respiration.
#' Knowing the minimal conductance is important because the conductance can become negative and lead to unrealistic values in photosynthesis models
#' @inheritParams f.make.param
#' @param RH Air relative humidity (0 to 100), only necessary for the BWB stomatal conductance model
#' @param cs CO2 concentration at the leaf surface in ppm
#' @param Tleaf Leaf temperature in Kelvin
#' @param ds Leaf to air vapor pressure deficit in Pa
#'
#' @return Minimum conductance in mol m-2 s-1
#' @export
#'
#' @examples f.gsmin(RdRef=	0.825,RdHa=	46390,RdHd=150650,RdS=490,Tleaf=300,cs=400,ds=1000,g0=0.02,g1=4.1,power=0.5,model="USO")
f.gsmin<-function(RdRef=	0.825,RdHa=	46390,RdHd=150650,RdS=490,Tleaf=300,cs=400,ds=1000,g0=0.02,g1=4.1,power=0.5,model="USO"){
  Rd=f.modified.arrhenius(PRef=RdRef,Ha = RdHa,Hd = RdHd,s = RdS,Tleaf = Tleaf)
  if(model=="USO"|model==0){
    gsmin=g0+1.6*(1-g1/(ds/1000)^power)*(Rd)/cs
  } else if(model=="USO_simpl"|model==1){
    gsmin=g0-1.6*g1*Rd/(cs*(ds/1000)^power)
  } else if(model=="BWB"|model==2){
    gsmin=g0-g1*(Rd*RH/100)/cs
  } else if(model=="Nonlinear"|model==3){
    gsmin=g0
  } 
  else{print(paste("Parameter model =",model,"is not in the list of implemented models"))}
  return(gsmin)
  
}

#' @title Maximum theoretical stomatal conductance
#' @details This function calculates the maximum theoretical conductance value according to morphological data and the physics of diffusion through pores. It follows the equation from Franks and Berling 2009.
#' @param Sarea Maximum area of the aperture of stomata when open (microm2).
#' @param Sdensity Number of stomata per mm2 of leaf.
#' @param Sdepth Stomatal pore depth (micro m).
#' @param Diffusivity Diffusivity of water vapor in air (m2 s-1).
#' @param mvair Molar volume of air (m3 mol-1).
#'
#' @return Maximum stomatal conductance to water vapour in mol m-2 s-1
#' @export
#' @references Franks PJ, Beerling DJ. Maximum leaf conductance driven by CO2 effects on stomatal size and density over geologic time. Proc Natl Acad Sci U S A. 2009;106(25):10343-10347. doi:10.1073/pnas.0904209106
#' @examples
#' ## The density of stomata is around 400 stomata.mm-2 in the tropical species.
#' ## The length of the stomata is around 20 micro m. Following Franks and Beerling 2009 we can estimate the
#' ## Sarea of the stomata: pi*(20/4*10^-6)^2 and the Sdepth: 20*10^-6/4.
#' f.gsmax(Sarea=0.78,Sdensity=400,Sdepth=5)
f.gsmax=function(Sarea=0.78,Sdensity=400,Sdepth=5,Diffusivity=0.282/1000,mvair=24.5/1000){
  Diffusivity/mvair*Sdensity*Sarea*10^-6/(Sdepth*10^-6+pi/2*(Sarea*10^-12/pi)^0.5)
}


#' @title Temperature dependence of Gamma star, Ko, Kc and Rd
#'
#' @param PRef Value of the parameter at the reference temperature.
#' @param TRef Reference temperature in Kelvin.
#' @param Ha Enthalpie of activation in J.mol-1.
#' @param Tleaf Temperature of the leaf in Kelvin.
#' @param R Ideal gas constant.
#'
#' @return Value of the parameter at the temperature of the leaf
#' @export
#' @references VON CAEMMERER, S. (2013), Steady state models of photosynthesis. Plant Cell Environ, 36: 1617-1630. doi:10.1111/pce.12098.
#'  Bernacchi, C.J., Singsaas, E.L., Pimentel, C., Portis Jr, A.R. and Long, S.P. (2001), Improved temperature response functions for models of Rubisco‐limited photosynthesis. Plant, Cell & Environment, 24: 253-259. doi:10.1111/j.1365-3040.2001.00668.x.
#' @examples plot(x=seq(25,35,0.1),y=f.arrhenius(PRef=1,Ha=46390,Tleaf=seq(273.15+25,273.15+35,0.1),R=8.314),xlab='Temperature degree C',ylab='Rd')

f.arrhenius<-function(PRef,Ha,Tleaf,TRef=298.16,R=8.314){
  P=PRef*exp(Ha/(R*TRef)-Ha/(R*Tleaf))
  return(P)
}

#' @title Temperature dependence of Gamma star, Ko, Kc and Rd
#' @details Retrieve the value of the parameter at TRef knowing its value at Tleaf
#' @inheritParams f.arrhenius
#' @param P Value of the parameter at Tleaf.
#' @return Value of the parameter at the reference temperature
#' @export
#'
#' @examples RdTleaf=f.arrhenius(PRef=1,Ha=46390,Tleaf=273.16+32,R=8.314)
#' f.arrhenius.inv(P=RdTleaf,Ha=46390,Tleaf=273.16+32,R=8.314)

f.arrhenius.inv<-function(P,Ha,Tleaf,TRef=298.16,R=8.314){
  PRef=P/exp(Ha/(R*TRef)-Ha/(R*Tleaf))
  return(PRef)
}

## @title Temperature dependence of photosynthetic parameters
## @details This equation is used in JULES TBM model
## @inheritParams f.make.param
##
## @return Value of the photosynthetic parameter at the specified leaf temperature
## @export
## @references Clark, D. B., Mercado, L. M., Sitch, S., Jones, C. D., Gedney, N., Best, M. J., . Cox, P. M. (2011). The Joint UK Land Environment Simulator (JULES), model description - Part 2: Carbon fluxes and vegetation dynamics. Geoscientific Model Development, 4(3), 701-722. doi:10.5194/gmd-4-701-2011
## @examples
#f.Q10=function(Pref,Q10,Tleaf,TRef){
#  P=Pref*Q10^(0.1*(Tleaf-TRef))
#  return(P)
#}
#
## @title Temperature dependence of photosynthetic parameters
## @details This equation is used in JULES TBM model
##
## @inheritParams f.make.param
##
## @return Value of the photosynthetic parameter at the specified leaf temperature
## @references Clark, D. B., Mercado, L. M., Sitch, S., Jones, C. D., Gedney, N., Best, M. J., . Cox, P. M. (2011). The Joint UK Land Environment Simulator (JULES), model description - Part 2: Carbon fluxes and vegetation dynamics. Geoscientific Model Development, 4(3), 701-722. doi:10.5194/gmd-4-701-2011
##
## @examples
#f.Q10.modified=function(Pref,Q10,Tleaf,TRef,Tlow,Tup){
#  P=Pref*Q10^(0.1*(Tleaf-TRef))/((1+exp(0.3*(Tleaf-Tup)))*(1+exp(0.3*(Tlow-Tleaf))))
#  return(P)
#}

#' @title Temperature dependence of Jmax and Vcmax
#' @description The temperature dependence of the photosynthetic parameters Vcmax, the maximum catalytic rate of the enzyme Rubisco, and Jmax, the maximum electron transport rate is modeled by a modified Arrhenius equation. It is modified to account for decreases in each parameter at high temperatures.
#' @param PRef Value of the parameter, here Vcmax or Jmax, at the reference temperature in micromol.m-2.s-1.
#' @param Ha Energy of activation in J.mol-1.
#' @param Hd Energy of desactivation in J.mol-1.
#' @param s Entropy term in J.mol-1.K-1.
#' @param Tleaf Temperature of the leaf in Kelvin.
#' @param TRef Reference temperature in Kelvin, usually 25 + 273.16 K.
#' @param R Ideal gas constant.
#'
#' @return Value of the parameter Jmax or Vcmax at a given temperature
#' @references Leuning, R. (2002), Temperature dependence of two parameters in a photosynthesis model. Plant, Cell & Environment, 25: 1205-1210. doi:10.1046/j.1365-3040.2002.00898.x.
#' @export
#'
#' @examples plot(x=seq(25,35,0.1),y=f.modified.arrhenius(PRef=50,Ha=73637,Hd=149252,s=486,Tleaf=seq(273.15+25,273.15+35,0.1)),xlab='Temperature degree C',ylab='Vcmax')
f.modified.arrhenius<-function(PRef,Ha,Hd,s,Tleaf,TRef=298.16,R=8.314){
  P=PRef*(1+exp((s*TRef-Hd)/(R*TRef)))*exp(Ha/(R*TRef)*(1-TRef/Tleaf))/(1+exp((s*Tleaf-Hd)/(R*Tleaf)))
  return(P)
}

#' @title Temperature dependence of Jmax and Vcmax
#' @description Retrieve the reference temperature value of a parameter knowing its value at Tleaf
#' @param P Value of the parameter, here Vcmax or Jmax, at the leaf temperature in micromol.m-2.s-1
#' @inheritParams f.modified.arrhenius
#'
#' @return Value of the parameter Jmax or Vcmax at the reference temperature
#' @export
#'
#' @examples VcmaxTleaf=f.modified.arrhenius(PRef=50,Ha=73637,Hd=149252,s=486,Tleaf=seq(273.15+25,273.15+35,0.1))
#' f.modified.arrhenius.inv(P=VcmaxTleaf,Ha=73637,Hd=149252,s=486,Tleaf=seq(273.15+25,273.15+35,0.1))
f.modified.arrhenius.inv<-function(P,Ha,Hd,s,Tleaf,TRef=298.16,R=8.314){
  PRef=P/(1+exp((s*TRef-Hd)/(R*TRef)))/exp(Ha/(R*TRef)*(1-TRef/Tleaf))*(1+exp((s*Tleaf-Hd)/(R*Tleaf)))
  return(PRef)
}

#' @title Coupled conductance photosynthesis model
#' @description Photosynthesis model at the leaf level using the farquhar equations coupled with a conductance model. The parameters can be defined by the
#' function f.make param and corresponds to the parameters implemented in FATES model. Note that the conductance model is selected in the param parameter.
#' @inheritParams f.make.param
#' @param PFD Photosynthetic light at the leaf surface in micro mol m-2 s-1.
#' @param cs CO2 concentration at the leaf surface.
#' @param Tleaf Leaf temperature in Kelvin.
#' @param Tair Air temperature in Kelvin.
#' @param RH Relative humidity at the leaf surface (between 0 and 100).
#' @param param List of parameters given by f.make.param().
#' @export
#' @references
#' Farquhar, G. D., von Caemmerer, S. V., & Berry, J. A. (1980). A biochemical model of photosynthetic CO 2 assimilation in leaves of C 3 species. planta, 149, 78-90.
#' Fates parameters: https://github.com/NGEET/fates/blob/main/parameter_files/fates_params_default.cdl.
#' @return 
#' List of different variables:
#'  - A: Raw assimilation of the leaf in micromol.m-2.s-1.
#'  
#'  - Ac: Rubisco limitation assimilation of the leaf in micromol.m-2.s-1.
#'  
#'  - Aj: Electron transport rate assimilation of the leaf in micromol.m-2.s-1.
#'  
#'  - Ap: TPU rate of the leaf in micromol.m-2.s-1.
#'  
#'  - Ag: Gross assimilation in micromol.m-2.s-1.
#'  
#'  - Rd: Respiration rate in micromol.m-2.s-1.
#'  
#'  - gs: Conductance of the leaf for water vapour in mol m-2 s-1.
#'  
#'  - ci: Intracellular CO2 concentration in micromol.mol-1.
#'  
#'  - ds: Leaf surface to air vapour pressure deficit in Pa.
#'  
#'  - Transp: Water transpiration in mL m-2 s-1.
#'
#' @examples f.A(PFD=2000,cs=400,Tleaf=273.16+29,Tair=273.16+28,RH=70,param=f.make.param())
f.A<-function(PFD,cs,Tleaf,Tair,RH,param=f.make.param()){
  #Calculation of temperature dependence of the parameters
  Kc=f.arrhenius(param[['KcRef']],param[['KcHa']],Tleaf)
  Ko=f.arrhenius(param[['KoRef']],param[['KoHa']],Tleaf)
  Gstar=f.arrhenius(param[['GstarRef']],param[['GstarHa']],Tleaf)
  Rd=f.modified.arrhenius(PRef=param[['RdRef']],param[['RdHa']],param[['RdHd']],param[['RdS']],Tleaf)
  Vcmax=f.modified.arrhenius(PRef=param[['VcmaxRef']],param[['VcmaxHa']],param[['VcmaxHd']],param[['VcmaxS']],Tleaf)
  Jmax=f.modified.arrhenius(PRef=param[['JmaxRef']],param[['JmaxHa']],param[['JmaxHd']],param[['JmaxS']],Tleaf)
  
  # Calculation of the electron transport rate
  I2=PFD*param[['abso']]*(param[['aQY']])
  J=(I2+Jmax-((I2+Jmax)^2-4*(param[['Theta']])*I2*Jmax)^0.5)/(2*(param[['Theta']]))
  
  # Calculation of the leaf to air vapor pressure deficit
  ds=f.ds(Tleaf,Tair,RH)
  
  # Analytical solution of the system of equations {E1 : A=f(ci), E2 : gs=f(A,cs) and ci=f(cs)}
  cic=f.solv(x=Vcmax,y=Kc*(1+param[['O2']]/Ko),cs=cs,Rd=Rd,Gstar=Gstar,g0=param[['g0']],g1=param[['g1']],power=param[['power']],ds=ds,RH=RH,model=param[["model.gs"]])
  Wc=(cic-Gstar)*Vcmax/(cic+Kc*(1+param[['O2']]/Ko))
    
  cij=f.solv(x=J/4,y=2*Gstar,cs=cs,Rd=Rd,Gstar=Gstar,g0=param[['g0']],g1=param[['g1']],param[['power']],ds=ds,RH=RH,model=param[["model.gs"]])
  #plot(cij[,1],ylim=c(-400,400))
  #points(cij[,2],col="blue")
  #points(cij[,3],col="red")
    
  Wj=(cij-Gstar)*J/(4*cij+8*Gstar)
    
  Tp=f.modified.arrhenius(PRef=param[['TpRef']],param[['TpHa']],param[['TpHd']],param[['TpS']],Tleaf)
  Wp=3*Tp
  cip=f.solv(x=3*Tp,y=-Gstar,cs=cs,Rd=Rd,Gstar=Gstar,g0=param[['g0']],g1=param[['g1']],power=param[['power']],ds=ds,RH=RH,model=param[["model.gs"]])
    
  ci=cij
  if(!is.null(which(Wc<Wj))&length(cic)==length(cij)){ci[which(Wc<Wj)]=cic[which(Wc<Wj)]}
  if(!is.null(which(Wc<Wj))&length(cic)!=length(cij)){ci[which(Wc<Wj)]=cic}
  W=pmin(Wc,Wj)
  if(!is.null(which(Wp<W))){ci[which(Wp<W)]=cip[which(Wp<W)]}
  Ai=f.smooth(A1 = Wc,A2 = Wj,theta=param[['thetacj']])
  A=f.smooth(A1=Ai,A2=Wp,theta=param[['thetaip']])-Rd
  gs=f.gs(A=A,cs=cs,ds=ds,Rd=Rd,RH=RH,g0=param[['g0']],g1=param[['g1']],power=param[['power']],model =param[['model.gs']])
  output=list(A=A,Ac=Wc-Rd,Aj=Wj-Rd,Ap=Wp-Rd,Ag=A+Rd,Rd=Rd,gs=gs,ci=ci,ds=ds,Transp=gs*ds/(param[['Patm']]*1000)*18)
  return(output)
}



#' @title Coupled conductance photosynthesis model and energy balance model
#' @details This function allows to calculate the photosynthesis from environmental variables PFD, RH, wind, cs and Tair.
#' The energy balance model is calculated using the package Tealeaves (see reference). The energy balance calculation involves the stomatal conductance and the cuticular conductance.
#' Here the cuticular conductance is considered to be equal to g0. 
#' @param PFD Photosynthetic light at the leaf surface in micro mol m-2 s-1.
#' @param ca CO2 concentration in the air in ppm.
#' @param RHa Relative humidity in the air (between 0 and 100).
#' @param Tair Air temperature in Kelvin.
#' @param param List of parameters given by f.make.param().
#' @param precision Precision of the leaf temperature prediction. The resolution of the energy balance coupled with the photosynthesis and stomatal conductance is numerical. The smaller the precision, the longer will be the resolution.
#' @param max_it Maximum number of iterations to find the solution.
#' @param wind Wind speed at the surface of the leaf in m.s-1.
#' @param abso_s absorptance of the leaves in shortwave, see documentation of tealeaves package.
#' @param leaf_size Dimension of the leaf in m (see Tealeaves package).
#' @param NIR NIR radiation in watt m-2, if not given, then by default the shortwave radiation is calculated as PFD/4.57+NIR = PFD/(4.57*0.45) (see Yun et al. 2020 for the constants 0.45 and 4.57).
#' @return List of different variables:
#' 
#'  - A: Raw assimilation of the leaf in micromol.m-2.s-1.
#'  
#'  - Ac: Rubisco limitation assimilation of the leaf in micromol.m-2.s-1.
#'  
#'  - Aj: Electron transport rate assimilation of the leaf in micromol.m-2.s-1.
#'  
#'  - Ap: TPU rate of the leaf in micromol.m-2.s-1.
#'  
#'  - Ag: Gross assimilation in micromol.m-2.s-1.
#'  
#'  - Rd: Respiration rate in micromol.m-2.s-1.
#'  
#'  - gs: Conductance of the leaf for water vapour in mol m-2 s-1.
#'  
#'  - ci: Intracellular CO2 concentration in micromol.mol-1.
#'  
#'  - ds: Leaf surface to air vapour pressure deficit in Pa.
#'  
#'  - Transp: Water transpiration in mL m-2 s-1.
#'  
#'  - gbw: boundary layer conductance in mol m-2 s-1.
#'  
#'  - cs: CO2 concentration at the leaf surface in ppm.
#'  
#'  - RHs: Relative humidity at the leaf surface (0 100).
#'  
#'  - Tleaf: Leaf temperature in K
#'  
#' @export
#' @keywords internal
#' @references 
#' tealeaves: an R package for modelling leaf temperature using energy budgets. Christopher. D. Muir. bioRxiv 529487; doi: https://doi.org/10.1101/529487 
#' 
#'Yun, S. H., Park, C. Y., Kim, E. S., & Lee, D. K. (2020). A Multi-Layer Model for Transpiration of Urban Trees Considering Vertical Structure. Forests, 11(11), 1164.
#' @examples 
#' f.ATnotvectorised(PFD=1500,ca=400,Tair=298,wind=2,RHa=70,param=f.make.param(g0=0.03))
f.ATnotvectorised<-function(PFD,ca,Tair,RHa,wind,precision=0.05,max_it=10,param,NIR=NA,abso_s=0.5,leaf_size=0.04){
  #Defining the initial conditions at the leaf surface
  cs=ca
  RHs=RHa
  Tleaf=Tair+1
  ea=0.6108*exp(17.27*(Tair-273.16)/(Tair-273.16+237.3))*RHa/100 ## water pressure in the air (using Tetens 1930 equation, see f.ds doc)
  n=1
  delta=precision+1
  while(delta>precision&n<max_it){
    Leaf_physio=f.A(PFD=PFD,Tleaf=Tleaf,Tair=Tair,cs = cs,RH = RHs,param=param)
    ds=f.ds(Tleaf,Tair,RHs)
    gs=Leaf_physio$gs-param[['g0']] ## The package Tealeaves separates the cuticular and the leaf conductance. Here I consider that the cuticular conductance is g0 which is a (wrong..) approximation
    if(gs<param[['g0']]){gs=param[['g0']]}
    g_sw <- set_units(gs, "mol/m^2/s")  ##Stomatal conductance
    g_uw<- set_units(param[['g0']], "mol/m^2/s")   ##Cuticular conductance
    g_sw=convert_conductance(g_sw,
                             Temp = set_units(Tair, "K"),
                             P = set_units(param[['Patm']], "kPa"))$`umol/m^2/s/Pa`
    g_uw=convert_conductance(g_uw,
                             Temp = set_units(Tair, "K"),
                             P = set_units(param[['Patm']], "kPa"))$`umol/m^2/s/Pa`
    
    leaf_par <- make_leafpar(replace = list(leafsize=set_units(c(leaf_size), "m"),g_sw=g_sw,g_uw=g_uw,abs_s=set_units((abso_s))))
    if(is.na(NIR))(SW=PFD/(4.57*0.45)) else (SW=PFD/4.57+NIR)
    enviro_par <- make_enviropar(replace=list(S_sw=set_units(SW,"W/m^2"),RH=set_units(RHa/100),T_air = set_units(Tair, "K"),wind=set_units(wind,"m/s")))## see A Multi-Layer Model for Transpiration of Urban Trees Considering Vertical Structure for 4.57 and 0.45
    constants <- make_constants()
    leaf_budget=tleaf(leaf_par, enviro_par, constants,quiet = TRUE,set_units = FALSE)
    Tleaf_mod <- as.numeric(leaf_budget$T_leaf)
    gbw=leaf_budget$g_bw*param[['Patm']]*1000/(param[['R']]*(Tair)) ## Converting the boundary layer value from m/s to mol/m-2/s-1
    Leaf_physio$gbw=gbw
    Leaf_physio$cs=ca-1.37/gbw*Leaf_physio$A ## Calculating the new CO2 concentration at the leaf surface using fick s law
    es=(gbw*ea+Leaf_physio$gs*0.6108*exp(17.27*(Tleaf_mod-273.16)/(Tleaf_mod-273.16+237.3)))/(gbw+Leaf_physio$gs)
    Leaf_physio$RHs=es/(0.6108*exp(17.27*(Tleaf_mod-273.16)/(Tleaf_mod-273.16+237.3)))*100
    delta=abs(Tleaf-Tleaf_mod)
    Tleaf=Tleaf_mod
    n=n+1
    if(n>max_it&delta>precision){
      Tleaf=NA
      print(paste("Temperature convergence limit"))}
    Leaf_physio$Tleaf=Tleaf
  }
  return(Leaf_physio)
}


#' @title Coupled conductance photosynthesis model and energy balance model
#' @details This function allows to calculate the photosynthesis from environmental variables PFD, RH, wind, cs and Tair.
#' The energy balance model is calculated using the package Tealeaves (see reference). The energy balance calculation involves the stomatal conductance and the cuticular conductance.
#' Here the cuticular conductance is considered to be equal to g0 as done in some TBMs even if it is probably a wrong representation.This choice was made to prevent unrealistic energy budgets when the conductance is too low (<= 0) for low light levels.
#' @inheritParams f.AT
#' @return 
#' List of different variables:
#' 
#'  - A: Raw assimilation of the leaf in micromol.m-2.s-1.
#'  
#'  - Ac: Rubisco limitation assimilation of the leaf in micromol.m-2.s-1.
#'  
#'  - Aj: Electron transport rate assimilation of the leaf in micromol.m-2.s-1.
#'  
#'  - Ap: TPU rate of the leaf in micromol.m-2.s-1.
#'  
#'  - Ag: Gross assimilation in micromol.m-2.s-1.
#'  
#'  - Rd: Respiration rate in micromol.m-2.s-1.
#'  
#'  - gs: Conductance of the leaf for water vapour in mol m-2 s-1.
#'  
#'  - ci: Intracellular CO2 concentration in micromol.mol-1.
#'  
#'  - ds: Leaf surface to air vapour pressure deficit in Pa.
#'  
#'  - Transp: Water transpiration in mL m-2 s-1.
#'  
#'  - gbw: boundary layer conductance in mol m-2 s-1.
#'  
#'  - cs: CO2 concentration at the leaf surface in ppm.
#'  
#'  - RHs: Relative humidity at the leaf surface (0 100).
#'  
#'  - Tleaf: Leaf temperature in K
#' @export
#' @references tealeaves: an R package for modelling leaf temperature using energy budgets. Christopher. D. Muir. bioRxiv 529487; doi: https://doi.org/10.1101/529487.
#' @examples leaf_physio=f.AT(PFD=seq(0,1500,50),ca=400,Tair=298,wind=2,RHa=70,param=f.make.param(g0=0.03))
#' plot(x=seq(0,1500,50),y=leaf_physio$A)

f.AT<-function(PFD,NIR=NA,ca,Tair,RHa,wind,precision=0.05,max_it=10,param,abso_s=0.5,leaf_size=0.04){
  input_variable=data.frame(PFD=PFD,ca=ca,Tair=Tair,RHa=RHa,wind=wind,NIR=NIR,abso_s=abso_s,leaf_size=leaf_size)
  result=apply(X = input_variable,1,FUN = function(x){
    f.ATnotvectorised(PFD = x['PFD'],ca = x['ca'],Tair = x['Tair'],RHa = x['RHa'],wind = x['wind'],NIR=x['NIR'],abso_s=x['abso_s'],
                      precision = precision,max_it = max_it,param = param, leaf_size=x['leaf_size'])})
  return(as.data.frame(matrix(unlist(result),ncol = length(result[[1]]),byrow = TRUE,dimnames = list(NULL,names(result[[1]])))))
}

#' @title Polynom solver
#' @param data 
#' @return
#' @export
#' @keywords internal
#' @examples
f.solv.poly<-function(data){
  solution=t(apply(X = data,MARGIN = 1,FUN=function(x){solve(polynomial(x))}))
  if(is.list(solution)){
    solution=t(sapply(solution, function(x) x[1:3],simplify = "array"))
  }
  solution[Im(solution)!=0]=NA
  solution=Re(solution)
  return(solution)
}

#' @title Analytical solution of the coupled photosynthesis and conductance models
#' @param x
#' @param y
#' @param cs
#' @param Rd
#' @param Gstar
#' @param g0
#' @param g1
#' @param ds
#' @param model
#'
#' @return
#' @export
#' @keywords internal
#'
#' @examples
f.solv<-function(x,y,cs,Rd,Gstar,g0,g1,power,ds,RH,model){
  if(model=="USO"|model==0){
    m=1.6*(1+g1/(ds/1000)^power)
    if(g1==-1){m=0}
  }else if(model=="USO_simpl"|model==1){
    m=1.6*(g1/(ds/1000)^power)
  }else if(model=="BWB"|model==2){
    m=(g1*RH/100)
  }
  if(model=='Nonlinear'|model==3){
    a1=5*g0*cs*sqrt(ds/1000)+8*g1*x^2
    a2=-16*Gstar*g1*x^2-5*sqrt(ds/1000)*cs^2*g0+10*sqrt(ds/1000)*cs*g0*y-8*cs*g1*x^2-8*sqrt(ds/1000)*Rd*cs+8*sqrt(ds/1000)*cs*x
    a3=8*Gstar^2*g1*x^2+16*Gstar*cs*g1*x^2-10*sqrt(ds/1000)*cs^2*g0*y+5*sqrt(ds/1000)*cs*g0*y^2-8*Gstar*sqrt(ds/1000)*cs*x-16*sqrt(ds/1000)*Rd*cs*y+8*sqrt(ds/1000)*cs*x*y
    a4=-8*Gstar^2*cs*g1*x^2-5*sqrt(ds/1000)*cs^2*g0*y^2-8*Gstar*sqrt(ds/1000)*cs*x*y-8*sqrt(ds/1000)*Rd*cs*y^2
    coef_pol=cbind(a4,a3,a2,a1)
    res=f.solv.poly(coef_pol)
    res[res<0]=NA
    if(dim(res)[2]==2)(result=pmax(res[,1],res[,2],na.rm = TRUE))else(result=pmax(res[,1],res[,2],res[,3],na.rm=TRUE)) ## 23/02/2023 I changed the solution as I had issues
    return(result)
    #return(res)
  } else {
    a=g0+m/cs*(x-Rd)
    b=y*g0+m/cs*(-Gstar*x-Rd*y)-cs*g0+(x-Rd)*(1.6-m)
    c=-y*cs*g0+(1.6-m)*(-Gstar*x-Rd*y)
    ci2=(-b-(b^2-4*a*c)^0.5)/(2*a)
    ci1=(-b+(b^2-4*a*c)^0.5)/(2*a)
    return(pmax(ci1,ci2)) ## Be careful that the max does not always correspond to the more sounding solution at low light
  }
}


#f.solv.Acc=function(x1,x2,g0,g1,power,ds,RH,model,gm,Rd,Gstar,cs){
#  if(model=="USO"|model==0){
#    m=1.6*(1+g1/(ds/1000)^power)
#  }else if(model=="USO_simpl"|model==1){
#    m=1.6*(g1/(ds/1000)^power)
#  }else if(model=="BWB"|model==2){
#    m=(g1*RH/100)
#  }else{
#    print(paste("model:",model,'is not in the list of implemented stomatal models'))
#  }
#  s=x1-Rd
#  t=-Gstar*x1-Rd*x2
#  u=g0+1.6*gm-cs*m*gm
#  a1=gm*(g0+m*s)
#  a2=m*s^2+u*s+t*m*gm+x2*m*gm*s+2*g0*gm*x2-g0*gm*cs
#  a3=2*m*t*s+u*t+x2*u*s+t*m*gm*x2+g0*gm*x2^2-2*cs*g0*gm*x2
#  a4=x2*u*t-cs*g0*gm*x2^2+m*t^2
#  p=a2/a1
#  q=a3/a1
#  r=a4/a1
#  U=(2*p^3-9*p*q+27*r)/54
#  Q=(p^2-3*q)/9
#  phi=acos(U/sqrt(Q^3))
#  #sol1=-2*sqrt(Q)*cos(phi/3)-p/3
#  sol2=-2*sqrt(Q)*cos((phi+2*pi)/3)-p/3
#  #sol3=-2*sqrt(Q)*cos((phi+4*pi)/3)-p/3
#  return(sol2)
#}


#' @title Photosynthesis and stomata model parameters
#' @description Function to create a list of parameters to be used in most of the functions of this package.
#' Depending on the function, all the parameters are not used. For example go and g1 are not used in f.Aci.
#' @details The call of this function is made using f.make.param(). If a parameter is modified, for example writing f.make.param(VcmaxRef=10), this function will return all the default parameters from FATES TBM with VcmaxRef = 10 instead of its default value
#' @param R Ideal gas constant.
#' @param O2 O2 concentration in ppm.
#' @param TRef Reference temperature for Kc, Ko, Rd,GammaStar Vcmax, Jmax in Kelvin.
#' @param Patm Atmospheric pressure in kPa.
#' @param JmaxRef Maximum electron transport rate in micromol.m-2.s-1.
#' @param JmaxHa Energy of activation for Jmax in J.mol-1.
#' @param JmaxHd Energy of desactivation for Jmax in J.mol-1.
#' @param JmaxS Entropy term for Jmax in J.mol-1.K-1.
#' @param VcmaxRef Maximum rate of Rubisco for carboxylation micromol.m-2.s-1.
#' @param VcmaxHa Energy of activation for Vcmax in J.mol-1.
#' @param VcmaxHd Energy of desactivation for Vcmax in J.mol-1.
#' @param VcmaxS Entropy term for Vcmax in J.mol-1.K-1.
#' @param TpRef Triose phosphate utilization  rate in micromol.m-2.s-1.
#' @param TpHa Activation energy for Tp in J.mol-1.
#' @param TpHd Energy of deactivation for Tp in J.mol-1.
#' @param TpS Entropy term for Tp in J.mol-1.K-1.
#' @param thetacj Collatz smoothing factor used to introduce a gradual transition from Ac to Aj (close to 0.999)
#' @param thetaip Collatz smoothing factor used to introduce a gradual transition from Aj to Ap (close to 0.999)
#' @param RdRef Respiration value at the reference temperature in micromol.m-2.s-1.
#' @param RdHa Energie of activation for Rd in J.mol-1.
#' @param KcRef Michaelis-Menten constant of Rubisco for CO2 at the reference temperature in micromol.mol-1.
#' @param KcHa Energy of activation for Kc in J.mol-1.
#' @param KoRef Michaelis-Menten constant of Rubisco for CO2 at the reference temperature in milimol.mol-1.
#' @param KoHa Energy of activation for Ko in J.mol-1.
#' @param GstarRef CO2 compensation point in absence of respiration in micromol.mol-1.
#' @param GstarHa Enthalpie of activation for Gstar in J.mol-1.
#' @param abso Absorptance of the leaf in the photosynthetic active radiation wavelenghts.
#' @param aQY Apparent quantum yield.
#' @param Theta Theta is the empirical curvature factor for the response of J to PFD. It takes its values between 0 and 1.
#' @param model.gs Type of conductance model (USO, USO_simpl,BWB or Nonlinear). See f.gs documentation for more information.
#' @param g0 Constant of the conductance model, representing the conductance when A is 0, in mol.m-2.s-1, usually around 0.01.
#' @param g1 Slope parameter, between 1.14 and 3.58 KPa^0.5 (Wu et al., 2019).
#' @param power Power of VPDl in USO model. By default power=0.5 as in Medlyn article.
#' 
#' @return List of parameters that can be used in other functions of the package such as f.A, f.Aci, and f.GPP
#' 
#' @references 
#' Bernacchi, C.J., Singsaas, E.L., Pimentel, C., Portis Jr, A.R. and Long, S.P. (2001), Improved temperature response functions for models of Rubisco‐limited photosynthesis. Plant, Cell & Environment, 24: 253-259. doi:10.1111/j.1365-3040.2001.00668.
#' 
#' FATES: https://fates-docs.readthedocs.io/en/latest/fates_tech_note.html.
#' 
#' Medlyn, B.E., Duursma, R.A., Eamus, D., Ellsworth, D.S., Colin Prentice, I., Barton, C.V.M., Crous, K.Y., de Angelis, P., Freeman, M. and Wingate, L. (2012), Reconciling the optimal and empirical approaches to modelling stomatal conductance. Glob Change Biol, 18: 3476-3476. doi:10.1111/j.1365-2486.2012.02790.x.
#' 
#' Leuning, R., Kelliher, F. M., De Pury, D. G. G., & Schulze, E. D. (1995). Leaf nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to canopies. Plant, Cell & Environment, 18(10), 1183-1200.
#' 
#' Ball, J. T., Woodrow, I. E., & Berry, J. A. (1987). A model predicting stomatal conductance and its contribution to the control of photosynthesis under different environmental conditions. In Progress in photosynthesis research (pp. 221-224). Springer, Dordrecht.

#' @export
#'
#' @examples 
#' param1=f.make.param(JmaxRef=100,VcmaxRef=60,RdRef=1,TpRef=10)
#' f.A(PFD=1500,cs=400,Tleaf=300,Tair=299,RH=70,param=param1)
f.make.param<-function(R=NA,O2=NA,TRef=NA,
                       Patm=NA,JmaxRef=	NA,
                       JmaxHa=	NA,
                       JmaxHd=	NA,
                       JmaxS= NA,
                       VcmaxRef=NA,
                       VcmaxHa	=NA,
                       VcmaxHd	=NA,
                       VcmaxS	=NA,
                       TpRef=NA,
                       TpHa=NA,
                       TpHd=NA,
                       TpS=NA,
                       thetacj=NA,
                       thetaip=NA,
                       RdRef=	NA,
                       RdHa=	NA,
                       RdHd=NA,
                       RdS=NA,
                       KcRef=	NA,
                       KcHa=	NA,
                       KoRef=	NA,
                       KoHa=	NA,
                       GstarRef= NA,
                       GstarHa	=NA,
                       abso=	NA,
                       aQY=NA,
                       Theta=NA,
                       model.gs=NA,
                       g0=NA,
                       g1=NA,
                       power=NA
                       ){
  if(!is.na(model.gs)&model.gs=="USO"){model.gs=0}else if(!is.na(model.gs)&model.gs=="USO_simpl"){model.gs=1}else if(!is.na(model.gs)&model.gs=="BWB"){model.gs=2}else if(!is.na(model.gs)&model.gs=="Nonlinear"){model.gs=3}else if(!is.na(model.gs)&!model.gs%in%c("USO","USO_simpl","BWB","Nonlinear")){print("Unknown model.gs")}
  param=list(R=8.314,O2=210,TRef=298.16,Patm=101,
               JmaxRef=	83.5,JmaxHa=	43540,JmaxHd=	152040,JmaxS	=495,
               VcmaxRef=	50,VcmaxHa	=65330,VcmaxHd	=149250,VcmaxS	=485,
               TpRef=1/6*50,TpHa=53100,TpHd=150650,TpS=490,
               thetacj=0.999,thetaip=0.999,
               RdRef=	1.43,RdHa=	46390,RdHd=150650,RdS=490,
               KcRef=	404.9,KcHa=	79430,KoRef=	278.4,KoHa=	36380,GstarRef=	42.75,GstarHa	=37830,
               abso=	0.85,aQY=	0.425,Theta=(0.7),g0=0.02,g1=4.1,model.gs=0,power=0.5)
  
  param_fun=list(R=R,O2=O2,TRef=TRef,Patm=Patm,JmaxRef=JmaxRef,JmaxHa=	JmaxHa,
                 JmaxHd=	JmaxHd,JmaxS	=JmaxS,VcmaxRef=VcmaxRef,VcmaxHa	= VcmaxHa,VcmaxHd	=VcmaxHd,
                 VcmaxS	=VcmaxS,
                 TpRef=TpRef,TpHa=TpHa,TpHd=TpHd,TpS=TpS,
                 thetacj=thetacj,thetaip=thetaip,
                 RdRef=RdRef,RdHa=RdHa, RdHd=RdHd,RdS=RdS,
                 KcRef= KcRef,KcHa=	KcHa,KoRef=KoRef,KoHa=	KoHa,GstarRef=	GstarRef,
                 GstarHa	=GstarHa,abso=	abso,aQY=aQY,Theta=Theta,g0=g0,
                 g1=g1,model.gs=model.gs,power=power)
  modified=which(lapply(X=param_fun,FUN = is.na)==FALSE)
  if(length(modified)>=1){
    for(i in 1: length(modified)){param[names(modified[i])]=param_fun[modified[i]]}
  }
  param_fun[names(param)]=param
  return(param_fun)
}

#' @title Photosynthesis model
#' @description Calculate the assimilation according to Farquhar equations. Contrary to f.A, this function uses inter-cellular CO2 and not ambient air CO2
#' @inheritParams f.make.param
#' @inheritParams f.A
#' @param ci Intercellular CO2 concentration in ppm 
#' @param param List of parameters, see f.make.param for details
#' 
#'
#' @return List of different variables:
#'  - A: Raw assimilation of the leaf in micromol.m-2.s-1.
#'  - Ac: Rubisco limitation assimilation of the leaf in micromol.m-2.s-1.
#'  - Aj: Electron transport rate assimilation of the leaf in micromol.m-2.s-1.
#'  - Ap: TPU rate of the leaf in micromol.m-2.s-1.
#'  - Ag: Gross assimilation in micromol.m-2.s-1.
#' @export
#'
#' @examples ci=seq(40,1500,10)
#' plot(x=ci,y=f.Aci(PFD=2000,ci=ci,Tleaf=300,param=f.make.param())$A)
f.Aci<-function(PFD,ci,Tleaf,param=f.make.param()){
  
  Kc=f.arrhenius(param[['KcRef']],param[['KcHa']],Tleaf)
  Ko=f.arrhenius(param[['KoRef']],param[['KoHa']],Tleaf)
  Gstar=f.arrhenius(param[['GstarRef']],param[['GstarHa']],Tleaf)
  
  Rd=f.modified.arrhenius(PRef=param[['RdRef']],param[['RdHa']],param[['RdHd']],param[['RdS']],Tleaf)
  Vcmax=f.modified.arrhenius(PRef=param[['VcmaxRef']],param[['VcmaxHa']],param[['VcmaxHd']],param[['VcmaxS']],Tleaf)
  Jmax=f.modified.arrhenius(PRef=param[['JmaxRef']],param[['JmaxHa']],param[['JmaxHd']],param[['JmaxS']],Tleaf)
  
  I2=PFD*param[['abso']]*(param[['aQY']])
  J=(I2+Jmax-((I2+Jmax)^2-4*(param[['Theta']])*I2*Jmax)^0.5)/(2*(param[['Theta']]))
  
  Tp=f.modified.arrhenius(PRef=param[['TpRef']],param[['TpHa']],param[['TpHd']],param[['TpS']],Tleaf)
  Wp=3*Tp
  Wc=Vcmax*(ci-Gstar)/(ci+Kc*(1+param[['O2']]/Ko))
  Wj=J/4*(ci-Gstar)/(ci+2*Gstar)
  Ai=f.smooth(A1 = Wc,A2 = Wj,theta=param[['thetacj']])*as.numeric(ci>Gstar)+f.smooth(A1 = Wc,A2 = Wj,theta=param[['thetacj']],root = 2)*as.numeric(ci<=Gstar)
  A=f.smooth(A1=Ai,A2=Wp,theta=param[['thetaip']])-Rd

  
  result=data.frame(A=A,Ac=Wc-Rd,Aj=Wj-Rd,Ap=Wp-Rd,Ag=A+Rd)
  return(result)
}


## @title Analytical solution for f.Aci with mesophyll conductance
## @param x1
## @param x2
## @param gm
## @param Rd
## @param Gstar
## @param ci
##
## @return
## @keywords internal
##
## @examples
#f.solv.cc=function(x1,x2,gm,Rd,Gstar,ci){
#  a=-gm
#  b=ci*gm-x2*gm-x1+Rd
#  c=x2*ci*gm+Gstar*x1+Rd*x2
#  sol1=(-b-sqrt(b^2-4*a*c))/(2*a)
#  sol2=(-b+sqrt(b^2-4*a*c))/(2*a)
#  return(sol1)
#}
#
#
#' @title Smoothing functions between photosynthesis limitations (for example between rubisco carboxylation and light limitation)
#' @param A1 Photosynthesis rate under limitation state 1 in micro mol m-2 s-1.
#' @param A2 Photosynthesis rate under limitation state 1 in micro mol m-2 s-1.
#' @param theta Smoothing factor between 0 and 1.
#' @return Smoothed photosynthesis value in micro mol m-2 s-1.
#' @export
#'
#' @examples A1= seq(0,20,1)
#' A2= seq(9,11,2/20)
#' Asmooth=f.smooth(A1=A1,A2=A2,theta=0.99)
#' plot(A1,type='l')
#' lines(A2)
#' lines(Asmooth,col='blue')
f.smooth=function(A1,A2,theta,root=1){
  if(root==1){sol=((A1+A2)-sqrt((A1+A2)^2-4*theta*A1*A2))/(2*theta)}else{sol=((A1+A2)+sqrt((A1+A2)^2-4*theta*A1*A2))/(2*theta)}
   return(sol)
}

#' @title Intercellular CO2 threshold between electron transport and carboxylation limitations
#' @inheritParams f.A
#' @return Intercellular CO2 such as Wc==Wj
#' @export
#'
#' @examples f.ci.treshold(PFD=2000,Tleaf=300,param=f.make.param(VcmaxRef=60,JmaxRef=85))
#' @examples f.ci.treshold(PFD=2000,Tleaf=300,param=f.make.param(VcmaxRef=70,JmaxRef=85))
f.ci.treshold<-function(PFD,Tleaf,param){
  
  Kc=f.arrhenius(param[['KcRef']],param[['KcHa']],Tleaf)
  Ko=f.arrhenius(param[['KoRef']],param[['KoHa']],Tleaf)
  Gstar=f.arrhenius(param[['GstarRef']],param[['GstarHa']],Tleaf)
  
  Vcmax=f.modified.arrhenius(PRef=param[['VcmaxRef']],param[['VcmaxHa']],param[['VcmaxHd']],param[['VcmaxS']],Tleaf)
  Jmax=f.modified.arrhenius(PRef=param[['JmaxRef']],param[['JmaxHa']],param[['JmaxHd']],param[['JmaxS']],Tleaf)
  
  I2=PFD*param[['abso']]*(param[['aQY']])
  J=(I2+Jmax-((I2+Jmax)^2-4*(param[['Theta']])*I2*Jmax)^0.5)/(2*(param[['Theta']]))
  ci_t=(J*Kc*(1+param[['O2']]/Ko)-8*Gstar*Vcmax)/(4*Vcmax-J)
  return(ci_t)
}

#' @title Plot data and model
#' @description Plot a generic graphic with observed data and predictions. Be careful to sort the data.frame beforehand.
#' @param measures Data frame obtained from CO2 or light curve with at least columns A, Ci, Qin and Tleaf.
#' @param type Type of the curve to plot (light curve: Aq or CO2 curve Aci).
#' @param list_legend Named list where the name and values will appear in the legend.
#' @inheritParams f.A
#' @param name Name of the curve to be displayed
#'
#' @return Plot a figure
#' @export
#'
#' @examples
#' param=f.make.param()
#' A=f.Aci(PFD=2000,Tleaf=300,ci=seq(40,1500,50),param=param)$A+rnorm(n = 30,mean = 0,sd = 0.5)
#' data=data.frame(Tleaf=rep(300,30),Ci=seq(40,1500,50),Qin=rep(2000,30),A=A)
#' f.plot(measures=data,param=param,list_legend=param['VcmaxRef'],name='Example 01',type='Aci')

f.plot<-function(measures=NULL,list_legend,param,name='',type='Aci',...){
  # Plot all data points
  if(type=='Aci'){x=measures$Ci
  xlab=expression(italic(C)[i]~ppm)}
  if(type%in%c('Aq','AQ')){x=measures$Qin
  xlab=expression(italic(Q)['in']~mu*mol~m^-2~s^-1)}
  if(!type%in%c('Aci','AQ','Aq')){print('type should be Aci or Aq')}
  plot(x=x,y=measures$A, main=name, xlab=xlab, ylab=expression(italic(A)~mu*mol~m^-2~s^-1),ylim=c(min(measures$A,na.rm = TRUE),1.15*max(measures$A,na.rm = TRUE)),...)
  if(!is.null(list_legend)){
    list_legend=list_legend[order(names(list_legend))]
    legend("bottomright",legend=mapply(FUN = function(x, i){paste(i,'=', round(x,2))}, list_legend, names(list_legend)),bty="n",cex=1)
  }
  legend("topleft",legend=c(expression(italic(A)[c]),expression(italic(A)[j]),expression(italic(A)[p]),expression(italic(A)),"Obs"),lty=c(2,2,2,1,0),
         pch=c(NA,NA,NA,NA,21),
         col=c("dark blue","dark red","dark green","dark grey","black"),bty="n",lwd=c(2,2,2,1,1),
         seg.len=2,cex=1,pt.cex=1)
  result=f.Aci(ci=measures$Ci,Tleaf=measures$Tleaf,PFD=measures$Qin,param=param)
  lines(x=x,y=result$A,col="dark grey",lwd=1)
  lines(x=x,y=result$Ac,lwd=2,col="dark blue",lty=2)
  lines(x=x,y=result$Aj,lwd=2,col="dark red",lty=2)
  lines(x=x,y=result$Ap,lwd=2,col="dark green",lty=2)
  box(lwd=1)
}


#' @title Compute the sum square of the difference between obervations and predictions
#' @description Function used to fit the parameters of a CO2 curve
#' @param x List of parameters to fit.
#' @param data Data frame obtained from CO2 curve with at least columns A, Ci, Qin and Tleaf.
#' @keywords internal
#' @return Sum square of the difference between predictions and observations
#' @export
#'
#' @examples
f.SumSq<-function(Fixed,data,Start){
  param=c(Fixed,Start)
  y<-data$A-f.Aci(ci=data$Ci,PFD=data$Qin,Tleaf=data$Tleaf,param=param)$A
  return(sum(y^2))
}

#' Title
#'
#' @inheritParams f.make.param
#' @param sigma Sigma value
#' @return
#' @export
#' @keywords internal
#'
#' @examples
f.MinusLogL<-function(data,sigma,R=0.75,O2=0.75,TRef=0.75,
                      Patm=0.75,JmaxRef=	0.75,
                      JmaxHa=	0.75,
                      JmaxHd=	0.75,
                      JmaxS= 0.75,
                      VcmaxRef=0.75,
                      VcmaxHa	=0.75,
                      VcmaxHd	=0.75,
                      VcmaxS	=0.75,
                      TpRef=0.75,
                      TpHa=0.75,
                      TpHd=0.75,
                      TpS=0.75,
                      thetacj=0.75,
                      thetaip=0.75,
                      RdRef=	0.75,
                      RdHa=	0.75,
                      RdHd=0.75,
                      RdS=0.75,
                      KcRef=	0.75,
                      KcHa=	0.75,
                      KcQ10=0.75,
                      KoRef=	0.75,
                      KoHa=	0.75,
                      GstarRef= 0.75,
                      GstarHa	=0.75,
                      abso=	0.75,
                      aQY=0.75,
                      Theta=0.75,
                      model.gs=NA,
                      g0=0.75,
                      g1=0.75,
                      power=0.75){
  
  param=list(R=R,O2=O2,TRef=TRef,Patm=Patm,JmaxRef=JmaxRef,JmaxHa=	JmaxHa,
             JmaxHd=	JmaxHd,JmaxS	=JmaxS,VcmaxRef=VcmaxRef,VcmaxHa	= VcmaxHa,VcmaxHd	=VcmaxHd,
             VcmaxS	=VcmaxS,
             TpRef=TpRef,TpHa=TpHa,TpHd=TpHd,TpS=TpS,
             thetacj=thetacj,thetaip=thetaip,
             RdRef=RdRef,RdHa=RdHa, RdHd=RdHd,RdS=RdS,
             KcRef= KcRef,KcHa=	KcHa,KoRef=KoRef,KoHa=	KoHa,GstarRef=	GstarRef,
             GstarHa	=GstarHa,abso=	abso,aQY=aQY,Theta=Theta,model.gs=model.gs,
             g0=g0,g1=g1,power=power)
  A_pred=f.Aci(ci=data$Ci,PFD=data$Qin,Tleaf=data$Tleaf,param=param)
  
  y<-dnorm(x=data$A,mean=A_pred$A,sd=(sigma),log=TRUE)
  return(-sum(y))
}

#' @title Fitting function for photosynthesis data (light curve or Aci curve)
#' @description Function to fit model to data. The parameters to fit have to be described in the list Start.
#' All the other parameters of the f.Aci functions have to be in param. If the parameters from Start are repeated in param, the later ones will be ignored.
#' This function uses two methods to fit the data. First by minimizing the residual sum-of-squares of the residuals and then by maximizing the likelihood function. The first method is more robust but the second one allows to calculate the confident interval of the parameters.
#' @param measures Data frame of measures obtained from gas exchange analyser with at least the columns A, Ci, Qin (the light at the leaf surface in micro mol m-2 s-1) and Tleaf (in K). If RHs, Tair, Patm, VPDleaf are also present, their mean will be added in the output, but those columns are not needed to estimate the parameters.
#' @param id.name Name of the column in the data frame measures with the name for the curve.
#' @param Start List of parameters to fit with their initial values.
#' @param param See f.make.param() for details.
#' @param modify.init TRUE or FALSE, allows to modify the Start values before fitting the data
#' @param do.plot TRUE or FALSE, plot data and fitted curves ?
#' @return Return a list with 3 components, 1 the result of the optim function which is used to estimate the parameters, 2 the output of the function bbmle, 3 the mean variable of the environment during the measurement
#' @export
#'
#' @examples ##Simulation of a CO2 curve
#' data=data.frame(Tleaf=rep(300,20),
#' Ci=seq(40,1500,75),Qin=rep(2000,20),Tair=300,RHs=70,VPDleaf=2,Patm=101,A=f.Aci(PFD=2000,Tleaf=300,ci=seq(40,1500,75),
#' param=f.make.param())$A+rnorm(n = 20,mean = 0,sd = 0.5))
#'
#' f.fitting(measures=data,id.name=NULL,Start=list(JmaxRef=90,VcmaxRef=70,RdRef=1),param=f.make.param())
f.fitting<-function(measures,id.name=NULL,Start=list(JmaxRef=90,VcmaxRef=70,RdRef=1),param=f.make.param(),modify.init=TRUE,do.plot=TRUE,type='Aci'){
  Fixed=param[!names(param)%in%names(Start)]
  if(modify.init){
    if('JmaxRef'%in%names(Start)){Start[['JmaxRef']]=f.modified.arrhenius.inv(P = 6*(max(measures$A,na.rm=TRUE)+1),Ha = param[['JmaxHa']],Hd = param[['JmaxHd']],s = param[['JmaxS']],Tleaf = mean(measures$Tleaf,na.rm=TRUE),TRef = param[['TRef']],R = param[['R']])}
    if('JmaxRef'%in%names(Start)&'VcmaxRef'%in%names(Start)){Start[['VcmaxRef']]=Start[['JmaxRef']]/2}
    grille=expand.grid(lapply(X = Start,FUN = function(x){x*c(0.2,1,2)}))
    grille.list=apply(X=grille,MARGIN = 1,FUN=as.list)
    value=9999999
    l_param=0
    for(l in 1:nrow(grille)){
      MoindresCarres=optim(par=grille.list[[l]],fn=f.SumSq,data=measures,Fixed=Fixed)
      if(!is.null(MoindresCarres)&MoindresCarres$value<value){value=MoindresCarres$value;l_param=l}
    }
    Start=grille.list[[l_param]]
  }
  
  if(is.null(id.name)){name=''}else{name=unique(measures[,id.name])}
  MoindresCarres<-Estimation2<-NULL
  try({
    MoindresCarres<-optim(par=Start,fn=f.SumSq,data=measures,Fixed=Fixed)
    print(MoindresCarres)
    print(paste('sd',sqrt(MoindresCarres$value/NROW(measures))))
    Start$sigma=sqrt(MoindresCarres$value/NROW(measures))
    for(l.name in names(MoindresCarres$par)){Start[l.name]=MoindresCarres$par[[l.name]]}
    for(l.name in names(MoindresCarres$par)){param[l.name]=MoindresCarres$par[[l.name]]}
    #if(do.plot){f.plot(measures=measures,name=name,param =param,list_legend = Start,type=type)}
  })
  
  try({
    Estimation2=mle2(minuslogl = f.MinusLogL,start = Start,fixed = Fixed,data = list(data=measures))
    print(summary(Estimation2))
    #conf=confint(Estimation2)
    #print(conf)
    for(i in names(Estimation2@coef[names(Estimation2@coef)%in%names(param)])){param[i]=Estimation2@coef[i]}
    if(do.plot){f.plot(measures=measures,name=name,param =param,list_legend = as.list(Estimation2@coef),type=type)}
  })
  Envir=NA
  try({
    Envir=c(Tair=mean(measures$Tair,na.rm=TRUE),Tleaf=mean(measures$Tleaf,na.rm=TRUE),RHs=mean(measures$RHs,na.rm=TRUE),VPDleaf=mean(measures$VPDleaf,na.rm=TRUE),Qin=mean(measures$Qin,na.rm=TRUE),Patm=mean(measures$Patm,na.rm=TRUE))
  })
   return(list(MoindresCarres,Estimation2,Envir))
}


## @title Function logit
## @description This function takes it values in 0;1- and returns values in Inf;+Inf. It is the inverse function of f.logistic
## @param x
## @return
## @export
### @encoding UTF-8
## @examples
## plot(x=seq(0,1,0.01),y=f.logit(x=seq(0,1,0.01)))
#f.logit<-function(x){ return(log(x/(1-x)))}



## @title Logistic function
## @description This function takes it values in -Inf;+Inf and returns values in 0;1. It is the inverse function of f.logit, ie f.logistic(logit(x))=x
## @details f.logistic(x)=1/(1+exp(-x)) if x<0, = exp(x)/(1+exp(x)) if x<=0
## @param x
## @return
## @export
#
## @examples
## plot(x=seq(-10,10,0.1),y=f.logistic(x=seq(-10,10,0.1)))
#f.logistic<-function(x){ return(ifelse(x>0,1/(1+exp(-x)),exp(x)/(1+exp(x))))}



###pack <- "LeafGasExchange"
#path <- find.package(pack)
#system(paste(shQuote(file.path(R.home("bin"), "R")),"CMD", "Rd2pdf", shQuote(path)))
