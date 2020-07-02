#' @title Leaf water vapour pressure deficit calculation
#' @description This function calculates the leaf water pressure deficit (VPDl or Ds) using the temperature of the leaf, the temperature of the air and its relative humidity
#' @param Tleaf Temperature of the leaf in Kelvin
#' @param Tair Temperature of the air in Kelvin
#' @param RH Humidity of the air (0 to 100)
#'
#' @return Ds in Pascal
#' @export
#' @examples f.ds(Tleaf=273.16 + 30, Tair=273.16+28, RH=70)
f.ds<-function(Tleaf,Tair,RH){
  ds=(0.6108*exp(17.27*(Tleaf-273.16)/(Tleaf-273.16+237.3))-0.6108*exp(17.27*(Tair-273.16)/(Tair-273.16+237.3))*RH/100)*1000
  return(ds)
}


#' @title Conductance model for stomatal conductance to water vapour
#' @description Semi-empirical model of the leaf conductance to water vapour
#' @param A Net assimilation in micromol.m-2.s-1, i-e, the assimilation in presence of respiration
#' @param cs CO2 at the surface of the leaf in ppm
#' @param ds Leaf surface to air vapour pressure deficit in Pa
#' @param g0 Constant of the USO model, representing the conductance when A is 0, in mol.m-2.s-1
#' @param g1 Slope parameter, between 1.14 and 3.58 KPa^0.5 (Wu et al., 2019)
#' @param power Power of the VPDl in USO model. By default is is 0.5 as in Medlin publication
#' @param model Stomatal model ("USO" or "USO_simpl")
#' @export
#' @return This function returns the optimal stomatal conductance to water vapour in mol.m-2.s-1
#' @references Medlyn, B.E., Duursma, R.A., Eamus, D., Ellsworth, D.S., Colin Prentice, I., Barton, C.V.M., Crous, K.Y., de Angelis, P., Freeman, M. and Wingate, L. (2012), Reconciling the optimal and empirical approaches to modelling stomatal conductance. Glob Change Biol, 18: 3476-3476. doi:10.1111/j.1365-2486.2012.02790.x
#'  Wu, J, Serbin, SP, Ely, KS, et al. The response of stomatal conductance to seasonal drought in tropical forests. Glob Change Biol. 2020; 26: 823– 839. https://doi.org/10.1111/gcb.14820
#'
#' @examples gs=f.gs(A=30,cs=400,ds=1500,g0=0.01,g1=2,power=0.5)
f.gs<-function(A,cs,ds,g0,g1,power=0.5,model="USO"){
  if(model=="USO"|model==0){
    gs=g0+1.6*(1+g1/(ds/1000)^power)*(A)/cs
  } else if(model=="USO_simpl"|model==1){
    gs=g0+1.6*(g1/(ds/1000)^power)*(A)/cs
  } else{print(paste("Parameter model =",model,"is not in the list of implemented models"))}
  return(gs)
}

#' @title Calculation of the minimal conductance given by a particular coupled conductance and photosynthesis model
#' @description The minimal conductance of a model depends on the parameters of the model (ie g0 and g1) but also on the minimum A value, which corresponds to the dark respiration.
#' Knowing the minimal conductance is important because the conductance can become negative and lead to unrealistic values in photosynthesis models
#' @inheritParams f.make.param
#'
#' @return Minimal conductance
#' @export
#'
#' @examples gs_min=f.gsmin(RdRef=	0.825,RdHa=	46390,RdHd=150650,RdS=490,Tleaf=300,cs=400,ds=1000,g0=0.02,g1=4.1,power=0.5,model="USO")
f.gsmin<-function(RdRef=	0.825,RdHa=	46390,RdHd=150650,RdS=490,Tleaf=300,cs=400,ds=1000,g0=0.02,g1=4.1,power=0.5,model="USO"){
  Rd=f.modified.arrhenius(PRef=RdRef,Ha = RdHa,Hd = RdHd,s = RdS,Tleaf = Tleaf)
  if(model=="USO"|model==0){
    gsmin=g0+1.6*(1-g1/(ds/1000)^power)*(Rd)/cs
  } else if(model=="USO_simpl"|model==1){
    gsmin=g0-1.6*g1*Rd/(cs*(ds/1000)^power)
  } else{print(paste("Parameter model =",model,"is not in the list of implemented models"))}
  return(gsmin)

}
#' @title Temperature dependence of Gamma star, Ko, Kc and Rd
#'
#' @param PRef Value of the parameter at the reference temperature
#' @param TRef Reference temperature
#' @param Ha Enthalpie of activation in J.mol-1
#' @param Tleaf Temperature of the leaf in Kelvin
#' @param R Ideal gas constant
#'
#' @return Value of the parameter at the temperature of the leaf
#' @export
#' @references VON CAEMMERER, S. (2013), Steady state models of photosynthesis. Plant Cell Environ, 36: 1617-1630. doi:10.1111/pce.12098
#'  Bernacchi, C.J., Singsaas, E.L., Pimentel, C., Portis Jr, A.R. and Long, S.P. (2001), Improved temperature response functions for models of Rubisco‐limited photosynthesis. Plant, Cell & Environment, 24: 253-259. doi:10.1111/j.1365-3040.2001.00668.x
#' @examples plot(x=seq(25,35,0.1),y=f.arrhenius(PRef=1,Ha=46390,Tleaf=seq(273.15+25,273.15+35,0.1),R=8.314),xlab='Temperature degree C',ylab='Rd')

f.arrhenius<-function(PRef,Ha,Tleaf,TRef=298.16,R=8.314){
  P=PRef*exp(Ha/(R*TRef)-Ha/(R*Tleaf))
  return(P)
}

#' @title Temperature dependence of Gamma star, Ko, Kc and Rd
#' @details Retrieve the value of the parameter at Tref knowing its value at Tleaf
#' @inheritParams f.arrhenius
#' @param P Value of the parameter at Tleaf
#' @return
#' @export
#'
#' @examples
f.arrhenius.inv<-function(P,Ha,Tleaf,TRef=298.16,R=8.314){
 PRef=P/exp(Ha/(R*TRef)-Ha/(R*Tleaf))
  return(PRef)
}

#' @title Temperature dependence of photosynthetic parameters
#' @details This equation is used in JULES TBM model
#' @inheritParams f.make.param
#'
#' @return Value of the photosynthetic parameter at the specified leaf temperature
#' @export
#' @references Clark, D. B., Mercado, L. M., Sitch, S., Jones, C. D., Gedney, N., Best, M. J., . Cox, P. M. (2011). The Joint UK Land Environment Simulator (JULES), model description - Part 2: Carbon fluxes and vegetation dynamics. Geoscientific Model Development, 4(3), 701-722. doi:10.5194/gmd-4-701-2011
#' @examples
f.Q10=function(Pref,Q10,Tleaf,TRef){
  P=Pref*Q10^(0.1*(Tleaf-TRef))
  return(P)
}

#' @title Temperature dependence of photosynthetic parameters
#' @details This equation is used in JULES TBM model
#'
#' @inheritParams f.make.param
#'
#' @return Value of the photosynthetic parameter at the specified leaf temperature
#' @references Clark, D. B., Mercado, L. M., Sitch, S., Jones, C. D., Gedney, N., Best, M. J., . Cox, P. M. (2011). The Joint UK Land Environment Simulator (JULES), model description - Part 2: Carbon fluxes and vegetation dynamics. Geoscientific Model Development, 4(3), 701-722. doi:10.5194/gmd-4-701-2011
#'
#' @examples
f.Q10.modified=function(Pref,Q10,Tleaf,TRef,Tlow,Tup){
  P=Pref*Q10^(0.1*(Tleaf-TRef))/((1+exp(0.3*(Tleaf-Tup)))*(1+exp(0.3*(Tlow-Tleaf))))
return(P)
  }

#' @title Temperature dependence of Jmax and Vcmax
#' @description The temperature dependence of the photosynthetic parameters Vcmax, the maximum catalytic rate of the enzyme Rubisco, and Jmax, the maximum electron transport rate is modelled by a modified Arrehenius equation. It is modified to account for decreases in each parameter at high temperatures.
#' @param PRef Value of the parameter, here Vcmax or Jmax, at the reference temperature in micromol.m-2.s-1
#' @param Ha Energy of activation in J.mol-1
#' @param Hd Energy of desactivation in J.mol-1
#' @param s Entropy term in J.mol-1.K-1
#' @param Tleaf Temperature of the leaf in Kelvin
#' @param TRef Reference temperature
#' @param R Ideal gas constant
#'
#' @return Value of the parameter Jmax or Vcmax at a given temperature
#' @references Leuning, R. (2002), Temperature dependence of two parameters in a photosynthesis model. Plant, Cell & Environment, 25: 1205-1210. doi:10.1046/j.1365-3040.2002.00898.x
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
#' @return
#' @export
#'
#' @examples
f.modified.arrhenius.inv<-function(P,Ha,Hd,s,Tleaf,TRef=298.16,R=8.314){
  PRef=P/(1+exp((s*TRef-Hd)/(R*TRef)))/exp(Ha/(R*TRef)*(1-TRef/Tleaf))*(1+exp((s*Tleaf-Hd)/(R*Tleaf)))
  return(PRef)
  }

#' @title Coupled conductance photosynthesis model
#' @description Photosynthesis model at the leaf level using the farquhar equations. The parameters can be defined by the
#' function f.make param and corresponds to the parameters inplemented in different Terrestrial Biosphere Modesl such as ORCHIDEE, JULES, CLM4.5 or FATES
#' @inheritParams f.make.param
#' @export
#' @return List of different variables:
#'  - A: Raw assimilation of the leaf in micromol.m-2.s-1
#'  - gs: Conductance of the leaf for water vapour
#'  - ci: Intracellular CO2 concentration in micromol.mol-1
#'  - cc: Mesophyll CO2 concentration in micromol.mol-1 (for the models using mesophyll conductance)
#'  - ds: Leaf surface to air vapour pressure deficit in Pa
#'
#' @examples f.A(PFD=2000,cs=400,Tleaf=273.16+29,Tair=273.16+28,RH=70,param=f.make.param())
f.A<-function(PFD,cs,Tleaf,Tair,RH,param=list()){
  #Calculation of temperature dependence of the parameters
  Kc=f.arrhenius(param[['KcRef']],param[['KcHa']],Tleaf)
  Ko=f.arrhenius(param[['KoRef']],param[['KoHa']],Tleaf)
  Gstar=f.arrhenius(param[['GstarRef']],param[['GstarHa']],Tleaf)
  Rd=f.modified.arrhenius(PRef=param[['RdRef']],param[['RdHa']],param[['RdHd']],param[['RdS']],Tleaf)
  Vcmax=f.modified.arrhenius(PRef=param[['VcmaxRef']],param[['VcmaxHa']],param[['VcmaxHd']],param[['VcmaxS']],Tleaf)
  Jmax=f.modified.arrhenius(PRef=param[['JmaxRef']],param[['JmaxHa']],param[['JmaxHd']],param[['JmaxS']],Tleaf)


  I2=PFD*param[['abso']]*(param[['aQY']])
  J=(I2+Jmax-((I2+Jmax)^2-4*(param[['Theta']])*I2*Jmax)^0.5)/(2*(param[['Theta']]))

  ds=f.ds(Tleaf,Tair,RH)
  cc=NA

  #Resolution for CLM4.5 and FATES
  if(param[['TBM']]%in%c(0,2)){

    # Analytical solution of the system of equations {E1 : A=f(ci), E2 : gs=f(A,cs) and ci=f(cs)}
    cic=f.solv(x=Vcmax,y=1,z=Kc*(1+param[['O2']]/Ko),cs=cs,Rd=Rd,Gstar=Gstar,g0=param[['g0']],g1=param[['g1']],power=param[['power']],ds=ds,model=param[["model.gs"]])
    Wc=Vcmax*cic/(cic+Kc*(1+param[['O2']]/Ko))

    cij=f.solv(x=J,y=4,z=8*Gstar,cs=cs,Rd=Rd,Gstar=Gstar,g0=param[['g0']],g1=param[['g1']],param[['power']],ds=ds,model=param[["model.gs"]])
    Wj=J*cij/(4*cij+8*Gstar)

    Tp=f.modified.arrhenius(PRef=param[['TpRef']],param[['TpHa']],param[['TpHd']],param[['TpS']],Tleaf)
    Wp=3*Tp
    cip=f.solv(x=3*Tp,y=1,z=-Gstar,cs=cs,Rd=Rd,Gstar=Gstar,g0=param[['g0']],g1=param[['g1']],param[['power']],ds=ds,model=param[["model.gs"]])

    ci=cij
    if(!is.null(which(Wc<Wj))&length(cic)==length(cij)){ci[which(Wc<Wj)]=cic[which(Wc<Wj)]}
    if(!is.null(which(Wc<Wj))&length(cic)!=length(cij)){ci[which(Wc<Wj)]=cic}
    W=pmin(Wc,Wj)
    if(!is.null(which(Wp<W))){ci[which(Wp<W)]=cip[which(Wp<W)]}
    Wc=Wc*(1-Gstar/ci)
    Wj=Wj*(1-Gstar/ci)
    Wp=Wp*(1-Gstar/ci)
    Ai=f.smooth(A1 = Wc,A2 = Wj,theta=param[['thetacj']])
    A=f.smooth(A1=Ai,A2=Wp,theta=param[['thetaip']])-Rd
    gs=f.gs(A=A,cs=cs,ds=ds,g0=param[['g0']],g1=param[['g1']],power=param[['power']],model =param[['model.gs']])
    output=list(A=A,gs=gs,ci=ci,ds=ds)
    return(output)
  }

  #Resolution for JULES
  if(param[['TBM']]==3){
    Kc=f.Q10(Pref = param[['KcRef']],Q10 = param[['KcQ10']],Tleaf,TRef=param[['TRef']])
    Ko=f.Q10(Pref = param[['KoRef']],Q10 = param[['KoQ10']],Tleaf,TRef=param[['TRef']])
    PFD2=PFD
    PFD2[PFD2<11]=exp(-10)
    Rd=f.Q10.modified(Pref = param[['RdRef']],Q10 = param[['VcmaxQ10']],Tlow = param[['Tlow']],Tup = param[['Tup']],Tleaf=Tleaf,TRef=param[['TRef']])
    #It is possible to include a light inhibition of the dark respiration here
    #Rd=(0.5-0.05*log(PFD2))*f.Q10.modified(Pref = param[['RdRef']],Q10 = param[['VcmaxQ10']],Tlow = param[['Tlow']],Tup = param[['Tup']],Tleaf=Tleaf,TRef=param[['TRef']])

    Vcmax=f.Q10.modified(Pref = param[['VcmaxRef']],Q10 = param[['VcmaxQ10']],Tlow = param[['Tlow']],Tup = param[['Tup']],Tleaf=Tleaf,TRef=param[['TRef']])
    Tau=f.arrhenius(param[['TauRef']],param[['TauQ10']],Tleaf,TRef=param[['TRef']])
    Gstar=param[['O2']]/(2*Tau)
    cic=f.solv(x=Vcmax,y=1,z=Kc*(1+param[['O2']]/Ko),cs=cs,Rd=Rd,Gstar=Gstar,g0=param[['g0']],g1=param[['g1']],power=param[['power']],ds=ds,model=param[["model.gs"]])
    Wc=Vcmax*cic/(cic+Kc*(1+param[['O2']]/Ko))
    J=param[['abso']]*param[['aQY']]*PFD

    cij=f.solv(x=J,y=1,z=2*Gstar,cs=cs,Rd=Rd,Gstar=Gstar,g0=param[['g0']],g1=param[['g1']],param[['power']],ds=ds,model=param[["model.gs"]])
    Wj=J*(cij-Gstar)/(cij+2*Gstar)

    ci=cij
    Tp=f.Q10.modified(Pref = param[['TpRef']],Q10 = param[['VcmaxQ10']],Tlow = param[['Tlow']],Tup = param[['Tup']],Tleaf=Tleaf,TRef=param[['TRef']])

    Wp=3*Tp
    cip=f.solv(x=Tp,y=1,z=-Gstar,cs=cs,Rd=Rd,Gstar=Gstar,g0=param[['g0']],g1=param[['g1']],param[['power']],ds=ds,model=param[["model.gs"]])

    if(!is.null(which(Wc<Wj))&length(cic)==length(cij)){ci[which(Wc<Wj)]=cic[which(Wc<Wj)]}
    if(!is.null(which(Wc<Wj))&length(cic)!=length(cij)){ci[which(Wc<Wj)]=cic}
    W=pmin(Wc,Wj)
    if(!is.null(which(Wp<W))){ci[which(Wp<W)]=cip[which(Wp<W)]}
    Ai=f.smooth(A1 = Wc,A2 = Wj,theta=param[['thetacj']])
    A=f.smooth(A1=Ai,A2=Wp,theta=param[['thetaip']])-Rd
    gs=f.gs(A=A,cs=cs,ds=ds,g0=param[['g0']],g1=param[['g1']],power=param[['power']],model =param[['model.gs']])
    output=list(A=A,gs=gs,ci=ci,ds=ds)
    return(output)

  }

    #Resolution for ORCHIDEE
    if(param[['TBM']]==1){
    ci=NA
    gm=f.modified.arrhenius(PRef=param[['gmRef']],param[['gmHa']],param[['gmHd']],param[['gmS']],Tleaf)
    ccc=f.solv.Acc(g0=param[['g0']],g1=param[['g1']],power=param[['power']],ds=ds,model=param[['model.gs']],gm=gm,Rd=Rd,Gstar=Gstar,x1=Vcmax,x2=Kc*(1+param[['O2']]/Ko),cs=cs)
    Wc=Vcmax*ccc/(ccc+Kc*(1+param[['O2']]/Ko))
    ccj=f.solv.Acc(g0=param[['g0']],g1=param[['g1']],power=param[['power']],ds=ds,model=param[['model.gs']],gm=gm,Rd=Rd,Gstar=Gstar,x1=J/4,x2=2*Gstar,cs=cs)
    Wj=J*ccj/(4*ccj+8*Gstar)
    cc=ccj
    if(!is.null(which(Wc<Wj))&length(ccc)==length(ccj)){cc[which(Wc<Wj)]=ccc[which(Wc<Wj)]}
    if(!is.null(which(Wc<Wj))&length(ccc)!=length(ccj)){cc[which(Wc<Wj)]=ccc}
    Wc=Wc*(1-Gstar/cc)
    Wj=Wj*(1-Gstar/cc)
    A=pmin(Wc,Wj)-Rd
    gs=f.gs(A=A,cs=cs,ds=ds,g0=param[['g0']],g1=param[['g1']],power=param[['power']],model =param[['model.gs']])
    output=list(A=A,gs=gs,ci=ci,cc=cc,ds=ds)
    return(output)
    }
}



#' @title Coupled conductance photosynthesis model with energy balance model
#' @details This function allo to calculate the photosynthesis from environmental variables PFD, RH, wind, cs, Tair.
#' (There is no boundary layer). The energy balance model comes from the package Tealeaves (see reference). The energy balance calculation involves the stomatal conductance and the cuticular conductance.
#' Here the cuticular conductance is considered to be equal to g0 even if it is wrong for the USO models. Most of the times, no precaution is taken on gs_min when fitting the conductance models so gs_min is often negative. This choice was made to prevent unrealistic energy budgets.
#' @inheritParams
#' @param param List of parameters given by f.make.param()
#' @param precision Precision of the leaf temperature prediction. The resolution of the energy balance coupled with the photosynthesis and stomatal conductance is numerical. The smaller the precision, the longer will be the resolution.
#' @param max_it Maximum number of iterations to find the solution
#' @return
#' @export
#' @references tealeaves: an R package for modelling leaf temperature using energy budgets. Christopher. D. Muir. bioRxiv 529487; doi: https://doi.org/10.1101/529487

#' @examples f.AT(PFD=1500,cs=400,Tair=299,wind=2,RH=70,param=f.make.param())
f.AT<-function(PFD,cs,Tair,RH,wind,precision=0.1,max_it=10,param=list()){
  Tleaf=Tair+1
  n=1
  delta=precision+1
  while(delta>0.1&n<10){
    Leaf_physio=f.A(PFD=PFD,Tleaf=Tleaf,Tair=Tair,cs = cs,RH = RH,param=param)
    ds=f.ds(Tleaf,Tair,RH)
    gs=Leaf_physio$gs-param[['g0']]
    if(gs<param[['g0']]){gs=param[['g0']]}
    g_sw <- set_units(gs, "mol/m^2/s")  ##Stomatal conductance
    g_uw<- set_units(param[['g0']], "mol/m^2/s")   ##Cuticular conductance
    g_sw=convert_conductance(g_sw,
                             Temp = set_units(Tair, "K"),
                             P = set_units(101.3246, "kPa"))$`umol/m^2/s/Pa`
    g_uw=convert_conductance(g_uw,
                             Temp = set_units(Tair, "K"),
                             P = set_units(101.3246, "kPa"))$`umol/m^2/s/Pa`

    leaf_par <- make_leafpar(replace = list(leafsize=set_units(c(0.04), "m"),g_sw=g_sw,g_uw=g_uw,abs_s=set_units((300*param[["abso"]]+3300*0.5242)/4000)))

    enviro_par <- make_enviropar(replace=list(S_sw=set_units(PFD/4.6*2,"W/m^2"),RH=set_units(RH/100),T_air = set_units(Tair, "K"),wind=set_units(wind,"m/s")))
    constants <- make_constants()
    Tleaf_mod <- as.numeric(tleaf(leaf_par, enviro_par, constants,quiet = TRUE)$T_leaf)
    delta=(Tleaf-Tleaf_mod)
    Tleaf=Tleaf_mod
    n=n+1
    if(n>10&delta>0.1){
      Tleaf=NA
      print(paste("Temperature convergence limit"))}
    Leaf_physio$Tleaf=Tleaf
  }
  return(Leaf_physio)
}

#' @title Analytical solution of the coupled photosynthesis and USO model
#' @param x
#' @param y
#' @param z
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
f.solv<-function(x,y,z,cs,Rd,Gstar,g0,g1,power,ds,model){
  if(model=="USO"|model==0){m=1.6*(1+g1/(ds/1000)^power)}else if(model=="USO_simpl"|model==1){m=1.6*(g1/(ds/1000)^power)}else{print(paste("model:",model,'is not in the list of implemented stomatal models'))}

  a=y*g0+m/cs*(x-Rd*y)
  b=z*g0+m/cs*(-Gstar*x-Rd*z)-cs*g0*y+(x-Rd*y)*(1.6-m)
  c=-z*cs*g0+(1.6-m)*(-Gstar*x-Rd*z)
  ci2=(-b-(b^2-4*a*c)^0.5)/(2*a)
  ci1=(-b+(b^2-4*a*c)^0.5)/(2*a)

  return(ci1)
}


#' @title Analytical solution of the coupled photosynthesis with mesophyll conductance and USO model
#' @param x1
#' @param x2
#' @param g0
#' @param g1
#' @param power
#' @param ds
#' @param model
#' @param gm
#' @param Rd
#' @param Gstar
#' @param cs
#'
#' @return
#' @export
#' @keywords internal
#' @examples
f.solv.Acc=function(x1,x2,g0,g1,power,ds,model,gm,Rd,Gstar,cs){
  if(model=="USO"|model==0){m=1.6*(1+g1/(ds/1000)^power)}else if(model=="USO_simpl"|model==1){m=1.6*(g1/(ds/1000)^power)}else{print(paste("model:",model,'is not in the list of implemented stomatal models'))}
  s=x1-Rd
  t=-Gstar*x1-Rd*x2
  u=g0+1.6*gm-cs*m*gm
  a1=gm*(g0+m*s)
  a2=m*s^2+u*s+t*m*gm+x2*m*gm*s+2*g0*gm*x2-g0*gm*cs
  a3=2*m*t*s+u*t+x2*u*s+t*m*gm*x2+g0*gm*x2^2-2*cs*g0*gm*x2
  a4=x2*u*t-cs*g0*gm*x2^2+m*t^2
  p=a2/a1
  q=a3/a1
  r=a4/a1
  U=(2*p^3-9*p*q+27*r)/54
  Q=(p^2-3*q)/9
  phi=acos(U/sqrt(Q^3))
 #sol1=-2*sqrt(Q)*cos(phi/3)-p/3
  sol2=-2*sqrt(Q)*cos((phi+2*pi)/3)-p/3
  #sol3=-2*sqrt(Q)*cos((phi+4*pi)/3)-p/3
  return(sol2)
}

#' @title Photosynthesis and stomata model parameters
#' @description Function to create a list of parameters to be used in most of the functions of this package.
#' Depending on the function, all the parameters are not used. For example go and g1 are not used in f.Aci.
#' The parameters from different TBM are implemented and can be chosen by selecting a TBM
#' @details The call of this function is made using f.make.param(). If a parameter is modified for example writting f.make.param(VcmaxRef=10), this function will return all the default parameters from FATES TBM with VcmaxRef = 10 instead of its default value
#' @param TBM Type of model (FATES, ORCHIDEE, CLM4.5 or JULES). Default is FATES
#' @param R Ideal gas constant
#' @param O2 O2 concentration in ppm
#' @param TRef Reference temperature for Kc, Ko, Rd,GammaStar Vcmax, Jmax
#' @param Patm Atmospheric pressure in Pa
#' @param JmaxRef Maximum electron transport rate in micromol.m-2.s-1
#' @param JmaxHa Energy of activation for Jmax in J.mol-1
#' @param JmaxHd Energy of desactivation for Jmax in J.mol-1
#' @param JmaxS Entropy term for Jmax in J.mol-1.K-1
#' @param VcmaxRef Maximum rate of Rubisco for carboxylation micromol.m-2.s-1
#' @param VcmaxHa Energy of activation for Vcmax in J.mol-1
#' @param VcmaxHd Energy of desactivation for Vcmax in J.mol-1
#' @param VcmaxS Entropy term for Vcmax in J.mol-1.K-1
#' @param RdRef Respiration value at the reference temperature
#' @param RdHa Energie of activation for Rd in J.mol-1
#' @param KcRef Michaelis-Menten constant of Rubisco for CO2 at the reference temperature in micromol.mol-1
#' @param KcHa Energy of activation for Kc in J.mol-1
#' @param KoRef ichaelis-Menten constant of Rubisco for CO2 at the reference temperature in milimol.mol-1
#' @param KoHa Energy of activation for Ko in J.mol-1
#' @param GstarRef CO2 compensation point in absence of respiration in micromol.mol-1
#' @param GstarHa Enthalpie of activation for Gstar in J.mol-1
#' @param abso Absorptance of the leaf in the photosynthetic active radiation wavelenghts
#' @param aQY Apparent quantum yield
#' @param Theta Theta is the empirical curvacture factor for the response of J to PFD. It takes its values between 0 and 1.
#' @param g0 Constant of the USO model, representing the conductance when A is 0, in mol.m-2.s-1
#' @param g1 Slope parameter, between 1.14 and 3.58 KPa^0.5 (Wu et al., 2019)
#' @param power Power of VPDl in USO model. By default power=0.5 as in Medlyn article
#' @param model.gs Type of conductance model (USO, USO_simpl)
#' @param gmRef Mesophyll conductance at Tref (25 deg C) mol m-2 s-1
#' @param gmS Entropy term for gm J K-1 mol-1
#' @param gmHa Energy of activation for gm in J.mol-1
#' @param gmHd Energy of deactivation for gm in J.mol-1
#' @return List of parameters that can be used in f.A
#' @references Bernacchi, C.J., Singsaas, E.L., Pimentel, C., Portis Jr, A.R. and Long, S.P. (2001), Improved temperature response functions for models of Rubisco‐limited photosynthesis. Plant, Cell & Environment, 24: 253-259. doi:10.1111/j.1365-3040.2001.00668.x
#' CLM4.5: http://www.cesm.ucar.edu/models/cesm2/land/CLM45_Tech_Note.pdf
#' ORCHIDEE: https://forge.ipsl.jussieu.fr/orchidee/wiki/Documentation/OrchideeParameters AND https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2003GB002199
#' JULES: https://www.geosci-model-dev.net/4/701/2011/gmd-4-701-2011.pdf
#' FATES: https://fates-docs.readthedocs.io/en/latest/fates_tech_note.html#
#' @export
#'
#' @examples param1=f.make.param(TBM='FATES',JmaxRef=100,VcmaxRef=60,RdRef=1,TpRef=10)
#'  param2=f.make.param(TBM='CLM4.5',JmaxRef=100,VcmaxRef=60,RdRef=1,TpRef=10)
#'  f.A(PFD=1500,cs=400,Tleaf=300,Tair=299,RH=70,param=param1)
#'  f.A(PFD=1500,cs=400,Tleaf=300,Tair=299,RH=70,param=param2)
f.make.param<-function(TBM='FATES',R=NA,O2=NA,TRef=NA,
                       Patm=NA,JmaxRef=	NA,
                       JmaxHa=	NA,
                       JmaxHd=	NA,
                       JmaxS= NA,
                       VcmaxRef=NA,
                       VcmaxHa	=NA,
                       VcmaxHd	=NA,
                       VcmaxS	=NA,
                       VcmaxQ10=NA,
                       Tlow=NA,
                       Tup=NA,
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
                       KcQ10=NA,
                       KoRef=	NA,
                       KoHa=	NA,
                       KoQ10=NA,
                       GstarRef= NA,
                       GstarHa	=NA,
                       TauRef=NA,
                       TauQ10=NA,
                       abso=	NA,
                       aQY=NA,
                       Theta=NA,
                       model.gs=NA,
                       g0=NA,
                       g1=NA,
                       power=NA,
                       gmRef=NA,
                       gmS=NA,
                       gmHa=NA,
                       gmHd=NA){
  if(!TBM%in%c('FATES','ORCHIDEE','CLM4.5','JULES')){print(paste('TBM',TBM,'is not in the list FATES, ORCHIDEE, CLM4.5, JULES'))}
  if(!is.na(model.gs)&model.gs=="USO"){model.gs=0}else if(!is.na(model.gs)&model.gs=="USO_simpl"){model.gs=1}else if(!is.na(model.gs)&!model.gs%in%c("USO","USO_simpl")){print("Unknown model.gs")}
   if(TBM=='FATES'){
    param=list(TBM=0,R=8.314,O2=210,TRef=298.16,Patm=101,
    JmaxRef=	83.5,JmaxHa=	43540,JmaxHd=	152040,JmaxS	=495,
    VcmaxRef=	50,VcmaxHa	=65330,VcmaxHd	=149250,VcmaxS	=485,
    TpRef=1/6*50,TpHa=53100,TpHd=150650,TpS=490,
    thetacj=0.999,thetaip=0.999,
    RdRef=	1.43,RdHa=	46390,RdHd=150650,RdS=490,
    KcRef=	404.9,KcHa=	79430,KoRef=	278.4,KoHa=	36380,GstarRef=	42.75,GstarHa	=37830,
    abso=	0.85,aQY=	0.425,Theta=(0.85),g0=0.02,g1=4.1,model.gs=0,power=0.5)
  }
  if(TBM=='ORCHIDEE'){
    param=list(TBM=1,R=8.314,O2=210,TRef=298.16,Patm=101,
               JmaxRef=	111.15,JmaxHa=	49884,JmaxHd=	20000,JmaxS	=640.57,
               VcmaxRef=	65,VcmaxHa	=71513,VcmaxHd	=200000,VcmaxS	=641.1,
               RdRef=0.01*65,RdHa=	46390,RdHd=150650,RdS=490,
               KcRef=	404.9,KcHa=	79430,KoRef=	278.4,KoHa=	36380,GstarRef=	42.75,GstarHa	=37830,
               abso=	0.84,aQY=	0.372,Theta=(0.7),g0=0.02,g1=4.1,model.gs=0,power=0.5,
               gmRef=0.4,gmS=1400,gmHa=49600,gmHd=437400)
  }
  if(TBM=='CLM4.5'){
    param=list(TBM=2,R=8.314,O2=210,TRef=298.16,Patm=101,
               JmaxRef=	94.05,JmaxHa=	50000,JmaxHd=	152040,JmaxS	=495,
               VcmaxRef=	55,VcmaxHa	=72000,VcmaxHd	=200000,VcmaxS	=641.1,
               TpRef=1/6*55,TpHa=72000,TpHd=200000,TpS=641.1,
               thetacj=0.98,thetaip=0.95,
               RdRef=	0.825,RdHa=	46390,RdHd=150650,RdS=490,
               KcRef=	404.9,KcHa=	79430,KoRef=	278.4,KoHa=	36380,GstarRef=	42.75,GstarHa	=37830,
               abso=	0.85,aQY=	0.425,Theta=(0.85),g0=0.02,g1=4.1,model.gs=0,power=0.5)
  }
  if(TBM=='JULES'){
    param=list(TBM=3,R=8.314,O2=210,TRef=298.16,Patm=101,
               VcmaxRef=36.8,VcmaxQ10=2,Tlow=0+273.15,Tup=36+273.15,
               TpRef=1/2*36.8,
               thetacj=0.83,thetaip=0.93,
               RdRef=	0.552,
               KcRef=	303,KcQ10=2.1,KoQ10=1.2,KoRef=	303,
               TauRef=2600/1010,TauQ10=0.57,
               abso=	0.85,aQY=	0.08,g0=0.02,g1=4.1,model.gs=0,power=0.5)
  }

   param_fun=list(TBM=TBM,R=R,O2=O2,TRef=TRef,Patm=Patm,JmaxRef=JmaxRef,JmaxHa=	JmaxHa,
             JmaxHd=	JmaxHd,JmaxS	=JmaxS,VcmaxRef=VcmaxRef,VcmaxHa	= VcmaxHa,VcmaxHd	=VcmaxHd,
             VcmaxS	=VcmaxS,VcmaxQ10=VcmaxQ10,Tlow=Tlow,Tup=Tup,
             TpRef=TpRef,TpHa=TpHa,TpHd=TpHd,TpS=TpS,
             thetacj=thetacj,thetaip=thetaip,
             RdRef=RdRef,RdHa=RdHa, RdHd=RdHd,RdS=RdS,
             KcRef= KcRef,KcHa=	KcHa,KcQ10=KcQ10,KoRef=KoRef,KoHa=	KoHa,KoQ10=KoQ10,GstarRef=	GstarRef,TauRef=TauRef,TauQ10=TauQ10,
             GstarHa	=GstarHa,abso=	abso,aQY=aQY,Theta=Theta,g0=g0,
             g1=g1,model.gs=model.gs,power=power,gmRef=gmRef,gmS=gmS,gmHa=gmHa,gmHd=gmHd)
  modified=which(lapply(X=param_fun,FUN = is.na)==FALSE)
  if(length(modified)>1){
    for(i in 2: length(modified)){param[names(modified[i])]=param_fun[modified[i]]}
  }
  param_fun[names(param)]=param
  return(param_fun)
}

#' @title Photosynthesis model
#' @description Calculate the assimilation according to Farquhar equations. Contrary to f.A, this function uses intracellular CO2 and not ambiant air CO2
#' @inheritParams f.make.param
#' @param param List of parameters, see f.make.param for details
#'
#' @return Assimilation in micromol.m-2.s-1
#' @export
#'
#' @examples ci=seq(40,1500,10)
#' plot(x=ci,y=f.Aci(PFD=2000,ci=ci,Tleaf=300,param=f.make.param())$A)
f.Aci=function(PFD,ci,Tleaf,param=list()){

  Kc=f.arrhenius(param[['KcRef']],param[['KcHa']],Tleaf)
  Ko=f.arrhenius(param[['KoRef']],param[['KoHa']],Tleaf)
  Gstar=f.arrhenius(param[['GstarRef']],param[['GstarHa']],Tleaf)

  Rd=f.modified.arrhenius(PRef=param[['RdRef']],param[['RdHa']],param[['RdHd']],param[['RdS']],Tleaf)
  Vcmax=f.modified.arrhenius(PRef=param[['VcmaxRef']],param[['VcmaxHa']],param[['VcmaxHd']],param[['VcmaxS']],Tleaf)
  Jmax=f.modified.arrhenius(PRef=param[['JmaxRef']],param[['JmaxHa']],param[['JmaxHd']],param[['JmaxS']],Tleaf)

  I2=PFD*param[['abso']]*(param[['aQY']])
  J=(I2+Jmax-((I2+Jmax)^2-4*(param[['Theta']])*I2*Jmax)^0.5)/(2*(param[['Theta']]))


  if(param[['TBM']]==1){
    gm=f.modified.arrhenius(PRef=param[['gmRef']],param[['gmHa']],param[['gmHd']],param[['gmS']],Tleaf)
    ccc=f.solv.cc(gm=gm,Rd=Rd,Gstar=Gstar,x1=Vcmax,x2=Kc*(1+param[['O2']]/Ko),ci=ci)
    ccj=f.solv.cc(gm=gm,Rd=Rd,Gstar=Gstar,x1=J/4,x2=2*Gstar,ci=ci)

    Wc=Vcmax*(ccc-Gstar)/(ccc+Kc*(1+param[['O2']]/Ko))
    Wj=J/4*(ccj-Gstar)/(ccj+2*Gstar)
    Wp=NA
    if(length(ccc)>length(ccj)){
      cc=ccc
      cc[Wj<Wc&!is.na(Wj)&!is.na(Wc)]=ccj
    }else{
      cc=ccj
      cc[Wc<Wj&!is.na(Wj)&!is.na(Wc)]=ccc[Wc<Wj&!is.na(Wj)&!is.na(Wc)]
      }
    cc[Wj<Wc&!is.na(Wj)&!is.na(Wc)]=ccj[Wj<Wc&!is.na(Wj)&!is.na(Wc)]
    A=pmin(Wc,Wj)-Rd
  }

  if(param[['TBM']]%in%c(0,2)){
    Tp=f.modified.arrhenius(PRef=param[['TpRef']],param[['TpHa']],param[['TpHd']],param[['TpS']],Tleaf)
    Wp=3*Tp
    Wc=Vcmax*(ci-Gstar)/(ci+Kc*(1+param[['O2']]/Ko))
    Wj=J/4*(ci-Gstar)/(ci+2*Gstar)
    Ai=f.smooth(A1 = Wc,A2 = Wj,theta=param[['thetacj']])
    A=f.smooth(A1=Ai,A2=Wp,theta=param[['thetaip']])-Rd
  }

  if(param[['TBM']]==3){
    Kc=f.Q10(Pref = param[['KcRef']],Q10 = param[['KcQ10']],Tleaf,TRef=param[['TRef']])
    Ko=f.Q10(Pref = param[['KoRef']],Q10 = param[['KoQ10']],Tleaf,TRef=param[['TRef']])
    PFD2=PFD
    PFD2[PFD2<11]=exp(-10)
    Rd=f.Q10.modified(Pref = param[['RdRef']],Q10 = param[['VcmaxQ10']],Tlow = param[['Tlow']],Tup = param[['Tup']],Tleaf=Tleaf,TRef=param[['TRef']])
    #Rd=(0.5-0.05*log(PFD2))*f.Q10.modified(Pref = param[['RdRef']],Q10 = param[['VcmaxQ10']],Tlow = param[['Tlow']],Tup = param[['Tup']],Tleaf=Tleaf,TRef=param[['TRef']])

    Vcmax=f.Q10.modified(Pref = param[['VcmaxRef']],Q10 = param[['VcmaxQ10']],Tlow = param[['Tlow']],Tup = param[['Tup']],Tleaf=Tleaf,TRef=param[['TRef']])
    Tp=f.Q10.modified(Pref = param[['TpRef']],Q10 = param[['VcmaxQ10']],Tlow = param[['Tlow']],Tup = param[['Tup']],Tleaf=Tleaf,TRef=param[['TRef']])
    Tau=f.arrhenius(param[['TauRef']],param[['TauQ10']],Tleaf,TRef=param[['TRef']])
    Gstar=param[['O2']]/(2*Tau)
    Wp=3*Tp
    Wc=Vcmax*(ci-Gstar)/(ci+Kc*(1+param[['O2']]/Ko))
    Wj=param[['abso']]*param[['aQY']]*PFD*(ci-Gstar)/(ci+2*Gstar)
    Ai=f.smooth(A1 = Wc,A2 = Wj,theta=param[['thetacj']])
    A=f.smooth(A1=Ai,A2=Wp,theta=param[['thetaip']])-Rd
  }
  result=data.frame(A=A,Ac=Wc-Rd,Aj=Wj-Rd,Ap=Wp-Rd)
  return(result)
}


#' @title Analytical solution for f.Aci with mesophyll conductance
#' @param x1
#' @param x2
#' @param gm
#' @param Rd
#' @param Gstar
#' @param ci
#'
#' @return
#' @keywords internal
#'
#' @examples
f.solv.cc=function(x1,x2,gm,Rd,Gstar,ci){
  a=-gm
  b=ci*gm-x2*gm-x1+Rd
  c=x2*ci*gm+Gstar*x1+Rd*x2
  sol1=(-b-sqrt(b^2-4*a*c))/(2*a)
  sol2=(-b+sqrt(b^2-4*a*c))/(2*a)
  return(sol1)
}


#' @title Smoothing functions between photosynthesis limitations (for example between rubisco carboxylation and light limitation)
#' @param A1
#' @param A2
#' @param theta Smoothing factor
#' @return Smoothed value
#' @export
#'
#' @examples A1= seq(0,20,1)
#' A2= seq(9,11,2/20)
#' Asmooth=f.smooth(A1=A1,A2=A2,theta=0.99)
#' plot(A1,type='l')
#' lines(A2)
#' lines(Asmooth,col='blue')
f.smooth=function(A1,A2,theta){
  return(((A1+A2)-sqrt((A1+A2)^2-4*theta*A1*A2))/(2*theta))
}

#' @title Intracellular CO2 threshold between electron transport and carboxylation limitations
#' @inheritParams f.A
#' @return Intracellular CO2 such as Wc==Wj
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
#' @param measures Data frame obtained from CO2 or light curve with at least columns Photo, Ci, PARi and Tleaf
#' Data frame obtained from CO2 or light curve with at least columns Photo, Ci, PARi and Tleaf
#' @param type Type of the curve to plot (light curve: Aq or CO2 curve Aci)
#' @param list_legend Named list where the name and values will appear in the legend
#' @inheritParams f.A
#' @param name Name of the curve to be displayed
#'
#' @return Plot a figure
#' @export
#'
#' @examples
#' param=f.make.param()
#' Photo=f.Aci(PFD=2000,Tleaf=300,ci=seq(40,1500,50),param=param)$A+rnorm(n = 30,mean = 0,sd = 0.5)
#' data=data.frame(Tleaf=rep(300,30),Ci=seq(40,1500,50),PARi=rep(2000,30),Photo=Photo)
#' f.plot(measures=data,param=param,list_legend=param['VcmaxRef'],name='Example 01',type='Aci')

f.plot<-function(measures=NULL,list_legend,param,name='',type='Aci'){
  # Plot all data points
  if(type=='Aci'){x=measures$Ci
  xlab="Ci in ppm"}
  if(type%in%c('Aq','AQ')){x=measures$PARi
  xlab="PARi"}
  if(!type%in%c('Aci','AQ','Aq')){print('type should be Aci or Aq')}
  plot(x=x,y=measures$Photo, main=name, xlab=xlab, ylab="Photo in micromol.m-2.s-1",ylim=c(min(measures$Photo,na.rm = TRUE),1.15*max(measures$Photo,na.rm = TRUE)))
  if(!is.null(list_legend)){
    list_legend=list_legend[order(names(list_legend))]
    legend("bottomright",legend=mapply(FUN = function(x, i){paste(i,'=', round(x,2))}, list_legend, names(list_legend)),bty="n",cex=1)
  }
  legend("topleft",legend=c("Rubisco","RuBP","TPU","Photo","Obs"),lty=c(2,2,2,1,0),
         pch=c(NA,NA,NA,NA,21),
         col=c("dark blue","dark red","dark green","dark grey","black"),bty="n",lwd=c(2,2,2,1,1),
         seg.len=2,cex=1,pt.cex=1)
  result=f.Aci(ci=measures$Ci,Tleaf=measures$Tleaf,PFD=measures$PARi,param=param)
  lines(x=x,y=result$A,col="dark grey",lwd=1)
  lines(x=x,y=result$Ac,lwd=2,col="dark blue",lty=2)
  lines(x=x,y=result$Aj,lwd=2,col="dark red",lty=2)
  lines(x=x,y=result$Ap,lwd=2,col="dark green",lty=2)
   box(lwd=1)
}


#' @title Compute the sum square of the difference between obervations and predictions
#' @description Function used to fit the parameters of a CO2 curve
#' @param x List of parameters to fit
#' @param data Data frame obtained from CO2 curve with at least columns Photo, Ci, PARi and Tleaf
#' @keywords internal
#' @return Sum square of the difference between predictions and observations
#' @export
#'
#' @examples
f.SumSq<-function(Fixed,data,Start){
  param=c(Fixed,Start)
  y<-data$Photo-f.Aci(ci=data$Ci,PFD=data$PARi,Tleaf=data$Tleaf,param=param)$A
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
f.MinusLogL<-function(data,sigma,TBM=0,R=0.75,O2=0.75,TRef=0.75,
                      Patm=0.75,JmaxRef=	0.75,
                      JmaxHa=	0.75,
                      JmaxHd=	0.75,
                      JmaxS= 0.75,
                      VcmaxRef=0.75,
                      VcmaxHa	=0.75,
                      VcmaxHd	=0.75,
                      VcmaxS	=0.75,
                      VcmaxQ10=0.75,
                      Tlow=0.75,
                      Tup=0.75,
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
                      KoQ10=0.75,
                      GstarRef= 0.75,
                      GstarHa	=0.75,
                      TauRef=0.75,
                      TauQ10=0.75,
                      abso=	0.75,
                      aQY=0.75,
                      Theta=0.75,
                      model.gs=NA,
                      g0=0.75,
                      g1=0.75,
                      power=0.75,
                      gmRef=0.75,
                      gmS=0.75,
                      gmHa=0.75,
                      gmHd=0.75){

  param=list(TBM=TBM,R=R,O2=O2,TRef=TRef,Patm=Patm,JmaxRef=JmaxRef,JmaxHa=	JmaxHa,
             JmaxHd=	JmaxHd,JmaxS	=JmaxS,VcmaxRef=VcmaxRef,VcmaxHa	= VcmaxHa,VcmaxHd	=VcmaxHd,
             VcmaxS	=VcmaxS,VcmaxQ10=VcmaxQ10,Tlow=Tlow,Tup=Tup,
             TpRef=TpRef,TpHa=TpHa,TpHd=TpHd,TpS=TpS,
             thetacj=thetacj,thetaip=thetaip,
             RdRef=RdRef,RdHa=RdHa, RdHd=RdHd,RdS=RdS,
             KcRef= KcRef,KcHa=	KcHa,KcQ10=KcQ10,KoRef=KoRef,KoHa=	KoHa,KoQ10=KoQ10,GstarRef=	GstarRef,TauRef=TauRef,TauQ10=TauQ10,
             GstarHa	=GstarHa,abso=	abso,aQY=aQY,Theta=Theta,model.gs=model.gs,
             g0=g0,g1=g1,power=power,gmRef=gmRef,gmS=gmS,gmHa=gmHa,gmHd=gmHd)
  A_pred=f.Aci(ci=data$Ci,PFD=data$PARi,Tleaf=data$Tleaf,param=param)$A
  y<-dnorm(x=data$Photo,mean=A_pred,sd=sigma,log=TRUE)
  return(-sum(y))
}

#' @title Fitting function for photosynthesis datadata (light curve or Aci curve)
#' @description Function to fit model to data. The parameters to fit have to be described in the list Start.
#' All the other parameters of the f.Aci functions have to be in param. If the parameters from Start are repeated in param, the later one will be ignored.
#' This function uses two methods to fit the data. First by minimizing the residual sum-of-squares of the residuals and then by maximizing the likelihood function. The first method is more robust but the second one allows to calculate the confident interval of the parameters.
#' @param measures Data frame of measures obtained from gas exchange analyser with at least the columns Photo, Ci, PARi and Tleaf (in K)
#' @param id.name Name of the colums in measures with the identifier for the curve.
#' @param Start List of parameters to fit with their initial values.
#' @param param See f.make.param() for details.
#' @param modify.init TRUE or FALSE, allows to modify the Start values before fitting the data
#' @param do.plot TRUE or FALSE, plot data and fitted curves ?
#' @return
#' @export
#'
#' @examples ##Simulation of a CO2 curve
#' data=data.frame(Tleaf=rep(300,20),
#' Ci=seq(40,1500,75),PARi=rep(2000,20),Photo=f.Aci(PFD=2000,Tleaf=300,ci=seq(40,1500,75),
#' param=f.make.param(TBM='FATES'))$A+rnorm(n = 20,mean = 0,sd = 0.5))
#'
#' f.fitting(measures=data,id.name=NULL,Start=list(JmaxRef=90,VcmaxRef=70,RdRef=1),param=f.make.param(TBM='FATES'))
f.fitting<-function(measures,id.name=NULL,Start=list(JmaxRef=90,VcmaxRef=70,RdRef=1),param=f.make.param(),modify.init=TRUE,do.plot=TRUE,type='Aci'){
  Fixed=param[!names(param)%in%names(Start)]
  if(modify.init){
      if('JmaxRef'%in%names(Start)){Start[['JmaxRef']]=f.modified.arrhenius.inv(P = 6*(max(measures$Photo,na.rm=TRUE)+1),Ha = param[['JmaxHa']],Hd = param[['JmaxHd']],s = param[['JmaxS']],Tleaf = mean(measures$Tleaf,na.rm=TRUE),TRef = param[['TRef']],R = param[['R']])}
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
    if(do.plot){f.plot(measures=measures,name=name,param =param,list_legend = Start,type=type)}
  })

  try({
    Estimation2=mle2(minuslogl = f.MinusLogL,start = Start,fixed = Fixed,data = list(data=measures))
    print(summary(Estimation2))
    #conf=confint(Estimation2)
    #print(conf)
    for(i in names(Estimation2@coef[names(Estimation2@coef)%in%names(param)])){param[i]=Estimation2@coef[i]}
    if(do.plot){f.plot(measures=measures,name=name,param =param,list_legend = as.list(Estimation2@coef),type=type)}
  })
  return(list(MoindresCarres,Estimation2))
}


#' @title Function logit
#' @description This function takes it values in 0;1- and returns values in Inf;+Inf. It is the inverse function of f.logistic
#' @param x
#' @return
#' @export

#' @encoding UTF-8
#' @examples
#' plot(x=seq(0,1,0.01),y=f.logit(x=seq(0,1,0.01)))
f.logit<-function(x){ return(log(x/(1-x)))}



#' @title Logistic function
#' @description This function takes it values in -Inf;+Inf and returns values in 0;1. It is the inverse function of f.logit, ie f.logistic(logit(x))=x
#' @details f.logistic(x)=1/(1+exp(-x)) if x<0, = exp(x)/(1+exp(x)) if x<=0
#' @param x
#' @return
#' @export

#' @examples
#' plot(x=seq(-10,10,0.1),y=f.logistic(x=seq(-10,10,0.1)))
f.logistic<-function(x){ return(ifelse(x>0,1/(1+exp(-x)),exp(x)/(1+exp(x))))}



###pack <- "TestGasEx"
#path <- find.package(pack)
#system(paste(shQuote(file.path(R.home("bin"), "R")),"CMD", "Rd2pdf", shQuote(path)))


