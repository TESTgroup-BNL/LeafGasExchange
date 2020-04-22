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


#' @title USO model for stomatal conductance to water vapour
#' @description Semi-empirical model of the leaf conductance to water vapour
#' @param A Raw assimilation in micromol.m-2.s-1, i-e, the assimilation in presence of respiration
#' @param cs CO2 at the surface of the leaf in ppm
#' @param ds Leaf surface to air vapour pressure deficit in Pa
#' @param g0 Constant of the USO model, representing the conductance when A is 0, in mol.m-2.s-1
#' @param g1 Slope parameter, between 1.14 and 3.58 KPa^0.5 (Wu et al., 2019)
#' @param power Power of the VPDl in USO model. By default is is 0.5 as in Medlin publication
#' @export
#' @return This function returns the optimal stomatal conductance to water vapour in mol.m-2.s-1
#' @references Medlyn, B.E., Duursma, R.A., Eamus, D., Ellsworth, D.S., Colin Prentice, I., Barton, C.V.M., Crous, K.Y., de Angelis, P., Freeman, M. and Wingate, L. (2012), Reconciling the optimal and empirical approaches to modelling stomatal conductance. Glob Change Biol, 18: 3476-3476. doi:10.1111/j.1365-2486.2012.02790.x
#'  Wu, J, Serbin, SP, Ely, KS, et al. The response of stomatal conductance to seasonal drought in tropical forests. Glob Change Biol. 2020; 26: 823– 839. https://doi.org/10.1111/gcb.14820
#'
#' @examples gs=f.gs(A=30,cs=400,ds=1500,g0=0.01,g1=2,power=0.5)
f.gs<-function(A,cs,ds,g0,g1,power=0.5){
  gs=g0+1.6*(1+g1/(ds/1000)^power)*(A)/cs
  return(gs)
}

#' Temperature dependence of Gamma star, Ko, Kc and Rd
#'
#' @param PRef Value of the parameter at the reference temperature
#' @param TRef Reference temperature
#' @param Ha Enthalpie of activation in J.mol-1
#' @param Tleaf Temperature of the leaf in Kelvin
#' @param R Ideal gas constant
#'
#' @return Value of the parameter at the temperature of the leaf
#' @export
#' @references VON CAEMMERER, S. (2013), Steady‐state models of photosynthesis. Plant Cell Environ, 36: 1617-1630. doi:10.1111/pce.12098
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
#' @description
#' @param PFD Photo Flux Density in micromol.m-2,s-1
#' @param cs CO2 at the surface of the leaf in ppm
#' @param Tleaf Temperature of the leaf in Kelvin
#' @param Tair Temperature of the air in Kelvin
#' @param RH Relative Humidity of the air, from 0 to 100
#' @param param List of parameters, see f.make.param for details
#' @export
#' @return List of different variables:
#'  - A: Raw assimilation of the leaf in micromol.m-2.s-1
#'  - gs: Conductance of the leaf for water vapour
#'  - ci: Intracellular CO2 concentration in micromol.mol-1
#'  - ds: Leaf surface to air vapour pressure deficit in Pa
#'
#' @examples f.A(PFD=2000,cs=400,Tleaf=273.16+29,Tair=273.16+28,RH=70)
f.A<-function(PFD,cs,Tleaf,Tair,RH,param=list(R=8.314,
                                              O2=210,
                                              TRef=298.16,
                                              Patm=101,
                                              JmaxRef=	160,
                                              JmaxHa=	50300,
                                              JmaxHd=	152044,
                                              JmaxS	=495,
                                              VcmaxRef=	120,
                                              VcmaxHa	=73637,
                                              VcmaxHd	=149252,
                                              VcmaxS	=486,
                                              RdRef=	1,
                                              RdHa=	46390,
                                              KcRef=	404.9,
                                              KcHa=	79430,
                                              KoRef=	278.4,
                                              KoHa=	36380,
                                              GstarRef=	42.75,
                                              GstarHa	=37830,
                                              abso=	0.85,
                                              f=	0.15,
                                              LogitTheta=f.logit(0.85),
                                              g0=0.01,
                                              g1=2)){

  Kc=f.arrhenius(param[['KcRef']],param[['KcHa']],Tleaf)
  Ko=f.arrhenius(param[['KoRef']],param[['KoHa']],Tleaf)
  Gstar=f.arrhenius(param[['GstarRef']],param[['GstarHa']],Tleaf)

  Rd=f.arrhenius(PRef=param[['RdRef']],param[['RdHa']],Tleaf)
  Vcmax=f.modified.arrhenius(PRef=param[['VcmaxRef']],param[['VcmaxHa']],param[['VcmaxHd']],param[['VcmaxS']],Tleaf)
  Jmax=f.modified.arrhenius(PRef=param[['JmaxRef']],param[['JmaxHa']],param[['JmaxHd']],param[['JmaxS']],Tleaf)

  ds=f.ds(Tleaf,Tair,RH)

  # Analytical solution of the system of equations {E1 : A=f(ci), E2 : gs=f(A,cs) and ci=f(cs)}
  cic=f.solv(x=Vcmax,y=1,z=Kc*(1+param[['O2']]/Ko),cs=cs,Rd=Rd,Gstar=Gstar,g0=param[['g0']],g1=param[['g1']],power=param[['power']],ds=ds)
  Wc=Vcmax*cic/(cic+Kc*(1+param[['O2']]/Ko))

  I2=PFD*param[['abso']]*(1-param[['f']])/2
  J=(I2+Jmax-((I2+Jmax)^2-4*f.logistic(param[['LogitTheta']])*I2*Jmax)^0.5)/(2*f.logistic(param[['LogitTheta']]))
  cij=f.solv(x=J,y=4,z=8*Gstar,cs=cs,Rd=Rd,Gstar=Gstar,g0=param[['g0']],g1=param[['g1']],param[['power']],ds=ds)
  Wj=J*cij/(4*cij+8*Gstar)

  ci=cij

  if(!is.null(which(Wc<Wj))&length(cic)==length(cij)){ci[which(Wc<Wj)]=cic[which(Wc<Wj)]}
  if(!is.null(which(Wc<Wj))&length(cic)!=length(cij)){ci[which(Wc<Wj)]=cic}
  A=pmin(Wc,Wj)*(1-Gstar/ci)-Rd
  gs=f.gs(A=A,cs=cs,ds=ds,g0=param[['g0']],g1=param[['g1']],power=param[['power']])
  output=list(A=A,gs=gs,ci=ci,ds=ds)
  return(output)
}


#' Analytical solution of the coupled photosynthesis and USO model
#'
#' @param x
#' @param y
#' @param z
#' @param cs
#' @param Rd
#' @param Gstar
#' @param g0
#' @param g1
#' @param ds
#'
#' @return
#' @export
#' @keywords internal
#'
#' @examples
f.solv<-function(x,y,z,cs,Rd,Gstar,g0,g1,power,ds){
  m=1.6*(1+g1/(ds/1000)^power)
  a=y*g0+m/cs*(x-Rd*y)
  b=z*g0+m/cs*(-Gstar*x-Rd*z)-cs*g0*y+(x-Rd*y)*(1.6-m)
  c=-z*cs*g0+(1.6-m)*(-Gstar*x-Rd*z)
  ci2=(-b-(b^2-4*a*c)^0.5)/(2*a)
  ci1=(-b+(b^2-4*a*c)^0.5)/(2*a)

  return(ci1)
}

#' @title Photosynthesis and stomata model parameters
#' @description Function to create a list of parameters to be used in f.A, f.Aci, f.AWc and f.AWj function.
#' Depending on the function, all the parameters are not used. For example go and g1 are not used in f.Aci
#'
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
#' @param f Correcting factor for the spectral quality of the light
#' @param LogitTheta Theta is the empirical curvacture factor for the response of J to PFD. It takes its values between 0 and 1. To avoid numerical issues when fitting data, this parameters is transformed in this model and called LogitTheta. The transformation between Theta and LogitTheta is: Theta=f.logistic(LogitTheta) and LogitTheta=f.logit(Theta)
#' @param g0 Constant of the USO model, representing the conductance when A is 0, in mol.m-2.s-1
#' @param g1 Slope parameter, between 1.14 and 3.58 KPa^0.5 (Wu et al., 2019)
#' @param power Power of VPDl in USO model. By default power=0.5 as in Medlyn article
#'
#' @return List of parameters that can be used in f.A
#' @references Bernacchi, C.J., Singsaas, E.L., Pimentel, C., Portis Jr, A.R. and Long, S.P. (2001), Improved temperature response functions for models of Rubisco‐limited photosynthesis. Plant, Cell & Environment, 24: 253-259. doi:10.1111/j.1365-3040.2001.00668.x
#' @export
#'
#' @examples param1=f.make.param(JmaxRef=100,VcmaxRef=60,RdRef=1)
#'  param2=f.make.param(JmaxRef=100,VcmaxRef=80,RdRef=1)
#'  f.A(PFD=1500,cs=400,Tleaf=300,Tair=299,RH=70,param=param1)
#'  f.A(PFD=1500,cs=400,Tleaf=300,Tair=299,RH=70,param=param2)
f.make.param<-function(R=8.314,
                       O2=210,
                       TRef=298.16,
                       Patm=101,
                       JmaxRef=	85,
                       JmaxHa=	50300,
                       JmaxHd=	152044,
                       JmaxS	=495,
                       VcmaxRef=	55,
                       VcmaxHa	=73637,
                       VcmaxHd	=149252,
                       VcmaxS	=486,
                       RdRef=	1,
                       RdHa=	46390,
                       KcRef=	404.9,
                       KcHa=	79430,
                       KoRef=	278.4,
                       KoHa=	36380,
                       GstarRef=	42.75,
                       GstarHa	=37830,
                       abso=	0.85,
                       f=	0.15,
                       LogitTheta=f.logit(0.85),
                       g0=0.01,
                       g1=2,
                       power=0.5){
  return(list(R=R,
              O2=O2,
              TRef=TRef,
              Patm=Patm,
              JmaxRef=JmaxRef,
              JmaxHa=	JmaxHa,
              JmaxHd=	JmaxHd,
              JmaxS	=JmaxS,
              VcmaxRef=VcmaxRef,
              VcmaxHa	= VcmaxHa,
              VcmaxHd	=VcmaxHd,
              VcmaxS	=VcmaxS,
              RdRef=RdRef,
              RdHa=	 RdHa,
              KcRef= KcRef,
              KcHa=	KcHa,
              KoRef=KoRef,
              KoHa=	KoHa,
              GstarRef=	GstarRef,
              GstarHa	=GstarHa,
              abso=	abso,
              f=	f,
              LogitTheta=LogitTheta,
              g0=g0,
              g1=g1,
              power=power))

}

#' @title Photosynthesis model
#' @description Calculate the assimilation according to Farquhar equations. Contrary to f.A, this function uses intracellular CO2 and not ambiant air CO2
#' @param PFD Photo Flux Density in micromol.m-2,s-1
#' @param ci Leaf intracellular CO2 in ppm
#' @param Tleaf Temperature of the leaf in Kelvin
#' @param param List of parameters, see f.make.param for details
#'
#' @return Assimilation in micromol.m-2.s-1
#' @export
#'
#' @examples ci=seq(1,1500,10)
#' plot(x=ci,y=f.Aci(PFD=2000,ci=ci,Tleaf=300))
f.Aci=function(PFD,ci,Tleaf,param=list(R=8.314,
                                       O2=210,
                                       TRef=298.16,
                                       Patm=101,
                                       JmaxRef=	160,
                                       JmaxHa=	50300,
                                       JmaxHd=	152044,
                                       JmaxS	=495,
                                       VcmaxRef=	120,
                                       VcmaxHa	=73637,
                                       VcmaxHd	=149252,
                                       VcmaxS	=486,
                                       RdRef=	1,
                                       RdHa=	46390,
                                       KcRef=	404.9,
                                       KcHa=	79430,
                                       KoRef=	278.4,
                                       KoHa=	36380,
                                       GstarRef=	42.75,
                                       GstarHa	=37830,
                                       abso=	0.85,
                                       f=	0.15,
                                       LogitTheta=f.logit(0.85))){

  Kc=f.arrhenius(param[['KcRef']],param[['KcHa']],Tleaf)
  Ko=f.arrhenius(param[['KoRef']],param[['KoHa']],Tleaf)
  Gstar=f.arrhenius(param[['GstarRef']],param[['GstarHa']],Tleaf)

  Rd=f.arrhenius(PRef=param[['RdRef']],param[['RdHa']],Tleaf)
  Vcmax=f.modified.arrhenius(PRef=param[['VcmaxRef']],param[['VcmaxHa']],param[['VcmaxHd']],param[['VcmaxS']],Tleaf)
  Jmax=f.modified.arrhenius(PRef=param[['JmaxRef']],param[['JmaxHa']],param[['JmaxHd']],param[['JmaxS']],Tleaf)

  Wc=Vcmax*ci/(ci+Kc*(1+param[['O2']]/Ko))

  I2=PFD*param[['abso']]*(1-param[['f']])/2
  J=(I2+Jmax-((I2+Jmax)^2-4*f.logistic(param[['LogitTheta']])*I2*Jmax)^0.5)/(2*f.logistic(param[['LogitTheta']]))
  Wj=J*ci/(4*ci+8*Gstar)

  A=pmin(Wc,Wj)*(1-Gstar/ci)-Rd
  return(A)
}


#' @title Carbon assimilation under electron transport limitation
#' @inheritParams f.A
#'
#' @return Assimilation under Jmax limitation
#' @export
#'
#' @examples plot(x=seq(0,1500,10),y=f.AWj(PFD=seq(0,1500,10),ci=270,Tleaf=300,param=f.make.param()),xlab='PFD',ylab='Awj in micromol.m-2.s-1')

f.AWj<-function(PFD,ci,Tleaf,param=list(RdHa,GstarRef,GstarHa,JmaxHa,JmaxHd,JmaxS,
                JmaxRef,RdRef,abso,f,LogitTheta)){

  Gstar=f.arrhenius(param[['GstarRef']],param[['GstarHa']],Tleaf)
  Rd=f.arrhenius(PRef=param[['RdRef']],param[['RdHa']],Tleaf)
  Jmax=f.modified.arrhenius(PRef=param[['JmaxRef']],param[['JmaxHa']],param[['JmaxHd']],param[['JmaxS']],Tleaf)

  I2=PFD*param[['abso']]*(1-param[['f']])/2
  J=(I2+Jmax-((I2+Jmax)^2-4*f.logistic(param[['LogitTheta']])*I2*Jmax)^0.5)/(2*f.logistic(param[['LogitTheta']]))
  Wj=J*ci/(4*ci+8*Gstar)

  A=Wj*(1-Gstar/ci)-Rd
  return(A)
}


#' @title Carbon assimilation under rubisco carboxylation limitation
#'
#' @inheritParams f.A
#' @return
#' @export
#'
#' @examples plot(x=seq(1,1500,10),y=f.AWc(ci=seq(1,1500,10),Tleaf=300,param=f.make.param()),xlab='Intracellular CO2 in ppm',ylab='Awc in micromol.m-2.s-1')
f.AWc<-function(ci,Tleaf,param=list(ci,Tleaf,RdHa,GstarRef,GstarHa,RdRef,VcmaxRef,VcmaxHa,VcmaxHd,VcmaxS,KoRef,KoHa,KcRef,KcHa)){
  Kc=f.arrhenius(param[['KcRef']],param[['KcHa']],Tleaf)
  Ko=f.arrhenius(param[['KoRef']],param[['KoHa']],Tleaf)
  Gstar=f.arrhenius(param[['GstarRef']],param[['GstarHa']],Tleaf)

  Rd=f.arrhenius(PRef=param[['RdRef']],param[['RdHa']],Tleaf)
  Vcmax=f.modified.arrhenius(PRef=param[['VcmaxRef']],param[['VcmaxHa']],param[['VcmaxHd']],param[['VcmaxS']],Tleaf)
  Wc=Vcmax*ci/(ci+Kc*(1+param[['O2']]/Ko))

  A=Wc*(1-Gstar/ci)-Rd
  return(A)
}


#' @title Intracellular CO2 threshold between electron transport and carboxylation limitations
#' @inheritParams f.A
#' @return Intracellular CO2 such as Wc==Wj
#' @export
#'
#' @examples f.ci.treshold(PFD=2000,Tleaf=300,param=f.make.param(VcmaxRef=60,JmaxRef=85))
#' @examples f.ci.treshold(PFD=2000,Tleaf=300,param=f.make.param(VcmaxRef=70,JmaxRef=85))
f.ci.treshold<-function(PFD,Tleaf,param=list(GstarRef,GstarHa,KoRef,
                       KoHa,KcRef,KcHa,VcmaxHa,VcmaxHd,VcmaxS,JmaxHa,JmaxHd,JmaxS,
                       VcmaxRef,JmaxRef,abso,f,LogitTheta)){

  Kc=f.arrhenius(param[['KcRef']],param[['KcHa']],Tleaf)
  Ko=f.arrhenius(param[['KoRef']],param[['KoHa']],Tleaf)
  Gstar=f.arrhenius(param[['GstarRef']],param[['GstarHa']],Tleaf)

  Vcmax=f.modified.arrhenius(PRef=param[['VcmaxRef']],param[['VcmaxHa']],param[['VcmaxHd']],param[['VcmaxS']],Tleaf)
  Jmax=f.modified.arrhenius(PRef=param[['JmaxRef']],param[['JmaxHa']],param[['JmaxHd']],param[['JmaxS']],Tleaf)

  I2=PFD*param[['abso']]*(1-param[['f']])/2
  J=(I2+Jmax-((I2+Jmax)^2-4*f.logistic(param[['LogitTheta']])*I2*Jmax)^0.5)/(2*f.logistic(param[['LogitTheta']]))
  ci_t=(J*Kc*(1+param[['O2']]/Ko)-8*Gstar*Vcmax)/(4*Vcmax-J)
  return(ci_t)
}

#' @title Plot Aci data and model
#' @description Plot a generic graphic with observed data and predictions. Be careful to sort the data.frame beforehand.
#' @param measures Data frame obtained from CO2 curve with at least columns Photo, Ci, PARi and Tleaf
#' @inheritParams f.A
#' @param name Name of the curve to be displayed
#'
#' @return Plot a figure
#' @export
#'
#' @examples
#' param=f.make.param()
#' Photo=f.Aci(PFD=2000,Tleaf=300,ci=seq(40,1500,50),param=param)+rnorm(n = 30,mean = 0,sd = 0.5)
#' data=data.frame(Tleaf=rep(300,30),Ci=seq(40,1500,50),PARi=rep(2000,30),Photo=Photo)
#' f.plot.Aci(measures=data,param=param,list_legend=param['VcmaxRef'],name='Example 01')

f.plot.Aci<-function(measures=NULL,list_legend,param,name=''){
  # Plot all data points

  plot(x=measures$Ci,y=measures$Photo, main=name, xlab="Ci in ppm", ylab="Photo in micromol.m-2.s-1",ylim=c(min(measures$Photo,na.rm = TRUE),1.15*max(measures$Photo,na.rm = TRUE)))
  if(!is.null(list_legend)){
    list_legend=list_legend[order(names(list_legend))]
    legend("bottomright",legend=mapply(FUN = function(x, i){paste(i,'=', round(x,2))}, list_legend, names(list_legend)),bty="n",cex=1)
  }
  legend("topleft",legend=c("Rubisco","RuBP","Photo","Obs"),lty=c(2,2,1,0),
         pch=c(NA,NA,NA,21),
         col=c("dark blue","dark red","dark grey","black"),bty="n",lwd=c(2,2,1,1),
         seg.len=2,cex=1,pt.cex=1)
  lines(x=measures$Ci,y=f.Aci(ci=measures$Ci,Tleaf=measures$Tleaf,PFD=measures$PARi,param=param),col="dark grey",lwd=1)
  lines(x=measures$Ci,y=f.AWc(ci=measures$Ci,Tleaf=measures$Tleaf,param=param),lwd=2,col="dark blue",lty=2)
  lines(x=measures$Ci,y=f.AWj(ci=measures$Ci,Tleaf=measures$Tleaf,PFD=measures$PARi,param=param),lwd=2,col="dark red",lty=2)
   box(lwd=1)
}

#' @title Plot AQ data and model
#' @description Plot a generic graphic with observed data and predictions. Be careful to sort the data.frame beforehand.
#' @inheritParams f.plot.Aci
#'
#' @return Plot a figure
#' @export
#'
#' @examples
#' param=f.make.param()
#' data=data.frame(Tleaf=seq(298,305,0.24),Ci=seq(300,271),PARi=seq(0,2000,67),Photo=f.AWj(PFD=seq(0,2000,67),Tleaf=seq(298,305,0.24),ci=seq(300,271),param=param)+rnorm(n = 30,mean = 0,sd = 0.5))
#' f.plot.AQ(measures=data,param=param,list_legend=param[c('VcmaxRef','RdRef','JmaxRef')],name='Example 01')

f.plot.AQ<-function(measures=NULL,param,list_legend=NULL,name=''){
  # Plot all data points

  plot(x=measures$PARi,y=measures$Photo, main=name, xlab="PFD in micromol.m-2.s-1", ylab="Photo in micromol.m-2.s-1",ylim=c(min(measures$Photo,na.rm = TRUE),1.15*max(measures$Photo,na.rm = TRUE)))
  if(!is.null(list_legend)){
      list_legend=list_legend[order(names(list_legend))]
      legend("bottomright",legend=mapply(FUN = function(x, i){paste(i,'=', round(x,2))}, list_legend, names(list_legend)),bty="n",cex=1)
  }
  legend("topleft",legend=c("RuBP","Obs"),lty=c(2,0),
         pch=c(NA,21),
         col=c("dark blue","black"),bty="n",lwd=c(2,1),
         seg.len=2,cex=1,pt.cex=1)
  lines(x=measures$PARi,y=f.AWj(ci=measures$Ci,Tleaf=measures$Tleaf,PFD=measures$PARi,param=param),lwd=2,col="dark blue",lty=2)
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
  y<-data$Photo-f.Aci(ci=data$Ci,PFD=data$PARi,Tleaf=data$Tleaf,param=param)
  return(sum(y^2))
}

#' @title Compute the sum square of the difference between obervations and predictions
#' @description Function used to fit the parameters of a CO2 curve
#' @param x List of parameters to fit
#' @param data Data frame obtained from CO2 curve with at least columns Photo, Ci, PARi and Tleaf
#'
#' @return Sum square of the difference between predictions and observations
#' @export
#' @keywords internal
#'
#' @examples
f.SumSq_Aq<-function(Fixed,data,Start){
  param=c(Fixed,Start)
  y<-data$Photo-f.AWj(ci=data$Ci,PFD=data$PARi,Tleaf=data$Tleaf,param=param)
  return(sum(y^2))
}

#' Title
#'
#' @inheritParams f.make.param
#' @inheritParams f.plot.Aci
#' @param sigma Sigma value
#' @return
#' @export
#' @keywords internal
#'
#' @examples
f.MinusLogL<-function(data,sigma,R=8.314,
                       O2=210,
                       TRef=298.16,
                       Patm=101,
                       JmaxRef=	85,
                       JmaxHa=	50300,
                       JmaxHd=	152044,
                       JmaxS	=495,
                       VcmaxRef=	55,
                       VcmaxHa	=73637,
                       VcmaxHd	=149252,
                       VcmaxS	=486,
                       RdRef=	1,
                       RdHa=	46390,
                       KcRef=	404.9,
                       KcHa=	79430,
                       KoRef=	278.4,
                       KoHa=	36380,
                       GstarRef=	42.75,
                       GstarHa	=37830,
                       abso=	0.85,
                       f=	0.15,
                       LogitTheta=0.85,
                       g0=0.01,
                       g1=2,
                      power=0.5){
  param=list(R=R,
             O2=O2,
             TRef=TRef,
             Patm=Patm,
             JmaxRef=JmaxRef,
             JmaxHa=	JmaxHa,
             JmaxHd=	JmaxHd,
             JmaxS	=JmaxS,
             VcmaxRef=VcmaxRef,
             VcmaxHa	= VcmaxHa,
             VcmaxHd	=VcmaxHd,
             VcmaxS	=VcmaxS,
             RdRef=RdRef,
             RdHa=	 RdHa,
             KcRef= KcRef,
             KcHa=	KcHa,
             KoRef=KoRef,
             KoHa=	KoHa,
             GstarRef=	GstarRef,
             GstarHa	=GstarHa,
             abso=	abso,
             f=	f,
             LogitTheta=LogitTheta,
             g0=g0,
             g1=g1,
             power=power)
  y<-dnorm(x=data$Photo,mean=f.Aci(ci=data$Ci,PFD=data$PARi,Tleaf=data$Tleaf,param=param),sd=sigma,log=TRUE)
  return(-sum(y))
}

#' Title
#'
#' @inheritParams f.make.param
#' @inheritParams f.plot.Aci
#' @param sigma Sigma value
#' @return
#' @export
#' @keywords internal
#'
#' @examples
f.MinusLogL_Aq<-function(data,sigma,R=8.314,
                      O2=210,
                      TRef=298.16,
                      Patm=101,
                      JmaxRef=	85,
                      JmaxHa=	50300,
                      JmaxHd=	152044,
                      JmaxS	=495,
                      VcmaxRef=	55,
                      VcmaxHa	=73637,
                      VcmaxHd	=149252,
                      VcmaxS	=486,
                      RdRef=	1,
                      RdHa=	46390,
                      KcRef=	404.9,
                      KcHa=	79430,
                      KoRef=	278.4,
                      KoHa=	36380,
                      GstarRef=	42.75,
                      GstarHa	=37830,
                      abso=	0.85,
                      f=	0.15,
                      LogitTheta=f.logit(0.85),
                      g0=0.01,
                      g1=2,
                      power=0.5){
  param=list(R=R,
             O2=O2,
             TRef=TRef,
             Patm=Patm,
             JmaxRef=JmaxRef,
             JmaxHa=	JmaxHa,
             JmaxHd=	JmaxHd,
             JmaxS	=JmaxS,
             VcmaxRef=VcmaxRef,
             VcmaxHa	= VcmaxHa,
             VcmaxHd	=VcmaxHd,
             VcmaxS	=VcmaxS,
             RdRef=RdRef,
             RdHa=	 RdHa,
             KcRef= KcRef,
             KcHa=	KcHa,
             KoRef=KoRef,
             KoHa=	KoHa,
             GstarRef=	GstarRef,
             GstarHa	=GstarHa,
             abso=	abso,
             f=	f,
             LogitTheta=LogitTheta,
             g0=g0,
             g1=g1,
             power=power)
  y<-dnorm(x=data$Photo,mean=f.AWj(ci=data$Ci,PFD=data$PARi,Tleaf=data$Tleaf,param=param),sd=sigma,log=TRUE)
  return(-sum(y))
}

#' @title Fitting function for Aci data
#' @description Function to fit f.Aci model to data. The parameters to fit have to be described in the list Start.
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
#' param=f.make.param(RdRef=1.25,VcmaxRef=57,JmaxRef=92))+rnorm(n = 20,mean = 0,sd = 0.5))
#'
#' f.CO2.fitting(measures=data,id.name=NULL,Start=list(JmaxRef=90,VcmaxRef=70,RdRef=1))
f.CO2.fitting<-function(measures,id.name=NULL,Start=list(JmaxRef=90,VcmaxRef=70,RdRef=1),param=f.make.param(),modify.init=TRUE,do.plot=TRUE){
  Fixed=param[!names(param)%in%names(Start)]
  if(modify.init){
      if('JmaxRef'%in%names(Start)){Start[['JmaxRef']]=f.modified.arrhenius.inv(P = 6*(max(measures$Photo,na.rm=TRUE)+1),Ha = param[['JmaxHa']],Hd = param[['JmaxHd']],s = param[['JmaxS']],Tleaf = mean(measures$Tleaf,na.rm=TRUE),TRef = param[['TRef']],R = param[['R']])}
      if('JmaxRef'%in%names(Start)&'VcmaxRef'%in%names(Start)){Start[['VcmaxRef']]=Start[['JmaxRef']]/2}
      grille=expand.grid(lapply(X = Start,FUN = function(x){x*c(0.2,0.6,1,1.5,2)}))
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
    if(do.plot){f.plot.Aci(measures=measures,name=name,param =param,list_legend = Start)}
  })

  try({
    Estimation2=mle2(minuslogl = f.MinusLogL,start = Start,fixed = Fixed,data = list(data=measures))
    print(summary(Estimation2))
    print(confint(Estimation2))
    if(do.plot){f.plot.Aci(measures=measures,name=name,param =param,list_legend = as.list(Estimation2@coef))}
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

#' Fitting function for AQ data
#'
#' @inheritParams f.CO2.fitting
#'
#' @return
#' @export
#'
#' @examples
#' data=data.frame(Tleaf=300,Ci=280,PARi=seq(0,2000,67),Photo=f.AWj(PFD=seq(0,2000,67),Tleaf=300,ci=280,param=f.make.param(LogitTheta=f.logit(0.89),JmaxRef=215,RdRef=1))+rnorm(n = 30,mean = 0,sd = 0.35))
#' Start=list(JmaxRef=50,LogitTheta=f.logit(0.8),RdRef=0.5)
#' f.light.fitting(measures=data,Start=Start)
f.light.fitting<-function(measures,id.name=NULL,Start=list(JmaxRef=90,LogitTheta=0.6,RdRef=1),param=f.make.param(),modify.init=TRUE,do.plot=TRUE){
  Fixed=param[!names(param)%in%names(Start)]
  param=c(Start,Fixed)
  if(modify.init){
     ##Find reasonable start values for RdRef and JmaxRef if they are in Start
      if('RdRef'%in%names(Start)&0%in%round(measures$PARi,-1)){
          Rd.Est=abs(mean(measures[round(measures$PARi,-1)==0,'Photo'],na.rm = TRUE))
          Tleaf.Est=mean(measures$Tleaf,na.rm=TRUE)
          Start[['RdRef']]=f.arrhenius.inv(P = Rd.Est,Ha = param[['RdHa']],Tleaf = Tleaf.Est,TRef = param[['TRef']],R = param[["R"]])}
      if('JmaxRef'%in%names(Start)&0%in%round(measures$PARi,-1)){Start['JmaxRef']=f.modified.arrhenius.inv(P = 6*(max(measures$Photo,na.rm=TRUE)+f.arrhenius(PRef = param[['RdRef']],Ha = param[['RdHa']],Tleaf = mean(measures$Tleaf,na.rm=TRUE),TRef = param[['TRef']],R = param[['R']])),Ha = param[['JmaxHa']],Hd = param[['JmaxHd']],s = param[['JmaxS']],Tleaf = mean(measures$Tleaf,na.rm=TRUE),TRef = param[['TRef']],R = param[['R']])}
      grille=expand.grid(lapply(X = Start,FUN = function(x){x*c(0.5,1,1.5,2)}))
      grille.list=apply(X=grille,MARGIN = 1,FUN=as.list)
      value=9999999
      l_param=0
      for(l in 1:nrow(grille)){
          MoindresCarres=optim(par=grille.list[[l]],fn=f.SumSq_Aq,data=measures,Fixed=Fixed,control = list(parscale=unlist(Start),maxit=1000))
        if(!is.null(MoindresCarres)&MoindresCarres$value<value){value=MoindresCarres$value;l_param=l}
      }
      Start=grille.list[[l]]
  }
  if(is.null(id.name)){name=''}else{name=unique(measures[,id.name])}
  MoindresCarres<-Estimation2<-NULL
  try({
    MoindresCarres<-optim(par=Start,fn=f.SumSq_Aq,data=measures,Fixed=Fixed,control = list(parscale=unlist(Start),maxit=2000))
    print(MoindresCarres)
    print(paste('sd',sqrt(MoindresCarres$value/NROW(measures))))
    Start$sigma=sqrt(MoindresCarres$value/NROW(measures))
    for(l.name in names(MoindresCarres$par)){Start[l.name]=MoindresCarres$par[[l.name]]}
    for(l.name in names(MoindresCarres$par)){param[l.name]=MoindresCarres$par[[l.name]]}
    if(do.plot){f.plot.AQ(measures=measures,name=name,param =param,list_legend = Start)}

  })

  try({
    Estimation2=mle2(minuslogl = f.MinusLogL_Aq,start = Start,fixed = Fixed,data = list(data=measures),control = list(parscale=unlist(Start),maxit=2000))
    print(summary(Estimation2))
    print(confint(Estimation2))
    param =as.list(Estimation2@fullcoef )
    if(do.plot){f.plot.AQ(measures=measures,name=name,param =param,list_legend = as.list(Estimation2@coef ))}
  })
return(list(MoindresCarres,Estimation2))
}

###pack <- "TestGasEx"
#path <- find.package(pack)
#system(paste(shQuote(file.path(R.home("bin"), "R")),"CMD", "Rd2pdf", shQuote(path)))


