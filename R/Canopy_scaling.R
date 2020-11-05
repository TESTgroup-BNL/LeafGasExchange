
#' @title Canopy scale GPP calculation
#' @description Generic function to calculate the GPP within a forest (Here GPP = sum of Anet at the canopy level, so it takes into account the leaf mitochondrial respiration)
#' @param TBM Specific TBM to use (ORCHIDEE, CLM4.5, FATES or JULES)
#' @param meteo_hourly Hourly weather data frame with at least the column at (air temperature in degree C) tl (leaf temperature in degree C) rh (humidity in pc) and sr the total PAR in micro mol m-2 s-1
#' @param Vcmax_Profile Vector of the values of Vcmax at the reference temperature at each layer of the canopy
#' @param Jmax_Profile Vector of the values of Jmax at the reference temperature at each layer of the canopy
#' @param Rd_Profile Vector of the values of Rd at the reference temperature at each layer of the canopy
#' @param Tp_Profile Vector of the values of Tp at the reference temperature at each layer of the canopy
#' @param g0_Profile Vector of the values of g0 at the reference temperature at each layer of the canopy
#' @param g1_Profile Vector of the values of g1 at the reference temperature at each layer of the canopy
#' @param gsmin Minimum stomatal conductance for water to consider. This value will be used as the minimum conductance value to avoid 0 and negative values obtained from the coupled assimilation and conductance models
#' @param canopy Description of the canopy interception (see canopy_interception function)
#' @param Patm Atmospheric pressure (used to calculate the transpiration)
#' @param ... Other parameters of the photosynthetic model, without gradients, for example curvature factor, quantum yield.. see the help of f.make.param()
#'
#' @return
#' @export
#'
#' @examples
#'## Simulation of photosynthetic gradients
#' LAI=seq(0,6.2,6.2/49)
#' Vcmax=f.VcmaxRef.LAI(kn=0.11,LAI=LAI,Vcmax0=70)
#' Jmax=1.7*Vcmax; Tp=1/5*Vcmax; Rd=0.03*Vcmax
#' ##Simulation of weather data
#' meteo_hourly=data.frame(time=0:23,rh=80,at=25,sr=sin(seq(0,pi,pi/23))*2000,tl=25)
#' meteo_hourly[!meteo_hourly$time%in%7:17,'sr']=0
#' ##Representation of the light interception inside the canopy
#' canopy=f.canopy.interception(meteo_hourly=meteo_hourly,lat = 9.2801048,t.d = 0:23,DOY = 60,n_layers = 50,Height = 26,LAI = 6)
#' GPP_sc1=f.GPP(TBM = "FATES",meteo_hourly = meteo_hourly,Vcmax_Profile = Vcmax,
#' Jmax_Profile =Jmax ,Rd_Profile =Rd ,Tp_Profile = Tp,
#' g0_Profile = rep(0.02,length(Vcmax)),g1_Profile = rep(4,length(Vcmax)),canopy=canopy,gsmin = 0.01)
f.GPP<-function(TBM,meteo_hourly,Vcmax_Profile,Jmax_Profile,Rd_Profile,Tp_Profile,g0_Profile,g1_Profile,gsmin,canopy,Patm=100,...){
  VpdL_dir=VpdL_dif=Photosynthesis_rate_dir=Photosynthesis_rate_dif=gs_dir=gs_dif=canopy$canopy_time_dir
  for(Layer in 1:nrow(canopy$canopy_time_dir)){
    res_dir=f.A(PFD = canopy$canopy_time_dir[Layer,],
                cs = 400,
                Tair = meteo_hourly[,"at"]+273.15,
                Tleaf= meteo_hourly[,"tl"]+273.15,
                RH = meteo_hourly[,"rh"],
                param = f.make.param(TBM=TBM,
                                     VcmaxRef =Vcmax_Profile[Layer],
                                     RdRef = Rd_Profile[Layer],
                                     JmaxRef=Jmax_Profile[Layer],
                                     TpRef=Tp_Profile[Layer],
                                     g0=g0_Profile[Layer],
                                     g1=g1_Profile[Layer],...
                ))
    ls.gs=which(res_dir$gs<gsmin)
    res_dir$gs[ls.gs]=gsmin
    res_dir$A[ls.gs]=f.A(PFD = 0,cs = 400,Tleaf = meteo_hourly[,"tl"]+273.15,Tair = meteo_hourly[,"at"]+273.15,RH = meteo_hourly[,"rh"],param = f.make.param(TBM=TBM,
                                                                                                                                                             VcmaxRef =Vcmax_Profile[Layer],
                                                                                                                                                             RdRef = Rd_Profile[Layer],
                                                                                                                                                             JmaxRef=Jmax_Profile[Layer],
                                                                                                                                                             TpRef=Tp_Profile[Layer],
                                                                                                                                                             g0=1,
                                                                                                                                                             g1=g1_Profile[Layer]
    ))$A[ls.gs]
    #((-g0_Profile[Layer])*400*sqrt(f.ds(Tleaf = meteo_hourly[,"tl"]+273.15,Tair = meteo_hourly[,"at"]+273.15,RH = meteo_hourly[,"rh"])/1000)/(1.6*g1_Profile[Layer]))[ls.gs]
    Photosynthesis_rate_dir[Layer,]=res_dir$A
    VpdL_dir[Layer,]=res_dir$ds/1000
    gs_dir[Layer,]=res_dir$gs
    res_dif=f.A(PFD = canopy$canopy_time_dif[Layer,],
                cs = 400,
                Tair = meteo_hourly[,"at"]+273.15,
                Tleaf= meteo_hourly[,"tl"]+273.15,
                RH = meteo_hourly[,"rh"],
                param = f.make.param(TBM=TBM,
                                     VcmaxRef =Vcmax_Profile[Layer],
                                     RdRef = Rd_Profile[Layer],
                                     JmaxRef=Jmax_Profile[Layer],
                                     TpRef=Tp_Profile[Layer],
                                     g0=g0_Profile[Layer],
                                     g1=g1_Profile[Layer],...
                ))
    ls.gs=which(res_dif$gs<gsmin)
    res_dif$gs[ls.gs]=gsmin
    res_dif$A[ls.gs]=f.A(PFD = 0,cs = 400,Tleaf = meteo_hourly[,"tl"]+273.15,Tair = meteo_hourly[,"at"]+273.15,RH = meteo_hourly[,"rh"],param = f.make.param(TBM=TBM,
                                                                                                                                                             VcmaxRef =Vcmax_Profile[Layer],
                                                                                                                                                             RdRef = Rd_Profile[Layer],
                                                                                                                                                             JmaxRef=Jmax_Profile[Layer],
                                                                                                                                                             TpRef=Tp_Profile[Layer],
                                                                                                                                                             g0=1,
                                                                                                                                                             g1=g1_Profile[Layer]
    ))$A[ls.gs]
    Photosynthesis_rate_dif[Layer,]=res_dif$A
    gs_dif[Layer,]=res_dif$gs
    VpdL_dif[Layer,]=res_dif$ds/1000
  }
  Photosynthesis_rate=(Photosynthesis_rate_dir*canopy$f_sun+Photosynthesis_rate_dif*(1-canopy$f_sun))
  figure_photosynthesis=melt(Photosynthesis_rate)
  a=(ggplot(data=figure_photosynthesis,aes(x=time,y=Layer,fill=value))+geom_raster()
     +scale_fill_distiller(palette = "Spectral", direction = -1) +scale_y_reverse()
     +xlab("Time in the day")
     +ylab("Vertical level (0= top, 50 = ground)")
     +labs(fill=expression(A~(mu~mol~m^-2~s^-1))))
  
  Conductance_rate=(gs_dir*canopy$f_sun+gs_dif*(1-canopy$f_sun))
  Trans=(gs_dir*VpdL_dir/Patm*canopy$f_sun+gs_dif*VpdL_dif/Patm*(1-canopy$f_sun))
  figure_conductance=melt(Conductance_rate)
  b=(ggplot(data=figure_conductance,aes(x=time,y=Layer,fill=value))+geom_raster()
     +scale_fill_distiller(palette = "Spectral", direction = -1) +scale_y_reverse()
     +xlab("Time in the day")
     +ylab("Vertical level (0= top, 50 = ground)")
     +labs(fill=expression(g[sw]~(mol~m^-2~s^-1))))
  print(a)
  print(b)
  totalGPP=mean(Photosynthesis_rate,na.rm=TRUE)*max(canopy$Light_Profile$LAI)*365*3600*24*44/10^6
  totalET= mean(Trans,na.rm=TRUE)*max(canopy$Light_Profile$LAI)*365*3600*24*18*10^-3
  print(paste("GPP = ",totalGPP,"g CO2 m-2 Ground Y-1"))
  print(paste("ET = ",totalET,"L H20 m-2 Ground Y-1"))
  return(list(A=Photosynthesis_rate,gs=Conductance_rate,A_dir=Photosynthesis_rate_dir,gs_dir=gs_dir,A_dif=Photosynthesis_rate_dif,gs_dif=gs_dif,GPP=totalGPP,ET=totalET,fig_A=a,fig_gs=b))
}

#' @title Canopy scale GPP calculation, with leaf energy budget
#' @description Generic function to calculate the GPP within a forest (Here GPP = sum of Anet at the canopy level, so it takes into account the leaf mitochondrial respiration)
#' @param TBM Specific TBM to use (ORCHIDEE, CLM4.5, FATES or JULES)
#' @param meteo_hourly Hourly weather data frame with at least the column at (air temperature in degree C) rh (humidity in pc) and sr the total PAR in micro mol m-2 s-1
#' @param Vcmax_Profile Vector of the values of Vcmax at the reference temperature at each layer of the canopy
#' @param Jmax_Profile Vector of the values of Jmax at the reference temperature at each layer of the canopy
#' @param Rd_Profile Vector of the values of Rd at the reference temperature at each layer of the canopy
#' @param Tp_Profile Vector of the values of Tp at the reference temperature at each layer of the canopy
#' @param g0_Profile Vector of the values of g0 at the reference temperature at each layer of the canopy
#' @param g1_Profile Vector of the values of g1 at the reference temperature at each layer of the canopy
#' @param gsmin Minimum stomatal conductance for water to consider. This value will be used as the minimum conductance value to avoid 0 and negative values obtained from the coupled assimilation and conductance models
#' @param canopy Description of the canopy interception (see canopy_interception function)
#' @param Patm Atmospheric pressure (used to calculate the transpiration)
#' @param ... Other parameters of the photosynthetic model, without gradients, for example curvature factor, quantum yield.. see the help of f.make.param()
#'
#' @return
#' @export
#'
#' @examples
#'## Simulation of photosynthetic gradients
#' LAI=seq(0,6.2,6.2/49)
#' Vcmax=f.VcmaxRef.LAI(kn=0.11,LAI=LAI,Vcmax0=70)
#' Jmax=1.7*Vcmax; Tp=1/5*Vcmax; Rd=0.03*Vcmax
#' ##Simulation of weather data
#' meteo_hourly=data.frame(time=0:23,rh=80,at=25,sr=sin(seq(0,pi,pi/23))*2000,tl=25,wind=2)
#' meteo_hourly[!meteo_hourly$time%in%7:17,'sr']=0
#' ##Representation of the light interception inside the canopy
#' canopy=f.canopy.interception(meteo_hourly=meteo_hourly,lat = 9.2801048,t.d = 0:23,DOY = 60,n_layers = 50,Height = 26,LAI = 6)
#' GPP_sc1=f.GPPT(TBM = "FATES",meteo_hourly = meteo_hourly,Vcmax_Profile = Vcmax,
#' Jmax_Profile =Jmax ,Rd_Profile =Rd ,Tp_Profile = Tp,
#' g0_Profile = rep(0.02,length(Vcmax)),g1_Profile = rep(4,length(Vcmax)),canopy=canopy,gsmin = 0.01)
f.GPPT<-function(TBM,meteo_hourly,Vcmax_Profile,Jmax_Profile,Rd_Profile,Tp_Profile,g0_Profile,g1_Profile,gsmin,canopy,Patm=100,...){
  VpdL_dir=VpdL_dif=Photosynthesis_rate_dir=Photosynthesis_rate_dif=gs_dir=gs_dif=Tleaf_dir=Tleaf_dif=canopy$canopy_time_dir
  for(Layer in 1:nrow(canopy$canopy_time_dir)){
    print(paste('Layer',Layer,'of', nrow(canopy$canopy_time_dir),'layers'))
    res_dir=f.AT(PFD = canopy$canopy_time_dir[Layer,],
                cs = 400,
                Tair = meteo_hourly[,"at"]+273.15,
                wind= meteo_hourly[,'wind'],
                RH = meteo_hourly[,"rh"],
                param = f.make.param(TBM=TBM,
                                     VcmaxRef =Vcmax_Profile[Layer],
                                     RdRef = Rd_Profile[Layer],
                                     JmaxRef=Jmax_Profile[Layer],
                                     TpRef=Tp_Profile[Layer],
                                     g0=g0_Profile[Layer],
                                     g1=g1_Profile[Layer],...
                ))
    ls.gs=which(res_dir$gs<gsmin)
    res_dir$gs[ls.gs]=gsmin
    res_dir$A[ls.gs]=f.AT(PFD = 0,cs = 400,Tair = meteo_hourly[,"at"]+273.15,RH = meteo_hourly[,"rh"],wind=meteo_hourly[,'wind'],param = f.make.param(TBM=TBM,
                                                                                                                                                             VcmaxRef =Vcmax_Profile[Layer],
                                                                                                                                                             RdRef = Rd_Profile[Layer],
                                                                                                                                                             JmaxRef=Jmax_Profile[Layer],
                                                                                                                                                             TpRef=Tp_Profile[Layer],
                                                                                                                                                             g0=1,
                                                                                                                                                             g1=g1_Profile[Layer]
    ))$A[ls.gs]
    #((-g0_Profile[Layer])*400*sqrt(f.ds(Tleaf = meteo_hourly[,"tl"]+273.15,Tair = meteo_hourly[,"at"]+273.15,RH = meteo_hourly[,"rh"])/1000)/(1.6*g1_Profile[Layer]))[ls.gs]
    Photosynthesis_rate_dir[Layer,]=res_dir$A
    VpdL_dir[Layer,]=res_dir$ds/1000
    gs_dir[Layer,]=res_dir$gs
    Tleaf_dir[Layer,]=res_dir$Tleaf
    res_dif=f.AT(PFD = canopy$canopy_time_dif[Layer,],
                cs = 400,
                Tair = meteo_hourly[,"at"]+273.15,
                wind=meteo_hourly[,'wind'],
                RH = meteo_hourly[,"rh"],
                param = f.make.param(TBM=TBM,
                                     VcmaxRef =Vcmax_Profile[Layer],
                                     RdRef = Rd_Profile[Layer],
                                     JmaxRef=Jmax_Profile[Layer],
                                     TpRef=Tp_Profile[Layer],
                                     g0=g0_Profile[Layer],
                                     g1=g1_Profile[Layer],...
                ))
    ls.gs=which(res_dif$gs<gsmin)
    res_dif$gs[ls.gs]=gsmin
    res_dif$A[ls.gs]=f.AT(PFD = 0,cs = 400,Tair = meteo_hourly[,"at"]+273.15,RH = meteo_hourly[,"rh"],wind=meteo_hourly[,'wind'],param = f.make.param(TBM=TBM,
                                                                                                                                                             VcmaxRef =Vcmax_Profile[Layer],
                                                                                                                                                             RdRef = Rd_Profile[Layer],
                                                                                                                                                             JmaxRef=Jmax_Profile[Layer],
                                                                                                                                                             TpRef=Tp_Profile[Layer],
                                                                                                                                                             g0=1,
                                                                                                                                                             g1=g1_Profile[Layer]
    ))$A[ls.gs]
    Photosynthesis_rate_dif[Layer,]=res_dif$A
    gs_dif[Layer,]=res_dif$gs
    VpdL_dif[Layer,]=res_dif$ds/1000
    Tleaf_dif[Layer,]=res_dif$Tleaf
  }
  Photosynthesis_rate=(Photosynthesis_rate_dir*canopy$f_sun+Photosynthesis_rate_dif*(1-canopy$f_sun))
  figure_photosynthesis=melt(Photosynthesis_rate)
  a=(ggplot(data=figure_photosynthesis,aes(x=time,y=Layer,fill=value))+geom_raster()
     +scale_fill_distiller(palette = "Spectral", direction = -1) +scale_y_reverse()
     +xlab("Time in the day")
     +ylab("Vertical level (0= top, 50 = ground)")
     +labs(fill=expression(A~(mu~mol~m^-2~s^-1))))
  
  Conductance_rate=(gs_dir*canopy$f_sun+gs_dif*(1-canopy$f_sun))
  Trans=(gs_dir*VpdL_dir/Patm*canopy$f_sun+gs_dif*VpdL_dif/Patm*(1-canopy$f_sun))
  Tleaf=(Tleaf_dir*canopy$f_sun+Tleaf_dif*(1-canopy$f_sun))
  figure_conductance=melt(Conductance_rate)
  b=(ggplot(data=figure_conductance,aes(x=time,y=Layer,fill=value))+geom_raster()
     +scale_fill_distiller(palette = "Spectral", direction = -1) +scale_y_reverse()
     +xlab("Time in the day")
     +ylab("Vertical level (0= top, 50 = ground)")
     +labs(fill=expression(g[sw]~(mol~m^-2~s^-1))))
  figure_Tleaf=melt(Tleaf)
  c=(ggplot(data=figure_Tleaf,aes(x=time,y=Layer,fill=value))+geom_raster()
     +scale_fill_distiller(palette = "Spectral", direction = -1) +scale_y_reverse()
     +xlab("Time in the day")
     +ylab("Vertical level (0= top, 50 = ground)")
     +labs(fill=expression(Tleaf~(K))))
  print(a)
  print(b)
  print(c)
  totalGPP=mean(Photosynthesis_rate,na.rm=TRUE)*max(canopy$Light_Profile$LAI)*365*3600*24*44/10^6
  totalET= mean(Trans,na.rm=TRUE)*max(canopy$Light_Profile$LAI)*365*3600*24*18*10^-3
  print(paste("GPP = ",totalGPP,"g CO2 m-2 Ground Y-1"))
  print(paste("ET = ",totalET,"L H20 m-2 Ground Y-1"))
  return(list(A=Photosynthesis_rate,gs=Conductance_rate,A_dir=Photosynthesis_rate_dir,gs_dir=gs_dir,A_dif=Photosynthesis_rate_dif,gs_dif=gs_dif,Tleaf_dir=Tleaf_dir,Tleaf_dif=Tleaf_dif,Tleaf=Tleaf,GPP=totalGPP,ET=totalET,fig_A=a,fig_gs=b,fig_Tleaf=c))
}



#' @title Gradients of photosynthetic parameters
#' @description Several versions of gradients can be found in the litterature, see for example Lloyd et al. 2010 (Fig. 10 and equation A2), but also the equation A14 from Krinner et al. 2005 and the equation 33 from Clark et al. 2011
#' The simpler model describing the gradients is Vcmax(LAI)=Vcmax0 x exp(-kn x LAI) with Vcmax0 Vcmax at the top of the canopy
#' kn can be also calculated as a function of Vcmax0: kn=exp(alpha x Vcmax0+beta)
#' If kn is NULL, then the function will use the default alpha and beta to calculate kn. If, on the contrary, kn is given, this specific one will be used to calculate the gradients.
#' Krinner et al use a slightly different version of this equation with the parameter lambda: Vcmax(LAI)=Vcmax0 x (1-lambda x (1-exp(-kn*LAI))). The previous equation is a particular case of this one for lambda = 1
#' @param alpha Slope of the relationship between Vcmax0 and log(kn), see Lloyd et al. 2010
#' @param beta Intercept of the relationship between Vcmax0 and log(kn), see Lloyd et al. 2010
#' @param Vcmax0 Vcmax at 25 degree C at the top of the canopy
#' @param LAI Vector of Leaf Area Index (or depth within the canopy see Clark et al. 2011)
#' @param kn Exponential decrease
#' @param lambda Asymptot of the decrease (see Krinner et al. 2005)
#' @references Krinner, G., Viovy, N., de Noblet-Ducoudr?, N., Og?e, J., Polcher, J., Friedlingstein, P., . Prentice, I. C. (2005). A dynamic global vegetation model for studies of the coupled atmosphere-biosphere system. Global Biogeochemical Cycles, 19(1). doi:10.1029/2003gb002199
#' Clark, D. B., Mercado, L. M., Sitch, S., Jones, C. D., Gedney, N., Best, M. J., . Cox, P. M. (2011). The Joint UK Land Environment Simulator (JULES), model description - Part 2: Carbon fluxes and vegetation dynamics. Geoscientific Model Development, 4(3), 701-722. doi:10.5194/gmd-4-701-2011
#' Lloyd, J., Pati?o, S., Paiva, R. Q., Nardoto, G. B., Quesada, C. A., Santos, A. J. B., . Mercado, L. M. (2010). Optimisation of photosynthetic carbon gain and within-canopy gradients of associated foliar traits for Amazon forest trees. Biogeosciences, 7(6), 1833-1859. doi:10.5194/bg-7-1833-2010
#' @return Vector of Vcmax at the different LAI specified in the call of the function
#' @export
#'
#' @examples
#' LAI=seq(0,6.2,6.2/49)
#' Vcmax=f.VcmaxRef.LAI(kn=0.11,LAI=LAI,Vcmax0=70)
#' Vcmax2=f.VcmaxRef.LAI(kn=0.11,LAI=LAI,Vcmax0=70,lambda=0.7)
#' plot(Vcmax)
#' lines(Vcmax2)
f.VcmaxRef.LAI=function(alpha=0.00963,beta=-2.43,Vcmax0=50,LAI=0:8,kn=NULL,lambda=1){
  if(is.null(kn)){kn=exp(alpha*Vcmax0+beta)}
  return(Vcmax0*(1-lambda*(1-exp(-kn*LAI))))
}


#' @title Wrapper of biocro lightME and sunML function to describe the light levels inside the canopy
#' @param meteo_hourly Hourly weather data frame with at least the column at (air temperature in degree C) tl (leaf temperature in degree C) rh (humidity in pc) and sr the total PAR in micro mol m-2 s-1
#' @param lat Latitude of the canopy to model (see lightME from biocro)
#' @param t.d time of the day (see lightME from biocro)
#' @param DOY Day of Year (see lightME from biocro)
#' @param n_layers Number of layers inside the canopy (max = 50, see sunML from biocro)
#' @param Height Total height of the canopy (see sunML from biocro)
#' @param LAI Total LAI of the cnaopy (see sunML from biocro)
#' @param chi.l Orientation of the leaves (see sunML from biocro)
#'
#' @return
#' @export
#'
#' @examples
#' ##Simulation of weather data
#' meteo_hourly=data.frame(time=0:23,rh=80,at=25,sr=sin(seq(0,pi,pi/23))*2000,tl=25)
#' meteo_hourly[!meteo_hourly$time%in%7:17,'sr']=0
#' ##Representation of the light interception inside the canopy
#' canopy=f.canopy.interception(meteo_hourly=meteo_hourly,lat = 9.2801048,t.d = 0:23,DOY = 60,n_layers = 50,Height = 26,LAI = 6)
f.canopy.interception=function(meteo_hourly,lat,t.d,DOY,n_layers,Height,LAI,chi.l=0.9){
  Light_carac=lightME(lat = lat,t.d = t.d,DOY = DOY)# This gives the proportion of diffuse light and direct light
  PFD_dir=meteo_hourly$sr*Light_carac$propIdir
  PFD_dif=meteo_hourly$sr*Light_carac$propIdiff
  cos.th=Light_carac$cos.th
  
  plot(x=t.d,y=PFD_dir,type="l",xlab="Time of the day",ylab=expression(Light~intensity~(mu~mol~m^-2~s^-1)))
  lines(x=t.d,y=PFD_dif,col="blue")
  legend('topleft',c('Direct light','Diffuse light'),col=c('black','blue'),lty=c(1,1))
  ### Creation of matrices with 50 vertical layers and 24 hours
  Canopy_time_dir=Canopy_time_dif=Photosynthesis_rate_dir=Photosynthesis_rate_dif=f_sun=Temp_leaf_dir=Temp_leaf_dif=gs_dir=gs_dif=matrix(data = NA,nrow = n_layers,ncol = length(t.d),dimnames = list(Layer=1:n_layers,time=t.d))
  for(i in 1:length(t.d)){
    Light_Profile=sunML(Idir = PFD_dir[i],Idiff = PFD_dif[i],LAI = LAI,nlayers = n_layers,cos.theta = cos.th[i],heightf = LAI/Height,chi.l=chi.l)
    Canopy_time_dir[,i]=(Light_Profile$layIdir)
    Canopy_time_dif[,i]=(Light_Profile$layIdiff)
    f_sun[,i]=(Light_Profile$layFsun)
  }
  
  Light=Canopy_time_dir*f_sun+Canopy_time_dif*(1-f_sun)
  Light_Profile$LAI=LAI-Light_Profile$layHeight/Height*LAI
  figure_light_dir=melt(Canopy_time_dir)
  #print(ggplot(data=figure_light_dir,aes(x=time,y=Layer,fill=value))
  #  +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
  # +scale_y_reverse()
  #  +labs(fill=expression(Direct~light~(mu~mol~m^-2~s^-1))))
  
  figure_light_dif=melt(Canopy_time_dif)
  #print(ggplot(data=figure_light_dif,aes(x=time,y=Layer,fill=value))
  #  +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
  # +scale_y_reverse()
  #  +labs(fill=expression(Diffuse~light~(mu~mol~m^-2~s^-1))))
  
  figure_f_sun=melt(f_sun)
  #print(ggplot(data=figure_f_sun,aes(x=time,y=Layer,fill=value))+geom_raster()
  # +scale_fill_distiller(palette = "Spectral", direction = -1) +scale_y_reverse()
  #  +ggtitle('Fraction of leaves in direct light')
  #+labs(fill='f_sun %'))
  
  figure_light_tot=melt(Light)
  print(ggplot(data=figure_light_tot,aes(x=time,y=Layer,fill=value))
    +ggtitle('Mean PFD of an average leaf')
    +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
    +scale_y_reverse()+labs(fill=expression(PFD~(mu~mol~m^-2~s^-1))))
  return(list(canopy_time_dir=Canopy_time_dir,canopy_time_dif=Canopy_time_dif,f_sun=f_sun,Light_Profile=Light_Profile))
}
