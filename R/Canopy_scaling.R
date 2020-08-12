#' @title Generic function to calculate the GPP within a forest
#' @description Here GPP = sum of Anet at the canopy level, so it takes into account the leaf mitochondrial respiration
#' @param TBM Specific TBM to use (ORCHIDEE, CLM4.5, FATES or JULES)
#' @param meteo_hourly Hourly weather data frame with at least the column at (air temperature in degree C) tl (leaf temperature in degree C) rh (humidity in %) and sr (solar radiation in micro mol photosynthetic photon m-2 s-1)
#' @param Vcmax_Profile Vector of the values of Vcmax at the reference temperature at each layer of the canopy
#' @param Jmax_Profile Vector of the values of Jmax at the reference temperature at each layer of the canopy
#' @param Rd_Profile Vector of the values of Rd at the reference temperature at each layer of the canopy
#' @param Tp_Profile Vector of the values of Tp at the reference temperature at each layer of the canopy
#' @param g0_Profile Vector of the values of g0 at the reference temperature at each layer of the canopy
#' @param g1_Profile Vector of the values of g1 at the reference temperature at each layer of the canopy
#' @param LAI LAI on the ground of the forest
#' @param Height Heigt of the canopy in m
#' @param gsmin Minimum stomatal conductance for water to consider. THis value will be used as the minimum conductance value to avoid 0 and negative values obtained from the coupled assimilation and conductance models
#' @param ... Other parameters of the photosynthetic model, without gradients, for example curvature factor, quantum yield.. see the help of f.make.param()
#' @examples
#' #Example for a simulated day in Panama (lat = 9,..). THe solar radiation is simulated for this example but not likely
#' meteo_hourly=data.frame(at=25,tl=26,rh=80,sr=sin(seq(0,pi,pi/23))*2000,t.d=0:23)
#' meteo_hourly[c(1:7,18:24),'sr']=0
#' canopy=canopy_interception(meteo_hourly=meteo_hourly,lat=9.2801048,t.d=0:23,DOY=20,Height=30,LAI=6,n_layers=50)
#' Vcmax_Profile=rep(65,length(canopy$LAI))
#' Rd_Profile=rep(1.2,length(canopy$LAI))
#' Jmax_Profile=Vcmax_Profile*1.67
#' Tp_Profile=Vcmax_Profile*1/6
#' g0_Profile=rep(0.021,length(canopy$LAI))
#' g1_Profile=rep(4.2,length(canopy$LAI))
#' GPP(TBM = "FATES",meteo_hourly = meteo_hourly,Vcmax_Profile = Vcmax_Profile,Jmax_Profile =Jmax_Profile ,Rd_Profile =Rd_Profile ,Tp_Profile = Tp_Profile,g0_Profile = g0_Profile,g1_Profile = g1_Profile,canopy=canopy,gsmin = 0.01)

#' @return
#' @export
#'
#' @examples
GPP<-function(TBM,meteo_hourly,Vcmax_Profile,Jmax_Profile,Rd_Profile,Tp_Profile,g0_Profile,g1_Profile,gsmin,canopy,Patm=100,...){
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
  figure_conductance_dir=melt(Conductance_rate)
  b=(ggplot(data=figure_conductance_dir,aes(x=time,y=Layer,fill=value))+geom_raster()
     +scale_fill_distiller(palette = "Spectral", direction = -1) +scale_y_reverse()
     +xlab("Time in the day")
     +ylab("Vertical level (0= top, 50 = ground)")
     +labs(fill=expression(g[sw]~(mol~m^-2~s^-1))))
  print(a)
  print(b)
  totalGPP=mean(Photosynthesis_rate,na.rm=TRUE)*max(canopy$LAI)*365*3600*24*44/10^6
  totalET= mean(Trans,na.rm=TRUE)*max(canopy$LAI)*365*3600*24*18*10^-3
  print(paste("GPP = ",totalGPP,"g CO2 m-2 Ground Y-1"))
  print(paste("ET = ",totalET,"L H20 m-2 Ground Y-1"))
  return(list(A=Photosynthesis_rate,gs=Conductance_rate,A_dir=Photosynthesis_rate_dir,gs_dir=gs_dir,A_dif=Photosynthesis_rate_dif,gs_dif=gs_dif,GPP=totalGPP,ET=totalET,fig_A=a,fig_gs=b))
}


#' @title Modeling the canopy structure and interception
#' @details This function is a wrapper of the functions lightME and sunML from bioCro package
#' @param lat Latitude of the modeled canopy (this parameter is important to model the light caracteristics, see lightME)
#' @param t.d Time of the day (0-23) see lightME
#' @param DOY Day Of Year see lightME
#' @param n_layers Number of layers to consider inside the canopy (max 50) see lightME
#' @param Height Total height of the canopy in m
#' @param LAI Total LAI of the canopy
#' @param chi.l Description of the leaves of the canopy (see sunML)
#' @describeIn GPP
#'
#' @return list
#' @export
#'
#' @examples
#' meteo_hourly=data.frame(at=25,tl=26,rh=80,sr=sin(seq(0,pi,pi/23))*2000,t.d=0:23)
#' meteo_hourly[c(1:7,18:24),'sr']=0
#' canopy_interception(meteo_hourly=meteo_hourly,lat=9.2801048,t.d=meteo_hourly$t.d,DOY=20,Height=30,LAI=6)
canopy_interception=function(meteo_hourly,lat,t.d,DOY,n_layers=50,Height,LAI,chi.l=0.9){
  Light_carac=lightME(lat = lat,t.d = t.d,DOY = DOY)# This gives the proportion of diffuse light and direct light
  PFD_dir=meteo_hourly$sr*Light_carac$propIdir
  PFD_dif=meteo_hourly$sr*Light_carac$propIdiff
  cos.th=Light_carac$cos.th

  plot(x=t.d,y=PFD_dir,type="l",xlab="Time of the day",ylab=expression(Light~intensity~(mu~mol~m^-2~s^-1)))
  lines(x=t.d,y=PFD_dif,col="blue")
  legend('topleft',c('Direct light','Diffuse light'),col=c('black','blue'),lty=c(1,1))

  Canopy_time_dir=Canopy_time_dif=f_sun=matrix(data = NA,nrow = n_layers,ncol = length(t.d),dimnames = list(Layer=1:n_layers,time=t.d))
  for(i in 1:length(t.d)){
    Light_Profile=sunML(Idir = PFD_dir[i],Idiff = PFD_dif[i],LAI = LAI,nlayers = n_layers,cos.theta = cos.th[i],heightf = LAI/Height,chi.l=chi.l)
    Canopy_time_dir[,i]=(Light_Profile$layIdir)
    Canopy_time_dif[,i]=(Light_Profile$layIdiff)
    f_sun[,i]=(Light_Profile$layFsun)
  }

  Light=Canopy_time_dir*f_sun+Canopy_time_dif*(1-f_sun)
  Light_Profile$LAI=LAI-Light_Profile$layHeight/Height*LAI
  figure_light_dir=melt(Canopy_time_dir)
  (ggplot(data=figure_light_dir,aes(x=time,y=Layer,fill=value))
    +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
    +scale_y_reverse()
    +labs(fill=expression(Direct~light~(mu~mol~m^-2~s^-1))))

  figure_light_dif=melt(Canopy_time_dif)
  (ggplot(data=figure_light_dif,aes(x=time,y=Layer,fill=value))
    +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
    +scale_y_reverse()
    +labs(fill=expression(Direct~light~(mu~mol~m^-2~s^-1))))

  figure_f_sun=melt(f_sun)
  (ggplot(data=figure_f_sun,aes(x=time,y=Layer,fill=value))+geom_raster()
    +scale_fill_distiller(palette = "Spectral", direction = -1) +scale_y_reverse()
    +ggtitle('Fraction of leaves in direct light')
    +labs(fill='f_sun %'))

  figure_light_tot=melt(Light)
  (ggplot(data=figure_light_tot,aes(x=time,y=Layer,fill=value))
    +ggtitle('Mean PFD of an average leaf')
    +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
    +scale_y_reverse()+labs(fill=expression(PFD~(mu~mol~m^-2~s^-1))))
  return(list(canopy_time_dir=Canopy_time_dir,canopy_time_dif=Canopy_time_dif,f_sun=f_sun,LAI=Light_Profile$LAI,Height=Light_Profile$layHeight))
}
