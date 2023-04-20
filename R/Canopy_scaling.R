#' Tridiagonal solver
#' Converted into a R code from the original code of Gordon Bonan: Bonan, G. (2019). Climate Change and Terrestrial Ecosystem Modeling. Cambridge: Cambridge University Press. doi:10.1017/9781107339217
#' @title Tridiagonal solver
#' @description 
#' % Solve for U given the set of equations R * U = D, where U is a vector
#' of length N, D is a vector of length N, and R is an N x N tridiagonal
#' matrix defined by the vectors A, B, C each of length N. A(1) and
#' C(N) are undefined and are not referenced.
#'
#'     |B(1) C(1) ...  ...  ...                     |
#'     |A(2) B(2) C(2) ...  ...                     |
#' R = |     A(3) B(3) C(3) ...                     |
#'     |                    ... A(N-1) B(N-1) C(N-1)|
#'     |                    ... ...    A(N)   B(N)  |
#'
#' The system of equations is written as:
#'
#'    A_i * U_i-1 + B_i * U_i + C_i * U_i+1 = D_i
#'
#' for i = 1 to N. The solution is found by rewriting the
#' equations so that:
#'
#'    U_i = F_i - E_i * U_i+1
  
#' @param a See description.
#' @param b See description.
#' @param c See description.
#' @param d See description.
#' @param n See description.
#'
#' @return Solution U
#' @export
#'
#' @examples
f.tridiagonal.solver=function(a,b,c,d,n){
  e=rep(NA,n-1)
  e[1] = c[1] / b[1]
  
  for (i in 2:(n-1)){
    e[i] = c[i] / (b[i] - a[i] * e[i-1])
  } 

  f=rep(NA,n)
  f[1] = d[1] / b[1]
  
  for (i in 2:(n)){
    f[i] = (d[i] - a[i] * f[i-1]) / (b[i] - a[i] * e[i-1])
  }
  
  u=rep(NA,n)
  u[n] = f[n];
  
  for (i in seq(n-1,1,-1)){
    u[i] = f[i] - e[i] * u[i+1];
  }
  return(u)
}

#' @title Norman 1979 Radiation interception model
#' Converted into a R code from the original code of Gordon Bonan: Bonan, G. (2019). Climate Change and Terrestrial Ecosystem Modeling. Cambridge: Cambridge University Press. doi:10.1017/9781107339217
#' @param Rho Leaf reflectance.
#' @param Tau Leaf transmittance.
#' @param Rho_soil_dir Direct beam albedo of ground (soil).
#' @param Rho_soil_dif Direct beam albedo of ground (soil).
#' @param cosz Cosinus of the solar zenith angle.
#' @param chil Index of departure of the leaf angles from a spherical distribution. -0.4 < chil < 0.6.
#' @param clumpfac Clumping factor, index of non random spatial distribution of leaves. = 1 for randomly spaced leaves, <1 for clumed leaves (Chen et al. 2012).
#' @param dLAI LAI of each one of the n layers of vegetation in the canopy, layer 1 is the top of canopy, layer n is the bottom.
#' @param nlayers Number of vegetation layers.
#' @param PARdir Atmospheric direct beam solar radiation (W/m2).
#' @param PARdif Atmospheric diffuse solar radiation (W/m2).
#'
#' @return list of output:
#' PARsun Absorbed PFD by the sunlit leaves
#' PARsha Absorbed PFD by the shaded leaves
#' fracsun Proportion of sunlit leaves
#' fracsha Proportion of shaded leaves
#' @export
#'
#' @examples
#' f.Norman.Radiation(Rho=0.1, Tau=0.05, PARdir=1000,PARdif=200,dLAI=c(rep(6/20,20)),nlayers=20,Rho_soil_dif = 0.1,Rho_soil_dir = 0.1,cosz = 0.88,chil = 0.1,clumpfac = 0.8)

f.Norman.Radiation=function(Rho=0.1, Tau=0.05, Rho_soil_dir=0.1,Rho_soil_dif=0.1,cosz,chil,clumpfac,dLAI,nlayers,PARdir=0.8,PARdif=0.2){
  if(length(dLAI)!=nlayers){print('Error: the input parameters nlayers does not correspond to the length of the input vector dLAI')}
  ## be careful, the original code by Bonan works with the first layer being the ground and the last being the top.
  ## so dlai has to be reverted.
  dLAI=rev(dLAI)
  lai=sum(dLAI)
  sumlai=c(NA,lai-cumsum(dLAI)+dLAI/2)
  dLAI=c(NA,dLAI)
  
  if(chil>0.6|chil<(-0.4)){print('Chil is not inside the interval -0.4, 0.6 and was changed')}
  chil = min(max(chil, -0.4), 0.6)
  phi1 = 0.5 - 0.633 * chil - 0.330 * chil^2
  phi2 = 0.877 * (1 - 2 * phi1)
  
  gdir = phi1 + phi2 * cosz
  
  # Direct beam extinction coefficient
  
  Kb = gdir / cosz;
  
  # Prevent large Kb at low sun angle
  
  Kb = min(Kb, 20)
  
  fracsun = clumpfac * exp(-Kb * sumlai *clumpfac);
  fracsha = 1 - fracsun
  
  laisun = (1 - exp(-Kb * lai * clumpfac)) / Kb
  laisha = lai - laisun
  
  tb = exp(-Kb * dLAI * clumpfac)
  td=rep(0,length(dLAI))
  for (j in 1:9){
    angle = (5 + (j - 1) * 10) * pi / 180
    gdirj = phi1 + phi2 * cos(angle)
    td = td+exp(-gdirj / cos(angle) *dLAI * clumpfac) * sin(angle) * cos(angle)
  }
  td = td * 2 * (10 * pi / 180)
  
  tbcum=rep(NA,nlayers+1)
  cumlai = 0
  iv = nlayers+1
  tbcum[iv] = 1
  for (iv in (nlayers+1):2){
    cumlai = cumlai + dLAI[iv]
    tbcum[iv-1] = exp(-Kb * cumlai * clumpfac);
  }
  
  
  print(paste('Radiation model for a total LAI of ', lai))
  
  swup=swdn=rep(0,nlayers+1)
  a=b=c=d=rep(0,nlayers+1)
  omega=Rho+Tau
  # Soil: upward flux
  m=1
  iv = 1
  a[m] = 0
  b[m] = 1
  c[m] = -Rho_soil_dif
  d[m] = PARdir * tbcum[m] * Rho_soil_dir
  
  # Soil: downward flux
  
  refld = (1 - td[iv+1]) * Rho
  trand = (1 - td[iv+1]) * Tau + td[iv+1]
  aiv = refld - trand * trand / refld;
  biv = trand / refld;
  
  m = 2
  a[m] = -aiv
  b[m] = 1
  c[m] = -biv
  d[m] = PARdir * tbcum[iv+1] * (1 - tb[iv+1]) * (Tau - Rho * biv)
  
  # Leaf layers, excluding top layer
  
  for (iv in 2:(nlayers)){
    # Upward flux
    
    refld = (1 - td[iv]) * Rho
    trand = (1 - td[iv]) * Tau + td[iv]
    fiv = refld - trand * trand / refld
    eiv = trand / refld;
    
    m = m + 1
    a[m] = -eiv
    b[m] = 1
    c[m] = -fiv
    d[m] = PARdir * tbcum[iv] * (1 - tb[iv]) * (Rho - Tau * eiv)
    
    # Downward flux
    
    refld = (1 - td[iv+1]) * Rho
    trand = (1 - td[iv+1]) * Tau + td[iv+1]
    aiv = refld - trand * trand / refld;
    biv = trand / refld;
    
    m = m + 1
    a[m] = -aiv
    b[m] = 1
    c[m] = -biv
    d[m] = PARdir * tbcum[iv+1] * (1 - tb[iv+1]) * (Tau - Rho * biv)
    
  }
  
  # Top canopy layer: upward flux
  
  iv = nlayers+1
  refld = (1 - td[iv]) * Rho
  trand = (1 - td[iv]) * Tau + td[iv]
  fiv = refld - trand * trand / refld
  eiv = trand / refld
  
  m = m + 1
  a[m] = -eiv
  b[m] = 1
  c[m] = -fiv
  d[m] = PARdir * tbcum[iv] * (1 - tb[iv]) * (Rho - Tau * eiv)
  
  # Top canopy layer: downward flux
  
  m = m + 1
  a[m] = 0
  b[m] = 1
  c[m] = 0
  d[m] = PARdif
  
  # Solve tridiagonal equations for fluxes
  
  u = f.tridiagonal.solver(a, b, c, d, m)
  
 # Now copy the solution (u) to the upward (swup) and downward (swdn) fluxes for each layer
  # swup - Upward diffuse solar flux above layer
  # swdn - Downward diffuse solar flux onto layer
  
  # Soil fluxes
  
  iv = 1
  m = 1
  swup[iv] = u[m]
  m = m + 1
  swdn[iv] = u[m]
  
  # Leaf layer fluxes
  
  for (iv in 2:(nlayers+1)){
    m = m + 1
    swup[iv] = u[m]
    m = m + 1
    swdn[iv] = u[m] 
  }

 # --- Compute flux densities
  
  # Absorbed direct beam and diffuse for ground (soil)
  
  iv = 1
  direct = PARdir * tbcum[iv] * (1 - Rho_soil_dir)
  diffuse = swdn[iv] * (1 - Rho_soil_dif)
  swsoi = direct + diffuse
  
  # Absorbed direct beam and diffuse for each leaf layer and sum
  # for all leaf layers
  
  swveg = 0
  swvegsun = 0
  swvegsha = 0
  swleafsun=swleafsha=rep(NA,nlayers+1)
  
  for (iv in 2:(nlayers+1)){
    # Per unit ground area (W/m2 ground)
    
    direct = PARdir * tbcum[iv] * (1 - tb[iv]) * (1 - omega)
    diffuse = (swdn[iv] + swup[iv-1]) * (1 - td[iv]) * (1 - omega)
    
    # Absorbed solar radiation for shaded and sunlit portions of leaf layer
    # per unit ground area (W/m2 ground)
    
    sun = diffuse * fracsun[iv] + direct
    shade = diffuse * fracsha[iv]
    
    # Convert to per unit sunlit and shaded leaf area (W/m2 leaf)
    
    swleafsun[iv] = sun / (fracsun[iv] * dLAI[iv])
    swleafsha[iv] = shade / (fracsha[iv] * dLAI[iv])
    
    # Sum fluxes over all leaf layers
    
    swveg = swveg + (direct + diffuse)
    swvegsun = swvegsun + sun
    swvegsha = swvegsha + shade
  }
  
  # --- Albedo
  
  incoming = PARdir + PARdif
  reflected = swup[nlayers+1]
  if (incoming > 0){
    albcan = reflected / incoming
  }  else {
    albcan = 0;
  }

  
  # --- Conservation check
  
  # Total radiation balance: absorbed = incoming - outgoing
  
  suminc = PARdir + PARdif
  sumref = albcan * (PARdir + PARdif)
  sumabs = suminc - sumref
  
  err = sumabs - (swveg + swsoi)
  if (abs(err) > 1e-03){
    print('err = %15.5f\n',err)
    error ('NormanRadiation: Total solar conservation error')
  }

    # Sunlit and shaded absorption
  
  err = (swvegsun + swvegsha) - swveg
  if (abs(err) > 1e-03){
    print('err = %15.5f\n',err)
    error ('NormanRadiation: Sunlit/shade solar conservation error')
    end
  }
  return(list(PARsun=rev(swleafsun[2:(nlayers+1)]),PARsha=rev(swleafsha[2:(nlayers+1)]),fracsha=rev(fracsha[2:(nlayers+1)]),fracsun=rev(fracsun[2:(nlayers+1)])))
}


#' @title Function to describe the light levels inside the canopy
#' @param meteo_hourly Hourly weather data frame with at least the column time (time in numeric, for example 0 to 23),Tair (air temperature in degree C), RH (humidity in pc from 0 to 100), cs the CO2 concentration and PFD the total PFD in micro mol m-2 s-1.
#' @param lat Latitude of the canopy to model.
#' @param DOY Day of Year .
#' @param nlayers Number of layers inside the canopy (max = 50).
#' @param dLAI LAI of each one of the n layers of vegetation in the canopy.
#' @param LAI Cumulated LAI in the midle of each layer.
#' @param Rho Leaf reflectance in the visible wavelengths.
#' @param Tau Leaf transmittance in the visible wavelengths.
#' @param Rho_soil Soil reflectance in the visible wavelengths.
#' @param Rho_NIR Leaf reflectance in the NIR wavelengths.
#' @param Tau_NIR Leaf transmittance in the NIR wavelengths.
#' @param Rho_soil_NIR Soil reflectance in the NIR wavelengths.
#' @param chil Index of departure of the leaf angles from a spherical distribution. -0.4 < chil < 0.6.
#' @param clumpfac Clumping factor, index of non random spatial distribution of leaves. = 1 for randomly spaced leaves, <1 for clumed leaves (Chen et al. 2012).

#' @param model Model for the radiation interception model, default is Norman (only Norman implemented so far).
#'
#' @return
#' @export
#'
#' @examples
#' ##Simulation of the vegetation
#' LAItot = 6
#' nlayers=20
#' dLAI=rep(6/nlayers,nlayers)
#' LAI=cumsum(dLAI)-dLAI/2 # LAI in the midle of each layer
#'##Simulation of the weather
#' meteo_hourly=data.frame(time=0:23,RH=80,Tair=25,cs=400,PFD=dnorm(x = seq(0,23,1),mean = 12,sd = 2.5)/0.16*2000,Tleaf=25)
#' ##Simulation of position and moment of the simulation
#' lat=9.2801048
#' DOY = 60
#' ##Representation of the light interception inside the canopy
#' canopy=f.canopy.interception(meteo_hourly=meteo_hourly,lat = lat,DOY = DOY,nlayers = nlayers,dLAI = dLAI,LAI=LAI)
f.canopy.interception=function(meteo_hourly,lat,DOY,nlayers,dLAI,LAI,Rho=0.11,Tau=0.06,Rho_soil=0.1,Rho_NIR=0.46,Tau_NIR=0.33,Rho_soil_NIR=0.33,chil=0.32,clumpfac=0.85,model='Norman'){
  if(is.null(meteo_hourly$NIR)){meteo_hourly$NIR=meteo_hourly$PFD/4.57}
  time=meteo_hourly$time
  ##Calculation of cosz according to Miguel et al. 2009
  phi=lat*pi/180
  delta=-23.5*cos(360*(DOY+10)/365*pi/180)*pi/180
  h=15*(time-12)*pi/180
  cosz=pmax(0.01,sin(phi)*sin(delta)+cos(phi)*cos(delta)*cos(h))
 
  ## Calculation of the proportion of direct and diffuse light using the empirical equation from CLM5 (33.7)
  a0=0.17639
  a1=0.00380
  a2=-9.0039*10^-6
  a3=8.1351*10^-9
  PFD_W=meteo_hourly$PFD/4.57 # Conversion from micro mol m-2 s-1 to W m-2 to be consistent with CLM5 eqn
  prop_dir=pmax(0.01,pmin(0.99,a0+a1*PFD_W+a2*(PFD_W)^2+a3*(PFD_W)^3))
  PFD_dir=prop_dir*meteo_hourly$PFD
  PFD_dif=meteo_hourly$PFD-PFD_dir
  
  b0=0.29548
  b1=0.00504
  b2=-1.4957*10-5
  b3=1.4881*10-8
  
  prop_NIR_dir=pmax(0.01,pmin(0.99,b0+b1*meteo_hourly$NIR+b2*(meteo_hourly$NIR)^2+b3*(meteo_hourly$NIR)^3))
  NIR_dir=prop_NIR_dir*meteo_hourly$NIR
  NIR_dif=meteo_hourly$NIR-NIR_dir
   
  #########################
  plot(x=time,y=meteo_hourly$PFD,type="l",xlab="Time of the day",ylab=expression(Light~intensity~(mu~mol~m^-2~s^-1)),col='black')
  lines(x=time,y=PFD_dif,col="blue")
  lines(x=time,y=PFD_dir,col='red')
  legend('topleft',c('Total','Diffuse light','Direct light'),col=c('black','blue','red'),lty=c(1,1,1))
  ### Creation of matrices with 50 vertical layers and 24 hours
  Canopy_time_dir=Canopy_time_dif=Canopy_time_NIR_dir=Canopy_time_NIR_dif=Canopy_time_tot=Photosynthesis_rate_dir=Photosynthesis_rate_dif=f_sun=f_shade=Temp_leaf_dir=Temp_leaf_dif=gs_dir=gs_dif=matrix(data = NA,nrow = nlayers,ncol = length(time),dimnames = list(Layer=1:nlayers,time=time))
  for(i in 1:length(time)){
    Light_Profile=f.Norman.Radiation(PARdir = PFD_dir[i],PARdif = PFD_dif[i], dLAI = dLAI,nlayers = nlayers,cosz = cosz[i],chil=chil,clumpfac = clumpfac,Rho = Rho,Tau = Tau,Rho_soil_dif =Rho_soil,Rho_soil_dir=Rho_soil)
    Light_Profile_NIR=f.Norman.Radiation(PARdir = NIR_dir[i],PARdif = NIR_dif[i], dLAI = dLAI,nlayers = nlayers,cosz = cosz[i],chil=chil,clumpfac = clumpfac,Rho =Rho_NIR,Tau=Tau_NIR,Rho_soil_dif =Rho_soil_NIR,Rho_soil_dir=Rho_soil_NIR)
    Canopy_time_dir[,i]=(Light_Profile$PARsun)
    Canopy_time_dif[,i]=(Light_Profile$PARsha)
    Canopy_time_NIR_dir[,i]=(Light_Profile_NIR$PARsun)
    Canopy_time_NIR_dif[,i]=(Light_Profile_NIR$PARsha)
    f_sun[,i]=(Light_Profile$fracsun)
    f_shade[,i]=(Light_Profile$fracsha)
  }
  
  ## VIS
  Light=Canopy_time_dir*f_sun+Canopy_time_dif*(1-f_sun)
  figure_light_dir=melt(Canopy_time_dir)
  print(ggplot(data=figure_light_dir,aes(x=time,y=Layer,fill=value))
       +scale_y_reverse()
       +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
       +labs(fill=expression(Direct~light~(mu~mol~m^-2~s^-1))))
  
  figure_light_dif=melt(Canopy_time_dif)
  print(ggplot(data=figure_light_dif,aes(x=time,y=Layer,fill=value))
        +scale_y_reverse()
        +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
        +labs(fill=expression(Diffuse~light~(mu~mol~m^-2~s^-1))))
  
  figure_f_sun=melt(f_sun)
  print(ggplot(data=figure_f_sun,aes(x=time,y=Layer,fill=value))+geom_raster()
   +scale_fill_distiller(palette = "Spectral", direction = -1) +scale_y_reverse()
    +ggtitle('Fraction of leaves in direct light')
  +labs(fill='f_sun %'))
  
  figure_light_tot=melt(Light)
  print(ggplot(data=figure_light_tot,aes(x=time,y=Layer,fill=value))
        +ggtitle('Mean PFD of an average leaf')
        +scale_y_reverse()
        +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
        +labs(fill=expression(PFD~(mu~mol~m^-2~s^-1))))
  
  ## NIR
  NIR=Canopy_time_NIR_dir*f_sun+Canopy_time_NIR_dif*(1-f_sun)
  figure_NIR_dir=melt(Canopy_time_NIR_dir)
  print(ggplot(data=figure_NIR_dir,aes(x=time,y=Layer,fill=value))
        +scale_y_reverse()
        +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
        +labs(fill=expression(Direct~NIR~(W~m^-2))))
  
  figure_NIR_dif=melt(Canopy_time_NIR_dif)
  print(ggplot(data=figure_NIR_dif,aes(x=time,y=Layer,fill=value))
        +scale_y_reverse()
        +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
        +labs(fill=expression(Diffuse~NIR~(W~m^-2))))
  
  figure_NIR_tot=melt(NIR)
  print(ggplot(data=figure_NIR_tot,aes(x=time,y=Layer,fill=value))
        +ggtitle('Mean NIR of an average leaf')
        +scale_y_reverse()
        +geom_raster()+scale_fill_distiller(palette = "Spectral", direction = -1)
        +labs(fill=expression(NIR~(W~m^-2))))
  
  return(list(Canopy_time_dir=Canopy_time_dir,Canopy_time_tot=Canopy_time_tot,Canopy_time_dif=Canopy_time_dif,Canopy_time_NIR_dir=Canopy_time_NIR_dir,Canopy_time_NIR_dif=Canopy_time_NIR_dif,f_sun=f_sun,f_shade=f_shade,Light_Profile=Light_Profile,LAItot=sum(dLAI),LAI=LAI))
}



#' @title Gradients of photosynthetic parameters
#' @description Several versions of gradients can be found in the litterature, see for example Lloyd et al. 2010 (Fig. 10 and equation A2), but also the equation A14 from Krinner et al. 2005 and the equation 33 from Clark et al. 2011
#' The simpler model describing the gradients is Vcmax(LAI)=Vcmax0 x exp(-kn x LAI) with Vcmax0 Vcmax at the top of the canopy
#' kn can be also calculated as a function of Vcmax0: kn=exp(alpha x Vcmax0+beta).
#' If kn is NULL, then the function will use the default alpha and beta to calculate kn. If, on the contrary, kn is given, this specific one will be used to calculate the gradients.
#' Krinner et al use a slightly different version of this equation with the parameter lambda: Vcmax(LAI)=Vcmax0 x (1-lambda x (1-exp(-kn*LAI))). The previous equation is a particular case of this one for lambda = 1.
#' @param alpha Slope of the relationship between Vcmax0 and log(kn), see Lloyd et al. 2010.
#' @param beta Intercept of the relationship between Vcmax0 and log(kn), see Lloyd et al. 2010.
#' @param Vcmax0 Vcmax at 25 degree C at the top of the canopy.
#' @param LAI Vector of Leaf Area Index (or depth within the canopy see Clark et al. 2011).
#' @param kn Exponential decrease.
#' @param lambda Asymptot of the decrease (see Krinner et al. 2005).
#' @references Krinner, G., Viovy, N., de Noblet-Ducoudr?, N., Og?e, J., Polcher, J., Friedlingstein, P., . Prentice, I. C. (2005). A dynamic global vegetation model for studies of the coupled atmosphere-biosphere system. Global Biogeochemical Cycles, 19(1). doi:10.1029/2003gb002199.
#' Clark, D. B., Mercado, L. M., Sitch, S., Jones, C. D., Gedney, N., Best, M. J., . Cox, P. M. (2011). The Joint UK Land Environment Simulator (JULES), model description - Part 2: Carbon fluxes and vegetation dynamics. Geoscientific Model Development, 4(3), 701-722. doi:10.5194/gmd-4-701-2011.
#' Lloyd, J., Pati?o, S., Paiva, R. Q., Nardoto, G. B., Quesada, C. A., Santos, A. J. B., . Mercado, L. M. (2010). Optimisation of photosynthetic carbon gain and within-canopy gradients of associated foliar traits for Amazon forest trees. Biogeosciences, 7(6), 1833-1859. doi:10.5194/bg-7-1833-2010.
#' @return Vector of Vcmax (or any other parameter) at the different LAI specified in the call of the function
#' @export
#'
#' @examples
#' LAI=seq(0,6.2,6.2/49)
#' Vcmax=f.VcmaxRef.LAI(kn=0.11,LAI=LAI,Vcmax0=70)
#' Vcmax2=f.VcmaxRef.LAI(kn=0.11,LAI=LAI,Vcmax0=70,lambda=0.7)
#' plot(Vcmax)
#' lines(Vcmax2)
f.VcmaxRef.LAI=function(alpha=0.00963,beta=-2.43,Vcmax0=50,LAI=0:8,kn=NULL,lambda=1){
  if(is.null(kn)){kn=exp(alpha*Vcmax0+beta)
  print(paste('kn is',kn))}
  return(Vcmax0*(1-lambda*(1-exp(-kn*LAI))))
}



#' @title Canopy scale GPP calculation
#' @description Generic function to calculate the GPP within a forest (Here GPP = sum of Anet at the canopy level, so it takes into account the leaf mitochondrial respiration)
#' @param meteo_hourly See f.canopy.interception doc. In addition to the requirement for f.canopy.interception, the leaf temperature has to be informed within the column Tleaf.
#' @param Vcmax_Profile Vector of the values of Vcmax at the reference temperature at each layer of the canopy.
#' @param Jmax_Profile Vector of the values of Jmax at the reference temperature at each layer of the canopy.
#' @param Rd_Profile Vector of the values of Rd at the reference temperature at each layer of the canopy.
#' @param Tp_Profile Vector of the values of Tp at the reference temperature at each layer of the canopy.
#' @param g0_Profile Vector of the values of g0 at the reference temperature at each layer of the canopy.
#' @param g1_Profile Vector of the values of g1 at the reference temperature at each layer of the canopy.
#' @param gsmin Minimum stomatal conductance for water to consider. This value will be used as the minimum conductance value to avoid 0 and negative values obtained from the coupled assimilation and conductance models.
#' @param canopy Description of the canopy interception (see canopy_interception function).
#' @param Patm Atmospheric pressure (used to calculate the transpiration).
#' @param ... Other parameters of the photosynthetic model, without gradients, for example curvature factor, quantum yield.. see the help of f.make.param().
#'
#' @return
#' @export
#'
#' @examples
#' #See vignettes on github
f.GPP<-function(meteo_hourly,Vcmax_Profile,Jmax_Profile,Rd_Profile,Tp_Profile,g0_Profile,g1_Profile,gsmin,canopy,Patm=100,...){
  if(length(Vcmax_Profile)!=nrow(canopy$Canopy_time_dir)){print(paste('Are you sure you want to use',length(Vcmax_Profile),'different Vcmax but ',nrow(canopy$Canopy_time_dir),'vertical canopy layers ?'))}
  VpdL_dir=VpdL_dif=Photosynthesis_rate_dir=Photosynthesis_rate_dif=gs_dir=gs_dif=canopy$Canopy_time_dir
  nlayer=nrow(canopy$Canopy_time_dir)
  g1_min=-1 ## This trick is used to fix gsw to gswmin. 
  for(Layer in 1:nlayer){
    res_dir=f.A(PFD = canopy$Canopy_time_dir[Layer,],
                cs = meteo_hourly[,"cs"],
                Tair = meteo_hourly[,"Tair"]+273.15,
                Tleaf= meteo_hourly[,"Tleaf"]+273.15,
                RH = meteo_hourly[,"RH"],
                param = f.make.param(VcmaxRef =Vcmax_Profile[Layer],
                                     RdRef = Rd_Profile[Layer],
                                     JmaxRef=Jmax_Profile[Layer],
                                     TpRef=Tp_Profile[Layer],
                                     g0=g0_Profile[Layer],
                                     g1=g1_Profile[Layer],
                                     abso=1,...
                ))
    ls.gs=which(res_dir$gs<gsmin)
    res_dir$gs[ls.gs]=gsmin
    
    res_dir$A[ls.gs]=f.A(PFD = canopy$Canopy_time_dir[Layer,],cs = meteo_hourly[,"cs"],Tleaf = meteo_hourly[,"Tleaf"]+273.15,Tair = meteo_hourly[,"Tair"]+273.15,RH = meteo_hourly[,"RH"],param = f.make.param(
                                                                                                                                                             VcmaxRef =Vcmax_Profile[Layer],
                                                                                                                                                             RdRef = Rd_Profile[Layer],
                                                                                                                                                             JmaxRef=Jmax_Profile[Layer],
                                                                                                                                                             TpRef=Tp_Profile[Layer],
                                                                                                                                                             g0=gsmin,
                                                                                                                                                             g1=g1_min,
                                                                                                                                                             abso=1,...
    ))$A[ls.gs]
    Photosynthesis_rate_dir[Layer,]=res_dir$A
    VpdL_dir[Layer,]=res_dir$ds/1000
    gs_dir[Layer,]=res_dir$gs
    res_dif=f.A(PFD = canopy$Canopy_time_dif[Layer,],
                cs = meteo_hourly[,"cs"],
                Tair = meteo_hourly[,"Tair"]+273.15,
                Tleaf= meteo_hourly[,"Tleaf"]+273.15,
                RH = meteo_hourly[,"RH"],
                param = f.make.param(VcmaxRef =Vcmax_Profile[Layer],
                                     RdRef = Rd_Profile[Layer],
                                     JmaxRef=Jmax_Profile[Layer],
                                     TpRef=Tp_Profile[Layer],
                                     g0=g0_Profile[Layer],
                                     g1=g1_Profile[Layer],
                                     abso=1,...
                ))
    ls.gs=which(res_dif$gs<gsmin)
    res_dif$gs[ls.gs]=gsmin
    res_dif$A[ls.gs]=f.A(PFD = canopy$Canopy_time_dif[Layer,],cs =meteo_hourly[,"cs"],Tleaf = meteo_hourly[,"Tleaf"]+273.15,Tair = meteo_hourly[,"Tair"]+273.15,RH = meteo_hourly[,"RH"],param = f.make.param(
                                                                                                                                                             VcmaxRef =Vcmax_Profile[Layer],
                                                                                                                                                             RdRef = Rd_Profile[Layer],
                                                                                                                                                             JmaxRef=Jmax_Profile[Layer],
                                                                                                                                                             TpRef=Tp_Profile[Layer],
                                                                                                                                                             g0=gsmin,
                                                                                                                                                             g1=g1_min,
                                                                                                                                                             abso=1,...
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
     +ylab(paste("Vertical level (0= top,",nlayer," = ground)"))
     +labs(fill=expression(A~(mu~mol~m^-2~s^-1))))
  
  Conductance_rate=(gs_dir*canopy$f_sun+gs_dif*(1-canopy$f_sun))
  Trans=(gs_dir*VpdL_dir/Patm*canopy$f_sun+gs_dif*VpdL_dif/Patm*(1-canopy$f_sun))
  figure_conductance=melt(Conductance_rate)
  b=(ggplot(data=figure_conductance,aes(x=time,y=Layer,fill=value))+geom_raster()
     +scale_fill_distiller(palette = "Spectral", direction = -1) +scale_y_reverse()
     +xlab("Time in the day")
     +ylab(paste("Vertical level (0= top,",nlayer," = ground)"))
     +labs(fill=expression(g[sw]~(mol~m^-2~s^-1))))
  print(a)
  print(b)
  totalGPP= sum(Photosynthesis_rate*dLAI,na.rm=TRUE)*365*3600*44/10^6
  totalET= sum(Trans*dLAI,na.rm=TRUE)*365*3600*18*10^-3
  print(paste("GPP = ",totalGPP,"g CO2 m-2 Ground Y-1"))
  print(paste("ET = ",totalET,"L H20 m-2 Ground Y-1"))
  return(list(A=Photosynthesis_rate,gs=Conductance_rate,A_dir=Photosynthesis_rate_dir,gs_dir=gs_dir,A_dif=Photosynthesis_rate_dif,gs_dif=gs_dif,Trans=Trans,GPP=totalGPP,ET=totalET,fig_A=a,fig_gs=b))
}

#' @title Canopy scale GPP calculation, with leaf energy budget
#' @description Generic function to calculate the GPP within a forest (Here GPP = sum of Anet at the canopy level, so it takes into account the leaf mitochondrial respiration)
#' @param meteo_hourly See f.canopy.interception
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
#' #See github vignettes
f.GPPT<-function(meteo_hourly,Vcmax_Profile,Jmax_Profile,Rd_Profile,Tp_Profile,g0_Profile,g1_Profile,gsmin,canopy,Patm=100,...){
  if(length(Vcmax_Profile)!=nrow(canopy$Canopy_time_dir)){print(paste('Are you sure you want to use',length(Vcmax_Profile),'different Vcmax but ',nrow(canopy$Canopy_time_dir),'vertical canopy layers ?'))}
  VpdL_dir=VpdL_dif=Photosynthesis_rate_dir=Photosynthesis_rate_dif=gs_dir=gs_dif=rd_dir=rd_dif=Tleaf_dir=Tleaf_dif=RHs_dir=RHs_dif=cs_dir=cs_dif=canopy$Canopy_time_dir
  nlayer=nrow(canopy$Canopy_time_dir)
  param=f.make.param()
 g1_min=-1  #This trick is used to fix gsw to gswmin
  for(Layer in 1:nlayer){
    print(paste('Layer',Layer,'of', nrow(canopy$Canopy_time_dir),'layers'))
    res_dir=f.AT(PFD = canopy$Canopy_time_dir[Layer,],
                 NIR= canopy$Canopy_time_NIR_dir[Layer,],
                 ca = meteo_hourly[,"cs"],
                 Tair = meteo_hourly[,"Tair"]+273.15,
                 wind= meteo_hourly[,'wind']*exp(-0.5*canopy$LAI[Layer]),
                 RHa = meteo_hourly[,"RH"],
                 abso_s=1,
                 param = f.make.param(VcmaxRef =Vcmax_Profile[Layer],
                                      RdRef = Rd_Profile[Layer],
                                      JmaxRef=Jmax_Profile[Layer],
                                      TpRef=Tp_Profile[Layer],
                                      g0=g0_Profile[Layer],
                                      g1=g1_Profile[Layer],abso=1,...
                 ))
    ls.gs=which(res_dir$gs<gsmin)
    res_dir$gs[ls.gs]=gsmin
    res_dir$A[ls.gs]=f.AT(PFD = canopy$Canopy_time_dir[Layer,],NIR = canopy$Canopy_time_NIR_dir[Layer,],ca = meteo_hourly[,"cs"],Tair = meteo_hourly[,"Tair"]+273.15,RHa = meteo_hourly[,"RH"],wind=meteo_hourly[,'wind']*exp(-0.5*canopy$LAI[Layer]),abso_s=1,param = f.make.param(
                                                                                                                                                      VcmaxRef =Vcmax_Profile[Layer],
                                                                                                                                                      RdRef = Rd_Profile[Layer],
                                                                                                                                                      JmaxRef=Jmax_Profile[Layer],
                                                                                                                                                      TpRef=Tp_Profile[Layer],
                                                                                                                                                      g0=gsmin,
                                                                                                                                                      g1=g1_min,abso=1,...
    ))$A[ls.gs]
    #((-g0_Profile[Layer])*400*sqrt(f.ds(Tleaf = meteo_hourly[,"tl"]+273.15,Tair = meteo_hourly[,"at"]+273.15,RH = meteo_hourly[,"RH"])/1000)/(1.6*g1_Profile[Layer]))[ls.gs]
    Photosynthesis_rate_dir[Layer,]=res_dir$A
    VpdL_dir[Layer,]=res_dir$ds/1000
    gs_dir[Layer,]=res_dir$gs
    rd_dir[Layer,]=res_dir$Rd
    Tleaf_dir[Layer,]=res_dir$Tleaf
    RHs_dir[Layer,]=res_dir$RHs
    cs_dir[Layer,]=res_dir$cs
    res_dif=f.AT(PFD = canopy$Canopy_time_dif[Layer,],
                 NIR= canopy$Canopy_time_NIR_dif[Layer,],
                 ca = meteo_hourly[,"cs"],
                 Tair = meteo_hourly[,"Tair"]+273.15,
                 wind=meteo_hourly[,'wind']*exp(-0.5*canopy$LAI[Layer]),
                 RHa = meteo_hourly[,"RH"],
                 abso_s=1,
                 param = f.make.param(VcmaxRef =Vcmax_Profile[Layer],
                                      RdRef = Rd_Profile[Layer],
                                      JmaxRef=Jmax_Profile[Layer],
                                      TpRef=Tp_Profile[Layer],
                                      g0=g0_Profile[Layer],
                                      g1=g1_Profile[Layer],abso=1,...
                 ))
    ls.gs=which(res_dif$gs<gsmin)
    res_dif$gs[ls.gs]=gsmin
    res_dif$A[ls.gs]=f.AT(PFD = canopy$Canopy_time_dif[Layer,],NIR = canopy$Canopy_time_NIR_dif[Layer,],ca = meteo_hourly[,"cs"],Tair = meteo_hourly[,"Tair"]+273.15,RHa = meteo_hourly[,"RH"],wind=meteo_hourly[,'wind']*exp(-0.5*canopy$LAI[Layer]),abso_s=1,param = f.make.param(
                                                                                                                                                      VcmaxRef =Vcmax_Profile[Layer],
                                                                                                                                                      RdRef = Rd_Profile[Layer],
                                                                                                                                                      JmaxRef=Jmax_Profile[Layer],
                                                                                                                                                      TpRef=Tp_Profile[Layer],
                                                                                                                                                      g0=gsmin,
                                                                                                                                                      g1=g1_min,abso=1,...
    ))$A[ls.gs]
    Photosynthesis_rate_dif[Layer,]=res_dif$A
    gs_dif[Layer,]=res_dif$gs
    rd_dif[Layer,]=res_dif$Rd
    VpdL_dif[Layer,]=res_dif$ds/1000
    Tleaf_dif[Layer,]=res_dif$Tleaf
    RHs_dif[Layer,]=res_dif$RHs
    cs_dif[Layer,]=res_dif$cs
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
     +ylab(paste("Vertical level (0= top,",nlayer," = ground)"))
     +labs(fill=expression(g[sw]~(mol~m^-2~s^-1))))
  figure_Tleaf=melt(Tleaf)
  c=(ggplot(data=figure_Tleaf,aes(x=time,y=Layer,fill=value))+geom_raster()
     +scale_fill_distiller(palette = "Spectral", direction = -1) +scale_y_reverse()
     +xlab("Time in the day")
     +ylab(paste("Vertical level (0= top,",nlayer," = ground)"))
     +labs(fill=expression(Tleaf~(K))))
  print(a)
  print(b)
  print(c)
  totalGPP= sum(Photosynthesis_rate*dLAI,na.rm=TRUE)*365*3600*44/10^6
  totalET= sum(Trans*dLAI,na.rm=TRUE)*365*3600*18*10^-3
  print(paste("GPP = ",totalGPP,"g CO2 m-2 Ground Y-1"))
  print(paste("ET = ",totalET,"L H20 m-2 Ground Y-1"))
  return(list(A=Photosynthesis_rate,gs=Conductance_rate,A_dir=Photosynthesis_rate_dir,gs_dir=gs_dir,A_dif=Photosynthesis_rate_dif,gs_dif=gs_dif,Tleaf_dir=Tleaf_dir,Tleaf_dif=Tleaf_dif,Tleaf=Tleaf,Rd_dir=rd_dir,Rd_dif=rd_dif,Trans=Trans,GPP=totalGPP,ET=totalET,VpdL_dif=VpdL_dif,VpdL_dir=VpdL_dir,RHs_dif=RHs_dif,RHs_dir=RHs_dir,cs_dif=cs_dif,cs_dir=cs_dir,fig_A=a,fig_gs=b,fig_Tleaf=c))
}


