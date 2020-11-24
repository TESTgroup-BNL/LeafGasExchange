## The quantity of light absorbed + reflected is different than the total light coming into the system...
## THe calcul of PFD_sun and PFD_sha should be verified

### Calculation of the direct beam extinction coefficient kdir
cosz= 0.8 ## Solar zenith angle (>0 and <=1)
chil= 0.1 ##Angle distribution parameter of the leaves
clumping_index= 0.85 ##
nlayers=100 ## number of vertical layers
LAItot=1  ## Total LAI of the canopy
LAIl= LAItot/nlayers ## Here we assume that the LAI is homogeneous in the different layers of the canopy
LAI=seq(LAIl,LAItot,LAIl) ## LAI at the bottom of each layer 
LAI_above=c(0,LAI[1:(length(LAI)-1)]) ## LAI at the top of each layer
if(nlayers==1){LAI_above=0}
LAI_all=c(0,LAI[1:(length(LAI))]) 
nitermax=1000 ## Number of light rebounds considered

rhol=0.1  ## Reflectance of the leaves (for direct and diffuse light)
taul=0.05 ## Transmittance of the leaves (for direct and diffuse light)
rhos_dir=0.1 ## Reflectance of the soil for direct light
rhos_dif=0.1 ## Reflectance of the soil for diffuse light

PFD_dir_top=1500 ## QUantity of direct light at the top of the canopy
PFD_dif_top=150 ## QUantity of diffuse light at the top of the canopy

## Calculation of the coefficient of extinction of the light inside the canopy
sb = (90- (acos(cosz)*180/pi)) * (pi/ 180)
phi1b = 0.5 - 0.633*chil - 0.330*chil^2
phi2b = 0.877* (1 - 2*phi1b)
gdir = phi1b + phi2b * sin(sb)
k_dir = clumping_index * gdir / sin(sb)




## Numerical calculation of the proportion of diffuse light which is not intercepted by the different layers
tr_dif_l=rep(0,nlayers) ## Proportion of diffuse light not intercepted by the layer l
tr_dif_top_l=rep(0,(nlayers+1)) 

n_angle=180 ## number of different angles to consider
angle=seq(180/(2*n_angle),180,180/n_angle)*pi/360
for(alpha in  angle){
  gdir_j = phi1b + phi2b * sin(alpha)
  k_dir_j=clumping_index * gdir_j / sin(alpha)
  tr_dif_l=tr_dif_l+exp(-k_dir_j*(LAI-LAI_above))
  tr_dif_top_l=tr_dif_top_l+exp(-k_dir_j*(LAI_all))
}
tr_dif_l=tr_dif_l/n_angle
tr_dif_top_l=tr_dif_top_l/n_angle

tr_dir_l = exp(-k_dir * (LAI-LAI_above)) ## Proportion of direct light crossing layer l without being intercepted

f_sun=exp(-k_dir*LAI) ## Proportion of sunlit leaves
f_sha=1-f_sun ## Proportion of shaded leaves

int_dir_l=(1-tr_dir_l) ## Proportion of the direct light which is intercepted when crossing layer l
trans_dir_l=taul*int_dir_l ## Proportion of the direct light which is transmitted when crossing layer l
refl_dir_l=rhol*int_dir_l ## Proportion of the direct light which is reflected when crossing layer l
abs_dir_l=(1-rhol-taul)*int_dir_l ## Proportion of the direct light which is absorbed when crossing layer l

tr_dif_l  ## Proportion of diffuse light crossing layer l without being intercepted
int_dif_l=(1-tr_dif_l) ## Proportion of the diffuse light which is intercepted when crossing layer l
trans_dif_l=taul*int_dif_l ## Proportion of diffuse light which is transmitted when crossing layer l
refl_dif_l=rhol*int_dif_l ## Proportion of the diffuse light which is reflected when crossing layer l
abs_dif_l=(1-rhol-taul)*int_dif_l ## Proportion of the diffuse light which is absorbed when crossing layer l

## Adding the coefficients of the soil, which is considered here as the nlayers +1 
tr_dir_l[(length(tr_dir_l)+1)]=0
tr_dif_l[(length(tr_dif_l)+1)]=0
refl_dir_l[(length(refl_dir_l)+1)]=rhos_dir
refl_dif_l[(length(refl_dif_l)+1)]=rhos_dif
abs_dir_l[(length(abs_dir_l)+1)]=1-rhos_dir
abs_dif_l[(length(abs_dif_l)+1)]=1-rhos_dif
trans_dir_l[(length(trans_dir_l)+1)]=0
trans_dif_l[(length(trans_dif_l)+1)]=0
refl_dir_l+abs_dir_l+trans_dir_l+tr_dir_l

top_rad_dir= exp(-k_dir * (LAI_all))*PFD_dir_top ## Calculation of the quantity of direct light entering each layer 
top_rad_dif=tr_dif_top_l*PFD_dif_top ## Calculation of the quantity of diffuse light entering each layer 
I_up_l=matrix(data=NA,nrow = nlayers+1,ncol = nitermax) ## Matrix of the radiation going up 
I_down_l=matrix(data=NA,nrow = nlayers+1,ncol = nitermax) ## Matrix of the radiation going down 
I_abs_l=matrix(data=NA,nrow = nlayers+1,ncol = nitermax) ## Matrix of the absorbed radiation

## Initialization of the matrix. We consider here the radiations going up, down and absorbed after the first rebound of light
I_up_l[,1]=refl_dir_l*top_rad_dir+refl_dif_l*top_rad_dif
I_down_l[,1]=trans_dir_l*top_rad_dir+trans_dif_l*top_rad_dif
I_abs_dir_l=abs_dir_l*top_rad_dir
I_abs_l[,1]=abs_dif_l*top_rad_dif

## We consider here the folowing rebounds
for (iter in 2:nitermax){
  if(nlayers>2){
    
    ## The light going up for the layer l corresponds to the light going down at the previous iteration for the layer above (l-1) which is reflected
    # + the light going up from layer l+1, which is not intercepted by the layer l
    # + the light going up from layer l+1, which is intercepted by the layer l and which is transmitted.
  I_up_l[(nlayers+1),iter]=refl_dif_l[(nlayers+1)]*I_down_l[nlayers,iter-1]
  I_up_l[2:(nlayers),iter]=I_up_l[3:(nlayers+1),iter-1]*tr_dif_l[2:nlayers]+I_up_l[3:(nlayers+1),iter-1]*trans_dif_l[2:nlayers]+refl_dif_l[2:nlayers]*I_down_l[1:(nlayers-1),iter-1]
  I_up_l[1,iter]=I_up_l[2,iter-1]*tr_dif_l[1]+I_up_l[2,iter-1]*trans_dif_l[1]
  
  I_abs_l[1,iter]=abs_dif_l[1]*I_up_l[2,iter-1]
  I_abs_l[2:nlayers,iter]=abs_dif_l[2:nlayers]*I_up_l[3:(nlayers+1),iter-1]+abs_dif_l[2:nlayers]*I_down_l[1:(nlayers-1),iter-1]
  I_abs_l[nlayers+1,iter]=abs_dif_l[nlayers+1]*I_down_l[(nlayers),iter-1]
  
  ## The light going down for the layer l corresponds to the light going down at the previous iteration for the layer above (l-1) which is not intercepted or which is intercepted and transmitted
  # + the light going up from layer l+1, which is not reflected by the layer l

  I_down_l[1,iter]=refl_dif_l[1]*I_up_l[2,iter-1]
  I_down_l[2:(nlayers),iter]=I_down_l[1:(nlayers-1),iter-1]*tr_dif_l[2:nlayers]+I_down_l[1:(nlayers-1),iter-1]*trans_dif_l[2:nlayers]+refl_dif_l[2:(nlayers)]*I_up_l[3:(nlayers+1),iter-1]
  I_down_l[(nlayers+1),iter]=0
  }
  
  if(nlayers==1){
    I_up_l[1,iter]=I_up_l[2,iter-1]*tr_dif_l[1]+I_up_l[2,iter-1]*trans_dif_l[1]
    I_up_l[2,iter]=I_down_l[1,iter-1]*refl_dif_l[2]
    
    I_abs_l[1,iter]=I_up_l[2,iter-1]*abs_dif_l[1]
    I_abs_l[2,iter]=I_down_l[1,iter-1]*abs_dif_l[2]
    
    I_down_l[1,iter]=I_up_l[2,iter-1]*refl_dif_l[1]
    I_down_l[2,iter]=0
  }
  
  if(nlayers==2){
    I_up_l[1,iter]=I_up_l[2,iter-1]*tr_dif_l[1]+I_up_l[2,iter-1]*trans_dif_l[1]
    I_up_l[2,iter]=I_down_l[1,iter-1]*refl_dif_l[2]+I_up_l[3,iter-1]*tr_dif_l[2]+I_up_l[3,iter-1]*trans_dif_l[2]
    I_up_l[3,iter]=I_down_l[2,iter-1]*refl_dif_l[3]
    
    I_abs_l[1,iter]=I_up_l[2,iter-1]*abs_dif_l[1]
    I_abs_l[2,iter]=I_down_l[1,iter-1]*abs_dif_l[2]+I_up_l[3,iter-1]*abs_dif_l[2]
    I_abs_l[3,iter]=I_down_l[2,iter-1]*abs_dif_l[3]
    
    I_down_l[1,iter]=I_up_l[2,iter-1]*refl_dif_l[1]
    I_down_l[2,iter]=I_up_l[3,iter-1]*refl_dif_l[2]+I_down_l[1,iter-1]*tr_dif_l[2]+I_down_l[1,iter-1]*trans_dif_l[2]
    I_down_l[3,iter]=0
  }
}


## Verification
# The quantity of light that reaches the system is:
total_light=PFD_dif_top+PFD_dir_top

# The quantity of light that goes out the system reflected by the top canopy layer is:
refl_canopy=sum(I_up_l[1,])

# The quantity of light that goes out the system absorbed by the soil is
# the sum of the direct light which is absorbed by the soil and the diffuse light
# which is absorbed by the soil

abs_soil=sum(I_abs_l[(nlayers+1),])+I_abs_dir_l[nlayers+1]

## The quantity of light absorbed by the canopy corresponds to the sum of the direct light and the diffuse light
abs_veg=sum(I_abs_l[1:(nlayers),]*(f_sha+f_sun))+sum(I_abs_dir_l[1:nlayers])


total_light
abs_soil+abs_veg+refl_canopy

## Here I consider that the PFD received by the sunlit leaves correspond to the direct radiation at the top of each layers + 
## the diffuse radiation and the diffuse radiation propagated into the canopy
sun_PFD_l=top_rad_dir[1:nlayers]+top_rad_dif[1:nlayers]+rowSums(I_up_l[2:(nlayers+1),])+c(0,rowSums(I_down_l)[1:(nlayers-1)])
sha_PFD_l=top_rad_dif[1:nlayers]+rowSums(I_up_l[2:(nlayers+1),])+c(0,rowSums(I_down_l)[1:(nlayers-1)])

print(sun_PFD_l)
print(sha_PFD_l)
