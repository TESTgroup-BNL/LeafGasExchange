######################################################################################################
#### The aim of this model is to describe the light conditions inside the canopy.                 ####
#### Particularly the amount of light that reaches the 'sunlit' leaves, and the amount of light   ####
#### that reaches the 'shaded' leaves.                                                            ####
#### Here we consider two sources of incident light. The direct light and the diffuse light. The  ####
#### direct  light is charachterised by the solar zenith angle whereas the diffuse light comes    ####
#### from all the directions. We consider that the light which is intercepted by a layer l and    ####
#### reflected or transmitted will travel to other layers of vegetation, soil or atmosphere.      ####
#### The multiple exchanges of light between layers are calculated iteratively.                   ####

#### Bug: the quantity of light absorbed + reflected is different than the total light coming into the system 
## for nlayers>1. I cant find why... See at the end of the code


##### Inputs of the model
cosz= 0.45            ## Solar zenith angle (>0 and <=1)
chil= 0.1           ## Angle distribution parameter of the leaves
clumping_index= 0.85  ## Vegetation dispersion parameter that quantifies the level of foliage distribution 
                      #non-randomness
nlayers=10            ## Number of vertical layers of vegetation
                      # !!!!!! Important, in the code we consider nlayers +1 and the last layers correspond to the soil
LAItot=6              ## Total LAI of the canopy
nitermax=1000         ## Number of light rebounds considered in the iterative loop
rhol=0.1              ## Reflectance of the leaves (for direct and diffuse light)
taul=0.05             ## Transmittance of the leaves (for direct and diffuse light)
rhos_dir=0.1          ## Reflectance of the soil for direct light
rhos_dif=0.1          ## Reflectance of the soil for diffuse light

PFD_dir_top=2000      ## Quantity of direct light at the top of the canopy (micro mol m-2 s-1)
PFD_dif_top=200       ## Quantity of diffuse light at the top of the canopy (micro mol m-2 s-1)


#### Construction of a LAI profile inside the canopy
LAIl= LAItot/nlayers  ## Here we assume that the LAI is homogeneous in the different layers of the canopy, 
                      #with a LAI of LAIl inside each layer

LAI_l=seq(LAIl,LAItot,LAIl) ## LAI at the bottom of each layer 
LAI_above_l=c(0,LAI_l[1:(length(LAI_l)-1)]) ## LAI at the top of each layer
if(nlayers==1){LAI_above_l=0}
LAI_all=c(0,LAI_l[1:(length(LAI_l))]) 


## Calculation of the coefficient of extinction of the light inside the canopy
sb = (90- (acos(cosz)*180/pi)) * (pi/ 180)
phi1b = 0.5 - 0.633*chil - 0.330*chil^2
phi2b = 0.877* (1 - 2*phi1b)
gdir = phi1b + phi2b * sin(sb)
k_dir = clumping_index * gdir / sin(sb)


## Numerical calculation of the proportion of diffuse light which is not intercepted by the different layers
tr_dif_l=rep(0,nlayers) ## Proportion of diffuse light not intercepted by the layer l
tr_dif_top_l=rep(0,(nlayers+1))  ## Proportion of diffuse light not intercepted at the top of each layer

n_angle=180 ## number of different angles to consider
angle=seq(180/(2*n_angle),180,180/n_angle)*pi/360
for(alpha in  angle){
  gdir_j = phi1b + phi2b * sin(alpha)
  k_dir_j=clumping_index * gdir_j / sin(alpha)
  tr_dif_l=tr_dif_l+exp(-k_dir_j*(LAI_l-LAI_above_l))
  tr_dif_top_l=tr_dif_top_l+exp(-k_dir_j*(LAI_all))
}
tr_dif_l=tr_dif_l/n_angle
tr_dif_top_l=tr_dif_top_l/n_angle

tr_dir_l = exp(-k_dir * (LAI_l-LAI_above_l)) ## Proportion of direct light crossing layer l without being intercepted

f_sun_l=exp(-k_dir*LAI_l) ## Proportion of sunlit leaves
f_sha_l=1-f_sun_l ## Proportion of shaded leaves

int_dir_l=(1-tr_dir_l) ## Proportion of the direct light which is intercepted when crossing layer l
trans_dir_l=taul*int_dir_l ## Proportion of the direct light which is transmitted when crossing layer l
refl_dir_l=rhol*int_dir_l ## Proportion of the direct light which is reflected when crossing layer l
abs_dir_l=(1-rhol-taul)*int_dir_l ## Proportion of the direct light which is absorbed when crossing layer l

tr_dif_l  ## Proportion of diffuse light crossing layer l without being intercepted
int_dif_l=(1-tr_dif_l) ## Proportion of the diffuse light which is intercepted when crossing layer l
trans_dif_l=taul*int_dif_l ## Proportion of diffuse light which is transmitted when crossing layer l
refl_dif_l=rhol*int_dif_l ## Proportion of the diffuse light which is reflected when crossing layer l
abs_dif_l=(1-rhol-taul)*int_dif_l ## Proportion of the diffuse light which is absorbed when crossing layer l

## Adding the coefficients of the soil, which is considered here as the layer nlayers +1 
tr_dir_l[(length(tr_dir_l)+1)]=0 #There is no light transmitted below the soil layer
tr_dif_l[(length(tr_dif_l)+1)]=0 #There is no light transmitted below the soil layer
refl_dir_l[(length(refl_dir_l)+1)]=rhos_dir
refl_dif_l[(length(refl_dif_l)+1)]=rhos_dif
abs_dir_l[(length(abs_dir_l)+1)]=1-rhos_dir
abs_dif_l[(length(abs_dif_l)+1)]=1-rhos_dif
trans_dir_l[(length(trans_dir_l)+1)]=0
trans_dif_l[(length(trans_dif_l)+1)]=0
## Verification, this line should give a vector of 1: refl_dir_l+abs_dir_l+trans_dir_l+tr_dir_l

top_rad_dir_l= exp(-k_dir * (LAI_all))*PFD_dir_top ## Calculation of the quantity of direct light at the top of 
                                                  #each layer 
top_rad_dif_l=tr_dif_top_l*PFD_dif_top ## Calculation of the quantity of diffuse light entering each layer 


I_up_l=matrix(data=NA,nrow = nlayers+1,ncol = nitermax) ## Matrix of the radiation going up 
I_down_l=matrix(data=NA,nrow = nlayers+1,ncol = nitermax) ## Matrix of the radiation going down 
I_abs_dif_l=matrix(data=NA,nrow = nlayers+1,ncol = nitermax) ## Matrix of the absorbed radiation

## Initialization of the matrix. We consider here the radiations going up, down and absorbed after the first 
# 'throw' of light
I_up_l[,1]=refl_dir_l*top_rad_dir_l+refl_dif_l*top_rad_dif_l
I_down_l[,1]=trans_dir_l*top_rad_dir_l+trans_dif_l*top_rad_dif_l
I_abs_dir_l=abs_dir_l*top_rad_dir_l
I_abs_dif_l[,1]=abs_dif_l*top_rad_dif_l

## We consider here the fate of the light which is not already absorbed or reflected out of the system
# by the top canopy layer
for (iter in 2:nitermax){
  if(nlayers>2){
    ## The light going up for the layer l corresponds to the light going down at the previous iteration for the layer above (l-1) which is reflected
    # + the light going up from layer l+1, which is not intercepted by the layer l
    # + the light going up from layer l+1, which is intercepted by the layer l and which is transmitted.

  I_up_l[2:(nlayers),iter]=I_up_l[3:(nlayers+1),iter-1]*tr_dif_l[2:nlayers]
                          +I_up_l[3:(nlayers+1),iter-1]*trans_dif_l[2:nlayers]
                          +refl_dif_l[2:nlayers]*I_down_l[1:(nlayers-1),iter-1]
    # Particular cases of the soil layer and top layer
  I_up_l[(nlayers+1),iter]=refl_dif_l[(nlayers+1)]*I_down_l[nlayers,iter-1]
  I_up_l[1,iter]=I_up_l[2,iter-1]*tr_dif_l[1]+I_up_l[2,iter-1]*trans_dif_l[1]
  
  
  ## The light absorbed by the layer l corresponds to:
  I_abs_dif_l[2:nlayers,iter]=abs_dif_l[2:nlayers]*I_up_l[3:(nlayers+1),iter-1]
                          +abs_dif_l[2:nlayers]*I_down_l[1:(nlayers-1),iter-1]
  
  I_abs_dif_l[nlayers+1,iter]=abs_dif_l[nlayers+1]*I_down_l[(nlayers),iter-1]
  I_abs_dif_l[1,iter]=abs_dif_l[1]*I_up_l[2,iter-1]
  
  ## The light going down for the layer l corresponds to the light going down at the previous iteration 
  # for the layer above (l-1) which is not intercepted or which is intercepted and transmitted
  # + the light going up from layer l+1, which is not reflected by the layer l
  
  
  I_down_l[2:(nlayers),iter]=I_down_l[1:(nlayers-1),iter-1]*tr_dif_l[2:nlayers]
                            +I_down_l[1:(nlayers-1),iter-1]*trans_dif_l[2:nlayers]
                            +refl_dif_l[2:(nlayers)]*I_up_l[3:(nlayers+1),iter-1]
  I_down_l[(nlayers+1),iter]=0
  I_down_l[1,iter]=refl_dif_l[1]*I_up_l[2,iter-1]
  }
  
  if(nlayers==1){
    I_up_l[1,iter]=I_up_l[2,iter-1]*tr_dif_l[1]+I_up_l[2,iter-1]*trans_dif_l[1]
    I_up_l[2,iter]=I_down_l[1,iter-1]*refl_dif_l[2]
    
    I_abs_dif_l[1,iter]=I_up_l[2,iter-1]*abs_dif_l[1]
    I_abs_dif_l[2,iter]=I_down_l[1,iter-1]*abs_dif_l[2]
    
    I_down_l[1,iter]=I_up_l[2,iter-1]*refl_dif_l[1]
    I_down_l[2,iter]=0
  }
  
  if(nlayers==2){
    I_up_l[1,iter]=I_up_l[2,iter-1]*tr_dif_l[1]+I_up_l[2,iter-1]*trans_dif_l[1]
    I_up_l[2,iter]=I_down_l[1,iter-1]*refl_dif_l[2]+I_up_l[3,iter-1]*tr_dif_l[2]
                  +I_up_l[3,iter-1]*trans_dif_l[2]
    I_up_l[3,iter]=I_down_l[2,iter-1]*refl_dif_l[3]
    
    I_abs_dif_l[1,iter]=I_up_l[2,iter-1]*abs_dif_l[1]
    I_abs_dif_l[2,iter]=I_down_l[1,iter-1]*abs_dif_l[2]
                    +I_up_l[3,iter-1]*abs_dif_l[2]
    I_abs_dif_l[3,iter]=I_down_l[2,iter-1]*abs_dif_l[3]
    
    I_down_l[1,iter]=I_up_l[2,iter-1]*refl_dif_l[1]
    I_down_l[2,iter]=I_up_l[3,iter-1]*refl_dif_l[2]
                    +I_down_l[1,iter-1]*tr_dif_l[2]
                    +I_down_l[1,iter-1]*trans_dif_l[2]
    I_down_l[3,iter]=0
  }
}



## Verification of the model. 
## !!! I have a mistake somewhere since the quantity of light that reaches the system is different from the quantity of
## light which is absorbed by the soil or absorbed by the vegetation or reflected in the atmosphere

# The quantity of light that reaches the system is:
total_light=PFD_dif_top+PFD_dir_top

# The quantity of light that goes out the system reflected by the top canopy layer is:
refl_canopy=sum(I_up_l[1,])

# The quantity of light that goes out the system absorbed by the soil is
# the sum of the direct light which is absorbed by the soil and the diffuse light
# which is absorbed by the soil

abs_soil=sum(I_abs_dif_l[(nlayers+1),])+I_abs_dir_l[nlayers+1]

## The quantity of light absorbed by the canopy corresponds to the sum of the direct light and the diffuse light
abs_sha_l=rowSums(I_abs_dif_l[1:(nlayers),])*(f_sha_l)  
abs_sha=sum(abs_sha_l)
abs_sun_l=rowSums(I_abs_dif_l[1:(nlayers),])*(f_sun_l)+I_abs_dir_l[1:nlayers]
abs_sun=sum(abs_sun_l)
abs_veg_l=abs_sha_l+abs_sun_l
abs_veg=sum(abs_veg_l)
total_light
abs_soil+abs_veg+refl_canopy

###Calculation of the PFD received by the sunlit leaves and the PFD received by the shaded leaves
## LAI of sunlit leaves
LAI_sun_l=(LAI_l-LAI_above_l)*f_sun_l
LAI_sha_l=(LAI_l-LAI_above_l)*f_sha_l
LAI_sun=sum(LAI_sun_l)
LAI_sha=sum(LAI_sha_l)
LAI_sun+LAI_sha

PFD_abs_sun_l=abs_sun_l/(LAI_sun_l)
PFD_abs_sha_l=abs_sha_l/(LAI_sha_l)

