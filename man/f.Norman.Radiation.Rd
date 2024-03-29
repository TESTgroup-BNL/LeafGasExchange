% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Canopy_scaling.R
\name{f.Norman.Radiation}
\alias{f.Norman.Radiation}
\title{Norman 1979 Radiation interception model
Converted into a R code from the original code of Gordon Bonan: Bonan, G. (2019). Climate Change and Terrestrial Ecosystem Modeling. Cambridge: Cambridge University Press. doi:10.1017/9781107339217}
\usage{
f.Norman.Radiation(
  Rho = 0.1,
  Tau = 0.05,
  Rho_soil_dir = 0.1,
  Rho_soil_dif = 0.1,
  cosz,
  chil,
  clumpfac,
  dLAI,
  nlayers,
  PARdir = 0.8,
  PARdif = 0.2
)
}
\arguments{
\item{Rho}{Leaf reflectance.}

\item{Tau}{Leaf transmittance.}

\item{Rho_soil_dir}{Direct beam albedo of ground (soil).}

\item{Rho_soil_dif}{Direct beam albedo of ground (soil).}

\item{cosz}{Cosinus of the solar zenith angle.}

\item{chil}{Index of departure of the leaf angles from a spherical distribution. -0.4 < chil < 0.6.}

\item{clumpfac}{Clumping factor, index of non random spatial distribution of leaves. = 1 for randomly spaced leaves, <1 for clumed leaves (Chen et al. 2012).}

\item{dLAI}{LAI of each one of the n layers of vegetation in the canopy, layer 1 is the top of canopy, layer n is the bottom.}

\item{nlayers}{Number of vegetation layers.}

\item{PARdir}{Atmospheric direct beam solar radiation (W/m2).}

\item{PARdif}{Atmospheric diffuse solar radiation (W/m2).}
}
\value{
list of output:
PARsun Absorbed PFD by the sunlit leaves
PARsha Absorbed PFD by the shaded leaves
fracsun Proportion of sunlit leaves
fracsha Proportion of shaded leaves
}
\description{
Norman 1979 Radiation interception model
Converted into a R code from the original code of Gordon Bonan: Bonan, G. (2019). Climate Change and Terrestrial Ecosystem Modeling. Cambridge: Cambridge University Press. doi:10.1017/9781107339217
}
\examples{
f.Norman.Radiation(Rho=0.1, Tau=0.05, PARdir=1000,PARdif=200,dLAI=c(rep(6/20,20)),nlayers=20,Rho_soil_dif = 0.1,Rho_soil_dir = 0.1,cosz = 0.88,chil = 0.1,clumpfac = 0.8)
}
