% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Photosynthesis_functions.R
\name{f.gsmin}
\alias{f.gsmin}
\title{Calculation of the minimal conductance given by a particular coupled conductance and photosynthesis model}
\usage{
f.gsmin(
  RdRef = 0.825,
  RdHa = 46390,
  RdHd = 150650,
  RdS = 490,
  Tleaf = 300,
  cs = 400,
  ds = 1000,
  g0 = 0.02,
  g1 = 4.1,
  power = 0.5,
  model = "USO"
)
}
\arguments{
\item{RdRef}{Respiration value at the reference temperature in micromol.m-2.s-1.}

\item{RdHa}{Energie of activation for Rd in J.mol-1.}

\item{Tleaf}{Leaf temperature in Kelvin}

\item{cs}{CO2 concentration at the leaf surface in ppm}

\item{ds}{Leaf to air vapor pressure deficit in Pa}

\item{g0}{Constant of the conductance model, representing the conductance when A is 0, in mol.m-2.s-1, usually around 0.01.}

\item{g1}{Slope parameter, between 1.14 and 3.58 KPa^0.5 (Wu et al., 2019).}

\item{power}{Power of VPDl in USO model. By default power=0.5 as in Medlyn article.}

\item{RH}{Air relative humidity (0 to 100), only necessary for the BWB stomatal conductance model}
}
\value{
Minimum conductance in mol m-2 s-1
}
\description{
The minimal conductance of a model depends on the parameters of the model (ie g0 and g1) but also on the minimum A value, which corresponds to the dark respiration.
Knowing the minimal conductance is important because the conductance can become negative and lead to unrealistic values in photosynthesis models
}
\examples{
f.gsmin(RdRef=	0.825,RdHa=	46390,RdHd=150650,RdS=490,Tleaf=300,cs=400,ds=1000,g0=0.02,g1=4.1,power=0.5,model="USO")
}
