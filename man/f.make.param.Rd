% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Photosynthesis_functions.R
\name{f.make.param}
\alias{f.make.param}
\title{Photosynthesis and stomata model parameters}
\usage{
f.make.param(
  R = NA,
  O2 = NA,
  TRef = NA,
  Patm = NA,
  JmaxRef = NA,
  JmaxHa = NA,
  JmaxHd = NA,
  JmaxS = NA,
  VcmaxRef = NA,
  VcmaxHa = NA,
  VcmaxHd = NA,
  VcmaxS = NA,
  TpRef = NA,
  TpHa = NA,
  TpHd = NA,
  TpS = NA,
  thetacj = NA,
  thetaip = NA,
  RdRef = NA,
  RdHa = NA,
  RdHd = NA,
  RdS = NA,
  KcRef = NA,
  KcHa = NA,
  KoRef = NA,
  KoHa = NA,
  GstarRef = NA,
  GstarHa = NA,
  abso = NA,
  aQY = NA,
  Theta = NA,
  model.gs = NA,
  g0 = NA,
  g1 = NA,
  power = NA
)
}
\arguments{
\item{R}{Ideal gas constant.}

\item{O2}{O2 concentration in ppm.}

\item{TRef}{Reference temperature for Kc, Ko, Rd,GammaStar Vcmax, Jmax in Kelvin.}

\item{Patm}{Atmospheric pressure in kPa.}

\item{JmaxRef}{Maximum electron transport rate in micromol.m-2.s-1.}

\item{JmaxHa}{Energy of activation for Jmax in J.mol-1.}

\item{JmaxHd}{Energy of desactivation for Jmax in J.mol-1.}

\item{JmaxS}{Entropy term for Jmax in J.mol-1.K-1.}

\item{VcmaxRef}{Maximum rate of Rubisco for carboxylation micromol.m-2.s-1.}

\item{VcmaxHa}{Energy of activation for Vcmax in J.mol-1.}

\item{VcmaxHd}{Energy of desactivation for Vcmax in J.mol-1.}

\item{VcmaxS}{Entropy term for Vcmax in J.mol-1.K-1.}

\item{TpRef}{Triose phosphate utilization  rate in micromol.m-2.s-1.}

\item{TpHa}{Activation energy for Tp in J.mol-1.}

\item{TpHd}{Energy of deactivation for Tp in J.mol-1.}

\item{TpS}{Entropy term for Tp in J.mol-1.K-1.}

\item{thetacj}{Collatz smoothing factor used to introduce a gradual transition from Ac to Aj (close to 0.999)}

\item{thetaip}{Collatz smoothing factor used to introduce a gradual transition from Aj to Ap (close to 0.999)}

\item{RdRef}{Respiration value at the reference temperature in micromol.m-2.s-1.}

\item{RdHa}{Energie of activation for Rd in J.mol-1.}

\item{KcRef}{Michaelis-Menten constant of Rubisco for CO2 at the reference temperature in micromol.mol-1.}

\item{KcHa}{Energy of activation for Kc in J.mol-1.}

\item{KoRef}{Michaelis-Menten constant of Rubisco for CO2 at the reference temperature in milimol.mol-1.}

\item{KoHa}{Energy of activation for Ko in J.mol-1.}

\item{GstarRef}{CO2 compensation point in absence of respiration in micromol.mol-1.}

\item{GstarHa}{Enthalpie of activation for Gstar in J.mol-1.}

\item{abso}{Absorptance of the leaf in the photosynthetic active radiation wavelenghts.}

\item{aQY}{Apparent quantum yield.}

\item{Theta}{Theta is the empirical curvature factor for the response of J to PFD. It takes its values between 0 and 1.}

\item{model.gs}{Type of conductance model (USO, USO_simpl,BWB or Nonlinear). See f.gs documentation for more information.}

\item{g0}{Constant of the conductance model, representing the conductance when A is 0, in mol.m-2.s-1, usually around 0.01.}

\item{g1}{Slope parameter, between 1.14 and 3.58 KPa^0.5 (Wu et al., 2019).}

\item{power}{Power of VPDl in USO model. By default power=0.5 as in Medlyn article.}
}
\value{
List of parameters that can be used in other functions of the package such as f.A, f.Aci, and f.GPP
}
\description{
Function to create a list of parameters to be used in most of the functions of this package.
Depending on the function, all the parameters are not used. For example go and g1 are not used in f.Aci.
}
\details{
The call of this function is made using f.make.param(). If a parameter is modified, for example writing f.make.param(VcmaxRef=10), this function will return all the default parameters from FATES TBM with VcmaxRef = 10 instead of its default value
}
\examples{
param1=f.make.param(JmaxRef=100,VcmaxRef=60,RdRef=1,TpRef=10)
f.A(PFD=1500,cs=400,Tleaf=300,Tair=299,RH=70,param=param1)
}
\references{
Bernacchi, C.J., Singsaas, E.L., Pimentel, C., Portis Jr, A.R. and Long, S.P. (2001), Improved temperature response functions for models of Rubisco‐limited photosynthesis. Plant, Cell & Environment, 24: 253-259. doi:10.1111/j.1365-3040.2001.00668.

FATES: https://fates-docs.readthedocs.io/en/latest/fates_tech_note.html.

Medlyn, B.E., Duursma, R.A., Eamus, D., Ellsworth, D.S., Colin Prentice, I., Barton, C.V.M., Crous, K.Y., de Angelis, P., Freeman, M. and Wingate, L. (2012), Reconciling the optimal and empirical approaches to modelling stomatal conductance. Glob Change Biol, 18: 3476-3476. doi:10.1111/j.1365-2486.2012.02790.x.

Leuning, R., Kelliher, F. M., De Pury, D. G. G., & Schulze, E. D. (1995). Leaf nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to canopies. Plant, Cell & Environment, 18(10), 1183-1200.

Ball, J. T., Woodrow, I. E., & Berry, J. A. (1987). A model predicting stomatal conductance and its contribution to the control of photosynthesis under different environmental conditions. In Progress in photosynthesis research (pp. 221-224). Springer, Dordrecht.
}
