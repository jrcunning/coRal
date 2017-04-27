#' Get default model parameters.
#'
#' \code{def_pars} takes no arguments and returns the default values for all
#' parameters used in the model.
#'
#' @return A list of 26 length-one numerics containing the default values for
#'   all parameters used in the model. Each list element is named as the
#'   parameter name.
#' @examples
#' def_pars()
#' pars <- def_pars()

def_pars <- function(nsym=1) {
  return(list(
    jHT0=0.03,  # Host specific biomass turnover rate (d^-1)
    nNH=0.18,  # N:C ratio in host biomass (-)
    nNX=0.2,  # N:C ratio in prey biomass (-)
    sigmaNH=0.9,  # Proportion of host nitrogen turnover recycled (-)
    sigmaCH=0.1,  # Proportion of host carbon turnover recycled (-)
    jXm=0.13,  # Maximum specific host feeding rate (molX/CmolH/d)
    jNm=0.035,  # Maximum specific host DIN uptake rate (molN/CmolH/d)
    jHGm=1,  # Maximum specific host growth rate (CmolH/CmolH/d)
    kCO2=10,  # Rate of host CCM's (molCO2/molC/d)
    KN=1.5e-6,  # Half-saturation constant for host DIN uptake (molN/L)
    KX=1e-6,  # Half-saturation constant for host feeding (CmolX/L)
    initH=1,  # Initial host biomass (CmolH)
    yC=0.8,
    jST0=rep(0.03, nsym),  # Symbiont specific biomass turnover rate (d^-1)
    nNS=rep(0.13, nsym),  # N:C ratio in symbiont biomass (-)
    yCL=rep(0.1, nsym),  # L:C ratio in fixed carbon (=quantum yield) (molC/mol ph)
    kNPQ=rep(112, nsym),  # capacity of non-photochemical quenching (mol ph/CmolS/d)
    # calculated as 4x max. photochemical quenching (Gorbunov et al. 2001)
    kROS=rep(80, nsym),  # amount of excess light beyond NPQ capacity (e.g., jeL-jNPQ) that doubles ROS production relative to baseline (mol ph/CmolS/d)
    k=rep(1, nsym),  # exponent on ROS production (-)
    astar=rep(1.34, nsym),  # Symbiont specific cross-sectional area (m^2/C-molS)
    sigmaNS=rep(0.9, nsym),  # Proportion of symbiont nitrogen turnover recylced (-)
    sigmaCS=rep(0.9, nsym),  # Proportion of symbiont carbon turnover recycled (-)
    jCPm=rep(2.8, nsym),  # Maximum specific photosynthate production rate (Cmol/CmolS/d)
    jSGm=rep(0.25, nsym),  # Maximum specific symbiont growth rate (CmolS/CmolS/d)
    initS=rep(1, nsym),  # Initial symbiont biomass (CmolS)
    b=rep(5, nsym)  # Scaling parameter for bleaching response
  ))
}

