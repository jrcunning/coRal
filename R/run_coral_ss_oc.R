#' Run coral-\emph{Symbiodinium} model to steady state.
#'
#' This function runs the coral-Symbiodinium bioenergetic model in a constant
#' environment until it reaches a steady state.
#'
#' @param env A list defining the constant environment in which to run model to steady state, with three named elements:
#' \describe{
#'     \item{L}{External light in mol photons m^-2 d^-1 (numeric).}
#'     \item{N}{External [DIN] in mol N L^-1 (numeric).}
#'     \item{X}{External [prey] in C-mol X L^-1 (numeric).}
#'   }
#' @param pars A list of (named) parameter values to use in the simulation.
#'   Parameter names must match those returned by \code{def_pars()}.
#' @param dt The time step (in days) to use in the numerical analysis. Values > 0.1 are not recommended.
#' @return A list with the following named elements:
#'   \describe{
#'     \item{env}{A list of the environment (L, N, X) used in the simulation}
#'     \item{pars}{A list of the parameter values used in the simulation}
#'     \item{H}{A tibble of all host biomass-specific model
#'     fluxes at each time step}
#'     \item{S}{A tibble of all symbiont biomass-specific
#'     model fluxes at each time step}
#'   }
#' @seealso \code{\link{init_env}}, \code{\link{def_pars}}
#' @examples
#' ssrun <- run_coral_ss(env=list(L=20, N=1e-7, X=0), pars=def_pars(), dt=0.1)

run_coral_ss_oc <- function(env, pars=def_pars(), dt=0.1, jOC = 0.0167) {

  # Initialize environment

  # Set initial values
  # ==================
  # Host
  H <- dplyr::data_frame(
    H=pars$initH,
    jX=(pars$jXm * env$X / (env$X + pars$KX)),
    jN=(pars$jNm * env$N / (env$N + pars$KN)),
    rhoN=jN,
    jeC=10,
    jCO2=pars$kCO2 * jeC,
    jHG=0.25,
    jHT=pars$jHT0,
    rNH=jHT * pars$nNH * pars$sigmaNH,
    rCH=jHT * pars$sigmaCH,
    dH.Hdt=pars$jHGm
  )

  # Symbiont
  S <- dplyr::data_frame(
    S=pars$initS,
    jL=env$L * pars$astar,
    jCP=max(0, synth(jL * pars$yCL, H$jCO2[1]*H$H/S, pars$jCPm), na.rm=T),
    jeL=max(jL - jCP/pars$yCL, 0),
    jNPQ=pars$kNPQ,
    jCO2w=H$jCO2*H$H/S - jCP,
    jSG=pars$jSGm/10,
    rhoC=jCP,
    jNw=0,
    jST=pars$jST0,
    rNS=jST * pars$nNS * pars$sigmaNS,
    rCS=jST * pars$sigmaCS,
    cROS=1,
    dS.Sdt=pars$jSGm)


  # Run simulation by updating
  # ==========================
  grss <- FALSE
  shss <- FALSE
  t <- 2

  while (grss==FALSE | shss==FALSE) {

    # Add row of NAs to be filled by next time step update
    S[t,] <- NA
    H[t,] <- NA

    # Photosynthesis
    # ==============
    # Light input flux
    S$jL[t] <- (1.256307 + 1.385969 * exp(-6.479055 * (S$S[t-1]/H$H[t-1]))) * env$L * pars$astar
    # CO2 input flux
    S$rCS[t] <- pars$sigmaCS * (pars$jST0 + (1-pars$yC)*S$jSG[t-1]/pars$yC)  # metabolic CO2 recycled from symbiont biomass turnover
    H$rCH[t] <- pars$sigmaCH * (H$jHT[t-1] + (1-pars$yC)*H$jHG[t-1]/pars$yC)  # metabolic CO2 recycled from host biomass turnover
    H$jCO2[t] <- pars$kCO2 * H$jeC[t-1]  # carbon not used in host biomass is used to activate CCM's that deliver CO2 to photosynthesis
    # Production flux (photosynthetic carbon fixation)
    S$jCP[t] <- synth(S$jL[t] * pars$yCL, (H$jCO2[t] + H$rCH[t])*H$H[t-1]/S$S[t-1] + S$rCS[t], pars$jCPm) / S$cROS[t-1]
    # Rejection flux: CO2 (wasted to the environment)
    S$jCO2w[t] <- max((H$jCO2[t] + H$rCH[t]) * H$H[t-1]/S$S[t-1] + S$rCS[t] - S$jCP[t], 0)
    # Rejection flux: excess light energy not quenched by carbon fixation
    S$jeL[t] <- max(S$jL[t] - S$jCP[t] / pars$yCL, 0)
    # Amount of excess light energy quenched by NPQ
    S$jNPQ[t] <- (pars$kNPQ^(-1) + S$jeL[t]^(-1))^(-1/1)
    # Scaled ROS production due to excess excitation energy (=not quenched by carbon fixation AND NPQ)
    S$cROS[t] <- 1 + (max(S$jeL[t] - S$jNPQ[t], 0) / pars$kROS)^pars$k

    # Symbiont biomass
    # ================
    # Nitrogen input flux
    S$rNS[t] <- pars$jST0 * pars$nNS * pars$sigmaNS  # Recylced N from symbiont biomass turnover.
    H$rhoN[t-1] <- H$rhoN[t-1]  # Nitrogen shared from the host (defined below, so previous time step used)
    # Carbon input flux
    S$jCP[t] <- S$jCP[t]  # Production of fixed carbon from photosynthesis SU
    # Production flux (symbiont biomass formation)
    S$jSG[t] <- synth(pars$yC*S$jCP[t], (H$rhoN[t-1]*H$H[t-1]/S$S[t-1] + S$rNS[t])/pars$nNS, pars$jSGm)
    # Rejection flux: carbon (surplus carbon shared with the host)
    S$rhoC[t] <- max(S$jCP[t] - S$jSG[t]/pars$yC, 0)
    # Rejection flux: nitrogen (surplus nitrogen wasted to the environment)
    S$jNw[t] <- max(H$rhoN[t-1]*H$H[t-1]/S$S[t-1] + S$rNS[t] - pars$nNS * S$jSG[t], 0)
    # Symbiont biomass loss (turnover)
    S$jST[t] <- pars$jST0 * (1 + pars$b * (S$cROS[t] - 1))

    # Host biomass
    # ============
    # Food input flux (prey=both carbon and nitrogen)
    H$jX[t] <- (pars$jXm * env$X / (env$X + pars$KX))  # Prey uptake from the environment
    # Nitrogen input flux
    H$jN[t] <- (pars$jNm * env$N / (env$N + pars$KN))  # N uptake from the environment
    H$rNH[t] <- H$jHT[t-1] * pars$nNH * pars$sigmaNH  # Recycled N from host biomass turnover
    # Carbon input flux
    S$rhoC[t] <- S$rhoC[t]  # Carbon shared by the symbiont (defined above)
    # Production flux (host biomass formation)
    H$jHG[t] <- synth(pars$yC*(S$rhoC[t]*S$S[t-1]/H$H[t-1] + H$jX[t] + jOC), (H$jN[t] + pars$nNX*H$jX[t] + H$rNH[t]) / pars$nNH, pars$jHGm)
    # Rejection flux: nitrogen (surplus nitrogen shared with the symbiont)
    H$rhoN[t] <- max(H$jN[t] + pars$nNX * H$jX[t] + H$rNH[t] - pars$nNH * H$jHG[t], 0)
    # Rejection flux: carbon -- given back to symbiont as CO2 input to photosynthesis
    H$jeC[t] <- max(H$jX[t] + jOC + S$rhoC[t]*S$S[t-1]/H$H[t-1] - H$jHG[t]/pars$yC, 0)
    # Host biomass loss
    H$jHT[t] <- pars$jHT0

    # State equations
    # ===============
    # Specific growth rates (Cmol/Cmol/d)
    S$dS.Sdt[t] <- S$jSG[t] - S$jST[t]
    H$dH.Hdt[t] <- H$jHG[t] - H$jHT[t]
    # Biomass (Cmol)
    S$S[t] <- S$S[t-1] + S$dS.Sdt[t] * S$S[t-1] * dt
    H$H[t] <- H$H[t-1] + H$dH.Hdt[t] * H$H[t-1] * dt

    # Test if steady state has been reached
    if (t > 20/dt) {
      grss <- ifelse(abs(H$dH.Hdt[t] - H$dH.Hdt[t-10/dt]) < 0.00001, T, F)
      shss <- ifelse(abs(S$S[t]/H$H[t] - S$S[t-10/dt]/H$H[t-10/dt]) < 0.00001, T, F)
    }

    # If system is oscillating after 1000 days, classify as steady state with no growth
    if (t > 1000/dt) {
     if (sum(diff(sign(diff(S$S/H$H)))!=0)>=3) {  # if slope of S/H has changed three times or more...then oscillating
       grss <- T; shss <- T
       H$dH.Hdt[t] <- 0
     }
    }

    # Increment time
    t <- t + 1
  }

  # Return results
  # ==============
  return(list(env=env, pars=pars, H=H, S=S, dt=dt))
}

