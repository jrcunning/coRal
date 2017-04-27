#' Run simulation of coral-\emph{Symbiodinium} model.
#'
#' This function runs a simulation of the coral-Symbiodinium bioenergetic model
#' using a specified time vector, environmental inputs, and parameters.
#'
#' @param time A vector of time steps at which the model should be evaluated
#'   (units=days) (e.g., seq(0, 365, 0.1))
#' @param env A list of numeric vectors named L, N, and X, defining the external
#'   light, [DIN], and [prey] for the model simulation. Each of these vectors
#'   should have length=length(time). The object returned from \code{init_env}
#'   can be used here.
#' @param pars A list of (named) parameter values to use in the simulation.
#'   Parameter names must match those returned by \code{def_pars()}.
#' @return A list with the following named elements:
#'   \describe{
#'     \item{time}{The time vector used in the simulation}
#'     \item{env}{The environment object used in the simulation}
#'     \item{pars}{The parameter values used in the simulation}
#'     \item{H}{A tibble of all host biomass-specific model
#'     fluxes at each time step}
#'     \item{S}{A tibble of all symbiont biomass-specific
#'     model fluxes at each time step}
#'   }
#' @seealso \code{\link{init_env}}, \code{\link{def_pars}}
#' @examples
#' time1 <- seq(0, 365, 0.1)
#' pars1 <- def_pars()
#' env1 <- init_env(time=time1, L=c(20,40,1), N=c(1e-7,1e-7,0), X=c(0,0,0))
#' run1 <- run_coral(time=time1, env=env1, pars=pars1)

run_coral <- function(time, env, pars) {

  # Get number of symbionts in run
  nsym <- length(pars$initS)

  # Set initial values
  # ==================
  # Host
  H <- dplyr::data_frame(
    time=time,
    H=pars$initH,
    jX=(pars$jXm * env$X[1] / (env$X[1] + pars$KX)),
    jN=(pars$jNm * env$N[1] / (env$N[1] + pars$KN)),
    rhoN=jN,
    jeC=10,
    jCO2=pars$kCO2 * jeC,
    jHG=0.25,
    jHT=pars$jHT0,
    rNH=jHT * pars$nNH * pars$sigmaNH,
    rCH=jHT * pars$sigmaCH,
    dH.Hdt=pars$jHGm
  )
  H[2:nrow(H), 2:ncol(H)] <- NA

  # Symbiont
  S <- list()
  for (i in 1:nsym) {
    S[[i]] <- dplyr::data_frame(
      time=time,
      S=pars$initS[i],
      jL=env$L[1] * pars$astar[i],
      jCP=max(0, synth(jL * pars$yCL[i], H$jCO2[1]*H$H/S, pars$jCPm[i]), na.rm=T),
      jeL=max(jL - jCP/pars$yCL[i], 0),
      jNPQ=pars$kNPQ[i],
      jCO2w=H$jCO2*H$H/S - jCP,
      jSG=pars$jSGm[i]/10,
      rhoC=jCP,
      jNw=0,
      jST=pars$jST0[i],
      rNS=jST * pars$nNS[i] * pars$sigmaNS[i],
      rCS=jST * pars$sigmaCS[i],
      cROS=1,
      dS.Sdt=pars$jSGm[i])
    S[[i]][2:nrow(S[[i]]), 2:ncol(S[[i]])] <- NA
  }


  # Run simulation by updating
  # ==========================
  for (t in 2:length(time)) {
    S.t <- rowSums(sapply(S, "[[", 2))[t-1]  # Get total symbiont abundance from prev time step
    S.p <- sapply(S, "[[", 2)[t-1, ] / S.t  # Get proportion of each symbiont

    for (i in 1:nsym) {
      # Photosynthesis
      # ==============
      # Light input flux
      S[[i]]$jL[t] <- (1.256307 + 1.385969 * exp(-6.479055 * (S.t/H$H[t-1]))) * env$L[t] * pars$astar[i]
      # CO2 input flux
      S[[i]]$rCS[t] <- pars$sigmaCS[i] * (pars$jST0[i] + (1-pars$yC)*S[[i]]$jSG[t-1]/pars$yC)  # metabolic CO2 recycled from symbiont biomass turnover
      H$rCH[t] <- pars$sigmaCH * (H$jHT[t-1] + (1-pars$yC)*H$jHG[t-1]/pars$yC)  # metabolic CO2 recycled from host biomass turnover
      H$jCO2[t] <- pars$kCO2 * H$jeC[t-1]  # carbon not used in host biomass is used to activate CCM's that deliver CO2 to photosynthesis
      # Production flux (photosynthetic carbon fixation)
      S[[i]]$jCP[t] <- synth(S[[i]]$jL[t] * pars$yCL[i], (H$jCO2[t] + H$rCH[t])*H$H[t-1]/S.t + S[[i]]$rCS[t], pars$jCPm[i]) / S[[i]]$cROS[t-1]
      # Rejection flux: CO2 (wasted to the environment)
      S[[i]]$jCO2w[t] <- max((H$jCO2[t] + H$rCH[t])*H$H[t-1]/S.t + S[[i]]$rCS[t] - S[[i]]$jCP[t], 0)
      # Rejection flux: excess light energy not quenched by carbon fixation
      S[[i]]$jeL[t] <- max(S[[i]]$jL[t] - S[[i]]$jCP[t]/pars$yCL[i], 0)
      # Amount of excess light energy quenched by NPQ
      S[[i]]$jNPQ[t] <- (pars$kNPQ[i]^(-1)+S[[i]]$jeL[t]^(-1))^(-1/1)  # single substrate SU
      #S$jNPQ[t] <- min(S$jeL[t], pars$kNPQ)  # minimum rule
      # Scaled ROS production due to excess excitation energy (=not quenched by carbon fixation AND NPQ)
      S[[i]]$cROS[t] <- 1 + (max(S[[i]]$jeL[t] - S[[i]]$jNPQ[t], 0) / pars$kROS[i])^pars$k[i]

      # Symbiont biomass
      # ================
      # Nitrogen input flux
      S[[i]]$rNS[t] <- pars$jST0[i] * pars$nNS[i] * pars$sigmaNS[i]  # Recylced N from symbiont biomass turnover.
      H$rhoN[t-1] <- H$rhoN[t-1]  # Nitrogen shared from the host (defined below, so previous time step used)
      # Carbon input flux
      S[[i]]$jCP[t] <- S[[i]]$jCP[t]  # Production of fixed carbon from photosynthesis SU
      # Production flux (symbiont biomass formation)
      S[[i]]$jSG[t] <- synth(pars$yC*S[[i]]$jCP[t], (H$rhoN[t-1]*H$H[t-1]/S.t + S[[i]]$rNS[t])/pars$nNS[i], pars$jSGm[i])
      # Rejection flux: carbon (surplus carbon shared with the host)
      S[[i]]$rhoC[t] <- max(S[[i]]$jCP[t] - S[[i]]$jSG[t]/pars$yC, 0)
      # Rejection flux: nitrogen (surplus nitrogen wasted to the environment)
      S[[i]]$jNw[t] <- max(H$rhoN[t-1]*H$H[t-1]/S.t + S[[i]]$rNS[t] - pars$nNS[i] * S[[i]]$jSG[t], 0)
      # Symbiont biomass loss (turnover)
      S[[i]]$jST[t] <- pars$jST0[i] * (1 + pars$b[i] * (S[[i]]$cROS[t] - 1))
    }

    # Host biomass
    # ============
    # Food input flux (prey=both carbon and nitrogen)
    H$jX[t] <- (pars$jXm * env$X[t] / (env$X[t] + pars$KX))  # Prey uptake from the environment
    # Nitrogen input flux
    H$jN[t] <- (pars$jNm * env$N[t] / (env$N[t] + pars$KN))  # N uptake from the environment
    H$rNH[t] <- H$jHT[t-1] * pars$nNH * pars$sigmaNH  # Recycled N from host biomass turnover
    # Carbon input flux
    rhoC.t <- sum(sapply(S, function(S) with(S[t-1, ], rhoC*S)))  # Total amount of carbon shared by all symbionts
    #S[[i]]$rhoC[t] = Carbon shared by the symbiont (defined above)
    # Production flux (host biomass formation)
    H$jHG[t] <- synth(pars$yC*(rhoC.t/H$H[t-1] + H$jX[t]), (H$jN[t] + pars$nNX*H$jX[t] + H$rNH[t]) / pars$nNH, pars$jHGm)
    # Rejection flux: nitrogen (surplus nitrogen shared with the symbiont)
    H$rhoN[t] <- max(H$jN[t] + pars$nNX * H$jX[t] + H$rNH[t] - pars$nNH * H$jHG[t], 0)
    # Rejection flux: carbon -- given back to symbiont as CO2 input to photosynthesis
    H$jeC[t] <- max(H$jX[t] + rhoC.t/H$H[t-1] - H$jHG[t]/pars$yC, 0)
    # Host biomass loss
    H$jHT[t] <- pars$jHT0

    # State equations
    # ===============
    # Specific growth rates (Cmol/Cmol/d)
    for (i in 1:nsym) S[[i]]$dS.Sdt[t] <- S[[i]]$jSG[t] - S[[i]]$jST[t]
    H$dH.Hdt[t] <- H$jHG[t] - H$jHT[t]
    # Biomass (Cmol)
    for (i in 1:nsym) S[[i]]$S[t] <- S[[i]]$S[t-1] + S[[i]]$dS.Sdt[t] * S[[i]]$S[t-1] * (time[2] - time[1])
    H$H[t] <- H$H[t-1] + H$dH.Hdt[t] * H$H[t-1] * (time[2] - time[1])

  }

  # Return results
  # ==============
  if (length(S)==1) S <- S[[1]]
  return(list(time=time, env=env, pars=pars, H=H, S=S))
}

