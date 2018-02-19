#' Simulate and plot a \emph{Symbiodinium} PI curve.
#'
#' Given a set of parameters defining \emph{Symbiodinium} photobiology (e.g.,
#' yCL, jCPm, kROS, kNPQ, k), \code{sym_PI} calculates and plots the rate of
#' photosynthesis (jCP), NPQ (jNPQ), and ROS production (cROS) as a function of
#' irradiance. This simulation is meant to show the photosynthetic performance
#' of a given symbiont independent from a host (i.e., in a state that is not
#' influenced by host carbon-concentrating mechanisms), so the "jCO2" parameter
#' is set at a very high value.
#'
#' @param pars A list of model parameters to use in the simulation, i.e. the
#'   object returned by \code{def_pars}.
#' @param draw Whether to plot PI curve (TRUE or FALSE)
#' @return A data frame with simulation results at each time step.
#' @examples
#' PI1 <- sym_PI(def_pars())
#'
sym_PI <- function(pars, draw=TRUE) {

  with(pars, {
    # Set light range for PI curve
    time <- seq(1,100,0.01)
    jL <- seq(0,250,along=time)
    # Set carbon delivery rate (mol C / CmolS / d)
    jCO2 <- 1000  # Extremely high carbon delivery - no carbon limitation
    # Initialize values...
    jCP <- c(0, rep(NA, length(time)-1))
    jeL <- c(0, rep(NA, length(time)-1))
    jNPQ <- c(0, rep(NA, length(time)-1))
    cROS <- c(1, rep(NA, length(time)-1))
    # Time-stepping solution...
    for (t in 2:length(time)) {
      # Calculate photosynthesis rate
      jCP[t] <- synth(jCO2, jL[t]*pars$yCL, pars$jCPm)/cROS[t-1]
      # Calculate light in excess of photosynthetic quenching
      jeL[t] <- max(0, jL[t] - jCP[t]/pars$yCL)
      # Calculate light energy quenched by NPQ
      jNPQ[t] <- (pars$kNPQ^(-1)+jeL[t]^(-1))^(-1/1)
      # Calculate ROS (cROS) produced due to excess light
      cROS[t] <- 1 + ((jeL[t] - jNPQ[t]) / pars$kROS)^pars$k
    }
    # Return ROS and photosynthesis rate for plotting
    par(mfrow=c(1,1), mar=c(3,3,3,7), mgp=c(1.2,0,0), cex=1, tck=0.025, xaxs="i")
    if (draw) {
      plot(jL, jCP/pars$yCL, xlab="Light (mol photons/C-molS/d)", ylab="",
           type="l", lwd=3, col="red")
      mtext(side=2, text="Photochemical quenching", cex=1, line=1.8)
      mtext(side=2, text="(mol photons/CmolS/d)", cex=0.8, line=1)
      par(new=T)
      plot(jL, jNPQ, type="l", lwd=3, axes=F, xlab="", ylab="")
      axis(side=4); mtext(side=4, text="Non-photochemical quenching", line=1, cex=1)
                    mtext(side=4, text="(mol photons/CmolS/d)", line=1.8, cex=0.8)
      par(new=T)
      plot(jL, cROS, type="l", lwd=3, axes=F, col="orange", xlab="", ylab="", ylim=c(1, max(cROS)*1.5))
      axis(side=4, line=4); mtext(side=4, line=5, text = "ROS production (relative)", cex=1)
      legend("bottomright", legend=c("Photo.", "NPQ", "ROS"), lwd=2, col=c("red", "black", "orange"),
            inset=c(0.1, 0.05), bty="n")
    }
    return(data.frame(time, jL, jCP, jNPQ, cROS))
  })

}


