#' Plot model simulations
#'
#' The \code{plot_run} function will plot a multipanel figure with all model
#' fluxes. Other functions can be used to plot single environmental factors or
#' model fluxes independently.
#'
#' @param run The model simulation to plot from. This can be the object returned
#'   from \code{run_coral} or \code{run_coral_ss}.
#'
plot_run <- function(run) {
  # If run was a steady state run, set up time and environment vectors for plotting
  if (is.null(run$time)) {
    run <- within(run, {
      time <- seq(1, nrow(S))
      env$L <- rep(env$L, last(time))
      env$N <- rep(env$N, last(time))
      env$X <- rep(env$X, last(time))
      time <- time*dt
    })
  }
  # Plotting
  with(run, {
    # Set up graphical parameters
    par(mar=c(2,2,1,1), oma=c(0,0,0,0), mfrow=c(5,2), mgp=c(1.2,0.2,0), tck=0.025, lwd=1, xaxs="i",
        cex.main=1.5, cex.axis=1, cex.lab=1)

    # External irradiance
    plot(time, env$L, type="l", col="yellow", ylim=c(0,60), lwd=2, xlab="", ylab="mol/m2/d")
    title("Light", adj=0.05, line=0)

    # External DIN concentration
    plot(time, env$N * 10^6, type="l", col="blue", ylim=c(0,4), ylab="µmol/L", xlab="", lwd=2)
    title("DIN", adj=0.05, line=0)

    # Specific growth rates of host and symbiont
    Hgrf <- H$dH.Hdt[length(H$dH.Hdt)]
    plot(time, H$dH.Hdt, type="l", ylim=c(min(0, min(c(H$dH.Hdt[-(1:10)], S$dS.Sdt[-(1:10)]))), max(c(H$dH.Hdt[-(1:10)], S$dS.Sdt[-(1:10)]))),
         xlab="", ylab="d-1", lwd=2, cex=1, cex.lab=1)
    if(any(c(H$dH.Hdt[-(1:10)], S$dS.Sdt[-(1:10)])<0)) abline(h=0, lwd=0.5, lty=3)
    text(time[0.95*length(time)], Hgrf, labels=as.character(round(Hgrf, 3)), pos=3, xpd=T, cex=0.75)
    #title("State variable dynamics", adj=0, cex.main=2, outer = T)
    title("Specific growth rate", adj=0.05, line=0)
    lines(time, S$dS.Sdt, col="black", lwd=1)
    legend("topright", legend=c("Host", "Sym"), lwd=c(3,1), col="black", bty="n", y.intersp=1)

    # Symbiont to host ratio
    totSH <- S$S / H$H
    totSHf <- totSH[length(totSH)]
    plot(time, totSH, type="l", col="black", ylim=c(0, max(totSH, na.rm=T)), ylab="C-molS/C-molH", xlab="", lwd=2)
    text(time[0.95*length(time)], totSHf, labels=as.character(round(totSHf, 3)), pos=3, xpd=T, cex=0.75)
    title("Symbiont to host ratio", adj=0.05, line=0)

    # Light quenching (carbon fixation, NPQ, and excess (=ROS producing))
    plot(NA, xlim=range(time), xlab="", ylim=c(0, max(c(S$jL, S$jeL), na.rm=T)), ylab="mol photons/C-molS/d")
    title("Light quenching", adj=0.05, line=0)
    lines(time, S$jCP/pars$yCL, col="yellow", lwd=1) # amt. quenched by photosynthesis
    lines(time, S$jCP/pars$yCL + S$jNPQ, col="yellow", lty=1, lwd=2) # amt. quenched by photo + NPQ
    lines(time, S$jL, col="yellow", lty=3, lwd=2) # total amt. absorbed
    legend("topright", legend=c("Excess", "NPQ", "Photo."), lty=c(3,1,1), lwd=c(3,3,1), col="yellow", bty="n")

    # Symbiont biomass SU dynamics (proportions of C and N rejected from SU)
    plot(NA, xlim=range(time), xlab="", ylim=c(0, 1), ylab="Proportion of uptake")
    title("Symbiont: substrate excess", adj=0.05, line=0)
    lines(time, S$rhoC/S$jCP, col="red", lwd=2)
    lines(time, S$jNw/(H$rhoN*H$H/S$S + S$rNS), col="blue", lwd=2)
    legend("topright", legend=c("C", "N"), lwd=2, bty="n", col=c("red", "blue"))

    # Photosynthesis SU dynamics (proportions of light and CO2 rejected from SU)
    pl <- log(pmin(((H$jCO2+H$rCH)*H$H/S$S + S$rCS), pars$jCPm) / pmin((S$jL*run$pars$yCL), pars$jCPm))
    maxabs <- max(abs(pl))
    plot(NA, xlim=range(time), xlab="", ylab="", ylim=c(-maxabs, maxabs))
    title("Photosynthesis C- vs. L-limitation", adj=0.05, line=0)
    mtext(side=4, line=-1, "C-lim.  L-lim.", cex=0.5)
    lines(time, pl, col="black", lty=1, lwd=2)
    abline(h=0, lty=3)

    # Photosynthesis rate (carbon fixation)
    plot(NA, xlim=range(time), xlab="", ylim=c(0, max(S$jCP)), ylab="C-mol/C-molS/d")
    title("Photosynthesis rate", adj=0.05, line=0)
    lines(time, S$jCP, col="red", lwd=2)

    # Relative ROS production due to excitation energy in excess of carbon fixation and NPQ
    plot(NA, xlim=range(time), xlab="", ylim=c(1, max(2, max(S$cROS))), ylab="Relative to baseline")
    title("ROS production", adj=0.05, line=0)
    lines(time, S$cROS, col="orange", lwd=2)

    # Host biomass SU dynamics (proportions of C and N rejected from SU)
    plot(NA, xlim=range(time), xlab="", ylim=c(0, 1), ylab="Proportion of uptake")
    title("Coral: substrate excess", adj=0.05, line=0)
    lines(time, H$jeC/(pars$yC*S$rhoC*S$S/H$H + H$jX), col="red", lwd=2)
    lines(time, H$rhoN/(H$jN + H$jX*pars$nNX + H$rNH), col="blue", lwd=2)
    legend("topright", legend=c("C", "N"), lwd=2, bty="n", col=c("red", "blue"))
  })
}

#' @describeIn plot_run Plot external irradiance
plot_L <- function(run) with(run, {
  plot(time, env$L, type="l", col="yellow", ylim=c(0,60), lwd=2, xlab="", ylab="mol/m2/d", main="Irradiance")
})

#' @describeIn plot_run Plot external DIN concentration
plot_DIN <- function(run) with(run, {
  plot(time, env$N * 10^6, type="l", col="blue", ylim=c(0,4), lwd=2, xlab="", ylab="µmol/L", main="DIN")
})

#' @describeIn plot_run Plot specific growth rates of host and symbiont
plot_gr <- function(run) with(run, {
  Hgrf <- H$dH.Hdt[length(H$dH.Hdt)]
  plot(time, H$dH.Hdt, type="l", ylim=c(min(0, min(c(H$dH.Hdt[-(1:10)], S$dS.Sdt[-(1:10)]))), max(c(H$dH.Hdt[-(1:10)], S$dS.Sdt[-(1:10)]))),
       xlab="", ylab="d-1", lwd=2, cex=1, cex.lab=1,
       main="Specific growth rate")
  if(any(c(H$dH.Hdt[-(1:10)], S$dS.Sdt[-(1:10)])<0)) abline(h=0, lwd=0.5, lty=3)
  text(time[0.95*length(time)], Hgrf, labels=as.character(round(Hgrf, 3)), pos=3, xpd=T, cex=0.75)
  lines(time, S$dS.Sdt, col="black", lwd=1)
  legend("topright", legend=c("Host", "Sym"), lwd=c(3,1), col="black", bty="n", y.intersp=1)
})

#' @describeIn plot_run Plot the symbiont to host biomass ratio
plot_sh <- function(run) with(run, {
  totSH <- S$S / H$H
  totSHf <- totSH[length(totSH)]
  plot(time, totSH, type="l", col="black", ylim=c(0, max(totSH, na.rm=T)), ylab="C-molS/C-molH", xlab="", lwd=2,
       main="Symbiont:host biomass")
  text(time[0.95*length(time)], totSHf, labels=as.character(round(totSHf, 3)), pos=3, xpd=T, cex=0.75)
})

#' @describeIn plot_run Plot light quenching (carbon fixation, NPQ, and excess (=ROS-producing))
plot_Lq <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(0, max(c(S$jL, S$jeL), na.rm=T)), ylab="mol photons/C-molS/d",
       main="Light quenching")
  lines(time, S$jCP/pars$yCL, col="yellow", lwd=1) # amt. quenched by photosynthesis
  lines(time, S$jCP/pars$yCL + S$jNPQ, col="yellow", lty=1, lwd=2) # amt. quenched by photo + NPQ
  lines(time, S$jL, col="yellow", lty=3, lwd=2) # total amt. absorbed
  legend("topright", legend=c("Excess", "NPQ", "Photo."), lty=c(3,1,1), lwd=c(3,3,1), col="yellow", bty="n")
})

#' @describeIn plot_run Plot photosynthesis rate (carbon fixation; jCP)
plot_photo <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(0, max(S$jCP)), ylab="C-mol/C-molS/d")
  title("Photosynthesis rate", adj=0.05, line=0)
  lines(time, S$jCP, col="red", lwd=2)
})

#' @describeIn plot_run Plot limitation coefficient for photosynthesis SU
plot_pl <- function(run) with(run, {
  pl <- log(   pmin((H$jCO2 + H$rCH)*H$H/S$S + S$rCS, pars$jCPm)   /   pmin(S$jL * pars$yCL, pars$jCPm)    )
  maxabs <- max(abs(pl))
  plot(NA, xlim=range(time), xlab="", ylab="", ylim=c(-maxabs, maxabs))
  title("Photosynthesis C- vs. L-limitation", adj=0.05, line=0)
  lines(time, pl, col="black", lty=1, lwd=2)
  abline(h=0, lty=3)
  mtext(side=4, line=-1, "C-lim.  L-lim.", cex=0.5)
})

#' @describeIn plot_run Plot photosynthesis SU dynamics (proportions of light and CO2 rejected from SU)
plot_PSU <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(0,1), ylab="Proportion of uptake")
  title("Photosynthesis: substrate excess", adj=0.05, line=0)
  lines(time, S$jeL/S$jL, col="yellow", lwd=2)
  lines(time, S$jCO2w/((H$jCO2 + H$rCH)*H$H/S$S + S$rCS), col="red", lwd=2)
  legend("topright", legend=c("Light", "DIC"), lty=c(1,1), lwd=2, bty="n", col=c("yellow", "red"))
})

#' @describeIn plot_run Plot relative ROS production due to excitation energy in excess of carbon fixation and NPQ
plot_ROS <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(1, max(2, max(S$cROS))), ylab="Relative to baseline")
  title("ROS production", adj=0.05, line=0)
  lines(time, S$cROS, col="orange", lwd=2)
})

#' @describeIn plot_run Plot symbiont biomass SU dynamics (proportions of C and N rejected from SU)
plot_symSU <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(0, 1), ylab="Proportion of uptake")
  title("Symbiont: substrate excess", adj=0.05, line=0)
  lines(time, S$rhoC/S$jCP, col="red", lwd=2)
  lines(time, S$jNw/(H$rhoN*H$H/S$S + S$rNS), col="blue", lwd=2)
  legend("topright", legend=c("C", "N"), lwd=2, bty="n", col=c("red", "blue"))
})

#' @describeIn plot_run Plot symbiont biomass SU limitation coefficient
plot_sl <- function(run) with(run, {
  sl <- log(  pmin(pars$yC*S$jCP, pars$jSGm)   /   pmin((H$rhoN*H$H/S$S + S$rNS)/pars$nNS, pars$jSGm)  )
  maxabs <- max(abs(sl))
  plot(NA, xlim=range(time), xlab="", ylab="", ylim=c(-maxabs, maxabs))
  title("Symbiont C- vs. N-limitation", adj=0.05, line=0)
  lines(time, sl, col="black", lty=1, lwd=2)
  abline(h=0, lty=3)
  mtext(side=4, line=-1, "C-lim.  N-lim.", cex=0.5)
})


#' @describeIn plot_run Plot host biomass SU dynamics (proportions of C and N rejected from SU)
plot_corSU <- function(run) with(run, {
  plot(NA, xlim=range(time), xlab="", ylim=c(0, 1), ylab="Proportion of uptake")
  title("Coral: substrate excess", adj=0.05, line=0)
  lines(time, H$jeC/(pars$yC*S$rhoC*S$S/H$H + H$jX), col="red", lwd=2)
  lines(time, H$rhoN/(H$jN + H$jX*pars$nNX + H$rNH), col="blue", lwd=2)
  legend("topright", legend=c("C", "N"), lwd=2, bty="n", col=c("red", "blue"))
})

#' @describeIn plot_run Plot host biomass SU limitation coefficient
plot_hl <- function(run) with(run, {
  hl <- log(  pmin(pars$yC*S$rhoC*S$S/H$H + H$jX, pars$jHGm) / pmin((H$jN + pars$nNX*H$jX + H$rNH) / pars$nNH, pars$jHGm)  )
  maxabs <- max(abs(hl))
  plot(NA, xlim=range(time), xlab="", ylab="", ylim=c(-maxabs, maxabs))
  title("Host C- vs. N-limitation", adj=0.05, line=0)
  lines(time, hl, col="black", lty=1, lwd=2)
  abline(h=0, lty=3)
  mtext(side=4, line=-1, "C-lim.  N-lim.", cex=0.5)
})

#' @describeIn plot_run Plot both symbiont and host biomass SU limitation coefficients
plot_bl <- function(run) with(run, {
  sl <- log(  pmin(pars$yC*S$jCP, pars$jSGm)   /   pmin((H$rhoN*H$H/S$S + S$rNS)/pars$nNS, pars$jSGm)  )
  hl <- log(  pmin(pars$yC*S$rhoC*S$S/H$H + H$jX, pars$jHGm) / pmin((H$jN + pars$nNX*H$jX + H$rNH) / pars$nNH, pars$jHGm)  )
  maxabs <- max(c(abs(hl), abs(sl)))
  range <- range(c(hl, sl))
  plot(NA, xlim=range(time), xlab="", ylab="", ylim=range)
  title("Host C- vs. N-limitation", adj=0.05, line=0)
  lines(time, hl, col="black", lty=1, lwd=2)
  lines(time, sl, col="black", lty=1, lwd=1)
  abline(h=0, lty=3)
  mtext(side=4, line=-1, "C-lim.  N-lim.", cex=0.5)
})

#' @describeIn plot_run Plot another version of the multipanel run output
plot_run_v2 <- function(run) {
  if (is.null(run$time)) time <- seq(1, nrow(run$S))
  with(run, {
    # Set up graphical parameters
    par(mar=c(2,2,1,1), oma=c(0,0,0,0), mfrow=c(6,1), mgp=c(1,0,0), tck=0.025, lwd=1, xaxs="i",
        cex.main=1.2, cex.axis=1, cex.lab=1)

    # External irradiance
    plot(time, env$L, type="l", col="gold", ylim=c(0,60), lwd=2, xlab="", ylab="mol/m2/d")
    title("A. External irradiance", adj=0, line=0.25)

    # Symbiont to host biomass ratio
    totSH <- S$S / H$H
    totSHf <- totSH[length(totSH)]
    plot(time, totSH, type="l", col="black", ylim=c(0, max(totSH, na.rm=T)), ylab="CmolS/CmolH", xlab="", lwd=2)
    text(time[0.95*length(time)], totSHf, labels=as.character(round(totSHf, 3)), pos=3, xpd=T, cex=0.75)
    title("B. Symbiont:host biomass", adj=0, line=0.25)

    # Specific growth rates of host and symbiont
    Hgrf <- H$dH.Hdt[length(H$dH.Hdt)]
    plot(time, H$dH.Hdt, type="l", ylim=c(min(0, min(c(H$dH.Hdt[-(1:10)], S$dS.Sdt[-(1:10)]))), max(c(H$dH.Hdt[-(1:10)], S$dS.Sdt[-(1:10)]))),
         xlab="", ylab="d-1", lwd=2, cex=1, cex.lab=1)
    if(any(c(H$dH.Hdt[-(1:10)], S$dS.Sdt[-(1:10)])<0)) abline(h=0, lwd=0.5, lty=3)
    text(time[0.95*length(time)], Hgrf, labels=as.character(round(Hgrf, 3)), pos=3, xpd=T, cex=0.75)
    #title("State variable dynamics", adj=0, cex.main=2, outer = T)
    title("C. Specific growth rates", adj=0, line=0.25)
    lines(time, S$dS.Sdt, col="black", lwd=1)
    legend("topright", legend=c("Host", "Sym"), lwd=c(2,1), col="black", bty="n", y.intersp=1, cex=0.75)

    # Light quenching (carbon fixation, NPQ, and excess (=ROS producing))
    plot(NA, xlim=range(time), xlab="", ylim=c(0, max(c(S$jL, S$jeL), na.rm=T)), ylab="mol ph/CmolS/d")
    title("D. Light quenching", adj=0, line=0.25)
    lines(time, S$jCP/pars$yCL, col="gold", lwd=1) # amt. quenched by photosynthesis
    lines(time, S$jCP/pars$yCL + S$jNPQ, col="gold", lty=1, lwd=2) # amt. quenched by photo + NPQ
    lines(time, S$jL, col="gold", lty=3, lwd=2) # total amt. absorbed
    legend("topright", legend=c("ROS", "NPQ", "Photo."), lty=c(3,1,1), lwd=c(2,2,1), col="gold", bty="n", cex=0.75)

    # Photosynthesis: substrate-limitation
    pl <- log(   pmin((H$jCO2 + H$rCH)*H$H/S$S + S$rCS, pars$jCPm)   /   pmin(S$jL * pars$yCL, pars$jCPm)    )
    maxabs <- max(abs(pl))
    plot(NA, xlim=range(time), xlab="", ylab="", ylim=c(-maxabs, maxabs))
    title("E. Photosynthesis: substrate-limitation", adj=0, line=0.25)
    lines(time, pl, col="black", lty=1, lwd=2)
    abline(h=0, lty=3)
    text(par("usr")[2], 0, labels="L-lim.\nCO2-lim.", cex=0.75, adj=1)

    # Biomass formation - substrate-limitation
    sl <- log(  pmin(pars$yC*S$jCP, pars$jSGm)   /   pmin((H$rhoN*H$H/S$S + S$rNS)/pars$nNS, pars$jSGm)  )
    hl <- log(  pmin(pars$yC*S$rhoC*S$S/H$H + H$jX, pars$jHGm) / pmin((H$jN + pars$nNX*H$jX + H$rNH) / pars$nNH, pars$jHGm)  )
    maxabs <- max(c(abs(hl), abs(sl)))
    range <- range(c(hl, sl))
    plot(NA, xlim=range(time), xlab="", ylab="", ylim=range)
    title("F. Biomass: substrate-limitation", adj=0, line=0.25)
    lines(time, hl, col="black", lty=1, lwd=2)
    lines(time, sl, col="black", lty=1, lwd=1)
    abline(h=0, lty=3)
    text(par("usr")[2], 0, labels="N-lim.\nC-lim.", cex=0.75, adj=1)
    legend("topright", legend=c("Host", "Sym"), lty=1, lwd=c(2,1), col="black", bty="n", cex=0.75)

    # Fixed carbon fate
    #plot(NA, xlim=range(time), xlab="", ylab="", ylim=c(0, max(c(S$jCP))))
    #lines(time, S$jCP, lwd=2, lty=2) # Total fixed carbon available to symbiont
    #lines(time, S$jSG, lwd=2) # Fixed carbon used by symbiont

    # nitrogen fate
    #plot(NA, xlim=range(time), xlab="", ylab="", ylim=c(0, max(c(H$jN+H$rNH))))
    #lines(time, H$jN + H$rNH + H$jX*pars$nNX, lwd=2, lty=2) # Total nitrogen available to host
    #lines(time, H$jHG*pars$nNH, lwd=2) # Nitrogen used by host
  })
}
