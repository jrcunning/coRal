#' Initialize abiotic environment.
#'
#' This function returns a list of vectors that define the external irradiance
#' (L), DIN concentration (N), and prey concentration (X) to be used as
#' subsequent inputs for a model simulation.
#'
#' The options for the functional form of each environmental factor are:
#' \itemize{
#'   \item 1=linear increase from min to max
#'   \item 2=linear decrease from max to min
#'   \item 3=sinusoid with period of one year and a range from min to max
#' }
#'
#' @param time A vector of time steps at which the model should be evaluated
#'   (units=days) (e.g., seq(0, 365, 0.1))
#' @param L A vector with length 3 defining external light (min, max, functional
#'   form (see details)). Units=mol photons m^-2 d^-1.
#' @param N A vector with length 3 defining external DIN concentration (min,
#'   max, functional form (see details)). Units=mol N L^-1.
#' @param X A vector with length 3 defining external prey concentration (min,
#'   max, functional form (see details)). Units=C-mol X L^-1.
#' @return A list of 3 numeric vectors (each with length=length(time))
#'   corresponding to light, DIN, and prey environmental forcing functions.
#' @examples
#' env1 <- init_env(time=seq(1,365,0.1), L=c(20,40,1), N=c(1e-7,1e-7,0), X=c(0,0,0))

init_env <- function(time, L, N, X) {

  # Light
  L <- if (L[3]==0) {
    seq(L[1], L[2], along=time)
  } else if (L[3]==1) {
    seq(L[2], L[1], along=time)
  } else if (L[3]==2) {
    0.5 * (L[2] - L[1]) * sin(0.0172*time) + L[1] + 0.5 * (L[2] - L[1])
  } else if (L[3]==3) {
    0.5 * (L[2] - L[1]) * sin(0.0172*(time-182)) + L[1] + 0.5 * (L[2] - L[1])
  } else if (L[3]==4) {
    f <- 0.5 * (L[2] - L[1]) * sin(0.0172*time)
    f <- scales::rescale(f, to=c(50, L[2]))
    ff <- scales::rescale(1 + (0.01 / (1 + exp(0.05*(time-200)))), to=c(0.6,1))
    f * ff
  }

  # DIN
  N <- if (N[3]==0) {
    seq(N[1], N[2], along=time)
  } else if (N[3]==1) {
    seq(N[2], N[1], along=time)
  } else {
    0.5 * (N[2] - N[1]) * sin(0.0172*time) + N[1] + 0.5 * (N[2] - N[1])
  }

  # Prey
  X <- if (X[3]==0) {
    seq(X[1], X[2], along=time)
  } else if (X[3]==1) {
    seq(X[2], X[1], along=time)
  } else {
    0.5 * (X[2] - X[1]) * sin(0.0172*time) + X[1] + 0.5 * (X[2] - X[1])
  }
  # Set environment specifications
  env <- list(L=L, N=N, X=X)
  return(env)
}


