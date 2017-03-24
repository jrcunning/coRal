#' Parallel complementary synthesizing unit.
#'
#' Given two substrate input fluxes, \code{synth} calculates the production flux
#' of a parallel complementary synthesizing unit. See page 105, Fig. 3.7, of:
#' Kooijman, S. A. L. M. (2010). Dynamic Energy Budget Theory for Metabolic
#' Organisation. Cambridge University Press. 3rd edition.
#'
#' @param x,y The two substrate input fluxes to the synthesizing unit as
#'   length-one numerics.
#' @param m The maximum production flux of the synthesizing unit.
#' @return The production flux as a length-one numeric.
#' @examples
#' synth(1, 1, 10)

synth <- function(x, y, m) {
  1 / ((1 / m) + (1 / x) + (1 / y) - (1 / (x + y)))
}
