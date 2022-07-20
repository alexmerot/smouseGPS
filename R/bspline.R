#' @title R6 class that creates b-splines.
#'
#' @description Creates b-splines given some set of knot points, and then
#' evaluates the splines given some set of coefficients for the splines. It has
#' 3 arguments initialisation : \code{f = BSpline(K, t_knot, m)}.
#'
#' @field t_knot Spline knot points.
#' @field K Order of polynomial.
#' @field m Spline coefficients (\eqn{m = M \cdot D}).
#' @field B
#' @field x_mean If set, these will be used to scale the output (default to 0).
#' @field x_std \eqn{x_{out} = x_{std} \cdot (X \cdot m) + x_{mean}}
#' @field t_pp pp break points (\code{size(t_pp) = length(t_knot) - 2*K + 1}).
#' @field C Piecewise polynomial coefficients (\code{size(C) = [length(t_pp)-1, K]}).
#'
#' @return
#' @importFrom R6 R6Class
#' @export
#'
#' @examples

bspline <- R6Class(
  classname = "bspline",

  public = list(
    t_knot = NULL,
    K = NULL,
    m = NULL,

    x_mean = 0,
    x_std = 1,
    B = list(),
    t_pp = NULL,
    C = NULL,

    #' @description Instantiate an object of the class \code{bspline}.
    #' @param t_knot Spline knot points.
    #' @param K Order of polynomial.
    #' @param m Spline coefficients (\eqn{m = M \cdot D}).
    initialize = function(t_knot, K, m) {
      stopifnot("t_knot must be numeric." = is.numeric(t_knot))
      stopifnot("K must be numeric." = is.numeric(K))
      stopifnot("m must be numeric." = is.numeric(m))

      self$K = K
      self$m = m
      self$t_knot = t_knot
    }
  ),

  private = list(

  )
)
