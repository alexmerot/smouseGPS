#' @title R6 class that creates b-splines.
#'
#' @description Creates b-splines given some set of knot points, and then
#' evaluates the splines given some set of coefficients for the splines. It has
#' 3 arguments initialisation : \code{f = BSpline(K, t_knot, m)}.
#'
#' @field t_knot Spline knot points.
#' @field K Order of polynomial.
#' @field m Spline coefficients (\eqn{m = M \cdot D}).
#' @field B Numeric
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
    B = vector("numeric"),
    t_pp = NULL,
    C = NULL,

    #' @description Instantiate an object of the class \code{bspline}.
    #' @param t_knot Spline knot points.
    #' @param K Order of polynomial.
    #' @param m Spline coefficients (\eqn{m = M \cdot D}).
    initialize = function(t_knot, K, m) {
      stopifnot("t_knot must be numeric." = is.numeric(t_knot) & !is.na(t_knot))
      stopifnot("K must be numeric." = is.numeric(K) & !is.na(K))
      stopifnot("m must be numeric." = is.numeric(m) & !is.na(m))

      self$K = K
      self$m = m
      self$t_knot = t_knot
    }
  ),

  private = list(
    domain = NULL,
    S = NULL
  ),

  active = list(
    #' @field get_S Get \eqn{S = K - 1}
    get_S = function() self$S = self$K - 1,

    #' @field get_domain Get the domain of \code{t_knot}
    get_domain = function() {
      self$domain = c(start = self$t_knot[1], end = self$t_knot[length(self$t_knot)])
    }
  )
)
