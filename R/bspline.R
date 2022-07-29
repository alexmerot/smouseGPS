#' @title R6 class that creates b-splines.
#'
#' @description Creates b-splines given some set of knot points, and then
#' evaluates the splines given some set of coefficients for the splines. It has
#' 3 arguments initialisation : \code{f = BSpline(K, t_knot, m)}.
#'
#' @field t_knot a Matrix. Spline knot points.
#' @field K an integer. Order of polynomial.
#' @field m a numeric. Spline coefficients (\eqn{m = M \cdot D}).
#' @field B a numeric array with 3 dimension.
#' @field x_mean a numeric. If set, these will be used to scale the output
#' (default to 0).
#' @field x_std a numeric. \eqn{x_{out} = x_{std} \cdot (X \cdot m) + x_{mean}}
#' @field t_pp a numeric columnwise matrix. pp break points (\code{size(t_pp) = nrow(t_knot) - 2*K + 1}).
#' @field C a vector. Piecewise polynomial coefficients
#' (\code{size(C) = [nrow(t_pp)-1, K]}).
#'
#' @return b-splines
#'
#' @importFrom R6 R6Class
#' @importFrom Rfast sort_mat colRanks
#' @importFrom pracma polyval mrdivide isempty numderiv numdiff
#' @export

bspline <- R6Class(
  classname = "bspline",
  public = list(
    t_knot = NULL,
    K = NULL,
    m = NULL,
    B = NULL,
    x_mean = 0,
    x_std = 1,
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

      self$t_knot <- t_knot
      self$K <- K
      self$m <- m
    },

    # Overload the subscript operator '['.
    # `[.self` <- function(x, i = NULL, j = NULL, ...) {
    #   if(nrow(idx) >= 1) t <- idx[[1]]
    #
    #   if(nrow(idx >= 2)) {
    #     num_derivatives <- idx[[2]]
    #   } else {
    #     num_derivatives <- 0
    #   }
    #
    #   return(self$value_at_points(t, num_derivatives))
    # },
    #
    # `[<-.self` <- function(x, i = NULL, j = NULL, value) {
    #
    # },

    #' @description Get the value at points.
    #' @param t a matrix. Points locations.
    #' @param num_derivatives an integer. Numero of the derivatives.
    value_at_points = function(t, num_derivatives) {
      if (!exists("num_derivatives")) num_derivatives <- 0

      x_out <- private$evaluate_from_pp_coeff(t, self$C, self$t_pp, num_derivatives)

      if (!isempty(self$x_std)) x_out <- self$x_std * x_out

      if (!isempty(self$x_mean) & num_derivatives == 0) x_out <- x_out + self$x_mean

      return(x_out)
    }
  ),
  private = list(
    domain = NULL,
    S = NULL
  ),
  active = list(
    #' @field get_S Get \eqn{S = K - 1}
    get_S = function() private$S <- self$K - 1,

    #' @field get_domain Get the domain of \code{t_knot}
    get_domain = function() {
      private$domain <- list(
        start = self$t_knot[1],
        end = self$t_knot[nrow(self$t_knot)]
      )
      attr(private$domain, "class") <- "domain"

      return(private$domain)
    }
  )
)



#' @section Static methods:
#'
#' Static methods for the \code{bspline} class can be use with the
#' syntax : \code{bs_statics$static_method_name}. See [smouseGPS::bs_statics()]
#' for the documentation.
