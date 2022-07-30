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
#' @examples
#' K <- 3
#' D <- K - 1
#' t <- Conj(t(matrix(0:10, nrow = 1)))
#'
#' # Increase the multiplicity of the end knots for higher order splines.
#' t_knot <- rbind(
#'    matrix(rep(t[1], K-1), ncol = 1),
#'    t,
#'    matrix(rep(t[length(t)], K-1))
#' )
#'
#' n_splines <-  length(t) + K - 2
#'
#' m <-  array(0, dim = c(n_splines, 1))
#'
#' B <- bspline$new(t_knot, K, m)
#'
#' m[3] <-  1
#'
#' B$m <- m
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
    B = matrix(),
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

    #' @description Overload the subscript operator `[` for our class.
    #' @param x object
    #' @param i rows
    #' @param j columns
    `[` = function(x, i, j) {

      if(length(i) >= 1) t <- i

      if(length(j) >= 1) {
        num_derivatives <- j
      } else {
        num_derivatives <- 0
      }

      self$value_at_points(t, num_derivatives)
    },

    #' @description Overload the subscript operator `[<-` for our class
    #' @param x object
    #' @param i an integer
    #' @param j an integer
    #' @param value value
    `[<-` = function(x, i, j, value) {

      if(length(i) >= 1) t <- i

      if(length(j >= 1)) {
        num_derivatives <- j
      } else {
        num_derivatives <- 0
      }
      x <- value
      invisible(self)
    },

    #' @description Get the value at points.
    #' @param t a matrix. Points locations.
    #' @param num_derivatives an integer. Number of the derivatives.
    value_at_points = function(t, num_derivatives) {
      if (!exists("num_derivatives")) num_derivatives <- 0

      x_out <- private$evaluate_from_pp_coeff(t, self$C, self$t_pp, num_derivatives)

      if (!isempty(self$x_std)) x_out <- self$x_std * x_out

      if (!isempty(self$x_mean) & num_derivatives == 0) x_out <- x_out + self$x_mean

      return(x_out)
    },

    #' @description Set m
    #' @param new_m a matrix
    set_m = function(new_m) {
      self$m <- new_m
      self$spline_coefficients_did_change()
    },

    #' @description Change matrix if spline coefficients have changed
    spline_coefficients_did_change = function() {
      out_matrix <- cbind(self$C, self$t_pp, self$B)
     out_matrix <-  bs_utils$pp_coeff_from_spline_coeff(
        self$m,
        self$t_knot,
        self$K,
        self$B
      )
    }
  ),
  private = list(
    domain = NULL,
    S = NULL
  ),
  active = list(
    #' @field get_S Get \eqn{S = K - 1}
    get_S = function(self) {
      private$S <- self$K - 1

      return(private$S)
    },

    #' @field get_domain Get the domain of \code{t_knot}
    get_domain = function(self) {
      private$domain <- list(
        start = self$t_knot[1],
        end = self$t_knot[nrow(self$t_knot)]
      )
      attr(private$domain, "class") <- "domain"

      return(private$domain)
    }
  )
)

#' @rdname bspline
#' @details Overloaded subscript operator
`[.bspline` <- function(obj, ...) obj$`[`(...)

#' @rdname bspline
#' @details Overloaded subscript operator assignment.
`[<-.bspline`  <- function(obj, ...) obj$`[<-`(...)

#' @section Static methods:
#'
#' Static methods for the \code{bspline} class can be use with the
#' syntax : \code{bs_statics$static_method_name}. See [smouseGPS::bs_statics()]
#' for the documentation.
