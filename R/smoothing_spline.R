#' @title R6 class that fits noisy data with a tensioned interpolating spline.
#'
#' @description This R6 class inherits from the \code{bspline} class.
#'
#' @return Fitted data
#' @importFrom R6 R6Class

#' @rdname smoothing_spline
#' @field t array of values for the independent axis.
#' @field x array of values for the dependent axis.
#' @field distribution distribution of the noise.
#' @field ... Other arguments, see details.
#'
#' @field T degree at which tension is applied.
#' @field lambda smoothing parameter, either pass a numeric value, or the
#' lambda enumeration. Default is lambda$optimal_iterated.
#' @field is_constrained indicates wether or not lambda was so big that the
#' solution is just a constrained solution.
#' @field mu mean value of the tension variable.
#' @field knot_dof knot dofs
#' @field covariance computed from the given distribution, this is the
#' covariance structure of the observations. It may be a scalar, vector, or
#' matrix.
#' @field variable_cache structure storing several cached variables, useful for
#' quick tension spline computation.
#' @field did_override_sigma
#' @field sigma initial weight (given as normal standard deviation).
#' @field constraints constraints =
#' @field outlier_distribution
#' @field alpha
#' @field lambda_at_full_tension what was the value of lambda at full tension.
#' @field sigma_at_full_tension what was the set of 'sigmas' produced from the
#' full tension solution.
#' @field outlier_indices
#' @field outlier_threshold set to a distance with < 1/10000 odds.
#'
#' @export
smoothing_spline <- R6Class(
  classname = "smoothing_spline",
  inherit = bspline,
  public = list(
    t = NULL,
    x = NULL,
    distribution = NULL,
    T = NULL,
    lambda = NULL,
    is_constrained = NULL,
    mu = NULL,
    knot_dof = NULL,
    covariance = NULL,
    variable_cache = NULL,
    did_override_sigma = NULL,
    sigma = NULL,
    constraints = NULL,
    outlier_distribution = NULL,
    alpha = NULL,
    lambda_at_full_tension = NULL,
    sigma_at_full_tension = NULL,
    outlier_indices = NULL,
    outlier_threshold = NULL,

    #' @description Instantiate an object of the class \code{smoothing_spline}.
    initialize = function(t, x, distribution, ...) {
      args <- list(...)

      for (i in 1:length(args)) {
        if (names(args)[i] == "K") {
          K = args[[i]]
        } else if (names(args)[i] == "K") {

        } else if (names(args)[i] == "K") {

        } else if (names(args)[i] == "K") {

        } else if (names(args)[i] == "K") {

        } else if (names(args)[i] == "K") {

        } else if (names(args)[i] == "K") {

        } else if (names(args)[i] == "K") {

        } else if (names(args)[i] == "K") {

        }
      }

      self$t = t
      self$x = x
      self$distribution = distribution
    }
  ),
  private = list(
    X = NULL,
    W = NULL,
    XWX = NULL,
    Cm = NULL,
    Cm_inv = NULL,
    non_outlier_indices = NULL,
    tension_value = NULL
  )
)


