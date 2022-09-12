#' @title R6 class that fits noisy data with a tensioned interpolating spline.
#'
#' @description This R6 class inherits from the \code{bspline} class.
#'
#' @return Fitted data
#' @importFrom R6 R6Class

#' @rdname smoothing_spline
#' @field t array of values for the independent axis.
#' @field x array of values for the dependent axis.
#' @field distribution distribution of the noise
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
    Non_outlier_indices = NULL,
    tension_value = NULL
  )
)


