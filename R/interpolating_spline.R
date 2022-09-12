#' @title R6 class that interpolates splines.
#'
#' @description This R6 class inherits from the \code{bspline} class.
#'
#' @return Fitted data
#' @importFrom R6 R6Class

#' @rdname interpolating_spline
#' @field t array of values for the independent axis.
#' @field x array of values for the dependent axis.
#' @field K Order of the spline.
#' @field fun Cubic spline interpolant.
#' @export
smoothing_spline <- R6Class(
  classname = "interpolating_spline",
  inherit = bspline,
  public = list(
    #' @description Instantiate an object of the class \code{interpolating_spline}.
    initialize = function(t, x, distribution, ...) {
      args <- list(...)
    }
  ),
  private = list(

  )
)
