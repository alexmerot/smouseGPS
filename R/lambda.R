#' @title Get lambda methods.
#'
#' @description Get a list of available lambda methods.
#'
#' @return A character which represents the lamba methods.
#'
#' @export

lambda <- structure(
  list(
    cross_validation = "cross_validation",
    optimal_iterated = "optimal_iterated",
    optimal_expected = "optimal_expected",
    full_tension_iterated = "full_tension_iterated",
    full_tension_expected = "full_tension_expected",
    optimal_ranged_iterated = "optimal_ranged_iterated",
    optimal_expected_robust = "optimal_expected_robust",
    full_tension_expected_robust = "full_tension_expected_robust"
  ),
  class = "lambda_methods"
)

