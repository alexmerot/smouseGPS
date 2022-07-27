#' @title Static methods for the \code{bspline} class
#'
#' @description Static methods for the [smouseGPS::bspline] class can be use with the
#' syntax : \code{bs_static$static_method_name()}.
#'
#' @return A list that contains all the static methods for the [smouseGPS::bspline] class.
#'
#' @importFrom R6 R6Class
#' @importFrom Rfast sort_mat colRanks
#' @importFrom pracma polyval mrdivide isempty numderiv
#' @export

# @field m a numeric. Spline coefficients (\eqn{m = M \cdot D}).
# @field B a numeric.
# @field x_mean a numeric. If set, these will be used to scale the output
# (default to 0).
# @field x_std a numeric. \eqn{x_{out} = x_{std} \cdot (X \cdot m) + x_{mean}}
# @field t_pp a matrix. pp break points (\code{size(t_pp) = nrow(t_knot) - 2*K + 1}).
# @field C a vector. Piecewise polynomial coefficients
# (\code{size(C) = [nrow(t_pp)-1, K]}).

#' @rdname bs_static
bs_static <- structure(
  list(
    #' @field point_of_support() This static method assumes that the splines are
    #' terminated at the boundary with repeat knot points.
    #'
    #' Arguments:
    #' * \code{t_knot} a Matrix. Spline knot points.
    #' * \code{K} an integer. Order of polynomial.
    #' * \code{D} an array.
    point_of_support = function(t_knot, K, D) {
      interior_knots <- t_knot[(K+1):(nrow(t_knot)-K)]

      if(isempty(interior_knots)) {
        if (K == 1) {
          t <- t_knot
        } else {
          dt = (t_knot[nrow(t_knot)] - t_knot[1]) / (K - 1)
          t <- t_knot[1] + dt * (0:(K-1))
        }

        return(t)
      }

      if(K %% 2 == 1) {
        interior_support <- interior_knots[1:(nrow(interior_knots)-1)] +
          numderiv(interior_knots) /
          2

        n <- K/2

        dt_start <- (interior_knots[1] - t_knot[1])/n
        dt_end <- (t_knot[nrow(t_knot)] - interior_knots[nrow(interior_knots)])/n

        n <- ceiling(n)

        t <- rbind(
          t_knot[1] + dt_start * Conj(t(matrix(0:(n-1), nrow = 1))),
          interior_support,
          t_knot[nrow(t_knot)] - dt_end * Conj(t(matrix(seq(n - 1, 0, -1), nrow = 1)))
        )
      } else {
        interior_support <- interior_knots

        n <- floor((K+1)/2)

        dt_start <- (interior_knots[1]-t_knot[1]) / n
        dt_end <- (t_knot[nrow(t_knot) - interior_knots[nrow(interior_knots)]]) / n

        t <- rbind(
          t_knot[1] + dt_start * Conj(t(matrix(0:(n-1), nrow = 1))),
          interior_support,
          t_knot[nrow(t_knot)] - dt_end * Conj(t(matrix(seq(n - 1, 0, -1), nrow = 1)))
        )
      }

      return(t)
    },

    #' @field pp_coeff_from_spline_coeff() Returns the piecewise polynomial coefficients in matrix C
    #' from spline coefficients in vector m.
    #'
    #' Arguments:
    #' * \code{m} Spline coefficients (\eqn{m = M \cdot D}).
    #' * \code{t_knot} a two-column matrix. Spline knot points.
    #' * \code{K} an integer. Order of polynomial.
    #' * \code{B} An array. Optional
    pp_coeff_from_spline_coeff = function(m, t_knot, K, B) {
      Nk <- nrow(t_knot)
      t_pp <- t_knot[K:(Nk-K+1)]

      if(!exists("B") | isempty(B)) B <- private$spline(t_pp, t_knot, K, K-1)

      C <- matrix(0, nrow = nrow(t_pp-1), ncol = K)

      for(i in 1:K) {
        C[,K-i+1] <- B[1:nrow(B)-1, , i] * m
      }

      return(c(C, t_pp, B))
    },

    # @rdname bs_static
    #' @field evaluate_from_pp_coeff() Returns the value of the function with derivative \code{D}
    #' represented by \code{PP} coefficients \code{C} at locations \code{t}.
    #' \code{t_pp} containes the intervals.
    #'
    #' Arguments:
    #' * \code{t} a matrix. Points locations.
    #' * \code{C} a vector. Piecewise polynomial coefficients.
    #' * \code{t_pp} a matrix. pp break points. (\code{size(t_pp) = nrow(t_knot) - 2*K + 1}).
    #' * \code{D}
    evaluate_from_pp_coeff = function(t, C, t_pp, D) {
      if(!is.unsorted(t[,1]) & t[1, 1] <= t[nrow(t), 1]) { # Is ascending order
        did_flip <- 0

      } else if(!is.unsorted(t[,1]) & t[1, 1] >= t[nrow(t), 1]) { # Is descending order
        t <- t[nrow(t):1,] # Flip the order
        did_flip <- 1

      } else {
        warning("t was not sorted.")
        return_indices <- colRanks(t, method = "first", stable = TRUE)
        t <- sort_mat[t]
        did_flip <- 2
      }

      K <-  dim(C)[2]
      f <- matrix(0, nrow = dim(t)[1], ncol = dim(t)[2])

      if(nargs() < 4) {
        D <-  0
      } else if (D > K-1)
        return(f) # # By construction the splines are zero for K or more derivs

      scale <- factorial(K-1-D)
      i <- 1
      while(scale > 0) {
        scale[i + 1] <- factorial(scale[i] - 1)
      }

      indices <- 1:(K-D)

      t_pp_bin <- findInterval(t, c(-Inf, tpp[2:(nrow(tpp)-1)], Inf))
      start_index <- 1

      while(start_index <= nrow(t)) {
        i_bin <- t_pp_bin(start_index)
        index <- which(t_pp_bin == i_bin)
        end_index <- index[nrow(index)]
        f[start_index:end_index] <- polyval(
          mrdivide(C[i_bin, indices], scale),
          t[start_index:end_index] - t_pp[i_bin]
        )

        start_index <- end_index + 1
      }

      # Include an extrapolated points past the end.
      if(start_index <= nrow(t)) {
        f[start_index:nrow(f)] <- polyval(
          mrdivide(C[i_bin, indices], scale),
          t[start_index:nrow(t)] - t_pp[i_bin]
        )
      }

      if(did_flip == 0){
        return(f)
      } else if(did_flip == 1) {
        f <- f[nrow(f):1,]
      } else {
        f <- f[return_indices]
      }

      return(f)
    },

    #' @field spline() Returns the basis splines of order K evaluated at point t,
    #' given knot points t_knot. If you optionally provide D, then D derivatives
    #' will be returned.
    #' \code{dim(B) = c(nrow(t), M, D)} where \code{M = nrow(t_knot) - K} is the
    #' number of splines.
    #'
    #' Arguments:
    #' * \code{t} a matrix. Points locations.
    #' * \code{t_knot} a Matrix. Spline knot points.
    #' * \code{K} an integer. Order of polynomial.
    #' * \code{D}
    spline = function(t, t_knot, K, D) {
      if(any(numdiff(t_knot) < 0)) {
        stop("t_knot must be non-decreasing.")
      }

      if(nargs() < 4) {
        D <- 0
      } else if(D > K-1) {
        D <- K-1
      }

      # Number of knots
      M <- nrow(t_knot)

      nl <- which(t_knot <= t_knot[1])
      nl <- n1[length(N1)]



      return(B)
    }
  ),
  class = "bspline_static_methods"
)
