#' @title Utils for the \code{bspline} class
#'
#' @description A list of utils for the [smouseGPS::bspline] class. They can be
#' used with the synthax \code{bs_utils$function_name()}. These functions are
#' equivalent to static methods for the \code{bspline} class.
#'
#' @return \code{bs_utils} is a list containing functions for the [smouseGPS::bspline] class.
#'
#' @importFrom Rfast sort_mat colRanks
#' @importFrom pracma polyval mrdivide isempty numderiv numdiff

#' @describeIn bs_utils This static method assumes that the splines are terminated at the boundary with repeat knot points.
#' @param t_knot a matrix. Spline knot points.
#' @param K an integer. Order of polynomial.
#' @param D an array.
point_of_support <-  function(t_knot, K, D) {
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
}

#' @describeIn bs_utils Returns the piecewise polynomial coefficients in matrix C from spline coefficients in vector m.
#' @param m a numeric. Spline coefficients (\eqn{m = M \cdot D}).
#' @param t_knot a matrix. Spline knot points.
#' @param K an integer. Order of polynomial.
#' @param B an array.
pp_coeff_from_spline_coeff <-  function(m, t_knot, K, B) {
  Nk <- nrow(t_knot)
  t_pp <- t_knot[K:(Nk-K+1)]
  if(!exists("B") | isempty(B)) B <- spline(t_pp, t_knot, K, K-1)
  C <- matrix(0, nrow = nrow(t_pp-1), ncol = K)
  for(i in 1:K) {
    C[,K-i+1] <- B[1:nrow(B)-1, , i] * m
  }
  return(c(C, t_pp, B))
}

#' @describeIn bs_utils Returns the value of the function with derivative \code{D} represented by \code{PP} coefficients \code{C} at locations \code{t}. \code{t_pp} containes the intervals.
#' @param t a matrix. Point locations.
#' @param C a vector. Piecewise polynomial coefficients (\code{size(C) = [nrow(t_pp)-1, K]})
#' @param t_pp a matrix. pp break points (\code{size(t_pp) = nrow(t_knot) - 2*K + 1}).
#' @param D an array.
evaluate_from_pp_coeff <-  function(t, C, t_pp, D) {
  if(!is.unsorted(t[,1]) & t[1, 1] <= t[nrow(t), 1]) { # Is ascending ord
    did_flip <- 0
  } else if(!is.unsorted(t[,1]) & t[1, 1] >= t[nrow(t), 1]) { # Is descen
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
    return(f) # # By construction the splines are zero for K or more deri
  scale <- factorial(K-1-D)
  i <- 1
  while(scale > 0) {
    scale[i + 1] <- factorial(scale[i] - 1)
  }
  indices <- 1:(K-D)
  t_pp_bin <- findInterval(t, c(-Inf, t_pp[2:(nrow(t_pp)-1)], Inf))
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
}


#' @describeIn bs_utils Returns the basis splines of order K evaluated at point t, given knot points t_knot. If you optionally provide D, then D derivatives will be returned. \code{dim(B) = c(nrow(t), M, D)} where \code{M = nrow(t_knot) - K} is the number of splines.
#' @param t a matrix. Point locations.
#' @param t_knot a matrix. Spline knot points.
#' @param K an integer. Order of polynomial.
#' @param D an array.
spline <-  function(t, t_knot, K, D) {
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
  nl <- nl[length(nl)]
  nr <- which(t_knot == t_knot[length(t_knot)])
  nr <- M - nr[1] + 1
  if(nl < K | nr < K) {
    stop(paste(
      "Your splines are not terminated. You need to have K repeat ",
      "knot points at the beginning and end."
    ))
  }
  N_splines = M - K
  N = nrow(t)
  B <- array(0, dim = c(N, N_splines, D+1))
  delta_r <- matrix(0, nrow = N, ncol = K)
  delta_l <- matrix(0, nrow = N, ncol = K)
  knot_indices <- findInterval(t, t_knot[1:(M-K+1)])

  # XB will contain all splines from (K-D) through order (K-1).
  # These are needed to compute the derivatives of the spline, if
  # requested by the user.
  #
  # The indexing is such that the spline of order m, is located
  # at index m-(K-D)+1

  if(D > 0) {
    XB <- array(0, dim = c(N, N_splines, D))

    # if we go tho the max derivative, we need to manually create the 0th order spline.
    if(D+1 == K) {
      to_indices <- (
        Conj(t(matrix(1:N, nrow = 1))) +
          dim(XB)[1] *
          ((knot_indices - 1) + dim(XB)[2] * (1-1))
      )

      XB[to_indices] <- 1
    }
  }

  b <- array(0, dim = c(N, K))
  b[,1] <- 1

  for(j in 1:(K-1)) { # Loop through splines of increasing order j+1.
    delta_r[,j] <- t_knot[knot_indices + j] - t
    delta_l[,j] <- t - t_knot[knot_indices + 1 - j]

    saved <- array(0, dim = c(N, 1))

    for(r in 1:j) { # Loop through the nonzero splines.
      term <- mrdivide(b[,r], delta_r[,r] + delta_l[,j+1-r])
      b[,r] <- saved + delta_r[,r] * term
      saved <- delta_l[,j+1-r] * term
    }

    b[,j+1] <- saved

    # Save this info for later use in computing the derivatives have to loop
    # through one index.
    if(j+1 >= K-D & j+1 <= K-1) {
      for(r in 1:j+1) {
        to_indices <- (
          Conj(t(matrix(1:N, nrow = 1))) +
            dim(XB)[1] *
            ((knot_indices - j+(r-1)-1) + dim(XB)[2] * (j+1-(K-1-D) - 1))
        )

        from_indices <- Conj(t(matrix(1:N, nrow = 1))) + dim(b)[1] * (r-1)

        XB[to_indices] <- b[from_indices]
      }
    }
  }

  for(r in 1:K) {
    to_indices <- (
      Conj(t(matrix(1:N, nrow = 1))) +
        dim(XB)[1] *
        ((knot_indices - (K-1) + (r-1)-1) + dim(B)[2] * (1-1))
    )

    from_indices <- Conj(t(matrix(1:N, nrow = 1))) + dim(b)[1] * (r-1)

    B[to_indices] <- b[from_indices]
  }

  diff_coeff <- function(a, r, m) (K-m) * (a[2]-a[1]) / (t_knot[r+K-m] - t_knot[r])

  if(D > 0) {
    for(r in 1:N_splines) {
      alpha <- array(0, dim = c(K+1, K+1))
      alpha[2,1] <- 1

      for(m in 1:D) { # loop over derivatives
        for(i in 1:(m+1)) { # Loop over coefficients
          a <- alpha[,m]
          alpha[i+1,m+1] <- diff_coeff(a[i:length(a)], r+i-1, m)

          if(is.infinite(alpha[i+1, m+1]) | is.nan(alpha[i+1, m+1])) {
            alpha[i+1, m+1] <- 0
          }

          if(r+i-1 > N_splines) {
            B0 <- array(0, dim = c(N, 1))
          } else {
            B0 <- XB[, r+i-1, D+1-m] # Want the k-m order spline, in position D+1-m
          }

          B[, r, m+1] <- B[, r, m+1] + alpha[i+1,m+1] * B0
        }
      }
    }
  }

  return(B)
}

#' @rdname bs_utils
#' @export
bs_utils <- structure(
  list(
    point_of_support = point_of_support,
    pp_coeff_from_spline_coeff = pp_coeff_from_spline_coeff,
    evaluate_from_pp_coeff = evaluate_from_pp_coeff,
    spline = spline
  ),
  class = "bspline_static_methods"
)
