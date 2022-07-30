library(smouseGPS)


# Test de bspline ---------------------------------------------------------


K <- 3 # order of spline
D <- K-1 # number of derivates to return
t <- Conj(t(matrix(0:10, nrow = 1))) # observation points

# increase the multiplicity of the end knots for higher order splines
t_knot <- rbind(
  array(t[1], dim = c(K-1, 1)),
  t,
  array(t[length(t)], dim = c(K-1, 1))
)

n_splines <- nrow(t) + K-2

# coefficients for the bsplines---set all of them to zero for now.
m <- array(0, dim = c(n_splines,1))

# initialize the BSpline class
B <- bspline$new(K,t_knot,m)

m[3] <-  1
B$m <-  m

B[t]
B
