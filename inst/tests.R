library(smouseGPS)


# Test de bspline ---------------------------------------------------------


K <- 4 # order of spline
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
B <- bspline$new(K = K, t_knot = t_knot, m = m)

B
m[3] <-  1

# Change m
B$set_m(m)
B

tq <- Conj(t(matrix(seq(min(t), max(t), length.out = 1000), nrow = 1)))

head(B[tq])

# Plot the spline
plot(tq, B[tq], type = "l")

# Plot the first derivative
plot(tq, B[tq, 1], type = "l")
