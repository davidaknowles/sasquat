##
# Module  which implements the fast update of the cholesky factor after a rank 1 update.
#
# Author: Jose Miguel Hernandez Lobato and Daniel Hernandez Lobato
# Date: June 2008
#
# Reference:
#
# P.E. Gill, G. H. Golub, W. Murray, M. A. Saunders. Methods for Modifying Matrix Factorizations.
# Mathematics of Computation. Vol. 28, No. 126, (Apr., 1976). pp. 505-535
#

dyn.load("cholUpdate.so")

##
# Let A = L %*% t(L) and /A = A + alpha * z %*% t(z) = L %*% /L %*% t(/L) %*% t(L).
# The function computes the updated cholesky factor L %*% /L.
# The function is based on the equality (/L)^-1 = I - (alpha / (1 + alpha * t(p) %*% p)) * p %*% t(p).
#
# Author: Jose Miguel Hernandez Lobato and Daniel Hernandez Lobato
# Date: June 2008
#
# Input:
# 	- alpha -> a constant.
# 	- z	-> a vector.
# 	- L	-> the previous cholesky factor.
#
# Output:
#	- A matrix with the updated cholesky factor.
#
# Reference:
#
# P.E. Gill, G. H. Golub, W. Murray, M. A. Saunders. Methods for Modifying Matrix Factorizations.
# Mathematics of Computation. Vol. 28, No. 126, (Apr., 1976). pp. 505-535
#

rankOneUpdate <- function(alpha, z, L) {

	# We call the c function which performs the updates

	ret <- .C("rankOneUpdate", z = as.double(z), p = as.double(z), alpha = as.double(alpha),
		rows = as.integer(nrow(L)), beta = as.double(z), dhat = as.double(z), L = as.double(L))

	matrix(ret$L, nrow(L), ncol(L))
}
