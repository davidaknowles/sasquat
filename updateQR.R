# Module: updateQR.R
#
# Module that updates the R  part of a QR decomposition of a matrix
# B when a column is appended to the right of B of when a column
# of B is removed.
# The functions in this module are based on the article on
#
#	Methods for Modifying Matrix Factorizations.
#	P. E. Gill; G. H. Golub; W. Murray; M. A. Saunders 
#	Mathematics of Computation, Vol. 28, No. 126. (Apr., 1974), pp. 505-535. 
# 
# and on the rutine "dchud" of the linpack library.
#

# We are going to use the linpack library...

dyn.load("dchud.so")
dyn.load("backsolve.so")

##
# Updates the R part of a QR decomposition of a matrix B
# when a column is appended to the right of B.
#
# R -> Old R part of the QR decomposition
# B -> Matrix to which a column is appended
# b -> Column appended to the right of B
#
# Returns: updated R part of the QR decomposition
#

addColumn <- function(R, B, b, alpha) {

	L <- t.default(R)
	v <- t.default(B) %*% b
	p <- v
	
	u <- .C("backsolve", p = p, nrow(L), ncol(L), L, v)$p
	gamma <- sqrt(sum(b^2) - sum(u^2) + alpha)

	rbind(cbind(R, u), c(rep(0, ncol(R)), gamma))
}

##
# Updates the R part of a QR decomposition of a matrix B
# when a column in B is deleted.
#
# R -> Old R part of the QR decomposition
# col -> Index of the column to be deleted
#
# Returns: updated R part of the QR decomposition
#

delColumn <- function(R, col) {

	if (col == ncol(R)) {

		# We do not have to update R when we take out the last column of B
	
		R <- R[ 1 : (ncol(R) - 1), 1 : (nrow(R) - 1) ]
	} else {
		T2 <- as.matrix(R[ col : nrow(R), (col + 1) : ncol(R) ])

		row <- as.matrix(T2[ 1, ])
		Raux <- as.matrix(T2[ 2 : nrow(T2), ])
	
		# We call the rutine dchud of linpack to update Raux when the row "row" is appended to the matrix
		# to which Raux is the R part of its QR decomposition.
	
		retornoDchud <- .Fortran("dchud", r = Raux, ldr = nrow(Raux), p = nrow(Raux), x = row,
			z = Raux, ldz = nrow(Raux), nz = 1, y = 0,
			rho = 1, c = rep(0, nrow(Raux)), s = rep(0, nrow(Raux)))
	
		Raux <- retornoDchud$r
		if (col == 1) {
			R <- as.matrix(R[ 1 : (ncol(R) - 1), 2 : ncol(R) ])
		} else {
			R <- as.matrix(R[ 1 : (ncol(R) - 1), c(1 : (col - 1), (col + 1) : ncol(R)) ])
		}
		R[ col : nrow(R) , col : ncol(R) ] <- Raux
	}

	as.matrix(R)
}
