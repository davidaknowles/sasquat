/**
 * Module  which implements the fast update of the cholesky factor after a rank 1 update.
 *
 * Author: Jose Miguel Hernandez Lobato and Daniel Hernandez Lobato
 * Date: June 2008
 *
 * Reference:
 *
 * P.E. Gill, G. H. Golub, W. Murray, M. A. Saunders. Methods for Modifying Matrix Factorizations.
 * Mathematics of Computation. Vol. 28, No. 126, (Apr., 1976). pp. 505-535
 *
 */

#include <R.h>
#include <Rmath.h>

/**
 * Let A = L %*% t(L) and /A = A + alpha * z %*% t(z) = L %*% /L %*% t(/L) %*% t(L).
 * This function computes the updated cholesky factor /L.
 * The function is based on the equality (/L)^-1^t(/L)^-1 = I - (alpha / (1 + alpha * t(p) %*% p)) * p %*% t(p).
 *
 * Author: Jose Miguel Hernandez Lobato and Daniel Hernandez Lobato
 * Date: June 2008
 *
 * Input:
 * 	- p     -> vector where to store the solution of L %*% p = z.
 * 	- z     -> a vector.
 * 	- alpha -> a constant.
 * 	- rows  -> number of rows in the matrix L.
 * 	- beta  -> a vector where to store the beta coefficients.
 * 	- dhat  -> vector where to store the square root of the dhat coefficients.
 * 	- L	-> current cholesky factor.
 *
 *
 * Reference:
 *
 * P.E. Gill, G. H. Golub, W. Murray, M. A. Saunders. Methods for Modifying Matrix Factorizations.
 * Mathematics of Computation. Vol. 28, No. 126, (Apr., 1976). pp. 505-535
 *
 */

void rankOneUpdate(double *z, double *p, double *alpha, int *rows, double *beta, double *dhat, double *L) {

	int i, j, n1;
	double aux1, aux2;
	
	/* We get the number of rows */

	n1 = *rows;

	/* We find the vector p such that L %*% p = z by back substitution. See the article in Wikipedia about the Triangular Matrix. */
	
	p[ 0 ] = z[ 0 ] / L[ 0 ];
	for (i = 1 ; i < n1 ; i++) {
		aux1 = 0;
		for (j = 0 ; j < i ; j++)
			aux1 = aux1 + p[ j ] * L[ i + j * n1 ];

		p[ i ] = (z[ i ] - aux1) / L[ i + i * n1 ];
	}

	/* We compute the recurrence relations for the betas and the dhats */

	aux1 = *alpha;
	for (i = 0 ; i < n1 ; i++) {
		dhat[ i ] = 1 + aux1 * p[ i ] * p[ i ];
		beta[ i ] = aux1 * p[ i ] / dhat[ i ];
		aux1 = aux1 / dhat[ i ];
		dhat[ i ] = sqrt(dhat[ i ]);
	}

	/* We update the Cholesky factor */

	for (i = 0 ; i < n1 ; i++) {
		aux1 = 0;
		aux2 = 0;
		for (j = i ; j >= 0 ; j--) {
			aux1 = aux1 + aux2;
			aux2 = L[ i + j * n1 ] * p[ j ];
			L[ i + j * n1 ] = (L[ i + j * n1 ] + aux1 * beta[ j ]) * dhat[ j ];
		}
	}
}
