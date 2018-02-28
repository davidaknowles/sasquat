#include <R.h>
#include <Rmath.h>

void backsolve(double *p, int *rows, int *columns, double *L, double *z) {

	int i, j, n1, n2;
	double aux1, aux2;
	
	/* We get the number of rows and columns in the matrix M. The number of rows in M is equal to the dimension of /L */

	n1 = *rows;
	n2 = *columns;

	/* We find the vector p such that L %*% p = z by back substitution. See the article in Wikipedia about the Triangular Matrix. */
	
	p[ 0 ] = z[ 0 ] / L[ 0 ];
	for (i = 1 ; i < n1 ; i++) {
		aux1 = 0;
		for (j = 0 ; j < i ; j++)
			aux1 = aux1 + p[ j ] * L[ i + j * n1 ];

		p[ i ] = (z[ i ] - aux1) / L[ i + i * n1 ];
	}
}

