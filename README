# Efficient Gibbs for Bayesian spike-and-spike regression

This code is adapted from code written by José Miguel Hernández-Lobato (https://jmhl.org/). The idea contributions over a naive implementation are
 - collapsing out the continuous component
 - rank-1 updates to the covariance of the conditional distribution. 

## Compilation

Run

R CMD SHLIB *.c *.f

to compile the C and fortran code. If you're having problems on Mac you probably need an up-to-date gfortran see e.g. 
https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks--lgfortran-and--lquadmath-error/

## Usage

See `testSampling.R` for a complete example with synthetic data. 

At some point I'll make this into an actual package! 