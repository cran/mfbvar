#ifndef _PROGUTILS_H_
#define _PROGUTILS_H_

/* Contains the following code modules:

   a) some helper functions such as progress bar tools and return value
      prettifier

   b) some functions related to the Cholesky decomposition used for
      sampling AWOL and efficiently solving the systems of linear
      equations

   c) function for inverse transform sampling

   d) a very basic Newton-Raphson algorithm for finding the root
      of dlogdnu (defined in densities.h)
*/

#include <RcppArmadillo.h>

// b)
// Cholesky factor for a tridiagonal matrix with constant off-diagonal
void cholTridiag(
    const arma::vec& omega_diag,
    double omega_offdiag,
    arma::vec& chol_diag,
    arma::vec& chol_offdiag);

// Solves Chol*x = covector ("forward algorithm")
void forwardAlg(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& covector,
    arma::vec& htmp);

// Solves (Chol')*x = htmp ("backward algorithm")
void backwardAlg(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& htmp,
    arma::vec& h);

// c)
// draws length(r) RVs, expects the non-normalized CDF mixprob
void invTransformSampling(
    const arma::vec& mixprob,
    arma::ivec& r,
    int T);

// d)
// truncated normal (stationary)
double rtruncnorm(double m, double v);
#endif
