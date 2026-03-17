// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
int is_positive(const NumericVector &x) {
  int n = x.size();
  for (int i = 0; i < n; i++) {
    if (x[i] <= 0) return 0;  // not positive
  }
  return 1;  // all positive
}

// [[Rcpp::export]]
NumericVector propose_step(const NumericVector &current,
                           const arma::mat &L,
                           double eta) {
  int p = current.size();
  NumericVector proposal(p);
  arma::vec base(p);

  for (int j = 0; j < p; ++j) base[j] = current[j];

  while (true) {
    arma::vec eps = arma::randn<arma::vec>(p);   // N(0,I)
    arma::vec prop = base + std::sqrt(eta) * (L * eps);
    for (int j = 0; j < p; ++j) proposal[j] = prop[j];

    if (is_positive(proposal)) {
      return proposal; // valid proposal found
    }
  }
}

// current: current state vector
// Sigma: covariance matrix (symmetric, positive-definite)
// eta: step-size


