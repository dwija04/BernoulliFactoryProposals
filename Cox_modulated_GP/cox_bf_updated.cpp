// cox_bf: Bernoulli factory MCMC for the Cox model.
//
// Proposal step uses propose_step() from cox_proposal_bf.cpp, which draws
// z = current + sqrt(eta) * L * eps and retries until all components are
// positive — replacing the old inline while(!accept) loop.
//
// Everything else (Bernoulli factory loop, auxiliary m1/m2 draws, accept
// bookkeeping) is a faithful translation of the R reference cox_bf().

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// ---------- Forward declarations ----------

// From cox_functions.cpp
double target(const NumericVector &chi_r, const List &x, const NumericVector &c_r,
              const NumericVector &t_knots, const NumericMatrix &cov_r,
              const IntegerVector &n_r, double delta_m);

// From cox_proposal_bf.cpp
int is_positive(const NumericVector &x);
NumericVector propose_step(const NumericVector &current, const arma::mat &L, double eta);

// [[Rcpp::export]]
List cox_bf(int N,
            const NumericVector &init,
            const IntegerVector &n_r,
            const List &x,
            const NumericVector &c_r,
            const NumericVector &t_knots,
            const NumericMatrix &cov_r,
            double eta,
            double delta_m) {

  const int p = init.size();
  NumericMatrix chi(N, p);
  NumericVector logpost(N);
  NumericVector bernoulli_loops(N);

  // L <- chol(cov)  — upper Cholesky in R, lower in Armadillo;
  // propose_step takes L (lower) and computes L * eps, matching R's t(L) * eps
  // because for a symmetric PD matrix: Sigma = L L^T = R^T R, so both give the
  // same distribution.
  arma::mat cov_mat(const_cast<double*>(cov_r.begin()),
                    cov_r.nrow(), cov_r.ncol(), /*copy=*/true);
  arma::mat L;
  if (!arma::chol(L, cov_mat, "lower")) {
    cov_mat += 1e-8 * arma::eye(cov_mat.n_rows, cov_mat.n_cols);
    if (!arma::chol(L, cov_mat, "lower"))
      Rcpp::stop("cox_bf: Cholesky decomposition failed.");
  }

  // Iteration 1 in R == index 0 in C++
  for (int j = 0; j < p; ++j) chi(0, j) = init[j];
  logpost[0] = target(init, x, c_r, t_knots, cov_r, n_r, delta_m);

  int accept_count = 0;
  const int print_stride = std::max(1, N / 10);

  // R loop: for(i in 2:N)
  for (int i = 1; i < N; ++i) {

    if (i % print_stride == 0)
      Rcpp::Rcout << "iter: " << i << "\n";

    // --- Proposal via propose_step (replaces old while(!accept) loop) ---
    // y <- propose_step(chi[i-1, ], L, eta)
    NumericVector curr_r(p);
    for (int j = 0; j < p; ++j) curr_r[j] = chi(i - 1, j);

    NumericVector y = propose_step(curr_r, L, eta);

    // --- Bernoulli factory ---
    // c1 <- target(y, ...);  c2 <- target(chi[i-1,], ...)
    // C  <- exp(c1-c2) / (1 + exp(c1-c2))
    double c1 = target(y,      x, c_r, t_knots, cov_r, n_r, delta_m);
    double c2 = logpost[i - 1];   // always in sync with target(chi[i-1,],...)
    double C  = std::exp(c1 - c2) / (1.0 + std::exp(c1 - c2));

    bool accepted = false;
    int bern_loops = 0;

    while (!accepted) {
      ++bern_loops;

      int C1 = R::rbinom(1, C);

      if (C1 == 1) {
        // m1 <- chi[i-1, ] + sqrt.cov %*% rnorm(p, sd = sqrt(eta))
        // if(is_positive(m1)) p_x <- 1  else p_x <- 0
        arma::vec eps1  = arma::randn<arma::vec>(p);
        arma::vec m1_v  = arma::vec(curr_r.begin(), p) + std::sqrt(eta) * (L * eps1);
        NumericVector m1(p);
        for (int j = 0; j < p; ++j) m1[j] = m1_v[j];

        int p_x = is_positive(m1);   // 1 or 0
        int C2  = R::rbinom(1, double(p_x));

        if (C2 == 1) {
          for (int j = 0; j < p; ++j) chi(i, j) = y[j];
          logpost[i] = c1;
          accepted = true;
        }

      } else {
        // m2 <- y + sqrt.cov %*% rnorm(p, sd = sqrt(eta))
        // if(is_positive(m2)) p_y <- 1  else p_y <- 0
        arma::vec eps2  = arma::randn<arma::vec>(p);
        arma::vec m2_v  = arma::vec(y.begin(), p) + std::sqrt(eta) * (L * eps2);
        NumericVector m2(p);
        for (int j = 0; j < p; ++j) m2[j] = m2_v[j];

        int p_y = is_positive(m2);   // 1 or 0
        int C2  = R::rbinom(1, double(p_y));

        if (C2 == 1) {
          for (int j = 0; j < p; ++j) chi(i, j) = chi(i - 1, j);
          logpost[i] = c2;
          accepted = true;
        }
      }
    }

    bernoulli_loops[i] = bern_loops;

    // if(chi[i, 1] == y[1]) accept_rate <- accept_rate + 1
    if (chi(i, 0) == y[0]) ++accept_count;
  }

  // accept_rate <- accept_rate / N  (matches R)
  double accept_rate = double(accept_count) / double(N);

  return List::create(Named("chi")             = chi,
                      Named("bernoulli_loops") = bernoulli_loops,
                      Named("accept_rate")     = accept_rate,
                      Named("logpost")         = logpost);
}
