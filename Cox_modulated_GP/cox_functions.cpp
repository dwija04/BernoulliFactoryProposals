// cox_mcmc_rcpp.cpp
// Rcpp/Armadillo translation of the provided R MCMC code (cox_bf, cox_mh, cox_rwmh, helpers)
// Compile with: Rcpp::sourceCpp('cox_mcmc_rcpp.cpp')

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// ---------- Helpers ----------



// phi: triangular kernel. del: scalar, x: scalar, t: numeric vector
// [[Rcpp::export]]
NumericVector phi(double del, double x, const NumericVector &t) {
  int p = t.size();
  NumericVector out(p);
  for (int j = 0; j < p; ++j) {
    double val = std::abs((x - t[j]) / del);
    out[j] = (val <= 1.0) ? (1.0 - val) : 0.0;
  }
  return out;
}

// is_positive: all entries > 0
// [[Rcpp::export]]
int is_positive(const NumericVector &v) {
  for (int i = 0; i < v.size(); ++i) if (v[i] <= 0.0) return 0;
  return 1;
}

// approx_trunc_prob: Monte Carlo estimate P(X >= 0 elementwise) for X ~ N(mu, cov)
// [[Rcpp::export]]
double approx_trunc_prob(const NumericVector &mu_r, const NumericMatrix &cov_r, int B = 200) {
  arma::vec mu(mu_r.begin(), mu_r.size());
  arma::mat cov(const_cast<double*>(cov_r.begin()), cov_r.nrow(), cov_r.ncol(), true);


  arma::mat L;
  bool ok = chol(L, cov, "lower");
  if (!ok) {
    cov += 1e-8 * arma::eye(cov.n_rows, cov.n_cols);
    chol(L, cov, "lower");
  }
  int count = 0;
  for (int b = 0; b < B; ++b) {
    arma::vec eps = arma::randn<arma::vec>(mu.n_elem);
    arma::vec samp = mu + L * eps;
    bool allpos = true;
    for (size_t j = 0; j < samp.n_elem; ++j) if (samp[j] < 0.0) { allpos = false; break; }
    if (allpos) ++count;
  }
  return double(count) / double(B);
}

// target: log-posterior
// [[Rcpp::export]]
double target(const NumericVector &chi_r, const List &x, const NumericVector &c_r,
              const NumericVector &t_knots, const NumericMatrix &cov_r, const IntegerVector &n_r,
              double delta_m) {
  int p = chi_r.size();
  for (int i = 0; i < p; ++i) if (chi_r[i] <= 0.0) return R_NegInf;

  arma::vec chi(chi_r.begin(), chi_r.size());
  arma::mat cov(const_cast<double*>(cov_r.begin()), cov_r.nrow(), cov_r.ncol(), true);

  arma::mat invcov = arma::inv_sympd(cov);

  int n0 = n_r.size();
  double quad = -0.5 * as_scalar(chi.t() * invcov * chi);
  arma::vec cvec(c_r.begin(), c_r.size());
  double linear = - double(n0) * as_scalar(cvec.t() * chi);
  double ret = quad + linear;

  for (int j = 0; j < n0; ++j) {
    int nj = n_r[j];
    NumericVector temp = x[j];
    for (int i = 0; i < nj; ++i) {
      NumericVector phi_v = phi(delta_m, temp[i], t_knots);
      double dot = 0.0;
      for (int k = 0; k < p; ++k) dot += phi_v[k] * chi_r[k];
      if (dot <= 0.0) return R_NegInf;
      ret += std::log(dot);
    }
  }
  return ret;
}

// smooth: t(phi) %*% samp
// [[Rcpp::export]]
double smooth(double delta_m, const NumericVector &t_knots, const NumericVector &samp, double val) {
  NumericVector phi_v = phi(delta_m, val, t_knots);
  double out = 0.0;
  for (int i = 0; i < samp.size(); ++i) out += phi_v[i] * samp[i];
  return out;
}

// ---------- MCMC samplers ----------

// cox_rwmh
// [[Rcpp::export]]
List cox_rwmh(int N, const NumericVector &init, const IntegerVector &n_r,
              const List &x, const NumericVector &c_r, const NumericVector &t_knots,
              const NumericMatrix &cov_r, double eta, double delta_m) {
  int p = init.size();
  NumericMatrix chi(N, p);
  NumericVector logpost(N);

  arma::mat cov(const_cast<double*>(cov_r.begin()), cov_r.nrow(), cov_r.ncol(), true);
  arma::mat L;
  bool ok = chol(L, cov, "lower");
  if (!ok) {
    cov += 1e-8 * arma::eye(cov.n_rows, cov.n_cols);
    chol(L, cov, "lower");
  }

  for (int j = 0; j < p; ++j) chi(0, j) = init[j];
  double lp0 = target(init, x, c_r, t_knots, cov_r, n_r, delta_m);
  logpost[0] = lp0;
  int accept_count = 0;

  for (int i = 1; i < N; ++i) {
    if (i % (N/10) == 0) Rcout << "iter: " << i << '\n';

    arma::vec eps = arma::randn<arma::vec>(p);
    arma::vec curr(p);
    for (int j = 0; j < p; ++j) curr[j] = chi(i-1, j);
    arma::vec prop = curr + std::sqrt(eta) * (L * eps);
    NumericVector z_r(p);
    for (int j = 0; j < p; ++j) z_r[j] = prop[j];

    bool okpos = true;
    for (int j = 0; j < p; ++j) if (z_r[j] <= 0.0) { okpos = false; break; }
    if (!okpos) {
      for (int j = 0; j < p; ++j) chi(i, j) = chi(i-1, j);
      logpost[i] = logpost[i-1];
      continue;
    }

    double num = target(z_r, x, c_r, t_knots, cov_r, n_r, delta_m);
    double denom = logpost[i-1];
    double ratio = std::exp(num - denom);
    if (R::runif(0.0, 1.0) < std::min(1.0, ratio)) {
      for (int j = 0; j < p; ++j) chi(i, j) = z_r[j];
      logpost[i] = num;
      ++accept_count;
    } else {
      for (int j = 0; j < p; ++j) chi(i, j) = chi(i-1, j);
      logpost[i] = denom;
    }
  }
  double accept_rate = double(accept_count) / double(N);
  return List::create(Named("chi") = chi,
                      Named("accept_rate") = accept_rate,
                      Named("logpost") = logpost);
}

// cox_mh
// [[Rcpp::export]]
List cox_mh(int N, const NumericVector &init, const IntegerVector &n_r,
            const List &x, const NumericVector &c_r, const NumericVector &t_knots,
            const NumericMatrix &cov_r, double eta, double delta_m, int B = 200) {
  int p = init.size();
  NumericMatrix chi(N, p);
  NumericVector logpost(N);

  arma::mat cov(const_cast<double*>(cov_r.begin()), cov_r.nrow(), cov_r.ncol(), true);
  arma::mat L;
  bool ok = chol(L, cov, "lower");
  if (!ok) {
    cov += 1e-8 * arma::eye(cov.n_rows, cov.n_cols);
    chol(L, cov, "lower");
  }

  for (int j = 0; j < p; ++j) chi(0, j) = init[j];
  double lp0 = target(init, x, c_r, t_knots, cov_r, n_r, delta_m);
  logpost[0] = lp0;
  int accept_count = 0;

  for (int i = 1; i < N; ++i) {
    if (i % (N/10) == 0) Rcout << "iter: " << i << '\n';
    arma::vec eps = arma::randn<arma::vec>(p);
    arma::vec curr(p);
    for (int j = 0; j < p; ++j) curr[j] = chi(i-1, j);
    arma::vec prop = curr + std::sqrt(eta) * (L * eps);
    NumericVector y_r(p);
    for (int j = 0; j < p; ++j) y_r[j] = prop[j];


    while(true) {
      arma::vec eps = arma::randn<arma::vec>(p);
      arma::vec prop = curr + std::sqrt(eta) * (L * eps);
      NumericVector y_r(p);
      for (int j = 0; j < p; ++j) y_r[j] = prop[j];
      if (is_positive(y_r)) break;
    }
    
    

    double num = target(y_r, x, c_r, t_knots, cov_r, n_r, delta_m);
    double denom = logpost[i-1];
    double log_alpha_k = num - denom;

    double beta_curr = approx_trunc_prob(NumericVector(curr.begin(), curr.end()), cov_r, B);
    double beta_prop = approx_trunc_prob(y_r, cov_r, B);
    // double beta_diff = beta_curr - beta_prop;  // raw difference, no logs

    beta_curr = std::max(beta_curr, 1e-300);
    beta_prop = std::max(beta_prop, 1e-300);
    double log_beta_ratio = std::log(beta_curr) - std::log(beta_prop);
    double log_ratio = log_alpha_k + log_beta_ratio;

    if (R::runif(0.0, 1.0) < std::min(1.0, std::exp(log_ratio))) {
      for (int j = 0; j < p; ++j) chi(i, j) = y_r[j];
      logpost[i] = num;
      ++accept_count;
    } else {
      for (int j = 0; j < p; ++j) chi(i, j) = chi(i-1, j);
      logpost[i] = denom;
    }
  }
  double accept_rate = double(accept_count) / double(N);
  return List::create(Named("chi") = chi,
                      Named("accept_rate") = accept_rate,
                      Named("logpost") = logpost);
}

// cox_bf (Bernoulli factory)
// [[Rcpp::export]]
List cox_bf(int N, const NumericVector &init, const IntegerVector &n_r,
            const List &x, const NumericVector &c_r, const NumericVector &t_knots,
            const NumericMatrix &cov_r, double eta, double delta_m) {
  int p = init.size();
  NumericMatrix chi(N, p);
  NumericVector logpost(N);
  NumericVector bernoulli_loops(N);

  arma::mat cov(const_cast<double*>(cov_r.begin()), cov_r.nrow(), cov_r.ncol(), true);
  arma::mat L;
  bool ok = chol(L, cov, "lower");
  if (!ok) {
    cov += 1e-8 * arma::eye(cov.n_rows, cov.n_cols);
    chol(L, cov, "lower");
  }

  for (int j = 0; j < p; ++j) chi(0, j) = init[j];
  double lp0 = target(init, x, c_r, t_knots, cov_r, n_r, delta_m);
  logpost[0] = lp0;
  int accept_count = 0;

  for (int i = 1; i < N; ++i) {
    if (i % (N/10) == 0) Rcout << "iter: " << i << '\n';
    // propose y
    NumericVector y_r(p);
    while (true) {
      arma::vec eps = arma::randn<arma::vec>(p);
      arma::vec prop = arma::vec(p);
      for (int j = 0; j < p; ++j) prop[j] = chi(i-1, j);
      prop += std::sqrt(eta) * (L * eps);
      for (int j = 0; j < p; ++j) y_r[j] = prop[j];
      if (is_positive(y_r)) break;
    }
    double c1 = target(y_r, x, c_r, t_knots, cov_r, n_r, delta_m);
    double c2 = logpost[i-1];
    double C = std::exp(c1 - c2) / (1.0 + std::exp(c1 - c2));

    bool accepted = false;
    int loops = 0;
    while (!accepted) {
      loops++;
      int C1 = R::rbinom(1, C);
      if (C1 == 1) {
        arma::vec eps = arma::randn<arma::vec>(p);
        arma::vec m1(p);
        for (int j = 0; j < p; ++j) m1[j] = chi(i-1, j);
        m1 += std::sqrt(eta) * (L * eps);
        NumericVector m1_r(p);
        for (int j = 0; j < p; ++j) m1_r[j] = m1[j];
        int p_x = is_positive(m1_r);
        int C2 = R::rbinom(1, p_x);
        if (C2 == 1) {
          for (int j = 0; j < p; ++j) chi(i, j) = y_r[j];
          logpost[i] = c1;
          accepted = true;
        }
      } else {
        arma::vec eps = arma::randn<arma::vec>(p);
        arma::vec m2(p);
        for (int j = 0; j < p; ++j) m2[j] = y_r[j];
        m2 += std::sqrt(eta) * (L * eps);
        NumericVector m2_r(p);
        for (int j = 0; j < p; ++j) m2_r[j] = m2[j];
        int p_y = is_positive(m2_r);
        int C2 = R::rbinom(1, p_y);
        if (C2 == 1) {
          for (int j = 0; j < p; ++j) chi(i, j) = chi(i-1, j);
          logpost[i] = c2;
          accepted = true;
        }
      }
    }
    bernoulli_loops[i] = loops;
    if (chi(i, 0) == y_r[0]) accept_count++;
  }
  double accept_rate = double(accept_count) / double(N);
  return List::create(Named("chi") = chi,
                      Named("bernoulli_loops") = bernoulli_loops,
                      Named("accept_rate") = accept_rate,
                      Named("logpost") = logpost);
}


