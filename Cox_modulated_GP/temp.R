// [[Rcpp::export]]
List cox_mh(int N, const NumericVector &init, const IntegerVector &n_r,
            const List &x, const NumericVector &c_r, const NumericVector &t_knots,
            const NumericMatrix &cov_r, double eta, double delta_m) {
  
  int p = init.size();
  NumericMatrix chi(N, p);
  NumericVector logpost(N);
  
  // Cholesky decomposition
  arma::mat cov(const_cast<double*>(cov_r.begin()), cov_r.nrow(), cov_r.ncol(), true);
  arma::mat L;
  bool ok = chol(L, cov, "lower");
  if (!ok) {
    cov += 1e-8 * arma::eye(cov.n_rows, cov.n_cols);
    chol(L, cov, "lower");
  }
  
  // Initialize
  for (int j = 0; j < p; ++j) chi(0, j) = init[j];
  logpost[0] = target(init, x, c_r, t_knots, cov_r, n_r, delta_m);
  int accept_count = 0;
  
  // R function pmvnorm
  Rcpp::Function pmvnorm("pmvnorm");
  NumericVector lower(p, 0.0);
  NumericVector upper(p, R_PosInf);
  
  for (int i = 1; i < N; ++i) {
    if (i % (N/10) == 0) Rcout << "iter: " << i << "\n";
    
    // Current vector
    arma::vec curr(p);
    for (int j = 0; j < p; ++j) curr[j] = chi(i-1, j);
    
    // Proposal in positive orthant
    NumericVector y_r(p);
    while (true) {
      arma::vec eps = arma::randn<arma::vec>(p);
      arma::vec prop = curr + std::sqrt(eta) * (L * eps);
      for (int j = 0; j < p; ++j) y_r[j] = prop[j];
      if (is_positive(y_r)) break;
    }
    
    // Compute log-posterior
    double num = target(y_r, x, c_r, t_knots, cov_r, n_r, delta_m);
    double denom = logpost[i-1];
    double log_alpha_k = num - denom;
    
    // Orthant probabilities using pmvnorm
    NumericVector curr_r(p);
    for (int j = 0; j < p; ++j) curr_r[j] = curr[j];
    
    double beta_curr = as<double>(pmvnorm(Named("mean") = curr_r,
                                          Named("sigma") = cov,
                                          Named("lower") = lower,
                                          Named("upper") = upper));
    double beta_prop = as<double>(pmvnorm(Named("mean") = y_r,
                                          Named("sigma") = cov,
                                          Named("lower") = lower,
                                          Named("upper") = upper));
    
    beta_curr = std::max(beta_curr, 1e-300);
    beta_prop = std::max(beta_prop, 1e-300);
    
    double log_beta_ratio = std::log(beta_curr) - std::log(beta_prop);
    double log_ratio = log_alpha_k + log_beta_ratio;
    
    // Accept/reject
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


