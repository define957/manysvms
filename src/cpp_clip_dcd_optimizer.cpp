#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
SEXP cpp_clip_dcd_optimizer(arma::mat H, arma::mat q,
                            arma::mat lb, arma::mat ub,
                            double eps, unsigned int max_steps){
  arma::mat u = (lb + ub) / 2;
  unsigned int n = H.n_rows;

  unsigned int i = 0;
  unsigned int j = 0;

  for(i = 0; i < max_steps; i++){

    arma::mat numerator = q - arma::trans(arma::trans(u) * H);
    arma::vec L_idx_val(numerator.n_rows);
    arma::vec L_val(numerator.n_rows);

    for(j = 0; j < n; j++){
      L_idx_val(j)= numerator(j) / H(j,j);
    }
    for(j = 0; j < n; j++){
      L_val(j)= std::pow(numerator(j), 2) / H(j,j);
    }

    if(as_scalar(max(L_val)) < eps){
      break;
    }
    arma::uvec idx1 = arma::find(u > lb && L_idx_val < 0);
    arma::uvec idx2 = arma::find(u < ub && L_idx_val > 0);
    arma::uvec unique_idx = arma::unique(arma::join_cols(idx1, idx2));

    if(unique_idx.n_elem == 0){
      break;
    }

    unsigned int max_idx = 0;
    double lambda_max = - datum::inf;
    double lambda_opt;
    double max_L_val = - datum::inf;

    for(j = 0; j < unique_idx.n_elem; j++){
      unsigned int temp_idx = as_scalar(unique_idx(j));
      if(L_val(temp_idx) > max_L_val){
        max_L_val = L_val(temp_idx);
        max_idx = temp_idx;
      }
    }

    lambda_max = L_idx_val(max_idx);
    lambda_opt = std::max(as_scalar(lb(max_idx) - u(max_idx)),
                          std::min(lambda_max, as_scalar(ub(max_idx) - u(max_idx))));
    u(max_idx) = u(max_idx) + lambda_opt;
  }
  return Rcpp::wrap(u);
}
