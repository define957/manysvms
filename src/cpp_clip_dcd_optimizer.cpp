#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
SEXP cpp_clip_dcd_optimizer(arma::mat H, arma::mat q,
                            arma::mat lb, arma::mat ub,
                            double eps, int max_steps){
  arma::mat u = (lb + ub) / 2;
  int n = H.n_rows;

  int i = 0;
  int j = 0;

  for(i = 0; i < max_steps; i++){

    arma::mat numerator = q - arma::trans(arma::trans(u) * H);
    arma::mat L_idx_val(numerator.n_rows, numerator.n_cols);
    arma::mat L_val(numerator.n_rows, numerator.n_cols);

    for(j = 0; j < n; j++){
      L_idx_val(j)= numerator(j) / H(j,j);
    }
    for(j = 0; j < n; j++){
      L_val(j)= pow(numerator(j), 2) / H(j,j);
    }

    if(as_scalar(max(L_val)) < eps){
      break;
    }

    // check index for optimization
    int max_idx = 0;
    double lambda_max =  - datum::inf;
    double lambda_opt;
    int cnt = 0;

    for(j = 0; j < n; j ++){
      if(((u(j) > lb(j)) && (L_idx_val(j) < 0)) ||
         ((u(j) < ub(j)) && (L_idx_val(j) > 0))){
        if(L_idx_val(j) > lambda_max){
          lambda_max = L_idx_val(j);
          max_idx = j;
        }
        cnt ++;
      }
    }
    if(cnt == 0){
      break;
    }


    lambda_opt = std::max(as_scalar(lb(max_idx) - u(max_idx)),
                          std::min(lambda_max, as_scalar(ub(max_idx) - u(max_idx))));
    u(max_idx) = u(max_idx) + lambda_opt;
  }
  return Rcpp::wrap(u);
}
