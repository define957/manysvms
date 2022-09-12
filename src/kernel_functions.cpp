#include <RcppArmadillo.h>

using namespace std;
using namespace arma;

//[[Rcpp::export]]
SEXP cpp_rbf_kernel(arma::mat x1, arma::mat x2, double gamma)
{
  int n = x1.n_rows;
  int m = x2.n_rows;
  int i = 0;
  int j = 0;
  arma::mat kernelx(n, m);
  for(i = 0; i < n; i ++){
    for(j = 0; j < n; j ++){
      kernelx(i, j) = exp(as_scalar(- gamma * pow(norm(x1.row(i) - x2.row(j),2), 2)));
    }
  }
  return Rcpp::wrap(kernelx);
}

//[[Rcpp::export]]
SEXP cpp_poly_kernel(arma::mat x1, arma::mat x2, double gamma, int degree, double coef0)
{
  arma::mat kernelx = arma::pow((gamma * x1 * x2.t() + coef0), degree);
  return Rcpp::wrap(kernelx);
}



