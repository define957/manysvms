#include <RcppEigen.h>
#include <Rcpp.h>

using namespace std;
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

//[[Rcpp::export]]
SEXP cpp_spd_solve (Eigen::MatrixXd A, Eigen::MatrixXd b) {
  Eigen::MatrixXd res = A.llt().solve(b);
  return Rcpp::wrap(res);
}
