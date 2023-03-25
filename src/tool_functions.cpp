#include <RcppEigen.h>
#include <Rcpp.h>

using namespace std;
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

//[[Rcpp::export]]
SEXP cpp_chol_solve (Eigen::MatrixXd A, Eigen::MatrixXd b) {
  Eigen::MatrixXd res = A.llt().solve(b);
  return Rcpp::wrap(res);
}

// Rcpp::List Eigen_cpp_clip_dcd_optimizer(Eigen::MatrixXd H, Eigen::MatrixXd q,
//                                   Eigen::MatrixXd lb, Eigen::MatrixXd ub,
//                                   double eps, unsigned int max_steps,
//                                   Eigen::MatrixXd u){
//   unsigned int n = H.rows();
//
//   unsigned int i = 0;
//   unsigned int j = 0;
//
//   for(i = 0; i < max_steps; i++){
//     Eigen::MatrixXd numerator = q - (u.transpose() * H).transpose;
//   }
// }
