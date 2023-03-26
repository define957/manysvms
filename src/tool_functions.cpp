// #include <RcppEigen.h>
// #include <Rcpp.h>
//
// using namespace std;
// using namespace Eigen;
// using namespace Rcpp;
// // [[Rcpp::depends(RcppEigen)]]
//
// //[[Rcpp::export]]
// SEXP cpp_chol_solve (Eigen::MatrixXd A, Eigen::MatrixXd b) {
//   Eigen::MatrixXd res = A.llt().solve(b);
//   return Rcpp::wrap(res);
// }
// //[[Rcpp::export]]
// Rcpp::List Eigen_cpp_clip_dcd_optimizer(Eigen::MatrixXd H, Eigen::MatrixXd q,
//                                         Eigen::MatrixXd lb, Eigen::MatrixXd ub,
//                                         double eps, unsigned int max_steps,
//                                         Eigen::MatrixXd u) {
//   unsigned int n = H.rows();
//   unsigned int i = 0;
//   unsigned int j = 0;
//   Eigen::MatrixXd L_idx_val(n, 1);
//   Eigen::MatrixXd L_val(n, 1);
//   Eigen::MatrixXd uH(1, n);
//   Eigen::MatrixXd numerator(n, 1);
//   unsigned int temp_idx;
//
//   for(i = 0; i < max_steps; i++){
//     uH << u.transpose() * H;
//     numerator << q - uH.transpose();
//     for(j = 0; j < n; j++){
//       L_idx_val(j, 0) = numerator(j, 0) / H(j, j);
//       L_val(j, 0) = numerator(j, 0)*numerator(j, 0) / H(j, j);
//     }
//     if(L_val.maxCoeff() < eps){
//       break;
//     }
//     Eigen::MatrixXi unique_idx = ((u.array() > lb.array() && \
//                                    u.array() < L_idx_val.array()) \
//                                   || \
//                                   (u.array() < lb.array() && \
//                                    u.array() > L_idx_val.array())).cast<int>();
//     if(unique_idx.rows() == 0){
//       break;
//     }
//
//     unsigned int max_idx = 0;
//     double lambda_max  = - Eigen::Infinity;
//     double lambda_opt;
//     double max_L_val = - Eigen::Infinity;
//
//     for(j = 0; j < unique_idx.rows(); j++){
//       temp_idx = int(unique_idx(j));
//       if(double(L_val(temp_idx, 0)) > max_L_val){
//         max_L_val = L_val(temp_idx, 0);
//         max_idx = temp_idx;
//       }
//     }
//
//     lambda_max = L_idx_val(max_idx, 0);
//     lambda_opt = std::max(double(lb(max_idx, 0) - u(max_idx, 0)),
//                           std::min(lambda_max,
//                                    double(ub(max_idx, 0) - u(max_idx, 0))));
//     u(max_idx, 0) = u(max_idx, 0) + lambda_opt;
//   }
//   // double obj_val = double(0.5 * u.transpose() * H * u - q.transpose() * u);
//
//   Rcpp::List res = Rcpp::List::create(Rcpp::Named("x") = u,
//                                       Rcpp::Named("iterations") = i);
//                                       //Rcpp::Named("objectiv.value") = obj_val);
//   return res;
// }
