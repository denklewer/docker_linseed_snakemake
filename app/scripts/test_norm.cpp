#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat jump_norm(arma::mat& X, const double r_const_X = 0) {
  arma::mat norm_(X.n_rows,X.n_cols);
  norm_.fill(1.0);
  arma::mat X_trunc(X.n_rows,X.n_cols-1);
  arma::uvec ids = arma::regspace<arma::uvec>(1, X.n_cols-1);
  X_trunc = X.cols(ids);
  for (int k = 0; k < X.n_rows; k++) {
    double row_norm = norm(X_trunc.row(k),2);
    for (int j = 1; j < X.n_cols; j++) {
      if (r_const_X>row_norm) {
        norm_.at(k,j) = r_const_X/row_norm;
      } else {
        norm_.at(k,j) = 1;
      }  
    }
  }
  return norm_;
}


/*** R
# toSave <- readRDS("/Users/aladyeva.e/Dropbox (ArtyomovLab)/Deconvolution/norm_test_params.rds")
# initX <- toSave$initX
# R_0_75_new <- toSave$R_0_75_new
jump_norm(initX,R_0_75_new)
*/
