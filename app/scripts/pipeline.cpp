#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List scaleDataset(const arma::mat& V, int iterations) {
  arma::mat V_row = V;
  arma::mat V_column = V;
  for (int i=0; i<=iterations; i++){
    V_row = diagmat(1/arma::sum(V_column,1)) * V_column;
    V_column = V_row * diagmat(1/arma::sum(V_row,0));
  }
  return List::create(Named("V_row") = V_row,
                      Named("V_column") = V_column);
}
// [[Rcpp::export]]
arma::mat getDoubleProjection(const arma::mat& V, const arma::mat& R, const arma::mat& S) {
  arma::mat T = S.t() * S * V * R.t() * R;
  return T;
}

uword getNegative(arma::mat X) {
  uvec q1 = find(X < 0);
  vec B = conv_to<vec>::from(q1);
  return B.n_elem;
}
// [[Rcpp::export]]
double getSum(arma::mat X, arma::mat M) {
  return accu(X) / M.n_rows;
}

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

// [[Rcpp::export]]
arma::mat correctByNorm(arma::mat& X) {
  arma::vec norm_(X.n_rows);
  for (int k = 0; k < X.n_rows; k++) {
    norm_[k] = norm(X.row(k),2);
  }
  mat B = diagmat(1/norm_) * X;
  B.elem( find_nonfinite(B) ).zeros();
  return B;
}

// [[Rcpp::export]]
// added cosine similarity for matrix
arma::rowvec find_cosine(const arma::mat& X) {
  arma::mat Y = arma::trans(X) * X;
  arma::mat res = Y / (arma::sqrt(arma::diagvec(Y)) * arma::trans(arma::sqrt(arma::diagvec(Y))));
  res.diag(0).fill(0);
  arma::rowvec X_dist = max(res,0);
  return X_dist;
}

// [[Rcpp::export]]
arma::uvec update_idx(const arma::mat& prev_X,const arma::mat& new_X, const double thresh=0.8) {
  arma::rowvec prev_values = find_cosine(prev_X);
  arma::rowvec new_values = find_cosine(new_X);
  arma::uvec idx2 = find(new_values >= thresh);
  arma::uvec new_idx = {};
  for (int i=0;i<idx2.n_elem;i++){
    if (new_values.at(idx2[i])>=prev_values.at(idx2[i])) {
      int sz = new_idx.size();
      new_idx.resize(sz+1);
      new_idx(sz) = idx2[i];
    }
  }
  return new_idx;
}

arma::vec nnls_col(const mat &A, const subview_col<double> &b, int max_iter = 500, double tol = 1e-6, bool verbose = false)
{
  
  vec mu = -A.t() * b;
  mat H = A.t() * A;
  vec x(A.n_cols), x0(A.n_cols);
  x.fill(0);
  x0.fill(-9999);
  
  int i = 0;
  double tmp;
  while(i < max_iter && max(abs(x - x0)) > tol) 
  {
    x0 = x;
    for (int k = 0; k < A.n_cols; k++) 
    {
      tmp = x[k] - mu[k] / H.at(k,k);
      if (tmp < 0) tmp = 0;
      if (tmp != x[k]) mu += (tmp - x[k]) * H.col(k);
      x[k] = tmp;
    }
    ++i;
  }
  
  return x;
}

// [[Rcpp::export]]
arma::mat nnls_C__(mat A, mat b, int max_iter = 500, double tol = 1e-6)
{
  // solving Ax = b, where x and b are both matrices
  if(A.n_rows != b.n_rows)
    throw std::invalid_argument("A and b must have the same number of rows.");
  mat x(A.n_cols, b.n_cols);
  for (int i = 0; i < b.n_cols; i++)
    x.col(i) = nnls_col(A, b.col(i), max_iter, tol);
  
  return x;
}


// [[Rcpp::export]]
arma::mat hinge_der_proportions_C__(const arma::mat& H,
                                    const arma::mat& R, double precision_ = 1e-10) {
  int m = H.n_rows;
  int n = H.n_cols;
  
  arma::mat TMP(n, m*m, fill::zeros);
  
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      if (H(i,j) < 0) {
        for (int k = 1; k < m; k++) {
          TMP(j,k+i*m) = -R(k,j);
        }
      }
    }
  }
  
  return reshape(arma::sum(TMP,0),m,m).t();
  
}
// [[Rcpp::export]]
arma::mat hinge_der_basis_C__(const arma::mat& W,
                              const arma::mat& S, 
                              double precision_ = 1e-10) {
  int n = W.n_cols;
  
  arma::mat res(n, n, fill::zeros);
  
  for (int j = 0; j < n; j++) {
    
    arma::vec t = W.col(j); 
    res.col(j) =  arma::sum(-S.cols(find(t < -precision_)),1);
  }
  res.row(0).zeros();
  
  return res;
  
}
// [[Rcpp::export]]
double hinge_C__(const arma::mat& X) {
  arma::mat X_(X.n_rows,X.n_cols,fill::zeros);
  double elem_ = 0;
  
  for (int i=0; i<X.n_rows; i++) {
    for (int j=0; j<X.n_cols; j++) {
      elem_ = X(i,j);
      if (elem_<0){
        X_(i,j) = -elem_;
      }
    }
  }
  return accu(X_);
}
// [[Rcpp::export]]
List calcErrors(const arma::mat& X,
                const arma::mat& Omega,
                const arma::mat& D_w,
                const arma::mat& D_h,
                const arma::mat& V_row,
                const arma::mat& R,
                const arma::mat& S,
                const double coef_,
                const double coef_der_X,
                const double coef_der_Omega,
                const double coef_hinge_H,
                const double coef_hinge_W,
                const double coef_pos_D_h,
                const double coef_pos_D_w) {
  arma::mat V__ = S * V_row * R.t();
  arma::mat D_w_diag = diagmat(D_w);
  
  double error_ = pow(norm(V__ - Omega * D_w_diag * X,"fro"),2.0);
  double orig_deconv_error = pow(norm(V_row - S.t() * Omega * D_w_diag * X * R,"fro"),2);
  double lambda_error = coef_ * coef_hinge_H * hinge_C__(X * R);
  double beta_error = coef_ * coef_hinge_W * hinge_C__(S.t() * Omega);
  arma::mat A = arma::sum(R,1);
  arma::mat B = arma::sum(S,1);
  double D_h_error = coef_pos_D_h * pow(norm(X.t() * D_h - A,"fro"),2);
  double D_w_error = coef_pos_D_w * pow(norm(Omega * D_w - B,"fro"),2);
  double new_error = error_ + lambda_error + beta_error + D_h_error + D_w_error;
  
  return List::create(Named("error_") = error_,
                      Named("lambda_error") = lambda_error,
                      Named("beta_error") = beta_error,
                      Named("D_h_error") = D_h_error,
                      Named("D_w_error") = D_w_error,
                      Named("new_error") = new_error,
                      Named("orig_deconv_error") = orig_deconv_error);
}


// [[Rcpp::export]]
field<mat> derivative_stage2(const arma::mat& X,
                             const arma::mat& Omega,
                             const arma::mat& D_w,
                             const arma::mat& V_row,
                             const arma::mat& R,
                             const arma::mat& S,
                             const double coef_der_X,
                             const double coef_der_Omega,
                             const double coef_hinge_H,
                             const double coef_hinge_W,
                             const double coef_pos_D_h,
                             const double coef_pos_D_w,
                             const int cell_types,
                             const double N,
                             const double M,
                             const int global_iterations,
                             arma::mat& errors_statistics,
                             const int start_idx,
                             arma::mat& points_statistics_X,
                             arma::mat& points_statistics_Omega,
                             const double mean_radius_X,
                             const double mean_radius_Omega,
                             const double r_const_X=0,
                             const double r_const_Omega=0,
                             const double thresh=0.8) {
  arma::mat new_X = X;
  arma::mat new_Omega = Omega;
  arma::mat new_D_w = D_w;
  arma::mat new_D_h = new_D_w * (N/M);
  arma::mat jump_X, jump_Omega;
  
  arma::mat V__ = S * V_row * R.t();
  mat B = join_cols(arma::vectorise(V__),coef_pos_D_w * arma::sum(S,1));
  mat C = join_cols(arma::vectorise(V__),coef_pos_D_h * arma::sum(R,1));
  arma::mat der_X,der_Omega;
  
  for (int itr_= start_idx; itr_ < global_iterations+start_idx; itr_++) {
    
    bool has_jump_X = false;
    bool has_jump_Omega = false;
    // derivative X
    der_X = -2 * (diagmat(new_D_w) * new_Omega.t() * (V__ - new_Omega * diagmat(new_D_w) * new_X));
    der_X += coef_hinge_H * hinge_der_proportions_C__(new_X * R, R);
    der_X += coef_pos_D_h * 2 * new_D_h * (new_X.t() * new_D_h - arma::sum(R,1)).t();
    der_X.col(0).zeros();
    der_X = correctByNorm(der_X) * mean_radius_X;
    
    if (thresh > 0) {
      arma::mat tmp_X = (new_X - coef_der_X * der_X).t();
      arma::mat tmp_X_2 = (new_X).t();
      arma::uvec idx = update_idx(tmp_X,tmp_X_2,thresh);
    
      if (idx.n_elem > 0) {
        der_X.rows(idx).zeros();
      }
    }
    
    
    new_X = new_X - coef_der_X * der_X;
    
    if (r_const_X > 0) {
      jump_X = jump_norm(new_X,r_const_X);
      new_X = new_X % jump_X;
    }
    
    arma::mat vec_mtx(cell_types*cell_types,cell_types,fill::zeros);
    for (int c=0; c<cell_types; c++) {
      vec_mtx.col(c) = arma::vectorise(new_Omega.col(c) * new_X.row(c));
    }
    arma::mat A = join_cols((M/N) * vec_mtx, coef_pos_D_h * new_X.t());
    
    
    new_D_h = nnls_C__(A, C);
    new_D_w = new_D_h * (M/N);
    
    // derivative Omega
    der_Omega = -2 * (V__ - new_Omega * diagmat(new_D_w) * new_X) * new_X.t() * diagmat(new_D_w);
    der_Omega += coef_hinge_W * hinge_der_basis_C__(S.t() * new_Omega, S);
    der_Omega += coef_pos_D_w * 2 * (new_Omega*new_D_w-arma::sum(S,1)) * new_D_w.t();
    der_Omega.row(0).zeros();
    der_Omega = correctByNorm(der_Omega) * mean_radius_Omega;
    

    if (thresh > 0) {
      arma::mat tmp_Omega = new_Omega - coef_der_Omega * der_Omega;
      arma::uvec idx2 = update_idx(tmp_Omega,new_Omega,thresh);
      
      if (idx2.n_elem > 0) {
        der_Omega.cols(idx2).zeros();
      }
    }
    
    new_Omega = new_Omega - coef_der_Omega * der_Omega;

    if (r_const_Omega > 0) {
      arma::mat t_Omega = new_Omega.t();
      jump_Omega = jump_norm(t_Omega,r_const_Omega);
      jump_Omega = jump_Omega.t();
      //has_jump_Omega = any(jump_Omega != 1);
      new_Omega = new_Omega % jump_Omega;
    }
    
    
    vec_mtx.fill(fill::zeros);
    A.fill(fill::zeros);
    
    for (int c=0; c<cell_types; c++) {
      vec_mtx.col(c) = arma::vectorise(new_Omega.col(c) * new_X.row(c));
    }
    A = join_cols(vec_mtx, coef_pos_D_w * new_Omega);
    
    
    new_D_w = nnls_C__(A, B);
    
    new_D_h = new_D_w * (N/M);
    
    uword neg_props = getNegative(new_X * R);
    uword neg_basis = getNegative(S.t() * new_Omega);
    double sum_ = accu(new_D_w) / V_row.n_rows;
    
    List err_ = calcErrors(new_X,new_Omega,new_D_w, new_D_h,
                           V_row, R, S, 1, coef_der_X,
                           coef_der_Omega, coef_hinge_H,
                           coef_hinge_W, coef_pos_D_h,
                           coef_pos_D_w);
    errors_statistics.row(itr_) = { err_["error_"],
                                    err_["lambda_error"],
                                        err_["beta_error"],
                                            err_["D_h_error"],
                                                err_["D_w_error"],
                                                    err_["new_error"],
                                                        err_["orig_deconv_error"],
                                                            neg_props,
                                                            neg_basis,
                                                            sum_};
    points_statistics_X.row(itr_) = new_X.as_row();
    points_statistics_Omega.row(itr_) = new_Omega.as_row();
    
    
  }
  
  field<mat> ret_(7,1);
  ret_(0,0) = new_X;
  ret_(1,0) = new_Omega;
  ret_(2,0) = new_D_w;
  ret_(3,0) = new_D_h;
  ret_(4,0) = errors_statistics;
  ret_(5,0) = points_statistics_X;
  ret_(6,0) = points_statistics_Omega;
  
  return ret_;
  
}




