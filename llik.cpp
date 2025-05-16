#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List llik_cpp(
  const mat& Wyt, const vec& vOmega, const vec& vA, const vec& vB, 
  const vec& f1) {

  // get the formats
  int dim_N = Wyt.n_cols;
  int dim_T = Wyt.n_rows;

  // initialize storage objects
  mat ft(dim_T, dim_N);
  vec likt(dim_T);
  vec f_now = f1;
  
  // loop filter over observations  
  for (int i = 0; i < dim_T; i++) {
    ft.row(i) = f_now.t();
    vec ee2 = Wyt.row(i).t(); ee2 = ee2 % ee2;
    likt(i) -= 0.5 * sum(log(f_now)) + 0.5 * sum( ee2 / f_now);
      f_now = (1 - vB) % vOmega + vB % f_now + vA % (ee2 - f_now);
  }
  likt -= 0.5 * dim_N * log(2 * M_PI);

  double llik = mean(likt) / dim_N;
  if (llik != llik) llik = -1e10;
  
  return List::create(Named("llik") = llik,
                      Named("lik_t") = likt,
                      Named("ft") = ft);
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List llik1_cpp(
    const vec& Wyt, const double& dOmega, const double& dlogisticA, 
    const double& dlogisticB, const double& f1) {
  
  // get the formats
  int dim_T = Wyt.n_elem;
  
  // initialize storage objects
  vec ft(dim_T);
  vec likt(dim_T);
  double f_now = f1;
  if ((fabs(dlogisticA) > 10) || (fabs(dlogisticB) > 10)) 
    return List::create(Named("llik") = -1e10);
  double dA = 1 / (1 + exp(-dlogisticA));
  double dB = 1 / (1 + exp(-dlogisticB));
  
  // loop filter over observations  
  for (int i = 0; i < dim_T; i++) {
    ft(i) = f_now;
    double ee2 = Wyt(i); ee2 *= ee2;
    likt(i) = -0.5 * log(f_now) - 0.5 * ee2 / f_now;
    f_now = (1 - dB) * dOmega + dB * f_now + dA * (ee2 - f_now);
  }
  likt -= 0.5 * log(2 * M_PI);
  
  double llik = mean(likt);
  if (llik != llik) llik = -1e10;
  
  return List::create(Named("llik") = llik,
                      Named("lik_t") = likt,
                      Named("ft") = ft);
}



