// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h" 

// [[Rcpp::export(name = ".conditionalMVNRcppCpp")]]
Rcpp::List conditionalMVNRcppCpp(const arma::uvec cdex, 
                                    const arma::uvec gdex, 
                                    const arma::mat xx, arma::mat mu, 
                                    const arma::mat sigma) {
  
  arma::mat sinv = inv(sigma.submat(gdex,gdex));
  arma::mat p1 = sigma.submat(cdex, gdex) * sinv;
  arma::mat mu1 = mu.cols(cdex) + trans(p1 * trans(xx.cols(gdex) - mu.cols(gdex)));
  arma::mat vr1 = sigma.submat(cdex, cdex) - p1 * sigma.submat(gdex,cdex);
  
  return Rcpp::List::create(Rcpp::Named("mu")=mu1,
                            Rcpp::Named("vr")=vr1);
}
