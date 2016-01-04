// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h" 

// [[Rcpp::export(name = ".byRccpCpp")]]
Rcpp::List byRccpCpp(const int nr, 
                      const arma::mat frommat,
                      arma::mat totmat, 
                      arma::mat summat, 
                      arma::mat minmat, 
                      arma::mat maxmat){
  int k, i, j;
  double s;
  
  for(k = 0; k < nr; k++){
    
    i = frommat(k,0) - 1;
    j = frommat(k,1) - 1;
    s = frommat(k,2);
    totmat(i,j) = totmat(i,j) + 1;
    summat(i,j) = summat(i,j) + s;
    
    if(s > maxmat(i,j))
      maxmat(i,j) = s;   
    
    if(s < minmat(i,j))
      minmat(i,j) = s;
  }
  
  return Rcpp::List::create(Rcpp::Named("total")=totmat,
                            Rcpp::Named("sum")=summat,
                            Rcpp::Named("min")=minmat,
                            Rcpp::Named("max")=maxmat);
}
