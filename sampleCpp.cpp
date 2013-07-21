#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// [[Rcpp::export]]

NumericVector sampleC(NumericVector num, int size, bool replace)
{
  NumericVector out(size);
  NumericVector out2(1);
  NumericVector num2 = NumericVector::create(0.0,1.0);
  /*out = RcppArmadillo::sample(num2, size, FALSE);
  
  out = rnorm(size,0,1);*/
  
  out = RcppArmadillo::sample(num, size ,TRUE);
  NumericVector out3 = NumericVector::create(1,2,3,4);
  NumericVector xx = NumericVector::create(1,2,3,4);
  
  
  
  return(out);
}

/**
Z <- rnorm(10000,0,1)
Mean34 <- mean((1-pnorm((rho^(-1)-1)^(-0.5)*Z,0,1))^2) 
*/