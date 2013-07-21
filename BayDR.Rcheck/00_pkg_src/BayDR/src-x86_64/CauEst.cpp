
/**
  Causal <- apply(CSample,1,function(x)
  {
    tem <- gamma0+gamma2/sqrt(K)*(sum(x[seq(1,((r-2)*K),by=r)])-0.5*K)+
      lambda2/sqrt(K)*(sum(x[seq(3,(r*K),by=r)]*x[seq(4,(r*K),by=r)])-Mean34*K)
    (inv.logit(tem+gamma1+lambda1/sqrt(K)*(sum(x[seq(2,(r*K),by=r)])-0.5*K)) - inv.logit(tem))
    
  }
  )
*/

/* #include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
*/

#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]
double CauEstC(double mcerror, double gamma0, double gamma1,double gamma2,
                      double lambda1,double lambda2, int r,int K,double mean34) {
                        
  
  double out=0;
  double N = pow(0.5/mcerror,2);
  NumericVector num = NumericVector::create(0.0,1.0);
  for (int i=0;i<N;++i)
  {
    double s1=0,s2=0,s34=0;
    NumericVector csample = RcppArmadillo::sample(num, r*K ,TRUE);
    for (int j=0;j<K;++j)
    {
      s1+=csample(0+j*r);
      s2+=csample(1+j*r);
      s34+=csample(2+j*r)*csample(3+j*r);
    }
    double tp,tp2;
    tp = gamma0+gamma2/sqrt(K)*(s1-0.5*K)+lambda2/sqrt(K)*(s34-mean34*K);
    tp2 = tp+gamma1+lambda1/sqrt(K)*(s2-0.5*K);
    out+=exp(tp2)/(exp(tp2)+1) - exp(tp)/(exp(tp)+1) ;
  }
  out = out/N;
  return out;
}

/**
  CausalEst <- function(K=2,r=6,rho=0.3,gamma0=0.2,gamma1=1,gamma2=1,lambda1=1,lambda2=0,CSample) {
    
    Z <- rnorm(10000,0,1)
    Mean34 <- mean((1-pnorm((rho^(-1)-1)^(-0.5)*Z,0,1))^2)
    K=2;r=6;rho=0.3;gamma0=0.2;gamma1=1;gamma2=1;lambda1=1;lambda2=0;mc.error=0.001
    causal2 <- CauEstC(CSample,gamma0,gamma1,gamma2,lambda1,lambda2,r,K,Mean34)
    
    mean(Causal)
    
  }

*/
  
  
  
  
