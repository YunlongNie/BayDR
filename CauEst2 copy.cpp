
/**
Causal <- apply(CSample,1,function(x)
{
  tem <- phi0+phi2/sqrt(K)*(sum(x[seq(1,((r-2)*K),by=r)])-0.5*K)+
    lambda2/sqrt(K)*(sum(x[seq(3,(r*K),by=r)]*x[seq(4,(r*K),by=r)])-Mean34*K)
  (inv.logit(tem+phi1+lambda1/sqrt(K)*(sum(x[seq(2,(r*K),by=r)])-0.5*K)) - inv.logit(tem))
  
}
)
*/

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector CauEstC2(NumericMatrix CSample, double phi0, double phi1,double phi2,
   double lambda1,double lambda2, int r,int K,double mean34) {
  int nrow=CSample.nrow();
  NumericVector out(nrow);
  for (int i=0;i<nrow;++i)
  {
   double s1=0,s2=0,s34=0;
     for (int j=0;j<K;++j)
    {
      s1+=CSample(i,0+j*r);
      s2+=CSample(i,1+j*r);
      s34+=CSample(i,2+j*r)*CSample(i,3+j*r);
    }
    double tp,tp2;
    tp = phi0+phi2/sqrt(K)*(s1-0.5*K)+lambda2/sqrt(K)*(s34-mean34*K);
    tp2 = tp+phi1+lambda1/sqrt(K)*(s2-0.5*K);
    out[i]=exp(tp2)/(exp(tp2)+1) - exp(tp)/(exp(tp)+1) ;
}

  return out;
  }
  
  /**
  CausalEst <- function(K=2,r=6,rho=0.3,phi0=0.2,phi1=1,phi2=1,lambda1=1,lambda2=0,CSample) {

Z <- rnorm(10000,0,1)
Mean34 <- mean((1-pnorm((rho^(-1)-1)^(-0.5)*Z,0,1))^2)

causal2 <- CauEstC(CSample,phi0,phi1,phi2,lambda1,lambda2,r,K,Mean34)

mean(Causal)

}

*/
    
    
  
  