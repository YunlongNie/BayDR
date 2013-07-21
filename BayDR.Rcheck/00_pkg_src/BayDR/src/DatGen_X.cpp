
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
double VSumC(NumericVector x, NumericVector y){
  double out=0; 
  int length=y.size();
  for (int i=0;i<length;++i)
  {
    out+=x[y[i]-1];
  }
  return out;
}

double logitC(double x){
  double out=log(x/(1-x)); 
  return out;
}


// [[Rcpp::export]]
List DatDen_X_C(NumericMatrix C, double phi0, double phi1,double phi2,
  double lambda1,double lambda2,int K,
  NumericVector No123,NumericVector No1,NumericVector No2,NumericVector No34, double Mean34) 
  {
  int nrow=C.nrow();
  List out(2);
  NumericVector X(nrow);
  NumericVector Y(nrow);
  double alpha0 = (logitC(0.5) + logitC(0.05))/2;
  double alpha1 = (logitC(0.5) - logitC(0.05))/4;
  
  for (int i=0;i<nrow;++i)
  {
   double s123=VSumC(C(i,_),No123);
   double temp=alpha0+alpha1/sqrt(K)*(s123-K*1.5);
   X[i] = rbinom(1,1,exp(temp)/(1+exp(temp)))[0];
   double tp =phi0+phi1*X[i]+phi2/sqrt(K)*(VSumC(C(i,_),No1)-0.5*K) + 
   lambda1*X[i]/sqrt(K)*(VSumC(C(i,_),No2)-0.5*K)+
   lambda2/sqrt(K)*(VSumC(C(i,_),No34)-Mean34*K);
   Y[i]= rbinom(1,1,exp(tp)/(exp(tp)+1))[0];
  }
  out[0]=X;
  out[1]=Y;
  return(out);
  }
   /**
   rbinom(1,size=1,prob=inv.logit(alpha0+alpha1/sqrt(K)*(
      sum(x[No123])-K*1.5
    )))
    */
    
    
  
  /**
  Dat$Y <- apply(Dat,1,function(x)
  {
    rbinom(1,size=1,
    prob=inv.logit(phi0+ phi1*x[r*K+1]+phi2/sqrt(K)*(sum(x[seq(1,(r*K),by=r)])-0.5*K)+
    lambda1*x[r*K+1]/sqrt(K)*(sum(x[seq(2,(r*K),by=r)])-0.5*K)+
      lambda2/sqrt(K)*(sum(x[c(seq(3,(r*K),by=r),seq(4,(r*K),by=r))])-Mean34*K)
    ))
  }
}

*/
    
    
  
  
