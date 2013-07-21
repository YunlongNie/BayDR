#include <Rcpp.h>
using namespace Rcpp;
 
// [[Rcpp::export]]
NumericVector Comb_C(NumericMatrix Comb, NumericMatrix Number) 
  {
  int mfull=Number.nrow(),nrowComb = Comb.nrow();
  NumericVector out(nrowComb);
  for (int i=0;i<nrowComb;++i)
  {
    double stemp=0;
    for (int j=0;j<mfull;++j)
    {
      
      stemp+=R::lbeta(Comb(i,0)+Number(j,0),Comb(i,1)+Number(j,1)-Number(j,0));
      /**sum(apply(N..Number,1,function(x)
        lbeta(xx[1]+x[1],xx[2]+x[2]-x[1])
           )) + log(xx[3]) - mfull*lbeta(xx[1],xx[2])*/
           }
      out[i]=stemp + log(Comb(i,2))-mfull*R::lbeta(Comb(i,0),Comb(i,1));   
      
  }
  double maxout=max(out);
  for (int i=0;i<nrowComb;++i)
  {
    out[i]=exp(out[i]-maxout);
  }
  double sumTemp=sum(out);
  for (int i=0;i<nrowComb;++i)
  {
    out[i]=out[i]/sumTemp;
  }
  return(out);
  }

/**
Comb <- function(comb1,N..Number){   ## r 
mfull <- nrow(N..Number)
ff1 <- function(xx) # xx = c(a0,b0,r(a0,b0)) # 
  {
    
    
      sum(apply(N..Number,1,function(x)
        lbeta(xx[1]+x[1],xx[2]+x[2]-x[1])
           )) + log(xx[3]) - mfull*lbeta(xx[1],xx[2])

  } 
  
  ss1 <- apply(comb1,1,function(x) 
  {
    ff1(as.numeric(x))  
  }
  )
  expss11<- exp(ss1-max(ss1))
  comb1$w <- expss11/sum(expss11)
  return(comb1)
}
*/
