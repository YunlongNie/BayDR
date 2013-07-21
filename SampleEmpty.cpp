#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// [[Rcpp::export]]

double sampleEmpty2(double mcerror, NumericMatrix UqC, NumericVector coef,NumericVector seq,int type)
{
  
  NumericVector ext(UqC.nrow());
  for (int k=0;k<UqC.nrow();k++)
  {
    ext(k) = sum(seq*UqC(k,_)); 
  }
  
  /*return(ext);*/
  
  double out;
  
  /*NumericVector out(seq.size());*/
  int N = pow(0.5/mcerror,2);
  int clength=coef.size();
  NumericVector num = NumericVector::create(0.0,1.0);
  for (int i=0;i<N;++i)
  { 
    
    NumericVector csample(UqC.ncol()+2);
      
   if (type==0)
   {
   int m = 1;
   while (m < ext.size()-1)
    {
      csample=RcppArmadillo::sample(num, UqC.ncol()+2 ,TRUE);
      NumericVector cstemp(seq.size());
      for (int q= 0;q<seq.size();q++)
      {
        cstemp(q) = csample(q+1);
      }
      double nw = sum(seq*cstemp);
      int l=0;
      while (nw!=ext(l))
      {
        l=l+1;
        if(l>=ext.size()-1) break;
      }
      m=l;
    } 
   } else 
   {
     csample=RcppArmadillo::sample(num, UqC.ncol()+2 ,TRUE);
     }
    csample(0)=1;
    csample(clength-1)=1;
    double s1 = exp(sum(csample*coef))/(1+exp(sum(csample*coef)));
    csample(clength-1)=0;
    double s0 = exp(sum(csample*coef))/(1+exp(sum(csample*coef)));
    out = out + s1-s0;
  }
  out = out/N ;
  return out;
}

/**
    CSample = sample.empty(N=(0.5/Mc.error)^2,UqC=as.matrix(UqC))}
    Hatm1C1full <- inv.logit(as.matrix(cbind(1,CSample,1))%*%as.matrix(coef(outModel)))
    Hatm0C0full <- inv.logit(as.matrix(cbind(1,CSample,0))%*%as.matrix(coef(outModel)))
    
    estimate.p <- mean(Hatm1C1full-Hatm0C0full)*k_q*n.mpty/(level.con*k_q+no.ob) + sum(PrC*(Hatm1C1-Hatm0C0))
    NumericVector out(1);
  NumericVector out2(1);
  NumericVector num2 = NumericVector::create(0.0,1.0);
 
  NumericVector out3 = NumericVector::create(1,2,3,4);
  NumericVector xx = NumericVector::create(1,2,3,4);
  out2 = xx*out3;
  return(out2);
*/