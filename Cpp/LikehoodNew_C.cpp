#include <Rcpp.h>
using namespace Rcpp;
 
// [[Rcpp::export]]
double Like_C(NumericMatrix Comb1,NumericMatrix Comb0, NumericMatrix Number) 
{
  int nrow1=Comb1.nrow(),nrow0=Comb0.nrow();
  int nrowC=Number.nrow();
  NumericVector p1(nrow1),p0(nrow0);
  
    for (int i=0;i<nrow1;++i)
    {
      double s1=0;
      for (int j=0;j<nrowC;++j)
      {
        s1+=R::lbeta(Comb1(i,0)+Number(j,0),Comb1(i,1)+Number(j,1)-Number(j,0))-
        R::lbeta(Comb1(i,0),Comb1(i,1)); 
      }
      p1[i] = exp(log(Comb1(i,2)) + s1);
    }
  for (int i=0;i<nrow0;++i)
    {
      double s0=0;
      for (int j=0;j<nrowC;++j)
      {
        s0+=R::lbeta(Comb0(i,0)+Number(j,2),Comb0(i,1)+Number(j,3)-Number(j,2))-
        R::lbeta(Comb0(i,0),Comb0(i,1)); 
      }
      p0[i] = exp(log(Comb0(i,2)) + s0);
    }
   double out=log(sum(p0)*sum(p1));
   return(out);
}


/**
 * 
x=as.matrix(t.comb1[1,])
y=as.matrix(t.N..Number[1,])
out <- Like_C(x,x,y)
out
out <- Like_C(as.matrix(t.comb1),as.matrix(t.comb0),
Number=as.matrix(t.N..Number))


likelihood.new <- log(sum(apply(t.comb1,1,function(x) gg(as.numeric(x),t.N..Number1))))+ # likehihood of estimate.s.new
    log(sum(apply(t.comb0,1,function(x) gg(as.numeric(x),t.N..Number0))))
 * 
 */

