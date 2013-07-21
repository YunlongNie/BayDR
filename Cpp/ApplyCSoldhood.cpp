/**
#   likelihood.old <- sum(apply(N..Number,1,function(x) {  # likehihood of estimate.s.old 
#     (
#       lgamma(x[1]+a1)+lgamma(x[2]-x[1]+b1)+lgamma(a1+b1)-
#         lgamma(x[2]+a1+b1)-lgamma(a1)-lgamma(b1)+
#         lgamma(x[3]+a0)+lgamma(x[4]-x[3]+b0)+lgamma(a0+b0)-
#         lgamma(x[4]+a0+b0)-lgamma(a0)-lgamma(b0)
#     )
#   }
#   )
#   )

*/


#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector ApplyCSoldhood(NumericMatrix Number, double a1, double b1, double a0, double b0) 
  {
  int nRow=Number.nrow();
  NumericVector out(nRow);
  for (int i=0;i<nRow;++i)
  {
    
//      out[i]= (kq + Number(i,4))/(levelCon*kq+noObs)*
//      ((a1+Number(i,0))/(a1+b1+Number(i,1))-(a0+Number(i,2))/(a0+b0+Number(i,3)));
       
    out[i]  = R::lbeta(Number(i,0)+a1,Number(i,1)-Number(i,0)+b1)-R::lbeta(a1,b1)+
         R::lbeta(Number(i,2)+a0,Number(i,3)-Number(i,2)+b0)-R::lbeta(a0,b0);
  }
  
  return(out);
  }