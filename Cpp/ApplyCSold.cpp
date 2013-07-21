 /**
 DeltaS <-  apply(N..Number,1,function(x) {
    
    (k_q+x[5])/(level.con*k_q+no.ob)*((a1+x[1])/(a1+b1+x[2])-(a0+x[3])/(a0+b0+x[4]))
    
  })
  */
  
  
  
#include <Rcpp.h>
using namespace Rcpp;
 // [[Rcpp::export]]
NumericVector ApplyCSold(NumericMatrix Number, double kq, 
double levelCon, double noObs, double a1, double b1, double a0, double b0) 
  {
  int nRow=Number.nrow();
  NumericVector out(nRow);
  for (int i=0;i<nRow;++i)
  {
    
      out[i]= (kq + Number(i,4))/(levelCon*kq+noObs)*
      ((a1+Number(i,0))/(a1+b1+Number(i,1))-(a0+Number(i,2))/(a0+b0+Number(i,3)));
      
  }
  
  return(out);
  }