


 
#include <Rcpp.h>
using namespace Rcpp;
 
// [[Rcpp::export]]
NumericVector EstS_C(NumericMatrix Comb1,NumericMatrix Comb0, NumericMatrix Number,
int kq,int levelC,int N) 
  {
    int nrow1=Comb1.nrow(),nrow0=Comb0.nrow();
    int nrowC=Number.nrow();
    NumericVector out(nrowC);
    
    
    for (int i=0;i<nrowC;i++)
    {
      double sc1=0,sc0=0;
      for (int j=0;j<nrow1;j++)
      {
        sc1+= (Comb1(j,0)+Number(i,0))/(Number(i,1) + Comb1(j,0)+Comb1(j,1))*Comb1(j,2);
      }
      for (int k=0;k<nrow0;k++)
      {
        sc0+= (Comb0(k,0)+Number(i,2))/(Number(i,3) + Comb0(k,0)+Comb0(k,1))*Comb0(k,2);
      }
      out[i] = (kq+Number(i,4))/(levelC*kq+N)*(sc1-sc0);
      /**out[i]=sc1;out[i] = (kq+Number(i,4))/(levelC*kq+N)*(sc1-sc0);*/
    }

return(out);
}


/** R
benchmark(
out <- EstS_C(as.matrix(t.comb1),as.matrix(t.comb0),
Number=as.matrix(t.N..Number),kq=t.k_q,levelC=t.level.con,N=t.no.ob),

temp.new<- (
    apply(t.N..Number,1,function(x) {
      
      (t.k_q+x[5])/(t.level.con*t.k_q+t.no.ob)*
      (Exp(as.numeric(x),t.comb1,1)-
      Exp(as.numeric(x),t.comb0,0))
      
    })
  ),replications=3
  )

*/
