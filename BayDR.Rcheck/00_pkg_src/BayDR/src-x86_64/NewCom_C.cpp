#include <Rcpp.h>
using namespace Rcpp;
NumericVector seqC(double start, double end, int length)
{
  NumericVector out(length);
  if (length > 1) {
  for (int i=0;i<length;i++)
  {
    out[i]=start + i*(end-start)/(length-1);
  }
  } else {
    out=start;
  }
  return out;
}

// [[Rcpp::export]]
List NewComb_C(NumericMatrix Add, NumericVector weight,int Addbin,double intMean,double intCon) 
  {
   int nrow=Add.nrow(); 
   NumericMatrix Addtemp(nrow,2);
   List out(nrow);
   Environment base("package:base");
   Function cbind = base["cbind"];
   Function expandGrid = base["expand.grid"];
   Function asMatrix = base["as.matrix"];
  
   for (int i=0;i<nrow;++i)
   {
     Addtemp(i,0) = Add(i,0)/(Add(i,0)+Add(i,1));
     Addtemp(i,1)= Add(i,0)+Add(i,1);
   }
   
   NumericMatrix Addbox = cbind(Add,Addtemp);
   
   for (int i=0;i<nrow;i++)
   {
     NumericVector addmean=seqC(Addbox(i,2)-intMean/2,Addbox(i,2)+intMean/2,Addbin);
     NumericVector addcon=seqC(Addbox(i,3)-intCon/2,Addbox(i,3)+intCon/2,Addbin);
     
      NumericMatrix temp = asMatrix(expandGrid(addmean,addcon));
      NumericMatrix temp2(temp.nrow(),3);
     for (int j=0;j< temp.nrow();++j)
     {
       temp2(j,0)= temp(j,0)*temp(j,1);
       temp2(j,1)= temp(j,1)-temp2(j,0);
       temp2(j,2)= weight[i]/temp.nrow();
     }
     out[i]=cbind(temp,temp2);
    
   }
  
    return(out);
  }







 
