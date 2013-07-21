#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector seqC(double start, double end, int length)
{
  NumericVector out(length);
  for (int i=0;i<length;i++)
  {
    out[i]=start + i*(end-start)/(length-1);
  }
  return out;
}