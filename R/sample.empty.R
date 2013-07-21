#' Sample the empty cells data
#'
#' This function generate the empty-cell samples
#' 
#' 
#' @param UqC 
#' @param N
#' @export
#' @examples
#' \dontrun{
#' sample.empty<- function(UqC,N)
#' }


sample.empty<- function(UqC,N,colnumber){
  if(missing(UqC)) {ncol=colnumber;ext=NULL} else 
  {
    ncol <- ncol(UqC)
    ext <- apply(UqC,1,function (x) 
    {
      sum(x*2^((ncol-1):0))
    }
    )
    
  }
  temp.N=0
  sample.now <- NULL
  out <- NULL
  while (temp.N < N)
  {
    temp <-  matrix(sample(2,10000*ncol,replace=TRUE)-1,ncol=ncol)
    
    if (is.null(ext)) out <- rbind(out,temp) else 
    {
      sample.now <-apply(temp,1,function(x) sum(x*2^((ncol-1):0)))
      out <- rbind(out,temp[!sample.now%in%ext,])
    }
    temp.N =nrow(out)
  }
  out
}

# /**
#   CSample = sample.empty(N=(0.5/Mc.error)^2,UqC=as.matrix(UqC))}
# Hatm1C1full <- inv.logit(as.matrix(cbind(1,CSample,1))%*%as.matrix(coef(outModel)))
# Hatm0C0full <- inv.logit(as.matrix(cbind(1,CSample,0))%*%as.matrix(coef(outModel)))
# 
# estimate.p <- mean(Hatm1C1full-Hatm0C0full)*k_q*n.mpty/(level.con*k_q+no.ob) + sum(PrC*(Hatm1C1-Hatm0C0))
# NumericVector out(1);
# NumericVector out2(1);
# NumericVector num2 = NumericVector::create(0.0,1.0);
# 
# NumericVector out3 = NumericVector::create(1,2,3,4);
# NumericVector xx = NumericVector::create(1,2,3,4);
# out2 = xx*out3;
# return(out2);
# */
