#' Compute the Sold.
#'
#' This function computes the old saturated estimators
#'
#' @param Y response varaible
#' @param X exposure
#' @param C confounder matrix
#' @seealso \code{\link{Allest_C.R}}
#' @param prior0,prior1 prior parameter of the saturated old estimate
#' @param kappa the prior weight when calculating the bayesian esimate
#' @param beta the prior of the paramatric estimate
#' @param Mc.error(0.001), controls the Mc.error when estimating the parametric estimate
#' @return a list of  saturated old estimate and its likelihood
#' @export
#' @importFrom plyr dlply
#' @importFrom boot inv.logit
#' @importFrom LaplacesDemon as.inverse
#' @importFrom mnormt dmnorm
#' @examples
#' \dontrun{
#' AllEstD_C <- function(Y,X,C,k_q=1,prior0=c(1,1),prior1=c(1,1))
#' }
sat.est.old <- function(Y,X,C,k_q,prior0,prior1,Dat) {
  if (missing(Y)|missing(X)|missing(C)) {Y=Dat$Y;X=Dat$X;C=Dat[,paste0("C",1:(ncol(Dat)-2))]}
  if (missing(Dat)&(missing(Y)|missing(X)|missing(C))) stop("Data entry wrong Y or X or C or Dat is missing") 
  if (missing(k_q)) k_q=1
  if (missing(prior0)) 
  {
    a0=b0=1
  }
  else {a0=prior0[1]*prior0[2];b0=prior0[2]-a0}
  
  if (missing(prior1)) 
  {
    a1=b1=1
  }
  else {a1=prior1[1]*prior1[2];b1=prior1[2]-a1}
  
  Dat <- as.data.frame(cbind(C,X,Y)) # combine the C confounder, X exposure and Y response into a dataset
  no.confounder = ncol(C) # number of confounders
  level.con <- 2^no.confounder # levels of confounder's combination
  no.ob <- nrow(Dat) # number of observations
  names(Dat) <- c(paste("C",1:no.confounder,sep=""),"X","Y")
  Dat$X <- factor(Dat$X)
  Dat$Y <- factor(Dat$Y)
  temp.list <- dlply(Dat,paste("C",1:no.confounder,sep=""),
                     function(x) {
                       temp.table <- table(x[,(no.confounder+1):(no.confounder+2)])
                       as.numeric(temp.table)
                     }
  )
  UqC <- attr(temp.list,"split_labels")
  N..Number <- as.data.frame(do.call(rbind,lapply(temp.list,function(x) c(x[4],x[4]+x[2],x[3],x[3]+x[1],sum(x)))))
  colnames(N..Number) <- c('C11','C1.','C01','C0.','Number') # matrix for n_cxx 
  rownames(N..Number) <- 1:nrow(N..Number)
  N..Number0 <- N..Number[,3:4];N..Number1 <- N..Number[,1:2] # for C=0 & for C=1
  n.full <- nrow(UqC)                 ###  # non-empty cells               
  n.mpty <- level.con-n.full               ###  # empty cells 
  

  DeltaS = ApplyCSold(Number=as.matrix(N..Number),kq=k_q,levelCon=level.con,
                 noObs=no.ob,a1=a1,a0=a0,b0=b0,b1=b1)
  estimate.s.old<- sum(DeltaS) + n.mpty*k_q/(level.con*k_q + no.ob)*(a1/(a1+b1)-a0/(a0+b0)) 

   
  likelihood.old = sum(ApplyCSoldhood(Number=as.matrix(N..Number),a1=a1,b1=b1,a0=a0,b0=b0))
  
  return(
    list(est.old=estimate.s.old, # without smoothing for the prior
         likelihood.old=likelihood.old 
         )
  )  
}
