#' Compute the P.
#'
#' This function computes paramatrix estimator.
#'
#' @param Y response varaible
#' @param X exposure
#' @param C confounder matrix
#' @param beta prior matrix for the paramatric model
#' @param Mc.error default 0.001
#' @importFrom plyr dlply
#' @importFrom boot inv.logit
#' @importFrom LaplacesDemon as.inverse
#' @importFrom mnormt dmnorm
#' @export
#' @examples
#' \dontrun{
#' para.est <- function(Y,X,C,beta,k_q=1,Mc.error=0.001)
#' }
para.est <- function(Y,X,C,beta,k_q=1,Mc.error=0.001,Dat) {
  
  
  if (missing(Y)|missing(X)|missing(C)) {Y=Dat$Y;X=Dat$X;C=Dat[,paste0("C",1:(ncol(Dat)-2))]}
  if (missing(Dat)&(missing(Y)|missing(X)|missing(C))) stop("Data entry wrong Y or X or C or Dat is missing") 
  if (missing(beta)) {beta=30^2*diag(rep(1,ncol(C)+2))}
  if (missing(k_q)) {k_q=1}
  if (missing(Mc.error)) {Mc.error=0.001}
  Dat <- cbind(C,X,Y) # combine the C confounder, X exposure and Y response into a dataset
  no.confounder = ncol(C) # number of confounders
  level.con <- 2^no.confounder # levels of confounder's combination
  no.ob <- nrow(Dat) # number of observations
  
  Dat <- as.data.frame(Dat)
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
  
  
  n.full <- nrow(UqC)                 ###  # non-empty cells               
  n.mpty <- level.con-n.full               ###  # empty cells 
  
  Dat$X <-X # otherwise Hatm1C1 might have some error. 
  Dat$Y <-Y
  outModel <- glm(Y~(.),data=Dat,family='binomial')
  Hatm1C1 <- inv.logit(as.matrix(cbind(1,UqC,1))%*%as.matrix(coef(outModel)))
  Hatm0C0 <- inv.logit(as.matrix(cbind(1,UqC,0))%*%as.matrix(coef(outModel)))
  
  PrC <- apply(N..Number,1,function(x) {
    
    (k_q+x[5])/(level.con*k_q+no.ob)
    
  })
  if (n.mpty==0) { estimate.p <- sum(PrC*(Hatm1C1-Hatm0C0)) } else { 
    CSample = sample.empty(N=(0.5/Mc.error)^2,UqC=as.matrix(UqC))
    Hatm1C1full <- inv.logit(as.matrix(cbind(1,CSample,1))%*%as.matrix(coef(outModel)))
    Hatm0C0full <- inv.logit(as.matrix(cbind(1,CSample,0))%*%as.matrix(coef(outModel)))
    estimate.p <- mean(Hatm1C1full-Hatm0C0full)*k_q*n.mpty/(level.con*k_q+no.ob) + sum(PrC*(Hatm1C1-Hatm0C0))
  } 
  
  
  
  theta <- matrix(nrow=no.confounder+2, as.numeric(coef(outModel))) # MLE \beta
  H <- as.inverse(vcov(outModel)) 
  Sigma_star <- as.inverse(as.inverse(beta)+H)
  beta_star <- Sigma_star%*%H%*%theta 
  phi <- dmnorm(x=as.numeric(beta_star),mean=rep(0,no.confounder+2),varcov= beta)/
    dmnorm(x=as.numeric(beta_star),mean=as.numeric(beta_star),varcov = Sigma_star)
  Dat$X <- as.numeric(Dat$X)
  temp= as.matrix(cbind(inter=1,(Dat[,-no.confounder-2])))
  fy <- inv.logit(temp%*%as.matrix(beta_star))
  likelihood.p <- sum(log(fy^Y*(1-fy)^(1-Y))) + log(phi) # likelihood of parametric estimate 
  output <- list(P=estimate.p,likelihood=likelihood.p)
  
  
  return(
    output 
  )  
}
