#' Compute the Sold P and Bold
#'
#' This function computes the old bayesian estimator
#'
#' @param Y response varaible
#' @param X exposure
#' @param C confounder matrix
#' @param k_q(1) prior parameter of the Dirichlet distribution for Pr(C=c) 
#' @seealso \code{\link{Allest_C}}
#' @param prior0(1),prior1(1) prior parameter of the saturated old estimate
#' @param kappa(0.5) the prior weight when calculating the bayesian esimate
#' @param beta the prior of the paramatric estimate
#' @param Mc.error(0.001), controls the Mc.error when estimating the parametric estimate
#' @return a list of paramatric estimate, bayesian estimate and saturated old estimate
#' @export
#' @importFrom plyr dlply
#' @importFrom boot inv.logit
#' @importFrom LaplacesDemon as.inverse
#' @importFrom mnormt dmnorm
#' @examples
#' \dontrun{
#'  bay.est.old <- function(Y,X,C,k_q=1,prior0,prior1,kappa,beta,Mc.error=0.001)
#' }



bay.est.old <- function(Y,X,C,k_q=1,prior0,prior1,kappa,beta,Mc.error=0.001,Dat)
  {
  
  if (missing(Y)|missing(X)|missing(C)) {Y=Dat$Y;X=Dat$X;C=Dat[,paste0("C",1:(ncol(Dat)-2))]}
  if (missing(Dat)&(missing(Y)|missing(X)|missing(C))) stop("Data entry wrong Y or X or C or Dat is missing") 
  if (missing(k_q)) k_q=1
  if (missing(prior0)) {prior0=c(0.5,2)}
  if (missing(prior1)) {prior1=c(0.5,2)}
  if (missing(kappa)) kappa=0.5
  if (missing(beta)) beta=30^2*diag(rep(1,ncol(C)+2))
  if (missing(Mc.error)) Mc.error=0.001

  temp.p = para.est(Y,X,C,beta=beta,k_q=k_q,Mc.error=Mc.error)
  estimate.p = temp.p$P
  likelihood.p=temp.p$likelihood
  temp.s.old = sat.est.old(Y,X,C,k_q=k_q,prior0=prior0,prior1=prior1)
  likelihood.old=temp.s.old$likelihood.old
  estimate.s.old = temp.s.old$est.old
  
  ratio.likelihood.old <- exp(likelihood.old - likelihood.p)
  w.old <- kappa/(kappa+(1-kappa)*ratio.likelihood.old)
  estimate.b.old <- w.old*estimate.p + (1-w.old)*estimate.s.old # bayesian estimate for simple version
  
  output=list(est=c(Sold=estimate.s.old, # without smoothing for the prior
                    P=estimate.p,
                    Bold=estimate.b.old),
              likelihood=c(l.old=likelihood.old,l.p=likelihood.p)
  )
  
  return(
    output
  )  
}
