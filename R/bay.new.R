#' Compute the Snew P and Bnew
#'
#' This function computes the new bayesian estimator
#'
#' @param Y binary repsonse vector length of n
#' @param X binary exposure vector same length with Y
#' @param C binary confounder matrix n*r, r the number of counfouders
#' @param k_q(1) prior parameter of the Dirichlet distribution for Pr(C=c)
#' @param mean(0.01,0.99),con(1,20) the range of the grid for the hyperparameter of the saturated new estimate
#' @param BinMean(20),BinCon(20) control the number of new-added points 
#' @param kappa(0.5) the prior weight when calculating the bayesian esimate
#' @param beta the prior of the paramatric estimate
#' @param Mc.error(0.001), controls the Mc.error when estimating the parametric estimate
#' @param addBin(10) n_2
#' @return a list of saturated, bayesian estimate and paramatric estimate
#' @seealso AllEst_C
#' @export
#' @importFrom plyr dlply
#' @importFrom boot inv.logit
#' @importFrom LaplacesDemon as.inverse
#' @importFrom mnormt dmnorm
#' @examples 
#' \dontrun{
#' data(sampleDat)
#' Y=sample.dataset$Y
#' X=sample.dataset$X
#' C=sample.dataset$C
#' bay.est.new(Y,X,C)
#' }
bay.est.new <- function(Y,X,C,k_q,mean,BinMean,con,BinCon,addBin,kappa,beta,Mc.error,liketype,Dat)
  {
  if (missing(Y)|missing(X)|missing(C)) {Y=Dat$Y;X=Dat$X;C=Dat[,paste0("C",1:(ncol(Dat)-2))]}
  if (missing(Dat)&(missing(Y)|missing(X)|missing(C))) stop("Data entry wrong Y or X or C or Dat is missing") 
  if (missing(k_q)) k_q=1
  if (missing(mean)) mean=c(0.01,0.99)
  if (missing(con)) con=c(1,20)
  if (missing(BinMean)) BinMean=20
  if (missing(BinCon)) BinCon=20
  if (missing(beta)) {beta=30^2*diag(rep(1,ncol(C)+2))}
  if (missing(Mc.error)) Mc.error=0.001
  if (missing(kappa)) kappa=0.5
  if (missing(addBin)) addBin=10
  if (missing(liketype)) liketype="bernoulli"
  
  temp.p = para.est(Y,X,C,beta=beta,k_q=k_q,Mc.error=Mc.error)
  estimate.p = temp.p$P
  likelihood.p=temp.p$likelihood
  temp.s.new <- sat.est.new(Y,X,C,k_q,mean,BinMean,con,BinCon,addBin,liketype=liketype)
  likelihood.new = temp.s.new$likelihood.new
  estimate.s.new = temp.s.new$est.new
  ratio.likelihood.new<- exp(likelihood.new - likelihood.p)
  w.new <- kappa/(kappa+(1-kappa)*ratio.likelihood.new)
  estimate.b.new <- w.new*estimate.p+(1-w.new)*estimate.s.new # bayesian estimate for hirarchical version
  
  
  return(
    list(est=c(Snew=estimate.s.new, 
         P=estimate.p,
         Bnew=estimate.b.new),
         likelihood=c(l.p=likelihood.p,l.Snew=likelihood.new)    
         )
  )  
}
