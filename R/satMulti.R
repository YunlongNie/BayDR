#' Compute the S1 P1 and B1 and S2 P2 B2 in fact P2=P1
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
#' @param addBin(10) can be a vector the default is 10
#' @return a list of saturated, bayesian estimate and paramatric estimate
#' @seealso AllEst_C
#' @export
#' @importFrom plyr dlply
#' @importFrom boot inv.logit
#' @importFrom LaplacesDemon as.inverse as.symmetric.matrix
#' @importFrom mnormt dmnorm
#' @importFrom speedglm speedglm
#' @examples 
#' \dontrun{
#' data(sampleDat)
#' Y=sample.dataset$Y
#' X=sample.dataset$X
#' C=sample.dataset$C
#' bay.est.new(Y,X,C)
#' }
sat2 <- function(Y,X,C,k_q,mean,BinMean,con,BinCon,addBin,kappa,beta,Mc.error,liketype,Dat)
{
  if (missing(Y)|missing(X)|missing(C)) {Y=Dat$Y;X=Dat$X;C=Dat[,paste0("C",1:(ncol(Dat)-2))]}
  if (missing(Dat)&(missing(Y)|missing(X)|missing(C))) stop("Data entry wrong Y or X or C or Dat is missing") 
  if (missing(k_q)) k_q=1
  if (missing(mean)) mean=c(0.01,0.99)
  if (missing(con)) con=c(1,20)
  if (missing(BinMean)) BinMean=20
  if (missing(BinCon)) BinCon=20
  if (missing(beta)) beta=30^2*diag(rep(1,ncol(C)+2))
  if (missing(Mc.error)) Mc.error=0.001
  if (missing(kappa)) kappa=0.5
  if (missing(addBin)) addBin=c(0,10)
  if (missing(liketype)) liketype="bernoulli"
  
  temp.p = para.est(Y,X,C,beta=beta,k_q=k_q,Mc.error=Mc.error)
  estimate.p = temp.p$P
  likelihood.p=temp.p$likelihood
  temp.s0 <- sat.est.new(Y,X,C,k_q,mean,BinMean,con,BinCon,addBin=addBin[1],liketype=liketype)
  likelihood0 = temp.s0$likelihood.new
  estimate.s0 = temp.s0$est.new
  
  ratio.likelihood0<- exp(likelihood0 - likelihood.p)
  w0 <- kappa/(kappa+(1-kappa)*ratio.likelihood0)
  estimate.b0 <- w0*estimate.p+(1-w0)*estimate.s0 # bayesian estimate for hirarchical version
  
  temp.s1 <- sat.est.new(Y,X,C,k_q,mean,BinMean,con,BinCon,addBin=addBin[2],liketype=liketype)
  likelihood1 = temp.s1$likelihood.new
  estimate.s1 = temp.s1$est.new
  
  ratio.likelihood1<- exp(likelihood1 - likelihood.p)
  w1 <- kappa/(kappa+(1-kappa)*ratio.likelihood1)
  estimate.b1 <- w1*estimate.p+(1-w1)*estimate.s1 # bayesian estimate for hirarchical version
  
  
  return(
    list(est=c(S0=estimate.s0, 
               P=estimate.p,
               B0=estimate.b0,
               S1=estimate.s1,
               B1=estimate.b1
                ),
         likelihood=c(l.p=likelihood.p,l.S0=likelihood0,l.S1=likelihood1)    
    )
  )  
}
