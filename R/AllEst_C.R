#' Five estimators 
#'
#' This function gives P snew sold bold bnew
#' 
#' 
#' @param Y binary repsonse vector length of n
#' @param X binary exposure vector same length with Y
#' @param C binary confounder matrix n*r, r the number of counfouders
#' @param k_q(1) prior parameter of the Dirichlet distribution for Pr(C=c)
#' @param a0,b0,a1,b1(1) prior parameter of the saturated old estimate
#' @param mean default (0.01,0.99) con (1,20) the range of the grid for the hyperparameter of the saturated new estimate
#' @param BinMean(20),BinCon(20) n_1 control the number of new-added points 
#' @param kappa(0.5) the prior weight when calculating the bayesian esimate
#' @param beta 30^2*diag(rep(1,ncol(C)+2)),the prior of the paramatric estimate
#' @param Mc.error(0.001) controls the Mc.error when estimating the parametric estimate
#' @param Nhead(10) t The number of points in the first where new points are attached to 
#' @param addBin(10) n_2 
#' @return a list of five estimates 
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
#' est=AllEs_C(Y,X,C,addBin=10)
#  est$est  #Sold        Snew           P        Bold        Bnew 
#' est$likelihood # l.old        l.p      l.new 
#' 
#' 
#' $est
#' Sold        Snew           P        Bold        Bnew 
#' 0.009647389 0.388687875 0.396289809 0.009647389 0.388687875 
#' 
#' $likelihood
#' l.old        l.p      l.new 
#' -699.9814 -1182.9841  -608.2210 
#' 
#' }

AllEs_C <- function(Y,X,C,k_q,a0,b0,a1,b1,mean,BinMean,
                      con,BinCon,addBin,
                      kappa,
                      beta,Mc.error,Nhead) {
  if (missing(k_q)) k_q=1
  if (missing(mean)) mean=c(0.01,0.99)
  if (missing(a0)) a0=1
  if(missing(b0)) b0=1
  if (missing(a1)) a1=1
  if(missing(b1)) b1=1
  if (missing(con)) con=c(1,20)
  if (missing(BinMean)) BinMean=20
  if (missing(BinCon)) BinCon=20
  if (missing(beta)) {beta=30^2*diag(rep(1,ncol(C)+2))}
  if (missing(Mc.error)) Mc.error=0.001
  if (missing(kappa)) kappa=0.5
  if (missing(addBin)) addBin=10
  if (missing(Nhead)) Nhead=10
  #require(boot);require(LaplacesDemon);require(Rcpp);require(mnormt)
  
  y = seq(con[1],con[2],length=BinCon+2)[c(-1,-BinCon-2)]  # the prior con sequence :  the prior value that a_x b_x can take
  
  x = seq(mean[1],mean[2],length=BinMean+2)[c(-1,-BinMean-2)] # the prior mean sequence
  
  IntMean <- x[2]-x[1]; IntCon <- y[2]-y[1]  # get the interval length of y and x
  
  comb <- expand.grid(x = x, y = y); comb$a <- comb$x*comb$y;comb$b <- comb$y-comb$a # calculate the value of a and b prior
  
  Dat <- cbind(C,X,Y) # combine the C confounder, X exposure and Y response into a dataset
  
  
  no.confounder = ncol(C) # number of confounders
  
  
  level.con <- 2^no.confounder # levels of confounder's combination
  
  no.ob <- nrow(Dat) # number of observations
  
  
  #iff <- function(x) {Iff(x,dat=Dat,C=C)} # define a function to calculate the N..Number the value of n_cxx
  
  #UqC <- unique(C) # the combination of confounders that can be observed from the data 
  
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
  N..Number0 <- N..Number[,3:4];N..Number1 <- N..Number[,1:2] # for C=0 & for C=1
  
  
  
  
  n.full <- nrow(UqC)                 ###  # non-empty cells               
  n.mpty <- level.con-n.full               ###  # empty cells 
  
  DeltaS <-  apply(N..Number,1,function(x) {
    
    (k_q+x[5])/(level.con*k_q+no.ob)*((a1+x[1])/(a1+b1+x[2])-(a0+x[3])/(a0+b0+x[4]))
    
  })
  estimate.s.old<- sum(DeltaS) + n.mpty*k_q/(level.con*k_q + no.ob)*(a1/(a1+b1)-a0/(a0+b0)) # saturated estimate simple version two parts: non empty cells and empty cells 
  #####
  likelihood.old <- sum(apply(N..Number,1,function(x) {  # likehihood of estimate.s.old 
    (
      lgamma(x[1]+a1)+lgamma(x[2]-x[1]+b1)+lgamma(a1+b1)-
        lgamma(x[2]+a1+b1)-lgamma(a1)-lgamma(b1)+
        lgamma(x[3]+a0)+lgamma(x[4]-x[3]+b0)+lgamma(a0+b0)-
        lgamma(x[4]+a0+b0)-lgamma(a0)-lgamma(b0)
    )
  }
  )
  )
  
  comb1 <- comb0 <- data.frame(a=comb$a,b=comb$b,w=1/nrow(comb)) # set up the prior value and distribution for a_x b_x
  
  comb0$w <- Comb_C(as.matrix(comb0),as.matrix(N..Number0))## calculate the posterior weights for a_x and b_x 
  if (addBin>0){
    comb0 <- NewCom_C(comb0,Head=Nhead,IntMean=IntMean,IntCon=IntCon,AddBin=addBin) # add points to the high weighted posterior points
    Newcomb0=comb0
    comb0$w  <-Comb_C(as.matrix(comb0),as.matrix(N..Number0))## calulated the posterior weights again 
  } else {
    Newcomb0=data.frame(a=comb$a,b=comb$b,w=1/nrow(comb))
  }
  
  comb1$w <- Comb_C(as.matrix(comb1),as.matrix(N..Number1))## calculate the posterior weights for a_x and b_x 
  if (addBin>0){
    comb1<- NewCom_C(comb1,Head=Nhead,IntMean=IntMean,IntCon=IntCon,AddBin=addBin) # add points to the high weighted posterior points
    Newcomb1=comb1
    comb1$w  <-Comb_C(as.matrix(comb1),as.matrix(N..Number1))
  } else {
    Newcomb1=data.frame(a=comb$a,b=comb$b,w=1/nrow(comb))
  }
  
  ####  
  
  
  temp.new <- sum(EstS_C(as.matrix(comb1),as.matrix(comb0),
                         Number=as.matrix(N..Number),kq=k_q,levelC=level.con,N=no.ob))
 
  mty0 <- sum(apply(comb0,1,function(x) {
    x[1]/(x[1]+x[2])*x[3]
  })) 
  
  mty1 <- sum(apply(comb1,1,function(x) {
    x[1]/(x[1]+x[2])*x[3]
  }))
  
  estimate.s.new <- temp.new + n.mpty*k_q/(level.con*k_q + no.ob)*(mty1-mty0) # saturated estimate hierarchical version also two parts
  
  likelihood.new <- Like_C(as.matrix(Newcomb1),as.matrix(Newcomb0),
                           Number=as.matrix(N..Number))
  
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
  likelihood.p <- sum(log( fy^Y*(1-fy)^(1-Y))) + log(phi) # likelihood of parametric estimate 
  
  ratio.likelihood.old <- exp(likelihood.old - likelihood.p)
  w.old <- kappa/(kappa+(1-kappa)*ratio.likelihood.old)
  estimate.b.old <- w.old*estimate.p + (1-w.old)*estimate.s.old # bayesian estimate for simple version
  
  ratio.likelihood.new<- exp(likelihood.new - likelihood.p)
  w.new <- kappa/(kappa+(1-kappa)*ratio.likelihood.new)
  estimate.b.new <- w.new*estimate.p+(1-w.new)*estimate.s.new # bayesian estimate for hirarchical version
  
  output=list(est=c(Sold=estimate.s.old, # without smoothing for the prior
             Snew=estimate.s.new, 
             P=estimate.p,
             Bold=estimate.b.old,
             Bnew=estimate.b.new),
       likelihood=c(l.old=likelihood.old,l.p=likelihood.p,l.new=likelihood.new)
  )
  
  return(
    output
  )  
}
