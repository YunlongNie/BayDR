#' Compute the Snew
#'
#' This function computes the new saturated estimator
#' @seealso \code{\link{Allest_C}}
#' @param Y response varaible
#' @param X exposure
#' @param C confounder matrix
#' @param mean default (0.01,0.99) con (1,20) the range of the grid for the hyperparameter of the saturated new estimate
#' @param BinMean(20),BinCon(20) control the number of new-added points 
#' @param addBin(10)
#' @param kappa the prior weight when calculating the bayesian esimate
#' @param beta the prior of the paramatric estimate
#' @param liketype type of likelihood for saturated model bernoulli or binomial
#' @return a list of saturated estimate and its likelihood
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
#'  sat.est.new <- function(Y,X,C,k_q,mean=c(0.01,0.99),BinMean = 20,
#' con = c(1,20),BinCon = 20)
#' }
sat.est.new <- function(Y,X,C,k_q,mean,BinMean,con,BinCon,addBin,Nhead,liketype,Dat){
  
  if (missing(Y)|missing(X)|missing(C)) {Y=Dat$Y;X=Dat$X;C=Dat[,paste0("C",1:(ncol(Dat)-2))]}
  if (missing(Dat)&(missing(Y)|missing(X)|missing(C))) stop("Data entry wrong Y or X or C or Dat is missing") 
  if (missing(k_q)) k_q=1
  if (missing(mean)) mean=c(0.01,0.99)
  if (missing(con)) con=c(1,20)
  if (missing(BinMean)) BinMean=20
  if (missing(BinCon)) BinCon=20
  if (missing(addBin)) addBin=10
  if (missing(Nhead)) Nhead=10
  if (missing(liketype)) liketype="bernoulli"
  y = seq(con[1],con[2],length=BinCon+2)[c(-1,-BinCon-2)]  # the prior con sequence :  the prior value that a_x b_x can take
  
  x = seq(mean[1],mean[2],length=BinMean+2)[c(-1,-BinMean-2)] # the prior mean sequence
  
  IntMean <- x[2]-x[1]; IntCon <- y[2]-y[1]  # get the interval length of y and x
  
  comb <- expand.grid(x = x, y = y); comb$a <- comb$x*comb$y;comb$b <- comb$y-comb$a # calculate the value of a and b prior
  
  Dat <-  as.data.frame(cbind(C,X,Y)) # combine the C confounder, X exposure and Y response into a dataset
  
  
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
  
  comb1 <- comb0 <- data.frame(a=comb$a,b=comb$b,w=1/nrow(comb)) # set up the prior value and distribution for a_x b_x
  
  comb0$w <- Comb_C(as.matrix(comb0),as.matrix(N..Number0))## calculate the posterior weights for a_x and b_x 
  if (addBin>0){
  comb0 <- NewCom_C(comb0,Head=Nhead,IntMean=IntMean,IntCon=IntCon,AddBin=addBin) # add points to the high weighted posterior points
  Newcomb0=comb0
  comb0$w  <- Comb_C(as.matrix(comb0),as.matrix(N..Number0))## calulated the posterior weights again 
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
  
  temp.new <- sum(EstS_C(as.matrix(comb1),as.matrix(comb0),
                         Number=as.matrix(N..Number),kq=k_q,levelC=level.con,N=no.ob))
  
  mty0 <- sum(apply(comb0,1,function(x) {
    x[1]/(x[1]+x[2])*x[3]
  })) 
  
  mty1 <- sum(apply(comb1,1,function(x) {
    x[1]/(x[1]+x[2])*x[3]
  }))
  
  estimate.s.new <- temp.new + n.mpty*k_q/(level.con*k_q + no.ob)*(mty1-mty0) # saturated estimate hierarchical version also two parts
  if (liketype=="bernoulli")
    {
    likelihood.new <- Like_C(as.matrix(Newcomb1),as.matrix(Newcomb0),
                           Number=as.matrix(N..Number))
    } else {
      likelihood.new <- Like_CBin(as.matrix(Newcomb1),as.matrix(Newcomb0),
                               Number=as.matrix(N..Number))
                           }
  return(
    list(est.new=estimate.s.new, 
         likelihood.new=likelihood.new)
  )  
}


