#' Simulating the data
#'
#' This function simulates the data
#'
#' @param n number of ob
#' @param K number of block
#' @param r number of confounders in each block
#' @param rho correlation coef
#' @importFrom mnormt rmnorm 
#' @export
#' @examples
#' \dontrun{
#'  GenDat2 <- function(n=1000,K=2,r=6,rho=0.3,gamma0=0.2,gamma1=0,gamma2=0,
#' lambda1=0,lambda2=0) 
#' }



GenDat2 <- function(n=1000,K=2,r=6,rho=0.3,gamma0=0.2,gamma1=0,gamma2=0,
                   lambda1=0,lambda2=0) {
 
  Sigma <- diag(rep(1,r))
  Sigma[Sigma==0] <- rho
  C<- do.call(cbind,lapply(1:K,function(x)
  {
    temp <-  rmnorm(n,mean=rep(0,r),varcov = Sigma)
    temp[temp>0] <- 1
    temp[temp!=1] <- 0
    temp
  })
  )
  Seq <- 1:(6*K)
  No123 <- Seq[((Seq)%%6)%in%1:3] # get the seq 1,2,3,7,8,9....
  No1 <- seq(1,(r*K),by=r);No2 =seq(2,(r*K),by=r);No34=c(seq(3,(r*K),by=r),seq(4,(r*K),by=r))
  Z <- rnorm(10000,0,1)
  Mean34 <- mean((1-pnorm((rho^(-1)-1)^(-0.5)*Z,0,1))^2)
  CDat <- do.call(cbind,DatDen_X_C(C,gamma0,gamma1,gamma2,
                                   lambda1,lambda2,K,No123,No1,No2,No34,Mean34))
  Dat <- as.data.frame(cbind(C,CDat))
  names(Dat) <- c(paste("C",1:(r*K),sep=""),"X","Y")
  return(Dat)
}
