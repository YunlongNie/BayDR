#' Estimate the true causal effect
#'
#' This function estimate the true causal effect
#'
#' @param K number of block
#' @param r number of confounders in each block
#' @param rho correlation coef
#' @param gamma0,gamma1,gamma2
#' @param lambda1,lambda2
#' @param CSample
#' @export
#' @examples
#' \dontrun{
#' CausalEst <- function(K=2,r=6,rho=0.3,gamma0=0.2,gamma1=1,gamma2=1,lambda1=1,lambda2=0,CSample)  
#' }



CausalEst <- function(K=2,r=6,rho=0.3,gamma0=0.2,gamma1=1,gamma2=1,lambda1=1,lambda2=0,mc.error=0.0001) {

Z <- rnorm(10000,0,1)
Mean34 <- mean((1-pnorm((rho^(-1)-1)^(-0.5)*Z,0,1))^2)

Causal <-CauEstC(mc.error,gamma0=gamma0,gamma1=gamma1,gamma2=gamma2,lambda1=lambda1,lambda2=lambda2,r=r,K=K,mean34=Mean34)
Causal
}

#CausalEst(mc.error=0.001,gamma0=0.2,gamma1=1,gamma2=1,lambda1=1,lambda2=0,r=6,K=2)