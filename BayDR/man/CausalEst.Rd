\name{CausalEst}
\alias{CausalEst}
\title{Estimate the true causal effect}
\usage{
  CausalEst(K = 2, r = 6, rho = 0.3, gamma0 = 0.2,
    gamma1 = 1, gamma2 = 1, lambda1 = 1, lambda2 = 0,
    mc.error = 1e-04)
}
\arguments{
  \item{K}{number of block}

  \item{r}{number of confounders in each block}

  \item{rho}{correlation coef}

  \item{gamma0,gamma1,gamma2}{}

  \item{lambda1,lambda2}{}

  \item{CSample}{}
}
\description{
  This function estimate the true causal effect
}
\examples{
\dontrun{
CausalEst <- function(K=2,r=6,rho=0.3,gamma0=0.2,gamma1=1,gamma2=1,lambda1=1,lambda2=0,CSample)
}
}

