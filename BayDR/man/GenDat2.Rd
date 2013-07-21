\name{GenDat2}
\alias{GenDat2}
\title{Simulating the data}
\usage{
  GenDat2(n = 1000, K = 2, r = 6, rho = 0.3, gamma0 = 0.2,
    gamma1 = 0, gamma2 = 0, lambda1 = 0, lambda2 = 0)
}
\arguments{
  \item{n}{number of ob}

  \item{K}{number of block}

  \item{r}{number of confounders in each block}

  \item{rho}{correlation coef}
}
\description{
  This function simulates the data
}
\examples{
\dontrun{
 GenDat2 <- function(n=1000,K=2,r=6,rho=0.3,gamma0=0.2,gamma1=0,gamma2=0,
lambda1=0,lambda2=0)
}
}

