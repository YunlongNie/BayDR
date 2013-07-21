\name{para.est}
\alias{para.est}
\title{Compute the P.}
\usage{
  para.est(Y, X, C, beta, k_q = 1, Mc.error = 0.001, Dat)
}
\arguments{
  \item{Y}{response varaible}

  \item{X}{exposure}

  \item{C}{confounder matrix}

  \item{beta}{prior matrix for the paramatric model}

  \item{Mc.error}{default 0.001}
}
\description{
  This function computes paramatrix estimator.
}
\examples{
\dontrun{
para.est <- function(Y,X,C,beta,k_q=1,Mc.error=0.001)
}
}

