\name{bay.est.old}
\alias{bay.est.old}
\title{Compute the Sold P and Bold}
\usage{
  bay.est.old(Y, X, C, k_q = 1, prior0, prior1, kappa,
    beta, Mc.error = 0.001, Dat)
}
\arguments{
  \item{Y}{response varaible}

  \item{X}{exposure}

  \item{C}{confounder matrix}

  \item{k_q(1)}{prior parameter of the Dirichlet
  distribution for Pr(C=c)}

  \item{prior0(1),prior1(1)}{prior parameter of the
  saturated old estimate}

  \item{kappa(0.5)}{the prior weight when calculating the
  bayesian esimate}

  \item{beta}{the prior of the paramatric estimate}

  \item{Mc.error(0.001),}{controls the Mc.error when
  estimating the parametric estimate}
}
\value{
  a list of paramatric estimate, bayesian estimate and
  saturated old estimate
}
\description{
  This function computes the old bayesian estimator
}
\examples{
\dontrun{
 bay.est.old <- function(Y,X,C,k_q=1,prior0,prior1,kappa,beta,Mc.error=0.001)
}
}
\seealso{
  \code{\link{Allest_C}}
}

