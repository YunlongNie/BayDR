\name{sat.est.old}
\alias{sat.est.old}
\title{Compute the Sold.}
\usage{
  sat.est.old(Y, X, C, k_q, prior0, prior1, Dat)
}
\arguments{
  \item{Y}{response varaible}

  \item{X}{exposure}

  \item{C}{confounder matrix}

  \item{prior0,prior1}{prior parameter of the saturated old
  estimate}

  \item{kappa}{the prior weight when calculating the
  bayesian esimate}

  \item{beta}{the prior of the paramatric estimate}

  \item{Mc.error(0.001),}{controls the Mc.error when
  estimating the parametric estimate}
}
\value{
  a list of saturated old estimate and its likelihood
}
\description{
  This function computes the old saturated estimators
}
\examples{
\dontrun{
AllEstD_C <- function(Y,X,C,k_q=1,prior0=c(1,1),prior1=c(1,1))
}
}
\seealso{
  \code{\link{Allest_C.R}}
}

