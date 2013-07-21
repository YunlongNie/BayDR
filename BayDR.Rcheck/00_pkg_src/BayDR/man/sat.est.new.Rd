\name{sat.est.new}
\alias{sat.est.new}
\title{Compute the Snew}
\usage{
  sat.est.new(Y, X, C, k_q, mean, BinMean, con, BinCon,
    addBin, Nhead, liketype, Dat)
}
\arguments{
  \item{Y}{response varaible}

  \item{X}{exposure}

  \item{C}{confounder matrix}

  \item{mean}{default (0.01,0.99) con (1,20) the range of
  the grid for the hyperparameter of the saturated new
  estimate}

  \item{BinMean(20),BinCon(20)}{control the number of
  new-added points}

  \item{addBin(10)}{}

  \item{kappa}{the prior weight when calculating the
  bayesian esimate}

  \item{beta}{the prior of the paramatric estimate}

  \item{liketype}{type of likelihood for saturated model
  bernoulli or binomial}
}
\value{
  a list of saturated estimate and its likelihood
}
\description{
  This function computes the new saturated estimator
}
\examples{
\dontrun{
data(sampleDat)
Y=sample.dataset$Y
X=sample.dataset$X
C=sample.dataset$C
 sat.est.new <- function(Y,X,C,k_q,mean=c(0.01,0.99),BinMean = 20,
con = c(1,20),BinCon = 20)
}
}
\seealso{
  \code{\link{Allest_C}}
}

