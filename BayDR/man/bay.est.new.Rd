\name{bay.est.new}
\alias{bay.est.new}
\title{Compute the Snew P and Bnew}
\usage{
  bay.est.new(Y, X, C, k_q, mean, BinMean, con, BinCon,
    addBin, kappa, beta, Mc.error, liketype, Dat)
}
\arguments{
  \item{Y}{binary repsonse vector length of n}

  \item{X}{binary exposure vector same length with Y}

  \item{C}{binary confounder matrix n*r, r the number of
  counfouders}

  \item{k_q(1)}{prior parameter of the Dirichlet
  distribution for Pr(C=c)}

  \item{mean(0.01,0.99),con(1,20)}{the range of the grid
  for the hyperparameter of the saturated new estimate}

  \item{BinMean(20),BinCon(20)}{control the number of
  new-added points}

  \item{kappa(0.5)}{the prior weight when calculating the
  bayesian esimate}

  \item{beta}{the prior of the paramatric estimate}

  \item{Mc.error(0.001),}{controls the Mc.error when
  estimating the parametric estimate}

  \item{addBin(10)}{n_2}
}
\value{
  a list of saturated, bayesian estimate and paramatric
  estimate
}
\description{
  This function computes the new bayesian estimator
}
\examples{
\dontrun{
data(sampleDat)
Y=sample.dataset$Y
X=sample.dataset$X
C=sample.dataset$C
bay.est.new(Y,X,C)
}
}
\seealso{
  AllEst_C
}

