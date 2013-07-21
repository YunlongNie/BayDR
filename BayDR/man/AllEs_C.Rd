\name{AllEs_C}
\alias{AllEs_C}
\title{Five estimators}
\usage{
  AllEs_C(Y, X, C, k_q, a0, b0, a1, b1, mean, BinMean, con,
    BinCon, addBin, kappa, beta, Mc.error, Nhead)
}
\arguments{
  \item{Y}{binary repsonse vector length of n}

  \item{X}{binary exposure vector same length with Y}

  \item{C}{binary confounder matrix n*r, r the number of
  counfouders}

  \item{k_q(1)}{prior parameter of the Dirichlet
  distribution for Pr(C=c)}

  \item{a0,b0,a1,b1(1)}{prior parameter of the saturated
  old estimate}

  \item{mean}{default (0.01,0.99) con (1,20) the range of
  the grid for the hyperparameter of the saturated new
  estimate}

  \item{BinMean(20),BinCon(20)}{n_1 control the number of
  new-added points}

  \item{kappa(0.5)}{the prior weight when calculating the
  bayesian esimate}

  \item{beta}{30^2*diag(rep(1,ncol(C)+2)),the prior of the
  paramatric estimate}

  \item{Mc.error(0.001)}{controls the Mc.error when
  estimating the parametric estimate}

  \item{Nhead(10)}{t The number of points in the first
  where new points are attached to}

  \item{addBin(10)}{n_2}
}
\value{
  a list of five estimates
}
\description{
  This function gives P snew sold bold bnew
}
\examples{
\dontrun{
data(sampleDat)
Y=sample.dataset$Y
X=sample.dataset$X
C=sample.dataset$C
est=AllEs_C(Y,X,C,addBin=10)
est$likelihood # l.old        l.p      l.new


$est
Sold        Snew           P        Bold        Bnew
0.009647389 0.388687875 0.396289809 0.009647389 0.388687875

$likelihood
l.old        l.p      l.new
-699.9814 -1182.9841  -608.2210

}
}

