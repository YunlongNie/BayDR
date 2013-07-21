\name{dcom}
\alias{dcom}
\title{Decompose the confounder matrix}
\usage{
  dcom(C, X, Y, Dat)
}
\arguments{
  \item{Y}{binary repsonse vector length of n}

  \item{X}{binary exposure vector same length with Y}

  \item{C}{binary confounder matrix n*r, r the number of
  counfouders}

  \item{Dat}{can be missing if given Y,X,C}
}
\value{
  a list UqC of the unique levels of C, N..Number,
  N..Number0, N..Number1
}
\description{
  This function decomposes the confounder matrix.
}
\examples{
\dontrun{
data(sampleDat)
Y=sample.dataset$Y
X=sample.dataset$X
C=sample.dataset$C
dcom(Y,X,C)
}
}

