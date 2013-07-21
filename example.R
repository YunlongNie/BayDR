# install the package
install.packages(file.choose(), repos=NULL)
# then choose the zip file


# a simple example
 
library(BayDR)

# sampleDat is a sample dataset in the package.
data(sampleDat)

Y=sample.dataset$Y
X=sample.dataset$X
C=sample.dataset$C
(Causal=sample.dataset$CausalEst)
## Five Estimates
AllEs_C(Y,X,C)

## bayes new
bay.est.new(Y,X,C)

## bayes old
bay.est.old(Y,X,C)


## saturated new
sat.est.new(Y,X,C)

## saturated old
sat.est.old(Y,X,C)

## function of the package
ls("package:BayDR")
?AllEs_C # for help

