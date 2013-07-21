require(BayDR)

dat=GenDat2(n=1000,K=5,r=6)
head(dat)

Y=dat$Y
X=dat$X
C=dat[,1:30]
UqC=as.matrix(dcom(Y,X,C)$UqC)
seq=0:29
coef = rep(0.2,32)
sampleEmpty(mcerror=0.001,UqC=UqC,coef=coef,seq=seq)

head(dcom(Y,X,C)$N..Number)
res0=sat.est.new(Y,X,C,addBin=0)

bay.est.new(Dat=dat)
para.est(Dat=dat)
res0$est.new
w=res0$comb0[[3]]
sum(-sort(-w)[1:30])
colSums(res0$comb1)
with(res0$comb0,order(w))

sat.est.new(Y,X,C,addBin=10)$est.new
sat.est.new(Y,X,C,addBin=40)$est.new
sat.est.old(Y,X,C)
library(BayDR)

dat=GenDat2(phi1=2,K=2)
Y = dat$Y
X= dat$X
C=dat[,1:6]


bay.est.old(Y,X,C)

sat.est.new(Y,X,C,addBin=0)
sat.est.old(Y,X,C)

UqC=dcom(Dat=dat)$UqC

para.est(Y,X,C)
dat.sp <- dat[1:100,]
UqC.sp=dcom(Dat=dat.sp)$UqC
Y = dat.sp$Y
X= dat.sp$X
C=dat.sp[,1:6]

sat.est.new(Y,X,C)
sat.est.old(Y,X,C)


AllEs_C(Y,X,C)

CSample=sample.empty(UqC=dcom(dat)$UqC,N=100000)
CausalEst(phi1=2,CSample=CSample)

sample.dataset=list(Y=Y,X=X,C=C,phi=2,CausalEst=CausalEst(phi1=2,CSample=CSample),CSample=CSample)

save(sample.dataset,file="sampleDat.rda")

#
# sampleDat is a sample dataset in the package.
require(BayDR)
data(sampleDat)

Y=sample.dataset$Y
X=sample.dataset$X
C=sample.dataset$C
para.est(Y,X,C)
(Causal=sample.dataset$CausalEst)
## Five Estimate
est=AllEs_C(Y,X,C)
est$est  #Sold        Snew           P        Bold        Bnew 
## bayes new
bay.est.new(Y,X,C)

## bayes old
bay.est.old(Y,X,C)

##para
para.est(Y,X,C)
> para.est(Y,X,C)
$P
[1] 0.3962461

$likelihood
[1] -669.6466
## saturated new
sat.est.new(Y,X,C)

## saturated old
sat.est.old(Y,X,C)

## function of the package
ls("package:BayDR")
?AllEs_C

## dataset in the package
try(data(package="BayDR"))
