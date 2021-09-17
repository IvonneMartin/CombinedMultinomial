## Data analysis including the information regarding the baseline helminth.

rm(list=ls())

setwd("H:/OtherStuff")

library(ecoreg)
library(mvtnorm)
library(devtools)
install_github("IvonneMartin/CombinedMultinomial")


dat.long <- as.data.frame(read.table("datlong.txt",header=TRUE))  ## metadata and bacteria count for 150 subject
phyla <- c(22,16,17,27,33,35)  ## column indices for bacterial phyla of interest, the first one is the reference (Firmicutes)


## Y : multivariate outcome which are the count of 6 bacteria categories.
## X : matrix of covariate values
## Q : number of bacteria categories
## N : number of subject
## Des : design matrix of covariate.
## dat1 : arranging the outcome into wide format.

Y <- as.matrix(dat.long[,phyla])
X <- cbind(dat.long$Treatment,dat.long$inf,dat.long$tp,dat.long$basehelm)+1
Q <- dim(Y)[2]
N <- dim(Y)[1]/2

n.cov <- ncol(X)
lvl.cov <- apply(X,2,function (x) length(levels(as.factor(x))))

FunDes <- DesMat(Q = Q, lvl.cov = lvl.cov) # The left-hand side of the matrix equation and the design matrix, saturated model
# alternatively, user can supply the design matrix and the left hand side of the matrix equation as follows :
  src_dir <- 'H:/OtherStuff/CombMult/DataAnalysis/'
  LHS <- as.matrix(read.table(paste0(src_dir,"LHS64.txt"), header=TRUE))
  Des <- as.matrix(read.table(paste0(src_dir,"DesMat64.txt"), header=TRUE))
  FunDes <- list(LHS, Des)


BactDat <- list("Y" = Y, "X" = X)
## The loglikelihood when we assume that all observations are independent
th <- 0.2
ptm <- proc.time()
pars <- c(rep(0.2,41),log(th))
ansCSFir <- try(optim(pars,loglikcs,dataset = BactDat, Des = FunDes, model = "DMM",method="BFGS",control=list(fnscale=-1,maxit = 5000,parscale=rep(1e-1,length(pars))),hessian=TRUE))
proc.time() - ptm

## takes 145.70s
## Longitudinal estimated using only one random effect:


u1 <- 0.5

ptm <- proc.time()
pars <- c(ansCSFir$par[1:41],log(c(u1)),ansCSFir$par[42])
ansLFir <- try(optim(pars,logliklong, dataset = BactDat, Des = FunDes, model = "DMM",method="BFGS",control=list(fnscale=-1,maxit = 5000,parscale=rep(1e-1,length(pars))),hessian=TRUE))
proc.time() - ptm

## user  system elapsed
## 1097.84    0.76 1098.99

## mixture of chisquare test to test whether the variance of random effect is significantly different than 0

D <- -2*(ansCSFir$value - ansLFir$value)
pchisq(D,df = 1,lower.tail = FALSE)/2

## [1] 0.0003738501

## Presenting the result : Regression terms

RegVar <- strsplit(colnames(FunDes[[2]])[-1],"_")
VarName <- VarNum <- NULL

for (i in 1: length(RegVar)){
  VarName[i] <- RegVar[[i]][1]
  VarNum[i] <- RegVar[[i]][2]
}

## In this loglikelihood, the following variables are estimated: intercept (l_0 and a_); infection (c_), time(d_), treatment at second time point (bd_),
## baseline helminth at the second time points (de_), interaction between treatment at the second time point with baseline helminth (bde_),
## the interaction between treatment at the second time point, baseline helminth and current infection status (bcde_).

Est.vec <- ansLFir$par[c(1:41,43)]
Var.Est <- solve(-ansLFir$hessian[c(1:41,43),c(1:41,43)])
s.u.est <- exp(ansLFir$par[42])
theta.est <- Est.vec[42]


## Intercept form
  idx <- which(VarName == "a")

  Icpt.vec <- Est.vec[which(VarName == "a")]

  SE.vec <- NULL
  for(i in 1:5){
    SE.vec[i] <- sqrt(Var.Est[i,i])
  }
  SE.1 <- SE.vec

## Infection

  idx <- which(VarName == "ac")

  Inf.vec <- Est.vec[idx]
  SE.vec <- NULL
  for(i in 1:5){
    SE.vec[i] <- sqrt(Var.Est[idx[i],idx[i]])
  }
  SE.2 <- SE.vec

## time
  idx <- which(VarName == "ad")

  time.vec <- c(Est.vec[idx])
  SE.vec <- NULL
  for(i in 1:5){
    SE.vec[i] <- sqrt(Var.Est[idx[i],idx[i]])
  }
  SE.3 <- SE.vec

## treatment at the second time point
  idx <- which(VarName == "abd")

  treatt1.vec <- c(Est.vec[idx])
  SE.vec <- NULL
  for(i in 1:5){
    SE.vec[i] <- sqrt(Var.Est[idx[i],idx[i]])
  }
  SE.4 <- SE.vec


## baseline helminth at the second time point
  idx <- which(VarName == "ade")

  bhelmt1.vec <- c(Est.vec[idx])
  SE.vec <- NULL
  for(i in 1:5){
    SE.vec[i] <- sqrt(Var.Est[idx[i],idx[i]])
  }
  SE.5 <- SE.vec

## interaction term between treatment at the second time point and baseline helminth

  idx <- which(VarName == "abde")
  inter1.vec <- c(Est.vec[idx])

  SE.vec <- NULL
  for(i in 1:5){
    SE.vec[i] <- sqrt(Var.Est[idx[i],idx[i]])
  }
  SE.6 <- SE.vec

## interaction term between treatment at the second time point, baseline helminth and current infection status

  idx <- which(VarName == "abcde")
  inter2.vec <- c(Est.vec[idx])

  SE.vec <- NULL
  for(i in 1:5){
    SE.vec[i] <- sqrt(Var.Est[idx[i],idx[i]])
  }
  SE.7 <- SE.vec



## Testing the interaction between baseline helminth, treatment and time. (df = 6)
  idxsel <- which(VarName == "bde" | VarName == "abde")
  Des1 <- FunDes[[2]][,-idxsel]
  FunDes1 <- list(LHS, Des1)

  pars <- ansLFir$par[-c(idxsel)]
  ansLFirinter2 <- try(optim(pars,logliklong, dataset = BactDat, Des = FunDes1, model = "DMM",method="BFGS",control=list(fnscale=-1,maxit = 5000,parscale=rep(1e-1,length(pars))),hessian=TRUE))

  D <- -2*(ansLFirinter2$value - ansLFir$value)
  pval.inter2 <- pchisq(D,df=6,lower.tail = FALSE)


## Testing the interaction between baseline helminth infection with time. (df = 12)
  idxsel <- which(VarName == "de" | VarName == "ade")
  Des1 <- FunDes[[2]][,-idxsel]
  FunDes1 <- list(LHS, Des1)

  pars <- ansLFir$par[-idxsel]
  ansLFirbhelmtime <- try(optim(pars,logliklong, dataset = BactDat, Des = FunDes1, model = "DMM",method="BFGS",control=list(fnscale=-1,maxit = 10000,parscale=rep(1e-1,length(pars))),hessian=TRUE))

  D <- -2*(ansLFirbhelmtime$value - ansLFir$value)
  pval.bhelmtime <- pchisq(D,df=6,lower.tail = FALSE)

## Testing the interaction between treatment and time. (df = 6)
  idxsel <- which(VarName == "bd" | VarName == "abd")
  Des1 <- FunDes[[2]][,-idxsel]
  FunDes1 <- list(LHS, Des1)

  pars <- ansLFir$par[-idxsel]
  ansLFirtrttime <- try(optim(pars,logliklong, dataset = BactDat, Des=FunDes1, model = "DMM",method="BFGS",control=list(fnscale=-1,maxit = 10000,parscale=rep(1e-1,length(pars))),hessian=TRUE))

## Testing the time covariates. (df = 12)

Des <- as.matrix(read.delim("clipboard"))
Des <- Des[,-c(3,17:21)]
pars <- ansLFir$par[c(2,16:20)]
ansLFirtime <- try(optim(pars,Loglik1,method="BFGS",control=list(fnscale=-1,maxit = 10000,parscale=rep(1e-1,length(pars))),hessian=TRUE))

## Testing the infection covariates. (df = 18)

Des <- as.matrix(read.delim("clipboard"))
Des <- Des[,-c(2,12:16)]
pars <- ansLFir$par[-c(1,11:15)]
ansLFirinf <- try(optim(pars,Loglik1,method="BFGS",control=list(fnscale=-1,maxit = 10000,parscale=rep(1e-1,length(pars))),hessian=TRUE))

## Testing the time covariates. (df = 20)

Des <- as.matrix(read.delim("clipboard"))
Des <- Des[,-c(2,4:7,13:17,23:42)]
pars <- ansLFir$par[c(2,7:11,17:21,42:43)]
ansLFirtime <- try(optim(pars,Loglik1,method="BFGS",control=list(fnscale=-1,maxit = 10000,parscale=rep(1e-1,length(pars))),hessian=TRUE))




