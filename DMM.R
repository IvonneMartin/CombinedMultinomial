rm(list=ls())

########################################################################################################
##	                Dirichlet - multinomial mixed regression model                                ##
##   						2018 - 08 - 20			                      ##
########################################################################################################

# libraries --------------------------------------------------------------------------------------------
library(HMP)
library(ecoreg)

# log-likelihood function ------------------------------------------------------------------------------
Loglik <- function(params){
 # N : number of samples
 # T : number of repeated measurements, defined by user
 # Q : number of bacterial categories 
 # P : the number of covariates per sample
 # Y : TN times Q matrix of bacterial counts 
 # X : TN times P matrix of covariates 
 # Des : design matrix for the covariates effect 

	theta <- exp(params[length(params)])
	s.u <- exp(params[(length(params) - 1)])
	lambdas <- c(0,params[1:(length(params) - 2)])

	l <- dim(Y)[1]
	N <- l/T
	Q <- dim(Y)[2]
	L <- NULL 
	G <- matrix(NA,nrow = l, ncol = Q)

	f <- match.fun(F)
	G <- t(apply(X,1,FixEf,lambdas))
	index <- seq(1,N,by = 1)
	L <- sapply(index, evalAGHQ, G, overdisp = theta, sigma.u = s.u,Y,f = f,l = l,N = N, T = T, Q= Q)	
	
res <- sum(log(L))	
return(res)}

# The marginal distribution for the multivariate count ------------ -------------------------------------

CNB <- function(Ct,g){
	gs <- sum(g)
	ys <- sum(Ct)
	val <- lgamma(ys+1) + lgamma(gs) - lgamma(ys+gs) + 
		 sum(lgamma(Ct+g) - lgamma(g) - lgamma(Ct+1)) 
return(val)
}

# The function to compute the distribution conditioned on the random effect----

F <- function(z,sigma.u,mus,overdisp,Yt,T,Q)
{	Eta <- (1/overdisp)*(exp(mus + rep(z,Q)))
	val <- NULL
      for (j in 1:T){
		val[j] <- CNB(Yt[j,],Eta[j,])
		}
 est <- sum(val) + dnorm(z,mean = 0,sd = sigma.u,log=TRUE)
return(est)
}


# Comnputing the Adaptive Gauss-Hermite Quadrature approximation ----------------------------------------------------------
evalAGHQ <- function(i,G,overdisp,sigma.u,Y,f,l,N,T,Q){
 rpt <- seq(i,l,by = N)
	
 Yt <- matrix(NA,nrow = T, ncol = Q)
 mus <- matrix(NA,nrow = T,ncol=Q)
	
	Yt <- Y[rpt,]
	mus <- G[rpt,]

	FF1 <- function(z) F(z,sigma.u,mus,overdisp,Yt = Yt,T,Q)
	opt <- try(optim(0.1,FF1,method="BFGS",control=list(fnscale=-1,maxit=5000),hessian=TRUE))
	est <- integrate.gh(function(z) exp(FF1(z)),mu = opt$par,scale=1/sqrt(-opt$hessian),points=20)
 
 return(est)
}

# Simulating dataset following the Dirichlet - multinomial mixed model 

FixEf <- function(x.vec,b){
 gt <- c(Des %*% b)
		if (x.vec[1] == 0) {
  	  		g = gt[seq(1,6,by=2)]
		} else
  			g = gt[seq(2,6,by=2)]
return(g)}

simulDMM <- function(N,Q,S,theta,s.u,T = 2,P = 1){
	U <- rnorm(N,mean=0,sd=s.u)
	X <- rbinom(N,size=1,prob = 0.5)
 	Des <- as.matrix(cbind(rep(1,6),rep(0:1,3),c(0,0,1,1,0,0),c(rep(0,4),rep(1,2)),c(0,0,0,1,0,0),c(rep(0,5),1)))

	XB1 <- t(sapply(X[,1],FixEf,b = beta))
	XB.tilde1 <- (1/th)*exp(cbind(XB1+rep(U,3)))
	C1 <- t(apply(XB.tilde1,1,Dirichlet.multinomial,Nrs = S))

	XB2 <- t(sapply(X[,2],FixEf,b = beta))
	XB.tilde2 <- (1/th)*exp(cbind(XB2+rep(U,3)))
	C2 <- t(apply(XB.tilde2,1,Dirichlet.multinomial,Nrs = S))

return(list(C1 = C1, C2 = C2, X = X))
}
