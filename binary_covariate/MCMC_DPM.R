#Draw samples from a truncated normal:
rtnorm<-function(n,mu,sigma,lower,upper){ 
  lp<-pnorm(lower,mu,sigma) 
  up<-pnorm(upper,mu,sigma)  
  quntile <- qnorm(runif(n,lp,up),mu,sigma) 
  quntile[quntile==Inf] <- upper[quntile==Inf]
  quntile[quntile==-Inf] <- lower[quntile==-Inf]
  quntile[quntile==Inf] <- 2
  quntile[quntile==-Inf] <- -1
  return (quntile)
}

########################################################
# 
# Inputs
#
#  X     := p-vector of covariates
#  Y     := nt+1-vector of outcomes
#  delta := nt-vector of times between visits
#  A     := nt-vector of recommended times between visits 
# 
# Model
#  
#  (X,Y[1]) ~ N(b0,s0)
#  delta[t] ~ N(X1[,t]%*%b1,s1)
#  Y[t]     ~ N(X2[,t]%*%b2,s2)
#
#  X1[,t]   = (X,A[t],Y[t])
#  X2[,t]   = (X,A[t-1],Y[t-1],delta[t-1]) 
#
#########################################################


MCMC_DPM <- function(base,nt,XX1,XY1,YY1,XX2,XY2,YY2,bi_ind,
                     L=5,DV=.1,eps=1,
                     iters=5000,burn=3000,update=100000){
  
  library(MCMCpack)
  
  n  <- nrow(base)
  p0 <- ncol(base)
  p1 <- ncol(XY1)
  p2 <- ncol(XY2)
  
  # Introduce latent variables for binary covariates
  # The latent variable is in (0,Inf) if the binary covariate is 1
  # The latent variable is in (-Inf, 0) if the binary covariate is 0
  nbi <- length(bi_ind)
  if(nbi>0){
    low_bi <- high_bi <- matrix(0,n,nbi)
    for (i in 1:nbi){
      ind <- bi_ind[i]
      binary <- base[,ind]
      low_bi[,i] <- ifelse(binary==1,0,-Inf)
      high_bi[,i] <- ifelse(binary==1,Inf,0)
    }
  }
  
  b0 <- matrix(0,L,p0)
  b1 <- matrix(0,L,p1)
  b2 <- matrix(0,L,p2)
  
  taub <- rep(1,3)
  
  Q0 <- diag(p0)
  Q1 <- 1 
  Q2 <- 1
  
  m0 <- rep(0,p0)
  m1 <- rep(0,p1)
  m2 <- rep(0,p2)
  
  Qb0 <- diag(p0)
  Qb1 <- diag(p1)
  Qb2 <- diag(p2)
  
  g  <- rep(1:L,n)[1:n]
  g  <- sample(1:L,n,replace=TRUE)
  pg <- rep(1/L,L)
  
  keep.b0 <- array(0,c(iters,L,p0))
  keep.b1 <- array(0,c(iters,L,p1))
  keep.b2 <- array(0,c(iters,L,p2))
  keep.s0 <- array(0,c(iters,p0,p0))
  keep.s1 <- rep(0,iters)
  keep.s2 <- rep(0,iters)
  keep.pg <- matrix(0,iters,L)
  
  for(iter in 1:iters){
    
    # baseline parameters
    
    SSS <- base - b0[g,]
    SSS <- t(SSS)%*%SSS
    Q0  <- rwish(n+p0+eps,solve(SSS+diag(p0)*(p0+eps)))    
    for(l in 1:L){
      id  <- which(g==l)
      ng  <- length(id)
      VVV <- solve(ng*Q0 + Qb0)
      if(ng==0){MMM <- Qb0%*%m0}
      if(ng==1){MMM <- Q0%*%base[id,] + Qb0%*%m0}
      if(ng>1){ MMM <- Q0%*%colSums(base[id,]) + Qb0%*%m0}
      b0[l,] <- VVV%*%MMM + t(chol(VVV))%*%rnorm(p0)   
    }
    
    #Force varaince of binary covariates to be 1
    if(nbi>0){
      Sigma0 <- solve(Q0)
      ifelse(nbi==1,var_bi<-Sigma0[bi_ind,bi_ind],var_bi<-diag(Sigma0[bi_ind,bi_ind]))
      for (i in 1:nbi){
        ind <- bi_ind[i]
        Sigma0[ind,]  <- Sigma0[ind,]/sqrt(var_bi[i])
        Sigma0[,ind]  <- Sigma0[,ind]/sqrt(var_bi[i])
      }
      Q0 <- solve(Sigma0)
      for (i in 1:nbi){
        ind <- bi_ind[i]
        b0[g,ind] <- b0[g,ind]/sqrt(var_bi[i])
      }
      
      # Update baseline latent variables
      for (i in 1:nbi){
        ind <- bi_ind[i]
        muZ <- b0[g,ind] + Sigma0[ind,-ind]%*%solve(Sigma0[-ind,-ind],t(base[,-ind]-b0[g,-ind]))
        varZ <- Sigma0[ind,ind]- Sigma0[ind,-ind]%*%solve(Sigma0[-ind,-ind], Sigma0[-ind,ind])
        base[,ind] <- rtnorm(n,muZ,sqrt(varZ),low_bi[,i], high_bi[,i])
      }
    }
    
    # delta parameters
    
    SS <- sum(YY1)
    for(l in 1:L){
      id  <- which(g==l)
      ng  <- length(id)
      if(ng==0){
        Z0 <- 0
        Z1 <- rep(0,p1)
        Z2 <- 0*diag(p1)
      }
      if(ng==1){
        Z0  <- YY1[id]
        Z1  <- XY1[id,]
        Z2  <- XX1[id,,]
      }
      if(ng>1){ 
        Z0  <- sum(YY1[id])
        Z1  <- apply(XY1[id,],2,sum)
        Z2  <- apply(XX1[id,,],2:3,sum)
      }
      VVV    <- solve(Q1*Z2 + Qb1)
      MMM    <- Q1*Z1 + Qb1%*%m1
      b1[l,] <- VVV%*%MMM + t(chol(VVV))%*%rnorm(p1)
      
      SS     <- SS - 2*sum(b1[l,]*Z1) + t(b1[l,])%*%Z2%*%b1[l,]
    }
    Q1  <- rgamma(1,sum(nt)/2+eps,SS/2+eps)
    
    # Y parameters
    
    SS <- sum(YY2)
    for(l in 1:L){
      id  <- which(g==l)
      ng  <- length(id)
      if(ng==0){
        Z0 <- 0
        Z1 <- rep(0,p2)
        Z2 <- 0*diag(p2)
      }
      if(ng==1){
        Z0  <- YY2[id]
        Z1  <- XY2[id,]
        Z2  <- XX2[id,,]
      }
      if(ng>1){ 
        Z0  <- sum(YY2[id])
        Z1  <- apply(XY2[id,],2,sum)
        Z2  <- apply(XX2[id,,],2:3,sum)
      }
      VVV    <- solve(Q2*Z2 + Qb2)
      MMM    <- Q2*Z1 + Qb2%*%m2
      b2[l,] <- VVV%*%MMM + t(chol(VVV))%*%rnorm(p2)
      
      SS     <- SS - 2*sum(b2[l,]*Z1) + t(b2[l,])%*%Z2%*%b2[l,]
    }
    Q2  <- rgamma(1,sum(nt)/2+eps,SS/2+eps)
    
    # Cluster labels:
    
    bbb1 <- base%*%Q0%*%t(b0)
    bbb2 <- rowSums((b0%*%Q0)*b0)
    for(i in 1:n){
      
      logp <- log(pg) 
      logp <- logp+bbb1[i,]-0.5*bbb2
      logp <- logp+as.vector(XY1[i,]%*%t(b1))*Q1 - 0.5*rowSums((b1%*%XX1[i,,])*b1)*Q1                   
      logp <- logp+as.vector(XY2[i,]%*%t(b2))*Q2 - 0.5*rowSums((b2%*%XX2[i,,])*b2)*Q2                    
      
      g[i] <- sample(1:L,1,prob=exp(logp-max(logp)))
    }
    
    # Hyperparameters
    
    VVV <- solve(L*Qb0 + eps*diag(p0))
    MMM <- Qb0%*%colSums(b0)
    m0  <- VVV%*%MMM + t(chol(VVV))%*%rnorm(p0)
    
    VVV <- solve(L*Qb1 + eps*diag(p1))
    MMM <- Qb1%*%colSums(b1)
    m1  <- VVV%*%MMM + t(chol(VVV))%*%rnorm(p1)
    
    VVV <- solve(L*Qb2 + eps*diag(p2))
    MMM <- Qb2%*%colSums(b2)
    m2  <- VVV%*%MMM + t(chol(VVV))%*%rnorm(p2)
    
    
    SSS <- sweep(b0,2,m0,"-") 
    SSS <- t(SSS)%*%SSS
    Qb0 <- rwish(L+p0+eps,solve(SSS+diag(p0)*(p0+eps))) 
    
    SSS <- sweep(b1,2,m1,"-") 
    SSS <- t(SSS)%*%SSS
    Qb1 <- rwish(L+p1+eps,solve(SSS+diag(p1)*(p1+eps))) 
    
    SSS <- sweep(b2,2,m2,"-") 
    SSS <- t(SSS)%*%SSS
    Qb2 <- rwish(L+p2+eps,solve(SSS+diag(p2)*(p2+eps))) 
    
    
    pg <- rep(1,L)
    for(l in 1:(L-1)){
      pg[l] <- rbeta(1,sum(g==l)+1,sum(g>l)+DV)
    }
    pg[2:L] <- pg[2:L]*cumprod(1-pg[2:L-1])
    
    # Keep track of stuff:
    
    keep.b0[iter,,] <- b0
    keep.b1[iter,,] <- b1
    keep.b2[iter,,] <- b2
    keep.s0[iter,,] <- solve(Q0)
    keep.s1[iter]   <- 1/sqrt(Q1)
    keep.s2[iter]   <- 1/sqrt(Q2)
    keep.pg[iter,]  <- pg
    
    if(iter%%update==0){
      print(paste("Done with",iter,"of",iters))
    }
    
  }
  
  
  
  out <- list(b0=keep.b0[burn:iters,,],
              b1=keep.b1[burn:iters,,],
              b2=keep.b2[burn:iters,,],
              s0=keep.s0[burn:iters,,],
              s1=keep.s1[burn:iters],
              s2=keep.s2[burn:iters],
              pg=keep.pg[burn:iters,]) 
  
  return(out)}

