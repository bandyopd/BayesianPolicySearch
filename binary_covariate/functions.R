#Simulate the disease progression within 60 months for one subject
gen.one.sub <- function(params,bi_ind,policy=riskscore,alpha=NULL){
  
  #Sample group assignment for the subject
  g    <- sample(1:length(params$prob),1,prob=params$prob) 
  p    <- nrow(params$tchols0)-1
  #Simulate the baseline information under assigned group
  XY   <- params$b0[g,] + params$tchols0%*%rnorm(p+1)
  X    <- XY[1:p]
  for(j in bi_ind){
    X[j] <- ifelse(X[j]>0,1,0)
  }
  curY <- XY[p+1]
  Y    <- curY
  d    <- NULL
  A    <- NULL
  t    <- 0
  curd <- 6
  curA <- 6
  lagY <- curY
  
  X1 <- NULL
  X2 <- NULL
  
  #Simulate the disease progression for the subject within 60 months
  while(t<60 | length(Y)<3){
    lagA  <- curA
    lagd  <- curd
    lag2Y <- lagY
    lagY  <- curY

    curA  <- policy(X,log(abs(lagd-lagA)+1),lagY,lagY-lag2Y,alpha)
    
    A     <- c(A,curA)
    
    Xd    <- c(1,X,lagY,log(curA),X*log(curA),lagY*log(curA))
    curd  <- exp(sum(Xd*params$b1[g,])+params$s1*rnorm(1))+.1
    curd  <- ifelse(is.na(curd),60,curd)  
    curd  <- ifelse(curd>60,60,curd)  
    d     <- c(d,curd)
    
    Xy    <- c(1,X,lagY,(curd-6),X*(curd-6),(curd-6)*lagY)
    curY  <- sum(Xy*params$b2[g,])+params$s2*rnorm(1)
    Y     <- c(Y,curY)
    t     <- t+curd
    
    X1    <- rbind(X1,Xd)
    X2    <- rbind(X2,Xy)
    
  }       
  
  out <- list(X=X,Y=Y,A=A,delta=d,X1=X1,X2=X2)
  
  return(out)}


# Simulate the disease progression within 60 months for m subjects
gen.m.subs <- function(m,params,bi_ind,alpha=NULL,policy=riskscore,utility="red",indX_policy){
  
  p    <- nrow(params$tchols0)-1
  #Sample group assignment for m subjects
  g    <- sample(1:length(params$prob),m,prob=params$prob,replace=TRUE)
  b0   <- params$b0[g,]
  b1   <- params$b1[g,]
  b2   <- params$b2[g,]
  
  #Simulate the baseline information for m subjects
  XY   <- matrix(rnorm(m*(p+1)),p+1,m)
  XY   <- b0 + t(params$tchols0%*%XY)
  X    <- XY[,1:p]
  for(j in bi_ind){
    X[,j] <- ifelse(X[,j]>0,1,0)
  }
  curY <- XY[,p+1]
  Y    <- curY
  d    <- NULL
  A    <- NULL
  t    <- 0
  curd <- 6
  curA <- 6
  lagY <- curY
  
  b1   <- params$b1[g,]
  b2   <- params$b2[g,]
  
  #Simulate the disease progression for m subjects within 60 months
  visit <- 1
  while(min(t)<60 | visit<3){
    lagA  <- curA
    lagd  <- curd
    lag2Y <- lagY
    lagY  <- curY
    
    curA  <- policy(X[,indX_policy],log(abs(lagd-lagA)+1),lagY,lagY-lag2Y,alpha)
    
    A     <- cbind(A,curA)
    
    int   <- sweep(X,1,log(curA),"*")
    Xd    <- cbind(1,X,lagY,log(curA),int,lagY*log(curA))
    curd  <- exp(rowSums(Xd*b1)+params$s1*rnorm(m))+.1
    curd  <- ifelse(is.na(curd),60,curd)  
    curd  <- ifelse(curd>60,60,curd)  
    d     <- cbind(d,curd)
    
    int1  <- sweep(X,1,(curd-6),"*")
    #int2  <- sweep(X,1,lagY,"*")
    Xy    <- cbind(1,X,lagY,(curd-6),int1,(curd-6)*lagY)
    curY  <- rowSums(Xy*b2)+params$s2*rnorm(m)
    Y     <- cbind(Y,curY)
    
    t     <- t+curd
    visit <- visit+1
  }       
  
  if(utility=="ave"){
    # Average utility function
    DY       <- Y[,-1]
    DY       <- ifelse(DY>0,DY-0,0)
    tt       <- t(apply(d,1,cumsum))<=60
    junk     <- is.na(DY+tt+A)
    DY[junk] <-0
    tt[junk] <-0
    A[junk]  <-0
    
    # Mean responses within 60 months as the value
    r    <- rowSums(DY*tt)/rowSums(tt) 
    # Average recommended recall time
    c    <- rowSums(A*tt)/rowSums(tt)  
    
  }else if(utility=="red"){
    # Reduction utility function
    DY       <- Y[,-1]
    DY       <- ifelse(DY>0,DY-0,0)
    tt       <- t(apply(d,1,cumsum))<=60
    cumd <- t(apply(d,1,cumsum))
    tt <- ifelse(tt==TRUE,1,0)
    # The number of visits within 60 months + 1
    last <- apply(tt,1,which.min) 
    # The response of the last visit within 60 months
    Y.lastbf <- DY[cbind(1:m,last-1)] 
    # The response of the first visit after 60 months
    Y.lastaf <- DY[cbind(1:m,last)] 
    # The time between first visit and the last visit within 60 months
    d.lastbf <- cumd[cbind(1:m,last-1)]
    # The time between first visit and the first visit after 60 months
    d.lastaf <- cumd[cbind(1:m,last)]
    # Using interpolation to estimate the response on the 60th month
    # and negative reduction as the value
    r <- (Y.lastaf-Y.lastbf)/(d.lastaf-d.lastbf)*(60-d.lastbf)+Y.lastbf-Y[,1]
    
    junk     <- is.na(DY+tt+A)
    DY[junk] <-0
    tt[junk] <-0
    A[junk]  <-0
    
    # The average recommended recall time
    c    <- rowSums(A*tt)/rowSums(tt) 
  }else{stop('The utility function is not defined')}
  
  out <- list(utility=r,cost=c)
  
  return(out)}




# Compute the sufficient statistics X'X, X'Y, Y'Y
get.sufficient.stats <- function(X,Y,delta,X1,X2){
  
  # Baseline data
  base <- c(X,Y[1])
  
  # Time between visits
  logY <- log(ifelse(delta<=0.1,0.001,delta-0.1))
  
  YY1 <- sum(logY^2)
  XX1 <- t(X1)%*%X1
  XY1 <- as.vector(t(X1)%*%logY)
  
  # Disease progression
  
  YY2 <- sum(Y[-1]^2)
  XX2 <- t(X2)%*%X2
  XY2 <- as.vector(t(X2)%*%Y[-1])
  
  # Output
  nt    <- length(delta)   
  stats <- list(nt=nt,base=base,
                XX1=XX1,XY1=XY1,YY1=YY1,
                XX2=XX2,XY2=XY2,YY2=YY2)
  
  return(stats)}

# Policy that recommends recall interval A[t] based on logit(Pr(A[t]=3))=Y[t-1]
randompolicy <- function(X,DA,Y,diffY,alpha=1){
  prob <- 1/(1+exp(-Y))
  A    <- rbinom(length(prob),1,prob)
  A    <- ifelse(A==1,3,9)     
  return(A)}

# Baseline policy that recommends A=6 months between visits for all subjects and t
policy6 <- function(X,DA,Y,diffY,alpha){
  6+0*Y
}

# Policy based on risk score
riskscore<- function(X,DA,Y,diffY,alpha){
  risk <- cbind(X,DA,Y)%*%alpha[-1]
  A    <- ifelse(risk>alpha[1],3,9)     
  return(A)}


# Compute the value of a policy defined by a

get.value <- function(a,params,bi_ind,m=5000,policy=riskscore,thresh=NULL,utility="red",indX_policy){
  
  # Estimate the threshold that gives average recommended recall time = 6
  if(is.null(thresh)){
    thresh <- get.thresh(a,params,bi_ind,policy=policy,indX_policy=indX_policy)
  }
  alpha   <- c(thresh,a)
  dat     <- gen.m.subs(m,params,bi_ind,alpha=alpha,policy=policy,utility=utility,indX_policy=indX_policy)  
  results <- c(mean(dat$utility),se(dat$utility),
               mean(dat$cost),se(dat$cost))
  names(results)<-c("Value","SE(Value)","Cost","SE(Cost)")
  
  return(results)}



get.thresh <- function(a,params,bi_ind,m=2000,policy=riskscore,npts=8,indX_policy){
  
  a  <- a/sqrt(sum(a^2))
  p  <- ncol(params$tchols0)-1
  
  # Set the range of threshold to search over
  
  g    <- sample(1:length(params$prob),m,prob=params$prob,replace=TRUE)
  base <- params$b0[g,]+t(params$tchols0%*%matrix(rnorm(m*(p+1)),p+1,m))
  risk <- as.vector(base[,c(indX_policy,p+1)]%*%a[c(1:length(indX_policy),3)])
  
  # First pass - broad search over entire range
  
  cant <- mean(risk) + seq(-5,5,length=npts)*sd(risk)
  ntau <- length(cant)
  ccc  <- NULL
  
  for(ttt in 1:ntau){
    dat      <- gen.m.subs(m=m,params=params,bi_ind,
                           alpha=c(cant[ttt],a),policy=policy,indX_policy=indX_policy)  
    ccc[ttt] <- mean(dat$cost)
  }
  if(ntau>6){ccc <- loess(ccc~cant)$fitted}
  
  # Second pass - targeted search on a narrow window
  cant <- c(cant[1]-5*sd(risk),
            cant,
            cant[ntau]+5*sd(risk))
  ccc  <- c(-Inf,ccc,Inf)
  id   <- sum(ccc<6)
  cant <- seq(cant[id],cant[id+1],length=npts)
  ntau <- length(cant)
  ccc  <- NULL
  for(ttt in 1:ntau){
    dat      <- gen.m.subs(m,params,bi_ind,
                           alpha=c(cant[ttt],a),policy=policy,indX_policy=indX_policy)  
    ccc[ttt] <- mean(dat$cost)
  }
  if(ntau>6){ccc <- loess(ccc~cant)$fitted}
  
  # Select the threshold by interpolation
  id  <- sum(ccc<6)
  if(id==0){thresh<-min(cant)}
  if(id==ntau){thresh<-max(cant)}
  if(id>0 & id<ntau){
    slope  <- (ccc[id+1] - ccc[id])/(cant[id+1]-cant[id])
    thresh <- (6-ccc[id])/slope + cant[id]
  }  
  
  return(thresh)}


se <- function(x){sd(x)/sqrt(length(x))}

# Central composite design on a sphere:
ccd_sphere <- function(width,p,m){
  grid <- seq(-width,width,length=m)
  d    <- NULL
  for(j in 1:p){
    d[[j]] <- grid
  }
  x    <- as.matrix(expand.grid(d))
  s    <- sqrt(rowSums(x^2))
  x    <- x[s>0,]
  s    <- s[s>0]
  x    <- sweep(x,1,s,"/") # constrain L2 norm of x to be 1
  x    <- unique(x)
  return(x)}



############################################################
########  Gaussian process interpolation functions  ########
############################################################

# Maximum likelihood estimates of parameters in the Gaussian process
GP_MLE <- function(Y,X){
  mu    <- mean(Y)
  sigma <- sd(Y)
  Z     <- (Y-mu)/sigma
  
  # Likelihood of the Guassian process
  GP_nll <- function(params,Y,X,r){
    phi     <- exp(params)
    X       <- sweep(X,2,phi,"/")
    d       <- as.matrix(dist(X))   
    S       <- r*exp(-d^2)
    diag(S) <- 1   
    C       <- t(chol(S))
    nll     <- sum(log(diag(C))) + 0.5*sum(solve(C,Y)^2)
    return(nll)}
  
  r <- 0.99
  # MLE of parameters in the Gaussian process
  logphi    <- optim(rep(-2,ncol(X)),GP_nll,Y=Z,X=X,r=0.99)$par
  
  return(list(mu=mu,sigma=sigma,logphi=logphi,r=r))}

make.Q  <- function(X,Y,params){
  Z       <- (Y-params$mu)/params$sigma
  phi     <- exp(params$logphi)
  X       <- sweep(X,2,phi,"/")
  d       <- as.matrix(dist(X))   
  S       <- params$r*exp(-d^2)
  diag(S) <- 1   
  Q       <- solve(S)
  return(Q)}

# Predict the values at new data points
GP_predict <- function(Xnew,X,Y,params){
  
  library(fields)
  
  Z      <- (Y-params$mu)/params$sigma
  S11inv <- make.Q(X,Y,params)
  phi    <- exp(params$logphi)
  
  x1     <- sweep(Xnew,2,phi,"/")
  x2     <- sweep(X,2,phi,"/")
  d      <- rdist(x1,x2)
  S12    <- params$r*exp(-d^2)
  S12S11 <- S12%*%S11inv
  EY     <- S12S11%*%Z
  VY     <- 1-rowSums(S12S11*S12)
  # Use kriging to predict values
  EY     <- params$mu + params$sigma*EY
  VY     <- params$sigma*params$sigma*VY
  
  return(list(mean=EY,sd=sqrt(VY)))}

# Compute the expected decrease
gain <- function(thresh,mu,sigma){
  z     <- (thresh-mu)/sigma
  cdf   <- pnorm(z)
  pdf   <- dnorm(z)
  Egain <- cdf*(thresh-mu)+sigma*pdf
  return(Egain)}


