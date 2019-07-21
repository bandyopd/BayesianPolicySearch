 # General two group model
  p         <- 2 #Number of baseline covariates

  # Set true values of coefficients
  b0_or     <- rep(0,p+1) # E(X,Y_base) 
  b1_or     <- rep(0,2*p+4)    # log(delta) | int, X, lagY, log(A) X*log(A) lagY*log(A)
  b2_or     <- rep(0,2*p+6)    #     Y | int, X, lagY, D-6, X*(D-6), lagY*(D-6), D^2,lagY^2
  b1_or[5]  <- 0.9  # Effect of A on D
  b1_or[6]  <- 0.1  # Effect of X[1]*A on D
  b2_or[1]  <- 0.1  # Intercept
  b2_or[3]  <- 0.2  # Effect of X[2] on Y
  b2_or[5]  <- -0.3  # Effect of D on Y
  b2_or[4]  <- 0.9  # Effect of Y[t-1] on Y
  b2_or[8] <- 0.02 # Effect of D*Y[t-1] on Y

  b0_or    <- rbind(b0_or,b0_or)
  b1_or    <- rbind(b1_or,b1_or)
  b2_or    <- rbind(b2_or,b2_or) 

  b0_or[2,1]  <- 1
  b1_or[2,]   <- 0
  b1_or[2,1]  <- log(5.3)
  b2_or[2,2]  <- 0.3
  b2_or[2,3]  <- 0
  b2_or[2,5]  <- -0.2
  b2_or[2,8] <- 0
  b2_or[1,1] <- -2
  b2_or[1,9] <- 0.06

  s0_or  <- 0.5+0.5*diag(p+1)  # Covariance matrix of (X,Y_base)
  tc_or  <- t(chol(s0_or))
  # Standard deviation of the error term in the compliance model
  # and progression model
  s1_or  <- .1   
  s2_or  <- .5


 # Define the four models

 prob_or <- c(1,0)   # Group assignment probability (single group)
 if(type==2){
   prob_or <- c(0.8,0.2)  # Group assignment probability (mixture groups)
 }

 n <- 1000   # Sample size

 # True value of model parameters
 params_or <- list(b0=b0_or,b1=b1_or,b2=b2_or,
                   tchols0=tc_or,s1=s1_or,s2=s2_or,
                   prob=prob_or)

