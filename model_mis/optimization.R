# First batch 

a1   <- ccd_sphere(1,q,5)   # Create a grid of feature weights alpha
junk <- rowSums(a1[,-q]^2)==0
a1   <- a1[!junk,] 
n1   <- nrow(a1)
a1   <- a1 + rnorm(prod(dim(a1)))/50 #Perturb candidat alphas
val1 <- rep(0,n1)
for(k in 1:n1){
  #Estimate values for each alpha in the grid
  if(k==1){print("Estimate values for a grid of alphas:")}
  val1[k] <- get.value(a1[k,],params_fit,m=m,utility=utility,indX_policy=indX_policy)[1]
  if(k%%10==0){
    print(paste("Done with",k,"of",n1))
  }
}
junk   <- is.na(val1)
a1     <- a1[!junk,]
val1   <- val1[!junk]


# Second batch

val <- val1
a   <- a1

# Estimate GP parameters
params_gp <- GP_MLE(val,a) 

for(step in 1:n2){  
  
  # Update current minimum
  
  pred      <- GP_predict(a,a,val,params_gp)$mean  
  thresh    <- min(pred)   
  
  # Pick a new candidate
  
  cana  <- matrix(rnorm(1000*q),1000,q)
  cana  <- sweep(cana,1,sqrt(rowSums(cana^2)),"/")
  pred  <- GP_predict(cana,a,val,params_gp)
  # Select the weight with the largest expected decrease
  g     <- gain(thresh,pred$mean,pred$sd)
  aaa   <- cana[which.max(g),] 
  
  # Refine the candidate weight
  cana  <- matrix(rnorm(1000*q,0,.1),1000,q)
  cana  <- sweep(cana,2,aaa,"+")
  cana  <- sweep(cana,1,sqrt(rowSums(cana^2)),"/")
  pred  <- GP_predict(cana,a,val,params_gp)
  # Select the weight with the largest expected decrease
  g     <- gain(thresh,pred$mean,pred$sd)
  cana  <- cana[which.max(g),] 
  
  # Compute the value of the candidate alpha
  a      <- rbind(a,cana)
  val2   <- get.value(cana,params_fit,m=m,utility=utility,indX_policy=indX_policy)[1]
  val    <- c(val,val2)
  
  if(step==1){print("Sequential optimization:")}
  
  if(step%%10==0){
    print(paste("Done with",step,"of",n2))
  }
}


# Compute the final estimate that minimizes the predictive mean from Gaussian process
cura      <- a[which.min(val),]
for(iter in 1:20){
  cana      <- matrix(rnorm(1000*q,0,1/iter^2),1000,q)
  cana      <- sweep(cana,2,cura,"+")
  cana      <- sweep(cana,1,sqrt(rowSums(cana^2)),"/")
  pred      <- GP_predict(cana,a,val,params_gp)
  cura      <- cana[which.min(pred$mean),]
}
alpha_est <- cura



