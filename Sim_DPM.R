# Those R packages need to be installed before running the demo code
library(fields)
library(MCMCpack)

# Set working directory where you save all the code files
path <- "~/Demo_code"
setwd(path)
source("functions.R")
source("MCMC_DPM.R")


###############################################
#####         SIMULATION DESIGN          ######
###############################################

type <- 2  # 1 for simulating a single group, 2 for simulating mixture groups
source("SimDesign.R")   # Set simulation parameters

m       <- 20000  # Number of MC samples per candiate alpha
n2     <- 100    # Number of steps in the sequential optimization
indX_policy <- c(1)  # Indicator of the elements in X that will be a feature in the policy
q      <- length(indX_policy)+2    # Number of policy features
utility <- "red"   # Choose utility function: "red"/"ave"


####################################################
#####          GENERATE A FAKE DATASET         #####
  
source("Sim_One_Dataset.R")
# This returns a dataset and also a list of sufficient stats X'X, X'Y and Y'Y
# for each subject
  
####################################################
#####             FIT THE MODEL                #####
####################################################
  
fit <- MCMC_DPM(base,nt,XX1,XY1,YY1,XX2,XY2,YY2,update=1000)
  
# Select a subset of the MCMC posterior samples for policy search
iters <- dim(fit$b0)[1]
samps <- sample(1:iters,100,replace=FALSE)
  
B0  <- NULL
B1  <- NULL
B2  <- NULL
PR  <- NULL
for(draw in samps){
  B0 <- rbind(B0,fit$b0[draw,,])
  B1 <- rbind(B1,fit$b1[draw,,])
  B2 <- rbind(B2,fit$b2[draw,,])
  PR <- c(PR,fit$pg[draw,])
}
keep <- PR>0.01
B0   <- B0[keep,]
B1   <- B1[keep,]
B2   <- B2[keep,]
PR   <- PR[keep]
  
params_fit <- list(
  b0      = B0,
  b1      = B1,
  b2      = B2,
  tchols0 = t(chol(apply(fit$s0,2:3,mean))),
  s1      = mean(fit$s1),
  s2      = mean(fit$s2),
  prob    = PR/sum(PR))
  
####################################################
#####             POLICY SEARCH                #####
####################################################
  
# Optimize the feature weights that minimize the opposite number of the value
# (equivalent to maximizing the value).
source("optimization.R") 
  
####################################################
#####            STORE THE OUTPUT              #####
####################################################
  
# Estimate the threshold corresponding to the estimated optimal alpha
thresh  <- get.thresh(alpha_est,params_fit,m=50000,policy=riskscore,indX_policy=indX_policy)
# Compute the value corresponding to the estimated optimal alpha using fitted model
v1      <- get.value(alpha_est,params_fit,m=1000000,thresh=thresh,policy=riskscore,utility=utility,indX_policy=indX_policy)
# Compute the value corresponding to the estimated optimal alpha using oracle model
v2      <- get.value(alpha_est,params_or,m=1000000,thresh=thresh,policy=riskscore,utility=utility,indX_policy=indX_policy)

cat("The estimated optimal feature weights alpha: ",alpha_est)
cat("The threshold: ",thresh)
cat("The value, standard error of the value, cost and standard error of the cost 
    corresponding to the estimated optimal alpha using fitted model:\n",c(-v1[1],v1[2:4]))
cat("The value, standard error of the value, cost and standard error of the cost
    corresponding to the estimated optimal alpha using oracle model:\n",c(-v2[1],v2[2:4]))
