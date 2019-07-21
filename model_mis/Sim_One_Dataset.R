# Simulate the disease progression within 60 months for one subject
 d    <- gen.one.sub(params_or,policy=randompolicy,alpha=1)
 p0   <- p+1
 p1   <- ncol(d$X1)
 p2   <- ncol(d$X2[,1:8])
 nt   <- rep(0,n)
 base <- matrix(0,n,p0)
 XX1  <- array(0,c(n,p1,p1))
 XY1  <- matrix(0,n,p1)
 YY1  <- matrix(0,n)
 XX2  <- array(0,c(n,p2,p2))
 XY2  <- matrix(0,n,p2)
 YY2  <- matrix(0,n)

 dat    <- NULL
 Ychange <- NULL
 delta <- NULL

 # Compute the sufficient statistics that will be used to fit the model
 for(i0 in 1:n){
    d        <- gen.one.sub(params_or,policy=randompolicy,alpha=1)

    dat[[i0]] <- d
    Z         <- get.sufficient.stats(d$X,d$Y,d$delta,d$X1,d$X2[,1:8])
    nt[i0]    <- Z$nt
    base[i0,] <- Z$base
    XX1[i0,,] <- Z$XX1
    XY1[i0,]  <- Z$XY1
    YY1[i0]   <- Z$YY1
    XX2[i0,,] <- Z$XX2
    XY2[i0,]  <- Z$XY2
    YY2[i0]   <- Z$YY2
    Ychange <- c(Ychange,(d$Y[3:(Z$nt+1)]-d$Y[1:(Z$nt-1)])/d$delta[2:Z$nt])
    delta <- c(delta, d$delta[-Z$nt])
 }

