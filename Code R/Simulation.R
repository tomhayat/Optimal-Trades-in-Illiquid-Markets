# Set directories
setwd("/Users/danconstantini/Downloads/Projet\ Mopsi/Data ")
# Parameters
lambda <- 4
mu <- 8
Time <- 1
Q <- t(matrix(data =c(-2,2,0,1,-4,3,0,2,-2),nrow = 3, ncol = 3))
Q <- t(Q)
lambda.Markov <- 4*c(3, 3, 1)
mu.Markov <- c(8,4,4)
Gamma <- diag(lambda.Markov)

# random generation
gener <- function(x){
  return (-log(runif(1))/x);
}
PoissonState <- function(lambda,Time){
  X <- gener(lambda)
  if( X < Time){
    return(X)
  }
  else
    return(Time)
}
proba <- function(Q,i){
  state <- seq(1:3)
  state <- state[-i]
  p <- Q[i,c(-i)]
  p <- p/sum(p)
  X <- runif(1)
  if ( X < p[1] )
    return (state[1])
  else
    return(state[2])
}

# Poisson Process
PoissonProcess <- function(Time, lambda ){
  sigma <- c(PoissonState(lambda,Time))
  while(sigma[length(sigma)] < Time){
    sigma <- c(sigma, sigma[length(sigma)] + PoissonState(lambda,Time)) 
  }
  if (length(sigma >1)){
    return (sigma[1:length(sigma)-1])
  }
  else{
    return (sigma)
  }
}
test.PoissonProcess <- function(){
  OrderBook <- PoissonProcess(1, lambda)
  print(OrderBook)
  stripchart(OrderBook)
}
Example.PoissonProcess <- function(k,a){
  #a <- result.read("a.txt")
  OrderBook <- PoissonProcess(1, lambda)
  stripchart(OrderBook)
  Opti <- OptimalStrategy.PoissonProcess(a,OrderBook, k)
  print(Opti)
  plot(Opti[2:length(Opti[,2]),2]~OrderBook, type="h", xlim=c(0,1), lwd = 3,
       xlab = "Time", ylab = "Sold shares", main = "Optimal strategy for Poisson Process simulation")
  print(-diff(Naive.Vs.Opti(Opti)))
  print(Mean.PoissonProcess(a,k))
}
# Compound Poisson Process
CompoundPoissonProcess <- function( Time, mu, lambda ){
  sigma <- c(PoissonState(lambda,Time))
  ordersize <- c(floor(rpois(1,mu)*(1-exp(-mu))))
  while(sigma[length(sigma)] < Time){
    sigma <- c(sigma, sigma[length(sigma)] + PoissonState(lambda,Time)) 
    ordersize <- c(ordersize, floor(rpois(1,mu)*(1-exp(-mu))))
  }
  if (length(sigma) >1){
    #return (ordersize)
    return(cbind(sigma[1:length(sigma)-1],ordersize[1:length(ordersize)-1]))
  }
  else{
    #return (ordersize)
    return(cbind(as.numeric(sigma),as.numeric(ordersize)))
  }
}
test.CompoundPoissonProcess <- function(){
  OrderBook <- CompoundPoissonProcess(1, mu, lambda)
  print(OrderBook)
  plot(OrderBook[,2]~OrderBook[,1], type ="h", 
        xlab = "Time to Maturity",
        ylab = "Quantity Placed",
        main = "Simulation of an Order Book",
        sub = "Compound Poisson Process"
        )
}
Example.CompoundPoissonProcess <- function(k,apois){
  #apois <- result.read("apois.txt")
  OrderBook <- CompoundPoissonProcess(1,mu, lambda)
  Opti <- OptimalStrategy.CompoundPoissonProcess(apois,OrderBook, k)
  #print(Opti)
  plot(OrderBook[,2]~OrderBook[,1], col = "red", xlim=c(0,1),
       xlab = "Time", ylab = "Sold shares", main = "Optimal strategy for Compound Poisson Process simulation")
  points(OrderBook[,1],Opti[2:length(Opti[,2]),2], type="h",lwd = 3)
  
  #print(-diff(Naive.Vs.Opti(Opti)))
  #print(Mean.CompoundPoissonProcess(apois,k))
}

# Markov Switch Regime
PoissonState.Markov <- function(Q,i,Time){
  X <- gener(abs(Q[i,i]))
  if( X < Time){
    return(X)
  }
  else
    return(Time)
}
GenerateStates <- function(Q,Time){
  regime <- c(1)
  OrderBook <- c(PoissonState.Markov(Q,regime[length(regime)],Time))
  while(OrderBook[length(OrderBook)] < Time ){
    regime <- c(regime,proba(Q,regime[length(regime)]))
    OrderBook <- c(OrderBook,OrderBook[length(OrderBook)] +
                     PoissonState.Markov(Q,regime[length(regime)],Time))
  }
  OrderBook[length(OrderBook)] <- Time
  return(cbind(OrderBook,regime))
}
GenerateTime <- function(States, Time, Gamma){
  regime <- 1
  R <- c(regime)
  sigma <- c(PoissonState.Markov(Gamma, regime,Time))
  ordersize <- c(rpois(1,mu.Markov[regime])/(1-exp(-mu.Markov[regime])))
  while(sigma[length(sigma)] < Time){
    regime <- States[which(States[,1]>sigma[length(sigma)])[1],2]
    R <- c(R, regime)
    sigma <- c(sigma,sigma[length(sigma)] + PoissonState.Markov(Gamma, regime,Time)) 
    ordersize <- c(ordersize, rpois(1,mu.Markov[regime])/(1-exp(-mu.Markov[regime])))
  }
  if (length(sigma) > 1){
    return(cbind(sigma[1:length(sigma)-1],ordersize[1:length(ordersize)-1],R[1:length(R)-1]))
  }
  else{
    return(cbind(sigma,ordersize,R))
  }
}
test.RegimeSwitching <- function(){
  Order <- GenerateStates(Q,1)
  print(Order)
  Sigma <- GenerateTime(Order, Time, Gamma)
  plot( Sigma[,2]~Sigma[,1], 
        xlim = c(0,Time), 
        ylim = c(0,max(Sigma[,2])+3), 
        type ="h",
        xlab = "Time to Maturity",
        ylab = "Quantity Placed",
        main = "Simulation of an Order Book",
        sub = "Regime Switching"
        
        )
  for( i in (1:length(Order[,1]))){
    abline(v = Order[i,1], col = 30*Order[i,2])
  }
  print(Sigma)
}
Example.RSS <- function(k, as){
#   as1 <- result.read("as1.txt")
#   as2 <- result.read("as2.txt")
#   as3 <- result.read("as3.txt")
#   as <- list(as1,as2,as3)
  Order <- GenerateStates(Q,1)
  Sigma <- GenerateTime(Order, Time, Gamma)
  Opti <- OptimalStrategy.RS(as, Sigma, k)
  plot( Sigma[,2]~Sigma[,1], 
        xlim = c(0,Time), 
        ylim = c(0,max(Sigma[,2])+3), 
        xlab = "Time to Maturity",
        ylab = "Quantity Placed",
        main = "Simulation of an Order Book",
        sub = "Regime Switching"
        
  )
  for( i in (1:length(Order[,1]))){
    abline(v = Order[i,1], col = 30*Order[i,2])
  }

  points(Sigma[,1],Opti[2:length(Opti[,2]),2], type="h",lwd = 3)
  #print(Mean.CompoundPoissonProcess(apois,k))
}

# Read Data
result.read <- function(filename){
  data <- read.csv(filename,check.names = FALSE,header = FALSE)
  data <- matrix(unlist(data), nrow = dim(data)[1],ncol = dim(data)[2])
  return(data)
}

# Poisson Process Optimal Strategy
a <- result.read("a.txt")
Segmentation <- function(x, p){
  return(floor(x*p))
}
OptimalStrategy.PoissonProcess <- function(a, OrderBook, amount){
  OB <- Segmentation(OrderBook,dim(a)[2])
  OB[which(OB == dim(a)[2])] <- dim(a)[2] -1
  x <- c(amount)
  k <- c(0)
  for (i in 1:length(OB)){
    k <- c(k,a[x[length(x)]+1,OB[i]+1])
    x<-c(x, x[length(x)] - k[length(k)])
  }
  return(cbind(x,k))
}
Cost.Function <- function(x){
  return (as.numeric(abs(x)^2/2))
}
Naive.Vs.Opti <- function(Optimal){
  N <- as.numeric(Optimal[1,1])
  return(c(Cost.Function(N), Cost.Function(N - sum(Optimal[,2]))+ sum(Cost.Function(Optimal[,2]))))
}
Mean.PoissonProcess <- function(a, k){
  X <- c()
  for (i in (1:1000)){
    OrderBook <- PoissonProcess(1, lambda)
    Opti <- OptimalStrategy.PoissonProcess(a, OrderBook, k)
    X <- c(X,Naive.Vs.Opti(Opti)[1]-Naive.Vs.Opti(Opti)[2])
  }
  X <- X[which(!is.na(X))]
  X <- X/Cost.Function(k)
  vari <- sqrt(var(X))
  beta <- 1.96
  return(c(sum(X)/length(X),beta*vari/sqrt(1000)))
}
Overview.PoissonProcess <- function(){
  a <- result.read("a.txt")
  Y <- c()
  for(i in (1:18)){
    OrderBook <- PoissonProcess(1, lambda)
    Y <- rbind(Y, Mean.PoissonProcess(a,i))
  }
  plot(Y[,1]+Y[,2], cex=0.01, type ="l")
  points(1:length(Y[,1]),Y[,1], cex= 0.01, type ="l",  col ="red")
  points(1:length(Y[,1]),Y[,1]-Y[,2], cex= 0.01, type ="l")
}

# Compound Poisson Process Optimal Strategy
apois <- result.read("apois.txt")
OptimalStrategy.CompoundPoissonProcess <- function(apois, OrderBook, amount){
  OB <- OrderBook
  OB[,1] <- Segmentation(OB[,1],dim(apois)[2])
  OB[which(OB == dim(apois)[2]),1] <- dim(apois)[2]-1
  x <- c(amount)
  k <- c(0)
  for (i in 1:length(OB[,1])){
    w <- min(apois[x[length(x)]+1,OB[i,1]+1], OB[i,2])
    k <- c(k,w)
    x<-c(x, x[length(x)] - k[length(k)])
  }
  return(cbind(x,k))
}
Mean.CompoundPoissonProcess <- function(a, k){
  X <- c()
  for (i in (1:1000)){
    OrderBook <- CompoundPoissonProcess(1, mu, lambda)
    Opti <- OptimalStrategy.CompoundPoissonProcess(apois, OrderBook, k)
    X <- c(X,Naive.Vs.Opti(Opti)[1]-Naive.Vs.Opti(Opti)[2])
  }
  X <- X[which(!is.na(X))]
  X <- X/Cost.Function(k)
  vari <- sqrt(var(X))
  beta <- 1.96  # quantile dordre 95%
  return(c(sum(X)/length(X),beta*vari/sqrt(1000)))
}
Overview.CompoundPoissonProcess <- function(){
  a <- result.read("apois.txt")
  Y <- c()
  for(i in (2:18)){
    Y <- rbind(Y, Mean.CompoundPoissonProcess(a,i))
  }
  plot(Y[,1]+Y[,2], cex=0.01, type ="l")
  points(1:length(Y[,1]),Y[,1], cex= 0.01, type ="l",  col ="red")
  points(1:length(Y[,1]),Y[,1]-Y[,2], cex= 0.01, type ="l")
}

# Regime Switching
as1 <- result.read("as1.txt")
as2 <- result.read("as2.txt")
as3 <- result.read("as3.txt")
as <- list(as1,as2,as3)
OptimalStrategy.RS <- function(as, OrderBook, amount){
  OB <- matrix(as.numeric(OrderBook), ncol = 3,)
  OB[,1] <- Segmentation(OB[,1],dim(as[[1]])[2])
  OB[which(OB == dim(as[[1]])[2]),1] <- dim(as[[1]])[2]-1
  x <- c(amount)
  k <- c(0)
  for (i in 1:length(OB[,1])){
    w <- min(as[[OB[i,3]]][x[length(x)]+1,OB[i,1]+1], OB[i,2])
    k <- c(k,w)
    x<-c(x, x[length(x)] - k[length(k)])
  }
  return(cbind(x,k))
}
Mean.RS <- function(as, k){
  X <- c()
  for (i in (1:1000)){
    Order <- GenerateStates(Q,1)
    OrderBook <- GenerateTime(Order, Time, Gamma)
    Opti <- OptimalStrategy.RS(as, OrderBook, k)
    X <- c(X,Naive.Vs.Opti(Opti)[1]-Naive.Vs.Opti(Opti)[2]) 
  }
  X <- X[which(!is.na(X))]
  X <- X/Cost.Function(k)
  vari <- sqrt(var(X))
  beta <- 1.96
  return(c((sum(X)/length(X)),beta*vari/sqrt(1000)))
}
Overview.RS <- function(){
  Y <- c()
  for(i in (2:18)){
    Y <- rbind(Y, Mean.RS(as,i))
  }
  plot(Y[,1]+Y[,2], cex=0.01, type ="l")
  points(1:length(Y[,1]),Y[,1], cex= 0.01, type ="l",  col ="red")
  points(1:length(Y[,1]),Y[,1]-Y[,2], cex= 0.01, type ="l")
}

comparison <- function(){
  a <- result.read("a.txt")
  apois <- result.read("apois.txt")
  ac <- as.numeric(unlist(read.table("acontinuous.txt")))
  v <- result.read("v.txt")
  vpois <- result.read("vpois.txt")
  
  par(mfrow=c(1,2))
  plot(apois[20,]+1~seq(0,0.995, by = 0.005),cex= 0.1,col = "red",     main = "Model Comparison for k = 19",
       xlab =  "Time",
       ylab = "Optimal Cost Value")
  points(seq(0,0.995, by = 0.005),a[20,]+1.0125,cex= 0.1, col ="blue")
  points(seq(0,0.995, by = 0.005),apois[20,]+1,cex= 0.1, col ="red")
  #points(seq(0,1, by = 0.005),20*ac,cex= 0.1, col ="green")
  
  legend("topright",  
         legend= c("Base Model", "Constrained Model"),cex=0.4,
         lty=c(1,1),
         lwd=c(2.5,2.5),
         col=c("blue","red")
  )
  
  
  plot(vpois[20,]~seq(0,0.995, by = 0.005),cex= 0.1,col = "red",
       main = "Model Comparison for k = 19",
       xlab =  "Time",
       ylab = "Optimal Execution Cost")
  points(seq(0,0.995, by = 0.005),v[20,],cex= 0.1, col ="blue")
  points(seq(0,0.995, by = 0.005),vpois[20,],cex= 0.1, col ="red")
  
  #points(seq(0,1, by = 0.005),20^2*(ac),cex= 0.1, col ="green")
  legend("topright",  
         legend= c("Base Model", "Constrained Model"),cex=0.4,
         lty=c(1,1),
         lwd=c(2.5,2.5),
         col=c("blue","red")
  )
}
 
###########################################
#Naive 2
Naive2.PoissonProcess <- function(Opti,k){
  L <- rep(1,length(Opti[,2])-1)
  I.poisprocess <- lambda
  L <- L*floor(k/I.poisprocess)
  L <- c(L, k- sum(L))
  L <- sum(Cost.Function(L))
  O <- sum(Cost.Function(Opti[,2]))+Cost.Function(Opti[length(Opti[,1]),1])

  return((L-O)/L)
}
Naive2.PoissonProcess.MonteCarlo <- function(k){
  a <- result.read("a.txt")
  X <- c()
  for(i in (1:1000)){
    OrderBook <- PoissonProcess(1, lambda)
    Opti <- OptimalStrategy.PoissonProcess(a, OrderBook, k)
    X <- c(X,Naive2.PoissonProcess(Opti,k))
  }
  X <- X[which(!is.na(X))]
  X <- X/Cost.Function(k)
  vari <- sqrt(var(X))
  beta <- 1.96 
  return(c((sum(X)/length(X)),beta*vari/sqrt(1000)))
}
Naive2.CompoundPoissonProcess <- function(Opti,k){
  L <- rep(1,length(Opti[,2])-1)
  I <- floor(lambda*mu/(1-exp(-mu)))
  L <- L*floor(k/I)+1
  L <- c(L, k- sum(L))
  L <- sum(Cost.Function(L))
  O <- sum(Cost.Function(Opti[,2]))+Cost.Function(Opti[length(Opti[,1]),1])
  return((L-O)/L)
}
Naive2.CompoundPoissonProcess.MonteCarlo <- function(k){
  a <- result.read("apois.txt")
  X <- c()
  for(i in (1:1000)){
    OrderBook <- CompoundPoissonProcess(1,mu, lambda)
    Opti <- OptimalStrategy.CompoundPoissonProcess(a, OrderBook, k)
    X <- c(X,Naive2.CompoundPoissonProcess(Opti,k))
  }
  X <- X[which(!is.na(X))]
  X <- X/Cost.Function(k)
  vari <- sqrt(var(X))
  beta <- 1.96 
  return(c((sum(X)/length(X)),beta*vari/sqrt(1000)))
}

# OrderBook <- PoissonProcess(1, lambda)
# Opti <- OptimalStrategy.PoissonProcess(a, OrderBook, k)
# Opti
# L <- rep(1,length(Opti[,2])-1)
# A <-seq(1,18)
# sapply(A,Naive2.CompoundPoissonProcess.MonteCarlo)


