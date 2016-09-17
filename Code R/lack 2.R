# init
library("expm")
lambda <- c(3, 5, 10)
mu <- c(8,4,4)
Gamma <- diag(lambda)
Q <- t(matrix(data =c(-2,2,0,1,-4,3,0,2,-2),nrow = 3, ncol = 3))
Q <- t(Q)


# random generation
gener <- function(x){
  return (-log(runif(1))/x);
}

# proba
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

PoissonState <- function(Q,i,Time){
  X <- gener(abs(Q[i,i]))
  if( X < Time){
    return(X)
  }
  else
    return(Time)
}


GenerateOrderBook <- function(Q,Time){
  regime <- c(1)
  OrderBook <- c(PoissonState(Q,regime[length(regime)],Time))
  while(OrderBook[length(OrderBook)] < Time ){
    regime <- c(regime,proba(Q,regime[length(regime)]))
    OrderBook <- c(OrderBook,OrderBook[length(OrderBook)] +
                     PoissonState(Q,regime[length(regime)],Time))
  }
  return(cbind(OrderBook,regime))
}

GenerateTime <- function(States, Time, Gamma){
  regime <- 1
  sigma <- c(PoissonState(Gamma, regime,Time))
  ordersize <- c(rpois(1,lambda[regime]))
  while(sigma[length(sigma)] < Time){
    regime <- States[which(States[,1]>sigma[length(sigma)])[1],2]
    sigma <- c(sigma,sigma[length(sigma)] + PoissonState(Gamma, regime,Time)) 
    ordersize <- c(ordersize, rpois(1,lambda[regime]))
  }
  return(cbind(sigma[1:length(sigma)-1],ordersize[1:length(ordersize)-1]))
}

Order <- GenerateOrderBook(Q,20)
Sigma <- GenerateTime(Order, 20, Gamma)


# Actualisation de linfo
m <- function(t, pi){
  return(pi%*%expm(t*(Q-Gamma)))
}
x <- function(t, pi){
  return((1/sum(m(t,pi)))*m(t,pi))
}

nu <- function(y,i,mu){
  return ((exp(-mu[i])*mu[i]^y)/((1-exp(-mu[i]))*factorial(y)))
}

jump <- function(X,sigma,j){
  Y <- nu(sigma[j,2],seq(1:length(lambda)),mu)*X
  Y <- lambda*Y
  Y <- Y/sum(Y)
  return(Y)
}

LackInfo.Test <- function(){
  pi = c(0.01,0.2,0.7)
  OB <- Sigma
  OB <- cbind(c(0,OB[,1]),c(0,OB[,2]))
  X <- c(pi)
  for (i in (2:length(OB[,1]))){
    pi <- x(OB[i]-OB[i-1],pi)
    X <- rbind(X,pi)
    pi <- jump(pi,OB,i)
    X <- rbind(X,pi)
  }
  
  plot(X[,1],type='l')
  points(X[,2],type='l', col ="red")
  points(X[,3],type='l', col = "blue")
}
