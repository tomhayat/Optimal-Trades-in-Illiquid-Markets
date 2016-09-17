setwd("/Users/tom/Google Drive/Projet Mopsi/Code C++/Build/Debug")
k <- 10
lambda <- 4
mu <- 8
Time <- 1
Q <- t(matrix(data =c(-2,2,0,1,-4,3,0,2,-2),nrow = 3, ncol = 3))
Q <- t(Q)
lambda.Markov <- 4*c(3, 3, 1)
mu.Markov <- c(8,4,4)
Gamma <- diag(lambda.Markov)
# Poisson Process
resultplot("aTEST.txt")
resultplot("vTEST.txt")
test.PoissonProcess()
a <- result.read("as1TEST.txt")
Example.PoissonProcess(k,a)

####
# Compound Poisson Process
resultplot("apoisTEST.txt")
resultplot("vpoisTEST.txt")
test.CompoundPoissonProcess()
apois <- result.read("apoisTEST.txt")
Example.CompoundPoissonProcess(k,apois)

# RSS
resultplot("Vs1TEST.txt")
resultplot("Vs2TEST.txt")
resultplot("Vs3TEST.txt")
resultplot("as1TEST.txt")
resultplot("as2TEST.txt")
resultplot("as3TEST.txt")
as1 <- result.read("as1TEST.txt")
as2 <- result.read("as2TEST.txt")
as3 <- result.read("as3TEST.txt")
as <- list(as1,as2,as3)
test.RegimeSwitching()
Example.RSS(k,as)




