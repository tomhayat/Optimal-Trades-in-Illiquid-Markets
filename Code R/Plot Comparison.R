Plot <- function(){
  #setwd("/Users/tom/Downloads/TP1/Build/Debug")
  list.files()
  a <- result.read("asc.txt")
  v <- result.read("vsc.txt")
  ap <- result.read("apoisc.txt")
  vp <- result.read("vpoisc.txt")
  ac <- result.read("acontinuous.txt")
  plot(ap[200,]+1~seq(0,0.995, by = 0.005),cex= 0.1,col = "red",     main = "Model Comparison for k = 19",
       xlab =  "Time",
       ylab = "Optimal Cost Value")
  points(seq(0,0.995, by = 0.005),a[200,]+1.0125,cex= 0.1, col ="blue")
  points(seq(0,0.995, by = 0.005),ap[200,]+1,cex= 0.1, col ="red")
  points(seq(0,1, by = 0.005),200*ac,cex= 0.1, col ="green")
  
  plot(vp[5,]~seq(0,0.995, by = 0.005),cex= 0.1,col = "red",
       main = "Model Comparison for k = 19",
       xlab =  "Time",
       ylab = "Optimal Execution Cost")
  points(seq(0,0.995, by = 0.005),v[5,],cex= 0.1, col ="blue")
  points(seq(0,0.995, by = 0.005),vp[5,],cex= 0.1, col ="red")
  points(seq(0,1, by = 0.005),5^4*(ac),cex= 0.1, col ="green")
  plot(3^4*(ac)~seq(0,1, by = 0.005))
  plot(v[5,]~seq(0,0.995, by = 0.005),cex= 0.1,col = "red",
       main = "Model Comparison for k = 19",
       xlab =  "Time",
       ylab = "Optimal Execution Cost")
  v[5,]-vp[5,]
  plot(vp[20,] - ac*20^4)
}


