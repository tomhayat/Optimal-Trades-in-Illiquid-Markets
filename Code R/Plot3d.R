# init directory
#setwd("/Users/danconstantini/Documents/IMI/MOPSI/Projet/TP1/build/Debug")
#read the table and adapt the type
library("plotly")
#read the table and adapt the type
result.read <- function(filename){
  data <- read.csv(filename,check.names = FALSE,header = FALSE)
  data <- matrix(unlist(data), nrow = dim(data)[1],ncol = dim(data)[2])
  return(data)
}

resultplot <- function(filename){
  data <- read.csv(filename,check.names = FALSE,header = FALSE)
  data <- matrix(unlist(data), nrow = dim(data)[1],ncol = dim(data)[2])
  z <- data[c(1:dim(data)[1]),]
  # plot 3d
  RemainingShares <- (1:nrow(z))   # 10 meter spacing (S to N)
  TimeToMaturity <-  (1:ncol(z))  # 10 meter spacing (E to W)
  ## Don't draw the grid lines :  border = NA
  #persp(RemainingShares, TimeToMaturity, z,theta = -140,phi = 25)
  plot_ly(z = z, type = "surface")
}










