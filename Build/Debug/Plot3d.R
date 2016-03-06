
# init directory
setwd(path)
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
  z <- data
  # plot 3d
  RemainingShares <- (1:nrow(z))   # 10 meter spacing (S to N)
  TimeToMaturity <-  (1:ncol(z))   # 10 meter spacing (E to W)
  ## Don't draw the grid lines :  border = NA
  #persp(RemainingShares, TimeToMaturity, z,theta = -140,phi = 25)
  plot_ly(z = z, type = "surface")
}
resultplot("a.txt")
resultplot('v.txt')
resultplot('vpoisscomp.txt')
resultplot("apoisscomp.txt")

v <- result.read('v.txt')
vp <- result.read('vpoisscomp.txt')
which(vp-v>0)
quartz()
resultplot("apoisscomp.txt")
