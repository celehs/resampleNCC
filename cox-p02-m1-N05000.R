library(MASS)
library(mvtnorm)
library(survival)
library(glmpath)
library(lars)

source("R/SIM.R")
source("R/FUN.R")

set.seed(123)
N <- 05000
p.y <- 02
gam.true <- c(1, 0.5)
myB.dist <- "exp"
yes.CLR <- TRUE
yes.DepCen <- TRUE
yes.Z <- FALSE
m.match <- 1
a0.match <- c(0, 1)
yes.regularize <- (p.y > length(gam.true))
myind.Zmatch <- 1:length(a0.match) + 3 + p.y
yes.betaonly <- FALSE
yes.pos.only <- FALSE
t0 <- seq(0.1, 1.5, by = 0.01) 
n.t0 <- length(t0)
yy0 <- matrix(rep(c(0, 1), rep(p.y, 2)), ncol = 2)

ind.remove <- TRUE
count0 <- 0
while (ind.remove) {
  mydata.cht <- SIM.FUN(nn = N, DepCen = yes.DepCen, Zmatch = yes.Z, 
                        a0 = a0.match, p.add.noise = p.y - length(gam.true))
  mydata.ncc <- mydata.cht[[1]]
  Vij.IND <- mydata.cht[[2]] 
  Iik0.mat <- mydata.cht[[3]]
  ind.remove <- sum(Vij.IND[, 2] == 0) > 0 
  count0 <- count0 + 1 
  if(count0 > 1) { print(count0 - 1) }
}
# Conditional Logistic Regression
if (yes.CLR) {
  junk1 <- NCC.Cox.FUN.clogit(
    mydata.ncc, Vij.IND, y0.cut = yy0, Zmatch.ind = myind.Zmatch)
  bhat1 <- junk1[-(1:(n.t0 * ncol(yy0))), ]
} else {
  bhat1 <- junk1 <- NULL
}
# IPW Resampling Approach
junk2 <- NCC.Cox.FUN.Regularize(
  mydata.ncc, Vij.IND, Iik0 = Iik0.mat, Zmatch.ind = myind.Zmatch, 
  rtn = "ALL", regularize = yes.regularize, y0.cut = yy0, 
  betaonly = yes.betaonly, B0 = 1000, pos.only = yes.pos.only)
if (yes.betaonly) {
  bhat2 <- matrix(junk2, ncol = 3)
} else { 
  bhat2 <- matrix(junk2[-(1:ncol(yy0))], ncol = 3)[-(1:(n.t0 * ncol(yy0))), ]
}
colnames(bhat2) <- c("Est", "SE1", "SE2")

# CLR Results
bhat1

# Resampling Results
bhat2
