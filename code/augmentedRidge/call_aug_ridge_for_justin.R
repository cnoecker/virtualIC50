
source("./aridge/aug_ridge_code_for_justin.R")


burnin <- 3000
my.ylim = c(-5, 5)

################################################################
## simulate data
################################################################

set.seed(12345)

n1 <- 50
true.beta <- c(rnorm(50), rep(0, 50))
sig <- 1
rho <- 0.5
dat1 <- SimulateData(2 * n1, true.beta, sig, rho, 100, 20)
p <- length(true.beta)

y1train <- dat1[[1]][1:n1]
X1train <- dat1[[2]][1:n1,]
y1test <- dat1[[1]][(n1+1):(2*n1)]
X1test <- dat1[[2]][(n1+1):(2*n1),]

X1train <- as.matrix(scale(X1train))
X1test <- as.matrix(scale(X1test))
y1train <- y1train - mean(y1train)
y1test <- y1test - mean(y1test) 



ed1 <- eigen(t(X1train) %*% X1train, symmetric = TRUE)
P1 <- ed1$vectors
d1 <- ed1$values
Z1s <- X1train %*% P1
Z1sty <- as.vector(t(Z1s) %*% y1train)


########################################################################
## fit Bayesian ridge restricted to the first data set (fixed lambda)
########################################################################


o2 <- GibbsRidgeEigen(y1train,
                      Z1s,
                      Z1sty, 
                      d1,
                      ite = 13000,
                      r.tau = 0.01, 
                      s.tau = 0.01,
                      lambda = 10,
                      thetaj.0 = rep(0, length(d1)),   
                      tau.0 = 1)
Betah2 <- P1 %*% o2$Thetaj[, -c(1:burnin)]
m2 <- apply(Betah2, 1, mean)
GetMse(X = X1test, Y = y1test, m2)
su2 <- GetPostProbInterval(Betah2, 0.99, true.beta)
su2$correct.selections
PlotPostInterval(su2$posterior.interval, ylim = my.ylim)


fsvd1 <- svd(X1train)
V1 <- fsvd1$v
ZZ1s <- X1train %*% V1
ZZ1sty <- t(ZZ1s) %*% y1train
dd1 <- fsvd1$d
oo2 <- GibbsRidgeSvd(y1train,
                     ZZ1s,
                     ZZ1sty, 
                     dd1,
                     ite = 13000,
                     r.tau = 0.01, 
                     s.tau = 0.01,
                     lambda = 10,
                     thetaj.0 = rep(0, length(dd1)),   
                     tau.0 = 1)
BBetah2 <- V1 %*% oo2$Thetaj[, -c(1:burnin)]
mm2 <- apply(BBetah2, 1, mean)
GetMse(X = X1test, Y = y1test, mm2)
ssu2 <- GetPostProbInterval(BBetah2, 0.99, true.beta)
ssu2$correct.selections
PlotPostInterval(ssu2$posterior.interval, ylim = my.ylim)



X1ty <- t(X1train) %*% y1train
X1tX1 <- t(X1train) %*% X1train
ooo2 <- GibbsRidge(y1train,
                   X1train,
                   X1ty, 
                   X1tX1,
                   ite = 13000,
                   r.tau = 0.01, 
                   s.tau = 0.01,
                   lambda = 10,
                   betaj.0 = rep(0, length(p)),   
                   tau.0 = 1)
BBBetah2 <- ooo2$Betaj[, -c(1:burnin)]
mmm2 <- apply(BBBetah2, 1, mean)
GetMse(X = X1test, Y = y1test, mmm2)
sssu2 <- GetPostProbInterval(BBBetah2, 0.99, true.beta)
sssu2$correct.selections
PlotPostInterval(sssu2$posterior.interval, ylim = my.ylim)

my.ylim <- c(-5, 5)
par(mfrow = c(1, 3))
PlotPostInterval(ssu2$posterior.interval, ylim = my.ylim, main = "SVD")
PlotPostInterval(su2$posterior.interval, ylim = my.ylim, main = "Eigen")
PlotPostInterval(sssu2$posterior.interval, ylim = my.ylim, main = "original")
par(mfrow = c(1, 1))




########################################################################
## fit Bayesian ridge restricted to the first data set (full model)
########################################################################


o2 <- GibbsRidgeEigenFull(y1train,
                          Z1s,
                          Z1sty, 
                          d1,
                          ite = 13000,
                          r.tau = 0.01, 
                          s.tau = 0.01,
                          r.lambda = 0.01,
                          s.lambda = 0.01,
                          thetaj.0 = rep(0, length(d1)),   
                          tau.0 = 1,
                          lambda.0 = 1)
Betah2 <- P1 %*% o2$Thetaj[, -c(1:burnin)]
m2 <- apply(Betah2, 1, mean)
GetMse(X = X1test, Y = y1test, m2)
su2 <- GetPostProbInterval(Betah2, 0.99, true.beta)
su2$correct.selections
PlotPostInterval(su2$posterior.interval, ylim = my.ylim)


fsvd1 <- svd(X1train)
V1 <- fsvd1$v
ZZ1s <- X1train %*% V1
ZZ1sty <- t(ZZ1s) %*% y1train
dd1 <- fsvd1$d
oo2 <- GibbsRidgeSvdFull(y1train,
                         ZZ1s,
                         ZZ1sty, 
                         dd1,
                         ite = 13000,
                         r.tau = 0.01, 
                         s.tau = 0.01,
                         r.lambda = 0.01,
                         s.lambda = 0.01,
                         thetaj.0 = rep(0, length(dd1)),   
                         tau.0 = 1,
                         lambda.0 = 1)
BBetah2 <- V1 %*% oo2$Thetaj[, -c(1:burnin)]
mm2 <- apply(BBetah2, 1, mean)
GetMse(X = X1test, Y = y1test, mm2)
ssu2 <- GetPostProbInterval(BBetah2, 0.99, true.beta)
ssu2$correct.selections
PlotPostInterval(ssu2$posterior.interval, ylim = my.ylim)



X1ty <- t(X1train) %*% y1train
X1tX1 <- t(X1train) %*% X1train
ooo2 <- GibbsRidgeFull(y1train,
                       X1train,
                       X1ty, 
                       X1tX1,
                       ite = 13000,
                       r.tau = 0.01, 
                       s.tau = 0.01,
                       r.lambda = 0.01,
                       s.lambda = 0.01,
                       betaj.0 = rep(0, length(p)),   
                       tau.0 = 1,
                       lambda.0 = 1)
BBBetah2 <- ooo2$Betaj[, -c(1:burnin)]
mmm2 <- apply(BBBetah2, 1, mean)
GetMse(X = X1test, Y = y1test, mmm2)
sssu2 <- GetPostProbInterval(BBBetah2, 0.99, true.beta)
sssu2$correct.selections
PlotPostInterval(sssu2$posterior.interval, ylim = my.ylim)

my.ylim <- c(-5, 5)
par(mfrow = c(1, 3))
PlotPostInterval(ssu2$posterior.interval, ylim = my.ylim, main = "SVD")
PlotPostInterval(su2$posterior.interval, ylim = my.ylim, main = "Eigen")
PlotPostInterval(sssu2$posterior.interval, ylim = my.ylim, main = "original")
par(mfrow = c(1, 1))



################################################################
## larger p case
################################################################


set.seed(12345)

n1 <- 300
true.beta <- c(rnorm(250), rep(0, 250))
sig <- 1
rho <- 0.5
dat1 <- SimulateData(2 * n1, true.beta, sig, rho, 100, 20)
p <- length(true.beta)

y1train <- dat1[[1]][1:n1]
X1train <- dat1[[2]][1:n1,]
y1test <- dat1[[1]][(n1+1):(2*n1)]
X1test <- dat1[[2]][(n1+1):(2*n1),]

X1train <- as.matrix(scale(X1train))
X1test <- as.matrix(scale(X1test))
y1train <- y1train - mean(y1train)
y1test <- y1test - mean(y1test) 


ed1 <- eigen(t(X1train) %*% X1train, symmetric = TRUE)
P1 <- ed1$vectors
d1 <- ed1$values
Z1s <- X1train %*% P1
Z1sty <- as.vector(t(Z1s) %*% y1train)
o2 <- GibbsRidgeEigenFull(y1train,
                          Z1s,
                          Z1sty, 
                          d1,
                          ite = 13000,
                          r.tau = 0.01, 
                          s.tau = 0.01,
                          r.lambda = 0.01,
                          s.lambda = 0.01,
                          thetaj.0 = rep(0, length(d1)),   
                          tau.0 = 1,
                          lambda.0 = 1)
Betah2 <- P1 %*% o2$Thetaj[, -c(1:burnin)]
m2 <- apply(Betah2, 1, mean)
GetMse(X = X1test, Y = y1test, m2)
cor(y1test, X1test %*% m2)
su2 <- GetPostProbInterval(Betah2, 0.99, true.beta)
su2$correct.selections
PlotPostInterval(su2$posterior.interval, ylim = my.ylim)

mean(o2$Tau[-c(1:burnin)])
mean(o2$Lambda[-c(1:burnin)])

fsvd1 <- fast.svd(X1train)
V1 <- fsvd1$v
ZZ1s <- X1train %*% V1
ZZ1sty <- t(ZZ1s) %*% y1train
dd1 <- fsvd1$d
oo2 <- GibbsRidgeSvdFull(y1train,
                         ZZ1s,
                         ZZ1sty, 
                         dd1,
                         ite = 13000,
                         r.tau = 0.01, 
                         s.tau = 0.01,
                         r.lambda = 0.01,
                         s.lambda = 0.01,
                         thetaj.0 = rep(0, length(dd1)),   
                         tau.0 = 1,
                         lambda.0 = 1)
BBetah2 <- V1 %*% oo2$Thetaj[, -c(1:burnin)]
mm2 <- apply(BBetah2, 1, mean)
GetMse(X = X1test, Y = y1test, mm2)
ssu2 <- GetPostProbInterval(BBetah2, 0.99, true.beta)
ssu2$correct.selections
PlotPostInterval(ssu2$posterior.interval, ylim = my.ylim)

mean(oo2$Tau[-c(1:burnin)])
mean(oo2$Lambda[-c(1:burnin)])

v2 <- apply(Betah2, 1, var)
vv2 <- apply(BBetah2, 1, var)
summary(v2)
summary(vv2)

plot(m2, mm2)

###########################################################################
## simulate second data-set
###########################################################################


set.seed(54321)
n2 <- 500
dat2 <- SimulateData(2 * n2, true.beta, sig, rho, 100, 20)
y2train <- dat2[[1]][1:n2]
X2train <- dat2[[2]][1:n2,]
y2train <- y2train - mean(y2train)
X2train <- as.matrix(scale(X2train))



###########################################################################
## fit augmented Bayesian ridge on y1, X1, X2, but where y2 is missing 
###########################################################################

Xtrain <- rbind(X1train, X2train)

ed <- eigen(t(Xtrain) %*% Xtrain, symmetric = TRUE)
P <- ed$vectors
d <- ed$values
Z <- Xtrain %*% P
Z1 <- Z[1:n1,]
Z2 <- Z[-c(1:n1),]
Z1ty <- as.vector(t(Z1) %*% y1train)
theta.hat <- apply(o2$Thetaj[, -c(1:burnin)], 1, mean)
tau.hat <- mean(o2$Tau[-c(1:burnin)])
o3 <- GibbsAugmentedRidgeEigenFull(y1train, 
                                   Z2,
                                   Z,
                                   Z1ty,
                                   P1,
                                   Pt = t(P),
                                   theta.hat,
                                   tau.hat,
                                   d,
                                   ite = 13000,
                                   r.tau = 0.01, 
                                   s.tau = 0.01,
                                   r.lambda = 0.01,
                                   s.lambda = 0.01,
                                   thetaj.0 = rep(0, length(d)),   
                                   tau.0 = 1,
                                   lambda.0 = 1,
                                   w.0 = rnorm(n2))
Betah3 <- P %*% o3$Thetaj[, -c(1:burnin)]
m3 <- apply(Betah3, 1, mean)
GetMse(X = X1test, Y = y1test, m3)
fsu3 <- GetPostProbInterval(Betah3, 0.9, true.beta)
fsu3$correct.selections
par(mfrow = c(1, 1))
PlotPostInterval(fsu3$posterior.interval, ylim = my.ylim)




###########################################################################
## augmented fits using the posterior predictive distribution
###########################################################################



vis <- GetVs(Z1, Z2, mean(o2$Lambda[-c(1:burnin)]))

s.hat <- (0.01 + t(y1test - Z1 %*% theta.hat) %*% y1test)/(0.01 + n1)

o4 <- GibbsAugmentedRidgeEigenPPD(y1train, 
                                  Z2,
                                  Z,
                                  Z1ty,
                                  P1,
                                  Pt = t(P),
                                  theta.hat = theta.hat,
                                  s.hat = s.hat,
                                  vis = rep(mean(vis), n2),
                                  d = d,
                                  ite = 13000,
                                  a.tau = 0.01, 
                                  b.tau = 0.01,
                                  lambda = mean(o2$Lambda[-c(1:burnin)]),
                                  thetaj.0 = rep(0, length(d)),   
                                  tau.0 = 1,
                                  w.0 = rnorm(n2),
                                  etai.0 = rep(1, n2))
Betah4 <- P %*% o4$Thetaj[, -c(1:burnin)]
m4 <- apply(Betah4, 1, mean)
GetMse(X = X1test, Y = y1test, m4)
cor(y1test, X1test %*% m4)
fsu4 <- GetPostProbInterval(Betah4, 0.9999, true.beta)
fsu4$correct.selections
par(mfrow = c(1, 1))
PlotPostInterval(fsu4$posterior.interval, ylim = my.ylim)



oo4 <- GibbsAugmentedRidgeEigenPPD2(y1train, 
                                    Z2,
                                    Z,
                                    Z1ty,
                                    P1,
                                    Pt = t(P),
                                    theta.hat = theta.hat,
                                    s.hat = s.hat,
                                    vis = rep(mean(vis), n2),
                                    d = d,
                                    ite = 13000,
                                    a.tau = 0.01, 
                                    b.tau = 0.01,
                                    lambda = mean(o2$Lambda[-c(1:burnin)]),
                                    thetaj.0 = rep(0, length(d)),   
                                    tau.0 = 1,
                                    w.0 = rnorm(n2),
                                    etai.0 = rep(1, n2),
                                    thin.window = 10)
BBetah4 <- P %*% oo4$Thetaj[, -c(1:burnin)]
mm4 <- apply(BBetah4, 1, mean)
GetMse(X = X1test, Y = y1test, mm4)
cor(y1test, X1test %*% mm4)
ffsu4 <- GetPostProbInterval(BBetah4, 0.9999, true.beta)
ffsu4$correct.selections
par(mfrow = c(1, 1))
PlotPostInterval(ffsu4$posterior.interval, ylim = my.ylim)


ooo4 <- GibbsAugmentedRidgeEigenPPD2(y1train, 
                                     Z2,
                                     Z,
                                     Z1ty,
                                     P1,
                                     Pt = t(P),
                                     theta.hat = theta.hat,
                                     s.hat = s.hat,
                                     vis = rep(mean(vis), n2),
                                     d = d,
                                     ite = 13000,
                                     a.tau = 0.01, 
                                     b.tau = 0.01,
                                     lambda = mean(o2$Lambda[-c(1:burnin)]),
                                     thetaj.0 = rep(0, length(d)),   
                                     tau.0 = 1,
                                     w.0 = rnorm(n2),
                                     etai.0 = rep(1, n2),
                                     thin.window = 20)
BBBetah4 <- P %*% ooo4$Thetaj[, -c(1:burnin)]
mmm4 <- apply(BBBetah4, 1, mean)
GetMse(X = X1test, Y = y1test, mmm4)
cor(y1test, X1test %*% mmm4)
fffsu4 <- GetPostProbInterval(BBBetah4, 0.9999, true.beta)
fffsu4$correct.selections
par(mfrow = c(1, 1))
PlotPostInterval(fffsu4$posterior.interval, ylim = my.ylim)


k <- 1
plot(Betah2[k,], type = "n")
lines(Betah2[k,])

plot(Betah4[k,], type = "n")
lines(Betah4[k,])

plot(BBetah4[k,], type = "n")
lines(BBetah4[k,])

plot(BBBetah4[k,], type = "n")
lines(BBBetah4[k,])


plot(m2, m4)
plot(m2, mm4)
plot(m2, mmm4)

