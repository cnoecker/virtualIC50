
library(MASS)
library(corpcor)

########################################################################
## functions to simulate data
########################################################################


CreateSigma <- function(rho, p) {
  aux1 <- matrix(rep(1:p, p), p, p)
  aux2 <- matrix(rep(1:p, each = p), p, p) 
  rho^abs(aux1 - aux2)
}



GetBlockSizes <- function(p, max.block, min.block) {
  max.block <- min(c(p, max.block))
  block.sizes <- sample(min.block:max.block, ceiling(p/min.block), 
                        replace = TRUE)
  cum.sizes <- cumsum(block.sizes)
  block.sizes <- block.sizes[which(cum.sizes <= p)]
  aux <- which.min(block.sizes)
  block.sizes[aux] <- block.sizes[aux] + p - sum(block.sizes)
  block.sizes 
}



SimulateData <- function(n, beta, sig, rho, max.block, min.block) {
  p <- length(beta)
  block.sizes <- GetBlockSizes(p, max.block, min.block)
  nblocks <- length(block.sizes)
  Sigma <- CreateSigma(rho, block.sizes[1])
  X <- mvrnorm(n, rep(0, block.sizes[1]), Sigma)
  if (nblocks > 1) {
    for (i in 2:nblocks) {
      Sigma <- CreateSigma(rho, block.sizes[i])
      X <- cbind(X, mvrnorm(n, rep(0, block.sizes[i]), Sigma))  
    }
  }
  y <- X %*% beta + sig * rnorm(n)
  list(y = y, X = X)
}


########################################################################
## functions to fit Bayesian ridge and augmented Bayesian ridge 
########################################################################


######################################
## fixed lambda
######################################

GibbsRidgeEigen <- function(y,
                            Z,
                            Zty, 
                            d,
                            ite,
                            r.tau, 
                            s.tau,
                            lambda,
                            thetaj.0,   
                            tau.0) {
  UpdateThetaj <- function(r, Zty, d, lambda, tau) {
    aux <- 1/(d + lambda)
    rnorm(r, aux * Zty, sqrt(aux/tau))
  }
  UpdateTau <- function(y, Z, r, n, r.tau, s.tau, lambda, theta) {
    rgamma(1, shape = r.tau + (r + n)/2, 
           rate = s.tau + 0.5 * (sum(theta^2 * lambda) + 
                                   sum((y - Z %*% theta)^2)))
  }
  r <- length(d)
  n <- length(y)
  Thetaj <- matrix(NA, r, ite)
  Tau <- rep(NA, ite)
  Thetaj[, 1] <- thetaj.0
  Tau[1] <- tau.0
  ##
  for(i in 2:ite) {
    ## update thetaj: r, Zty, d, lambda, tau
    Thetaj[,i] <- UpdateThetaj(r, Zty, d, lambda, Tau[i-1])
    ## update tau: y, Z, r, n, r.tau, s.tau, lambda, theta
    Tau[i] <- UpdateTau(y, Z, r, n, r.tau, s.tau, lambda, Thetaj[,i])
    cat("ite", i, "\n")
  }
  list(Thetaj = Thetaj, Tau = Tau)
}



GibbsRidgeSvd <- function(y,
                          Z,
                          Zty, 
                          d,
                          ite,
                          r.tau, 
                          s.tau,
                          lambda,
                          thetaj.0,   
                          tau.0) {
  UpdateThetaj <- function(r, Zty, d, lambda, tau) {
    aux <- 1/(d^2 + lambda)
    rnorm(r, aux * Zty, sqrt(aux/tau))
  }
  UpdateTau <- function(y, Z, r, n, r.tau, s.tau, lambda, theta) {
    rgamma(1, shape = r.tau + (r + n)/2, 
           rate = s.tau + 0.5 * (sum(theta^2 * lambda) + 
                                   sum((y - Z %*% theta)^2)))
  }
  r <- length(d)
  n <- length(y)
  Thetaj <- matrix(NA, r, ite)
  Tau <- rep(NA, ite)
  Thetaj[, 1] <- thetaj.0
  Tau[1] <- tau.0
  ##
  for(i in 2:ite) {
    ## update thetaj: r, Zty, d, lambda, tau
    Thetaj[,i] <- UpdateThetaj(r, Zty, d, lambda, Tau[i-1])
    ## update tau: y, Z, r, n, r.tau, s.tau, lambda, theta
    Tau[i] <- UpdateTau(y, Z, r, n, r.tau, s.tau, lambda, Thetaj[,i])
    cat("ite", i, "\n")
  }
  list(Thetaj = Thetaj, Tau = Tau)
}



GibbsRidge <- function(y,
                       X,
                       Xty,
                       XtX,
                       ite,
                       r.tau, 
                       s.tau,
                       lambda,
                       betaj.0,   
                       tau.0) {
  UpdateBetaj <- function(r, XtX, Xty, lambda, tau) {
    aux <- solve(XtX + lambda * diag(r))
    mvrnorm(1, aux %*% Xty, aux/tau)
  }
  UpdateTau <- function(y, X, r, n, r.tau, s.tau, lambda, beta) {
    rgamma(1, shape = r.tau + (r + n)/2, 
           rate = s.tau + 0.5 * (sum(beta^2 * lambda) + 
                                   sum((y - X %*% beta)^2)))
  }
  r <- ncol(X)
  n <- length(y)
  Betaj <- matrix(NA, r, ite)
  Tau <- rep(NA, ite)
  Betaj[, 1] <- betaj.0
  Tau[1] <- tau.0
  ##
  for(i in 2:ite) {
    ## update betaj:
    Betaj[,i] <- UpdateBetaj(r, XtX, Xty, lambda, Tau[i-1])
    ## update tau:
    Tau[i] <- UpdateTau(y, X, r, n, r.tau, s.tau, lambda, Betaj[,i])
    cat("ite", i, "\n")
  }
  list(Betaj = Betaj, Tau = Tau)
}



######################################
## full model (including lambda)
######################################

GibbsRidgeEigenFull <- function(y,
                                Z,
                                Zty, 
                                d,
                                ite,
                                r.tau, 
                                s.tau,
                                r.lambda,
                                s.lambda,
                                thetaj.0,   
                                tau.0,
                                lambda.0) {
  UpdateThetaj <- function(r, Zty, d, lambda, tau) {
    aux <- 1/(d + lambda)
    rnorm(r, aux * Zty, sqrt(aux/tau))
  }
  UpdateTau <- function(y, Z, r, n, r.tau, s.tau, lambda, theta) {
    rgamma(1, shape = r.tau + (r + n)/2, 
           rate = s.tau + 0.5 * (sum(theta^2 * lambda) + 
                                   sum((y - Z %*% theta)^2)))
  }
  UpdateLambda <- function(r, r.lambda, s.lambda, tau, theta) {
    rgamma(1, shape = r.lambda + r/2, 
           rate = s.lambda + 0.5 * tau * sum(theta^2))
  }
  r <- length(d)
  n <- length(y)
  Thetaj <- matrix(NA, r, ite)
  Tau <- rep(NA, ite)
  Lambda <- rep(NA, ite)
  Thetaj[, 1] <- thetaj.0
  Tau[1] <- tau.0
  Lambda[1] <- lambda.0
  ##
  for(i in 2:ite) {
    ## update thetaj: r, Zty, d, lambda, tau
    Thetaj[,i] <- UpdateThetaj(r, Zty, d, Lambda[i-1], Tau[i-1])
    ## update tau: y, Z, r, n, r.tau, s.tau, lambda, theta
    Tau[i] <- UpdateTau(y, Z, r, n, r.tau, s.tau, Lambda[i-1], Thetaj[,i])
    Lambda[i] <- UpdateLambda(r, r.lambda, s.lambda, Tau[i], Thetaj[,i])
    cat("ite", i, "\n")
  }
  list(Thetaj = Thetaj, Tau = Tau, Lambda = Lambda)
}



GibbsRidgeSvdFull <- function(y,
                              Z,
                              Zty, 
                              d,
                              ite,
                              r.tau, 
                              s.tau,
                              r.lambda,
                              s.lambda,
                              thetaj.0,   
                              tau.0,
                              lambda.0) {
  UpdateThetaj <- function(r, Zty, d, lambda, tau) {
    aux <- 1/(d^2 + lambda)
    rnorm(r, aux * Zty, sqrt(aux/tau))
  }
  UpdateTau <- function(y, Z, r, n, r.tau, s.tau, lambda, theta) {
    rgamma(1, shape = r.tau + (r + n)/2, 
           rate = s.tau + 0.5 * (sum(theta^2 * lambda) + 
                                   sum((y - Z %*% theta)^2)))
  }
  UpdateLambda <- function(r, r.lambda, s.lambda, tau, theta) {
    rgamma(1, shape = r.lambda + r/2, 
           rate = s.lambda + 0.5 * tau * sum(theta^2))
  }
  r <- length(d)
  n <- length(y)
  Thetaj <- matrix(NA, r, ite)
  Tau <- rep(NA, ite)
  Lambda <- rep(NA, ite)
  Thetaj[, 1] <- thetaj.0
  Tau[1] <- tau.0
  Lambda[1] <- lambda.0
  ##
  for(i in 2:ite) {
    ## update thetaj: r, Zty, d, lambda, tau
    Thetaj[,i] <- UpdateThetaj(r, Zty, d, Lambda[i-1], Tau[i-1])
    ## update tau: y, Z, r, n, r.tau, s.tau, lambda, theta
    Tau[i] <- UpdateTau(y, Z, r, n, r.tau, s.tau, Lambda[i-1], Thetaj[,i])
    ## update lambda: r, r.lambda, s.lambda, tau, thetaj
    Lambda[i] <- UpdateLambda(r, r.lambda, s.lambda, Tau[i], Thetaj[,i])
    cat("ite", i, "\n")
  }
  list(Thetaj = Thetaj, Tau = Tau, Lambda = Lambda)
}



GibbsRidgeFull <- function(y,
                           X,
                           Xty,
                           XtX,
                           ite,
                           r.tau, 
                           s.tau,
                           r.lambda,
                           s.lambda,
                           betaj.0,   
                           tau.0,
                           lambda.0) {
  UpdateBetaj <- function(r, XtX, Xty, lambda, tau) {
    aux <- solve(XtX + lambda * diag(r))
    mvrnorm(1, aux %*% Xty, aux/tau)
  }
  UpdateTau <- function(y, X, r, n, r.tau, s.tau, lambda, beta) {
    rgamma(1, shape = r.tau + (r + n)/2, 
           rate = s.tau + 0.5 * (sum(beta^2 * lambda) + 
                                   sum((y - X %*% beta)^2)))
  }
  UpdateLambda <- function(r, r.lambda, s.lambda, tau, beta) {
    rgamma(1, shape = r.lambda + r/2, 
           rate = s.lambda + 0.5 * tau * sum(beta^2))
  }
  r <- ncol(X)
  n <- length(y)
  Betaj <- matrix(NA, r, ite)
  Tau <- rep(NA, ite)
  Lambda <- rep(NA, ite)
  Betaj[, 1] <- betaj.0
  Tau[1] <- tau.0
  Lambda[1] <- lambda.0
  ##
  for(i in 2:ite) {
    ## update betaj:
    Betaj[,i] <- UpdateBetaj(r, XtX, Xty, Lambda[i-1], Tau[i-1])
    ## update tau:
    Tau[i] <- UpdateTau(y, X, r, n, r.tau, s.tau, Lambda[i-1], Betaj[,i])
    ##
    Lambda[i] <- UpdateLambda(r, r.lambda, s.lambda, Tau[i], Betaj[,i])
    cat("ite", i, "\n")
  }
  list(Betaj = Betaj, Tau = Tau, Lambda = Lambda)
}



GibbsAugmentedRidgeEigenFull <- function(y, 
                                         Z2,
                                         Z,
                                         Z1ty,
                                         P1,
                                         Pt,
                                         theta.hat,
                                         tau.hat,
                                         d, 
                                         ite,
                                         r.tau, 
                                         s.tau,
                                         r.lambda, 
                                         s.lambda,
                                         thetaj.0,   
                                         tau.0,
                                         lambda.0,
                                         w.0) {
  UpdateThetaj <- function(r, Z1ty, Z2, d, lambda, tau, w) {
    Z2tw <- as.vector(crossprod(Z2, w))
    Ztu <- Z1ty + Z2tw
    aux <- 1/(d + lambda)
    rnorm(r, aux * Ztu, sqrt(aux/tau))
  }
  UpdateTau <- function(y, Z, r, n, r.tau, s.tau, lambda, w, theta) {
    u <- c(y, w)
    rgamma(1, shape = r.tau + (r + n)/2, 
           rate = s.tau + 0.5 * (sum(theta^2 * lambda) + sum((u - Z %*% theta)^2)))
  }
  UpdateLambda <- function(r, r.lambda, s.lambda, tau, theta) {
    rgamma(1, shape = r.lambda + r/2, 
           rate = s.lambda + 0.5 * tau * sum(theta^2))
  }
  UpdateW <- function(Z2, n2, PtP1, theta.hat, tau.hat, theta, tau) {
    aux <- PtP1 %*% theta.hat
    aux <- theta + (tau.hat/tau) * aux
    aux <- as.vector(Z2 %*% aux)
    rnorm(n2, (tau/(tau + tau.hat)) * aux, sqrt(1/(tau + tau.hat)))
  }
  r <- ncol(Z2)
  n2 <- nrow(Z2)
  n <- nrow(Z)
  Thetaj <- matrix(NA, r, ite)
  Tau <- rep(NA, ite)
  Lambda <- rep(NA, ite)
  W <- matrix(NA, n2, ite)
  Thetaj[, 1] <- thetaj.0
  Tau[1] <- tau.0
  Lambda[1] <- lambda.0
  W[, 1] <- w.0
  PtP1 <- Pt %*% P1
  ##
  for(i in 2:ite) {
    ## update thetaj: r, Z1ty, Z2, d, lambda, tau, w
    Thetaj[,i] <- UpdateThetaj(r, Z1ty, Z2, d, Lambda[i-1], Tau[i-1], W[,i-1])
    ## update tau: y, Z, r, n, r.tau, s.tau, lambda, w, theta
    Tau[i] <- UpdateTau(y, Z, r, n, r.tau, s.tau, Lambda[i-1], W[,i-1], Thetaj[,i])
    ##
    Lambda[i] <- UpdateLambda(r, r.lambda, s.lambda, Tau[i], Thetaj[,i])
    ## update w: Z2, n2, P1, P2t, theta.hat, tau.hat, theta, tau
    W[,i] <- UpdateW(Z2, n2, PtP1, theta.hat, tau.hat, Thetaj[,i], Tau[i])
    cat("ite = ", i, "\n")
  }
  list(Thetaj = Thetaj, Tau = Tau, Lambda = Lambda, W = W)
}



GibbsAugmentedRidgeSvdFull <- function(y, 
                                       Z2,
                                       Z,
                                       Z1ty,
                                       V1,
                                       Vt,
                                       theta.hat,
                                       tau.hat,
                                       d, 
                                       ite,
                                       r.tau, 
                                       s.tau,
                                       r.lambda, 
                                       s.lambda,
                                       thetaj.0,   
                                       tau.0,
                                       lambda.0,
                                       w.0) {
  UpdateThetaj <- function(r, Z1ty, Z2, d, lambda, tau, w) {
    Z2tw <- as.vector(crossprod(Z2, w))
    Ztu <- Z1ty + Z2tw
    aux <- 1/(d^2 + lambda)
    rnorm(r, aux * Ztu, sqrt(aux/tau))
  }
  UpdateTau <- function(y, Z, r, n, r.tau, s.tau, lambda, w, theta) {
    u <- c(y, w)
    rgamma(1, shape = r.tau + (r + n)/2, 
           rate = s.tau + 0.5 * (sum(theta^2 * lambda) + sum((u - Z %*% theta)^2)))
  }
  UpdateLambda <- function(r, r.lambda, s.lambda, tau, theta) {
    rgamma(1, shape = r.lambda + r/2, 
           rate = s.lambda + 0.5 * tau * sum(theta^2))
  }
  UpdateW <- function(Z2, n2, VtV1, theta.hat, tau.hat, theta, tau) {
    aux <- VtV1 %*% theta.hat
    aux <- theta + (tau.hat/tau) * aux
    aux <- as.vector(Z2 %*% aux)
    rnorm(n2, (tau/(tau + tau.hat)) * aux, sqrt(1/(tau + tau.hat)))
  }
  r <- ncol(Z2)
  n2 <- nrow(Z2)
  n <- nrow(Z)
  Thetaj <- matrix(NA, r, ite)
  Tau <- rep(NA, ite)
  Lambda <- rep(NA, ite)
  W <- matrix(NA, n2, ite)
  Thetaj[, 1] <- thetaj.0
  Tau[1] <- tau.0
  Lambda[1] <- lambda.0
  W[, 1] <- w.0
  VtV1 <- Vt %*% V1
  ##
  for(i in 2:ite) {
    ## update thetaj: r, Z1ty, Z2, d, lambda, tau, w
    Thetaj[,i] <- UpdateThetaj(r, Z1ty, Z2, d, Lambda[i-1], Tau[i-1], W[,i-1])
    ## update tau: y, Z, r, n, r.tau, s.tau, lambda, w, theta
    Tau[i] <- UpdateTau(y, Z, r, n, r.tau, s.tau, Lambda[i-1], W[,i-1], Thetaj[,i])
    ##
    Lambda[i] <- UpdateLambda(r, r.lambda, s.lambda, Tau[i], Thetaj[,i])
    ## update w: Z2, n2, V1, Vt, theta.hat, tau.hat, theta, tau
    W[,i] <- UpdateW(Z2, n2, VtV1, theta.hat, tau.hat, Thetaj[,i], Tau[i])
    cat("ite = ", i, "\n")
  }
  list(Thetaj = Thetaj, Tau = Tau, Lambda = Lambda, W = W)
}


#############################################################
## functions based on the posterior predictive distribution 
#############################################################


GetVs <- function(Z1, Z2, lambda) {
  n2 <- nrow(Z2)
  vs <- rep(NA, n2)
  for (i in 1:n2) {
    cat("vi = ", i, "\n")
    z2i <- Z2[i,, drop =F]
    Zstar <- rbind(z2i, Z1)
    fsvd <- fast.svd(Zstar)
    dstar <- fsvd$d
    Vstar <- fsvd$v
    aux <- t(Vstar) %*% t(z2i)
    aux <- aux/(dstar^2 + lambda)
    aux <- Vstar %*% aux
    aux <- z2i %*% aux
    vs[i] <- 1 - as.vector(aux) 
  }
  vs
}


GibbsAugmentedRidgeEigenPPD <- function(y, 
                                        Z2,
                                        Z,
                                        Z1ty,
                                        P1,
                                        Pt,
                                        theta.hat,
                                        s.hat,
                                        vis,
                                        d, 
                                        ite,
                                        a.tau, 
                                        b.tau,
                                        lambda, 
                                        thetaj.0,   
                                        tau.0,
                                        w.0,
                                        etai.0) {
  UpdateThetaj <- function(p, Z1ty, Z2, d, lambda, tau, w) {
    Z2tw <- as.vector(crossprod(Z2, w))
    Ztu <- Z1ty + Z2tw
    aux <- 1/(d + lambda)
    rnorm(p, aux * Ztu, sqrt(aux/tau))
  }
  UpdateTau <- function(y, Z, p, n, a.tau, b.tau, lambda, w, theta) {
    u <- c(y, w)
    rgamma(1, shape = a.tau + (p + n)/2, 
           rate = b.tau + 0.5 * (lambda * sum(theta^2) + sum((u - Z %*% theta)^2)))
  }
  UpdateEtai <- function(a.tau, n1, s.hat, vis, w, Z2PtP1ThetaHat) {
    rgamma(n2, shape = 0.5 * (2 * a.tau + n1 + 1), 
           rate = 0.5 * (2 * a.tau + n1 + (vis/s.hat) * (w - Z2PtP1ThetaHat)^2))
  }
  UpdateW <- function(Z2, n2, Z2PtP1ThetaHat, vis, theta.hat, s.hat, theta, tau, lambda, 
                      eta) {
    M <- vis * eta
    aux1 <- 1/(M + tau * s.hat)
    aux2 <- as.vector(Z2 %*% theta)
    aux2 <- tau * s.hat * aux2
    aux3 <- M * Z2PtP1ThetaHat
    aux3 <- aux3 + aux2
    rnorm(n2, aux1 * aux3, sqrt(s.hat * aux1))
  }
  p <- ncol(Z2)
  n2 <- nrow(Z2)
  n <- nrow(Z)
  Thetaj <- matrix(NA, p, ite)
  Tau <- rep(NA, ite)
  W <- matrix(NA, n2, ite)
  Etai <- matrix(NA, n2, ite)
  Thetaj[, 1] <- thetaj.0
  Tau[1] <- tau.0
  W[, 1] <- w.0
  Etai[, 1] <- etai.0
  Z2PtP1ThetaHat <- P1 %*% theta.hat
  Z2PtP1ThetaHat <- Pt %*% Z2PtP1ThetaHat
  Z2PtP1ThetaHat <- Z2 %*% Z2PtP1ThetaHat
  Z2PtP1ThetaHat <- as.vector(Z2PtP1ThetaHat)
  ##
  for(i in 2:ite) {
    ## update thetaj: p, Z1ty, Z2, d, lambda, tau, w
    Thetaj[,i] <- UpdateThetaj(p, Z1ty, Z2, d, lambda, Tau[i-1], W[,i-1])
    ## update tau: y, Z, p, n, a.tau, b.tau, lambda, w, theta
    Tau[i] <- UpdateTau(y, Z, p, n, a.tau, b.tau, lambda, W[,i-1], Thetaj[,i])
    ## update w: Z2, n2, Z2PtP1ThetaHat, vis, theta.hat, s.hat, theta, tau, lambda, eta
    W[,i] <- UpdateW(Z2, n2, Z2PtP1ThetaHat, vis, theta.hat, s.hat, Thetaj[,i], Tau[i], 
                     lambda, Etai[,i-1])
    ## update etaj: a.tau, n1, s.hat, vis, w, Z2PtP1ThetaHat
    Etai[,i] <- UpdateEtai(a.tau, n1, s.hat, vis, W[,i], Z2PtP1ThetaHat)
    cat("ite = ", i, "\n")
  }
  list(Thetaj = Thetaj, Tau = Tau, W = W, Etai = Etai)
}




GibbsAugmentedRidgeEigenPPD2 <- function(y, 
                                         Z2,
                                         Z,
                                         Z1ty,
                                         P1,
                                         Pt,
                                         theta.hat,
                                         s.hat,
                                         vis,
                                         d, 
                                         ite,
                                         a.tau, 
                                         b.tau,
                                         lambda, 
                                         thetaj.0,   
                                         tau.0,
                                         w.0,
                                         etai.0,
                                         thin.window = 10) {
  UpdateThetaj <- function(p, Z1ty, Z2, d, lambda, tau, w) {
    Z2tw <- as.vector(crossprod(Z2, w))
    Ztu <- Z1ty + Z2tw
    aux <- 1/(d + lambda)
    rnorm(p, aux * Ztu, sqrt(aux/tau))
  }
  UpdateTau <- function(y, Z, p, n, a.tau, b.tau, lambda, w, theta) {
    u <- c(y, w)
    rgamma(1, shape = a.tau + (p + n)/2, 
           rate = b.tau + 0.5 * (lambda * sum(theta^2) + sum((u - Z %*% theta)^2)))
  }
  UpdateEtai <- function(a.tau, n1, s.hat, vis, w, Z2PtP1ThetaHat) {
    rgamma(n2, shape = 0.5 * (2 * a.tau + n1 + 1), 
           rate = 0.5 * (2 * a.tau + n1 + (vis/s.hat) * (w - Z2PtP1ThetaHat)^2))
  }
  UpdateW <- function(Z2, n2, Z2PtP1ThetaHat, vis, theta.hat, s.hat, theta, tau, lambda, 
                      eta) {
    M <- vis * eta
    aux1 <- 1/(M + tau * s.hat)
    aux2 <- as.vector(Z2 %*% theta)
    aux2 <- tau * s.hat * aux2
    aux3 <- M * Z2PtP1ThetaHat
    aux3 <- aux3 + aux2
    rnorm(n2, aux1 * aux3, sqrt(s.hat * aux1))
  }
  ThinCycle <- function(y, 
                        Z2,
                        Z,
                        Z1ty,
                        theta.hat,
                        s.hat,
                        vis,
                        d, 
                        a.tau, 
                        b.tau,
                        lambda,
                        Z2PtP1ThetaHat,
                        ite = thin.window,
                        Thetaj,
                        Tau,
                        W,
                        Etai) {
    for (i in 1:ite) {
      Thetaj <- UpdateThetaj(p, Z1ty, Z2, d, lambda, Tau, W)
      Tau <- UpdateTau(y, Z, p, n, a.tau, b.tau, lambda, W, Thetaj)
      W <- UpdateW(Z2, n2, Z2PtP1ThetaHat, vis, theta.hat, s.hat, Thetaj, Tau, 
                   lambda, Etai)
      Etai <- UpdateEtai(a.tau, n1, s.hat, vis, W, Z2PtP1ThetaHat)  
    }
    list(Thetaj = Thetaj, Tau = Tau, W = W, Etai = Etai)
  }
  p <- ncol(Z2)
  n2 <- nrow(Z2)
  n <- nrow(Z)
  Thetaj <- matrix(NA, p, ite)
  Tau <- rep(NA, ite)
  W <- matrix(NA, n2, ite)
  Etai <- matrix(NA, n2, ite)
  Thetaj[, 1] <- thetaj.0
  Tau[1] <- tau.0
  W[, 1] <- w.0
  Etai[, 1] <- etai.0
  Z2PtP1ThetaHat <- P1 %*% theta.hat
  Z2PtP1ThetaHat <- Pt %*% Z2PtP1ThetaHat
  Z2PtP1ThetaHat <- Z2 %*% Z2PtP1ThetaHat
  Z2PtP1ThetaHat <- as.vector(Z2PtP1ThetaHat)
  ##
  for(i in 2:ite) {
    aux <-   ThinCycle(y, Z2, Z, Z1ty, theta.hat, s.hat, vis, d, a.tau, b.tau,
                       lambda, Z2PtP1ThetaHat, ite = thin.window, 
                       Thetaj[,i-1], Tau[i-1], W[,i-1], Etai[,i-1])
    Thetaj[,i] <- aux$Thetaj
    Tau[i] <- aux$Tau
    W[,i] <- aux$W
    Etai[,i] <- aux$Etai
    cat("ite = ", i, "\n")
  }
  list(Thetaj = Thetaj, Tau = Tau, W = W, Etai = Etai)
}





########################################################################
## miscellanea
########################################################################


GetMse <- function(X, Y, betah) {
  t(Y - X %*% betah) %*% (Y - X %*% betah)/nrow(X)
}



GetPostProbInterval <- function(Beta, prob, true.beta) {
  beta.mean <- apply(Beta, 1, mean)
  lower.bound <- apply(Beta, 1, quantile, (1 - prob)/2)
  upper.bound <- apply(Beta, 1, quantile, (1 - prob)/2 + prob)
  contain.zero <- (lower.bound < 0) & (upper.bound > 0)
  post.int <- data.frame(true.beta, beta.mean, lower.bound, upper.bound, contain.zero)
  aux.true.beta <- true.beta == 0
  correct.selections <- sum(contain.zero == aux.true.beta)
  aux <- order(abs(true.beta), decreasing = TRUE)
  list(posterior.interval = post.int[aux,], correct.selections = correct.selections, 
       index = aux)
}



PlotPostInterval <- function(x, cex.lab = 1, cex.axis = 1, ylim = NULL, main = "") {
  my.lim <- c(min(x$lower.bound) - 1, max(x$upper.bound) + 1)
  plot(x$beta.mean, ylab = expression(beta), 
       xlab = "sorted regression coefficient index", cex.lab = cex.lab, 
       cex.axis = cex.axis, type = "n", ylim = ylim, main = main)
  abline(h = 0, col = "black")
  for (i in 1:nrow(x)) {
    segments(x0 = i, x1 = i, y0 = x$lower.bound[i], y1 = x$upper.bound[i], col = "red")
    points(i, x$beta.mean[i], col = "black", cex = 1, pch = 19)
    points(i, x$true.beta[i], col = "blue", cex = 1, pch = 20)
  }  
}





