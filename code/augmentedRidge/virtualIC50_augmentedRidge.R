source("./aridge/aug_ridge_code_for_justin.R")

virtual_ic50_augmentedRidge <- function(cellLineEset, drugName, tcga.dat, seed=2013, reverseDrug=FALSE,burnin=3000,svd.mode=TRUE){
  
  # find common genes beween celllines and tCGA 
  common <- intersect(rownames(tcga.dat$geneExpr), featureNames(cellLineEset))
  tcgaExpr.m <- tcga.dat$geneExpr[match(common, rownames(tcga.dat$geneExpr)),]
  CL.m <- cellLineEset[match(common,featureNames(cellLineEset)),]
  
  set.seed(seed)
  drug_vec <- pData(CL.m)[,drugName]
  if(reverseDrug){
    drug_vec <- -drug_vec
  }
  mask <- !is.na(drug_vec)
  X <- scale(t(exprs(CL.m[,mask])))
  y = drug_vec[mask] - mean(drug_vec[mask])
  
  ## apply models on new data
 
  tcgaExpr.scaled <- scale(t(tcgaExpr.m))
  n1 <- nrow(X)
  n2 <- nrow(tcgaExpr.scaled)
  
  lambda_opt <- best.lambda(t(X), y)
  b_hat <- ridgesvd(t(X), y, lambda_opt)
  s1 <- (y - t(X %*% b_hat) %*% y) / n1
  Vs <- GetVs
}

mse <- function(y1, y2){
  mean((y1-y2)^2)
}

#########
# find optimal lambda that minimizes mse in cross validation
best.lambda <- function(X, y, lambda_grid=.1^(seq(-5,5,by=.2)), cv=10){
  mses <- sapply(lambda_grid, function(lambda){
    mse(y, cv.ridge(X, y, lambda, cv))
  })
  lambda_grid[order(mses)[1]]
}

cv.ridge <- function(X,y,lambda,cv=10){
  n <- ncol(X)
  cv.groups <- rep(1:cv, n / cv + 1)[sample(n)]
  yhat <- rep(NA,n)
  for(i in 1:cv){
    mask <- i %in% cv.groups
    bhat <- ridgesvd(X[, !mask], y[!mask], lambda)
    y_hat[mask] <- t(X[, mask]) %*% bhat
  }
  yhat
}

##############
# Ridge regresion using SVD trick
# X = row normlized p-by-n dimensional matrix
# y = mean centered, n-dimensional vector
# TODO optimized to take multiple lambdas
ridgesvd <- function(X,y,lambda=1){
  p <- nrow(X)
  
  pc <- fast.svd(X)
  div <- diag(pc$d^2 + lambda)
  rhs <- t(pc$v) %*% y
  
  Bhats <- pc$u %*% drop(pc$d * rhs) / div
  Bhats
}

GetVs <- function(X1, X2, lambda) {
  n2 <- nrow(X2)
  vs <- rep(NA, n2)
  for (i in 1:n2) {
    cat("vi = ", i, "\n")
    z2i <- X2[i,, drop =F]
    Zstar <- rbind(z2i, X1)
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


