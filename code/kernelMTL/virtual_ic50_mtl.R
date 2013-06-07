
library(abind)
source('code/kernelMTL/helper.R')  
source('code/kernelMTL/bemkl_multitask_supervised_regression_variational_train.R')
source('code/kernelMTL/bemkl_multitask_supervised_regression_variational_test.R')

defaultParameters <- function(){
  parameters <- list()
  parameters$alpha_lambda <- 1
  parameters$beta_lambda <- 1
  parameters$alpha_upsilon <- 1
  parameters$beta_upsilon <- 1
  parameters$alpha_gamma <- 1
  parameters$beta_gamma <- 1
  parameters$alpha_omega <- 1
  parameters$beta_omega <- 1
  parameters$alpha_epsilon <- 1
  parameters$beta_epsilon <- 1
  parameters$iteration <- 50 #number of iterations, 200 iterations are enough usually
  parameters$seed <- 2012
  parameters$progress <- 0 #if 1, calculates the lower bound at each iteration
  return (parameters)
}

combineKernels <- function(...){
  input <- list(...)
  Ktrain <- lapply(input, function(x) x$train)
  Ktest <- lapply(input, function(x) x$test)
 
  Ktrain.final <- abind(Ktrain,along=3)
  Ktest.final <- lapply(1:length(Ktest[[1]]), function(i) { 
    tmp <- lapply(Ktest, function(x) { x[[i]]})
    abind(tmp,along=3)
  })
  
  return(list(train=Ktrain.final, test=Ktest.final))
}

factorKernel <- function(trainLabeledFactorList, trainUnlabeledFactorList, testFactorList){
  allLevels <- c()
  for(f in c(trainLabeledFactorList, trainLabeledFactorList, testFactorList)){
    allLevels <- unique(union(allLevels, f))
  }
  buildFactorMatrix <- function(f){
    sapply(allLevels, function(x){ f %in% x})
  }
  
  trainLabeledLMatrix <- do.call("rbind",lapply(trainLabeledFactorList, buildFactorMatrix))
  trainUnlabeledLMatrix <- do.call("rbind",lapply(trainUnlabeledFactorList, buildFactorMatrix))
  testLMatrixList <- lapply(testFactorList, buildFactorMatrix)
  
  Ktrain <- logicalORKernel(rbind(trainLabeledLMatrix, trainUnlabeledLMatrix), trainLabeledLMatrix)
  Ktests <- lapply(testLMatrixList, function(L){
    logicalORKernel(rbind(trainLabeledLMatrix, trainUnlabeledLMatrix), L)
  })
  
  return (list(train=Ktrain,test=Ktests))
}



rbfGeneExprKernel <- function(trainLabeledExprList, trainUnlabeledExprList, testExprList, genesOfInterest=NULL, gaussian_factors=c(16, 8, 4, 2, 1, .5, .25, .125)){
  data <- c(trainLabeledExprList,trainUnlabeledExprList,testExprList)
  common <- rownames(data[[1]])
  for(expr in data[-1]){
    common <- intersect(common, rownames(expr))
  }
  
  if(!is.null(genesOfInterest)){
    common <- common[common %in% genesOfInterest]
  }

  trainLabeledExprMatrix <- do.call("rbind",lapply(trainLabeledExprList, function(expr){
    scale(t(expr[match(common, rownames(expr)),]))
  }))
  trainUnlabeledExprMatrix <- do.call("rbind",lapply(trainUnlabeledExprList, function(expr){
    scale(t(expr[match(common, rownames(expr)),]))
  }))
  testExprList <- lapply(testExprList, function(expr){
    scale(t(expr[match(common, rownames(expr)),]))
  })
  
 
  gaussian_sigs <- sqrt((nrow(trainLabeledExprMatrix) + nrow(trainUnlabeledExprMatrix)) * gaussian_factors)
  # Nbasis x Ntrain X P kernel matrices (Nbasis and Ntrain can be different for each task)
  Ktrain <- gaussianKernel(rbind(trainLabeledExprMatrix, trainUnlabeledExprMatrix), trainLabeledExprMatrix, gaussian_sigs)
  Ktests <- lapply(testExprList, function(expr){
                  gaussianKernel(rbind(trainLabeledExprMatrix, trainUnlabeledExprMatrix), expr, gaussian_sigs)
  })

  return (list(train=Ktrain,test=Ktests))
}

virtual_ic50_mtl_v2 <- function(trainKernel, y, testKernels, parameters){
  
  stopifnot(length(y)==ncol(trainKernel))
  
  if(any(is.na(y))){
    cat("Found NAs in y.  Removing...")
    trainKernel <- trainKernel[,!is.na(y),]
    y <- y[!is.na(y)]
  }
  
  Kmtrainlist <- list(trainKernel)
  ytrainlist <- list(y)
  
  state <- bemkl_multitask_supervised_regression_variational_train(Kmtrainlist, ytrainlist, parameters)
  
  predictions <- lapply(testKernels, function(Ktest){
    bemkl_multitask_supervised_regression_variational_test(list(Ktest), state)
  })
  y_hats <- lapply(predictions, function(x){x$f[[1]]$mean})
  
  return(list(y_hats=y_hats,state=state,predictions=predictions))
}


virtual_ic50_mtl <- function(cellLineEset, drugName, tcga.dat, seed=2013, reverseDrug=FALSE, tcga.disease=NULL, other.tests=list()){
  common <- intersect(rownames(tcga.dat$geneExpr), featureNames(cellLineEset))
  for(i in 1:length(other.tests)){
    common <- intersect(common, rownames(other.tests[[i]]))
  }
  
  tcga.m <- tcga.dat$geneExpr[match(common, rownames(tcga.dat$geneExpr)),]
  CL.m <- cellLineEset[match(common,featureNames(cellLineEset)),]
  for(i in 1:length(other.tests)){
    tmp <- other.tests[[i]]
    other.tests[[i]] <- scale(t(tmp[match(common, rownames(tmp)),]))
  }
  
  
  set.seed(seed)
  drug_vec <- pData(CL.m)[,drugName]
  if(reverseDrug){
    drug_vec <- -drug_vec
  }
  mask <- !is.na(drug_vec)
  cl.X <- scale(t(exprs(CL.m[,mask])))
  tcga.X <- scale(t(tcga.m))
  y = drug_vec[mask]
  
  gaussian_factors <- c(16, 8, 4, 2, 1, .5, .25, .125)
  gaussian_sigs <- sqrt(ncol(tcga.X) * gaussian_factors)
  # Nbasis x Ntrain X P kernel matrices (Nbasis and Ntrain can be different for each task)
  Ktrain <- gaussianKernel(rbind(cl.X,tcga.X), cl.X, gaussian_sigs)
  Ktest <- gaussianKernel(rbind(cl.X,tcga.X), tcga.X, gaussian_sigs)
  browser()
  
  # build kernels for other datasets of interest
  Kothertests <- lapply(other.tests, function(X){
    gaussianKernel(rbind(cl.X,tcga.X), X, gaussian_sigs)
  })
  
  if(!is.null(tcga.disease)){
    diseaseLbls <- gsub(".*?_(.*)","\\1",rownames(cl.X))
    diseases <- sort(unique(diseaseLbls))
    if(!(tcga.disease %in% diseases)){
      stop(tcga.disease," is not a valid tumor type.")
    }
    D.cl <- matrix(NA, nrow=length(diseaseLbls),ncol=length(diseases))
    for(i in 1:(length(diseases))){
      D.cl[,i] <- as.numeric(diseaseLbls %in% diseases[i])
    }
    
    D.tcga <- matrix(0, nrow=nrow(tcga.X), ncol=ncol(D.cl))
    D.tcga[, which(diseases %in% tcga.disease)] <- 1
    Ktrain.disease <- logicalORKernel(rbind(D.cl, D.tcga), D.cl)
    Ktest.disease <- logicalORKernel(rbind(D.cl, D.tcga), D.tcga)
    Ktrain <- abind(Ktrain, Ktrain.disease,along=3)
    Ktest <- abind(Ktest, Ktest.disease,along=3)
    
    #TODO add 'otherTest' kernel
    for(i in 1:length(other.tests)){
      # assume target disease
      Ktmp <- Kothertests[[i]]
      X <- other.tests[[i]]
      D <- matrix(0, nrow=nrow(X), ncol=ncol(D.cl))
      D[, which(diseases %in% tcga.disease)] <- 1
      Ktest.disease <- logicalORKernel(rbind(D.cl, D.tcga), D)
      Kothertests[[i]] <- abind(Ktmp, Ktest.disease,along=3)
    }
  }
  Kmtrainlist <- list(Ktrain)
  ytrainlist <- list(y)
  Kmtestlist <- list(Ktest)

  parameters <- list()
  parameters$alpha_lambda <- 1
  parameters$beta_lambda <- 1
  parameters$alpha_upsilon <- 1
  parameters$beta_upsilon <- 1
  parameters$alpha_gamma <- 1
  parameters$beta_gamma <- 1
  parameters$alpha_omega <- 1
  parameters$beta_omega <- 1
  parameters$alpha_epsilon <- 1
  parameters$beta_epsilon <- 1
  parameters$iteration <- 200 #number of iterations, 200 iterations are enough usually
  parameters$seed <- seed #random number seed for different replications
  parameters$progress <- 0 #if 1, calculates the lower bound at each iteration
  
  
  state <- bemkl_multitask_supervised_regression_variational_train(Kmtrainlist, ytrainlist, parameters)
  prediction <- bemkl_multitask_supervised_regression_variational_test(Kmtestlist, state)
  y_hat <- prediction$f[[1]]$mean
  names(y_hat) <- colnames(tcga.m)
  
  yhats_other <- lapply(1:length(Kothertests), function(i){
    pred <- bemkl_multitask_supervised_regression_variational_test(list(Kothertests[[i]]), state)
    yhat <- pred$f[[1]]$mean
    names(yhat) <- rownames(other.tests[[i]])
    yhat
  })
  
  return(list(y_hat=prediction$f[[1]]$mean,state=state,prediction=prediction, yhats_other=yhats_other))
}

##########
# helper functions

# n-sample by p-gene matrix
gaussianKernel <- function(X1, X2, sigs){
  require(kernlab)
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  n3 <- length(sigs)
  #X1 <- t(X1)
  #X2 <- t(X2)
  K <- array(NA,dim=c(n1,n2,n3))
  for(k in 1:n3){
    rbf <- rbfdot(sigma = 1/sigs[k]^2)
    K[,,k]  <- kernelMatrix(rbf, X1, X2)
  }
  return (K)
}

# n-sample by p-factor matrix
logicalORKernel <- function(X1, X2){
  return(X1 %*% t(X2))
}

jaccardKernel <- function(X1, X2){
  require(fpc)
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  K <- matrix(NA,nrow=n1,ncol=n2)
  for(i in 1:n1){
    for(j in 1:n2){
      K[i,j] <- clujaccard(X1[i,], X2[i,])
    }
  }
  return (K)
}