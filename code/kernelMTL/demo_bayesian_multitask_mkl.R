# Mehmet Gonen (mehmet.gonen@aalto.fi)
# Helsinki Institute for Information Technology HIIT
# Department of Information and Computer Science
# Aalto University School of Science
require('R.matlab')

demo_bayesian_multitask_mkl <- function() {
  source('./helper.R')  
  source('./bemkl_multitask_supervised_regression_variational_train.R')
  source('./bemkl_multitask_supervised_regression_variational_test.R')
  parameters <- list()
  parameters$alpha_lambda <- 1e-10
  parameters$beta_lambda <- 1e-10
  parameters$alpha_upsilon <- 1
  parameters$beta_upsilon <- 1
  parameters$alpha_gamma <- 1e-10
  parameters$beta_gamma <- 1e-10
  parameters$alpha_omega <- 1e-10
  parameters$beta_omega <- 1e-10
  parameters$alpha_epsilon <- 1
  parameters$beta_epsilon <- 1
  parameters$iteration <- 200 #number of iterations, 200 iterations are enough usually
  parameters$seed <- 1606 #random number seed for different replications
  parameters$progress <- 0 #if 1, calculates the lower bound at each iteration

  Ktrain <- readMat('Ktrain.mat')
  ytrain <- readMat('ytrain.mat')
  Ktest <- readMat('Ktest.mat')

  #the number of tasks is assumed to be T = 4 and the number of kernels is assumed to be P = 3
  T <- 4
  P <- 3
  Kmtrainlist <- vector("list", T) #Nbasis x Ntrain X P kernel matrices (Nbasis and Ntrain can be different for each task)
  ytrainlist <- vector("list", T) #Ntrain x T target outputs (Ntrain can be different for each task)
  Kmtrainlist[[1]] <- Ktrain$Ktrain
  ytrainlist[[1]] <- ytrain$ytrain
  Kmtrainlist[[2]] <- Ktrain$Ktrain
  ytrainlist[[2]] <- ytrain$ytrain
  Kmtrainlist[[3]] <- Ktrain$Ktrain
  ytrainlist[[3]] <- ytrain$ytrain
  Kmtrainlist[[4]] <- Ktrain$Ktrain
  ytrainlist[[4]] <- ytrain$ytrain
  
  Kmtestlist <- vector("list", T) #Nbasis x Ntest X P kernel matrices (Nbasis and Ntest can be different for each task)
  Kmtestlist[[1]] <- Ktest$Ktest
  Kmtestlist[[2]] <- Ktest$Ktest
  Kmtestlist[[3]] <- Ktest$Ktest
  Kmtestlist[[4]] <- Ktest$Ktest

  state <- bemkl_multitask_supervised_regression_variational_train(Kmtrainlist, ytrainlist, parameters)
  prediction <- bemkl_multitask_supervised_regression_variational_test(Kmtestlist, state)
  
  par(mfrow = c(2, 2))
  plot(prediction$f[[1]]$mean)
  plot(prediction$f[[2]]$mean)
  plot(prediction$f[[3]]$mean)
  plot(prediction$f[[4]]$mean)
}