



# Sample data (replace with your actual data)
#set.seed(123)
#X <- matrix(rnorm(5*5), 5, 5)
#Y <- matrix(rnorm(5*3), 5, 3)

#grid1=seq(0.1, 1, length.out = 3)
#grid2=seq(0.1, 1, length.out = 3)

# a = tune_rcc_parallel(X, Y, grid1, grid2)






# Parallelized tuning using foreach

## function for parallelly cross-validation per set of lambda parameters
source("scripts/func_Mfold-parallel.R")

## use the methods validation == "loo"
tune_rcc_parallel = function(X, Y, grid1, grid2){
  
  # prep for paralleling tuning using foreach{}
  M = nrow(X) 
  folds = split(1:M, 1:M)
  grid <- expand.grid(lambda1 = grid1, lambda2 = grid2)
  
  ## summary cross-validation of lambdas in parallel
  cv_score_summary = foreach(i = 1:nrow(grid), .combine = rbind, .packages = 'mixOmics') %do% {
    lam1 <- grid$lambda1[i]
    lam2 <- grid$lambda2[i]
    results = Mfold(X,Y, lam1, lam2, folds)
  }
  
  
  ## selection of best cross-validated parameters
  opt = cv_score_summary[cv_score_summary[, 3] == max(cv_score_summary[, 3]), ]
  
  out = list(opt.lambda1 = opt[[1]], opt.lambda2 = opt[[2]], 
             opt.cv.score = opt[[3]], grid1 = grid1, grid2 = grid2, cv_score = cv_score_summary)
  
  out$call = match.call()
  
  class(out) = "tune.rcc"
  
  return(invisible(out))
}


























