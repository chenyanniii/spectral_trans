
# set.seed(123)
# X <- matrix(rnorm(5*5), 5, 5)
# Y <- matrix(rnorm(5*3), 5, 3)

# lambda1=0.1
# lambda2=0.1

# M = nrow(X) 
# folds = split(1:M, 1:M)

#a = Mfold(X,Y,lambda1,lambda2,folds)




Mfold = function(X, Y, lam1, lam2, folds){
  results <- foreach(m = 1:length(folds), .combine = rbind, .packages = 'mixOmics') %dopar% {

  # Perform the cross-validation for tuning    
  omit=folds[[m]]

  result = rcc(X[-omit, , drop = FALSE], Y[-omit, , drop = FALSE], 
               lambda1 = lam1, lambda2 = lam2)
  X[omit, ][is.na(X[omit, ])] = 0
  Y[omit, ][is.na(Y[omit, ])] = 0
  xscore = X[omit, , drop=FALSE] %*% result$loadings$X[,1]
  yscore = Y[omit, , drop=FALSE] %*% result$loadings$Y[,1]
  xy_scores <- data.frame(
    m = m,
    lambda1 = lam1,
    lambda2 = lam2,
    x_score = as.vector(xscore),
    y_score = as.vector(yscore)
  )
  return(xy_scores)
  }
  cv_score = data.frame(lambda1 = lam1, lambda2 = lam2,
                        cv_score = cor(results$x_score, results$y_score, use = "pairwise"))
}

