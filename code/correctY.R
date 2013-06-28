correctY = function(Y, X, shrink = TRUE) {
  K = ncol(X)
  g = ncol(Y)
  #       X = scale(X, scale=FALSE) #remove column means Y = scale(Y,scale=scaleY)
  #       #remove column means, and scale if scale=T
  SX2 = colSums(X * X)  #a K vector of sums of squares
  EY2 = colMeans(Y * Y)  #variance of Y
  betahat = SX2^{-1} * (t(X) %*% Y)  # a K by g matrix
  betahat.shrunk = matrix(0, nrow = K, ncol = g)
  sebetahat.list = c()
  betahat.centered.list = c()
  betahat.shrunk.centered.list = c()
  for (k in 1:K) {
    #       for (k in 1:2) {
    resid = Y - outer(X[, k], betahat[k, ])  # residuals, n by g
    sigmahat = sqrt(colMeans(resid * resid))  # g vector, estimated residual variances
    sebetahat = SX2[k]^{-1} * sigmahat
    sebetahat.list <- cbind(sebetahat.list, sebetahat)
    betahat.mean <- mean(betahat[k,])
    betahat.centered <- betahat[k,] - betahat.mean
    betahat.centered.list <- rbind(betahat.centered.list, betahat.centered)
    #                       betahat.shrunk[k, ] = ash(betahat[k, ], sebetahat, prior = "uniform")$PosteriorMean
    betahat.shrunk.centered.list <- rbind(betahat.shrunk.centered.list, ash(betahat.centered, sebetahat, prior = "uniform")$PosteriorMean)
    betahat.shrunk[k, ] = ash(betahat.centered, sebetahat, prior = "uniform")$PosteriorMean + betahat.mean
  }
  Ycorr = Y - X %*% betahat.shrunk  # Y corrected for the effects of X
  
  ###Edits###
  betahat.shrunk.rowSums <- rowSums(betahat.shrunk)
  betahat.rowSums <- rowSums(betahat)
  sebetahat.rowMeans <- rowMeans(sebetahat.list)
  
  ###########
  
  return(list(Ycorr = Ycorr, betahat = betahat, betahat.centered = betahat.centered.list, betahat.shrunk = betahat.shrunk, betahat.shrunk.centered = betahat.shrunk.centered.list, sebetahat = sebetahat.list, sebetahat.rowMeans = sebetahat.rowMeans, shrunk.rowSums = betahat.shrunk.rowSums, betahat.rowSums = betahat.rowSums))
}