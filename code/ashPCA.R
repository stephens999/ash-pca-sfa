#input: Y an n by p data matrix and niter, the number of shrinkage iterations to perform

#output: the PCs (n-vectors) and the shrunken loadings (p vectors) and the residuals (n by p)
ashPCA = function(Y,niter){
  Y.svd = svd(Y)
  betahat = Y.svd$d * t(Y.svd$v)  #kth row of betahat are the coeff of the kth pc
  betahat.shrunk = betahat
  sebetahat = 1
  n = nrow(Y)
  for(i in 1:niter){
      for(k in 1:n){
          betahat.mean = mean(betahat[k,])
          betahat.ash = ash(betahat[k,]-betahat.mean,sebetahat,prior="uniform")
          betahat.shrunk[k,] = betahat.ash$PosteriorMean+betahat.mean
      }
      resid = Y - Y.svd$u %*% betahat.shrunk
      sebetahat = sqrt(colMeans(resid*resid))
  }
  return(list(betahat=betahat,sebetahat=sebetahat,Ycorr=resid))
}
    
