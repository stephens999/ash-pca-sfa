# Background

Currently when correcting for PCs $(v_1,\dots,v_K)$, and
testing a SNP $G$ we
fit $$Y_g = a_{1g} v_1 + .... + a_{Kg} v_K + G \beta + e$$
and estimate $a_{1g},\dots, a_{Kg}$ separately for each gene $g$
while testing $H_0:\beta=0$. 

The LMM instead assumes a prior for $a_{1g},\dots, a_{Kg}$ based
on the eigenvalues. Specifically, $a_{kg} \sim N(0, \sigma_g^2 \lambda_k)$ where $\sigma_g^2$ is estimated as part of the LMM fitting. 

Here's the idea: instead, how about fitting the same model as the LMM, but estimating the prior on $a_{kg}$ from the data, combining information
across genes, instead of fixing to the eigenvalues.
The idea is that some PCs will correlated strongly with multiple phenotype, and others (eg higher $k$) correlate with none, and we can learn that from the data.

Here's a suggested first thing to try.

1. estimate the effect $a_{kg}$ of PC $k$ on gene $g$, one at a time, by doing simple linear regression $Y_g = a_{kg} v_k + e$.
This gives $\hat{a}_{kg}$ and a standard error, we'll denote $s_{kg}$. (Note that because we are doing the PCs one at a time, rather than 
all at once, the residual variance will tend to be overestimated, so $s_{kg}$ will tend to overestimate the standard error).

2. Run the $g$ vectors $ahat_k, s_k$ through ASH to shrink them.
You'll have to download ASH from my github repo ash, and look
at the readme to see what it is doing.
But bottom line is that it will provide "shrunk" estimates for $ahat_{kg}$.

3. Compute the residuals $r_g = Y_g - \sum_k a_{kg} v_k$.
(3b. later - we could possibly use $r_g$ to estimate $\sigma_g$, the residual variance for $Y_g$,
and then this could be used to get a revised estimate of the standard error $s_{kg}$ )

4. Test the residuals for association with genotype.

Here is a start at writing R code for this:
First read in ash functions

```r
setwd("~/Documents/git/ash-pca-sfa/code")
ash.repodir = scan("../.ash.repodir.txt", what = character())  #say where on your computer you have the ash repo
source(file.path(ash.repodir, "/Rcode/ash.R"))
```





```r
source("correctY.R")
```


Here we'll simulate some data based on latent factors $Z$

```r
set.seed(124)
n = 100
G = 1000
K = 14
lambda = (0.9)^(1:K)
L = matrix(rnorm(K * n), nrow = n) %*% diag(lambda)  #n by k matrix of loadings
FF = matrix(rnorm(K * G), nrow = K)  #K by g matrix of factors
E = matrix(rnorm(n * G), nrow = n)
Y = L %*% FF + E
# Y = scale(Y) #centers and scales each gene

Y.svd = svd(Y)
plot(Y.svd$d, ylab = "singular values")  #singular values
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 

```r
X = Y.svd$u  # columns of X are PCs (n n-vectors)
Y.corr = correctY(Y, X)
```


Note that the early eigenvectors are shrunk only a little, whereas
the 15 and higher are shrunk very strongly

```r
plot(Y.corr$betahat[1, ], Y.corr$betahat.shrunk[1, ], ylim = c(-4, 4))
abline(a = 0, b = 1)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-41.png) 

```r
plot(Y.corr$betahat[7, ], Y.corr$betahat.shrunk[7, ], ylim = c(-4, 4))
abline(a = 0, b = 1)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-42.png) 

```r
plot(Y.corr$betahat[14, ], Y.corr$betahat.shrunk[14, ], ylim = c(-4, 4))
abline(a = 0, b = 1)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-43.png) 

```r
plot(Y.corr$betahat[15, ], Y.corr$betahat.shrunk[15, ], ylim = c(-4, 4))
abline(a = 0, b = 1)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-44.png) 


Ideally the corrected Y should look something like E, so check this:

```r
plot(E, Y.corr$Ycorr)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 


Let's compare with just correcting on the basis of the first k PCs.

```r
# this function returns estimate of $Y$ from truncating the svd, which is
# like correcting for the first k PCs
trunc.svd = function(Y.svd, k) {
    Y.svd$u[, 1:k, drop = F] %*% diag(Y.svd$d[1:k, drop = F], nrow = length(Y.svd$d[1:k, 
        drop = F])) %*% t(Y.svd$v[, 1:k, drop = F])
}
plot(E, Y - trunc.svd(Y.svd, 14))
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-61.png) 

```r
plot(Y.corr$Ycorr, Y - trunc.svd(Y.svd, 14))
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-62.png) 

So it is somewhat similar.

We could use correlation between E and the corrected Y as a measure of 
how well we are doing. (Maybe would be better to do this column by column and look at the average correlation, rather than pooling everything together.)

```r
cc = rep(0, 30)
for (i in 1:30) {
    cc[i] = cor(as.vector(Y - trunc.svd(Y.svd, i)), as.vector(E))
}
plot(cc)
abline(h = cor(as.vector(Y.corr$Ycorr), as.vector(E)))
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 

So we can see that by this measure correcting for 14 PCs does best
(which is kind of expected given the way the data were simulated).
The shrinkage does worse than this, but better than only
doing 12 PCs, and better than doing 20 PCs.

## Cross validation

One could worry about overfitting: the
PCs are selected from the data, so the effects of each PC will be systematically overestimated on average, and this
could cause overfitting/undershrinkage.

To asses whether this is a problem we try cross validation, correcting genes 1 to 500 with PCs learned from other genes, and vice versa.



```r
Y1 = Y[, 1:500]
Y2 = Y[, 501:1000]
Y1.svd = svd(Y1)
Y2.svd = svd(Y2)

X1 = Y1.svd$u  # columns of X are PCs (n n-vectors)
X2 = Y2.svd$u
Y1.corr = correctY(Y1, X2)
Y2.corr = correctY(Y2, X1)
Ycv.corr = cbind(Y1.corr$Ycorr, Y2.corr$Ycorr)
plot(E, Ycv.corr)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 

```r
cor(as.vector(Ycv.corr), as.vector(E))
```

```
## [1] 0.8636
```


### The Stranger Data


```r
setwd("~/Documents/git/ash-pca-sfa/")
Data1 <- read.table("data/beforeCorrection_onlyExpGenes.txt.pedOrder.Transposed.allPopNorm", 
    header = FALSE)
xx = matrix(as.numeric(as.matrix(Data1[-1, -(1:2)])), nrow = 210)
```



## An alternative idea

Another idea would be to ``correct" each gene by all the other genes, rather
than by the PCs. This might avoid the overfitting problem, because it avoids first choosing the PCs by looking at the data.
To investigate we'll try correcting the first J genes using all the others.

```r
J = 10
Y1.corr = correctY(Y[, 1:J], Y[, -(1:J)])
plot(Y1.corr$Ycorr, E[, 1:J])
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 

That doesn't work at all! This is because the Y[,-(1:J)] are not at all
orthogonal to one another, and the correctY function really assumes they
are orthogonal.

# Miscellenous

The remainder of this are just notes for future reference,
and should probably be ignored.

## An alternative approach to controlling for PCs

We could assume $a_{kg} \sim N(0, \sigma_g^2 \lambda_k)$
and estimate $\sigma_g$ and $\lambda_k$ from data.
This might involve
$\hat{a}_{kg} \sim N(0, \sigma_g^2 \lambda_k + s_{kg})$.


## Possible sparse linear regression this way?

Maybe we can do something similar for sparse regression?
So $Y=X\beta + \epsilon$.
Estimate $\beta$ separately for each $\beta$, and
then shrink the $\hat{beta}$ values. 
Then reestimate the standard errors?

Might want to be able to do something like this based
on summary statistics alone?

Might want to think about approximating posterior
on $X\beta$ this way. Then posterior on any given
$\beta$ can be obtained?


### Possible SFA this way?

Another thought along these lines is perhaps
we can do a version of sparse factor analysis this way.
Recall that the idea is to fit $Y= \Lambda F + E$
where $Y$ is $n$ by $p$, $\Lambda$ is $n \times k$ and $F$ is
$k \times p$, and both $\Lambda$ and/or $F$ may be small.
Given $F$, we can use regression to estimate the mean and sd of each element of $\Lambda$.
Then we can apply ASH to shrink the lambda values for each factor. 
Then given $\Lambda$ we can use regression to estimate the mean and sd of each $F$, and use ASH to adaptively shrink the $F$ values?
Since the $k$ $F$s and $k$ $\Lambdas$ will in general not be orthogonal we might want to worry about this, but we might also ignore it as a first try.



