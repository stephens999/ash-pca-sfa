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
You'll have to download ASH from my github repo BayesFDR, and look
at the readme to see what it is doing.
But bottom line is that it will provide "shrunk" estimates for $ahat_{kg}$.

3. Compute the residuals $r_g = Y_g - \sum_k a_{kg} v_k$.
(3b. later - we could possibly use $r_g$ to estimate $\sigma_g$, the residual variance for $Y_g$,
and then this could be used to get a revised estimate of the standard error $s_{kg}$ )

4. Test the residuals for association with genotype.

Here is a start at writing R code for this:

```r
# INPUT: Y an n by g matrix of phenotypes X an n by K matrix of covariates
# (Typically columns of X will be the PCs, and K=n since there are n PCs)
# scaleY a bool saying whether to scale Y to have unit variance
correctY = function(Y, X, scaleY = FALSE) {
    X = scale(X, scale = FALSE)  #remove column means
    Y = scale(Y, scale = scaleY)  #remove column means, and
    SX2 = colSums(X * X)  #a K vector of sums of squares
    EY2 = colMeans(Y * Y)  #variance of Y
    betahat = SX2^{
        -1
    } * (t(X) %*% Y)  # a K by g matrix
    resid = Y - X %*% betahat  # residuals, K by g
    sigmahat = colMeans(resid * resid)  # K vector, estimated residual variances
    sebetahat = SX2^{
        -1
    } * sigmahat
    return(list(betahat, sebetahat))
}
```


Here we'll simulate some data based on latent factors $Z$

```r
n = 100
G = 1000
K = 14
L = matrix(rnorm(K * n), nrow = n)  #n by k matrix of loadings
F = matrix(rnorm(K * G), nrow = K)  #K by g matrix of factors
Y = L %*% F + matrix(rnorm(n * G), nrow = n)
Y.pr = prcomp(Y)
```



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



