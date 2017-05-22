# Cai, Jian-Feng, Emmanuel J. Candès, and Zuowei Shen. "A singular value thresholding algorithm for matrix completion." SIAM Journal on Optimization 20.4 (2010): 1956-1982.

## WARNING sketchy prototype ahead...

library(Matrix)
library(irlba)
svt = function(M,  delta=1.2, epsilon=1e-3, tau=5 * nrow(M), kmax=200)
{
  F = norm(M, type="F")
  dp = diff(M@p)
  idx = cbind(M@i + 1, rep(seq_along(dp), dp)) # the non-zero indices
  k0 = max(floor(tau / (delta * irlba(M, 1)$d)), 1)
  Y = k0 * delta * M
  r = 0
  err = sqrt(drop(crossprod(M[idx])))
  for(k in seq(1, kmax))
  {
    S = NULL
    s = r + 1
    dmin = tau + 1
    while(dmin > tau)
    {
      S = irlba(Y, nv=s, v=S)
      dmin = min(S$d)
      s = s + 5
    }
    r = max(which(S$d > tau), 1)
cat(k, s, r, err, "\n") # debug
    X = S$u %*% ((S$d - tau) * t(S$v))
    err = sqrt(drop(crossprod(M[idx] - X[idx]))) / F
    if(err < epsilon) break
    Y[idx] = Y[idx] + delta * (M[idx] - X[idx])
  }
  list(X=X, err=err, k=k)
}

# EXAMPLE
set.seed(1)
r = 10
n = 1000
Mtrue = matrix(rnorm(n * r), n) %*% matrix(rnorm(n * r), r)
p = 0.1 * n * n
i = sample(n, p, replace=TRUE)
j = sample(n, p, replace=TRUE)
M = sparseMatrix(i=i, j=j, x=Mtrue[cbind(i, j)], dims=c(n, n), use.last.ij=TRUE)

x = svt(M)
cat("Relative error ||X - Mtrue||_F / ||Mtrue||_F :\n")
print(norm(Mtrue - x$X, "F") / norm(Mtrue, "F"))
sd = svd(Mtrue)$d
xd = svd(x$X)$d

plot(sd)
lines(xd)
legend("topright", legend=c("true singular values", "imputed values"), pch=c("o", "lines"))
