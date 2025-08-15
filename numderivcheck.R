library(numDeriv)
library(foreach)

set.seed(123)
k = 3

Rr = cor(matrix(rnorm((k+2)*k), nrow = k+2))
Rri = solve(Rr)

e.out = eigen(Rr)
lambda = e.out$values
W = e.out$vectors

lambda = sort(rchisq(k, df = 14), decreasing = TRUE)
Rq = W %*% diag(lambda) %*% t(W)
Rqi = W %*% diag(1/lambda) %*% t(W)
RqD = diag(sqrt(diag(Rq)))

makeRmatrix = function(loglambda, W, corr = FALSE) {
  Rrcc = W %*% diag(exp(loglambda)) %*% t(W)
  Rrcci = W %*% diag(exp(-loglambda)) %*% t(W)
  Rrc = diag(sqrt(diag(Rrcc))) %*% Rrcci %*% diag(sqrt(diag(Rrcc)))
  if (corr) Rrc = solve(Rrc)
  return(c(Rrc))
}

jb0 = t(jacobian(makeRmatrix, log(lambda), W = W))



jb1 = foreach(i1=1:k, .combine = rbind) %do% {
  DDS = diag(0.5 * lambda[i1] * W[ , i1]^2 / sqrt(diag(Rq)))
  c(
    -RqD %*% outer(W[,i1] , W[,i1] / lambda[i1]) %*% RqD +
      DDS %*% Rqi %*% RqD + RqD %*% Rqi %*% DDS
  )
}


jb2= foreach(i1=1:k, .combine = rbind) %do% {
  DDS = diag(0.5 * lambda[i1] * W[ , i1]^2 / sqrt(diag(Rq)))
  c(
    -RqD %*% outer(W[,i1] , W[,i1] / lambda[i1]) %*% RqD +
      2 * DDS %*% Rqi %*% RqD
  )
}

  

print(jb0)
print(jb1)
print(jb2)
