
rm(list=ls())
gc()
library(mvtnorm)
library(glmnet)
library(BDcocolasso)
# https://github.com/celiaescribe/BDcocolasso
library(pqr)
# https://github.com/xliusufe/pqr/

K <- 5
nK <- c(2, 3, 4, 5, 6) * 1e4
n <- sum(nK)
p <- 600
meanx <- rep(0, p)
corrx <- 0.5
sigmax <- matrix(0, p, p)
for (i in 1:p) {
  for (j in 1:p) {
    sigmax[i, j] <- corrx ^ (abs(i - j))
  }
}
d <- 3
s <- 8
bt <- rep(2, d)
gm <- c(2, rep(2, s - d), rep(0, p - s))
meanv <- rep(0, d)
sigmav <- diag(rep(0.5, d))
r1.set <- 1000
itr.max <- 500
bhat.dunif <- bhat.unif <- array(NA, c(itr.max, length(r1.set), d))

# data generation
UK <- array(NA, c(max(nK), K, p + 1))
XK <- array(NA, c(max(nK), K, d))
VK <- array(NA, c(max(nK), K, d))
WK <- array(NA, c(max(nK), K, d))
ZK <- array(NA, c(max(nK), K, p - d + 1))
eK <- matrix(NA, max(nK), K)
yK <- matrix(NA, max(nK), K)
for (k in 1:K) {
  UK[1:nK[k], k, ] <- cbind(1, rmvnorm(nK[k], meanx, sigmax))
  XK[1:nK[k], k, ] <- UK[1:nK[k], k, 2:(d + 1)]
  VK[1:nK[k], k, ] <- rmvnorm(nK[k], meanv, sigmav)
  ZK[1:nK[k], k, ] <- UK[1:nK[k], k, -c(2:(d + 1))]
  eK[1:nK[k], k] <- rnorm(nK[k], mean=0, sd=3)
  yK[1:nK[k], k] <- XK[1:nK[k], k, ] %*% bt + ZK[1:nK[k], k, ] %*% gm + eK[1:nK[k], k]
  WK[1:nK[k], k, ] <- XK[1:nK[k], k, ] + VK[1:nK[k], k, ]
}

time1 <- Sys.time()
for (itr in 1:itr.max) {
  
  cat(itr)
  for (r1.idx in 1:length(r1.set)) {
    
    r1 <- r1.set[r1.idx]
    # first step subsampling
    WK.simp <- c()
    ZK.simp <- c()
    yK.simp <- c()
    r1.simp <- round(r1 * nK / n)
    idx.simp <- matrix(NA, max(r1.simp), K)
    for (k in 1:K) {
      idx.simp[1:r1.simp[k], k] <- sample(1:nK[k], r1.simp[k], T)
      WK.simp <- rbind(WK.simp, WK[idx.simp[1:r1.simp[k], k], k, ])
      ZK.simp <- rbind(ZK.simp, ZK[idx.simp[1:r1.simp[k], k], k, ])
      yK.simp <- c(yK.simp, yK[idx.simp[1:r1.simp[k], k], k])
    }

    # cocolasso using the cross validation lambda
    cv.fit <- coco(cbind(ZK.simp, WK.simp), y=as.matrix(yK.simp), n=length(yK.simp),
                   p=p + 1, p1=p - d + 1, p2=d, center.Z=F, scale.Z=F, center.y=F, scale.y=F,
                   K=4, tau=sigmav, earlyStopping_max=3, noise="additive", block=T, penalty="lasso", mode="HM")
    theta.ini <- c(cv.fit$beta.opt[-c(1:(p - d + 1))], cv.fit$beta.opt[1:(p - d + 1)])
    # theta.ini[abs(theta.ini) < 0.5] <- 0
    # plot(theta.ini)
    
    # refit
    # if (sum(theta.ini[1:d]!=0)==d) {
    #   index <- theta.ini!=0
    #   mat <- cbind(WK.simp, ZK.simp)[, index]
    #   sigmahat <- t(mat) %*% mat / r1 - diag(c(diag(sigmav), rep(0, dim(mat)[2] - d)))
    #   rhohat <- t(mat) %*% yK.simp / r1
    #   if (min(eigen(sigmahat)$values) > 0) {
    #     theta.refit <- solve(sigmahat) %*% rhohat
    #   } else {
    #     refit <- lasso_covariance_block(n=r1, p1=dim(mat)[2] - d, p2=d,
    #                                     X1=mat[, -c(1:d)], Z2=mat[, 1:d], y=yK.simp, sigma1=t(mat[, -c(1:d)]) %*% mat[, -c(1:d)] / r1,
    #                                     sigma2=HM_proj(t(mat[, 1:d]) %*% mat[, 1:d] / r1 - sigmav),
    #                                     lambda=cv.fit$lambda.opt,
    #                                     noise="additive", beta1.start=rep(0, dim(mat)[2] - d),
    #                                     beta2.start=rep(0, d), penalty="lasso")
    #     theta.refit <- c(refit$coefficients.beta2, refit$coefficients.beta1)
    #   }
    #   theta.refit[abs(theta.refit) < 0.5] <- 0
    #   theta.ini[index] <- theta.refit
    # }
    # plot(theta.ini)
    bttilde <- theta.ini[c(1:d)]
    gmtilde <- theta.ini[-c(1:d)]
    bhat.unif[itr, r1.idx, ] <- theta.ini[c(1:d)]

    # estimate of projection matrix
    Hfit <- mvr(WK.simp, ZK.simp[, -1])
    Hhat <- rbind(as.vector(Hfit$muhat), Hfit$Bhat)
    Htilde <- t(Hhat)
    
    tmp1.dunif <- tmp2.dunif <- 0
    for (i in 1:length(yK.simp)) {
      tmp1.dunif <- tmp1.dunif + 
        ((WK.simp[i, ] - Htilde %*% ZK.simp[i, ]) %*% t(WK.simp[i, ]) - sigmav)
      tmp2.dunif <- tmp2.dunif +
        ((WK.simp[i, ] - Htilde %*% ZK.simp[i, ]) * c(yK.simp[i] - gmtilde %*% ZK.simp[i, ]))
    }
    bhat.dunif[itr, r1.idx, ] <- c(solve(tmp1.dunif) %*% tmp2.dunif)
  }
}
time2 <- Sys.time()

mse.dunif <- mse.unif <- absbias.dunif <- absbias.unif <- rep(NA, length(r1.set))
for (r1.idx in 1:length(r1.set)) {
  mse.dunif[r1.idx] <- 
    mean(diag((bhat.dunif[, r1.idx, ] - rep(1, itr.max) %*% t(bt)) %*% t(bhat.dunif[, r1.idx, ] - rep(1, itr.max) %*% t(bt))))
  mse.unif[r1.idx] <- 
    mean(diag((bhat.unif[, r1.idx, ] - rep(1, itr.max) %*% t(bt)) %*% t(bhat.unif[, r1.idx, ] - rep(1, itr.max) %*% t(bt))))
  absbias.dunif[r1.idx] <- mean(abs(apply(bhat.dunif[, r1.idx, ], 2, mean) - bt))
  absbias.unif[r1.idx] <- mean(abs(apply(bhat.unif[, r1.idx, ], 2, mean) - bt))
}
time2 - time1
c(mse.dunif[1], mse.unif[1])
c(absbias.dunif[1], absbias.unif[1])

