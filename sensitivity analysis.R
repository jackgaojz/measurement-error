
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
r1 <- 400

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

# initial estimate of beta
# cocolasso using the selected lambda
# fit <- lasso_covariance_block(n=r1, p1=p - d + 1, p2=d,
#   X1=z.simp, Z2=w.simp, y=y.simp, sigma1=t(z.simp) %*% z.simp / r1,
#   sigma2=HM_proj(t(w.simp) %*% w.simp / r1 - sigmav),
#   lambda=sqrt(log(p) / r1) * 2,
#   noise="additive", beta1.start=rep(0, p - d + 1),
#   beta2.start=rep(0, d), penalty="lasso")
# theta.ini <- c(fit$coefficients.beta2, fit$coefficients.beta1)
# theta.ini[abs(theta.ini) < 0.5] <- 0
# plot(theta.ini)

# cocolasso using the cross validation lambda
cv.fit <- coco(cbind(ZK.simp, WK.simp), y=as.matrix(yK.simp), n=length(yK.simp),
               p=p + 1, p1=p - d + 1, p2=d, center.Z=F, scale.Z=F, center.y=F, scale.y=F,
               K=4, tau=sigmav, earlyStopping_max=3, noise="additive", block=T, penalty="lasso", mode="HM")
theta.ini <- c(cv.fit$beta.opt[-c(1:(p - d + 1))], cv.fit$beta.opt[1:(p - d + 1)])
theta.ini[abs(theta.ini) < 0.5] <- 0
# plot(theta.ini)

# refit
if (sum(theta.ini[1:d]!=0)==d) {
  index <- theta.ini!=0
  mat <- cbind(WK.simp, ZK.simp)[, index]
  sigmahat <- t(mat) %*% mat / r1 - diag(c(diag(sigmav), rep(0, dim(mat)[2] - d)))
  rhohat <- t(mat) %*% yK.simp / r1
  if (min(eigen(sigmahat)$values) > 0) {
    theta.refit <- solve(sigmahat) %*% rhohat
  } else {
    refit <- lasso_covariance_block(n=r1, p1=dim(mat)[2] - d, p2=d,
                                    X1=mat[, -c(1:d)], Z2=mat[, 1:d], y=y.simp, sigma1=t(mat[, -c(1:d)]) %*% mat[, -c(1:d)] / r1,
                                    sigma2=HM_proj(t(mat[, 1:d]) %*% mat[, 1:d] / r1 - sigmav),
                                    lambda=cv.fit$lambda.opt,
                                    noise="additive", beta1.start=rep(0, dim(mat)[2] - d),
                                    beta2.start=rep(0, d), penalty="lasso")
    theta.refit <- c(refit$coefficients.beta2, refit$coefficients.beta1)
  }
  theta.refit[abs(theta.refit) < 0.5] <- 0
  theta.ini[index] <- theta.refit
}
# plot(theta.ini)
bttilde <- theta.ini[c(1:d)]
gmtilde <- theta.ini[-c(1:d)]

# lasso
# cv.fit.ini <- cv.glmnet(cbind(XK.simp, ZK.simp[, -1]), yK.simp, alpha=1, family="gaussian")
# fit.ini <- glmnet(cbind(XK.simp, ZK.simp[, -1]), yK.simp, alpha=1, lambda=cv.fit.ini$lambda.min)
# theta.ini <- as.vector(c(fit.ini$beta[1:d], fit.ini$a0, fit.ini$beta[-(1:d)]))
# theta.ini[abs(theta.ini) < 0.1] <- 0
# plot(theta.ini)

# estimate of projection matrix
Hfit <- mvr(WK.simp, ZK.simp[, -1])
Hhat <- rbind(as.vector(Hfit$muhat), Hfit$Bhat)
Htilde <- t(Hhat)

r2.set <- 1000
itr.max <- 500

lambda.theta <- c(cv.fit$lambda.opt - 0.10, cv.fit$lambda.opt - 0.05, cv.fit$lambda.opt,
                  cv.fit$lambda.opt + 0.05, cv.fit$lambda.opt + 0.10)
lambda.h <- c(Hfit$lambda_opt - 0.10, Hfit$lambda_opt - 0.05, Hfit$lambda_opt,
              Hfit$lambda_opt + 0.05, Hfit$lambda_opt + 0.10)
bhat.dunif <- bhat.dmV <- bhat.dmVc <- array(NA, c(itr.max, d, length(r2.set), length(lambda.theta), length(lambda.h)))

time1 <- Sys.time()
for (itr in 1:itr.max) {
  
  cat(itr)
  for (r2.idx in 1:length(r2.set)) {
    r2 <- r2.set[r2.idx]
    for (idx.lambda.theta in 1:length(lambda.theta)) {
      for (idx.lambda.h in 1:length(lambda.h)) {
        
        # cocolasso using the selected lambda
        fit <- lasso_covariance_block(n=r1, p1=p - d + 1, p2=d,
                                      X1=ZK.simp, Z2=WK.simp, y=yK.simp, sigma1=t(ZK.simp) %*% ZK.simp / r1,
                                      sigma2=HM_proj(t(WK.simp) %*% WK.simp / r1 - sigmav),
                                      lambda=lambda.theta[idx.lambda.theta],
                                      noise="additive", beta1.start=rep(0, p - d + 1),
                                      beta2.start=rep(0, d), penalty="lasso")
        theta.ini <- c(fit$coefficients.beta2, fit$coefficients.beta1)
        # theta.ini[abs(theta.ini) < 0.5] <- 0
        # plot(theta.ini)
        
        # refit
        if (sum(theta.ini[1:d]!=0)==d) {
          index <- theta.ini!=0
          mat <- cbind(WK.simp, ZK.simp)[, index]
          sigmahat <- t(mat) %*% mat / r1 - diag(c(diag(sigmav), rep(0, dim(mat)[2] - d)))
          rhohat <- t(mat) %*% yK.simp / r1
          if (min(eigen(sigmahat)$values) > 0) {
            theta.refit <- solve(sigmahat) %*% rhohat
          } else {
            refit <- lasso_covariance_block(n=r1, p1=dim(mat)[2] - d, p2=d,
                                            X1=mat[, -c(1:d)], Z2=mat[, 1:d], y=yK.simp, sigma1=t(mat[, -c(1:d)]) %*% mat[, -c(1:d)] / r1,
                                            sigma2=HM_proj(t(mat[, 1:d]) %*% mat[, 1:d] / r1 - sigmav),
                                            lambda=lambda.theta[itr, idx.lambda.theta],
                                            noise="additive", beta1.start=rep(0, dim(mat)[2] - d),
                                            beta2.start=rep(0, d), penalty="lasso")
            theta.refit <- c(refit$coefficients.beta2, refit$coefficients.beta1)
          }
          theta.refit[abs(theta.refit) < 0.5] <- 0
          theta.ini[index] <- theta.refit
        }
        # plot(theta.ini)
        
        bttilde <- theta.ini[c(1:d)]
        gmtilde <- theta.ini[-c(1:d)]
        
        # estimate of projection matrix
        Hfit <- mvr(WK.simp, ZK.simp[, -1], lambda=lambda.h[idx.lambda.h])
        Hhat <- rbind(as.vector(Hfit$muhat), Hfit$Bhat)
        Htilde <- t(Hhat)
        
        # method 1: dunif
        PI.opt.dunif <- matrix(NA, max(nK), K)
        r2.dunif <- rep(NA, K)
        for (k in 1:K) {
          PI.opt.dunif[1:nK[k], k] <- 1 / nK[k]
          r2.dunif[k] <- round(r2 * nK[k] / n)
        }
        WK.dunif <- c()
        ZK.dunif <- c()
        yK.dunif <- c()
        PI.dunif <- c()
        rK.dunif <- c()
        idx.dunif <- matrix(NA, max(r2.dunif), K)
        for (k in 1:K) {
          idx.dunif[1:r2.dunif[k], k] <- sample(1:nK[k], r2.dunif[k], T)
          WK.dunif <- rbind(WK.dunif, WK[idx.dunif[1:r2.dunif[k], k], k, ])
          ZK.dunif <- rbind(ZK.dunif, ZK[idx.dunif[1:r2.dunif[k], k], k, ])
          yK.dunif <- c(yK.dunif, yK[idx.dunif[1:r2.dunif[k], k], k])
          PI.dunif <- c(PI.dunif, PI.opt.dunif[idx.dunif[1:r2.dunif[k], k], k])
          rK.dunif <- c(rK.dunif, rep(r2.dunif[k], r2.dunif[k]))
        }
        
        tmp1.dunif <- tmp2.dunif <- 0
        for (i in 1:length(yK.dunif)) {
          tmp1.dunif <- tmp1.dunif + 
            ((WK.dunif[i, ] - Htilde %*% ZK.dunif[i, ]) %*% t(WK.dunif[i, ]) - sigmav) / (rK.dunif[i] * PI.dunif[i])
          tmp2.dunif <- tmp2.dunif +
            ((WK.dunif[i, ] - Htilde %*% ZK.dunif[i, ]) * c(yK.dunif[i] - gmtilde %*% ZK.dunif[i, ])) / (rK.dunif[i] * PI.dunif[i])
        }
        bhat.dunif[itr, , r2.idx, idx.lambda.theta, idx.lambda.h] <- c(solve(tmp1.dunif) %*% tmp2.dunif)
        
        # method 2: dmV
        Jtilde <- t(WK.simp) %*% (WK.simp - ZK.simp %*% t(Htilde)) / r1 - sigmav
        PI.opt.dmV <- matrix(NA, max(nK), K)
        r.opt.dmV <- rep(NA, K)
        
        tmp.dmV <- matrix(NA, max(nK), K)
        for (k in 1:K) {
          tmp.dmV[1:nK[k], k] <- sqrt(rowSums(((WK[1:nK[k], k, ] - ZK[1:nK[k], k, ] %*% t(Htilde)) %*% solve(Jtilde) * 
                                                 c(WK[1:nK[k], k, ] %*% bttilde + ZK[1:nK[k], k, ] %*% gmtilde - yK[1:nK[k], k]) - 
                                                 t(matrix(1, d, nK[k]) * c(solve(Jtilde) %*% sigmav %*% bttilde))) ^ 2))
          PI.opt.dmV[1:nK[k], k] <- tmp.dmV[1:nK[k], k] / sum(tmp.dmV[1:nK[k], k])
          r.opt.dmV[k] <- sum(tmp.dmV[1:nK[k], k])
        }
        r.opt.dmV <- round(r2 * r.opt.dmV / sum(r.opt.dmV))
        
        WK.dmV <- c()
        ZK.dmV <- c()
        yK.dmV <- c()
        PI.dmV <- c()
        rK.dmV <- c()
        idx.dmV <- matrix(NA, max(r.opt.dmV), K)
        for (k in 1:K) {
          idx.dmV[1:r.opt.dmV[k], k] <- sample(1:nK[k], r.opt.dmV[k], T, PI.opt.dmV[1:nK[k], k])
          WK.dmV <- rbind(WK.dmV, WK[idx.dmV[1:r.opt.dmV[k], k], k, ])
          ZK.dmV <- rbind(ZK.dmV, ZK[idx.dmV[1:r.opt.dmV[k], k], k, ])
          yK.dmV <- c(yK.dmV, yK[idx.dmV[1:r.opt.dmV[k], k], k])
          PI.dmV <- c(PI.dmV, PI.opt.dmV[idx.dmV[1:r.opt.dmV[k], k], k])
          rK.dmV <- c(rK.dmV, rep(r.opt.dmV[k], r.opt.dmV[k]))
        }
        
        tmp1.dmV <- tmp2.dmV <- 0
        for (i in 1:length(yK.dmV)) {
          tmp1.dmV <- tmp1.dmV + 
            ((WK.dmV[i, ] - Htilde %*% ZK.dmV[i, ]) %*% t(WK.dmV[i, ]) - sigmav) / (rK.dmV[i] * PI.dmV[i])
          tmp2.dmV <- tmp2.dmV +
            ((WK.dmV[i, ] - Htilde %*% ZK.dmV[i, ]) * c(yK.dmV[i] - gmtilde %*% ZK.dmV[i, ])) / (rK.dmV[i] * PI.dmV[i])
        }
        bhat.dmV[itr, , r2.idx, idx.lambda.theta, idx.lambda.h] <- c(solve(tmp1.dmV) %*% tmp2.dmV)
        
        # method 3: dmVc
        PI.opt.dmVc <- matrix(NA, max(nK), K)
        r.opt.dmVc <- rep(NA, K) 
        
        tmp.dmVc <- matrix(NA, max(nK), K)
        for (k in 1:K) {
          tmp.dmVc[1:nK[k], k] <- sqrt(rowSums(((WK[1:nK[k], k, ] - ZK[1:nK[k], k, ] %*% t(Htilde)) * 
                                                  c(WK[1:nK[k], k, ] %*% bttilde + ZK[1:nK[k], k, ] %*% gmtilde - yK[1:nK[k], k]) - 
                                                  t(matrix(1, d, nK[k]) * c(sigmav %*% bttilde))) ^ 2))
          PI.opt.dmVc[1:nK[k], k] <- tmp.dmVc[1:nK[k], k] / sum(tmp.dmVc[1:nK[k], k])
          r.opt.dmVc[k] <- sum(tmp.dmVc[1:nK[k], k])
        }
        r.opt.dmVc <- round(r2 * r.opt.dmVc / sum(r.opt.dmVc))
        
        WK.dmVc <- c()
        ZK.dmVc <- c()
        yK.dmVc <- c()
        PI.dmVc <- c()
        rK.dmVc <- c()
        idx.dmVc <- matrix(NA, max(r.opt.dmVc), K)
        for (k in 1:K) {
          idx.dmVc[1:r.opt.dmVc[k], k] <- sample(1:nK[k], r.opt.dmVc[k], T, PI.opt.dmVc[1:nK[k], k])
          WK.dmVc <- rbind(WK.dmVc, WK[idx.dmVc[1:r.opt.dmVc[k], k], k, ])
          ZK.dmVc <- rbind(ZK.dmVc, ZK[idx.dmVc[1:r.opt.dmVc[k], k], k, ])
          yK.dmVc <- c(yK.dmVc, yK[idx.dmVc[1:r.opt.dmVc[k], k], k])
          PI.dmVc <- c(PI.dmVc, PI.opt.dmVc[idx.dmVc[1:r.opt.dmVc[k], k], k])
          rK.dmVc <- c(rK.dmVc, rep(r.opt.dmVc[k], r.opt.dmVc[k]))
        }
        
        tmp1.dmVc <- tmp2.dmVc <- 0
        for (i in 1:length(yK.dmVc)) {
          tmp1.dmVc <- tmp1.dmVc + 
            ((WK.dmVc[i, ] - Htilde %*% ZK.dmVc[i, ]) %*% t(WK.dmVc[i, ]) - sigmav) / (rK.dmVc[i] * PI.dmVc[i])
          tmp2.dmVc <- tmp2.dmVc +
            ((WK.dmVc[i, ] - Htilde %*% ZK.dmVc[i, ]) * c(yK.dmVc[i] - gmtilde %*% ZK.dmVc[i, ])) / (rK.dmVc[i] * PI.dmVc[i])
        }
        bhat.dmVc[itr, , r2.idx, idx.lambda.theta, idx.lambda.h] <- c(solve(tmp1.dmVc) %*% tmp2.dmVc)
      }
    }
  }
}
time2 <- Sys.time()

mse.dunif <- mse.dmV <- mse.dmVc <- array(NA, c(length(r2.set), length(lambda.theta), length(lambda.h)))
for (idx.r2 in 1:length(r2.set)) {
  for (idx.lambda.theta in 1:length(lambda.theta)) {
    for (idx.lambda.h in 1:length(lambda.h)) {
      mse.dunif[idx.r2, idx.lambda.theta, idx.lambda.h] <- 
        mean(diag((bhat.dunif[, , idx.r2, idx.lambda.theta, idx.lambda.h] - rep(1, itr.max) %*% t(bt)) %*% t(bhat.dunif[, , idx.r2, idx.lambda.theta, idx.lambda.h] - rep(1, itr.max) %*% t(bt))))
      mse.dmV[idx.r2, idx.lambda.theta, idx.lambda.h] <- 
        mean(diag((bhat.dmV[, , idx.r2, idx.lambda.theta, idx.lambda.h] - rep(1, itr.max) %*% t(bt)) %*% t(bhat.dmV[, , idx.r2, idx.lambda.theta, idx.lambda.h] - rep(1, itr.max) %*% t(bt))))
      mse.dmVc[idx.r2, idx.lambda.theta, idx.lambda.h] <- 
        mean(diag((bhat.dmVc[, , idx.r2, idx.lambda.theta, idx.lambda.h] - rep(1, itr.max) %*% t(bt)) %*% t(bhat.dmVc[, , idx.r2, idx.lambda.theta, idx.lambda.h] - rep(1, itr.max) %*% t(bt))))
    }
  }
}

for (idx.r2 in 1:length(r2.set)) {
  print(mse.dunif[idx.r2, , ])
  print(mse.dmV[idx.r2, , ])
  print(mse.dmVc[idx.r2, , ])
}

time2 - time1

