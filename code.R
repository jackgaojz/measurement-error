
rm(list=ls())
gc()
library(mvtnorm)
library(glmnet)
library(BDcocolasso)
# https://github.com/celiaescribe/BDcocolasso
library(pqr)
# https://github.com/xliusufe/pqr/

########## If BDcocolasso cannot be installed ##########
l1proj <- function(v, b) {
  
  # Efficient projection onto L1 ball of specified radius (i.e. b), used by the admm algo
  # Ref. Duchi et al. (2008). Efficient Projections onto the L1-Ball for Learning in High Dimensions, ICML
  
  stopifnot(b>0)
  
  u <- sort(abs(v),decreasing=TRUE)
  sv <- cumsum(u)
  rho <- max(which(u>(sv-b)/1:length(u)))
  theta <- max(0, (sv[rho]-b)/rho)
  w <-sign(v) * pmax(abs(v)-theta,0)
  
  return(w)
}

ADMM_proj <- function(mat, epsilon=1e-4, mu=10, it.max=1e3, etol=1e-4, etol_distance=1e-4) {
  
  
  
  p<-nrow(mat)
  
  # Initialization
  R<-diag(mat)
  S<-matrix(0,p,p)
  L<-matrix(0,p,p)
  
  itr<-0
  iteration <- eps_R <- eps_S <- eps_primal <- time <- distance <- NULL
  while (itr<it.max) {
    #print(itr)
    Rp<-R
    Sp<-S
    start <- Sys.time()
    # Subproblem I: R step
    W<-mat+S+mu*L
    W.eigdec<-eigen(W, symmetric=TRUE) 
    W.V<-W.eigdec$vectors
    W.D<-W.eigdec$values
    R<-W.V%*%diag(pmax(W.D,epsilon))%*%t(W.V)
    
    # Subproblem II: S step
    M<-R-mat-mu*L     
    S[lower.tri(S, diag = TRUE)]<-M[lower.tri(M, diag = TRUE)]-l1proj(v=M[lower.tri(M, diag = TRUE)],b=mu/2)    
    for (i in 2:p){
      for (j in 1:(i-1)){
        S[j,i]<-S[i,j]
      }
    }
    
    # L step: update the Lagrange parameter
    L<-L-(R-S-mat)/mu
    end <- Sys.time()
    #Stocking the values of different parameters with the number of iterations
    iteration <- c(iteration, itr)
    eps_R <- c(eps_R,max(abs(R-Rp)))
    eps_S <- c(eps_S,max(abs(S-Sp)))
    eps_primal <- c(eps_primal, max(abs(R-S-mat)))
    time <- c(time, end - start)
    distance <- c(distance,max(abs(R-mat)))
    
    # Stopping Rule                        
    #cat("check the stopping criterion:",max(abs(R-S-mat)),"\n")
    if (((max(abs(R-Rp))<etol) && (max(abs(S-Sp))<etol) && (max(abs(R-S-mat))<etol)) || (abs(max(abs(Rp-mat)) - max(abs(R-mat)))<etol_distance)){
      itr<-it.max
    } else {
      itr<-itr+1
    }
    
    if (itr%%20==0) {
      mu<-mu/2
    }
  }
  df_ADMM <- data.frame(iteration = iteration, eps_R = eps_R, eps_S=eps_S, eps_primal=eps_primal, time=time, distance=distance)
  return(list(mat=R,df_ADMM=df_ADMM))
  
}

getRmat <- function(X, y=NA) {
  m=dim(X)[1]; n=dim(X)[2]
  R=matrix(0,n,n); rho_paired=c()
  for (j in (1:n)){ #col
    for (i in (1:j)){#row
      colProd=X[,i]*X[,j]; nij=sum(!is.na(colProd))
      R[i,j]=(nij/m)
      R[j,i]=R[i,j]
    }
    rho_j=(sum(X[,j]*y,na.rm=T))/nij
    rho_paired=c(rho_paired,rho_j)
  }
  return(list(rho_paired,R))
}

updateB_maxNorm_hmLasso <- function(Akp1, Lk, S_paired, mu, W) {
  C=Akp1-S_paired-mu*Lk
  C_vec=as.vector(C); W_vec=as.vector(W); n=dim(C)[1]; m=dim(C)[2]
  WC=abs(C_vec)*W_vec
  WC_sort_inx=order(WC,decreasing = T);
  W_sort=W_vec[WC_sort_inx]; C_sort=C_vec[WC_sort_inx]
  l=1; frac=0
  while (W_sort[l]*abs(C_sort[l])> frac && l<=length(C_vec)){
    frac=(sum(abs(C_sort[1:l]))-mu/2)/(sum(1/(W_sort[1:l])))
    l=l+1
  }
  d=frac
  
  b_vec=mapply(FUN=function(t) ifelse(W_vec[t]*abs(C_vec[t])>d, d*sign(C_vec[t])/W_vec[t],
                                      C_vec[t]),1:length(W_sort))
  Bkp1=matrix(b_vec,n,m)
  return(Bkp1)
}

HM_proj <- function(sigmaHat, R=NULL, a=1, iter_max=1000, epsilon=1e-4, mu=10, tolerance=1e-4, norm="F") {
  iter=0
  S_paired=sigmaHat; n=nrow(S_paired)
  if (is.null(R)){
    W=matrix(1,n,n)^a
  } else{ W=R^a}
  
  Ak=S_paired
  Bk=matrix(0,n,n); Lk=matrix(0,n,n)
  while (iter<iter_max){
    A=Bk+S_paired+mu*Lk
    A_eigdec=eigen(A, symmetric=TRUE) 
    A_V=A_eigdec$vectors
    A_D=A_eigdec$values
    Akp1=A_V%*%diag(pmax(A_D,epsilon))%*%t(A_V)
    
    if (norm =='F'){
      Bkp1 = (Akp1-S_paired-mu*Lk)/(mu*W*W+matrix(1,n,n)) 
    } else {
      Bkp1 = updateB_maxNorm_hmLasso(Akp1,Lk,S_paired,mu,W)
    }
    Lkp1=Lk-(Akp1-Bkp1-S_paired)/mu
    
    if (max(max(abs(Akp1-Ak)),max(abs(Bkp1-Bk)),max(abs(Lkp1-Lk)))<tolerance ){
      iter=iter_max
    } else { iter=iter+1}
    Ak=Akp1; Bk=Bkp1; Lk=Lkp1
  }
  return(Ak)
}

lambda_max.coordinate_descent <- function(Z, y, n, p, p1, p2, ratio_matrix=NULL, noise=c("additive", "missing")) {
  start <- p1 + 1
  X1 <- Z[,1:p1]
  Z2 <- Z[,start:p]
  if (noise=="additive"){
    rho_tilde <- 1/n*t(Z)%*%y
    lambda <- max(abs(rho_tilde))
  }else{
    rho_tilde1 <- 1/n*t(X1)%*%y
    rho_tilde2 <- 1/n*t(Z2)%*%y/ diag(ratio_matrix)
    lambda = max(max(abs(rho_tilde1)),max(abs(rho_tilde2)))
  }
  lambda
}

scale_manual <- function(j, Z) {
  sd <- stats::sd(Z[,j])
  if (sd != 0){
    return(Z[,j]/sd)
  }else{
    return (Z[,j])
  }
  
}

rescale_without_NA_block <- function(j, Z, p1) {
  m <- mean(Z[which(!is.na(Z[,p1 + j]), arr.ind = TRUE),p1+j])
  Z[,p1 + j] - m
}

mean_without_NA <- function(j, Z) {
  m <- mean(Z[which(!is.na(Z[,j]), arr.ind = TRUE),j])
  m
}

sd_without_NA_block <- function(j, Z) {
  sd <- stats::sd(Z[which(!is.na(Z[,j]), arr.ind = TRUE),j])
  sd
}

change_NA_value_block <- function(j, Z, p1) {
  Z[which(is.na(Z[,p1 + j]), arr.ind = TRUE),j + p1] <- 0
  Z[,p1 + j]
}

scale_manual_with_sd <- function(j, Z, v) {
  sd <- v[j]
  if (sd != 0){
    return(Z[,j]/sd)
  }else{
    return (Z[,j])
  }
}

cross_validation_function.block_descent <- function(k, Z, y, n, n_one_fold, n_without_fold, p, p1, p2, folds, lambda, 
                                                    list_PSD_lasso, list_sigma_lasso, list_PSD_error, list_sigma_error,
                                                    ratio_matrix=NULL, beta1.start, beta2.start, noise, penalty=penalty) {
  
  ### Calculating the error for the design matrix without the kth fold
  start = p1 + 1
  #Solving the lasso problem without the kth fold
  sigma_corrupted_train <- list_PSD_lasso[[k]]
  sigma_uncorrupted_train <- list_sigma_lasso[[k]] 
  index <- which(folds==k,arr.ind = TRUE)
  X1_cv_train <- Z[-index,1:p1]
  Z2_cv_train <- Z[-index,start:p]
  y_cv_train <- y[-index]
  out = lasso_covariance_block(n=n_without_fold,p1=p1,p2=p2,X1=X1_cv_train,Z2=Z2_cv_train,y=y_cv_train,
                               sigma1=sigma_uncorrupted_train,sigma2=sigma_corrupted_train,lambda=lambda,
                               noise=noise,ratio_matrix = ratio_matrix,beta1.start = beta1.start, beta2.start = beta2.start,
                               penalty=penalty)
  beta1.lambda <- out$coefficients.beta1
  beta2.lambda <- out$coefficients.beta2
  
  #Calculating the error on the remaining fold
  sigma_corrupted_test <- list_PSD_error[[k]]
  sigma_uncorrupted_test <- list_sigma_error[[k]]
  X1_cv_test <- Z[index,1:p1]
  Z2_cv_test <- Z[index,start:p]
  y_cv_test <- y[index]
  
  rho_1 <- 1/n_one_fold * t(X1_cv_test) %*% y_cv_test
  if (noise == "additive"){
    rho_2 <- 1/ n_one_fold * t(Z2_cv_test) %*% y_cv_test 
    sigma3 <- 1/ n_one_fold * t(Z2_cv_test) %*% X1_cv_test %*% beta1.lambda
  }else{
    Z2_tilde <- sapply(1:p2,function(j)Z2_cv_test[,j] / diag(ratio_matrix)[j])
    rho_2 <- 1/ n_one_fold * (t(Z2_tilde) %*% y_cv_test) 
    sigma3 <- 1/ n_one_fold * (t(Z2_tilde) %*% X1_cv_test %*% beta1.lambda) 
  }
  error <- t(beta1.lambda)%*%sigma_uncorrupted_test%*%beta1.lambda + t(beta2.lambda)%*%sigma_corrupted_test%*%beta2.lambda - 2*t(rho_1)%*%beta1.lambda - 2*t(rho_2)%*%beta2.lambda + 2*t(beta2.lambda)%*%sigma3
  
  error
}

blockwise_coordinate_descent <- function(Z, y, n, p, p1, p2, center.Z=TRUE, scale.Z=TRUE, center.y=TRUE, scale.y=TRUE,
                                         lambda.factor=ifelse(dim(Z)[1] < dim(Z)[2], 0.01, 0.001), step=100, K=4, mu=10, tau=NULL, etol=1e-4, optTol=1e-5,
                                         earlyStopping_max=10, noise=c("additive", "missing"), penalty=c("lasso", "SCAD"), mode="ADMM") {
  
  nrows = nrow(Z)
  ncols = ncol(Z)
  vnames = colnames(Z)
  
  if(!(is.matrix(Z))){
    stop("Z has to be a matrix")
  }
  if(!(is.matrix(y))){
    stop("y has to be a matrix")
  }
  if(!is.null(tau) & !is.numeric(tau)){
    stop("tau must be numeric")
  }
  if(n != nrows){
    stop(paste("Number of rows in Z (", nrows, ") different from n(", n, ")"),sep="")
  }
  if(p != ncols){
    stop(paste("Number of columns in Z (", ncols, ") different from p (", p, ")"),sep="")
  }
  if (nrows != dim(y)[1]){
    stop(paste("Number of rows in Z (", nrows, ") different from number of rows in y (", dim(y)[1], ") "),sep="")
  }
  if (!is.numeric(y)) {
    stop("The response y must be numeric. Factors must be converted to numeric")
  }
  if(lambda.factor >= 1){
    stop("lambda factor should be smaller than 1")
  }
  if(p1 + p2 != p){
    stop(paste("Sum of p1 and p2 (", p1 + p2, ") should be equal to p (", p, ")"),sep="")
  }
  if(mu>500 || mu<1){
    warning(paste("Mu value (", mu, ") is not in the usual range (10-500)"))
  }
  if (n %% K != 0){
    stop("K should be a divider of n")
  }
  if(noise=="missing" && center.Z == FALSE){
    stop("When noise is equal to missing, it is required to center matrix Z. Use center.Z=TRUE.")
  }
  if(scale.Z == FALSE && noise=='missing'){
    warning("When noise is equal to missing, it is recommended to use scale.Z equal to TRUE in order 
            to obtain trustworthy results.")
  }
  if(scale.Z == TRUE && noise=='additive'){
    warning("When noise is equal to additive, it is recommended to use scale.Z equal to FALSE in order 
            to obtain trustworthy results. Otherwise, the scaling should be taken into account
            when introducing the error parameter as a function parameter.")
  }
  
  #General variables we are going to use in the function
  start = 1 + p1
  n_without_fold = n - floor(n/K)
  n_one_fold = floor(n/K)
  earlyStopping = step
  
  ratio_matrix = NULL
  if (noise=="missing"){
    ratio_matrix = matrix(0,p2,p2)
    
    for (i in 1:p2){
      for (j in i:p2){
        n_ij = length(intersect(which(!is.na(Z[,p1 + i])),which(!is.na(Z[,p1 + j]))))
        ratio_matrix[i,j] = n_ij
        ratio_matrix[j,i] = n_ij
      }
    }
    ratio_matrix = ratio_matrix/n
  }
  
  mean.Z = sapply(1:p, function(j)mean_without_NA(j,Z))
  sd.Z <- sapply(1:p, function(j)sd_without_NA_block(j,Z))
  
  
  
  if (center.Z == TRUE){
    if (scale.Z == TRUE){
      Z[,start:p] = sapply(1:p2, function(j)rescale_without_NA_block(j,Z,p1))
      Z[,start:p] = sapply(1:p2, function(j)change_NA_value_block(j,Z,p1))
      Z[,1:p1] = scale(Z[,1:p1], center = TRUE, scale = FALSE)
      Z = sapply(1:p, function(j)scale_manual_with_sd(j,Z,sd.Z))
    }else{
      Z[,start:p] = sapply(1:p2, function(j)rescale_without_NA_block(j,Z,p1))
      Z[,start:p] = sapply(1:p2, function(j)change_NA_value_block(j,Z,p1))
      Z[,1:p1] = scale(Z[,1:p1], center = TRUE, scale = FALSE)
    }
  }else{
    if (scale.Z == TRUE){
      Z = sapply(1:p, function(j)scale_manual_with_sd(j,Z,sd.Z))
    }
  }
  
  mean.y = mean(y)
  sd.y = stats::sd(y)
  
  if (center.y == TRUE){
    if (scale.y == TRUE){
      y = scale(y, center = TRUE, scale = TRUE)
    }else{
      y = scale(y, center = TRUE, scale = FALSE)
    }
  }
  
  lambda_max <- lambda_max.coordinate_descent(Z=Z,y=y,n=n,p=p,p1=p1,p2=p2,ratio_matrix=ratio_matrix,noise=noise)
  lambda_min <- lambda.factor*lambda_max
  lambda_list <- emdbook::lseq(lambda_max,lambda_min,step)
  beta1.start <- rep(0,p1)
  beta2.start <- rep(0,p2)
  beta.start <- c(beta1.start,beta2.start)
  best.lambda <- lambda_max
  beta.opt <- beta.start
  best.error <- 10000
  error_list <- matrix(0,step,4)
  error <- 0
  earlyStopping_high = 0
  
  matrix_beta <- matrix(0,step,p)
  
  ### Creating the K matrices we are going to use for cross validation
  output = cv_covariance_matrices_block_descent(K=K, mat=Z, y=y, p=p, p1=p1, p2=p2, mu=mu, 
                                                tau=tau, ratio_matrix = ratio_matrix, etol=etol, noise = noise, mode=mode)
  list_PSD_lasso = output$list_PSD_lasso
  list_PSD_error = output$list_PSD_error
  list_sigma_lasso = output$list_sigma_lasso
  list_sigma_error = output$list_sigma_error
  sigma1 = output$sigma_global_uncorrupted
  sigma2 = output$sigma_global_corrupted
  folds = output$folds
  X1 = Z[,1:p1]
  Z2 = Z[,start:p]
  
  for (i in 1:step){
    
    lambda_step <- lambda_list[i]
    error_old <- error
    
    out = sapply(1:K, function(k)cross_validation_function.block_descent(k,
                                                                         Z,
                                                                         y,
                                                                         n,
                                                                         n_one_fold,
                                                                         n_without_fold,
                                                                         p,
                                                                         p1,
                                                                         p2,
                                                                         folds,
                                                                         lambda_step,
                                                                         list_PSD_lasso,
                                                                         list_sigma_lasso,
                                                                         list_PSD_error,
                                                                         list_sigma_error,
                                                                         ratio_matrix = ratio_matrix,
                                                                         beta1.start,
                                                                         beta2.start,
                                                                         noise=noise,
                                                                         penalty=penalty))
    error = mean(out)
    sd_low = stats::quantile(out, probs = c(0.1))
    sd_high = stats::quantile(out, probs = c(0.9))
    error_list[i,1] <- error
    error_list[i,2] <- sd_low
    error_list[i,3] <- sd_high
    error_list[i,4] <- stats::sd(out)
    out = lasso_covariance_block(n=n, p1=p1, p2=p2, X1=X1, Z2=Z2, y=y, sigma1=sigma1, sigma2=sigma2, lambda=lambda_step, 
                                 noise=noise, ratio_matrix = ratio_matrix, beta1.start = beta1.start, beta2.start = beta2.start,
                                 penalty=penalty)
    beta1 <- out$coefficients.beta1
    beta2 <- out$coefficients.beta2
    beta <- c(beta1,beta2)
    
    beta1.start <- beta1
    beta2.start <- beta2
    matrix_beta[i,] <- beta
    
    ### Checking for optimal parameters
    if (error <= best.error){
      best.error <- error
      best.lambda <- lambda_step
      beta.opt <- beta
    }
    
    ## Early stopping
    if (abs(error - error_old) < optTol){
      print("Early Stopping because of convergence of the error")
      earlyStopping = i
      break
    }
    
    # Trying to avoid getting to too high lambda values, leading to too time consuming steps
    if (i>= step/2 && error >= error_list[1]){
      print("Value of lambda yielding too high error : we exit the loop")
      earlyStopping = i
      break
    }
    if (error > best.error){
      earlyStopping_high = earlyStopping_high +1
      if (earlyStopping_high >= earlyStopping_max){
        print("Early stopping because of error getting too high")
        earlyStopping = i
        break
      }
    }
  }
  df <- data.frame(lambda=lambda_list[1:earlyStopping], error=error_list[1:earlyStopping,1],error.inf=error_list[1:earlyStopping,2],error.sup=error_list[1:earlyStopping,3],error.sd=error_list[1:earlyStopping,4])
  step.min <- which(df[,"error"] == best.error)
  sd.best <- df[step.min,"error.sd"]
  step.sd <- max(which(df[,"error"] > best.error + sd.best & df[,"lambda"] > df[step.min,"lambda"]))
  lambda.sd <- df[step.sd,"lambda"]
  beta.sd <- matrix_beta[step.sd,]
  
  
  data_intermediate <- data.frame(matrix_beta[1:earlyStopping,])
  names(data_intermediate) <- sapply(1:p, function(i)paste0("beta",i))
  
  data_beta <- data.frame(lambda = lambda_list[1:earlyStopping])
  
  data_beta <- cbind(data_beta,data_intermediate)
  
  fit = list(
    lambda.opt = best.lambda,
    lambda.sd = lambda.sd,
    beta.opt = beta.opt,
    beta.sd = beta.sd,
    data_error = df,
    data_beta = data_beta,
    earlyStopping = earlyStopping,
    vnames = vnames,
    mean.Z = mean.Z,
    sd.Z = sd.Z,
    mean.y = mean.y,
    sd.y = sd.y                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
  )
  return(fit)
}

lambda_max.coordinate_descent.general <- function(Z, y, n, p, p1, p2, p3, ratio_matrix=NULL) {
  start <- p1 + p2 + 1
  X1 <- Z[,1:(p1+p2)]
  Z2 <- Z[,start:p]
  rho_tilde1 <- 1/n*t(X1)%*%y
  rho_tilde2 <- 1/n*t(Z2)%*%y/ diag(ratio_matrix)
  lambda = max(max(abs(rho_tilde1)),max(abs(rho_tilde2)))
  lambda
}

scale_manual <- function(j, Z) {
  sd <- stats::sd(Z[,j])
  if (sd != 0){
    return(Z[,j]/sd)
  }else{
    return (Z[,j])
  }
  
}

rescale_without_NA_block <- function(j, Z, p1) {
  m <- mean(Z[which(!is.na(Z[,p1 + j]), arr.ind = TRUE),p1+j])
  Z[,p1 + j] - m
}

mean_without_NA <- function(j, Z) {
  m <- mean(Z[which(!is.na(Z[,j]), arr.ind = TRUE),j])
  m
}

sd_without_NA_block <- function(j, Z) {
  sd <- stats::sd(Z[which(!is.na(Z[,j]), arr.ind = TRUE),j])
  sd
}

change_NA_value_block <- function(j, Z, p1) {
  Z[which(is.na(Z[,p1 + j]), arr.ind = TRUE),j + p1] <- 0
  Z[,p1 + j]
}

scale_manual_with_sd <- function(j, Z, v) {
  sd <- v[j]
  if (sd != 0){
    return(Z[,j]/sd)
  }else{
    return (Z[,j])
  }
}

cross_validation_function.block_descent_general <- function(k, Z, y, n, n_one_fold, n_without_fold, p, p1, p2, p3, folds, lambda,
                                                            list_PSD_lasso_additive, list_PSD_lasso_missing, list_sigma_lasso, list_PSD_error_additive, list_PSD_error_missing,
                                                            list_sigma_error, ratio_matrix=NULL, beta1.start, beta2.start, beta3.start, penalty=penalty) {
  
  ### Calculating the error for the design matrix without the kth fold
  start = p1 + 1
  #Solving the lasso problem without the kth fold
  if (p1 != 0) {
    sigma_uncorrupted_train <- list_sigma_lasso[[k]] 
  }
  sigma_corrupted_train_additive <- list_PSD_lasso_additive[[k]]
  sigma_corrupted_train_missing <- list_PSD_lasso_missing[[k]]
  
  index <- which(folds==k,arr.ind = TRUE)
  X1_cv_train <- Z[-index,1:p1]
  Z2_cv_train <- Z[-index,start:(p1+p2)]
  Z3_cv_train <- Z[-index,(p1+p2+1):p]
  y_cv_train <- y[-index]
  out = lasso_covariance_block_general(n=n_without_fold,p1=p1,p2=p2,p3=p3,X1=X1_cv_train,Z2=Z2_cv_train,Z3=Z3_cv_train,y=y_cv_train,
                                       sigma1=sigma_uncorrupted_train,sigma2=sigma_corrupted_train_additive,sigma3=sigma_corrupted_train_missing,lambda=lambda,
                                       ratio_matrix = ratio_matrix,beta1.start = beta1.start, beta2.start = beta2.start, beta3.start = beta3.start,
                                       penalty=penalty)
  if (p1 != 0) {
    beta1.lambda <- out$coefficients.beta1
  }
  beta2.lambda <- out$coefficients.beta2
  beta3.lambda <- out$coefficients.beta3
  
  #Calculating the error on the remaining fold
  if (p1 != 0) {
    sigma_uncorrupted_test <- list_sigma_error[[k]]
    sigma_corrupted_test_additive <- list_PSD_error_additive[[k]]
    sigma_corrupted_test_missing <- list_PSD_error_missing[[k]]
    X1_cv_test <- Z[index,1:p1]
    Z2_cv_test <- Z[index,start:(p1+p2)]
    Z3_cv_test <- Z[index,(p1+p2+1):p]
    y_cv_test <- y[index]
    
    rho_1 <- 1/n_one_fold * t(X1_cv_test) %*% y_cv_test
    rho_2 <- 1/n_one_fold * t(Z2_cv_test) %*% y_cv_test
    Z3_tilde <- sapply(1:p3,function(j)Z3_cv_test[,j] / diag(ratio_matrix)[j])
    rho_3 <- 1/n_one_fold * t(Z3_tilde) %*% y_cv_test
    
    sigma21 <- 1/n_one_fold * t(Z2_cv_test) %*% X1_cv_test %*% beta1.lambda
    sigma31 <- 1/n_one_fold * t(Z3_tilde) %*% X1_cv_test %*% beta1.lambda
    sigma32 <- 1/n_one_fold * t(Z3_tilde) %*% Z2_cv_test %*% beta2.lambda
    
    error <- t(beta1.lambda)%*%sigma_uncorrupted_test%*%beta1.lambda + 
      t(beta2.lambda)%*%sigma_corrupted_test_additive%*%beta2.lambda +
      t(beta3.lambda)%*%sigma_corrupted_test_missing%*%beta3.lambda - 
      2*t(rho_1)%*%beta1.lambda - 
      2*t(rho_2)%*%beta2.lambda -
      2*t(rho_3)%*%beta3.lambda + 
      2*t(beta2.lambda)%*%sigma21 +
      2*t(beta3.lambda)%*%sigma31 +
      2*t(beta3.lambda)%*%sigma32
    
    return(error)
  }
  if (p1 == 0) {
    sigma_corrupted_test_additive <- list_PSD_error_additive[[k]]
    sigma_corrupted_test_missing <- list_PSD_error_missing[[k]]
    Z2_cv_test <- Z[index,start:(p1+p2)]
    Z3_cv_test <- Z[index,(p1+p2+1):p]
    y_cv_test <- y[index]
    
    rho_2 <- 1/n_one_fold * t(Z2_cv_test) %*% y_cv_test
    Z3_tilde <- sapply(1:p3,function(j)Z3_cv_test[,j] / diag(ratio_matrix)[j])
    rho_3 <- 1/n_one_fold * t(Z3_tilde) %*% y_cv_test
    
    sigma32 <- 1/n_one_fold * t(Z3_tilde) %*% Z2_cv_test %*% beta2.lambda
    
    error <- t(beta2.lambda)%*%sigma_corrupted_test_additive%*%beta2.lambda +
      t(beta3.lambda)%*%sigma_corrupted_test_missing%*%beta3.lambda - 
      2*t(rho_2)%*%beta2.lambda -
      2*t(rho_3)%*%beta3.lambda + 
      2*t(beta3.lambda)%*%sigma32
    
    return(error)
  }
}

blockwise_coordinate_descent_general <- function(Z, y, n, p, p1, p2, p3,center.Z=TRUE, scale.Z=TRUE, center.y=TRUE, scale.y=TRUE,
                                                 lambda.factor=ifelse(dim(Z)[1] < dim(Z)[2], 0.01, 0.001), step=100, K=4, mu=10, tau=NULL, etol=1e-4,
                                                 optTol=1e-5, earlyStopping_max=10, penalty=c("lasso", "SCAD"), mode="ADMM") {
  
  nrows = nrow(Z)
  ncols = ncol(Z)
  vnames = colnames(Z)
  
  if(!(is.matrix(Z))){
    stop("Z has to be a matrix")
  }
  if(!(is.matrix(y))){
    stop("y has to be a matrix")
  }
  if(!is.null(tau) & !is.numeric(tau)){
    stop("tau must be numeric")
  }
  if(n != nrows){
    stop(paste("Number of rows in Z (", nrows, ") different from n(", n, ")"),sep="")
  }
  if(p != ncols){
    stop(paste("Number of columns in Z (", ncols, ") different from p (", p, ")"),sep="")
  }
  if (nrows != dim(y)[1]){
    stop(paste("Number of rows in Z (", nrows, ") different from number of rows in y (", dim(y)[1], ") "),sep="")
  }
  if (!is.numeric(y)) {
    stop("The response y must be numeric. Factors must be converted to numeric")
  }
  if(lambda.factor >= 1){
    stop("lambda factor should be smaller than 1")
  }
  if(p1 + p2 + p3 != p){
    stop(paste("Sum of p1, p2 and p3 (", p1 + p2 + p3, ") should be equal to p (", p, ")"),sep="")
  }
  if(mu>500 || mu<1){
    warning(paste("Mu value (", mu, ") is not in the usual range (10-500)"))
  }
  if (n %% K != 0){
    stop("K should be a divider of n")
  }
  if(center.Z == FALSE){
    stop("When the data contain missingness, it is required to center matrix Z. Use center.Z=TRUE.")
  }
  if(scale.Z == FALSE){
    warning("When the data contain missingness, it is recommended to use scale.Z equal to TRUE in order 
            to obtain trustworthy results.")
  }
  
  #General variables we are going to use in the function
  start = 1 + p1 + p2
  n_without_fold = n - floor(n/K)
  n_one_fold = floor(n/K)
  earlyStopping = step
  
  ratio_matrix = NULL
  ratio_matrix = matrix(0,p3,p3)
  
  for (i in 1:p3){
    for (j in i:p3){
      n_ij = length(intersect(which(!is.na(Z[,p1 + p2 + i])),which(!is.na(Z[,p1 + p2 + j]))))
      ratio_matrix[i,j] = n_ij
      ratio_matrix[j,i] = n_ij
    }
  }
  ratio_matrix = ratio_matrix/n
  
  
  mean.Z = sapply(1:p, function(j)mean_without_NA(j,Z))
  sd.Z <- sapply(1:p, function(j)sd_without_NA_block(j,Z))
  
  
  
  if (center.Z == TRUE){
    if (scale.Z == TRUE){
      Z[,start:p] = sapply(1:p3, function(j)rescale_without_NA_block(j,Z,p1+p2))
      Z[,start:p] = sapply(1:p3, function(j)change_NA_value_block(j,Z,p1+p2))
      Z[,1:(p1+p2)] = scale(Z[,1:(p1+p2)], center = TRUE, scale = FALSE)
      if (p1 == 0) {
        Z[,c((p1+p2+1):p)] = sapply(1:(p - p2), function(j)scale_manual_with_sd(j,Z[,c((p1+p2+1):p)],sd.Z[c((p1+p2+1):p)]))
      }
      if (p1 != 0) {
        Z[,c(1:p1,(p1+p2+1):p)] = sapply(1:(p - p2), function(j)scale_manual_with_sd(j,Z[,c(1:p1,(p1+p2+1):p)],sd.Z[c(1:p1,(p1+p2+1):p)]))
      }
    }else{
      Z[,start:p] = sapply(1:p3, function(j)rescale_without_NA_block(j,Z,p1+p2))
      Z[,start:p] = sapply(1:p3, function(j)change_NA_value_block(j,Z,p1+p2))
      Z[,1:(p1+p2)] = scale(Z[,1:(p1+p2)], center = TRUE, scale = FALSE)
    }
  }else{
    if (scale.Z == TRUE){
      Z = sapply(1:p, function(j)scale_manual_with_sd(j,Z,sd.Z))
    }
  }
  
  mean.y = mean(y)
  sd.y = stats::sd(y)
  
  if (center.y == TRUE){
    if (scale.y == TRUE){
      y = scale(y, center = TRUE, scale = TRUE)
    }else{
      y = scale(y, center = TRUE, scale = FALSE)
    }
  }
  
  lambda_max <- lambda_max.coordinate_descent.general(Z=Z,y=y,n=n,p=p,p1=p1,p2=p2,p3=p3,ratio_matrix=ratio_matrix)
  lambda_min <- lambda.factor*lambda_max
  lambda_list <- emdbook::lseq(lambda_max,lambda_min,step)
  beta1.start <- rep(0,p1)
  beta2.start <- rep(0,p2)
  beta3.start <- rep(0,p3)
  beta.start <- c(beta1.start,beta2.start,beta3.start)
  best.lambda <- lambda_max
  beta.opt <- beta.start
  best.error <- 10000
  error_list <- matrix(0,step,4)
  error <- 0
  earlyStopping_high = 0
  
  matrix_beta <- matrix(0,step,p)
  
  ### Creating the K matrices we are going to use for cross validation
  output = cv_covariance_matrices_block_descent_general(K=K, mat=Z, y=y, p=p, p1=p1, p2=p2, p3=p3, mu=mu, 
                                                        tau=tau, ratio_matrix = ratio_matrix, etol=etol, mode=mode)
  list_PSD_lasso_additive = output$list_PSD_lasso_additive
  list_PSD_error_additive = output$list_PSD_error_additive
  list_PSD_lasso_missing = output$list_PSD_lasso_missing
  list_PSD_error_missing = output$list_PSD_error_missing
  if (p1!=0) {
    list_sigma_lasso = output$list_sigma_lasso
    list_sigma_error = output$list_sigma_error
    sigma1 = output$sigma_global_uncorrupted
  }
  if (p1==0) {
    list_sigma_lasso = NULL
    list_sigma_error = NULL
    sigma1 = NULL
  }
  sigma2 = output$sigma_global_corrupted_additive
  sigma3 = output$sigma_global_corrupted_missing
  folds = output$folds
  X1 = Z[,1:p1]
  Z2 = Z[,(p1+1):(p1+p2)]
  Z3 = Z[,(p1+p2+1):p]
  
  for (i in 1:step){
    print(paste0("iteration: ",i))
    
    lambda_step <- lambda_list[i]
    error_old <- error
    
    out = sapply(1:K, function(k)cross_validation_function.block_descent_general(k,
                                                                                 Z,
                                                                                 y,
                                                                                 n,
                                                                                 n_one_fold,
                                                                                 n_without_fold,
                                                                                 p,
                                                                                 p1,
                                                                                 p2,
                                                                                 p3,
                                                                                 folds,
                                                                                 lambda_step,
                                                                                 list_PSD_lasso_additive,
                                                                                 list_PSD_lasso_missing,
                                                                                 list_sigma_lasso,
                                                                                 list_PSD_error_additive,
                                                                                 list_PSD_error_missing,
                                                                                 list_sigma_error,
                                                                                 ratio_matrix = ratio_matrix,
                                                                                 beta1.start,
                                                                                 beta2.start,
                                                                                 beta3.start,
                                                                                 penalty=penalty))
    error = mean(out)
    sd_low = stats::quantile(out, probs = c(0.1))
    sd_high = stats::quantile(out, probs = c(0.9))
    error_list[i,1] <- error
    error_list[i,2] <- sd_low
    error_list[i,3] <- sd_high
    error_list[i,4] <- stats::sd(out)
    out = lasso_covariance_block_general(n=n, p1=p1, p2=p2, p3=p3, X1=X1, Z2=Z2, Z3=Z3, y=y, sigma1=sigma1, sigma2=sigma2, sigma3=sigma3, lambda=lambda_step, 
                                         ratio_matrix = ratio_matrix, beta1.start = beta1.start, beta2.start = beta2.start, beta3.start = beta3.start,
                                         penalty=penalty)
    if (p1 != 0) {
      beta1 <- out$coefficients.beta1
    }
    beta2 <- out$coefficients.beta2
    beta3 <- out$coefficients.beta3
    if (p1 != 0) {
      beta <- c(beta1,beta2,beta3)
    }
    if (p1 == 0) {
      beta <- c(beta2,beta3)
    }
    
    if (p1 != 0) {
      beta1.start <- beta1
    }
    if (p1 == 0) {
      beta1.start <- rep(0,p1)
    }
    beta2.start <- beta2
    beta3.start <- beta3
    matrix_beta[i,] <- beta
    
    ### Checking for optimal parameters
    if (error <= best.error){
      best.error <- error
      best.lambda <- lambda_step
      beta.opt <- beta
    }
    
    ## Early stopping
    if (abs(error - error_old) < optTol){
      print("Early Stopping because of convergence of the error")
      earlyStopping = i
      break
    }
    
    # Trying to avoid getting to too high lambda values, leading to too time consuming steps
    if (i>= step/2 && error >= error_list[1]){
      print("Value of lambda yielding too high error : we exit the loop")
      earlyStopping = i
      break
    }
    if (error > best.error){
      earlyStopping_high = earlyStopping_high +1
      if (earlyStopping_high >= earlyStopping_max){
        print("Early stopping because of error getting too high")
        earlyStopping = i
        break
      }
    }
  }
  df <- data.frame(lambda=lambda_list[1:earlyStopping], error=error_list[1:earlyStopping,1],error.inf=error_list[1:earlyStopping,2],error.sup=error_list[1:earlyStopping,3],error.sd=error_list[1:earlyStopping,4])
  step.min <- which(df[,"error"] == best.error)
  sd.best <- df[step.min,"error.sd"]
  step.sd <- max(which(df[,"error"] > best.error + sd.best & df[,"lambda"] > df[step.min,"lambda"]))
  lambda.sd <- df[step.sd,"lambda"]
  beta.sd <- matrix_beta[step.sd,]
  
  
  data_intermediate <- data.frame(matrix_beta[1:earlyStopping,])
  names(data_intermediate) <- sapply(1:p, function(i)paste0("beta",i))
  
  data_beta <- data.frame(lambda = lambda_list[1:earlyStopping])
  
  data_beta <- cbind(data_beta,data_intermediate)
  
  fit = list(
    lambda.opt = best.lambda,
    lambda.sd = lambda.sd,
    beta.opt = beta.opt,
    beta.sd = beta.sd,
    data_error = df,
    data_beta = data_beta,
    earlyStopping = earlyStopping,
    vnames = vnames,
    mean.Z = mean.Z,
    sd.Z = sd.Z,
    mean.y = mean.y,
    sd.y = sd.y                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
  )
  return(fit)
}

coco <- function(Z, y, n, p, p1=NULL, p2=NULL, center.Z=TRUE, scale.Z=TRUE, center.y=TRUE, scale.y=TRUE, 
                 lambda.factor=ifelse(dim(Z)[1] < dim(Z)[2], 0.01, 0.001), step=100, K=4, mu=10, tau=NULL, etol=1e-4, optTol=1e-5,
                 earlyStopping_max=10, noise=c("additive", "missing"), block=TRUE, penalty=c("lasso", "SCAD"), mode="ADMM") {
  
  this.call <- match.call()
  if(block){
    fit <- blockwise_coordinate_descent(Z=Z,
                                        y=y,
                                        n=n,
                                        p=p,
                                        p1=p1,
                                        p2=p2,
                                        center.Z = center.Z,
                                        scale.Z = scale.Z,
                                        center.y = center.y,
                                        scale.y = scale.y,
                                        lambda.factor = lambda.factor,
                                        step = step,
                                        K = K,
                                        mu = mu,
                                        tau = tau,
                                        etol = etol,
                                        optTol = optTol,
                                        earlyStopping_max = earlyStopping_max,
                                        noise = noise,
                                        penalty=penalty,
                                        mode=mode)
  }else{
    fit <- pathwise_coordinate_descent(Z=Z,
                                       y=y,
                                       n=n,
                                       p=p,
                                       center.Z = center.Z,
                                       scale.Z = scale.Z,
                                       center.y = center.y,
                                       scale.y = scale.y,
                                       lambda.factor = lambda.factor,
                                       step = step,
                                       K = K,
                                       mu = mu,
                                       tau = tau,
                                       etol = etol,
                                       optTol = optTol,
                                       earlyStopping_max = earlyStopping_max,
                                       noise = noise,
                                       penalty=penalty,
                                       mode=mode)
  }
  fit$call <- this.call
  class(fit) <- "coco"
  return(fit)
}

cv_covariance_matrices <- function(K, mat, y, p, mu, tau=NULL, ratio_matrix=NULL,
                                   etol=1e-4, noise=c("additive", "missing"), mode="ADMM") {
  
  
  # calculate the K nearest PSD covariance matrices 
  #      in the cross validation process
  
  # @param tau is the sd of the covariance matrix of the additive error. We study the simple case where 
  # there is no correlation between the error for different features.
  # @param ratio_matrix is the observation matrix. Its j,k term calculates the number of rows where we observe j and
  # k feature at the same time if j!=k, and the number of rows where we observe the j feature if j==k
  
  n = nrow(mat)
  p = ncol(mat)
  n_without_fold = n - floor(n/K)
  n_one_fold = floor(n/K)
  
  folds = sample(cut(seq(1,n),breaks=K,labels=FALSE))
  list_matrices_lasso <- list()
  list_matrices_error <- list()
  list_rho_lasso <- list()
  list_rho_error <- list()
  
  ### Case where there is additive noise
  if (noise == "additive"){
    
    #We calculate the global nearest PSD cov matrix and the global surrogate rho, when we take into account the whole data set
    
    print("Doing the global data")
    cov_modified <- 1/n*t(mat)%*%mat - tau**2*diag(p)
    if (mode=="ADMM") {
      sigma_global <- ADMM_proj(cov_modified,mu=mu, etol = etol)$mat
    }
    if (mode=="HM") {
      sigma_global <- HM_proj(sigmaHat = cov_modified,R=ratio_matrix,mu=mu, tolerance = etol)
    }
    rho_global <- t(mat)%*%y/n
    
    for (i in 1:K){
      # We calculate the necessary matrices for the cross validation
      
      print(paste("Doing the",i,"fold"))
      index <- which(folds==i, arr.ind= TRUE)
      
      #Calculating the nearest PSD cov matrix when we remove the kth fold, to resolve lasso problem during cross validation
      mat_train <- mat[-index,]
      cov_modified_train <- 1/n_without_fold*t(mat_train)%*%mat_train - tau**2*diag(p)
      if (mode=="ADMM") {
        mat_cov_train <- ADMM_proj(cov_modified_train,mu=mu, etol = etol)$mat
      }
      if (mode=="HM") {
        mat_cov_train <- HM_proj(sigmaHat = cov_modified_train,R=ratio_matrix,mu=mu, tolerance = etol)
      }
      list_matrices_lasso <- rlist::list.append(list_matrices_lasso,mat_cov_train)
      
      #Calculating the nearest PSD cov matrix for the kth fold, to calculate the error on the problem solved without the kth fold
      mat_test <- mat[index,]
      cov_modified_test <- 1/n_one_fold*t(mat_test)%*%mat_test - tau**2*diag(p)
      if (mode=="ADMM") {
        mat_cov_test <- ADMM_proj(cov_modified_test,mu=mu, etol = etol)$mat
      }
      if (mode=="HM") {
        mat_cov_test <- HM_proj(sigmaHat = cov_modified_test,R=ratio_matrix,mu=mu, tolerance = etol)
      }
      list_matrices_error <- rlist::list.append(list_matrices_error,mat_cov_test)
      
      #Calculating the surrogate rho when we remove the kth fold, to resolve lasso problem during cross validation
      y_train <- y[-index,]
      rho_train <- 1/n_without_fold * t(mat_train)%*%y_train
      list_rho_lasso <- rlist::list.append(list_rho_lasso,rho_train)
      
      #Calculating the surrogate rho for the kth fold, to calculate the error on the problem solved without the kth fold
      y_test <- y[index,]
      rho_test <- 1/n_one_fold * t(mat_test)%*%y_test
      list_rho_error <- rlist::list.append(list_rho_error,rho_test)
      
    }
  }
  
  else if (noise == "missing"){
    #We calculate the global nearest PSD cov matrix and the global surrogate rho, when we take into account the whole data set
    #mat_for_adjustment <- diag(probs*(1-probs),p,p) + (1-probs) %*% t(1 - probs)
    print("Doing the global data")
    cov_modified <- 1/n*t(mat)%*%mat / ratio_matrix
    if (mode=="ADMM") {
      sigma_global <- ADMM_proj(cov_modified,mu=mu, etol = etol)$mat
    }
    if (mode=="HM") {
      sigma_global <- HM_proj(sigmaHat = cov_modified,R=ratio_matrix,mu=mu, tolerance = etol)
    }
    rho_global <- 1/n*t(mat)%*%y/diag(ratio_matrix)
    
    for (i in 1:K){
      # We calculate the necessary matrices for the cross validation
      
      print(paste("Doing the",i,"fold"))
      index <- which(folds==i, arr.ind= TRUE)
      
      #Calculating the nearest PSD cov matrix when we remove the kth fold, to resolve lasso problem during cross validation
      mat_train <- mat[-index,]
      # cov_modified_train <- 1/n_without_fold*t(train_mat)%*%train_mat / mat_for_adjustment
      cov_modified_train <- 1/n_without_fold*t(mat_train)%*%mat_train / ratio_matrix
      if (mode=="ADMM") {
        mat_cov_train <- ADMM_proj(cov_modified_train,mu=mu, etol = etol)$mat
      }
      if (mode=="HM") {
        mat_cov_train <- HM_proj(sigmaHat = cov_modified_train,R=ratio_matrix,mu=mu, tolerance = etol)
      }
      list_matrices_lasso <- rlist::list.append(list_matrices_lasso,mat_cov_train)
      
      #Calculating the nearest PSD cov matrix for the kth fold, to calculate the error on the problem solved without the kth fold
      mat_test <- mat[index,]
      # cov_modified_test <- 1/n_one_fold*t(test_mat)%*%test_mat / mat_for_adjustment
      cov_modified_test <- 1/n_one_fold*t(mat_test)%*%mat_test / ratio_matrix
      if (mode=="ADMM") {
        mat_cov_test <- ADMM_proj(cov_modified_test,mu=mu, etol = etol)$mat
      }
      if (mode=="HM") {
        mat_cov_test <- HM_proj(sigmaHat = cov_modified_test,R=ratio_matrix,mu=mu, tolerance = etol)
      }
      list_matrices_error <- rlist::list.append(list_matrices_error,mat_cov_test)
      
      #Calculating the surrogate rho when we remove the kth fold, to resolve lasso problem during cross validation
      y_train <- y[-index,]
      # rho_train <- 1/n_without_fold * t(train_mat)%*%y_train / (1-probs)
      rho_train <- 1/n_without_fold * t(mat_train)%*%y_train / diag(ratio_matrix)
      list_rho_lasso <- rlist::list.append(list_rho_lasso,rho_train)
      
      #Calculating the surrogate rho for the kth fold, to calculate the error on the problem solved without the kth fold
      y_test <- y[index,]
      # rho_test <- 1/n_one_fold * t(test_mat)%*%y_test / (1-probs)
      rho_test <- 1/n_one_fold * t(mat_test)%*%y_test / diag(ratio_matrix)
      list_rho_error <- rlist::list.append(list_rho_error,rho_test)
      
    }
  }
  
  
  return(list(sigma_global = sigma_global, rho_global = rho_global, list_matrices_lasso = list_matrices_lasso, list_matrices_error = list_matrices_error, 
              list_rho_lasso=list_rho_lasso, list_rho_error=list_rho_error))
}

cv_covariance_matrices_block_descent <- function(K, mat, y, p, p1, p2, mu, tau=NULL,
                                                 ratio_matrix=NULL, etol=1e-4, noise=c("additive","missing"), mode="ADMM") {
  
  n = nrow(mat)
  p = ncol(mat)
  start = p1 + 1
  n_without_fold = n - floor(n/K)
  n_one_fold = floor(n/K)
  
  folds = sample(cut(seq(1,n),breaks=K,labels=FALSE))
  list_PSD_lasso <- list()
  list_PSD_error <- list()
  list_sigma_lasso <- list()
  list_sigma_error <- list()
  
  ### Case where there is additive noise
  if (noise == "additive"){
    
    #We calculate the global nearest PSD cov matrix and the global surrogate rho, when we take into account the whole data set
    
    print("Processing the global data")
    mat_corrupted <- mat[,start:p]
    mat_uncorrupted <- mat[,1:p1]
    cov_modified <- 1/n*t(mat_corrupted)%*%mat_corrupted - tau**2*diag(p2)
    if (mode=="ADMM") {
      sigma_global_corrupted <- ADMM_proj(cov_modified,mu=mu, etol = etol)$mat
    }
    if (mode=="HM") {
      sigma_global_corrupted <- HM_proj(sigmaHat = cov_modified,R=ratio_matrix,mu=mu, tolerance = etol)
    }
    sigma_global_uncorrupted <- 1/n * t(mat_uncorrupted)%*%mat_uncorrupted
    
    for (i in 1:K){
      # We calculate the necessary matrices for the cross validation
      
      print(paste("Processing the",i,"fold"))
      index <- which(folds==i, arr.ind= TRUE)
      
      #Calculating the nearest PSD cov matrix when we remove the kth fold, to resolve lasso problem during cross validation
      mat_train <- mat[-index,start:p]
      cov_modified_train <- 1/n_without_fold*t(mat_train)%*%mat_train - tau**2*diag(p2)
      if (mode=="ADMM") {
        mat_cov_train <- ADMM_proj(cov_modified_train,mu=mu, etol = etol)$mat
      }
      if (mode=="HM") {
        mat_cov_train <- HM_proj(sigmaHat = cov_modified_train,R=ratio_matrix,mu=mu, tolerance = etol)
      }
      list_PSD_lasso <- rlist::list.append(list_PSD_lasso,mat_cov_train)
      
      #Calculating the nearest PSD cov matrix for the kth fold, to calculate the error on the problem solved without the kth fold
      mat_test <- mat[index,start:p]
      cov_modified_test <- 1/n_one_fold*t(mat_test)%*%mat_test - tau**2*diag(p2)
      if (mode=="ADMM") {
        mat_cov_test <- ADMM_proj(cov_modified_test,mu=mu, etol = etol)$mat
      }
      if (mode=="HM") {
        mat_cov_test <- HM_proj(sigmaHat = cov_modified_test,R=ratio_matrix,mu=mu, tolerance = etol)
      }
      list_PSD_error <- rlist::list.append(list_PSD_error,mat_cov_test)
      
      #Calculating the cov matrix when we remove the kth fold, to resolve lasso problem during cross validation
      mat_train <- mat[-index,1:p1]
      cov_train <- 1/n_without_fold*t(mat_train)%*%mat_train
      list_sigma_lasso <- rlist::list.append(list_sigma_lasso,cov_train)
      
      #Calculating the cov matrix  for the kth fold, to calculate the error on the problem solved without the kth fold
      mat_test <- mat[index,1:p1]
      cov_test <- 1/n_one_fold*t(mat_test)%*%mat_test
      list_sigma_error <- rlist::list.append(list_sigma_error,cov_test)
      
      
      
    }
  }
  
  else if (noise == "missing"){
    #We calculate the global nearest PSD cov matrix and the global surrogate rho, when we take into account the whole data set
    #mat_for_adjustment <- diag(probs*(1-probs),p,p) + (1-probs) %*% t(1 - probs)
    print("Processing the global data")
    mat_corrupted <- mat[,start:p]
    mat_uncorrupted <- mat[,1:p1]
    cov_modified <- 1/n*t(mat_corrupted)%*%mat_corrupted / ratio_matrix
    if (mode=="ADMM") {
      sigma_global_corrupted <- ADMM_proj(cov_modified,mu=mu, etol = etol)$mat
    }
    if (mode=="HM") {
      sigma_global_corrupted <- HM_proj(sigmaHat = cov_modified,R = ratio_matrix,mu=mu, tolerance = etol)
    }
    sigma_global_uncorrupted <- 1/n * t(mat_uncorrupted)%*%mat_uncorrupted
    
    for (i in 1:K){
      # We calculate the necessary matrices for the cross validation
      
      print(paste("Processing the",i,"fold"))
      index <- which(folds==i, arr.ind= TRUE)
      
      #Calculating the nearest PSD cov matrix when we remove the kth fold, to resolve lasso problem during cross validation
      mat_train <- mat[-index,start:p]
      # cov_modified_train <- 1/n_without_fold*t(train_mat)%*%train_mat / mat_for_adjustment
      cov_modified_train <- 1/n_without_fold*t(mat_train)%*%mat_train / ratio_matrix
      if (mode=="ADMM") {
        mat_cov_train <- ADMM_proj(cov_modified_train,mu=mu, etol = etol)$mat
      }
      if (mode=="HM") {
        mat_cov_train <- HM_proj(sigmaHat = cov_modified_train,R=ratio_matrix,mu=mu, tolerance = etol)
      }
      list_PSD_lasso <- rlist::list.append(list_PSD_lasso,mat_cov_train)
      
      #Calculating the nearest PSD cov matrix for the kth fold, to calculate the error on the problem solved without the kth fold
      mat_test <- mat[index,start:p]
      # cov_modified_test <- 1/n_one_fold*t(test_mat)%*%test_mat / mat_for_adjustment
      cov_modified_test <- 1/n_one_fold*t(mat_test)%*%mat_test / ratio_matrix
      if (mode=="ADMM") {
        mat_cov_test <- ADMM_proj(cov_modified_test,mu=mu, etol = etol)$mat
      }
      if (mode=="HM") {
        mat_cov_test <- HM_proj(sigmaHat = cov_modified_test,R=ratio_matrix,mu=mu, tolerance = etol)
      }
      list_PSD_error <- rlist::list.append(list_PSD_error,mat_cov_test)
      
      #Calculating the cov matrix when we remove the kth fold, to resolve lasso problem during cross validation
      mat_train <- mat[-index,1:p1]
      cov_train <- 1/n_without_fold*t(mat_train)%*%mat_train
      list_sigma_lasso <- rlist::list.append(list_sigma_lasso,cov_train)
      
      #Calculating the cov matrix  for the kth fold, to calculate the error on the problem solved without the kth fold
      mat_test <- mat[index,1:p1]
      cov_test <- 1/n_one_fold*t(mat_test)%*%mat_test
      list_sigma_error <- rlist::list.append(list_sigma_error,cov_test)
      
    }
  }
  
  
  return(list(sigma_global_corrupted = sigma_global_corrupted, sigma_global_uncorrupted = sigma_global_uncorrupted, 
              list_PSD_lasso = list_PSD_lasso, list_PSD_error = list_PSD_error, 
              list_sigma_lasso=list_sigma_lasso, list_sigma_error=list_sigma_error, folds = folds))
}

cv_covariance_matrices_block_descent_general <- function(K, mat, y, p, p1, p2, p3, mu,
                                                         tau=NULL, ratio_matrix=NULL, etol=1e-4, mode="ADMM") {
  
  n = nrow(mat)
  p = ncol(mat)
  start = p1 + 1
  n_without_fold = n - floor(n/K)
  n_one_fold = floor(n/K)
  
  folds = sample(cut(seq(1,n),breaks=K,labels=FALSE))
  list_PSD_lasso_additive <- list()
  list_PSD_error_additive <- list()
  list_PSD_lasso_missing <- list()
  list_PSD_error_missing <- list()
  list_sigma_lasso <- list()
  list_sigma_error <- list()
  
  print("Processing the global data")
  if (p1 != 0) {
    mat_uncorrupted <- mat[,1:p1]
    sigma_global_uncorrupted <- 1/n * t(mat_uncorrupted)%*%mat_uncorrupted
  }
  
  mat_corrupted_additive <- mat[,start:(p1+p2)]
  mat_corrupted_missing <- mat[,(p1+p2+1):p]
  cov_modified_additive <- 1/n*t(mat_corrupted_additive)%*%mat_corrupted_additive - tau**2*diag(p2)
  cov_modified_missing <- 1/n*t(mat_corrupted_missing)%*%mat_corrupted_missing / ratio_matrix
  if (mode=="ADMM") {
    sigma_global_corrupted_additive <- ADMM_proj(cov_modified_additive,mu=mu, etol = etol)$mat
    sigma_global_corrupted_missing <- ADMM_proj(cov_modified_missing,mu=mu, etol = etol)$mat
  }
  if (mode=="HM") {
    sigma_global_corrupted_additive <- HM_proj(sigmaHat = cov_modified_additive,mu=mu, tolerance  = etol)
    sigma_global_corrupted_missing <- HM_proj(sigmaHat = cov_modified_missing,R=ratio_matrix,mu=mu, tolerance = etol)
  }
  
  for (i in 1:K){
    # We calculate the necessary matrices for the cross validation
    
    print(paste("Processing the",i,"fold"))
    index <- which(folds==i, arr.ind= TRUE)
    
    #Calculating the nearest PSD cov matrix when we remove the kth fold, to resolve lasso problem during cross validation
    mat_train <- mat[-index,start:(p1+p2)]
    cov_modified_train <- 1/n_without_fold*t(mat_train)%*%mat_train - tau**2*diag(p2)
    if (mode=="ADMM") {
      mat_cov_train <- ADMM_proj(cov_modified_train,mu=mu, etol = etol)$mat
    }
    if (mode=="HM") {
      mat_cov_train <- HM_proj(sigmaHat = cov_modified_train,mu = mu, tolerance = etol)
    }
    list_PSD_lasso_additive <- rlist::list.append(list_PSD_lasso_additive,mat_cov_train)
    
    mat_train <- mat[-index,(p1+p2+1):p]
    cov_modified_train <- 1/n_without_fold*t(mat_train)%*%mat_train / ratio_matrix
    if (mode=="ADMM") {
      mat_cov_train <- ADMM_proj(cov_modified_train,mu=mu, etol = etol)$mat
    }
    if (mode=="HM") {
      mat_cov_train <- HM_proj(sigmaHat = cov_modified_train,R=ratio_matrix,mu = mu, tolerance = etol)
    }
    list_PSD_lasso_missing <- rlist::list.append(list_PSD_lasso_missing,mat_cov_train)
    
    #Calculating the nearest PSD cov matrix for the kth fold, to calculate the error on the problem solved without the kth fold
    mat_test <- mat[index,start:(p1+p2)]
    cov_modified_test <- 1/n_one_fold*t(mat_test)%*%mat_test - tau**2*diag(p2)
    if (mode=="ADMM") {
      mat_cov_test <- ADMM_proj(cov_modified_test,mu=mu, etol = etol)$mat
    }
    if (mode=="HM") {
      mat_cov_test <- HM_proj(sigmaHat = cov_modified_test,mu = mu, tolerance = etol)
    }
    list_PSD_error_additive <- rlist::list.append(list_PSD_error_additive,mat_cov_test)
    
    mat_test <- mat[index,(p1+p2+1):p]
    cov_modified_test <- 1/n_one_fold*t(mat_test)%*%mat_test / ratio_matrix
    mat_cov_test <- ADMM_proj(cov_modified_test,mu=mu, etol = etol)$mat
    if (mode=="ADMM") {
      mat_cov_test <- ADMM_proj(cov_modified_test,mu=mu, etol = etol)$mat
    }
    if (mode=="HM") {
      mat_cov_test <- HM_proj(sigmaHat = cov_modified_test,R=ratio_matrix,mu = mu, tolerance = etol)
    }
    list_PSD_error_missing <- rlist::list.append(list_PSD_error_missing,mat_cov_test)
    
    #Calculating the cov matrix when we remove the kth fold, to resolve lasso problem during cross validation
    mat_train <- mat[-index,1:p1]
    cov_train <- 1/n_without_fold*t(mat_train)%*%mat_train
    list_sigma_lasso <- rlist::list.append(list_sigma_lasso,cov_train)
    
    #Calculating the cov matrix  for the kth fold, to calculate the error on the problem solved without the kth fold
    mat_test <- mat[index,1:p1]
    cov_test <- 1/n_one_fold*t(mat_test)%*%mat_test
    list_sigma_error <- rlist::list.append(list_sigma_error,cov_test)
  }
  if (p1 != 0) {
    return(list(sigma_global_uncorrupted = sigma_global_uncorrupted,
                sigma_global_corrupted_additive = sigma_global_corrupted_additive, 
                sigma_global_corrupted_missing = sigma_global_corrupted_missing,
                list_PSD_lasso_additive = list_PSD_lasso_additive, list_PSD_error_additive = list_PSD_error_additive,
                list_PSD_lasso_missing =list_PSD_lasso_missing, list_PSD_error_missing = list_PSD_error_missing,
                list_sigma_lasso=list_sigma_lasso, list_sigma_error=list_sigma_error, folds = folds))
  }
  if (p1 == 0) {
    return(list(sigma_global_corrupted_additive = sigma_global_corrupted_additive, 
                sigma_global_corrupted_missing = sigma_global_corrupted_missing,
                list_PSD_lasso_additive = list_PSD_lasso_additive, list_PSD_error_additive = list_PSD_error_additive,
                list_PSD_lasso_missing =list_PSD_lasso_missing, list_PSD_error_missing = list_PSD_error_missing,
                folds = folds))
  }
}

generalcoco <- function(Z, y, n, p, p1=NULL, p2=NULL, p3=NULL, center.Z=TRUE, scale.Z=TRUE, center.y=TRUE, scale.y=TRUE,
                        lambda.factor=ifelse(dim(Z)[1] < dim(Z)[2], 0.01, 0.001), step=100, K=4, mu=10, tau=NULL,
                        etol=1e-4, optTol=1e-5, earlyStopping_max=10, penalty=c("lasso","SCAD"), mode="ADMM") {
  
  this.call <- match.call()
  fit <- blockwise_coordinate_descent_general(Z=Z,
                                              y=y,
                                              n=n,
                                              p=p,
                                              p1=p1,
                                              p2=p2,
                                              p3=p3,
                                              center.Z = center.Z,
                                              scale.Z = scale.Z,
                                              center.y = center.y,
                                              scale.y = scale.y,
                                              lambda.factor = lambda.factor,
                                              step = step,
                                              K = K,
                                              mu = mu,
                                              tau = tau,
                                              etol = etol,
                                              optTol = optTol,
                                              earlyStopping_max = earlyStopping_max,
                                              penalty=penalty,
                                              mode=mode)
  fit$call <- this.call
  class(fit) <- "coco"
  return(fit)
}

lasso_covariance <- function(n, p, lambda, control=list(maxIter=1000, optTol=10 ^ (-5), zeroThreshold=10 ^ (-6)), 
                             XX, Xy, beta.start, penalty=c("lasso", "SCAD")) {
  
  beta <- beta.start
  wp <- beta
  m <- 1
  
  ## compute the product of XX with beta
  s <- XX %*% beta
  lambda0 <- lambda
  while (m < control$maxIter) {
    beta_old <- beta
    for (j in 1:p) {
      # Compute the Shoot and Update the variable
      S0 <- s[j] - XX[j,j]*beta_old[j] - Xy[j]
      if (sum(is.na(S0)) >= 1) {
        beta[j] <- 0
        next
      }
      
      w_j <- 1
      if(penalty == "SCAD"){
        a <- 3.7
        if(abs(beta[j]) <= lambda){
          w_j <- 1
        }else if(abs(beta[j]) <= a*lambda){
          w_j <- (a*lambda - abs(beta[j]))/(lambda*(a-1))
        }else{
          w_j <- 0
        }
      }
      
      lambda <- w_j * lambda0
      if (S0 > lambda){
        beta[j] <- (lambda - S0)/XX[j, j]
        s <- s + XX[,j]*( beta[j] - beta_old[j])
      } 
      if (S0 < -1 * lambda){
        beta[j] <- (-1 * lambda - S0)/XX[j, j]
        s <- s + XX[,j]*( beta[j] - beta_old[j])
      }
      if (abs(S0) <= lambda) 
        beta[j] <- 0
    }
    # Update
    wp <- cbind(wp, beta)
    # Check termination for early stopping
    if (sum(abs(beta - beta_old), na.rm = TRUE) < control$optTol) {
      break
    }
    m <- m + 1
  }
  w <- beta
  #We impose very small coefficients to be equal to zero
  w[abs(w) < control$zeroThreshold] <- 0
  return(list(coefficients = w, coef.list = wp, num.it = m))
}

lasso_covariance_block <- function(n,p1, p2, X1, Z2, y, sigma1, sigma2, lambda, noise=c("additive", "missing", "HM"),
                                   ratio_matrix=NULL, control=list(maxIter=1000, optTol=10 ^ (-5), zeroThreshold=10 ^ (-6)),
                                   beta1.start, beta2.start, penalty=c("lasso","SCAD")) {
  
  
  beta1 <- beta1.start
  beta2 <- beta2.start
  objective_function <- 1000
  m <- 1
  error_beta1 <- matrix(0,control$maxIter,1)
  error_beta2 <- matrix(0,control$maxIter,1)
  if (noise=="additive"){
    Z2_tilde <- Z2
    rho1 <- 1/n * (t(X1) %*% y)
    rho2 <- 1/n * (t(Z2) %*% y) 
  }
  else if (noise=="missing"){
    Z2_tilde <- sapply(1:p2,function(j)Z2[,j] / diag(ratio_matrix)[j])
    rho1 <- 1/n * (t(X1) %*% y)
    rho2 <- 1/n * (t(Z2) %*% y) / diag(ratio_matrix)
  }
  
  
  while (m < control$maxIter){
    beta1.old <- beta1
    beta2.old <- beta2
    objective_function.old <- objective_function
    
    # First step of the block descent : dealing with uncorrupted predictors
    if (noise == "additive"){
      Xy1 <- 1/n * t(X1) %*% (y - Z2 %*% beta2)
    } else if ((noise == "missing") || (noise == "HM")){
      Z2_tilde <- sapply(1:p2,function(j)Z2[,j] / diag(ratio_matrix)[j])
      Xy1 <- 1/n * t(X1) %*% (y - Z2_tilde %*% beta2 )
    }
    beta1 <- lasso_covariance(n=n, p=p1, lambda=lambda, XX=sigma1, Xy = Xy1, beta.start = beta1.old, penalty=penalty)$coefficients
    
    
    # Second step of the block descent : dealing with corrupted predictors
    if (noise == "additive"){
      Xy2 <- 1/n * t(Z2) %*% (y - X1 %*% beta1)
    } else if ((noise == "missing") || (noise == "HM")){
      Xy2 <- 1/n * (t(Z2) %*% (y - X1 %*% beta1)) / diag(ratio_matrix)
    }
    beta2 <- lasso_covariance(n=n, p=p2, lambda=lambda, XX=sigma2, Xy = Xy2, beta.start = beta2.old, penalty=penalty)$coefficients
    
    error_beta1[m,1] = sum(abs(beta1 - beta1.old))
    error_beta2[m,1] = sum(abs(beta2 - beta2.old))
    if ((sum(abs(beta1 - beta1.old), na.rm = TRUE) < control$optTol) && (sum(abs(beta2 - beta2.old), na.rm = TRUE) < control$optTol)) {
      break
    }
    
    
    objective_function <- t(beta1) %*% sigma1 %*% beta1 + t(beta2) %*% sigma2 %*% beta2 - 2 * t(rho1) %*% beta1 - 2 * t(rho2) %*% beta2 + 2 * t(beta2) %*% t(Z2_tilde) %*% X1 %*% beta1
    # if (objective_function > objective_function.old){
    #   print(paste0("Warning at step ",m))
    # }
    m <- m + 1
  }
  #We impose very small coefficients to be equal to zero
  beta1[abs(beta1) < control$zeroThreshold] <- 0
  beta2[abs(beta2) < control$zeroThreshold] <- 0
  return(list(coefficients.beta1 = beta1, coefficients.beta2 = beta2, num.it = m))
}

lasso_covariance_block_general <- function(n, p1, p2, p3, X1, Z2, Z3, y, sigma1, sigma2, sigma3, lambda, ratio_matrix=NULL,
                                           control=list(maxIter=1000, optTol=10 ^ (-5), zeroThreshold=10 ^ (-6)),
                                           beta1.start, beta2.start, beta3.start, penalty=c("lasso","SCAD")) {
  
  if (p1 != 0) {
    beta1 <- beta1.start
    beta2 <- beta2.start
    beta3 <- beta3.start
    objective_function <- 1000
    m <- 1
    error_beta1 <- matrix(0,control$maxIter,1)
    error_beta2 <- matrix(0,control$maxIter,1)
    error_beta3 <- matrix(0,control$maxIter,1)
    
    rho1 <- 1/n * (t(X1) %*% y)
    rho2 <- 1/n * (t(Z2) %*% y)
    Z3_tilde <- sapply(1:p3,function(j)Z3[,j] / diag(ratio_matrix)[j])
    rho3 <- 1/n * (t(Z3) %*% y) / diag(ratio_matrix)
    
    while (m < control$maxIter){
      beta1.old <- beta1
      beta2.old <- beta2
      beta3.old <- beta3
      objective_function.old <- objective_function
      
      # First step of the block descent : dealing with uncorrupted predictors
      Z3_tilde <- sapply(1:p3,function(j)Z3[,j] / diag(ratio_matrix)[j])
      Xy1 <- 1/n * t(X1) %*% (y - Z2 %*% beta2 - Z3_tilde %*% beta3)
      beta1 <- lasso_covariance(n=n, p=p1, lambda=lambda, XX=sigma1, Xy = Xy1, beta.start = beta1.old, penalty=penalty)$coefficients
      
      
      # Second step of the block descent : dealing with corrupted predictors - additive
      Xy2 <- 1/n * t(Z2) %*% (y - X1 %*% beta1 - Z3_tilde %*% beta3)
      beta2 <- lasso_covariance(n=n, p=p2, lambda=lambda, XX=sigma2, Xy = Xy2, beta.start = beta2.old, penalty=penalty)$coefficients
      
      # Third step of the block descent : dealing with corrupted predictors - missing
      Xy3 <- 1/n * (t(Z3) %*% (y - X1 %*% beta1 - Z2 %*% beta2)) / diag(ratio_matrix)
      beta3 <- lasso_covariance(n=n, p=p3, lambda=lambda, XX=sigma3, Xy = Xy3, beta.start = beta3.old, penalty=penalty)$coefficients
      
      error_beta1[m,1] = sum(abs(beta1 - beta1.old))
      error_beta2[m,1] = sum(abs(beta2 - beta2.old))
      error_beta3[m,1] = sum(abs(beta3 - beta3.old))
      if ((sum(abs(beta1 - beta1.old), na.rm = TRUE) < control$optTol) && (sum(abs(beta2 - beta2.old), na.rm = TRUE) < control$optTol) && (sum(abs(beta3 - beta3.old), na.rm = TRUE) < control$optTol)) {
        break
      }
      
      
      objective_function <- t(beta1) %*% sigma1 %*% beta1 + 
        t(beta2) %*% sigma2 %*% beta2 +
        t(beta3) %*% sigma3 %*% beta3 -
        2 * t(rho1) %*% beta1 - 
        2 * t(rho2) %*% beta2 -
        2 * t(rho3) %*% beta3 + 
        2 * t(beta2) %*% t(Z2) %*% X1 %*% beta1 +
        2 * t(beta3) %*% t(Z3_tilde) %*% X1 %*% beta1 +
        2 * t(beta3) %*% t(Z3_tilde) %*% Z2 %*% beta2
      # if (objective_function > objective_function.old){
      #   print(paste0("Warning at step ",m))
      # }
      m <- m + 1
    }
    #We impose very small coefficients to be equal to zero
    beta1[abs(beta1) < control$zeroThreshold] <- 0
    beta2[abs(beta2) < control$zeroThreshold] <- 0
    beta3[abs(beta3) < control$zeroThreshold] <- 0
    return(list(coefficients.beta1 = beta1, coefficients.beta2 = beta2, coefficients.beta3 = beta3, num.it = m))
  }
  if (p1 == 0) {
    beta2 <- beta2.start
    beta3 <- beta3.start
    objective_function <- 1000
    m <- 1
    error_beta2 <- matrix(0,control$maxIter,1)
    error_beta3 <- matrix(0,control$maxIter,1)
    
    rho2 <- 1/n * (t(Z2) %*% y)
    Z3_tilde <- sapply(1:p3,function(j)Z3[,j] / diag(ratio_matrix)[j])
    rho3 <- 1/n * (t(Z3) %*% y) / diag(ratio_matrix)
    
    while (m < control$maxIter){
      beta2.old <- beta2
      beta3.old <- beta3
      objective_function.old <- objective_function
      
      # First step of the block descent : dealing with uncorrupted predictors
      Z3_tilde <- sapply(1:p3,function(j)Z3[,j] / diag(ratio_matrix)[j])
      
      # Second step of the block descent : dealing with corrupted predictors - additive
      Xy2 <- 1/n * t(Z2) %*% (y - Z3_tilde %*% beta3)
      beta2 <- lasso_covariance(n=n, p=p2, lambda=lambda, XX=sigma2, Xy = Xy2, beta.start = beta2.old, penalty=penalty)$coefficients
      
      # Third step of the block descent : dealing with corrupted predictors - missing
      Xy3 <- 1/n * (t(Z3) %*% (y - Z2 %*% beta2)) / diag(ratio_matrix)
      beta3 <- lasso_covariance(n=n, p=p3, lambda=lambda, XX=sigma3, Xy = Xy3, beta.start = beta3.old, penalty=penalty)$coefficients
      
      error_beta2[m,1] = sum(abs(beta2 - beta2.old))
      error_beta3[m,1] = sum(abs(beta3 - beta3.old))
      if ((sum(abs(beta2 - beta2.old), na.rm = TRUE) < control$optTol) && (sum(abs(beta3 - beta3.old), na.rm = TRUE) < control$optTol)) {
        break
      }
      
      
      objective_function <- t(beta2) %*% sigma2 %*% beta2 +
        t(beta3) %*% sigma3 %*% beta3 -
        2 * t(rho2) %*% beta2 -
        2 * t(rho3) %*% beta3 + 
        2 * t(beta3) %*% t(Z3_tilde) %*% Z2 %*% beta2
      # if (objective_function > objective_function.old){
      #   print(paste0("Warning at step ",m))
      # }
      m <- m + 1
    }
    #We impose very small coefficients to be equal to zero
    beta2[abs(beta2) < control$zeroThreshold] <- 0
    beta3[abs(beta3) < control$zeroThreshold] <- 0
    return(list(coefficients.beta2 = beta2, coefficients.beta3 = beta3, num.it = m))
  }
}

print.coco <- function(x, ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(x$data_error[,c("lambda","error")])
}

predict.coco <- function(object, newx, s=NULL, lambda.pred=NULL, type=c("response", "coefficients"), ...) {
  
  if (type == "response"){
    if (missing(newx)){
      stop("newx is missing. Please supply the vector of covariates for which response is required.")
    }
  }
  lambda.seq = object$data_beta[,"lambda"]
  nbeta = object$data_beta
  if(!is.null(object$vnames)){
    
    colnames(nbeta) <- c("lambda", object$vnames)
  }
  if (!is.null(s)){
    index <- match(s,lambda.seq)
    nbeta <- nbeta[index,]
  }
  if(type == "coefficients"){
    coef = t(nbeta[,-1])
    if (!is.null(s)){
      colnames(coef) <- c("Coefficient")
    }
    
    return(coef)
  }else{
    #browser()
    mean.Z = object$mean.Z
    mean.y = object$mean.y
    sd.Z = object$sd.Z
    sd.y = object$sd.y
    matcoef = object$beta.sd
    if(!is.null(lambda.pred)){
      index <- match(lambda.pred,lambda.seq)
      matcoef <- t(object$data_beta[index,2:dim(object$data_beta)[2]])
    }
    newx <- (newx - mean.Z)/sd.Z * sd.y
    y <- newx %*% matcoef + mean.y
    return(y)
  }
}

coef.coco <- function(object, s=NULL, ...) {
  stats::predict(object, s = s, type = "coefficients", ...)
}

lambda_max <- function(Z, y, n, ratio_matrix=NULL, noise=c("additive", "missing")) {
  
  if (noise == "additive"){
    rho_tilde <- 1/n*t(Z)%*%y
  } else if ((noise=="missing") || (noise=="HM")){
    #rho_tilde <- 1/n*t(Z)%*%y/ (1 - probs)
    rho_tilde <- 1/n*t(Z)%*%y/ diag(ratio_matrix)
  }
  max(abs(rho_tilde))
}

rescale_without_NA <- function(j, Z) {
  m <- mean(Z[which(!is.na(Z[,j]), arr.ind = TRUE),j])
  Z[,j] - m
}

mean_without_NA <- function(j, Z) {
  m <- mean(Z[which(!is.na(Z[,j]), arr.ind = TRUE),j])
  m
}

sd_without_NA_block <- function(j, Z) {
  sd <- stats::sd(Z[which(!is.na(Z[,j]), arr.ind = TRUE),j])
  sd
}

change_NA_value <- function(j, Z) {
  Z[which(is.na(Z[,j]), arr.ind = TRUE),j] <- 0
  Z[,j]
}

scale_manual_with_sd <- function(j, Z, v) {
  sd <- v[j]
  if (sd != 0){
    return(Z[,j]/sd)
  }else{
    return (Z[,j])
  }
}

cross_validation_function <- function(k, n, p, lambda_step, list_matrices_lasso, list_rho_lasso,
                                      list_matrices_error, list_rho_error, beta_start, penalty=c("lasso","SCAD")) {
  
  #Solving the lasso problem without the kth fold
  sigma_train <- list_matrices_lasso[[k]]
  rho_train <- list_rho_lasso[[k]] 
  coef_lambda = lasso_covariance(n=n, p=p,lambda=lambda_step, XX=sigma_train, Xy=rho_train, beta.start=beta_start, penalty=penalty)$coefficients
  
  #Calculating the error on the remaining fold
  sigma_test <- list_matrices_error[[k]]
  rho_test <- list_rho_error[[k]]
  error <- t(coef_lambda)%*%sigma_test%*%coef_lambda - 2*t(rho_test)%*%coef_lambda
  
  error
}

pathwise_coordinate_descent <- function(Z, y, n, p, center.Z=TRUE, scale.Z=TRUE, center.y=TRUE,
                                        scale.y=TRUE, lambda.factor=ifelse(dim(Z)[1] < dim(Z)[2], 0.01, 0.001),
                                        step=100, K=4, mu=10, tau=NULL, etol=1e-4, optTol=1e-10, earlyStopping_max=10,
                                        noise=c("additive", "missing"), penalty=c("lasso","SCAD"), mode="ADMM") {
  
  
  nrows = nrow(Z)
  ncols = ncol(Z)
  vnames = colnames(Z)
  
  if(!(is.matrix(Z))){
    stop("Z has to be a matrix")
  }
  if(!(is.matrix(y))){
    stop("y has to be a matrix")
  }
  if(!is.null(tau) & !is.numeric(tau)){
    stop("tau must be numeric")
  }
  if(n != nrows){
    stop(paste("Number of rows in Z (", nrows, ") different from n(", n, ")"),sep="")
  }
  if(p != ncols){
    stop(paste("Number of columns in Z (", ncols, ") different from p (", p, ")"),sep="")
  }
  if (nrows != dim(y)[1]){
    stop(paste("Number of rows in Z (", nrows, ") different from number of rows in y (", dim(y)[1], ") "),sep="")
  }
  if (!is.numeric(y)) {
    stop("The response y must be numeric. Factors must be converted to numeric")
  }
  if (any(is.na(y))){
    stop("The response contains NA values. Remove NA values before calling the function.")
  }
  if(lambda.factor >= 1){
    stop("lambda factor should be smaller than 1")
  }
  if (n %% K != 0){
    stop("K should be a divider of n")
  }
  if(mu>500 || mu<1){
    warning(paste("Mu value (", mu, ") is not in the usual range (10-500)"))
  }
  if(noise=="missing" && center.Z == FALSE){
    stop("When noise is equal to missing, it is required to center matrix Z. Use center.Z=TRUE.")
  }
  if(scale.Z == FALSE && noise=='missing'){
    warning("When noise is equal to missing, it is recommended to use scale.Z equal to TRUE in order 
            to obtain trustworthy results.")
  }
  if(scale.Z == TRUE && noise=='additive'){
    warning("When noise is equal to additive, it is recommended to use scale.Z equal to FALSE in order 
            to obtain trustworthy results. Otherwise, the scaling should be taken into account
            when introducing the error parameter as a function parameter.")
  }
  
  ratio_matrix = NULL
  if (noise=="missing"){
    ratio_matrix = matrix(0,p,p)
    
    for (i in 1:p){
      for (j in i:p){
        n_ij = length(intersect(which(!is.na(Z[,i])),which(!is.na(Z[,j]))))
        ratio_matrix[i,j] = n_ij
        ratio_matrix[j,i] = n_ij
      }
    }
    ratio_matrix = ratio_matrix/n
  }
  
  mean.Z = sapply(1:p, function(j)mean_without_NA(j,Z))
  sd.Z <- sapply(1:p, function(j)sd_without_NA_block(j,Z))
  
  if (center.Z == TRUE){
    if (scale.Z == TRUE){
      Z = sapply(1:p, function(j)rescale_without_NA(j,Z))
      Z = sapply(1:p, function(j)change_NA_value(j,Z))
      Z = sapply(1:p, function(j)scale_manual_with_sd(j,Z,sd.Z))
    }else{
      Z = sapply(1:p, function(j)rescale_without_NA(j,Z))
      Z = sapply(1:p, function(j)change_NA_value(j,Z))
    }
  }else{
    if(scale.Z == TRUE){
      Z = sapply(1:p, function(j)scale_manual_with_sd(j,Z,sd.Z))
    }
  }
  
  
  
  mean.y = mean(y)
  sd.y = stats::sd(y)
  
  if (center.y == TRUE){
    if (scale.y == TRUE){
      y = scale(y, center = TRUE, scale = TRUE)
    }else{
      y = scale(y, center = TRUE, scale = FALSE)
    }
  }else{
    if (scale.y == TRUE){
      y = scale(y, center = FALSE, scale = TRUE)
    }
  }
  
  
  
  n_without_fold = n - floor(n/K)
  n_one_fold = floor(n/K)
  earlyStopping = step
  
  lambda_max <- lambda_max(Z=Z,y=y,n=n,ratio_matrix=ratio_matrix,noise=noise)
  lambda_min <- lambda.factor*lambda_max
  lambda_list <- emdbook::lseq(lambda_max,lambda_min,step)
  beta_start <- rep(0,p)
  best.lambda <- lambda_max
  beta.opt <- beta_start
  best.error <- 1000
  error_list <- matrix(0,step,4)
  error <- 1000
  earlyStopping_high = 0
  
  matrix_beta <- matrix(0,step,p)
  
  ### Creating the K matrices we are going to use for cross validation
  output = cv_covariance_matrices(K=K, mat=Z, y=y, p=p, mu=mu, tau=tau, ratio_matrix = ratio_matrix, etol=etol, noise = noise, mode=mode)
  list_matrices_lasso = output$list_matrices_lasso
  list_matrices_error = output$list_matrices_error
  list_rho_lasso = output$list_rho_lasso
  list_rho_error = output$list_rho_error
  ZZ = output$sigma_global
  Zy = output$rho_global
  for (i in 1:step){
    lambda_step <- lambda_list[i]
    error_old <- error
    error <- 0
    
    out = sapply(1:K, function(k)cross_validation_function(k,
                                                           n,
                                                           p,
                                                           lambda_step,
                                                           list_matrices_lasso,
                                                           list_rho_lasso,
                                                           list_matrices_error,
                                                           list_rho_error,
                                                           beta_start,
                                                           penalty=penalty))
    error = mean(out)
    sd_low = stats::quantile(out, probs = c(0.1))
    sd_high = stats::quantile(out, probs = c(0.9))
    error_list[i,1] <- error
    error_list[i,2] <- sd_low
    error_list[i,3] <- sd_high
    error_list[i,4] <- stats::sd(out)
    coef_tot = lasso_covariance(n=n, p=p, lambda=lambda_step, XX=ZZ, Xy=Zy, beta.start = beta_start, penalty=penalty)$coefficients
    beta_start <- coef_tot
    matrix_beta[i,] <- beta_start
    
    ### Checking for optimal parameters
    if (error <= best.error){
      best.error <- error
      best.lambda <- lambda_step
      beta.opt <- coef_tot
    }
    
    ## Early stopping
    if (abs(error - error_old) < optTol){
      print("Early Stopping because of convergence of the error")
      earlyStopping = i
      break
    }
    
    
    if (error > best.error){
      earlyStopping_high = earlyStopping_high +1
      if (earlyStopping_high >= earlyStopping_max){
        print("Early stopping because of error getting too high")
        earlyStopping = i
        break
      }
    }
  }
  df <- data.frame(lambda=lambda_list[1:earlyStopping], error=error_list[1:earlyStopping,1],error.inf=error_list[1:earlyStopping,2],error.sup=error_list[1:earlyStopping,3],error.sd=error_list[1:earlyStopping,4])
  step.min <- which(df[,"error"] == best.error)
  sd.best <- df[step.min,"error.sd"]
  step.sd <- max(which(df[,"error"] > best.error + sd.best & df[,"lambda"] > df[step.min,"lambda"]))
  lambda.sd <- df[step.sd,"lambda"]
  beta.sd <- matrix_beta[step.sd,]
  
  data_intermediate <- data.frame(matrix_beta[1:earlyStopping,])
  names(data_intermediate) <- sapply(1:p, function(i)paste0("beta",i))
  
  data_beta <- data.frame(lambda = lambda_list[1:earlyStopping])
  
  data_beta <- cbind(data_beta,data_intermediate)
  fit <- list(
    lambda.opt = best.lambda,
    lambda.sd = lambda.sd,
    beta.opt = beta.opt,
    beta.sd = beta.sd,
    data_error = df,
    data_beta = data_beta,
    earlyStopping = earlyStopping,
    vnames = vnames,
    mean.Z = mean.Z,
    sd.Z = sd.Z,
    mean.y = mean.y,
    sd.y = sd.y
  )
  return(fit)
}

plotCoef <- function(object, linetype="dashed", col="black") {
  data_error = object$data_error
  data_beta = object$data_beta
  beta <- object$beta.opt
  beta <- as.matrix(beta)
  best.lambda <- object$lambda.opt
  lambda.sd <- object$lambda.sd
  beta.sd <- object$beta.sd
  beta.sd <- as.matrix(beta.sd)
  
  data_beta.bis <- reshape2::melt(data_beta, id="lambda", value.name="value" ,variable.name="beta")
  ggplot2::ggplot(data=data_beta.bis) + ggplot2::geom_line(data=data_beta.bis,ggplot2::aes(x=log(.data$lambda),y=.data$value,colour=beta)) + ggplot2::geom_vline(xintercept = log(best.lambda), linetype=linetype, colour=col) + ggplot2::geom_vline(xintercept = log(lambda.sd), linetype=linetype, colour=col) + ggplot2::theme(legend.position = "none")+ ggplot2::xlab("Log Lambda") + ggplot2::ylab("Coefficients")
  
}

plotError <- function(object, linetype="dashed", col="black", colLine="red") {
  data_error = object$data_error
  data_beta = object$data_beta
  beta <- object$beta.opt
  beta <- as.matrix(beta)
  best.lambda <- object$lambda.opt
  lambda.sd <- object$lambda.sd
  beta.sd <- object$beta.sd
  beta.sd <- as.matrix(beta.sd)
  
  data_beta.bis <- melt(data_beta, id="lambda", value.name="value" ,variable.name="beta")
  ggplot2::ggplot(data=data_error) + ggplot2::geom_line(data=data_error,ggplot2::aes(log(.data$lambda),.data$error),colour=colLine,size=1) + ggplot2::xlab("Log Lambda") + ggplot2::ylab("Error") + ggplot2::geom_errorbar(data=data_error,ggplot2::aes(x=log(.data$lambda), ymin = .data$error.inf, ymax=.data$error.sup), colour=col) + ggplot2::geom_vline(xintercept = log(best.lambda), linetype=linetype) + ggplot2::geom_vline(xintercept = log(lambda.sd), linetype=linetype)
  
}

cov_autoregressive <- function(p) {
  cov <- matrix(0,p,p)
  for (i in 1:p){
    for (j in 1:p){
      cov[i,j] = 0.5**(abs(i-j))
    }
  }
  cov
}

########################################################

res <- function(itr, WK, XK, ZK, yK, K, n, nK, p, d, s, r1, r2, sigmav) {
  
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
  bhat.dunif <- c(solve(tmp1.dunif) %*% tmp2.dunif)
  
  Jhat.dunif <- Vhat.dunif <- 0
  for (i in 1:length(yK.dunif)) {
    Jhat.dunif <- Jhat.dunif + 
      (WK.dunif[i, ] %*% t(WK.dunif[i, ] - Htilde %*% ZK.dunif[i, ]) - sigmav) / (n * rK.dunif[i] * PI.dunif[i])
    Vhat.dunif <- Vhat.dunif +
      ((WK.dunif[i, ] - Htilde %*% ZK.dunif[i, ]) * c(gmtilde %*% ZK.dunif[i, ] + bhat.dunif %*% WK.dunif[i, ] - yK.dunif[i]) -
         sigmav %*% bhat.dunif) %*% 
      t(((WK.dunif[i, ] - Htilde %*% ZK.dunif[i, ]) * c(gmtilde %*% ZK.dunif[i, ] + bhat.dunif %*% WK.dunif[i, ] - yK.dunif[i]) -
           sigmav %*% bhat.dunif)) / (n * rK.dunif[i] * PI.dunif[i]) ^ 2
  }
  se.dunif <- sqrt(diag(solve(Jhat.dunif) %*% Vhat.dunif %*% solve(Jhat.dunif)))
  
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
  bhat.dmV <- c(solve(tmp1.dmV) %*% tmp2.dmV)
  
  Jhat.dmV <- Vhat.dmV <- 0
  for (i in 1:length(yK.dmV)) {
    Jhat.dmV <- Jhat.dmV + 
      (WK.dmV[i, ] %*% t(WK.dmV[i, ] - Htilde %*% ZK.dmV[i, ]) - sigmav) / (n * rK.dmV[i] * PI.dmV[i])
    Vhat.dmV <- Vhat.dmV +
      ((WK.dmV[i, ] - Htilde %*% ZK.dmV[i, ]) * c(gmtilde %*% ZK.dmV[i, ] + bhat.dmV %*% WK.dmV[i, ] - yK.dmV[i]) -
         sigmav %*% bhat.dmV) %*% 
      t(((WK.dmV[i, ] - Htilde %*% ZK.dmV[i, ]) * c(gmtilde %*% ZK.dmV[i, ] + bhat.dmV %*% WK.dmV[i, ] - yK.dmV[i]) -
           sigmav %*% bhat.dmV)) / (n * rK.dmV[i] * PI.dmV[i]) ^ 2
  }
  se.dmV <- sqrt(diag(solve(Jhat.dmV) %*% Vhat.dmV %*% solve(Jhat.dmV)))
  
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
  bhat.dmVc <- c(solve(tmp1.dmVc) %*% tmp2.dmVc)
  
  Jhat.dmVc <- Vhat.dmVc <- 0
  for (i in 1:length(yK.dmVc)) {
    Jhat.dmVc <- Jhat.dmVc + 
      (WK.dmVc[i, ] %*% t(WK.dmVc[i, ] - Htilde %*% ZK.dmVc[i, ]) - sigmav) / (n * rK.dmVc[i] * PI.dmVc[i])
    Vhat.dmVc <- Vhat.dmVc +
      ((WK.dmVc[i, ] - Htilde %*% ZK.dmVc[i, ]) * c(gmtilde %*% ZK.dmVc[i, ] + bhat.dmVc %*% WK.dmVc[i, ] - yK.dmVc[i]) -
         sigmav %*% bhat.dmVc) %*% 
      t(((WK.dmVc[i, ] - Htilde %*% ZK.dmVc[i, ]) * c(gmtilde %*% ZK.dmVc[i, ] + bhat.dmVc %*% WK.dmVc[i, ] - yK.dmVc[i]) -
           sigmav %*% bhat.dmVc)) / (n * rK.dmVc[i] * PI.dmVc[i]) ^ 2
  }
  se.dmVc <- sqrt(diag(solve(Jhat.dmVc) %*% Vhat.dmVc %*% solve(Jhat.dmVc)))
  
  # method 4: oradmV
  XK.ora <- XK
  ZK.ora <- ZK[, , 1:(s - d + 1)]
  yK.ora <- yK
  
  XK.ora.simp <- c()
  ZK.ora.simp <- c()
  yK.ora.simp <- c()
  r1.ora.simp <- round(r1 * nK / n)
  idx.ora.simp <- matrix(NA, max(r1.ora.simp), K)
  for (k in 1:K) {
    idx.ora.simp[1:r1.ora.simp[k], k] <- sample(1:nK[k], r1.ora.simp[k], T)
    XK.ora.simp <- rbind(XK.ora.simp, XK.ora[idx.ora.simp[1:r1.ora.simp[k], k], k, ])
    ZK.ora.simp <- rbind(ZK.ora.simp, ZK.ora[idx.ora.simp[1:r1.ora.simp[k], k], k, ])
    yK.ora.simp <- c(yK.ora.simp, yK.ora[idx.ora.simp[1:r1.ora.simp[k], k], k])
  }
  tmp1.ora.simp <- t(cbind(XK.ora.simp, ZK.ora.simp)) %*% cbind(XK.ora.simp, ZK.ora.simp) / r1
  tmp2.ora.simp <- t(cbind(XK.ora.simp, ZK.ora.simp)) %*% yK.ora.simp / r1
  theta.ora.simp <- c(solve(tmp1.ora.simp) %*% tmp2.ora.simp)
  
  Jtilde <- t(cbind(XK.ora.simp, ZK.ora.simp)) %*% cbind(XK.ora.simp, ZK.ora.simp) / r1
  PI.opt.ora.dmV <- matrix(NA, max(nK), K)
  r.opt.ora.dmV <- rep(NA, K) 
  
  tmp.ora.dmV <- matrix(NA, max(nK), K)
  for (k in 1:K) {
    tmp.ora.dmV[1:nK[k], k] <- sqrt(rowSums((cbind(XK.ora[1:nK[k], k, ], ZK.ora[1:nK[k], k, ]) %*% solve(Jtilde) * 
                                               c(cbind(XK.ora[1:nK[k], k, ], ZK.ora[1:nK[k], k, ]) %*% theta.ora.simp - yK.ora[1:nK[k], k]))) ^ 2)
    PI.opt.ora.dmV[1:nK[k], k] <- tmp.ora.dmV[1:nK[k], k] / sum(tmp.ora.dmV[1:nK[k], k])
    r.opt.ora.dmV[k] <- sum(tmp.ora.dmV[1:nK[k], k])
  }
  r.opt.ora.dmV <- round(r2 * r.opt.ora.dmV / sum(r.opt.ora.dmV))
  
  XK.ora.dmV <- c()
  ZK.ora.dmV <- c()
  yK.ora.dmV <- c()
  PI.ora.dmV <- c()
  rK.ora.dmV <- c()
  idx.ora.dmV <- matrix(NA, max(r.opt.ora.dmV), K)
  for (k in 1:K) {
    idx.ora.dmV[1:r.opt.ora.dmV[k], k] <- sample(1:nK[k], r.opt.ora.dmV[k], T, PI.opt.ora.dmV[1:nK[k], k])
    XK.ora.dmV <- rbind(XK.ora.dmV, XK.ora[idx.ora.dmV[1:r.opt.ora.dmV[k], k], k, ])
    ZK.ora.dmV <- rbind(ZK.ora.dmV, ZK.ora[idx.ora.dmV[1:r.opt.ora.dmV[k], k], k, ])
    yK.ora.dmV <- c(yK.ora.dmV, yK.ora[idx.ora.dmV[1:r.opt.ora.dmV[k], k], k])
    PI.ora.dmV <- c(PI.ora.dmV, PI.opt.ora.dmV[idx.ora.dmV[1:r.opt.ora.dmV[k], k], k])
    rK.ora.dmV <- c(rK.ora.dmV, rep(r.opt.ora.dmV[k], r.opt.ora.dmV[k]))
  }
  
  tmp1.ora.dmV <- tmp2.ora.dmV <- 0
  for (i in 1:length(yK.ora.dmV)) {
    tmp1.ora.dmV <- tmp1.ora.dmV + 
      c(XK.ora.dmV[i, ], ZK.ora.dmV[i, ]) %*% t(c(XK.ora.dmV[i, ], ZK.ora.dmV[i, ])) / (rK.ora.dmV[i] * PI.ora.dmV[i])
    tmp2.ora.dmV <- tmp2.ora.dmV +
      c(XK.ora.dmV[i, ], ZK.ora.dmV[i, ]) * yK.ora.dmV[i] / (rK.ora.dmV[i] * PI.ora.dmV[i])
  }
  that.oradmV <- c(solve(tmp1.ora.dmV) %*% tmp2.ora.dmV)
  bhat.oradmV <- c(that.oradmV[1:d])
  
  # estimated standard error
  Jhat.ora.dmV <- Vhat.ora.dmV <- 0
  for (i in 1:length(yK.ora.dmV)) {
    Jhat.ora.dmV <- Jhat.ora.dmV + 
      (c(XK.ora.dmV[i, ], ZK.ora.dmV[i, ]) %*% t(c(XK.ora.dmV[i, ], ZK.ora.dmV[i, ]))) / (n * rK.ora.dmV[i] * PI.ora.dmV[i])
    Vhat.ora.dmV <- Vhat.ora.dmV +
      (c(XK.ora.dmV[i, ], ZK.ora.dmV[i, ]) * c(c(XK.ora.dmV[i, ], ZK.ora.dmV[i, ]) %*% that.oradmV - yK.ora.dmV[i])) %*% 
      t(c(XK.ora.dmV[i, ], ZK.ora.dmV[i, ]) * c(c(XK.ora.dmV[i, ], ZK.ora.dmV[i, ]) %*% that.oradmV - yK.ora.dmV[i])) / (n * rK.ora.dmV[i] * PI.ora.dmV[i]) ^ 2
  }
  se.oradmV <- sqrt(diag(solve(Jhat.ora.dmV) %*% Vhat.ora.dmV %*% solve(Jhat.ora.dmV)))[1:d]
  
  # method 5: undmV 
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
  cv.fit.ini <- cv.glmnet(cbind(WK.simp, ZK.simp[, -1]), yK.simp, alpha=1, family="gaussian")
  fit.ini <- glmnet(cbind(WK.simp, ZK.simp[, -1]), yK.simp, alpha=1, lambda=cv.fit.ini$lambda.min)
  theta.ini <- as.vector(c(fit.ini$beta[1:d], fit.ini$a0, fit.ini$beta[-(1:d)]))
  theta.ini[abs(theta.ini) < 0.5] <- 0
  if (sum(theta.ini[1:d]!=0)==d) {
    index <- theta.ini!=0
    mat <- cbind(WK.simp, ZK.simp)[, index]
    sigmahat <- t(mat) %*% mat / r1
    rhohat <- t(mat) %*% yK.simp / r1
    if (min(eigen(sigmahat)$values) > 0) {
      theta.refit <- solve(sigmahat) %*% rhohat
    } 
    theta.refit[abs(theta.refit) < 0.5] <- 0
    theta.ini[index] <- theta.refit
  }
  
  # estimate of projection matrix H
  Hfit <- mvr(WK.simp, ZK.simp[, -1])
  Hhat <- rbind(as.vector(Hfit$muhat), Hfit$Bhat)
  
  # optimal probability
  Jtilde <- t(WK.simp) %*% (WK.simp - ZK.simp %*% t(Htilde)) / r1
  PI.opt.undmV <- matrix(NA, max(nK), K)
  r.opt.undmV <- rep(NA, K) 
  
  tmp.undmV <- matrix(NA, max(nK), K)
  for (k in 1:K) {
    tmp.undmV[1:nK[k], k] <- sqrt(rowSums(((WK[1:nK[k], k, ] - ZK[1:nK[k], k, ] %*% t(Htilde)) %*% solve(Jtilde) * 
                                             c(WK[1:nK[k], k, ] %*% bttilde + ZK[1:nK[k], k, ] %*% gmtilde - yK[1:nK[k], k])) ^ 2))
    PI.opt.undmV[1:nK[k], k] <- tmp.undmV[1:nK[k], k] / sum(tmp.undmV[1:nK[k], k])
    r.opt.undmV[k] <- sum(tmp.undmV[1:nK[k], k])
  }
  r.opt.undmV <- round(r2 * r.opt.undmV / sum(r.opt.undmV))
  
  WK.undmV <- c()
  ZK.undmV <- c()
  yK.undmV <- c()
  PI.undmV <- c()
  rK.undmV <- c()
  idx.undmV <- matrix(NA, max(r.opt.undmV), K)
  for (k in 1:K) {
    idx.undmV[1:r.opt.undmV[k], k] <- sample(1:nK[k], r.opt.undmV[k], T, PI.opt.undmV[1:nK[k], k])
    WK.undmV <- rbind(WK.undmV, WK[idx.undmV[1:r.opt.undmV[k], k], k, ])
    ZK.undmV <- rbind(ZK.undmV, ZK[idx.undmV[1:r.opt.undmV[k], k], k, ])
    yK.undmV <- c(yK.undmV, yK[idx.undmV[1:r.opt.undmV[k], k], k])
    PI.undmV <- c(PI.undmV, PI.opt.undmV[idx.undmV[1:r.opt.undmV[k], k], k])
    rK.undmV <- c(rK.undmV, rep(r.opt.undmV[k], r.opt.undmV[k]))
  }
  
  tmp1.undmV <- tmp2.undmV <- 0
  for (i in 1:length(yK.undmV)) {
    tmp1.undmV <- tmp1.undmV + 
      (WK.undmV[i, ] - Htilde %*% ZK.undmV[i, ]) %*% t(WK.undmV[i, ]) / (rK.undmV[i] * PI.undmV[i])
    tmp2.undmV <- tmp2.undmV +
      ((WK.undmV[i, ] - Htilde %*% ZK.undmV[i, ]) * c(yK.undmV[i] - gmtilde %*% ZK.undmV[i, ])) / (rK.undmV[i] * PI.undmV[i])
  }
  bhat.undmV <- c(solve(tmp1.undmV) %*% tmp2.undmV)
  
  # estimated standard error
  Jhat.undmV <- Vhat.undmV <- 0
  for (i in 1:length(yK.undmV)) {
    Jhat.undmV <- Jhat.undmV + 
      WK.undmV[i, ] %*% t(WK.undmV[i, ] - Htilde %*% ZK.undmV[i, ]) / (n * rK.undmV[i] * PI.undmV[i])
    Vhat.undmV <- Vhat.undmV +
      ((WK.undmV[i, ] - Htilde %*% ZK.undmV[i, ]) * c(yK.undmV[i] - gmtilde %*% ZK.undmV[i, ] - bhat.undmV %*% WK.undmV[i, ])) %*% 
      t((WK.undmV[i, ] - Htilde %*% ZK.undmV[i, ]) * c(yK.undmV[i] - gmtilde %*% ZK.undmV[i, ] - bhat.undmV %*% WK.undmV[i, ])) / 
      (n * rK.undmV[i] * PI.undmV[i]) ^ 2
  }
  se.undmV <- sqrt(diag(solve(Jhat.undmV) %*% Vhat.undmV %*% solve(Jhat.undmV)))
  return(list(bhat.dmV=bhat.dmV, bhat.dmVc=bhat.dmVc, bhat.dunif=bhat.dunif, 
              bhat.oradmV=bhat.oradmV, bhat.undmV=bhat.undmV,
              se.dmV=se.dmV, se.dmVc=se.dmVc, se.dunif=se.dunif, 
              se.oradmV=se.oradmV, se.undmV=se.undmV))
}

case <- c(1, 1)
# case <- c(2, 1)
# case <- c(3, 1)
# case <- c(4, 1)
# case <- c(1, 2)
# case <- c(2, 2)
# case <- c(3, 2)
# case <- c(4, 2)
if (case[2]==1) {
  K <- 5
  nK <- c(2, 3, 4, 5, 6) * 1e4
} else {
  K <- 25
  nK <- rep(8000, K)
}
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
r21 <- 400
r22 <- 600
r23 <- 800
r24 <- 1000
itr.max <- 500
bhat.dmV <- bhat.dmVc <- bhat.dunif <- bhat.oradmV <- bhat.undmV <- array(NA, c(itr.max, d, 4))
se.dmV <- se.dmVc <- se.dunif <- se.oradmV <- se.undmV <- array(NA, c(itr.max, d, 4))

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
  if (case[1]==1) {
    eK[1:nK[k], k] <- rnorm(nK[k], mean=0, sd=3)
  } else if (case[1]==2) {
    eK[1:nK[k], k] <- rt(nK[k], df=3)
  } else if (case[1]==3) {
    eK[1:nK[k], k] <- abs(XK[1:nK[k], k, 1]) * rnorm(nK[k], 0, sd=3)
  } else {
    eK[1:nK[k], k] <- rexp(nK[k], rate=1)
  }
  yK[1:nK[k], k] <- XK[1:nK[k], k, ] %*% bt + ZK[1:nK[k], k, ] %*% gm + eK[1:nK[k], k]
  WK[1:nK[k], k, ] <- XK[1:nK[k], k, ] + VK[1:nK[k], k, ]
}

time1 <- Sys.time()
for (itr in 1:itr.max) {
  
  cat(itr)
  
  res1 <- res(itr, WK, XK, ZK, yK, K, n, nK, p, d, s, r1, r2=r21, sigmav)
  bhat.dmV[itr, , 1] <- res1$bhat.dmV
  bhat.dmVc[itr, , 1] <- res1$bhat.dmVc
  bhat.dunif[itr, , 1] <- res1$bhat.dunif
  bhat.oradmV[itr, , 1] <- res1$bhat.oradmV
  bhat.undmV[itr, , 1] <- res1$bhat.undmV
  se.dmV[itr, , 1] <- res1$se.dmV
  se.dmVc[itr, , 1] <- res1$se.dmVc
  se.dunif[itr, , 1] <- res1$se.dunif
  se.oradmV[itr, , 1] <- res1$se.oradmV
  se.undmV[itr, , 1] <- res1$se.undmV
  
  res2 <- res(itr, WK, XK, ZK, yK, K, n, nK, p, d, s, r1, r2=r22, sigmav)
  bhat.dmV[itr, , 2] <- res2$bhat.dmV
  bhat.dmVc[itr, , 2] <- res2$bhat.dmVc
  bhat.dunif[itr, , 2] <- res2$bhat.dunif
  bhat.oradmV[itr, , 2] <- res2$bhat.oradmV
  bhat.undmV[itr, , 2] <- res2$bhat.undmV
  se.dmV[itr, , 2] <- res2$se.dmV
  se.dmVc[itr, , 2] <- res2$se.dmVc
  se.dunif[itr, , 2] <- res2$se.dunif
  se.oradmV[itr, , 2] <- res2$se.oradmV
  se.undmV[itr, , 2] <- res2$se.undmV
  
  res3 <- res(itr, WK, XK, ZK, yK, K, n, nK, p, d, s, r1, r2=r23, sigmav)
  bhat.dmV[itr, , 3] <- res3$bhat.dmV
  bhat.dmVc[itr, , 3] <- res3$bhat.dmVc
  bhat.dunif[itr, , 3] <- res3$bhat.dunif
  bhat.oradmV[itr, , 3] <- res3$bhat.oradmV
  bhat.undmV[itr, , 3] <- res3$bhat.undmV
  se.dmV[itr, , 3] <- res3$se.dmV
  se.dmVc[itr, , 3] <- res3$se.dmVc
  se.dunif[itr, , 3] <- res3$se.dunif
  se.oradmV[itr, , 3] <- res3$se.oradmV
  se.undmV[itr, , 3] <- res3$se.undmV
  
  res4 <- res(itr, WK, XK, ZK, yK, K, n, nK, p, d, s, r1, r2=r24, sigmav)
  bhat.dmV[itr, , 4] <- res4$bhat.dmV
  bhat.dmVc[itr, , 4] <- res4$bhat.dmVc
  bhat.dunif[itr, , 4] <- res4$bhat.dunif
  bhat.oradmV[itr, , 4] <- res4$bhat.oradmV
  bhat.undmV[itr, , 4] <- res4$bhat.undmV
  se.dmV[itr, , 4] <- res4$se.dmV
  se.dmVc[itr, , 4] <- res4$se.dmVc
  se.dunif[itr, , 4] <- res4$se.dunif
  se.oradmV[itr, , 4] <- res4$se.oradmV
  se.undmV[itr, , 4] <- res4$se.undmV
}
time2 <- Sys.time()

eva <- function(bhat.dmV, bhat.dmVc, bhat.dunif, bhat.oradmV, bhat.undmV,
                se.dmV, se.dmVc, se.dunif, se.oradmV, se.undmV) {
  alpha <- 0.05
  cp.dmV <- cp.dmVc <- cp.dunif <- cp.oradmV <- cp.undmV <- rep(NA, d)
  for (i in 1:d) {
    cp.dmV[i] <- mean(bhat.dmV[, i] - se.dmV[, i] * qnorm(1 - alpha / 2) < bt[i] & bt[i] < bhat.dmV[, i] + se.dmV[, i] * qnorm(1 - alpha / 2))
    cp.dmVc[i] <- mean(bhat.dmVc[, i] - se.dmVc[, i] * qnorm(1 - alpha / 2) < bt[i] & bt[i] < bhat.dmVc[, i] + se.dmVc[, i] * qnorm(1 - alpha / 2))
    cp.dunif[i] <- mean(bhat.dunif[, i] - se.dunif[, i] * qnorm(1 - alpha / 2) < bt[i] & bt[i] < bhat.dunif[, i] + se.dunif[, i] * qnorm(1 - alpha / 2))
    cp.oradmV[i] <- mean(bhat.oradmV[, i] - se.oradmV[, i] * qnorm(1 - alpha / 2) < bt[i] & bt[i] < bhat.oradmV[, i] + se.oradmV[, i] * qnorm(1 - alpha / 2))
    cp.undmV[i] <- mean(bhat.undmV[, i] - se.undmV[, i] * qnorm(1 - alpha / 2) < bt[i] & bt[i] < bhat.undmV[, i] + se.undmV[, i] * qnorm(1 - alpha / 2))
  }
  cp.mean.dmV <- mean(cp.dmV)
  cp.mean.dmVc <- mean(cp.dmVc)
  cp.mean.dunif <- mean(cp.dunif)
  cp.mean.oradmV <- mean(cp.oradmV)
  cp.mean.undmV <- mean(cp.undmV)
  
  length.mean.dmV <- 1.96 * mean(apply(se.dmV, 2, mean))
  length.mean.dmVc <- 1.96 * mean(apply(se.dmVc, 2, mean))
  length.mean.dunif <- 1.96 * mean(apply(se.dunif, 2, mean))
  length.mean.oradmV <- 1.96 * mean(apply(se.oradmV, 2, mean))
  length.mean.undmV <- 1.96 * mean(apply(se.undmV, 2, mean))
  
  mse.dmV <- mean(diag((bhat.dmV - rep(1, itr.max) %*% t(bt)) %*% t(bhat.dmV - rep(1, itr.max) %*% t(bt))))
  mse.dmVc <- mean(diag((bhat.dmVc - rep(1, itr.max) %*% t(bt)) %*% t(bhat.dmVc - rep(1, itr.max) %*% t(bt))))
  mse.dunif <- mean(diag((bhat.dunif - rep(1, itr.max) %*% t(bt)) %*% t(bhat.dunif - rep(1, itr.max) %*% t(bt))))
  mse.oradmV <- mean(diag((bhat.oradmV - rep(1, itr.max) %*% t(bt)) %*% t(bhat.oradmV - rep(1, itr.max) %*% t(bt))))
  mse.undmV <- mean(diag((bhat.undmV - rep(1, itr.max) %*% t(bt)) %*% t(bhat.undmV - rep(1, itr.max) %*% t(bt))))
  
  absbias.dmV <- mean(abs(apply(bhat.dmV, 2, mean) - bt))
  absbias.dmVc <- mean(abs(apply(bhat.dmVc, 2, mean) - bt))
  absbias.dunif <- mean(abs(apply(bhat.dunif, 2, mean) - bt))
  absbias.oradmV <- mean(abs(apply(bhat.oradmV, 2, mean) - bt))
  absbias.undmV <- mean(abs(apply(bhat.undmV, 2, mean) - bt))
  
  return(rbind(round(c(absbias.dmV, absbias.dmVc, absbias.dunif, absbias.oradmV, absbias.undmV), 3),
               round(c(mse.dmV, mse.dmVc, mse.dunif, mse.oradmV, mse.undmV), 3),
               round(c(cp.mean.dmV, cp.mean.dmVc, cp.mean.dunif, cp.mean.oradmV, cp.mean.undmV), 3),
               round(c(length.mean.dmV, length.mean.dmVc, length.mean.dunif, length.mean.oradmV, length.mean.undmV), 3)))
}

eva(bhat.dmV[, , 1], bhat.dmVc[, , 1], bhat.dunif[, , 1], bhat.oradmV[, , 1], bhat.undmV[, , 1],
     se.dmV[, , 1], se.dmVc[, , 1], se.dunif[, , 1], se.oradmV[, , 1], se.undmV[, , 1])
eva(bhat.dmV[, , 2], bhat.dmVc[, , 2], bhat.dunif[, , 2], bhat.oradmV[, , 2], bhat.undmV[, , 2],
    se.dmV[, , 2], se.dmVc[, , 2], se.dunif[, , 2], se.oradmV[, , 2], se.undmV[, , 2])
eva(bhat.dmV[, , 3], bhat.dmVc[, , 3], bhat.dunif[, , 3], bhat.oradmV[, , 3], bhat.undmV[, , 3],
    se.dmV[, , 3], se.dmVc[, , 3], se.dunif[, , 3], se.oradmV[, , 3], se.undmV[, , 3])
eva(bhat.dmV[, , 4], bhat.dmVc[, , 4], bhat.dunif[, , 4], bhat.oradmV[, , 4], bhat.undmV[, , 4],
    se.dmV[, , 4], se.dmVc[, , 4], se.dunif[, , 4], se.oradmV[, , 4], se.undmV[, , 4])

time2 - time1

