#####################Functions##################################################
csdar <- function(x, y, alpha=0, size=NULL, beta.int=NULL, tau=0.5, intercept=TRUE){
  x <- as.matrix(x)
  y <- as.matrix(y)
  np <- dim(x)
  n <- np[1]
  p <- np[2]
  
  if(intercept){
    meanx <- colMeans(x)
    meany <- mean(y)
  }else{
    meanx <- rep(0, p)
    meany <- 0
  }
  x <- scale(x, meanx, FALSE)
  y <- y - meany
  
  
  if(is.null(size)) {
    size <- seq(1, ceiling(n/log(log(n))/log(p)), 1)
  } 
  tx <- t(x)
  gk <- colSums(x^2) + alpha
  
  sdar <- function(ms){
    if(is.null(beta.int)){
      bk <- rep(0, p)
    }else{
      bk <- beta.int
    }
    dk <- (tx %*% (y - x %*% bk) - alpha*bk)/gk
    Hk <- sqrt(gk)*abs(tau*bk+(1-tau)*dk)
    Ak <- which(Hk >= sort(Hk, decreasing =TRUE)[ms])
    maxiter <- min(2*ms, 20)
    for (k in 1:maxiter) {
      Ik <- setdiff(1:p, Ak)
      bk[Ik] <- 0
      bk[Ak] <- solve(tx[Ak, ] %*% x[, Ak] + diag(alpha, ms)) %*% tx[Ak, ] %*% y
      ek <- y - as.matrix(x[, Ak]) %*% bk[Ak]
      dk[Ak] <- 0
      dk[Ik] <- (tx[Ik, ] %*% ek - alpha * bk[Ik])/gk[Ik]
      Hk <- sqrt(gk)*abs(tau*bk+(1-tau)*dk)
      Anew <- which(Hk >= sort(Hk, decreasing =TRUE)[ms])
      if(setequal(Ak, Anew)){
        return(c(k, bk))
      }
      Ak <- Anew
    }
    return(c(k, bk))
  }
  
  para <- sapply(size, function(ms){sdar(ms)})
  iter <- para[1, ]
  b <- t(para[-1, ])
  esti <- rbind(meany - meanx %*% t(b), t(b))
  
  return(list(beta=esti, size=size, iter=iter))
}

frbsr_ls <- function(x, y, alpha=0, size=NULL, h=length(y), tau=0.5, detla=0.0125,
                    intercept=TRUE, method=c("LTS", "REWLS", "RWLS")){
  x <- as.matrix(x)
  y <- as.matrix(y)
  np <- dim(x)
  n <- np[1]
  p <- np[2]
  
  if(is.null(size)) {
    size <- seq(1, ceiling(n/log(n)), 1)
  } 
  
  tx <- t(x)
  gk <- colSums(x^2) + alpha
  b0 <- ifelse(intercept, mean(y), 0)
  
  k_alpha <- function(gamma){
    a <- qnorm((gamma+1)/2)
    f_fun <- function(u){u^2 * exp(-u^2/2)/sqrt(2*pi)}
    sqrt(gamma /integrate(f_fun, -a, a)$value)
  }
  hatvalue <- function(x){
    x <- as.matrix(x)
    tx <- t(x)
    diag(x %*% solve(tx %*% x) %*%tx)
  }
  ridge <- function(x, y, alpha, intercept=TRUE){
    x <- as.matrix(x)
    y <- as.matrix(y)
    np <- dim(x)
    n <- np[1]
    p <- np[2]
    
    if(intercept){
      meanx <- colMeans(x)
      meany <- sum(y)/n
    }else{
      meanx <- rep(0, p)
      meany <- 0
    }
    
    x <- scale(x, meanx, FALSE)
    y <- y - meany
    tx <- t(x)
    
    para <- solve(tx %*% x  + diag(alpha, p)) %*% tx %*% y
    
    esti <- rbind(meany - meanx %*% para, para)
    esti
  }
  
  lts <- function(ms) {
    weight <- rep(0, n) 
    bk <- rep(0, p)
    dk <- (tx %*% (y - x %*% bk - b0) - alpha*bk)/gk
    Hk <- sqrt(gk)*abs(tau*bk+(1-tau)*dk)
    Ak <- which(Hk >= sort(Hk, decreasing =TRUE)[ms])
    fit <- ridge(x[, Ak], y, alpha, intercept)
    b0 <- fit[1]
    bk[Ak] <- fit[-1]
    res <- y - as.matrix(x[, Ak]) %*% bk[Ak] - b0
    # pos <- order(res^2)[1:h]
    pos <- order(hatvalue(x[, Ak]))[1:h]
    iter <- 0
    for (i in 1:100) {
      fit <- csdar(x[pos, ], y[pos], alpha, size=ms, beta.int=bk, tau, intercept)
      b0 <- fit$beta[1]
      bk <- fit$beta[-1]
      Ak <- which(bk!=0)
      iter <- fit$iter + iter
      res <- y - as.matrix(x[, Ak]) %*% bk[Ak] - b0
      posnew <- order(res^2)[1:h]
      if(setequal(pos, posnew)) break
      pos <- posnew
    }
    weight[pos] <- 1
    mu <- 0
    # mu <- mean(res[pos])
    gamma <- h/n
    kalpha <- k_alpha(gamma)
    sigma <- kalpha * sqrt(mean((res[pos]-mu)^2))
    hbic <- log(sigma) + log(log(n))*log(p)/n*ms
    return(c(i, iter, hbic, sigma, b0, bk, weight))
  }
  
  rewls <- function(ms) {
    weight <- rep(0, n) 
    bk <- rep(0, p)
    dk <- (tx %*% (y - x %*% bk - b0) - alpha*bk)/gk
    Hk <- sqrt(gk)*abs(tau*bk+(1-tau)*dk)
    Ak <- which(Hk >= sort(Hk, decreasing =TRUE)[ms])
    fit <- ridge(x[, Ak], y, alpha, intercept)
    b0 <- fit[1]
    bk[Ak] <- fit[-1]
    res <- y - as.matrix(x[, Ak]) %*% bk[Ak] - b0
    # pos <- order(res^2)[1:h]
    pos <- order(hatvalue(x[, Ak]))[1:h]
    iter <- 0
    for (i in 1:100) {
      fit <- csdar(x[pos, ], y[pos], alpha, size=ms, beta.int=bk, tau, intercept)
      b0 <- fit$beta[1]
      bk <- fit$beta[-1]
      Ak <- which(bk!=0)
      iter <- fit$iter + iter
      res <- y - as.matrix(x[, Ak]) %*% bk[Ak] - b0
      posnew <- order(res^2)[1:h]
      if(setequal(pos, posnew)) break
      pos <- posnew
    }
    
    mu <- 0
    # mu <- mean(res[pos])
    gamma <- h/n
    kalpha <- k_alpha(gamma)
    sigma <- kalpha * sqrt(mean((res[pos]-mu)^2))
    posnew <- which(abs(res - mu)/sigma < qnorm(1-detla))
    fit <- csdar(x[posnew, ], y[posnew], alpha, size=ms, beta.int=bk, tau, intercept)
    b0 <- fit$beta[1]
    bk <- fit$beta[-1]
    iter <- fit$iter + iter 
    weight[posnew] <- 1
    
    Aknew <- which(bk!=0)
    resnew <- y - as.matrix(x[, Aknew]) %*% bk[Aknew] - b0
    munew <- 0
    # munew <- mean(resnew[posnew])
    kalphanew <- k_alpha(length(posnew)/n)
    sigmanew <- kalphanew * sqrt(mean((resnew[posnew]-munew)^2))
    hbic <- log(sigmanew)+ log(log(n))*log(p)/n*ms
    return(c(i, iter, hbic, sigmanew, b0, bk, weight))
  }
  
  rwls <- function(ms) {
    bk <- rep(0, p)
    dk <- (tx %*% (y - x %*% bk - b0) - alpha*bk)/gk
    Hk <- sqrt(gk)*abs(tau*bk+(1-tau)*dk)
    Ak <- which(Hk >= sort(Hk, decreasing =TRUE)[ms])
    fit <- ridge(x[, Ak], y, alpha, intercept)
    b0 <- fit[1]
    bk[Ak] <- fit[-1]
    res <- y - as.matrix(x[, Ak]) %*% bk[Ak] - b0
    # pos <- order(res^2)[1:h]
    pos <- order(hatvalue(x[, Ak]))[1:h]
    iter <- 0
    for (i in 1:100) {
      fit <- csdar(x[pos, ], y[pos], alpha, size=ms, beta.int=bk, tau, intercept)
      b0 <- fit$beta[1]
      bk <- fit$beta[-1]
      Ak <- which(bk!=0)
      iter <- fit$iter + iter
      res <- y - as.matrix(x[, Ak]) %*% bk[Ak] - b0
      posnew <- order(res^2)[1:h]
      if(setequal(pos, posnew)) break
      pos <- posnew
    }
    
    mu <- 0
    # mu <- mean(res[pos])
    gamma <- h/n
    kalpha <- k_alpha(gamma)
    sigma <- kalpha * sqrt(mean((res[pos]-mu)^2))
    weight <- exp(-0.5*res^2/sigma^2)
    xnew <- as.vector(sqrt(weight)) * x 
    ynew <- sqrt(weight) * y
    fit <- csdar(xnew, ynew, alpha, size=ms, beta.int = bk, tau, intercept)
    b0 <- fit$beta[1]
    bk <- fit$beta[-1]
    iter <- fit$iter + iter 
    Aknew <- which(bk!=0)
    sigmanew <- sqrt(mean((ynew - as.matrix(xnew[, Aknew]) %*% bk[Aknew]-b0)^2))
    hbic <- log(sigmanew) + log(log(n))*log(p)/n*ms
    return(c(i, iter, hbic, sigma, fit$beta, weight))
  }
  
  solvefun <- switch(method, LTS=lts, REWLS=rewls, RWLS=rwls)
  
  para <- sapply(size, function(ms){solvefun(ms)})
  h_iter <- para[1, ]
  tol_iter <- para[2, ]
  hbic <- para[3, ]
  sigma <- para[4, ]
  esti <-  para[5:(p+5), ] 
  weight <- para[-c(1:(p+5)), ]
  
  return(list(beta=esti, size=size, hbic=hbic, sigma=sigma, weight=weight, tol_iter=tol_iter, h_iter=h_iter))
  
  
}

####################Data generation#############################################
DGPfun1_LM <- function(n=200, q=12, p=500, rho=0.5, sigma=1, R=3, pi_x=0, pi_y=0, error="d1"){
  ntrain <- n-100
  minbeta <- sigma*sqrt(2*log(p)/ntrain)
  beta <- rep(0, p)
  poi <- sample(1:p, q)
  beta[poi] <- runif(q, minbeta, R*minbeta)
  X <- matrix(, n, p)
  detla_x <- rbinom(n, 1, pi_x)
  X_normal <- matrix(rnorm(n * p, 0, 1), n, p)
  X_abnormal <- matrix(rnorm(n * p, 10, 1), n, p)
  if(rho!=0){
    corrmat <- toeplitz(rho^(0:(p-1)))
    X_normal <- X_normal  %*% chol(corrmat)
  }else{
    X_normal <- X_normal
  }
  
  X <- (1-detla_x) * X_normal + detla_x * X_abnormal
  
  detla_y <- rbinom(n, 1, pi_y)
  ee <- switch(error, d1=rcauchy(n), d2=rt(n,3), 
               d3=(1-detla_y) * rnorm(n, 0, sigma) + detla_y * rnorm(n, 10*sigma, sigma))
  
  y <- X %*% beta + ee
  Xtrain <- X[1:ntrain, ]
  Xtest <-  X[(ntrain+1):n, ]
  ytrain <- y[1:ntrain]
  ytest <- y[(ntrain+1):n]
  outlier_xtrain <- which(detla_x[1:ntrain]==1)
  outlier_ytrain <- which(detla_y[1:ntrain]==1)  
  
  return(list(xtrain=Xtrain, ytrain=ytrain, xtest=Xtest, ytest=ytest, beta=beta, 
              outlier_xtrain=outlier_xtrain, outlier_ytrain=outlier_ytrain))
}

####################Data merger#############################################
statfun.parallel_col <- function(sim.res){
  REE_mat = t(sapply(sim.res, function(u){u[, 1]}))
  APE_mat = t(sapply(sim.res, function(u){u[, 2]}))
  PDR_mat = t(sapply(sim.res, function(u){u[, 3]}))
  FDR_mat = t(sapply(sim.res, function(u){u[, 4]}))
  Time_mat = t(sapply(sim.res, function(u){u[, 5]}))

  MEAN = rbind(colMeans(REE_mat), colMeans(APE_mat), 
               colMeans(PDR_mat), colMeans(FDR_mat),
               colMeans(Time_mat))
  SD = rbind(colSds(REE_mat), colSds(APE_mat), 
             colSds(PDR_mat), colSds(FDR_mat),
              colSds(Time_mat))
  rownames(MEAN)=rownames(SD) = c("AEE", "APE", "PDR", "FDR", "Time")
  res_summary <- list(Mean=MEAN, Sd=SD)
  return(res_summary)
}

####################Simulation#############################################
library(foreach)
library(doParallel)
library(matrixStats)
N <- 100
## begin parallel
begin <- proc.time()
cores <- detectCores(logical = FALSE)
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
sim.res <- foreach(aa=1:N, .packages =c("matrixStats"),.errorhandling="remove") %dopar%
  {
  dat <- DGPfun1_LM(n=200, q=5, p=1000, rho=0.5, sigma=1, R=10, pi_x=0.1, pi_y=0.1, error="d3")
  beta <- c(0, dat$beta)
  xtrain <- dat$xtrain
  xtest <- dat$xtest
  ytrain <- dat$ytrain
  ytest <- dat$ytest
  
  n <- nrow(xtrain)
  p <- ncol(xtrain)
  size <- 1:floor(n/(log(log(n))*log(p)))
  
  intercept = TRUE
  h = floor(n*0.75)
  
  
  t1.csdar <- proc.time()
  fit.csdar <- frbsr_ls(xtrain, ytrain, size=size, h=n, tau=0.5, intercept =intercept, method="LTS")
  beta.csdar <- fit.csdar$beta[, which.min(fit.csdar$hbic)]
  t2.csdar <- proc.time()
  time.csdar <- (t2.csdar-t1.csdar)[3]
  
  
  t1.rbsr.lts <- proc.time()
  fit.rbsr.lts <- frbsr_ls(xtrain, ytrain, size=size, h=h, tau=0.5, intercept =intercept, method="LTS")
  beta.lts <- fit.rbsr.lts$beta[, which.min(fit.rbsr.lts$hbic)]
  t2.rbsr.lts <- proc.time()
  time.rbsr.lts <- (t2.rbsr.lts-t1.rbsr.lts)[3]
  
  
  t1.rbsr.rewls <- proc.time()
  fit.rbsr.rewls <- frbsr_ls(xtrain, ytrain, size=size, h=h, tau=0.5, intercept =intercept, method="REWLS")
  beta.rewls <- fit.rbsr.rewls$beta[, which.min(fit.rbsr.rewls$hbic)]
  t2.rbsr.rewls <- proc.time()
  time.rbsr.rewls <- (t2.rbsr.rewls-t1.rbsr.rewls)[3]
  
  
  t1.rbsr.rwls <- proc.time()
  fit.rbsr.rwls <- frbsr_ls(xtrain, ytrain, size=size, h=h, tau=0.5, intercept =intercept, method="RWLS")
  beta.rwls <- fit.rbsr.rwls$beta[, which.min(fit.rbsr.rwls$hbic)]
  t2.rbsr.rwls <- proc.time()
  time.rbsr.rwls <- (t2.rbsr.rwls-t1.rbsr.rwls)[3]
  
  
  pos <- which(beta!=0)
  Time <- c(time.csdar, time.rbsr.lts, time.rbsr.rewls, time.rbsr.rwls)
  
  REE <- c(norm(beta.csdar-beta,"2"),  norm(beta.lts-beta,"2"), 
           norm(beta.rewls-beta,"2"), norm(beta.rwls-beta,"2"))/norm(beta, "2")
  
  APE <- c(mean(abs(cbind(1, xtest)%*%beta.csdar-ytest)), 
           mean(abs(cbind(1, xtest)%*%beta.lts-ytest)),
           mean(abs(cbind(1, xtest)%*%beta.rewls-ytest)),
           mean(abs(cbind(1, xtest)%*%beta.rwls-ytest)))
  
  PDR <- c(length(which(beta.csdar[pos]!=0)), length(which(beta.lts[pos]!=0)),
           length(which(beta.rewls[pos]!=0)), length(which(beta.rwls[pos]!=0)))/length(pos)
  
  FDR <- c(length(which(beta.csdar[-c(1,pos)]==0)), length(which(beta.lts[-c(1,pos)]==0)),
           length(which(beta.rewls[-c(1,pos)]==0)), length(which(beta.rwls[-c(1,pos)]==0)))/(length(beta)-length(pos)-1)
  
  result <- cbind(REE, APE, PDR, FDR, Time)
  rownames(result)<-c("CSDAR", "LTS", "REWLS", "RWLS")
  result
}
stopImplicitCluster()
stopCluster(cl)

end <- proc.time()
end-begin

res_summary = statfun.parallel_col(sim.res)
res_summary$Mean
