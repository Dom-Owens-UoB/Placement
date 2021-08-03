## mosumvar_factor sims
library(mosumvar); library(factorcpt); set.seed(222)
N_sim <- 500
A_sim1 <- fcpt_ar$ar[,,] 
A_sim2 <- fcpt_ar2$ar[,,]
A_sim3 <- fcpt_ar3$ar[,,]

par(mfrow = c(1,3))
image(A_sim1); image(A_sim2); image(A_sim3);
par(mfrow = c(1,1))

###
###  ---------------------------------------------------
###

mosumvar_factor_sim <- function(N_sim=100,r=6,d=150, r_known = T){
  if(r_known) r_ <- r else r_ <- NULL
  out <- matrix(NA, nrow = N_sim, ncol = 7)
  colnames(out) <- c("nouniv.H","univ.H","factorcpt.H","nouniv.q","univ.q","factorcpt.q","r")
  true_cps <- c(100,200)
  for (sim in 1:N_sim) {
    data1 <- mosumvar::VAR.sim(100, coeffs = A_sim1[1:r,1:r])  
    data2 <- mosumvar::VAR.sim(100, coeffs = A_sim2[1:r,1:r])
    data3 <- mosumvar::VAR.sim(100, coeffs = A_sim3[1:r,1:r])
    factors_ <- rbind(data1,data2,data3)
    loadings <- matrix(rnorm(d*r, 0,1) , nrow = r, ncol = d)
    data_ <- factors_ %*% loadings
    nouniv <- mosumvar_factor(data_,p=1,G=50,r=r_, method = "Score", alpha = .2, criterion = "eta")
    univ <- mosumvar_factor(data_,p=1,G=50,r=r_, univ = T, method = "Score", alpha = .2, criterion = "eta")
    fcpt <- factorcpt::factor.seg.alg(t(data_), r= r_, q.seq = 4:8, B = 10, rule = 2, sig.lev = .2)
    if(!is.null(nouniv$cps)) out[sim,1] <- pracma::hausdorff_dist(nouniv$cps, true_cps) / 300 else out[sim,1] <- 1
    if(!is.null(univ$cps)) out[sim,2] <- pracma::hausdorff_dist(univ$cps, true_cps) / 300 else out[sim,2] <-1
    if(!is.null(fcpt$common.est.cps)) out[sim,3] <- pracma::hausdorff_dist(fcpt$common.est.cps, true_cps) / 300 else out[sim,3] <-1
    out[sim,4] <- length(nouniv$cps)
    out[sim,5] <- length(univ$cps)
    out[sim,6] <- length(fcpt$common.est.cps)
    out[sim,7] <- nouniv$r
  }
  return(out)
}
sim1_r_unknown <- mosumvar_factor_sim(100, r_known = F)
sim1_r_unknown[sim1_r_unknown[,3] == Inf, 3] <- NA
colMeans(sim1_r_unknown, T)
sim1_r_known <- mosumvar_factor_sim(100, r_known = T)
sim1_r_known[sim1_r_known[,3] == Inf, 3] <- NA
colMeans(sim1_r_known, T)


# Pesaran/Timmerman 2007

mosumvar_factor_sim_forecast <- function(N_sim=100,r=6,d=150, break_point = 75, A_2 = A_sim3, r_known = T){
  if(r_known) r_ <- r else r_ <- NULL
  out <- matrix(NA, nrow = N_sim, ncol = 5)
  colnames(out) <- c("full","split_oracle","split","pooled_oracle","pooled")
  for (sim in 1:N_sim) {
    data1 <- mosumvar::VAR.sim(break_point, coeffs = A_sim1[1:r,1:r])  
    #data2 <- mosumvar::VAR.sim(100, coeffs = A_sim2[1:r,1:r])
    data3 <- mosumvar::VAR.sim(150 - break_point, coeffs = A_2[1:r,1:r])
    loadings <- matrix(rnorm(d*r, 0,1) , nrow = r, ncol = d)
    onestep <- t(A_2[1:r,1:r] %*% data3[nrow(data3),]) %*% loadings
    factors_ <- rbind(data1,data3)
    data_ <- factors_ %*% loadings
    nouniv <- mosumvar_factor(data_,p=1,G=40,r=r_, method = "Score", nu=.1, alpha = .95)
    cp <- max(nouniv$cps, 1)
    fm <- get.factor.model(t(data_), q = r_)
    factors_ <- t(fm$f[1:r_,])
    
    full <- ar(factors_, aic = F, order.max = 1) ##full window
    split_oracle <- ar(factors_[break_point:150,], aic = F, order.max = 1) ##post-break window
    split <- ar(factors_[cp:150,], aic = F, order.max = 1) ##post-break window
    split_window_oracle <- pooled_forecast(factors_, cp = break_point, p =1) ## pooled
    split_window <- pooled_forecast(factors_, cp = cp, p =1) ## pooled
    
    pred_full <- t(full$ar[,,] %*%  factors_[nrow(data3),]) %*% t(fm$lam[,1:r_])
    pred_oracle <-t(split_oracle$ar[,,] %*%   factors_[nrow(data3),]) %*% t(fm$lam[,1:r_])
    pred_split <-t(split$ar[,,] %*%   factors_[nrow(data3),]) %*% t(fm$lam[,1:r_])
    pred_window_oracle <-  as.numeric(split_window_oracle %*% t(fm$lam[,1:r_]) )
    pred_window <-  as.numeric(split_window %*% t(fm$lam[,1:r_]) )
    
    sum_onestep <- sum(onestep^2)
    out[sim,1] <- sum( (pred_full-onestep)^2 )/sum_onestep
    out[sim,2] <- sum( (pred_oracle-onestep)^2 )/sum_onestep
    out[sim,3] <- sum( (pred_split-onestep)^2 )/sum_onestep
    out[sim,4] <- sum( (pred_window_oracle-onestep)^2 )/sum_onestep
    out[sim,5] <- sum( (pred_window-onestep)^2 )/sum_onestep
    #   out[sim,3] <- pracma::hausdorff_dist(fcpt$common.est.cps, true_cps) / 300
  }
  return(out)
}
sim1_forecast_A3 <- mosumvar_factor_sim_forecast(100, r =6)
sim1_forecast_A2 <- mosumvar_factor_sim_forecast(100, r =6, A_2 = A_sim2)

sim1_forecast_A3_cp30 <- mosumvar_factor_sim_forecast(100, r =6, break_point = 30)
sim1_forecast_A2_cp30 <- mosumvar_factor_sim_forecast(100, r =6, A_2 = A_sim2, break_point = 30)

sim1_forecast_A3_cp100 <- mosumvar_factor_sim_forecast(100, r =6, break_point = 100)
sim1_forecast_A2_cp100 <- mosumvar_factor_sim_forecast(100, r =6, A_2 = A_sim2, break_point = 100)

###
### New factor ---------------------------------------------------
###

mosumvar_factor_sim_newfactor <- function(N_sim=100, r2 =4,  d=150, r_known = T){
  if(r_known) r_ <- r2 else r_ <- NULL
  #r <- max(r1,r2)
  out <- matrix(NA, nrow = N_sim, ncol = 3)
  true_cps <- c(100,200)
  A_sim1_r <- A_sim1[1:r2,1:r2]
  A_sim1_r[r2,] <- A_sim1_r[,r2] <- 0 
  A_sim3_r <- A_sim3[1:r2,1:r2]
  A_sim3_r[r2,] <- A_sim3_r[,r2] <- 0 
  for (sim in 1:N_sim) {
    data1 <- mosumvar::VAR.sim(100, coeffs = A_sim1_r)  
    data2 <- mosumvar::VAR.sim(100, coeffs = A_sim2[1:r2,1:r2])
    data3 <- mosumvar::VAR.sim(100, coeffs = A_sim3_r)
    factors_ <- rbind(data1,data2,data3)
    loadings <- matrix(rnorm(d*r2, 0,1) , nrow = r2, ncol = d)
    data_ <- factors_ %*% loadings
    nouniv <- mosumvar_factor(data_,p=1,G=50,r=r_, method = "Score")
    univ <- mosumvar_factor(data_,p=1,G=50,r=r_, univ = T, method = "Score")
    #    fcpt <- factorcpt::factor.seg.alg(t(data_), r= r_)
    if(!is.null(nouniv$cps)) out[sim,1] <- pracma::hausdorff_dist(nouniv$cps, true_cps) / 300 else out[sim,1] <- 1
    if(!is.null(univ$cps)) out[sim,2] <- pracma::hausdorff_dist(univ$cps, true_cps) / 300 else out[sim,2] <-1
    #   out[sim,3] <- pracma::hausdorff_dist(fcpt$common.est.cps, true_cps) / 300
  }
  return(out)
}
sim2_r_unknown <- mosumvar_factor_sim_newfactor(500, r_known = F)
sim2_r_known <- mosumvar_factor_sim_newfactor(500, r_known = T)
