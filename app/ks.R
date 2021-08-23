## kalman smoother
library("dlm")
library("FKF")
library("factorcpt")

fkf_data <- cbind(mosumvar_gdp_window,data_mosumvar)

for (t in (sample_start+1):(sample_end-6) ) { #to end of regime
  t = sample_end-6
  window_start <- max(cps_3[t-18>cps_3])
  mosumvar_gdp_window_ks <- window(mosumvar_gdp ,time(data_mosumvar)[1],time(data_mosumvar)[t] ) 
  fm_ks <- get.factor.model(t(data_mosumvar[1:t,]), q = r_)
  #data_window_ks <- window(fkf_data ,time(fkf_data)[1],time(fkf_data)[t] )
  factors_ks <- t(fm_ks$f[1:r_,]) 
  gdp_lm_ks <- lm(mosumvar_gdp_window_ks~ factors_ks)
  window_robust <- pooled_forecast(factors_ks[,], cp = window_start-1, window_size = 30, p =1, n.ahead = 6, weights = "robust", return_params = T)
  #dlm::dlmModARMA(ar= list( as.matrix(window_robust$params[1,,]) ) , ma = list(diag(0, r_) ),  dV = 0.0001) #diag(window_robust$var)
  #, m0 = rep(0, r_), C0 = diag(0, r_)
  a0 <- colMeans(factors_ks)
  P0 <- cov(factors_ks)
  lm_resids <- rep(NA,t)
  lm_resids[as.numeric(names(gdp_lm_ks$residuals))] <- gdp_lm_ks$residuals
  ks_var <- #cov( cbind(gdp_lm_ks$residuals, (data_mosumvar[1:t,]) - t(fm_ks$lam[,1:r_] %*% fm_ks$f[1:r_,])), use = "complete.obs")
    diag(c(var(gdp_lm_ks$residuals, na.rm = T), diag(cov( (data_mosumvar[1:t,]) - t(fm_ks$lam[,1:r_] %*% fm_ks$f[1:r_,]) ) ) ) )
  Zt <- t(as.matrix(cbind( gdp_lm_ks$coefficients[-1], t(fm$lam[,1:r_]) ) ) )
  ct <- as.matrix(0*Zt[,1])
  ct[1,] <- gdp_lm_ks$coefficients[1]
  fkf_t <- FKF::fkf(a0 = a0, P0 = P0, dt = matrix(0, nrow=r_), ct = ct, Tt =window_robust$params[1,,], #gdp_lm_ks$coefficients[1]
           Zt = Zt, HHt = (window_robust$var), GGt = ks_var,#(diag(1,148)),
           #as.matrix(var(gdp_lm_ks$residuals, na.rm = T)), #matrix(gdp_lm_ks$coefficients[-1], nrow = 1, ncol = r_)
           yt =  as.matrix( t(fkf_data[1:t,])) )#matrix(mosumvar_gdp_window_ks, nrow = 1) )
  
}
  plot.ts(mosumvar_gdp_window[!is.na(mosumvar_gdp_window)])
  ts.plot(fkf_t$att[1,])
  