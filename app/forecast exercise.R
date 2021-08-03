## forecasting exercise
library(mosumvar); library(factorcpt)


##
## Most recent change point -------------



forecast_mosumvar <- mosumvar_factor(data_mosumvar[1:200,], p=1, G=24, nu=.3, criterion = "eta", method = "Score")
time(data_mosumvar)[forecast_mosumvar$cps]
cp <- max(forecast_mosumvar$cps, 1)

# MTS::MTSplot( scale(data_mosumvar)[1:200,])
# abline(v = forecast_mosumvar$cps, col = "red")
sample_end <- 190
forecast_out <- matrix(0, nrow = 10, ncol = 4)
for (t in 1:10) {
  fm <- get.factor.model(t(data_mosumvar[1:(sample_end+t),]), q = r_)
  factors_ <- t(fm$f[1:r_,])
  
  full <- ar(factors_, aic = F, order.max = 1) ##full window
  #split_oracle <- ar(factors_[break_point:150,], aic = F, order.max = 1) ##post-break window
  split <- ar(factors_[cp:(sample_end+t),], aic = F, order.max = 1) ##post-break window
  #split_window_oracle <- pooled_forecast(factors_, cp = break_point, p =1) ## pooled
  split_window <- pooled_forecast(factors_, cp = cp, p =1, window_size = 50) ## pooled
  split_window_all <- pooled_forecast(factors_, cp = cp, p =1) ## pooled
  
  pred_full <- t(full$ar[,,] %*%  factors_[nrow(factors_),]) %*% t(fm$lam[,1:r_])
  #pred_oracle <-t(split_oracle$ar[,,] %*%   factors_[nrow(data3),]) %*% t(fm$lam[,1:r_])
  pred_split <-t(split$ar[,,] %*%   factors_[nrow(factors_),]) %*% t(fm$lam[,1:r_])
  #pred_window_oracle <-  as.numeric(split_window_oracle %*% t(fm$lam[,1:r_]) )
  pred_window <-  as.numeric(split_window %*% t(fm$lam[,1:r_]) )
  pred_window_all <-  as.numeric(split_window_all %*% t(fm$lam[,1:r_]) )
  
  onestep <- data_mosumvar[190+t,]
  sum_onestep <- sum(onestep^2)
  forecast_out[t,1] <- sum( (pred_full-onestep)^2 )/sum_onestep
  #forecast_out[t,2] <- sum( (pred_oracle-onestep)^2 )/sum_onestep
  forecast_out[t,2] <- sum( (pred_split-onestep)^2 )/sum_onestep
  #forecast_out[t,4] <- sum( (pred_window_oracle-onestep)^2 )/sum_onestep
  forecast_out[t,3] <- sum( (pred_window-onestep)^2 )/sum_onestep
  forecast_out[t,4] <- sum( (pred_window_all-onestep)^2 )/sum_onestep
}
colMeans(forecast_out) ##marginal improvement using adaptive window
apply(forecast_out, 2, sd)


##
## Second change point -------------

cp2 <- forecast_mosumvar$cps[2]

# MTS::MTSplot( scale(data_mosumvar)[1:200,])
# abline(v = forecast_mosumvar$cps, col = "red")
sample_end2 <- 140
forecast_out2 <- matrix(0, nrow = 17, ncol = 4)
for (t in 1:17) {
  fm <- get.factor.model(t(data_mosumvar[1:(sample_end2+t),]), q = r_)
  factors_ <- t(fm$f[1:r_,])
  
  full <- ar(factors_, aic = F, order.max = 1) ##full window
  #split_oracle <- ar(factors_[break_point:150,], aic = F, order.max = 1) ##post-break window
  split <- ar(factors_[cp2:(sample_end2+t),], aic = F, order.max = 1) ##post-break window
  #split_window_oracle <- pooled_forecast(factors_, cp = break_point, p =1) ## pooled
  split_window <- pooled_forecast(factors_, cp = cp2, p =1, window_size = 50) ## pooled
  split_window_all <- pooled_forecast(factors_, cp = cp2, p =1) ## pooled
  
  pred_full <- t(full$ar[,,] %*%  factors_[nrow(factors_),]) %*% t(fm$lam[,1:r_])
  #pred_oracle <-t(split_oracle$ar[,,] %*%   factors_[nrow(data3),]) %*% t(fm$lam[,1:r_])
  pred_split <-t(split$ar[,,] %*%   factors_[nrow(factors_),]) %*% t(fm$lam[,1:r_])
  #pred_window_oracle <-  as.numeric(split_window_oracle %*% t(fm$lam[,1:r_]) )
  pred_window <-  as.numeric(split_window %*% t(fm$lam[,1:r_]) )
  pred_window_all <-  as.numeric(split_window_all %*% t(fm$lam[,1:r_]) )
  
  onestep <- data_mosumvar[sample_end2+t,]
  sum_onestep <- sum(onestep^2)
  forecast_out2[t,1] <- sum( (pred_full-onestep)^2 )/sum_onestep
  #forecast_out[t,2] <- sum( (pred_oracle-onestep)^2 )/sum_onestep
  forecast_out2[t,2] <- sum( (pred_split-onestep)^2 )/sum_onestep
  #forecast_ou2t[t,4] <- sum( (pred_window_oracle-onestep)^2 )/sum_onestep
  forecast_out2[t,3] <- sum( (pred_window-onestep)^2 )/sum_onestep
  forecast_out2[t,4] <- sum( (pred_window_all-onestep)^2 )/sum_onestep
}
colMeans(forecast_out2) ##marginal improvement using adaptive window
apply(forecast_out2, 2, sd)
