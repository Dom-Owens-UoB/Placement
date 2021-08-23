## forecast quarterly data only
## t = 1, 4, 7, ... 
library(factorcpt)

getF.5 <- function(pred, true, thr){
  tp <- sum( (pred > thr) & true,na.rm = T )
  fn <- sum( !(pred > thr) & true,na.rm = T )
 fd <- sum( (pred > thr) & !true , na.rm =T)
  tn <- sum( !(pred > thr) & !true, na.rm =T )
 
  # precision <- tp/(tp+fd)
  # recall <- tp/(tp+fn)
#  1.25 * (precision*recall)/(0.25 * precision + recall)
  return( 1.25*tp^2/(.25*tp*(tp+fn) + tp*(tp+fd)) )
}
getprecision <- function(pred, true, thr){
  tp <- sum( (pred > thr) & true,na.rm = T )
  fn <- sum( !(pred > thr) & true,na.rm = T )
  fd <- sum( (pred > thr) & !true , na.rm =T)
  tn <- sum( !(pred > thr) & !true, na.rm =T )
  
   precision <- tp/(tp+fd)
  # recall <- tp/(tp+fn)
  #  1.25 * (precision*recall)/(0.25 * precision + recall)
  return(precision )
}
getinvprecision <- function(pred, true, thr){
  tp <- sum( (pred > thr) & true,na.rm = T )
  fn <- sum( !(pred > thr) & true,na.rm = T )
  fd <- sum( (pred > thr) & !true , na.rm =T)
  tn <- sum( !(pred > thr) & !true, na.rm =T )
  
  invprecision <- tn/(tn+fn)
  # recall <- tp/(tp+fn)
  #  1.25 * (precision*recall)/(0.25 * precision + recall)
  return(invprecision )
}

##
## forecast up to/after first change point -------------


# 
# forecast_mosumvar <- mosumvar_factor(data_mosumvar[1:200,], p=1, G=24, nu=.3, criterion = "eta", method = "Score")
# time(data_mosumvar)[forecast_mosumvar$cps]
# cp <- max(forecast_mosumvar$cps, 1) ## 173
cp1 <- forecast_mosumvar$cps[1] ## 53 
cp2 <- forecast_mosumvar$cps[2] ##104
h_step <- matrix(0, 6, r_)


fm <- get.factor.model(t(data_mosumvar[1:200,]), q = r_)#[1:(30+t),]
factors_ <- t(fm$f[1:r_,])

## gdp lm -------------
time(data_mosumvar)[1]
mosumvar_gdp_window <- window(mosumvar_gdp ,time(data_mosumvar)[1],time(data_mosumvar)[200] )
ts.plot(mosumvar_gdp_window[!is.na(mosumvar_gdp_window)])
gdp_lm <- lm(mosumvar_gdp_window~ factors_)



pred_factor <- t(factors_[30+(t):(t),])
 
# gdp_fill_in <- cbind(1,(factors_)) %*%  (gdp_lm$coefficients)#predict.lm(gdp_lm, h_step) 
# gdp_filled_in <- mosumvar_gdp_window
# gdp_filled_in[is.na(gdp_filled_in)] <- gdp_fill_in[is.na(gdp_filled_in)]

#gdp_filled_diff <- diff(gdp_filled_in)

## gdp glm ----------
bin_response <- (mosumvar_gdp_window<0) ## is this recession?
bin_response_trim <- bin_response[!is.na(bin_response)]
plot.ts(bin_response_trim)
#wt <- (1-bin_response)*9 + 1
gdp_glm <- glm( bin_response ~ factors_[,], family = "binomial", weights = NULL)
plot.ts(gdp_glm$fitted.values)



bin_score <- 1:100 * 0
accuracy <- 1:100 * 0
F.5 <- precision <- recall <- tpr <- fdr <- 1:100 * 0
for (ii in 1:100 ) {
  bin_score[ii] <- sum(log(gdp_glm$fitted.values[(gdp_glm$fitted.values > ii* 0.01) == bin_response_trim ])) #log score rule
  accuracy[ii] <- mean( (gdp_glm$fitted.values > ii* 0.01) == bin_response_trim )
  tp<- sum( (gdp_glm$fitted.values > ii* 0.01) & bin_response_trim )
  fn <- sum( !(gdp_glm$fitted.values > ii* 0.01) & bin_response_trim )
  fd <- sum( (gdp_glm$fitted.values > ii* 0.01) & !bin_response_trim )
  tn <- sum( !(gdp_glm$fitted.values > ii* 0.01) & !bin_response_trim )
 tpr[ii] <-  tp/(tp+fn)
 fdr[ii] <-  fd/(tp+tn)
 precision[ii] <- tp/(tp+fd)
 recall[ii] <- tp/(tp+fn)
 invprecision[ii] <- tn/(tn+fn)
 F.5[ii] <- 1.25 * (precision[ii]*recall[ii])/(0.25 * precision[ii] + recall[ii])
}
plot(bin_score)
plot(fdr, tpr, type = "l"); abline(0,1, lty = 3) #roc
plot(1:100 * .01 , F.5, type = "l", xlab = "Threshold"); abline(v = which.max(F.5)*0.01, col = "red") ##34
class_thresh <- which.max(F.5)*0.01#which.max(precision)*0.01#(F.5)*0.01
class_thresh_invprecision <- which.max(invprecision)*0.01
# which.max(bin_score) * 0.01
# which.max(accuracy) * 0.01 
##







####################################
## first segment -------------------
####################################
q_forecast_out <- matrix(0, nrow = 30, ncol = 6)
q_gdp_lm_out <- q_gdp_glm_out <- q_gdp_glm_true <- matrix(0, nrow = 30, ncol = 7)#q_gdp_glm_out <-  q_forecast_out
for (t in 1:30 ) { #to end of regime
  full <- ar(factors_[1:(29+t),], aic = F, order.max = 1) ##full window
  # gdp_lm <- lm(mosumvar_gdp_window[1:(30+t)] ~ factors_[1:(30+t),])
  # gdp_glm <- glm( (mosumvar_gdp_window[1:(30+t)]>0) ~ factors_[1:(30+t),], family = "binomial")
  
  pred_data <- data_mosumvar[30+(t):(t+5),] 
  pred_factor <- t(factors_[30+(t):(t),])
  
  h_step <- predict(full, pred_factor, n.ahead = 6, se.fit = F) 
  pred_full <-  h_step %*% t(fm$lam[,1:r_])
  lm_pred <- cbind(1,h_step) %*%  (gdp_lm$coefficients)#predict.lm(gdp_lm, h_step) 
  glm_pred <- cbind(1,h_step) %*%  (gdp_glm$coefficients)#predict.glm(gdp_glm, h_step)
  
  sum_pred_data <- rowSums(pred_data^2)
  
  for (h in 1:6) {
    q_forecast_out[t,h] <- sum( (pred_full[h,] - pred_data[h,] )^2 )/sum_pred_data[h]
  }
  q_gdp_lm_out[t,] <- c( gdp_lm$residuals[paste(29+t)] , (lm_pred - mosumvar_gdp_window[30+(t):(t+5)]) )^2
  q_gdp_glm_out[t,] <- pracma::sigmoid(c(gdp_lm$fitted.values[paste(29+t)],glm_pred)) #abs(c( pracma::sigmoid(gdp_lm$fitted.values[paste(29+t)]) - (mosumvar_gdp_window[(29+t)]>class_thresh),
                      #        pracma::sigmoid(glm_pred) - (mosumvar_gdp_window[30+(t):(t+5)]< class_thresh) ) ) ## using thresh
  q_gdp_glm_true[t,] <- mosumvar_gdp_window[29+(t):(t+6)]>0
  
}
colMeans(q_forecast_out) 
colMeans(q_forecast_out[1:17,]) #up to cp 
colMeans(q_forecast_out[18:30,]) #after cp 
# apply(q_forecast_out, 2, sd)
plot(colMeans(q_forecast_out), type = "b", xlab = "Horizon", ylab = "X Error (normalised)", col = "blue" )
plot(colMeans(q_gdp_lm_out,T), type = "b", xaxt = "n", xlab = "Horizon", ylab = "Y Error", col = "red" ); axis(1, at = 1:7, labels = 0:6)


colMeans(q_gdp_lm_out ,T)
colMeans(q_gdp_glm_out ,T)
q_F.5 <- rep(0,7)
for (h in 1:7) {
  q_F.5[h] <- getF.5(q_gdp_glm_out[,h], q_gdp_glm_true[,h], class_thresh)
}
plot(q_F.5, type = "b", xaxt = "n", xlab = "Horizon", ylab = "F0.5 Score", col = "purple" ); axis(1, at = 1:7, labels = 0:6)

q_precision <- rep(0,7)
for (h in 1:7) {
  q_precision[h] <- getprecision(q_gdp_glm_out[,h], q_gdp_glm_true[,h], class_thresh)
}
q_precision


#############################################
## up to/ after second change point ------
#############################################
cp3 <- forecast_mosumvar$cps[3] ## 157



q_forecast_out2_full <- q_forecast_out2_window <-q_forecast_out2_window_all <- matrix(0, nrow = 30, ncol = 6)

q_gdp_lm_out2 <- q_gdp_lm_out2_window  <- q_gdp_lm_out2_window_all  <-   matrix(0, nrow = 30, ncol = 7)
q_gdp_glm_out2   <- q_gdp_glm_out2_window <- q_gdp_glm_out2_window_all  <- q_gdp_glm_true2 <- q_gdp_lm_out2
# fm <- get.factor.model(t(data_mosumvar[1:(cp2+t),]), q = r_)
# factors_ <- t(fm$f[1:r_,])
for (t in 1:(30) ) { #to end of regime

  full <- ar(factors_[(1):(cp1+29+t),], aic = F, order.max = 1) ## all
  window <- ar(factors_[cp1+(1):(29+t),], aic = F, order.max = 1) ##this segment
  window_all <- pooled_forecast(factors_[(1):(cp1+29+t),], cp = cp1, p =1, n.ahead = 6) ## pooled


  pred_data <- data_mosumvar[cp1+30+(t):(t+5),] 
  pred_factor <- t(factors_[cp1+30+(t):(t),])
  
  h_step <- predict(full, pred_factor, n.ahead = 6, se.fit = F) 
  pred_full <-  h_step %*% t(fm$lam[,1:r_])
  h_step_window <- predict(window, pred_factor, n.ahead = 6, se.fit = F) 
  pred_window <-  h_step_window %*% t(fm$lam[,1:r_])
  pred_window_all <- window_all %*% t(fm$lam[,1:r_])
  
  sum_pred_data <- rowSums(pred_data^2)
  
  for (h in 1:6) {
    q_forecast_out2_full[t,h] <- sum( (pred_full[h,] - pred_data[h,] )^2 )/sum_pred_data[h]
    q_forecast_out2_window[t,h] <- sum( (pred_window[h,] - pred_data[h,] )^2 )/sum_pred_data[h]
    q_forecast_out2_window_all[t,h] <- sum( (pred_window_all[h,] - pred_data[h,] )^2 )/sum_pred_data[h]
  }
  
  lm_pred <- cbind(1,h_step) %*%  (gdp_lm$coefficients)#predict.lm(gdp_lm, h_step) 
  glm_pred <- cbind(1,h_step) %*%  (gdp_glm$coefficients)#predict.glm(gdp_glm, h_step)
  
  lm_pred_window <- cbind(1,h_step_window ) %*%  (gdp_lm$coefficients)#predict.lm(gdp_lm, h_step) 
  glm_pred_window <- cbind(1,h_step_window) %*%  (gdp_glm$coefficients)#predict.glm(gdp_glm, h_step)
  
  lm_pred_window_window_all <- cbind(1,window_all) %*%  (gdp_lm$coefficients)
  glm_pred_window_all <- cbind(1,window_all) %*%  (gdp_glm$coefficients)#predict.glm(gdp_glm, h_step)
  
  q_gdp_lm_out2[t,] <- c( gdp_lm$residuals[paste(cp1+29+t)] , (lm_pred - mosumvar_gdp_window[cp1+30+(t):(t+5)]))^2
  q_gdp_glm_out2[t,] <- pracma::sigmoid(c(gdp_lm$fitted.values[paste(cp1+29+t)],glm_pred )) #abs(c( pracma::sigmoid(gdp_lm$fitted.values[paste(29+t)]) - (mosumvar_gdp_window[(29+t)]>class_thresh),
  #        pracma::sigmoid(glm_pred) - (mosumvar_gdp_window[30+(t):(t+5)]< class_thresh) ) ) ## using thresh
  q_gdp_glm_true2[t,] <- mosumvar_gdp_window[cp1+29+(t):(t+6)]>0
  
  
  q_gdp_lm_out2_window[t,] <- c( gdp_lm$residuals[paste(cp1+29+t)] ,(lm_pred_window - mosumvar_gdp_window[cp1+30+(t):(t+5)]) )^2
  q_gdp_glm_out2_window[t,] <- pracma::sigmoid(c(gdp_lm$fitted.values[paste(cp1+29+t)],glm_pred_window))
    #abs(c( pracma::sigmoid(gdp_lm$fitted.values[paste(cp1+29+t)]) - (mosumvar_gdp_window[(cp1+29+t)]>class_thresh),
     #                                 pracma::sigmoid(glm_pred_window) - (mosumvar_gdp_window[cp1+30+(t):(t+5)]<class_thresh) ) ) 
  #abs(pracma::sigmoid(glm_pred_window) - (mosumvar_gdp_window[cp1+30+(t):(t+5)]>0) )
  
  q_gdp_lm_out2_window_all[t,] <- c( gdp_lm$residuals[paste(cp1+29+t)] ,(lm_pred_window_window_all - mosumvar_gdp_window[cp1+30+(t):(t+5)]) )^2
  q_gdp_glm_out2_window_all[t,] <- pracma::sigmoid(c(gdp_lm$fitted.values[paste(cp1+29+t)],glm_pred_window_all))
}
colMeans(q_forecast_out2_full) 
colMeans(q_forecast_out2_window) 
colMeans(q_forecast_out2_window_all) 


q_F.5_2_full <- q_F.5_2_window <- q_F.5_2_window_all <- rep(0,7)
for (h in 1:7) {
  q_F.5_2_full[h] <- getF.5(q_gdp_glm_out2[,h], q_gdp_glm_true2[,h], class_thresh)
  q_F.5_2_window[h] <- getF.5(q_gdp_glm_out2_window[,h], q_gdp_glm_true2[,h], class_thresh)
  q_F.5_2_window_all[h] <- getF.5(q_gdp_glm_out2_window_all[,h], q_gdp_glm_true2[,h], class_thresh)
}


q_precision_2_full <- q_precision_2_window <- q_precision_2_window_all <- rep(0,7)
for (h in 1:7) {
  q_precision_2_full[h] <- getprecision(q_gdp_glm_out2[,h], q_gdp_glm_true2[,h], class_thresh)
  q_precision_2_window[h] <- getprecision(q_gdp_glm_out2_window[,h], q_gdp_glm_true2[,h], class_thresh)
  q_precision_2_window_all[h] <- getprecision(q_gdp_glm_out2_window_all[,h], q_gdp_glm_true2[,h], class_thresh)
}



## plots #############################################
plot(colMeans(q_forecast_out2_full), type = "b", xlab = "Horizon", ylab = "X Error (normalised)", ylim = c(1,1.004), col = "blue" );
lines(colMeans(q_forecast_out2_window) , type = "b", col = "lightblue");
lines(colMeans(q_forecast_out2_window_all) , type = "b", col = "darkblue");
legend(x="topleft", legend = c("Full","Current","Pooled"), fill = c("blue","lightblue","darkblue"))
 
plot(colMeans(q_gdp_lm_out2,T), type = "b", xaxt = "n", xlab = "Horizon", ylim = range(colMeans(q_gdp_lm_out2_window,T)), ylab = "Y Error", col = "red" ); axis(1, at = 1:7, labels = 0:6);
lines(colMeans(q_gdp_lm_out2_window,T) , type = "b", col = "pink");lines(colMeans(q_gdp_lm_out2_window_all,T) , type = "b", col = "darkred");
legend(x="topleft", legend = c("Full","Current","Pooled"), fill = c("red","pink","darkred"))

plot( (q_F.5_2_full ), type = "b", xaxt = "n", xlab = "Horizon",ylim = c(0.0,1),   ylab = "F0.5 Score", col = "yellow" ); axis(1, at = 1:7, labels = 0:6)
lines( (q_F.5_2_window) , type = "b", col = "green");
lines( (q_F.5_2_window_all) , type = "b", col = "darkgreen");
legend(x="bottomright", legend = c("Full","Current","Pooled"), fill = c("yellow","green","darkgreen"))

plot( (q_precision_2_full ), type = "b", xaxt = "n", xlab = "Horizon",ylim = c(0.0,1),   ylab = "Precision", col = "yellow" ); axis(1, at = 1:7, labels = 0:6)
lines( (q_precision_2_window) , type = "b", col = "green");
lines( (q_precision_2_window_all) , type = "b", col = "darkgreen");
legend(x="bottomright", legend = c("Full","Current","Pooled"), fill = c("yellow","green","darkgreen"))


#############################################
## entire sample ----------------------------
#############################################
sample_start <- 31
sample_end <- 200

q_forecast_out3_full <- q_forecast_out3_window <-q_forecast_out3_window_all <-q_forecast_out3_window_robust <- q_forecast_out3_roll_window <- matrix(0, nrow = sample_end-sample_start-6, ncol = 6)
q_square_full <- q_square_window <- q_square_window_all <- q_square_window_robust <- q_square_window_robust <- q_square_roll_window <- q_forecast_out3_full

q_square_roll50 <- q_abs_roll50 <- q_square_roll18 <- q_abs_roll18 <- q_square_ensemble <- q_abs_ensemble <- q_square_full

q_gdp_lm_out3 <- q_gdp_lm_out3_window  <- q_gdp_lm_out3_window_all  <- q_gdp_lm_out3_roll_window <- q_gdp_lm_out3_window_robust <-   matrix(0, nrow = sample_end-sample_start-6, ncol = 7)
q_gdp_glm_out3   <- q_gdp_glm_out3_window <- q_gdp_glm_out3_window_all  <- q_gdp_glm_true3 <- q_gdp_glm_out3_roll_window <-  q_gdp_lm_out3 


# fm <- get.factor.model(t(data_mosumvar[1:(cp2+t),]), q = r_)
# factors_ <- t(fm$f[1:r_,])
lm_pred_1step <- lm_pred_window_1step <- lm_pred_window_all_1step <- lm_pred_window_robust_1step <- lm_pred_roll_window_1step<- matrix(0, nrow= sample_end-sample_start-6, ncol = 6)

cps_3 <- c(1,forecast_mosumvar$cps)#
equal_res_store <- matrix(0, nrow = sample_end-sample_start-5, ncol = 147)

for (t in (sample_start+1):(sample_end-6) ) { #to end of regime
  fm_3 <- get.factor.model(t(data_mosumvar[1:t,]), q = r_)
  factors_3 <- t(fm_3$f[1:r_,])
  full <- ar(factors_3, aic = F, order.max = 1) ## all
  window_start <- max(cps_3[t-18>cps_3]) ## at least 18 ahead
  window <- ar(factors_3[window_start:t,], aic = F, order.max = 1) ##this segment 
  window_all <- pooled_forecast(factors_3[,], cp = window_start-1, window_size = 30, p =1, n.ahead = 6, weights = "equal") ## pooled equal
  window_robust <- pooled_forecast(factors_3[,], cp = window_start-1, window_size = 30, p =1, n.ahead = 6, weights = "robust") ## pooled robust
  roll_window <- ar(factors_3[max(t-30,1):t,], aic = F, order.max = 1) ##this segment
  roll_window50 <- ar(factors_3[max(t-50,1):t,], aic = F, order.max = 1) ##this segment
  roll_window18 <- ar(factors_3[max(t-18,1):t,], aic = F, order.max = 1) ##this segment
  #vars::SVEC()
  # vecm_roll <- tsDyn::VECM(factors_3[t:max(t-30,1),], lag=1, include = "none") ##reversed data; VECM
  # predict.nlar(vecm_roll, newdata = as.matrix(pred_factor), n.ahead = 6, se.fit =F, nroll = 1)
  
  pred_data <- data_mosumvar[1+ (t):(t+5),] 
  pred_factor <- t(factors_3[ (t):(t),])
  
  h_step <- predict(full, pred_factor, n.ahead = 6, se.fit = F) 
  pred_full <-  h_step %*% t(fm_3$lam[,1:r_])
  h_step_window <- predict(window, pred_factor, n.ahead = 6, se.fit = F) 
  pred_window <-  h_step_window %*% t(fm_3$lam[,1:r_])
  pred_window_all <- window_all %*% t(fm_3$lam[,1:r_])
  pred_window_robust <- window_robust %*% t(fm_3$lam[,1:r_])
  h_step_roll_window <- predict(roll_window, pred_factor, n.ahead = 6, se.fit = F) 
  pred_roll_window <- h_step_roll_window %*% t(fm_3$lam[,1:r_])
  
  h_step_roll_window50 <- predict(roll_window50, pred_factor, n.ahead = 6, se.fit = F) 
  pred_roll_window50 <- h_step_roll_window50 %*% t(fm_3$lam[,1:r_])
  h_step_roll_window18 <- predict(roll_window18, pred_factor, n.ahead = 6, se.fit = F) 
  pred_roll_window18 <- h_step_roll_window18 %*% t(fm_3$lam[,1:r_])
  
  sum_pred_data <- rowSums(pred_data)#^2)#rep(1,6)#
  sum_pred_data_2 <- rowSums(pred_data^2)#
  for (h in 1:6) {
    q_forecast_out3_full[t-sample_start,h] <- sum( abs(pred_full[h,] - pred_data[h,] ) )/sum_pred_data[h] #^2
    q_square_full[t-sample_start,h] <- sum( abs(pred_full[h,] - pred_data[h,] )^2 )/sum_pred_data_2[h] #
    q_forecast_out3_window[t-sample_start,h] <- sum( abs(pred_window[h,] - pred_data[h,] )  )/sum_pred_data[h]
    q_square_window[t-sample_start,h] <- sum( abs(pred_window[h,] - pred_data[h,] )^2  )/sum_pred_data_2[h]
    q_forecast_out3_window_all[t-sample_start,h] <- sum( abs(pred_window_all[h,] - pred_data[h,] )  )/sum_pred_data[h]
    q_square_window_all[t-sample_start,h] <- sum( abs(pred_window_all[h,] - pred_data[h,] )^2  )/sum_pred_data_2[h]
    q_forecast_out3_window_robust[t-sample_start,h] <- sum( abs(pred_window_robust[h,] - pred_data[h,] )  )/sum_pred_data[h] 
    q_square_window_robust[t-sample_start,h] <- sum( abs(pred_window_robust[h,] - pred_data[h,] )^2  )/sum_pred_data_2[h] 
    q_forecast_out3_roll_window[t-sample_start,h] <- sum( abs(pred_roll_window[h,] - pred_data[h,] ) )/sum_pred_data[h]
    q_square_roll_window[t-sample_start,h] <- sum( abs(pred_roll_window[h,] - pred_data[h,] )^2 )/sum_pred_data_2[h]
    ## ensemble of full / roll
    ensemble_h <- (pred_full[h,] + pred_roll_window[h,])/2
    q_abs_ensemble[t-sample_start,h] <- sum( abs( ensemble_h - pred_data[h,] ) )/sum_pred_data[h] 
    q_square_ensemble[t-sample_start,h] <- sum( abs( ensemble_h - pred_data[h,] )^2 )/sum_pred_data_2[h]
    ## other roll window sizes
    q_abs_roll50[t-sample_start,h] <- sum( abs(pred_roll_window50[h,] - pred_data[h,] ) )/sum_pred_data[h]
    q_square_roll50[t-sample_start,h] <- sum( abs(pred_roll_window50[h,] - pred_data[h,] )^2 )/sum_pred_data_2[h]
    q_abs_roll18[t-sample_start,h] <- sum( abs(pred_roll_window18[h,] - pred_data[h,] ) )/sum_pred_data[h]
    q_square_roll18[t-sample_start,h] <- sum( abs(pred_roll_window18[h,] - pred_data[h,] )^2 )/sum_pred_data_2[h]
    
  }
 
  equal_res <- ( (pred_window_all[1,] - pred_data[1,] )^2 )/ pred_data[1,]^2
  equal_res[equal_res==Inf] <- 1
  equal_res_store[t-sample_start,] <- equal_res
  
  
  ##lm 
  mosumvar_gdp_window_3 <- window(mosumvar_gdp ,time(data_mosumvar)[1],time(data_mosumvar)[t] ) 
  mosumvar_gdp_window_3_validate <- window(mosumvar_gdp ,time(data_mosumvar)[t+1],time(data_mosumvar)[t+6] ) 
  gdp_lm_3 <- lm(mosumvar_gdp_window_3~ factors_3)
  
  bin_response_3 <- (mosumvar_gdp_window_3<0) ## is this recession?
 #bin_response_trim <- bin_response[!is.na(bin_response)] 
  gdp_glm_3 <- glm( bin_response_3 ~ factors_3[,], family = "binomial", weights = NULL)
  
  lm_pred <- cbind(1,h_step) %*%  (gdp_lm_3$coefficients)#predict.lm(gdp_lm, h_step) 
  glm_pred <- cbind(1,h_step) %*%  (gdp_lm_3$coefficients)#predict.glm(gdp_glm, h_step)
  
  lm_pred_1step[t-sample_start,] <- lm_pred 
  
  lm_pred_window <- cbind(1,h_step_window ) %*%  (gdp_lm_3$coefficients)#predict.lm(gdp_lm, h_step) 
  glm_pred_window <- cbind(1,h_step_window) %*%  (gdp_glm_3$coefficients)#predict.glm(gdp_glm, h_step)
  
  lm_pred_window_1step[t-sample_start,] <- lm_pred_window 
  
  lm_pred_window_all <- cbind(1,window_all) %*%  (gdp_lm_3$coefficients)
  glm_pred_window_all <- cbind(1,window_all) %*%  (gdp_glm_3$coefficients)#predict.glm(gdp_glm, h_step)
  
  lm_pred_window_all_1step[t-sample_start,] <- lm_pred_window_all
  
  lm_pred_window_robust <- cbind(1,window_robust) %*%  (gdp_lm_3$coefficients)
  #glm_pred_window_all <- cbind(1,window_all) %*%  (gdp_glm_3$coefficients)#predict.glm(gdp_glm, h_step)
  
  lm_pred_window_robust_1step[t-sample_start,] <- lm_pred_window_robust
  
  lm_pred_roll_window <- cbind(1,h_step_roll_window ) %*%  (gdp_lm_3$coefficients)
  glm_pred_roll_window <- cbind(1,h_step_roll_window ) %*%  (gdp_glm_3$coefficients)#predict.glm(gdp_glm, h_step)
  
  lm_pred_roll_window_1step[t-sample_start,] <- lm_pred_roll_window
  
  ## truth
  q_gdp_lm_out3[t-sample_start,] <- abs( c( gdp_lm_3$residuals[paste( t-sample_start)] , (lm_pred - mosumvar_gdp_window_3_validate )) )#^2
 # q_gdp_glm_out3[t-sample_start,] <- pracma::sigmoid(c(gdp_lm$fitted.values[paste(t-sample_start)],glm_pred )) 
  #q_gdp_glm_true3[t-sample_start,] <- mosumvar_gdp_window_3[ (t):(t+6)]>0
  
  ## store residuals
  q_gdp_lm_out3_window[t-sample_start,] <- abs( c( gdp_lm_3$residuals[paste( t-sample_start)] ,(lm_pred_window - mosumvar_gdp_window_3_validate ) ) )#^2
  #q_gdp_glm_out3_window[t-sample_start,] <- pracma::sigmoid(c(gdp_lm_3$fitted.values[paste(t-sample_start)],glm_pred_window))
  #abs(c( pracma::sigmoid(gdp_lm$fitted.values[paste(cp1+29+t)]) - (mosumvar_gdp_window[(cp1+29+t)]>class_thresh),
  #                                 pracma::sigmoid(glm_pred_window) - (mosumvar_gdp_window[cp1+30+(t):(t+5)]<class_thresh) ) ) 
  #abs(pracma::sigmoid(glm_pred_window) - (mosumvar_gdp_window[cp1+sample_start+(t):(t+5)]>0) )
  q_gdp_lm_out3_window_all[t-sample_start,] <- abs( c( gdp_lm_3$residuals[paste( t-sample_start)] ,(lm_pred_window_all - mosumvar_gdp_window_3_validate ) ) )#^2
  #q_gdp_glm_out3_window_all[t-sample_start,] <- pracma::sigmoid(c(gdp_lm_3$fitted.values[paste( t-sample_start)],glm_pred_window_all))
  
  q_gdp_lm_out3_window_robust[t-sample_start,] <- abs( c( gdp_lm_3$residuals[paste( t-sample_start)] ,(lm_pred_window_robust - mosumvar_gdp_window_3_validate ) ) )#^2
  
  q_gdp_lm_out3_roll_window[t-sample_start,] <- abs( c( gdp_lm_3$residuals[paste( t-sample_start)] ,(lm_pred_roll_window - mosumvar_gdp_window_3_validate ) ) ) #^2
  #q_gdp_glm_out3_roll_window[t-sample_start,] <- pracma::sigmoid(c(gdp_lm_3$fitted.values[paste( t-sample_start)],glm_pred_roll_window))
}

# q_F.5_3_full <- q_F.5_3_window <- q_F.5_3_window_all <- rep(0,7)
# for (h in 1:7) {
#   q_F.5_3_full[h] <- getF.5(q_gdp_glm_out3[,h], q_gdp_glm_true3[,h], class_thresh)
#   q_F.5_3_window[h] <- getF.5(q_gdp_glm_out3_window[,h], q_gdp_glm_true3[,h], class_thresh)
#   q_F.5_3_window_all[h] <- getF.5(q_gdp_glm_out3_window_all[,h], q_gdp_glm_true3[,h], class_thresh)
# }
# 
# 
# q_precision_3_full <- q_precision_3_window <- q_precision_3_window_all <- rep(0,7)
# for (h in 1:7) {
#   q_precision_3_full[h] <- getprecision(q_gdp_glm_out3[,h], q_gdp_glm_true3[,h], class_thresh)
#   q_precision_3_window[h] <- getprecision(q_gdp_glm_out3_window[,h], q_gdp_glm_true3[,h], class_thresh)
#   q_precision_3_window_all[h] <- getprecision(q_gdp_glm_out3_window_all[,h], q_gdp_glm_true3[,h], class_thresh)
# }

## plots #############################################
  plot.ts(mosumvar_gdp[!is.na(mosumvar_gdp)], ylab = "Smoothed GDP")
  ts.plot(cbind(q_forecast_out3_full[,1],q_forecast_out3_window[,1],q_forecast_out3_window_all[,1],q_forecast_out3_window_robust[,1], q_forecast_out3_roll_window[,1]), col = 1:4)
    # x cast
  par(mfrow = c(1,2))
    plot(colMeans(q_forecast_out3_full), type = "b", xlab = "Horizon", ylab = "X Absolute Error (normalised)", ylim = c(2.2,2.7), col = 1 );#c(0.995,1.001) 2.4,2.7
    lines(colMeans(q_forecast_out3_window) , type = "b", col = 2);
    lines(colMeans(q_forecast_out3_window_all,T) , type = "b", col = 3);
    lines(colMeans(q_forecast_out3_window_robust,T) , type = "b", col =4);
    lines(colMeans(q_forecast_out3_roll_window,T) , type = "b", col = 5);
    lines(colMeans(q_abs_ensemble,T) , type = "b", col = 6);
    legend(x="bottomright", legend = c("Full","Current","Equal","Robust","Rolling"), fill = 1:5, cex = 0.65)#fill = c("blue","lightblue","darkblue","grey","black")
    
    plot(colMeans(q_square_full), type = "b", xlab = "Horizon", ylab = "X Square Error (normalised)", ylim = c(0.998,1.00), col = 1 );#c(0.995,1.001) 2.4,2.7
    lines(colMeans(q_square_window) , type = "b", col = 2);
    lines(colMeans(q_square_window_all,T) , type = "b", col = 3);
    lines(colMeans(q_square_window_robust,T) , type = "b", col =4);
    lines(colMeans(q_square_roll_window,T) , type = "b", col = 5);
    lines(colMeans(q_square_ensemble,T) , type = "b", col = 6);
  par(mfrow = c(1,1))
  # roll window comparison
  par(mfrow = c(1,2))
    plot(colMeans(q_forecast_out3_roll_window), type = "b", xlab = "Horizon", ylab = "X Absolute Error (normalised)", ylim = c(2.5,2.8), col = 1 );#c(0.995,1.001) 2.4,2.7
    lines(colMeans(q_abs_roll50) , type = "b", col = 2);
    lines(colMeans(q_abs_roll18,T) , type = "b", col = 3);
    legend(x="bottomright", title = "Window size", legend = c("30","50","18"), fill = 1:3, cex = 0.65)#fill = c("blue","lightblue","darkblue","grey","black")
    
    plot(colMeans(q_square_roll_window), type = "b", xlab = "Horizon", ylab = "X Square Error (normalised)", ylim = c(0.998,0.999), col = 1 );#c(0.995,1.001) 2.4,2.7
    lines(colMeans(q_square_roll50) , type = "b", col = 2);
    lines(colMeans(q_square_roll18,T) , type = "b", col = 3); 
  par(mfrow = c(1,1))  
    
#  legend(x="bottomleft", legend = c("Full","Current","Equal","Robust","Rolling"), fill = 1:5 )
  # plot(apply(q_forecast_out3_full,2,sd, T), type = "b", xlab = "Horizon", ylab = "X Error Variance", ylim =  c(0, 1), col = 1 );# c(0,.005)
  # lines(apply(q_forecast_out3_window,2,sd,T) , type = "b", col = 2);
  # lines(apply(q_forecast_out3_window_all,2,sd,T) , type = "b", col = 3);
  # lines(apply(q_forecast_out3_window_robust,2,sd,T) , type = "b", col = 4);
  # lines(apply(q_forecast_out3_roll_window,2,sd,T) , type = "b", col = 5);
  #legend(x="bottomright", legend = c("Full","Current","Equal","Robust","Rolling"), fill = 1:5)
  # y cast
  plot(colMeans(q_gdp_lm_out3,T), type = "b", xaxt = "n", xlab = "Horizon", ylim = range(colMeans(q_gdp_lm_out3_window_all,T)), ylab = "Y Error", col = "red" ); axis(1, at = 1:7, labels = 0:6);
  lines(colMeans(q_gdp_lm_out3_window,T) , type = "b", col = "pink");lines(colMeans(q_gdp_lm_out3_window_all,T) , type = "b", col = "darkred");
  lines(colMeans(q_gdp_lm_out3_roll_window,T) , type = "b", col = "black");
  legend(x="bottomright", legend = c("Full","Current","Pooled","Rolling"), fill = c("red","pink","darkred","black"))
  
  plot(apply(q_gdp_lm_out3,2,sd,T), type = "b", xaxt = "n", xlab = "Horizon", ylim = c(-500,500), ylab = "Y Error", col = "red" ); axis(1, at = 1:7, labels = 0:6);
  lines(apply(q_gdp_lm_out3_window,2,sd,T) , type = "b", col = "pink");
  lines(apply(q_gdp_lm_out3_window_all,2,sd,T) , type = "b", col = "darkred");
  lines(apply(q_gdp_lm_out3_roll_window,2,sd,T) , type = "b", col = "black");
  legend(x="bottomright", legend = c("Full","Current","Pooled","Rolling"), fill = c("red","pink","darkred","black"))
  #F score
  plot( (q_F.5_3_full ), type = "b", xaxt = "n", xlab = "Horizon",ylim = c(0.0,1),   ylab = "F0.5 Score", col = "yellow" ); axis(1, at = 1:7, labels = 0:6)
  lines( (q_F.5_3_window) , type = "b", col = "green");
  lines( (q_F.5_3_window_all) , type = "b", col = "darkgreen");
  legend(x="bottomright", legend = c("Full","Current","Pooled"), fill = c("yellow","green","darkgreen"))
  # precision
  plot( (q_precision_3_full ), type = "b", xaxt = "n", xlab = "Horizon",ylim = c(0.0,1),   ylab = "Precision", col = "yellow" ); axis(1, at = 1:7, labels = 0:6)
  lines( (q_precision_3_window) , type = "b", col = "green");
  lines( (q_precision_3_window_all) , type = "b", col = "darkgreen");
  legend(x="bottomright", legend = c("Full","Current","Pooled"), fill = c("yellow","green","darkgreen"))
  # probabilistic classifier
  ts.plot(pnorm(gdp_lm_3$fitted.values, sd= mad(gdp_lm_3$residuals[])), ylab = "P(Positive)");abline(h=0.5,col="red") #nowcast prob under gaussian
  lines(1-bin_response_3[!is.na(bin_response_3)], lty = 3)
  ts.plot(pnorm(lm_pred_1step[,1], sd= sqrt(mean(q_gdp_lm_out3[,2], na.rm=T))), ylab = "P(Positive)");abline(h=0.5,col="red") #one step cast prob under gaussian
  ts.plot(pt(lm_pred_1step[,1]/sqrt(mean(q_gdp_lm_out3[,2], na.rm=T)), df=2), ylab = "P(Positive)");abline(h=0.5,col="red") #one step cast prob under t dist
  ts.plot(pnorm(lm_pred_1step[,2], sd= sqrt(mean(q_gdp_lm_out3[,3], na.rm=T))), ylab = "P(Positive)");abline(h=0.5,col="red") #one step cast prob under gaussian
  ts.plot(pnorm(lm_pred_1step[,6], sd= sqrt(mean(q_gdp_lm_out3[,7], na.rm=T))), ylab = "P(Positive)");abline(h=0.5,col="red") #one step cast prob under gaussian
  
  colMeans(cbind(q_gdp_lm_out3[,2],q_gdp_lm_out3_window[,2], q_gdp_lm_out3_window_all[,2], q_gdp_lm_out3_roll_window[,2]),T    )
  # nowcast with prediction intervals
  plot(gdp_lm_3$fitted.values, type = "l", ylim = range(predict(gdp_lm_3,mosumvar_gdp_window_3, interval = "prediction")), xlab = "Time",ylab = "Y" )
  lines(predict(gdp_lm_3,mosumvar_gdp_window_3, interval = "prediction")[!is.na(mosumvar_gdp_window_3),2] ,lty = "dotted")#gdp_lm_3$fitted.values + 1.96*summary(gdp_lm_3)$sigma/sqrt(length(gdp_lm_3$fitted.values )
  lines(predict(gdp_lm_3,mosumvar_gdp_window_3, interval = "prediction")[!is.na(mosumvar_gdp_window_3),3]  ,lty = "dotted")
  points(mosumvar_gdp_window_3[!is.na(mosumvar_gdp_window_3)])
  # factor coeffs - this is only useful in so far as conf int contains zero
  plot(gdp_lm_3$coefficients, axes = F, ylim= range(confint(gdp_lm_3)));axis(1, at=1:7, labels=names(gdp_lm_3$coefficients),srt=180);axis(2); abline(h=0)
  points(confint(gdp_lm_3)[,1], pch = 24, bg = "blue")
  points(confint(gdp_lm_3 )[,2], pch = 25, bg = "blue")
  plot.ts(factors_3, main = "") #intercept positive, f1 market, f2 volatility?
  # panel errors
  which(colMeans(equal_res_store,T)>1000, arr.ind = T)
  barplot(colMeans(equal_res_store,T))

## score ######
logscore <- function(r,true) sum(true*log(r) + (1-true)*log(1-r), na.rm = T)


## SCORE THESE OVER ALL HORIZONS
nowcast_probs <-pnorm(gdp_lm_3$fitted.values, sd= mad(gdp_lm$residuals[]))
logscore(nowcast_probs, (mosumvar_gdp_window_3>0)[!is.na(mosumvar_gdp_window_3)] ) 
#1step
nowcast_scores <- matrix(0, 6, 5)
for (h in 1:6) {
  step_truth <- (mosumvar_gdp_window_3>0)[(sample_start+1):(nrow(factors_)-5)+h]
  
  nowcast_probs_1step <-pnorm(lm_pred_1step[,h], sd= sqrt(mean(q_gdp_lm_out3[,h+1]^2, na.rm=T))) 
  nowcast_scores[h,1] <- logscore(nowcast_probs_1step,  na.omit(step_truth))
  
  nowcast_probs_window_1step <-pnorm(lm_pred_window_1step[,h], sd= sqrt(mean(q_gdp_lm_out3_window[,h+1]^2, na.rm=T)) )
  nowcast_scores[h,2] <-logscore(nowcast_probs_window_1step,  na.omit(step_truth))
  
  nowcast_probs_window_all_1step <-pnorm(lm_pred_window_all_1step[,h], sd= sqrt(mean(q_gdp_lm_out3_window_all[,h+1]^2, na.rm=T)))
  nowcast_scores[h,3] <-logscore(nowcast_probs_window_all_1step,  na.omit(step_truth))
  
  nowcast_probs_window_robust_1step <-pnorm(lm_pred_window_robust_1step[,h], sd= sqrt(mean(q_gdp_lm_out3_window_robust[,h+1]^2, na.rm=T)))
  nowcast_scores[h,4] <-logscore(nowcast_probs_window_robust_1step,  na.omit(step_truth))
  
  nowcast_probs_roll_window_1step <-pnorm(lm_pred_roll_window_1step[,h], sd= sqrt(mean(q_gdp_lm_out3_roll_window[,h+1]^2, na.rm=T)))
  nowcast_scores[h,5] <-logscore(nowcast_probs_roll_window_1step,  na.omit(step_truth))
}
plot(nowcast_scores[,1], type="b", xlab= "Horizon", ylab = "Log Score", col = 1, ylim =range(nowcast_scores) )#c(-104.5,-98)
lines( nowcast_scores[,2], type = "b", col = 2);
lines(nowcast_scores[,3] , type = "b", col = 3);
lines(nowcast_scores[,4] , type = "b", col = 4);
lines(nowcast_scores[,5] , type = "b", col = 5);
legend(x="right", legend = c("Full","Current","Equal","Robust","Rolling"), fill = 1:5)#c("yellow","green","darkgreen","black")

# factor plot
plot.ts( factors_3, main = "" )
plot.ts(apply(factors_3,2,cumsum), main = "" )
rowscale_loadings <-  t(fm_3$lam[,1:r_]/rowSums( abs(fm_3$lam[,1:r_]))   )
image(rowscale_loadings)
cbind(apply(rowscale_loadings, 2, which.max), colnames(data_mosumvar))


