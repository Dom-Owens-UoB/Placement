## turning points as classifying sign of rate of change
library("nowcasting")

#true_roc <- window(survey_ts[,"GDPC1"],  start = time(survey_ts)[(nrow(survey_ts) - 100)])
#true_roc_diff <- diff(true_roc, lag = 3)/
true_roc <- data_survey[(nrow(survey_ts) - 99 - 5):(nrow(survey_ts)- 5),"data_alfred.GDPC1"]#diff(data_survey[(nrow(survey_ts) - 103):nrow(survey_ts),"GDPC1"], 3)/data_survey[(nrow(survey_ts) - 103):(nrow(survey_ts)-3),"GDPC1"] ##quarterly ROC
true_roc_diff <- diff(true_roc, lag = 3) #is my rate of change going up or down?
plot.ts(true_roc) #classify on this

#true_roc_diff_na <- true_roc_diff[!is.na(true_roc_diff)]
survey_trans_yr <- survey_trans
#survey_trans_yr[which(survey_trans_yr)==7] <- 6

out_survey <- rep(NA, 100)
for (t in 1:100) {
  survey_ts_trim <- ts(data_survey[0:(nrow(survey_ts) - 100) + t,], start = time(survey_ts)[t], end = time(survey_ts)[(nrow(survey_ts) - 100 + t)], frequency = 12) #[- (1:29) ,]
  data_survey_trim <- Bpanel(base = survey_ts_trim,
                        trans = (survey_trans_yr),
                        aggregate = FALSE, na.prop = .75) 
  #  time(data_trim) <- time(data)[1:(ncol(data) - 51 + t)]
  nowcast_survey_trim <- nowcast(formula = data_alfred.GDPC1 ~ ., data = data_survey_trim, r = 2, p = 2, q = 2, 
                               method = '2s_agg', frequency = freq_survey)
  # I want monthly not quarterly gdp
  pred <- nowcast_survey_trim$reg$coefficients %*% c(1,nowcast_survey_trim$factors$dynamic_factors[nrow(nowcast_survey_trim$factors$dynamic_factors) -12,])
  # true <- data_trim[nrow(data_trim)-12,"RGDPGR"]
  
   out_survey[t] <- pred#sign(pred) == sign(true) #consider nnon-zero threshold on pred classifier: -.1, -1? 
}
ts.plot(out_survey)
ts.plot(true_roc_diff)
roc_pred <-  diff(out_survey, lag = 3)
ts.plot(roc_pred)

ts.plot( cbind( scale(roc_pred), scale(true_roc_diff)) )
ts.plot( sign(cbind( roc_pred, true_roc_diff) ))
#plot(sign(roc_pred), sign(true_roc_diff) )

mean(sign(roc_pred[true_roc_diff!= 0 ])  == sign(true_roc_diff[true_roc_diff!= 0 ]) )

###
### For FREDMD DATA ----
###

out_fredmd <- rep(NA, 100)
for (t in 1:100) {
  fredmd_ts_trim <- ts(data_big[0:(nrow(data_big) - 100) + t,], start = time(data_big)[t], end = time(data_big)[(nrow(data_big) - 100   + t)], frequency = 12) #[- (1:29) ,]
  
 # data_fredmd_trim <- ts(data_big
  data_fredmd_trim <- Bpanel(base = fredmd_ts_trim,
                             trans = rep(0,ncol(fredmd_ts_trim)),
                             aggregate = FALSE, na.prop = .75, h = 12) 
  
  #  time(data_trim) <- time(data)[1:(ncol(data) - 51 + t)]
  nowcast_fredmd_trim <- nowcast(formula = data_survey.data_alfred.GDPC1 ~ ., data = data_fredmd_trim , r = 6, p = 3, q = 6, 
                            method = '2s_agg', frequency = freq_fredmd)#[-1] ) 
  # I want monthly not quarterly gdp
  pred <- nowcast_fredmd_trim$reg$coefficients %*% c(1,nowcast_fredmd_trim$factors$dynamic_factors[nrow(nowcast_fredmd_trim$factors$dynamic_factors) -12,])
  # true <- data_trim[nrow(data_trim)-12,"RGDPGR"]
  
  out_fredmd[t] <- pred#sign(pred) == sign(true) #consider nnon-zero threshold on pred classifier: -.1, -1? 
}
nowcast.plot(nowcast_fredmd_trim)
ts.plot(out_fredmd)

roc_pred_fredmd <-  diff(out_fredmd, lag = 3)
ts.plot(roc_pred_fredmd)

ts.plot( cbind( scale(roc_pred_fredmd), scale(true_roc_diff)) )
ts.plot( sign(cbind( roc_pred_fredmd, true_roc_diff) ))
#plot(sign(roc_pred), sign(true_roc_diff) )

mean(sign(roc_pred_fredmd[   ])  == sign(true_roc_diff[ ]) )

## start at estimated change point
out_fredmd_cp <- rep(NA, 50)
for (t in 1:50) {
  fredmd_ts_trim <- ts(data_big[run_factorcpt$idio.est.cps[2]:(nrow(data_big) - 50+ t) ,], start = time(data_big)[run_factorcpt$idio.est.cps[2]], end = time(data_big)[(nrow(data_big) - 50  + t)], frequency = 12) #[- (1:29) ,]
  
  # data_fredmd_trim <- ts(data_big
  data_fredmd_trim <- Bpanel(base = fredmd_ts_trim,
                             trans = rep(0,ncol(fredmd_ts_trim)),
                             aggregate = FALSE, na.prop = 1, h = 12) 
  
  #  time(data_trim) <- time(data)[1:(ncol(data) - 51 + t)]
  nowcast_fredmd_trim_cp <- nowcast(formula = data_survey.data_alfred.GDPC1 ~ ., data = data_fredmd_trim , r = 6, p = 2, q = 6, 
                                 method = '2s_agg', frequency = freq_fredmd)#[-1] ) 
  # I want monthly not quarterly gdp
  pred <- nowcast_fredmd_trim_cp$reg$coefficients %*% c(1,nowcast_fredmd_trim_cp$factors$dynamic_factors[nrow(nowcast_fredmd_trim_cp$factors$dynamic_factors) -12,])
  # true <- data_trim[nrow(data_trim)-12,"RGDPGR"]
  
  out_fredmd_cp[t] <- pred#sign(pred) == sign(true) #consider nnon-zero threshold on pred classifier: -.1, -1? 
}
nowcast.plot(nowcast_fredmd_trim_cp)
ts.plot(out_fredmd_cp)

roc_pred_fredmd_cp <-  diff(out_fredmd_cp, lag = 3)
ts.plot(roc_pred_fredmd_cp)

ts.plot( cbind( scale(roc_pred_fredmd_cp), scale(true_roc_diff[54:100])) )
ts.plot( sign(cbind( roc_pred_fredmd_cp, true_roc_diff[54:100]) ))
#plot(sign(roc_pred), sign(true_roc_diff) )

addrange <- seq(from= min(roc_pred_fredmd_cp), to = max(roc_pred_fredmd_cp), length.out = 20)
addclass <- addrange*0
for (ii in 1:length(addrange) ) {
  addclass[ii] <-   mean(sign(roc_pred_fredmd_cp[   ] + addrange[ii])  == sign(true_roc_diff[51:97 ]) )
          ##slightly higher...
} 
max(addclass)##0.638


## try bin reg for scoring
bin_response <- sign(true_roc_diff[51:97 ])
bin_response[bin_response == -1] <- 0
bin_try <- glm( bin_response~ roc_pred_fredmd_cp, family = "binomial")
sum(log(bin_try$fitted.values[(bin_try$fitted.values > .5) == bin_response ])) 
bin_score <- 1:100 * 0
for (ii in 1:100 ) {
  bin_score[ii] <- sum(log(bin_try$fitted.values[(bin_try$fitted.values > ii* 0.01) == bin_response ])) #log score rule
}
which.max(bin_score) * 0.01

accuracy_ <- sign(roc_pred_fredmd_cp[   ] + addrange[which.max(addclass)])  == sign(true_roc_diff[51:97 ]) 
hist(roc_pred_fredmd_cp[accuracy_])
hist(roc_pred_fredmd_cp[!accuracy_]) ## missclassified observations are sometimes of large order

mad(nowcast_fredmd_trim_cp$reg$residuals) ##classifying with se is similar to data-driven threshold
