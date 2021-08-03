library("factorcpt"); library("nowcasting")
run_factorcpt <- factor.seg.alg( t(na.omit(data_big[,])) , r = 6, rule = NULL, dw = NULL, sig.lev = 0.05)
  #factor.seg.alg( t(na.omit(data[,-190])), r = 2, q.seq = 2, rule = 1, dw = 20) # -25
run_factorcpt$cs.list
run_factorcpt$common.est.cps
run_factorcpt$idio.est.cps
time(data_big)[run_factorcpt$idio.est.cps]  #2012.833 i.e. 11/2012


ts.plot(data_big); abline(v = time(survey_ts)[run_factorcpt$common.est.cps], col = "red")
## fit model
cp_ts <- ts(survey_series[,-1], start = c(2012,8), end = c(2021, 6), frequency = 12) #[- (1:29) ,]
#survey_trans <- c(NYFED$legend$Transformation, 0,0)
data_cp <- Bpanel(base = cp_ts,
                      trans = survey_trans,
                      aggregate = FALSE, na.prop = .75)
#freq_survey <- rep(12,27); freq_survey[c(20,25)] <- 4
ts.plot(data_cp)

nowcast_cp <- nowcast(formula = GDPC1 ~ ., data = data_cp, r = run_factorcpt$r, p = , q = run_factorcpt$gfm$q.hat, 
                          method = '2s_agg', frequency = freq_survey)

nowcast.plot(nowcast_cp, type = "fcst")
nowcast.plot(nowcast_cp, type = "factors")  



###
### forecast comparison ----
###
cp_test <- ts(data_survey[,"GDPC1"] , start = c(2021,1),  end = c(2021, 6), frequency =12 )
## estimate with full sample
out_nyfed <- rep(NA, 6)
for (t in 1:6) {
  nyfed_trim <-   Bpanel(base =  ts(survey_series[,-1], start = c(2004,6),  end = c(2021, 6), frequency =12 ),
                        trans = survey_trans, na.prop = .8,
                        aggregate = FALSE, 
                        h = 12) # ts(data[1:(ncol(data) - 51 + t),], start =time(data)[1], end = time(data)[(ncol(data) - 51 + t)], frequency =12 )
  
  #  time(data_trim) <- time(data)[1:(ncol(data) - 51 + t)]
  nowcast_nyfed_trim <- nowcast(formula = GDPC1 ~ ., data = nyfed_trim, r =run_factorcpt$r, p = 2, q = run_factorcpt$gfm$q.hat, 
                               method = '2s_agg', frequency = freq_survey )
  # I want monthly not quarterly gdp
  pred <- nowcast_nyfed_trim$reg$coefficients %*% c(1, nowcast_nyfed_trim$factors$dynamic_factors[nrow(nowcast_nyfed_trim$factors$dynamic_factors) - 12,])
  # true <- data_trim[nrow(data_trim)-12,"RGDPGR"]
  
  #col2 <- na.omit(nowcastUSGDP_trim$yfcst[,2])
  #col3 <- na.omit(nowcastUSGDP_trim$yfcst[,3])
  #change <-  diff(c(col2[length(col2)], col3 ))
  #accuracy <- sign(change) == sign( window(qtr_gdp, time(col3)[1], time(col3)[4])  )
  #sign(qtr_gdp[(nrow(data) - 37 +1:12 + t)])  
  #predict.lm(nowcastUSGDP_trim$reg, USGDP$base[1:(nrow(data) - 37 + t),-c(drop, ncol(USGDP$base)) ])
  #nowcastUSGDP_trim$reg$coefficients %*% cbind(1, nowcastUSGDP_trim$factors$dynamic_factors[(nrow(data) - 37 + t) +1:12,])
  out_nyfed[t] <- pred#sign(pred) == sign(true) #consider nnon-zero threshold on pred classifier: -.1, -1? 
  if(t==6) print( summary(nowcast_nyfed_trim$reg))
}
sum( (cp_test -out_nyfed)[c(3,6)]^2 ) 

## estimate with full sample
out_nyfed_cp <- rep(NA, 6)
for (t in 1:6) {
  nyfed_trim <-   Bpanel(base =  ts(survey_series[,-1], start = c(2012,8),  end = c(2021, t), frequency =12 ),
                         trans = survey_trans, na.prop = .8,
                         aggregate = FALSE, 
                         h = 12) # ts(data[1:(ncol(data) - 51 + t),], start =time(data)[1], end = time(data)[(ncol(data) - 51 + t)], frequency =12 )
  
  #  time(data_trim) <- time(data)[1:(ncol(data) - 51 + t)]
  nowcast_nyfed_trim <- nowcast(formula = GDPC1 ~ ., data = nyfed_trim,  r =run_factorcpt$r, p = 2, q = run_factorcpt$gfm$q.hat, 
                                method = '2s_agg', frequency = freq_survey )
  # I want monthly not quarterly gdp
  pred <- nowcast_nyfed_trim$reg$coefficients %*% c(1, nowcast_nyfed_trim$factors$dynamic_factors[nrow(nowcast_nyfed_trim$factors$dynamic_factors) - 12,])
  # true <- data_trim[nrow(data_trim)-12,"RGDPGR"]
  
  #col2 <- na.omit(nowcastUSGDP_trim$yfcst[,2])
  #col3 <- na.omit(nowcastUSGDP_trim$yfcst[,3])
  #change <-  diff(c(col2[length(col2)], col3 ))
  #accuracy <- sign(change) == sign( window(qtr_gdp, time(col3)[1], time(col3)[4])  )
  #sign(qtr_gdp[(nrow(data) - 37 +1:12 + t)])  
  #predict.lm(nowcastUSGDP_trim$reg, USGDP$base[1:(nrow(data) - 37 + t),-c(drop, ncol(USGDP$base)) ])
  #nowcastUSGDP_trim$reg$coefficients %*% cbind(1, nowcastUSGDP_trim$factors$dynamic_factors[(nrow(data) - 37 + t) +1:12,])
  out_nyfed_cp[t] <- pred#sign(pred) == sign(true) #consider nnon-zero threshold on pred classifier: -.1, -1? 
  if(t==6) print( summary(nowcast_nyfed_trim$reg))
}
sum( (cp_test -out_nyfed_cp)[c(3,6)]^2 ) ## order of magnitude improvement

nowcast.plot(nowcast_nyfed_trim, "factors")

# more relevant comparison is to (1) rolling windows (2) turning points 
