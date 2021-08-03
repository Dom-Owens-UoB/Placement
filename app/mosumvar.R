## try change point analysis using mosumvar on estimated factor series
library("nowcasting"); library("mosumvar")
# mosumvar_ts <- ts(survey_series[,-1], start = c(2004,6), end = c(2021, 6), frequency = 12) #[- (1:29) ,]
# #survey_trans <- c(NYFED$legend$Transformation, 0,0)
# data_mosumvar <- Bpanel(base = mosumvar_ts,
#                   trans = survey_trans,
#                   aggregate = FALSE, na.prop = .75)
mosumvar_ts <- ts(data_big[,-146], start = time(data_big)[1], end = time(data_big)[nrow(data_big)], frequency = 12) #[- (1:29) ,]
mosumvar_gdp <- ts(data_big[,146], start = time(data_big)[1], end = time(data_big)[nrow(data_big)], frequency = 4)

#-146
# data_fredmd_trim <- ts(data_big
data_mosumvar <- Bpanel(base = mosumvar_ts,
                           trans = rep(0,ncol(mosumvar_ts)),
                           aggregate = FALSE, na.prop = 1, h = 12) 

#data_alfred[,25]
mosumvar_gdp <- diff(alfred_ts[,25],3)  


nowcast_mosumvar <- nowcast(formula = data_survey.data_alfred.GDPC1 ~ ., data = data_mosumvar, r = 6, p = 1, q = 4, 
                      method = '2s_agg', frequency = freq_fredmd)


## ar model
fcpt_mosumvar <- factorcpt::get.factor.model(t(data_mosumvar[1:200,]))
fcpt_ar <- ar(t(fcpt_mosumvar$f[1:6,1:53]), aic = T, order.max = 1)
fcpt_ar2 <- ar(t(fcpt_mosumvar$f[1:6,54:104]), aic = T, order.max = 1)
fcpt_ar3 <- ar(t(fcpt_mosumvar$f[1:6,105:157]), aic = T, order.max = 1)
#fcpt_mosumvar$lam[,1:6]

cp_mosumvar <- mosumvar::mosumvar( t(fcpt_mosumvar$f[1:6,]), p = 1, G= 24, nu = .1, method = "Score")
              #mosumvar::mosumvar(nowcast_mosumvar$factors$dynamic_factors, p = 2, G= 18, nu = .1)
cp_mosumvar_univ <- mosumvar::mosum_univ( t(fcpt_mosumvar$f[1:6,]), p = 1, G= 18, nu = .1, method = "Wald", do_bootstrap = T)


mosumvar_factor(data_mosumvar[1:200,], p=1, G=24, nu=.1, method = "Score")
mosumvar_factor(data_mosumvar[1:200,], p=1, G=18, nu=.01, method = "Score", univ = T)
###
### test for loadings of a single variable ----------- 
###
#data_reorder <- data_mosumvar[, c(25, (1:27)[-25])]
lm_data <- cbind(data_mosumvar[3:217,"GDPC1"], nowcast_mosumvar$factors$dynamic_factors[3:217,], nowcast_mosumvar$factors$dynamic_factors[2:216,],nowcast_mosumvar$factors$dynamic_factors[1:215,])
cp_lm_mosumvar <- mosumvar::mosum_lm(lm_data, G = 20, method = "Wald")
time(data_mosumvar)[cp_lm_mosumvar$cps + 2]
