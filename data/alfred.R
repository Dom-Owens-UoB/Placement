## get data from alfred
library("alfred")
library("dplyr")


try_series <- alfred::get_fred_series("PAYEMS", series_name = "PAYEMS", observation_start = "2002-01-01", observation_end = "2021-06-01")

for (name in colnames(NYFED$base)[-1]) {
  new_series <- alfred::get_fred_series(name, series_name = name, observation_start = "2002-01-01", observation_end = "2021-06-01")
  try_series <- full_join(try_series, new_series, by = "date")
}

#try_series[,26] <- c(rep(NA,12), (try_series[-(1:12),26] - try_series[-(223:234),26])/try_series[-(223:234),26])
  
nyfed_trans <- NYFED$legend$Transformation
nyfed_trans[25] <- 6
alfred_ts <- ts(try_series[,-1], start = c(2002,1), end = c(2021, 6), frequency = 12)
data_alfred <- Bpanel(base = alfred_ts,
               trans = nyfed_trans,#NYFED$legend$Transformation,
               aggregate = FALSE, na.prop = .75)
#data_alfred[,25] <- c(NA,NA,NA, timeSeries::diff(data_alfred[,25],3))
freq_alfred <- rep(12,25); freq_alfred[c(20,25)] <- 4

ts.plot(data_alfred)

nowcast_alfred <- nowcast(formula = GDPC1 ~ ., data = data_alfred, r = 2, p = 2, q = 2, 
                        method = '2s_agg', frequency = freq_alfred)

nowcast.plot(nowcast_alfred, type = "fcst")
nowcast.plot(nowcast_alfred, type = "factors") ## this would indicate structural breaks




###
### add more series ------------------
###
FRED_vars <- data.frame( colnames(NYFED$base), NYFED$legend$Transformation, NYFED$legend$SeriesName)
FRED_vars

## which surveys to add?
nowcasting::USGDP$legend$mnemonic[181:189]
nowcasting::USGDP$legend$Source.for.data[181:189]
nowcasting::USGDP$legend$Transformation[181:189]
# add texas and chicago diffusion indices to encode surveys

surveys <- c("CFNAIDIFF","BACTSAMFRBDAL")

survey_series <- full_join(alfred::get_fred_series("CFNAIDIFF", series_name = "CFNAIDIFF", observation_start = "2004-06-01", observation_end = "2021-06-01"),
                           alfred::get_fred_series("BACTSAMFRBDAL", series_name = "BACTSAMFRBDAL", observation_start = "2004-06-01", observation_end = "2021-06-01"))
#try_series
# for (name in surveys ) {
#   new_series <- alfred::get_fred_series(name, series_name = name, observation_start = "2004-06-01", observation_end = "2021-06-01")
#   survey_series <- full_join(survey_series, new_series, by = "date")
# }

survey_ts <- ts(survey_series[,-1], start = c(2004,6), end = c(2021, 6), frequency = 12) #[- (1:29) ,]
survey_trans <- c(nyfed_trans, 0,0) #NYFED$legend$Transformation
data_survey <-  ts.intersect(data_alfred, survey_ts)
          # Bpanel(base = survey_ts,
          #             trans = survey_trans,
          #             aggregate = FALSE, na.prop = .75)
freq_survey <- rep(12,27); freq_survey[c(20,25)] <- 4

ts.plot(data_survey)

nowcast_survey <- nowcast(formula = data_alfred.GDPC1 ~ ., data = data_survey, r = 2, p = 2, q = 2, 
                          method = '2s_agg', frequency = freq_survey)

nowcast.plot(nowcast_survey, type = "fcst")
nowcast.plot(nowcast_survey, type = "factors") ## this would indicate structural breaks



###
### add FREDMD --------------
###

#FREDMD_names <- colnames(fredmd_trans[,-1]) 

#fredmd_series <- inner_join(fredmd_4, survey_series, by = "date")
fredmd_ts <- ts(fredmd_4[,-1], start = c(2004,6), end = c(2021, 6), frequency = 12)
fredmd_trans <- c(rep(0, ncol(fredmd_4)-2), survey_trans)
data_fredmd <- Bpanel(base = fredmd_ts,
                      trans = rep(0, ncol(fredmd_4)-1),
                      aggregate = FALSE, na.prop = .75, NA.replace = T)
data_big <- ts.intersect(data_fredmd, data_survey)
#which(colnames(data_fredmd) == "ULCNFB")
freq_fredmd <- c( rep(12, ncol(fredmd_4)-1),freq_survey)

nowcast_fredmd <- nowcast(formula = data_survey.data_alfred.GDPC1 ~ ., data = data_big , r = 6, p = 2, q = 6, 
                          method = '2s_agg', frequency = freq_fredmd ) 

nowcast.plot(nowcast_fredmd, type = "fcst")
nowcast.plot(nowcast_fredmd, type = "factors")
