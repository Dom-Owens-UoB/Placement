## FRED-MD
library("nowcasting")
# 1 none, 2 monthly change, 3 monthly 2nd deriv, 4 log, 5 monthly change log, 6  monthly 2nd deriv log, 
transmap <- function(x){
  out1 <- x[-1,-1]
  trans_vec <- x[1,-1]
  for (trans in 1:length(trans_vec)) {
    if( trans_vec[trans] == 2) out1[,trans] <- c(NA,diff(out1[,trans]))
    if( trans_vec[trans] == 3) out1[,trans] <- c(NA,NA,diff(out1[,trans], differences =  2))
    if( trans_vec[trans] == 4) out1[,trans] <- log(out1[,trans])
    if( trans_vec[trans] == 5) out1[,trans] <- c(NA,diff( log(out1[,trans]) ))
    if( trans_vec[trans] == 6) out1[,trans] <- c(NA,NA,diff( log(out1[,trans]), differences = 2 ))
    if( trans_vec[trans] == 7) out1[,trans] <- c(NA, NA,diff(out1[-1,trans]/(out1[-nrow(out1),trans]) ))
  }
  out <-  data.frame(date = x[-1,1] ,out1)
  return(out)
}
  
  #c(0, 2, 0, 2, 2, 1)

fredmd <- read.csv("current.csv")
head(fredmd)
fredmd_1 <- fredmd[-1,]
ts.plot(fredmd_1[,2:5])
fredmd_2 <- transmap(fredmd[-(2:500),])
ts.plot(fredmd_2[-1,-1])

fredmd_3 <- fredmd_2[-(1:17),]
fredmd_3[,1] <- as.Date(fredmd_3[,1], format = "%m/%d/%y")
duplicates <- which(colnames(fredmd_3) %in% c(colnames(data_alfred), colnames(survey_series)))
fredmd_4 <- fredmd_3[,-duplicates[-1]]


#
fredmd_ts <- ts(fredmd_4[,-1], start = c(2004,6), end = c(2021, 6), frequency = 12)
fredmd_trans <- c(rep(0, ncol(fredmd_4)-2), survey_trans)
data_fredmd <- Bpanel(base = fredmd_ts,
                      trans = rep(0, ncol(fredmd_4)-1),
                      aggregate = FALSE, na.prop = .75, NA.replace = T)
data_only_fredmd <- ts.intersect(data_fredmd, data_survey.GDPC1 = data_survey[,25])
#which(colnames(data_fredmd) == "ULCNFB")
freq_only_fredmd <- c( rep(12, ncol(fredmd_4)-1),4)

factormodel_only_fredmd <-  factorcpt::get.factor.model(x=t(data_fredmd[1:205,]))
factormodel_only_fredmd$q.hat


nowcast_only_fredmd <- nowcast(formula = data_survey.GDPC1 ~ ., data = data_only_fredmd , r = 6, p = 2, q = 4, 
                          method = '2s_agg', frequency = freq_only_fredmd ) 

nowcast.plot(nowcast_only_fredmd, type = "fcst")
