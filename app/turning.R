## turning point exercise 
library(nowcasting)

nowcast.plot(nowcastUSGDP, type = "fcst")
plot.ts(data$RGDPR)

qtr_gdp <- timeSeries::diff(month2qtr(USGDP$base[,"RGDPGR"])) #diff(as.vector(USGDP$base[,"RGDPGR"]), 3, pad=1) #month2qtr(USGDP$base[,"RGDPGR"])
 

drop <- which(colnames(USGDP$base) %in% c("AVGWKLYCL" ,"CPOUT" ,    "M3"))

out_full <- rep(NA, 100)
for (t in 1:100) {
  data_trim <-   Bpanel(base =  ts(USGDP$base[1:(nrow(data) - 101 + t),-drop], start =time(data)[1], end = time(data)[(nrow(data) - 101 + t)], frequency =12 ),
                   trans = USGDP$legend$Transformation[-drop],
                   aggregate = FALSE) # ts(data[1:(ncol(data) - 51 + t),], start =time(data)[1], end = time(data)[(ncol(data) - 51 + t)], frequency =12 )
  
    #  time(data_trim) <- time(data)[1:(ncol(data) - 51 + t)]
  nowcastUSGDP_trim <- nowcast(formula = RGDPGR ~ ., data = data_trim, r = 2, p = 2, q = 2, 
                          method = '2s_agg', frequency = frequency)
  # I want monthly not quarterly gdp
  pred <- nowcastUSGDP_trim$reg$coefficients %*% c(1,nowcastUSGDP_trim$factors$dynamic_factors[nrow(nowcastUSGDP_trim$factors$dynamic_factors) -12,])
 # true <- data_trim[nrow(data_trim)-12,"RGDPGR"]
  
  #col2 <- na.omit(nowcastUSGDP_trim$yfcst[,2])
  #col3 <- na.omit(nowcastUSGDP_trim$yfcst[,3])
  #change <-  diff(c(col2[length(col2)], col3 ))
  #accuracy <- sign(change) == sign( window(qtr_gdp, time(col3)[1], time(col3)[4])  )
                        #sign(qtr_gdp[(nrow(data) - 37 +1:12 + t)])  
  #predict.lm(nowcastUSGDP_trim$reg, USGDP$base[1:(nrow(data) - 37 + t),-c(drop, ncol(USGDP$base)) ])
  #nowcastUSGDP_trim$reg$coefficients %*% cbind(1, nowcastUSGDP_trim$factors$dynamic_factors[(nrow(data) - 37 + t) +1:12,])
  out_full[t] <- pred#sign(pred) == sign(true) #consider nnon-zero threshold on pred classifier: -.1, -1? 
}
scale(out_full)
true <- window(data_trim[,"RGDPGR"],  start = time(data)[(nrow(data) - 99)])#data_trim[nrow(data_trim)-12 - 100:1,"RGDPGR"]
plot( scale(out_full), type = "l"); #lines(as.vector(na.omit(sign(true))), col = "blue");
lines( as.vector(true), col = "blue") 



## test monotonicity
# data_cor <- cor(data , use = "complete.obs")
# image(data_cor)
# ord <- order(data_cor[,ncol(data_cor)], decreasing = T) #which are most highly correlated with response?
# #ord <- ord[-(drop)]
# #drop2 <- c("AVGWKLYCL", "CPOUT")
# ord_set <-  ord[1:50]#[-which(colnames(USGDP$base) %in% drop2)] #c(which(ord %in% 1:50))
set1 <- c(1:50, 193)
#set1 <- set1[!(set1 %in% drop)]
frequency1 <- c(rep(12,  49), 4) 
out_1 <- rep(NA, 100)
for (t in 1:100) {
  data_trim <-   Bpanel(base =  ts(USGDP$base[1:(nrow(data) - 101 + t),set1], start =time(data)[1], end = time(data)[(nrow(data) - 101 + t)], frequency =12 ),
                        trans = USGDP$legend$Transformation[set1],
                        aggregate = FALSE) # ts(data[1:(ncol(data) - 51 + t),], start =time(data)[1], end = time(data)[(ncol(data) - 51 + t)], frequency =12 )
  
  #  time(data_trim) <- time(data)[1:(ncol(data) - 51 + t)]
  nowcastUSGDP_trim <- nowcast(formula = RGDPGR ~ ., data = data_trim, r = 2, p = 2, q = 2, 
                               method = '2s_agg', frequency = frequency1 )
  # I want monthly not quarterly gdp
  pred <- nowcastUSGDP_trim$reg$coefficients %*% c(1,nowcastUSGDP_trim$factors$dynamic_factors[nrow(nowcastUSGDP_trim$factors$dynamic_factors) -12,])
  # true <- data_trim[nrow(data_trim)-12,"RGDPGR"]
  
  #col2 <- na.omit(nowcastUSGDP_trim$yfcst[,2])
  #col3 <- na.omit(nowcastUSGDP_trim$yfcst[,3])
  #change <-  diff(c(col2[length(col2)], col3 ))
  #accuracy <- sign(change) == sign( window(qtr_gdp, time(col3)[1], time(col3)[4])  )
  #sign(qtr_gdp[(nrow(data) - 37 +1:12 + t)])  
  #predict.lm(nowcastUSGDP_trim$reg, USGDP$base[1:(nrow(data) - 37 + t),-c(drop, ncol(USGDP$base)) ])
  #nowcastUSGDP_trim$reg$coefficients %*% cbind(1, nowcastUSGDP_trim$factors$dynamic_factors[(nrow(data) - 37 + t) +1:12,])
  out_1[t] <- pred#sign(pred) == sign(true) #consider nnon-zero threshold on pred classifier: -.1, -1? 
}
scale(out_1)
 plot( scale(out_1), type = "l"); #lines(as.vector(na.omit(sign(true))), col = "blue");
lines( as.vector(true), col = "blue") 
 accuracy_out(out_1)
 
 
 
 
 
 
 
 
 set2 <- c(1:100, 193)
 #set1 <- set1[!(set1 %in% drop)]
 frequency2 <- c(rep(12,  99), 4) 
 out_2 <- rep(NA, 100)
 for (t in 1:100) {
   data_trim <-   Bpanel(base =  ts(USGDP$base[1:(nrow(data) - 101 + t),set2], start =time(data)[1], end = time(data)[(nrow(data) - 101 + t)], frequency =12 ),
                         trans = USGDP$legend$Transformation[set2],
                         aggregate = FALSE) # ts(data[1:(ncol(data) - 51 + t),], start =time(data)[1], end = time(data)[(ncol(data) - 51 + t)], frequency =12 )
   
   #  time(data_trim) <- time(data)[1:(ncol(data) - 51 + t)]
   nowcastUSGDP_trim <- nowcast(formula = RGDPGR ~ ., data = data_trim, r = 2, p = 2, q = 2, 
                                method = '2s_agg', frequency = frequency2 )
   # I want monthly not quarterly gdp
   pred <- nowcastUSGDP_trim$reg$coefficients %*% c(1,nowcastUSGDP_trim$factors$dynamic_factors[nrow(nowcastUSGDP_trim$factors$dynamic_factors) -12,])
   # true <- data_trim[nrow(data_trim)-12,"RGDPGR"]
   
   #col2 <- na.omit(nowcastUSGDP_trim$yfcst[,2])
   #col3 <- na.omit(nowcastUSGDP_trim$yfcst[,3])
   #change <-  diff(c(col2[length(col2)], col3 ))
   #accuracy <- sign(change) == sign( window(qtr_gdp, time(col3)[1], time(col3)[4])  )
   #sign(qtr_gdp[(nrow(data) - 37 +1:12 + t)])  
   #predict.lm(nowcastUSGDP_trim$reg, USGDP$base[1:(nrow(data) - 37 + t),-c(drop, ncol(USGDP$base)) ])
   #nowcastUSGDP_trim$reg$coefficients %*% cbind(1, nowcastUSGDP_trim$factors$dynamic_factors[(nrow(data) - 37 + t) +1:12,])
   out_2[t] <- pred#sign(pred) == sign(true) #consider nnon-zero threshold on pred classifier: -.1, -1? 
 }
 scale(out_2)
 plot( scale(out_2), type = "l"); #lines(as.vector(na.omit(sign(true))), col = "blue");
 lines( as.vector(true), col = "blue") 
 accuracy_out(out_2)
 
 
 
 
 set3 <- c(1:150, 193)
 #set1 <- set1[!(set1 %in% drop)]
 frequency3 <- c(rep(12,  147), 4) 
 out_3 <- rep(NA, 100)
 for (t in 1:100) {
   data_trim <-   Bpanel(base =  ts(USGDP$base[1:(nrow(data) - 101 + t),set3], start =time(data)[1], end = time(data)[(nrow(data) - 101 + t)], frequency =12 ),
                         trans = USGDP$legend$Transformation[set3],
                         aggregate = FALSE) # ts(data[1:(ncol(data) - 51 + t),], start =time(data)[1], end = time(data)[(ncol(data) - 51 + t)], frequency =12 )
   
   #  time(data_trim) <- time(data)[1:(ncol(data) - 51 + t)]
   nowcastUSGDP_trim <- nowcast(formula = RGDPGR ~ ., data = data_trim, r = 2, p = 2, q = 2, 
                                method = '2s_agg', frequency = frequency3 )
   # I want monthly not quarterly gdp
   pred <- nowcastUSGDP_trim$reg$coefficients %*% c(1,nowcastUSGDP_trim$factors$dynamic_factors[nrow(nowcastUSGDP_trim$factors$dynamic_factors) -12,])
   # true <- data_trim[nrow(data_trim)-12,"RGDPGR"]
   
   #col2 <- na.omit(nowcastUSGDP_trim$yfcst[,2])
   #col3 <- na.omit(nowcastUSGDP_trim$yfcst[,3])
   #change <-  diff(c(col2[length(col2)], col3 ))
   #accuracy <- sign(change) == sign( window(qtr_gdp, time(col3)[1], time(col3)[4])  )
   #sign(qtr_gdp[(nrow(data) - 37 +1:12 + t)])  
   #predict.lm(nowcastUSGDP_trim$reg, USGDP$base[1:(nrow(data) - 37 + t),-c(drop, ncol(USGDP$base)) ])
   #nowcastUSGDP_trim$reg$coefficients %*% cbind(1, nowcastUSGDP_trim$factors$dynamic_factors[(nrow(data) - 37 + t) +1:12,])
   out_3[t] <- pred#sign(pred) == sign(true) #consider nnon-zero threshold on pred classifier: -.1, -1? 
 }
 scale(out_3)
 plot(  (out_3), type = "l"); #lines(as.vector(na.omit(sign(true))), col = "blue");
 lines( as.vector(true), col = "blue") 
 accuracy_out(out_3)
 
 
 
 
 
 
 
 
 ##result 
 
 # accuracy_out <- function(out){
 #   out <- scale(out)
 #   ret <- 0 * (1:30)
 #   thresh <- -1* seq(0, 3, length.out = 30)
 #   for (i in 1:30) {
 #     ret[i] <- mean( na.omit((out < thresh[i]) == (true < 0)))
 #   }
 #   return(ret)
 # }
 # 