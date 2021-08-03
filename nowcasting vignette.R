library("nowcasting")
data <- Bpanel(base = USGDP$base,
               trans = USGDP$legend$Transformation,
               aggregate = FALSE)
frequency <- c(rep(12, ncol(data) -1), 4) #all monthly but gdp
nowcastUSGDP <- nowcast(formula = RGDPGR ~ ., data = data, r = 2, p = 2, q = 2, 
                        method = '2s_agg', frequency = frequency)
res <- ts(nowcastUSGDP$reg$residuals, start = start(data), frequency = 4)
acf(window(res, start = c(1985,1), end = c(2004,4)))

# y fcst
nowcast.plot(nowcastUSGDP, type = "fcst")

# factors
nowcast.plot(nowcastUSGDP, type = "factors") 

# how much of the variability in the dataset is explained by each factor 
nowcast.plot(nowcastUSGDP, type = "eigenvalues")

# importance of each variable in the first factor
nowcast.plot(nowcastUSGDP, type = "eigenvectors") 
