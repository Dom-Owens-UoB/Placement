## mosumvar::mosumvar_factor
# mosumvar extension: MOSUM change point analysis for a factor model with low-dimensional VAR dynamics 


mosumvar_factor <- function(x, p, G, method = c("Wald", "Score")[1], estim = c("DiagC","DiagH")[1], varEstim = c("Local", "Global")[1], alpha = 0.05, 
                            criterion = c("eps", "eta")[1], nu = 0.25, do_bootstrap = FALSE, M = 1000, 
                            univ = FALSE, rm_cross_terms = TRUE, global_resids = TRUE,
                            r = NULL, r_max = NA) {
  fcpt_mosumvar <- factorcpt::get.factor.model(t(x))
  #fcpt_ar <- ar(t(fcpt_mosumvar$f[1:6,1:50]), aic = T, order.max = 4)
  if(is.null(r)) r <- max(fcpt_mosumvar$q.hat, 1)
  r <- min(r, r_max, na.rm = T)
  if(!univ) {cp_mosumvar <- mosumvar::mosumvar( x=t(fcpt_mosumvar$f[1:r,]), p=p, G=G, method = method, estim =estim, 
                                     varEstim = varEstim, alpha = alpha, 
                                     criterion = criterion, nu = nu, do_bootstrap = do_bootstrap, M = M)
  #mosumvar::mosumvar(nowcast_mosumvar$factors$dynamic_factors, p = 2, G= 18, nu = .1)
  } else {cp_mosumvar <- mosumvar::mosum_univ( x=t(fcpt_mosumvar$f[1:r,]), p=p, G=G, method = method, estim =estim, 
                                            varEstim = varEstim, alpha = alpha, rm_cross_terms = rm_cross_terms, global_resids = global_resids,
                                            criterion = criterion, nu = nu, do_bootstrap = do_bootstrap, M = M)
  }
  cp_mosumvar$r <- r
  return(cp_mosumvar)
}
#mosumvar_factor()

## pooled forecast

pooled_forecast <- function(x, cp = 0, window_size = NULL, p = NULL,  n.ahead = 1, weights = "equal"){
  if(is.null(window_size)) window_size <- cp
  estim_start <- max(1, cp - window_size + 1)
  out <- matrix(0, n.ahead, ncol(x))
  if(weights == "equal"){
    weight_vec <- rep(1, cp+2 - estim_start)
  } else if(weights == "exp") {
    weight_vec <- exp(seq.int(cp+2 - estim_start) )
    weight_vec <- weight_vec * (cp+1 - estim_start) / sum(weight_vec) #normalise
  } else warning("incorrect weight specified")
  
  for (t in estim_start:(cp+1) ) {
    fit <- ar(x[t:nrow(x),], aic = is.null(p), order.max = p) 
    pred_data <- as.matrix(x[nrow(x) - (fit$order-1):0,])
    if(p==1) pred_data <- t(pred_data)
    out <- out + weight_vec[t-estim_start+1] * as.matrix(predict(fit, pred_data, n.ahead = n.ahead, se.fit=F))
  }
  out <- (out/(cp-estim_start))
  colnames(out) <- 1:ncol(x)
  return(out)
}
pooled_forecast(factors_, cp = 75,p=1)
