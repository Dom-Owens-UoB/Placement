## mosumvar::mosumvar_factor
# mosumvar extension: MOSUM change point analysis for a factor model with low-dimensional VAR dynamics


#' Change point detection for a factor model with VAR dynamics
#'
#' @param x matrix of data with series as columns
#' @param p integer VAR order
#' @param G integer MOSUM bandwidth
#' @param method String, which of "Wald" or "Score" statistics to use
#' @param estim String, estimation method
#' @param varEstim String, variance estimation method
#' @param alpha numeric, significance level
#' @param criterion String, location criterion
#' @param nu numeric, location criterion parameter
#' @param do_bootstrap Boolean, perform threshold bootstrap
#' @param M integer, bootstrap replications
#' @param univ Boolean, perform univariate dimension reduction
#' @param rm_cross_terms Boolean, remove cross terms under univ
#' @param global_resids Boolean, use global residuals under univ
#' @param r integer, factor number
#' @param r_max integer, maximum factor number to consider
#'
#' @return mosumvar object
#' @export
#'
#' @examples mosumvar_factor(panel$panel, p = 1, G = 24, method = "Score")
mosumvar_factor <- function(x, p, G, method = c("Wald", "Score")[1], estim = c("DiagC","DiagH")[1], varEstim = c("Local", "Global")[1], alpha = 0.05,
                            criterion = c("eps", "eta")[1], nu = 0.25, do_bootstrap = FALSE, M = 1000,
                            univ = FALSE, rm_cross_terms = TRUE, global_resids = TRUE,
                            r = NULL, r_max = NA) {
  fcpt_mosumvar <- factor_model(x) #factorcpt::get.factor.model(t(x))
  #fcpt_ar <- ar(t(fcpt_mosumvar$f[1:6,1:50]), aic = T, order.max = 4)
  if(is.null(r)) r <- max(fcpt_mosumvar$q.hat, 1)
  r <- min(r, r_max, na.rm = T)
  if(!univ) {cp_mosumvar <- mosumvar::mosumvar( x=(fcpt_mosumvar$f[,1:r]), p=p, G=G, method = method, estim =estim,
                                     varEstim = varEstim, alpha = alpha,
                                     criterion = criterion, nu = nu, do_bootstrap = do_bootstrap, M = M)
  #mosumvar::mosumvar(nowcast_mosumvar$factors$dynamic_factors, p = 2, G= 18, nu = .1)
  } else {cp_mosumvar <- mosumvar::mosum_univ( x=(fcpt_mosumvar$f[,1:r]), p=p, G=G, method = method, estim =estim,
                                            varEstim = varEstim, alpha = alpha, rm_cross_terms = rm_cross_terms, global_resids = global_resids,
                                            criterion = criterion, nu = nu, do_bootstrap = do_bootstrap, M = M)
  }
  cp_mosumvar$r <- r
  return(cp_mosumvar)
}



## pooled forecast

#' Change point based forecasting methods for factor series
#'
#' @param x matrix of data with series as columns
#' @param cp integer (vector) of estimated change points
#' @param window_size integer size of window to estimate back over
#' @param p integer VAR order
#' @param n.ahead integer number of steps ahead to forecast
#' @param weights String of weight method to use
#' @param return_model Boolean, return `ar` model
#'
#' @return list of class "pooled_forecast" containing matrix of forecasts, weight vector, and model
#' @export
#'
#' @examples fm <- factor_model(panel$panel)
#' pooled_forecast(fm$f.q)
pooled_forecast <- function(x, cp = 0, window_size = NULL, p = NULL,  n.ahead = 1, weights = c("equal","exp","robust","rolling")[1], return_model =F){
  cp <- max(cp)
  if(is.null(window_size)){
    if(weights != "rolling") window_size <- cp else {
      window_size <- max(p, 1) * ncol(x)^2
      warning("window_size set to ", window_size)
    }
  }

  #n <- nrow(x)
  estim_start <- max(1, cp - window_size + 1)
#  window_size <- min(window_size, n)
  fc <- matrix(0, n.ahead, ncol(x))
  if(cp != 0){
    n_size <- cp+2 - estim_start
    if(weights == "equal"){
    weight_vec <- rep(1, n_size)
    } else if(weights == "exp") {
      weight_vec <- exp(seq.int(n_size) )
      weight_vec <- weight_vec * (n_size-1) / sum(weight_vec) #normalise
    } else if (weights == "robust"){
      weight_vec <- -log(1 - seq.int(n_size)/n_size )#/(n_size-1)
      weight_vec[n_size] <- weight_vec[n_size-1]#log(n_size)/(n_size-1)##
      weight_vec <- weight_vec * (n_size-1)/ sum(weight_vec) #* (cp+1 - estim_start)
    }   else  warning("incorrect weight specified")
  } else {
    n_size <- 1
    weight_vec <- 1
  }

  if(weights != "rolling"){
  for (t in estim_start:(cp+1) ) {
    fit <- ar(x[t:nrow(x),], aic = is.null(p), order.max = p)
    if(return_model){
      if(t==estim_start) {
        out_params <- weight_vec[t-estim_start+1] *fit$ar;
        out_var <-  weight_vec[t-estim_start+1] *fit$var.pred
        out_mean <- fit$x.mean}  else {
        out_params <- out_params + weight_vec[t-estim_start+1]*fit$ar;
        out_var <- weight_vec[t-estim_start+1]*out_var+ fit$var.pred
        out_mean <- weight_vec[t-estim_start+1]*out_mean + fit$x.mean}
    }
    p <- fit$order
    pred_data <- as.matrix(x[nrow(x) - (p-1):0,])
    if(p==1) pred_data <- t(pred_data)
    fc <- fc + weight_vec[t-estim_start+1] * as.matrix(predict(fit, pred_data, n.ahead = n.ahead, se.fit=F))
    }
  fc <- fc/n_size
  } else {
    n_size <- 1
    weight_vec <- 1
    fit <- ar(x[(nrow(x) - window_size + 1):nrow(x),], aic = is.null(p), order.max = sqrt(window_size))
    pred_data <- as.matrix(x[nrow(x) - (fit$order-1):0,])
    p <- fit$order
    if(p==1) pred_data <- t(pred_data)
    fc <- as.matrix(predict(fit, pred_data, n.ahead = n.ahead, se.fit=F))
    out_params <- fit$ar
    out_var <- fit$var.pred
    out_mean <- fit$x.mean
  }
  colnames(fc) <- 1:ncol(x)


  if(return_model){
    model <- fit#list()
    model$ar <- out_params/n_size
    model$var.pred <- out_var/n_size
    model$x.mean <- fit$x.mean
    model$partialacf <- NULL
    resid <- x[-(1:p),]
    for (ii in 1:p) {
      resid <- resid - (x[-c(0:(p-ii), nrow(x)+(1-ii):0) ,]) %*% t(model$ar[ii,,])
    }
    model$resid <- resid
  } else model <- list()

  out <- list(forecast = fc, weights = weight_vec,  model = model)
  attr(out, "class") <- "pooled_forecast"
  return(out)
}
#pooled_forecast(factors_, cp = 0,p=1, return_params = T)



