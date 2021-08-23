## mosumvar::mosumvar_factor
# mosumvar extension: MOSUM change point analysis for a factor model with low-dimensional VAR dynamics


#' Title
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
#' @examples
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
#mosumvar_factor()

## pooled forecast

#' Title
#'
#' @param x matrix of data with series as columns
#' @param cp integer (vector) of estimated change points
#' @param window_size integer size of window to estimate back over
#' @param p integer VAR order
#' @param n.ahead integer number of steps ahead to forecast
#' @param weights String of weight method to use
#' @param return_params Boolean, return model parameters and forecast variance
#'
#' @return matrix of forecasts, or list also containing model parameters and forecast variance
#' @export
#'
#' @examples
pooled_forecast <- function(x, cp = 0, window_size = NULL, p = NULL,  n.ahead = 1, weights = c("equal","exp","robust")[1], return_params =F){
  cp <- max(cp)
  if(is.null(window_size)) window_size <- cp

  #n <- nrow(x)
  estim_start <- max(1, cp - window_size + 1)
#  window_size <- min(window_size, n)
  out <- matrix(0, n.ahead, ncol(x))
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
    }  else warning("incorrect weight specified")
  } else weight_vec <- 1


  for (t in estim_start:(cp+1) ) {
    fit <- ar(x[t:nrow(x),], aic = is.null(p), order.max = p)
    if(return_params){
      if(t==estim_start) {out_params <- weight_vec[t-estim_start+1] *fit$ar;out_var <-  weight_vec[t-estim_start+1] *fit$var.pred}  else {out_params <- out_params + weight_vec[t-estim_start+1]*fit$ar; out_var <- weight_vec[t-estim_start+1]*out_var+ fit$var.pred}

    }
    pred_data <- as.matrix(x[nrow(x) - (fit$order-1):0,])
    if(p==1) pred_data <- t(pred_data)
    out <- out + weight_vec[t-estim_start+1] * as.matrix(predict(fit, pred_data, n.ahead = n.ahead, se.fit=F))
  }
  out <- out/n_size
  colnames(out) <- 1:ncol(x)

  if(return_params){
    out_params <- out_params/n_size
    out_var <- out_var/n_size
    out <- list(forecast = out, params = out_params, var = out_var)
  }

  return(out)
}
#pooled_forecast(factors_, cp = 0,p=1, return_params = T)



