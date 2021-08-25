# factorcpt code
##########################################

#' Fit a factor_model
#'
#' Extracts factors, loadings, and estimates the factor number
#'
#' @param x data matrix, with series as columns
#' @param max.q integer maximum factor number to consider
#' @param q integer factor number
#' @param bn Boolean, use information criteria of Bai/Ng (2002) for estimating q
#' @param bn.op integer from 1 to 5 determining which information criterion to use (see details)
#' @param normalisation Boolean, centre and scale x
#'
#' @return factor_model object
#' @export
#'
#' @examples factor_model(panel$panel)
factor_model <- function(x, max.q=NULL, q=NULL, bn=TRUE, bn.op=2, normalisation=TRUE){
  x <- t(x)
	T <- ncol(x); n <- nrow(x)
	cnt <- min(n, T)
	if(is.null(max.q)) max.q <- round(sqrt(cnt))
	if(bn.op > 5 || bn.op < 1 ) bn.op <- 2

	if(normalisation){
		mx <- matrix(rep(apply(x, 1, mean), each=T), byrow=TRUE, nrow=n)
		x <- x-mx
		sdx <- apply(x, 1, sd)
		x <- x/sdx
	  # x <- t(scale(t(x)) )
	} else{
		mx <- rep(0, n); sdx <- rep(1, n)
	}
	xx <- x%*%t(x)/T
	eig <- eigen(xx, symmetric=TRUE)
	lam <- eig$vectors[, 1:(cnt-1), drop=FALSE]*sqrt(n)
	f <- t(eig$vectors[, 1:(cnt-1), drop=FALSE])%*%x/sqrt(n)

	if(bn){
		ic <- rep(0, 1+max.q)
		ic[1] <- (bn.op <= 4)*log(mean(x^2)) + (bn.op==5)*mean(x^2)
		l <- 1
		while(l<=max.q){
			hchi <- lam[, 1:l, drop=FALSE]%*%f[1:l, , drop=FALSE]
			ic[l+1] <- (bn.op <= 4)*log(mean((x-hchi)^2)) +
				(bn.op==1)*l*(n+T)/(n*T)*log(n*T/(n+T)) +
				(bn.op==2)*l*(n+T)/(n*T)*log(cnt) +
				(bn.op==3)*l*log(cnt)/cnt +
				(bn.op==4)*l*((n+T-l)*log(n*T)/(n*T) + (n+T)/(n*T)*log(cnt))/2 +
				(bn.op==5)*(mean((x-hchi)^2)+l*mean((x-hchi)^2)*(n+T-l)*log(n*T)/(n*T))
			l <- l+1
		}
		q.hat <- which(ic==min(ic))-1
	} else{
		ic <- rep(0, max.q)
		q.hat <- q
	}
	lam.q <- t(lam)[1:q.hat,]
	f.q <- t(f)[,1:q.hat]
	residuals <- t(x) - (f.q) %*% lam.q
  out <- list(lam = t(lam), f = t(f), lam.q = lam.q, f.q = f.q, residuals = residuals,
              var_explained = sum(eig$values[1:q.hat])/sum(eig$values),
	            norm.x=t(x), q.hat=q.hat, max.q=max.q, ic=ic)
  attr(out, "class") <- "factor_model"
	return(out)
}

#' Predict from a factor_model and forecasted factors
#'
#' @param fm factor_model object
#' @param newdata matrix of (factor) data to predict from
#'
#' @return matrix of predictions
#' @export
#'
#' @examples fm <- factor_model(panel$panel)
#' ar_fm <- ar(fm$f.q, p=1)
#' ar_pred <- predict(ar_fm, fm$f.q, n.ahead = 3)
#' fm_pred <- predict(fm, ar_pred$pred)
predict.factor_model <- function(fm, newdata) {
  out <- (newdata) %*% fm$lam.q
  out <- as.ts(out)
  return(out)
}


#' Plot a factor_model object
#'
#' @param fm factor_model object
#'
#' @return a plot of the information criterion against q, and a plot of the factor series
#' @export
#'
#' @examples fm <- factor_model(panel$panel)
#' plot(fm)
plot.factor_model <- function(fm){
  par(mfrow=c(1,1))
  plot(fm$ic, ylab = "IC", xlab = "factor number"); abline(v = fm$q.hat, col = "red")
  if(fm$q.hat > 6) ts.plot(fm$f.q, main = "factors") else plot.ts(fm$f.q, main = "factors")
}

#' Summarise a factor_model object
#'
#' @param fm factor_model object
#'
#' @return summary of factor_model object
#' @export
#'
#' @examples fm <- factor_model(panel$panel)
#' summary(fm)
summary.factor_model <- function(fm){
  cat("Information criterion: ", fm$ic, "\n")
  cat("Factor number: ", fm$q.hat, "\n")
  cat("Maximum factor number: ", fm$max.q, "\n")
  cat("Variance explained: ", fm$var_explained)
}
