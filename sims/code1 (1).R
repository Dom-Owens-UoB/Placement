power.sim.data <- function(T, n, r, chp, burnin=NULL,
	ar=TRUE, alpha=.4, chi.type=0, sigma=1, chi.prop=0, 
	cor=TRUE, H=NULL, vep.type=0, vep.prop=0, rho=.5, beta=.2, phi=1){

	q <- r
	if(is.null(burnin)) burnin <- max(50, round(T*.1))
	if(is.null(H)) H <- round(min(n/20, 10))
	lam <- matrix(rnorm(n*q), nrow=n)
	f <- matrix(rnorm(q*(T+burnin)), nrow=q)	
	B <- length(chp); brks <- c(0, chp+burnin, T+burnin)
	
	if(ar){
		ar.coef <- alpha-((1:q)-1)*.05
		for(t in (brks[1]+2):brks[2]) f[, t] <- ar.coef*f[, t-1] + f[, t]	
		for(b in 1:B){
			if(chi.type[b]==2) ar.coef <- -ar.coef
			for(t in (brks[b+1]+1):brks[b+2]) f[, t] <- ar.coef*f[, t-1] + f[, t]	
		} 
	}
	chi <- lam%*%f
	if(sum(chi.type==1)+sum(chi.type==3) > 0){
		for(b in 1:B){
			if(chi.type[b]==1){
				del <- matrix(0, nrow=n, ncol=q)
				del[sample(n, round(n*chi.prop[b])), 1:q] <- matrix(rnorm(round(n*chi.prop[b])*q)*sigma, ncol=q)
				lam <- lam+del
				chi[, (brks[b+1]+1):(T+burnin)] <- chi[, (brks[b+1]+1):(T+burnin)]+del%*%f[, (brks[b+1]+1):(T+burnin), drop=FALSE]
			}
			if(chi.type[b]==3){
				del <- matrix(0, nrow=n, ncol=1)
				del[sample(n, round(n*chi.prop[b])), 1] <- rnorm(round(n*chi.prop[b]))*sigma
				lam <- cbind(lam, del)
				f <- rbind(f, 0); q <- q+1; ar.coef <- c(ar.coef, alpha)
				for(t in (brks[b+1]+1):(T+burnin)) f[q, t] <- ar.coef[q]*f[q, t-1]+rnorm(1)
				chi[, (brks[b+1]+1):(T+burnin)] <- chi[, (brks[b+1]+1):(T+burnin)]+del%*%f[q, (brks[b+1]+1):(T+burnin), drop=FALSE]
			}			
		} 
	} 
	chi <- chi[, (1+burnin):(T+burnin)]
			
	v <- vep <- matrix(rnorm((T+burnin)*n), nrow=n)
	rng.mat <- matrix(0, n, n)
	if(cor){
		rho.seq <- runif(n, -1, 1)*rho
		beta.seq <- beta*sample(c(-1, 1), n, replace=TRUE)
		for(i in 1:n) rng.mat[i, setdiff(rep(1:n, 3)[n+(i-H):(i+H)], i)] <- 2
	} else{
		beta.seq <- rho.seq <- rep(0, n)
	}
	for(b in 0:B){
		if(b>0 && vep.type[b]==1){
			ind <- sample(n, round(n*vep.prop[b]))
			rho.seq[ind] <- -rho.seq[ind]
		}
		if(b>0 && vep.type[b]==2){
			ind <- sample(n, round(n*vep.prop[b]))
#			beta.seq[ind] <- beta.seq[ind]*sample(c(sigma, 1/sigma), length(ind), replace=TRUE)
			for(i in ind) rng.mat[i, setdiff(rep(1:n, 3)[n+(i-2*H):(i+2*H)], i)] <- 1
		}
		for(i in 1:n){
			rng <- which(rng.mat[i,]==1) #setdiff(rep(1:n, 3)[n+(i-H):(i+H)], i)
			vep[i, (brks[b+1]+1):brks[b+2]] <- v[i, (brks[b+1]+1):brks[b+2]]+beta.seq[i]*apply(v[rng, (brks[b+1]+1):brks[b+2]], 2, sum)
		}
		for(t in max(2, brks[b+1]+1):brks[b+2]) vep[, t] <- rho.seq*vep[, t-1] + vep[, t]	
	}
	vep <- vep*phi; vep <- vep[, (1+burnin):(T+burnin)]
	x <- chi + vep
	
	return(list(x=x, chi=chi, vep=vep, phi=phi, lam=lam, f=f[,-(1:burnin)]))	
}

#####

get.factor.model <- function(x, max.q=NULL, q=NULL, bn=TRUE, bn.op=2, normalisation=TRUE){
	T <- ncol(x); n <- nrow(x)
	cnt <- min(n, T)
	if(is.null(max.q)) max.q <- round(sqrt(cnt))
	if(is.null(q)) q <- max.q

	if(normalisation){
		mx <- matrix(rep(apply(x, 1, mean), each=T), byrow=TRUE, nrow=n)
		x <- x-mx
		sdx <- apply(x, 1, sd)
		x <- x/sdx
	} else{
		mx <- rep(0, n); sdx <- rep(1, n)
	}
	xx <- x%*%t(x)/T
	eig <- eigen(xx, symmetric=TRUE)
	lam <- eig$vectors[, 1:(cnt-1), drop=FALSE]*sqrt(n)
	f <- t(eig$vectors[, 1:(cnt-1), drop=FALSE])%*%x/sqrt(n)
	sigma.hat <- mean((x-lam%*%f)^2)

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
	
	return(list(lam = lam, f = f, norm.x=x, q.hat=q.hat, max.q=max.q, ic=ic))
}


make.tree <- function(y, op, dw, rule=NULL){
	d <- nrow(y); len <- ncol(y)
	if(is.null(rule)) rule <- round(log(len, 2)/2)
	tree <- list(matrix(0, 6, 1))
	mat <- c()

	fd <- func_density(y, 1)
	stat <- fd$res[, op]
	test.stat <- max(stat[-c((1:dw), (len-dw):len)])
	hat.chp <- min(which(stat==test.stat))
	
	tree[[1]][1, 1] <- 1
	tree[[1]][2, 1] <- 1
	tree[[1]][3, 1] <- hat.chp
	tree[[1]][4, 1] <- len
	tree[[1]][5, 1] <- 0
	tree[[1]][6, 1] <- test.stat
	mat <- cbind(mat, c(tree[[1]][-5, ], 1, 1))
		
	j <- 1
	while(length(tree)==j & j < rule){
		npc <- dim(tree[[j]])[2]
		if(sum(tree[[j]][4,]-tree[[j]][2,]-rep(4*dw, npc)>0)){
			ncc <- 0; i <- 1
			while(i <= npc){
				if(tree[[j]][3, i]-tree[[j]][2, i]+1>4*dw){
					s <- tree[[j]][2, i]; e <- tree[[j]][3, i]
					fd <- func_density(y[, s:e], 1)
					stat <- fd$res[, op]
					test.stat <- max(stat[-c((1:dw), (e-s+1-dw):(e-s+1))])
					hat.chp <- s+min(which(stat==test.stat))-1
						
					if(length(tree)==j) tree <- c(tree, list(matrix(0, 6, 0)))
					ncc <- ncc+1
					tree[[j+1]] <- matrix(c(tree[[j+1]], matrix(0, 6, 1)), 6, ncc)
					tree[[j+1]][1, ncc] <- 2*tree[[j]][1, i]-1
					tree[[j+1]][2, ncc] <- s
					tree[[j+1]][3, ncc] <- hat.chp
					tree[[j+1]][4, ncc] <- e
					tree[[j+1]][5, ncc] <- 0
					tree[[j+1]][6, ncc] <- test.stat
					mat <- cbind(mat, c(tree[[j+1]][-5, ncc], j+1, ncc))
				}
				if(tree[[j]][4, i]-tree[[j]][3, i]>4*dw){
					s <- tree[[j]][3, i]+1; e <- tree[[j]][4, i]
					fd <- func_density(y[, s:e], 1)
					stat <- fd$res[, op]
					test.stat <- max(stat[-c((1:dw), (e-s+1-dw):(e-s+1))])
					hat.chp <- s+min(which(stat==test.stat))-1

					if(length(tree)==j) tree <- c(tree, list(matrix(0, 6, 0)))
					ncc <- ncc+1
					tree[[j+1]] <- matrix(c(tree[[j+1]], matrix(0, 6, 1)), 6, ncc)
					tree[[j+1]][1, ncc] <- 2*tree[[j]][1, i]
					tree[[j+1]][2, ncc] <- s
					tree[[j+1]][3, ncc] <- hat.chp
					tree[[j+1]][4, ncc] <- e
					tree[[j+1]][5, ncc] <- 0
					tree[[j+1]][6, ncc] <- test.stat
					mat <- cbind(mat, c(tree[[j+1]][-5, ncc], j+1, ncc))
				}
				i <- i+1
			}
			j <- j+1
		} else{ 
			break
		}
	}
	list(tree=tree, mat=mat)
}

######

chi.seg <- function(gfm, q, scales=NULL, p=NULL, op=2, alpha=0.05, rule=NULL, B=200, M=NULL, dw=NULL, mby=NULL, tby=NULL, seg.op=1){
	lam <- gfm$lam; f <- gfm$f; nx <- gfm$norm.x
	T <- dim(nx)[2]
	if(is.null(p) || (p<=0 | p>=1)) p.seq <-  apply(f[1:q,,drop=FALSE], 1, function(z){g <- get.gg(z, M=M); min(.5, ((g[2]/g[1])^2)^(-1/3)*T^(-1/5))}) else p.seq <- rep(p, q)
	if(is.null(scales)) scales <- -(1:floor(log(log(T, 2), 2)))
	
	y <- c()
	hat.chi <- lam[, 1:q, drop=FALSE]%*%f[1:q, , drop=FALSE]
	for(sc in scales){
		cc <- func_coef(hat.chi, sc)
		if(sc>min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop=FALSE]
		chi.input <- t(func_input_on(cc, 0, 1)$input)
		y <- rbind(y, chi.input)
	}

	d <- nrow(y); len <- ncol(y)	
	ls <- length(scales)
	if(is.null(dw)) dw <- floor(2*log(len))
	if(is.null(mby)) mby <- ceiling(log(d))
	if(is.null(tby)) tby <- ceiling(log(len))
	if(is.null(rule)) rule <- round(log(len, 2)/2)

	mat <- make.tree(y, op, dw, rule)$mat
	null.stat <- chi.null.stat(mat, op, p.seq, B, q, lam, f, scales, mby, tby)

	k <- 1; pos.seq <- c()
	while(k <= ncol(mat)){
		pos <- mean(mat[5, k]>null.stat[k, ])
		pos.seq <- c(pos.seq, pos)
		if(pos < 1-alpha-(seg.op==1)*alpha){
			l <- mat[6, k]+1; rm.ind <- mat[1, k]
			while(l <= max(mat[6, ])){
				ind <- which(mat[6, ]==l & is.element(mat[1, ], c(rm.ind*2-1, rm.ind*2)))
				if(length(ind) > 0){
					rm.ind <- mat[1, ind]
					mat <- mat[, -ind, drop=FALSE]
					l <- l+1
				} else{
					break
				}
			}
		} 
		k <- k+1
	}
	tree <- list(matrix(0, 6, 1))
	mat <- rbind(mat, pos.seq); mat <- mat[, pos.seq >= 1-alpha, drop=FALSE]
	if(dim(mat)[2] > 0){
		for(l in 1:length(unique(mat[6, ]))){
			j <- unique(mat[6, ])[l]
			for(ncc in 1:sum(mat[6, ]==j)){
				k <- sum(mat[6, ] < j) + ncc
				if(length(tree)<j) tree <- c(tree, list(matrix(0, 6, 0)))
				tree[[l]] <- matrix(c(tree[[l]], matrix(0, 6, 1)), 6, ncc)
				tree[[l]][1, ncc] <- mat[1, k]
				tree[[l]][2, ncc] <- mat[2, k]
				tree[[l]][3, ncc] <- mat[3, k]
				tree[[l]][4, ncc] <- mat[4, k]					
				tree[[l]][5, ncc] <- mat[8, k]
				tree[[l]][6, ncc] <- mat[5, k]
			}
		}
	}
	est.cps <- sort(mat[3, ])

	ls <- list(op=op, tree=tree, est.cps=est.cps+2^(-min(scales))-1, dw=dw, p=p.seq)
	return(ls)
}

chi.null.stat <- function(mat, op, p.seq, B, q, lam, f, scales, mby, tby){
	len <- ncol(f)
	null.stat <- foreach(l=iter(1:B), .combine=cbind, .packages=c("Rcpp", "RcppArmadillo", "dcbsby")) %dopar% {
		boot.chi <- 0
		for(qq in 1:q){
			ind <- c()
			while(length(ind)<len){ L <- rgeom(1, p.seq[qq]); I <- sample(len, 1); ind <- c(ind, rep(1:len, 1+ceiling(L/len))[I:(I+L-1)]) }
			ind <- ind[1:len]
			boot.chi <- boot.chi + lam[, qq, drop=FALSE]%*%f[qq, ind, drop=FALSE]				
		}
		by <- c()
		for(sc in scales){
			cc <- func_coef(boot.chi, sc)
			if(sc>min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop=FALSE]
			by <- rbind(by, t(func_input_on(cc, 0, 1)$input))
		}
		tmp <- c()
		for(i in 1:ncol(mat)){
			s <- mat[2, i]; e <- mat[4, i]
			bfd <- func_density_by(by[, s:e], mby, tby, 1)
			tmp <- c(tmp, max(bfd$res[, op]))
			rm(bfd)
		}
		rm(by)
		tmp
	}
	null.stat
}


#####

vep.seg <- function(gfm, q, scales=NULL, diag=TRUE, p=NULL, op=2, alpha=0.05, rule=NULL, B=200, M=NULL, dw=NULL, mby=NULL, tby=NULL, seg.op=1){
	lam <- gfm$lam[,1:q,drop=FALSE]; f <- gfm$f[1:q,,drop=FALSE]; nx <- gfm$norm.x
	hat.vep <- nx - lam%*%f
	T <- dim(nx)[2]
	if(is.null(p) || (p<=0 | p>=1)) p <- min(.5, 1/mean(apply(hat.vep, 1, function(z){g <- get.gg(z, M=M); ((g[2]/g[1])^2)^(1/3)*T^(1/5)})))
	if(is.null(scales)) scales <- -(1:floor(log(log(T, 2), 2)))
	
	y <- NULL
	for(sc in scales){
		cc <- func_coef(hat.vep, sc)
		if(sc>min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop=FALSE]
		if(diag){
			y <- rbind(y, t(func_input_on(cc, 0, 1)$input))
		} else{
			sgn <- sign(cc%*%t(cc))
			y <- rbind(y, t(func_input(cc, sgn, 0, 1)))
		}
	}	
	d <- nrow(y); len <- ncol(y)	
	ls <- length(scales)
	if(is.null(dw)) dw <- floor(2*log(len))
	if(is.null(mby)) mby <- ceiling(log(d))
	if(is.null(tby)) tby <- ceiling(log(len))
	if(is.null(rule)) rule <- round(log(len, 2)/2)

	mat <- make.tree(y, op, dw, rule)$mat

	null.stat <- vep.null.stat(mat, op, p, B, q, hat.vep, diag, scales, mby, tby)

	k <- 1; pos.seq <- c()
	while(k <= ncol(mat)){
		pos <- mean(mat[5, k]>null.stat[k, ])
		pos.seq <- c(pos.seq, pos)
		if(pos < 1-alpha-(seg.op==1)*alpha){
			l <- mat[6, k]+1; rm.ind <- mat[1, k]
			while(l <= max(mat[6, ])){
				ind <- which(mat[6, ]==l & is.element(mat[1, ], c(rm.ind*2-1, rm.ind*2)))
				if(length(ind) > 0){
					rm.ind <- mat[1, ind]
					mat <- mat[, -ind, drop=FALSE]
					l <- l+1
				} else{
					break
				}
			}
		} 
		k <- k+1
	}
	tree <- list(matrix(0, 6, 1))
	mat <- rbind(mat, pos.seq); mat <- mat[, pos.seq >= 1-alpha, drop=FALSE]
	if(dim(mat)[2] > 0){
		for(l in 1:length(unique(mat[6, ]))){
			j <- unique(mat[6, ])[l]
			for(ncc in 1:sum(mat[6, ]==j)){
				k <- sum(mat[6, ] < j) + ncc
				if(length(tree)<j) tree <- c(tree, list(matrix(0, 6, 0)))
				tree[[l]] <- matrix(c(tree[[l]], matrix(0, 6, 1)), 6, ncc)
				tree[[l]][1, ncc] <- mat[1, k]
				tree[[l]][2, ncc] <- mat[2, k]
				tree[[l]][3, ncc] <- mat[3, k]
				tree[[l]][4, ncc] <- mat[4, k]					
				tree[[l]][5, ncc] <- mat[8, k]
				tree[[l]][6, ncc] <- mat[5, k]
			}
		}
	}
	est.cps <- sort(mat[3, ])

	ls <- list(op=op, tree=tree, est.cps=est.cps+2^(-min(scales))-1, dw=dw, p=p)
	return(ls)
}

vep.null.stat <- function(mat, op, p, B, q, hat.vep, diag, scales, mby, tby){
	len <- ncol(hat.vep)
	null.stat <- foreach(l=iter(1:B), .combine=cbind, .packages=c("Rcpp", "RcppArmadillo", "dcbsby")) %dopar% {
		ind <- c()
		while(length(ind)<len){
			L <- rgeom(1, p); I <- sample(len, 1)
			ind <- c(ind, rep(1:len, 1+ceiling(L/len))[I:(I+L-1)])
		}
		ind <- ind[1:len]
		boot.vep <- hat.vep[, ind]
		by <- NULL
		for(sc in scales){
			cc <- func_coef(boot.vep, sc)
			if(sc>min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop=FALSE]
			if(diag){
				by <- rbind(by, t(func_input_on(cc, 0, 1)$input))
			} else{
				sgn <- sign(cc%*%t(cc))
				by <- rbind(by, t(func_input(cc, sgn, 0, 1)))
			}
		}
		tmp <- c()
		for(i in 1:ncol(mat)){
			s <- mat[2, i]; e <- mat[4, i]; k <- 1
			bfd <- func_density_by(by[, s:e], mby, tby, 1)
			tmp <- c(tmp, max(bfd$res[, op]))
			rm(bfd)
		}
		rm(by)
		tmp
	}
	null.stat
}

#####

chi.detect <- function(gfm, q.seq, scales=NULL, p=NULL, B=200, M=NULL, dw=NULL, mby=NULL, tby=NULL){
	lam <- gfm$lam; f <- gfm$f
	T <- dim(f)[2]
	if(is.null(p) || (p<=0 | p>=1)) p.seq <-  apply(f[1:max(q.seq),,drop=FALSE], 1, function(z){g <- get.gg(z, M=M); min(.5, ((g[2]/g[1])^2)^(-1/3)*T^(-1/5))}) else p.seq <- rep(p, max(q.seq))
	if(is.null(scales)) scales <- -(1:floor(log(log(T, 2), 2)))
	
	y <- list()
	for(r in 1:length(q.seq)){
		q <- q.seq[r]
		z <- c()
		hat.chi <- lam[, 1:q, drop=FALSE]%*%f[1:q, , drop=FALSE]
		for(sc in scales){
			cc <- func_coef(hat.chi, sc)
			if(sc > min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop=FALSE]
			chi.input <- t(func_input_on(cc, 0, 1)$input)
			z <- rbind(z, chi.input)
		}
		y <- c(y, list(z))
	}	
	d <- nrow(y[[1]]); len <- ncol(y[[1]])	
	ls <- length(scales)
	if(is.null(dw)) dw <- floor(2*log(len))
	if(is.null(mby)) mby <- ceiling(log(d))
	if(is.null(tby)) tby <- ceiling(log(len))

	mat.list <- list()
	for(r in 1:length(q.seq)){
		fd <- func_density(y[[r]], 1)
		res <- fd$res[-c((1:dw), (len-dw):len), ]
		mat <- rbind(apply(res, 2, max), dw+apply(res, 2, which.max))
		mat.list <- c(mat.list, list(mat))
	}
	null.stat <- chi.detect.null.stat(p.seq, B, q.seq, lam, f, scales, mby, tby)
	
	res.list <- list()
	for(r in 1:length(q.seq)){
		mat <- mat.list[[r]]
		ns <- null.stat[(r-1)*ncol(mat)+1:ncol(mat), ]
		pos <- apply(matrix(mat[1, ], ncol=ncol(ns), nrow=ncol(mat), byrow=FALSE)>ns, 1, mean)
		res.list <- c(res.list, list(rbind(mat, pos)))
	}
	return(list(res.list=res.list, p=p.seq))
}

chi.detect.null.stat <- function(p.seq, B, q.seq, lam, f, scales, mby, tby){
	len <- ncol(f)
	null.stat <- foreach(l=iter(1:B), .combine=cbind, .packages=c("Rcpp", "RcppArmadillo", "dcbsby")) %dopar% {
		boot.chi <- 0
		tmp <- c()
		for(q in 1:max(q.seq)){
			ind <- c()
			while(length(ind)<len){ 
				L <- rgeom(1, p.seq[q]); I <- sample(len, 1); ind <- c(ind, rep(1:len, 1+ceiling(L/len))[I:(I+L-1)]) 
			}
			ind <- ind[1:len]
			boot.chi <- boot.chi + lam[, q, drop=FALSE]%*%f[q, ind, drop=FALSE]				
			if(is.element(q, q.seq)){
				by <- c()
				for(sc in scales){
					cc <- func_coef(boot.chi, sc)
					if(sc > min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop=FALSE]
					by <- rbind(by, t(func_input_on(cc, 0, 1)$input))
				}
				bfd <- func_density_by(by, mby, tby, 1)
				tmp <- c(tmp, apply(bfd$res, 2, max))
				rm(by, bfd)
			}
		}
		tmp
	}
	null.stat
}

vep.detect <- function(gfm, q.seq0, scales=NULL, p=NULL, diag=TRUE, B=200, M=NULL, dw=NULL,mby=NULL, tby=NULL){
	lam <- gfm$lam[,1:max(q.seq0),drop=FALSE]; f <- gfm$f[1:max(q.seq0),,drop=FALSE]; nx <- gfm$norm.x
	T <- dim(nx)[2]
	if(is.null(p) || (p<=0 | p>=1)) p <- min(1/2, 1/mean(apply(nx - lam[,1:max(q.seq0),drop=FALSE]%*%f[1:max(q.seq0),,drop=FALSE], 1, function(z){g <- get.gg(z, M=M); ((g[2]/g[1])^2)^(1/3)*T^(1/5)})))
	if(is.null(scales)) scales <- -(1:floor(log(log(T, 2), 2)))	
	q.seq <- sort(unique(q.seq0))

	y <- list()
	for(r in 1:length(q.seq)){
		q <- q.seq[r]; hat.vep <- nx - lam[,1:q,drop=FALSE]%*%f[1:q,,drop=FALSE]
		z <- c()
		for(sc in scales){
			cc <- func_coef(hat.vep, sc)
			if(sc > min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop=FALSE]
			if(diag){
				z <- rbind(z, t(func_input_on(cc, 0, 1)$input))
			} else{
				sgn <- sign(cc%*%t(cc))
				z <- rbind(z, t(func_input(cc, sgn, 0, 1)))
			}
		}
		y <- c(y, list(z))
	}	
	d <- nrow(y[[1]]); len <- ncol(y[[1]])	
	ls <- length(scales)
	if(is.null(dw)) dw <- floor(2*log(len))
	if(is.null(mby)) mby <- ceiling(log(d))
	if(is.null(tby)) tby <- ceiling(log(len))

	mat <- matrix(0, nrow=2, ncol=5)
	for(r in 1:length(q.seq)){
		op <- which(q.seq0==q.seq[r])
		fd <- func_density(y[[r]], 1)
		res <- fd$res[-c((1:dw), (len-dw):len), c(1:2, 4:6)[op], drop=FALSE]
		mat[, op] <- rbind(apply(res, 2, max), dw+apply(res, 2, which.max))
	}
	null.stat <- vep.detect.null.stat(p, B, q.seq, lam, f, nx, diag, scales, mby, tby)	

	res.list <- list()
	for(op in 1:5){
		r <- (1:length(q.seq))[q.seq0[op]==q.seq]
		ns <- null.stat[length(q.seq0)*(r-1)+op, ]
		pos <- mean(mat[1, op]>ns)
		res.list <- c(res.list, list(c(mat[, op], pos)))
	}

	return(list(res.list=res.list, p=p))
}

vep.detect.null.stat <- function(p, B, q.seq, lam, f, nx, diag, scales, mby, tby){
	len <- ncol(f)
	null.stat <- foreach(l=iter(1:B), .combine=cbind, .packages=c("Rcpp", "RcppArmadillo", "dcbsby")) %dopar% {
		tmp <- c()
		for(q in q.seq){
			hat.vep <- nx - lam[, 1:q, drop=FALSE]%*%f[1:q, , drop=FALSE]
			ind <- c()
			while(length(ind)<len){ 
				L <- rgeom(1, p); I <- sample(len, 1); ind <- c(ind, rep(1:len, 1+ceiling(L/len))[I:(I+L-1)]) 
			}
			ind <- ind[1:len]
			boot.vep <- hat.vep[, ind]	
			by <- NULL
			for(sc in scales){
				cc <- func_coef(boot.vep, sc)
				if(sc > min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop=FALSE]
				if(diag){
					by <- rbind(by, t(func_input_on(cc, 0, 1)$input))
				} else{
					sgn <- sign(cc%*%t(cc))
					by <- rbind(by, t(func_input(cc, sgn, 0, 1)))
				}
			}
			bfd <- func_density_by(by, mby, tby, 1)
			tmp <- c(tmp, apply(bfd$res[, -3], 2, max))
			rm(bfd, by)	
		}
		tmp
	}
	null.stat
}


calc.bias <- function(cps, hcps, len){
	a <- length(cps); b <- length(hcps)
	if(b==0) return(rep(len, a))
	mat <- abs(matrix(rep(cps, b), nrow=a, ncol=b, byrow=FALSE)-matrix(rep(hcps, a), nrow=a, ncol=b, byrow=TRUE))
	ind <- apply(mat, 1, which.min)
	res <- apply(mat, 1, min)
	if(sum(duplicated(ind))>0){
		for(l in 1:b){
			if(sum(ind==l)==1) next
			if(sum(ind==l)>1) res[setdiff(which(ind==l), which(ind==l)[which.min(res[ind==l])])] <- len
		}
	}
	return(res)
}

calc.gof.chi <- function(bx, bn.op=2, gfm, q, hcps){
	hat.chi <- gfm$lam[, 1:gfm$max.q, drop=FALSE]%*%gfm$f[1:gfm$max.q, , drop=FALSE]
	brks <- c(0, hcps, T)
	new.chi <- hat.chi*0
	for(b in 1:(length(brks)-1)){
		int <- (brks[b]+1):brks[b+1]
		ngfm <- get.factor.model(t(gfm$norm.x[, int]), bn.op = bn.op, max.q = gfm$q.hat, normalisation=FALSE)
		new.chi[, int] <- t(ngfm$lam[, 1:ngfm$q.hat, drop=FALSE]%*%ngfm$f[1:ngfm$q.hat, , drop=FALSE])
	}
	mean((hat.chi-new.chi)^2)
}

tri.kern <- function(len, h){
    filter <- rep(0, h+1)
    i <- 0
    while (i <= h) {
        u <- i/h
        if (u < 1/2) 
            filter[i+1] <- 1
        if (u >= 1/2 & u < 1) 
            filter[i+1] <- 2 * (1 - u)
        if (u > 1) 
            break
        i <- i + 1
    }
    filter
}

get.gg <- function(z, M=NULL, C=2, max.K=5){
	len <- length(z)
	max.K <- max(max.K, sqrt(log(len)))
	acv <- acf(z, type="covariance", lag.max=len-1, plot=FALSE)$acf[,,1]
	if(is.null(M)){
		l <- 1; ind <- 0
		while(l < sqrt(len)){
			if(abs(acv[l+1])/acv[1] < C*sqrt(log(len)/len)){
				ind <- ind+1
			} else{
				if(ind>0) ind <- 0
			}
			if(ind==max.K) break
			l <- l+1
		}
		lam <- max(1, l-max.K); M <- 2*lam
	}
	k <- tri.kern(len, M)	
	c(acv[1]+2*sum(k[-1]*acv[2:(M+1)]), 2*sum(k[-1]*(1:M)*acv[2:(M+1)]))
}