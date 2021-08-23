memory.limit(size=8142*10)

library(Rcpp); library(RcppArmadillo)
library(factorcpt)

library(doParallel)
cl <- makeCluster(2); registerDoParallel(cl)
# stopCluster(cl)

source('/users/h/dropbox/bcf/codes/haeran/tmp.R')
source('/users/h/dropbox/packages/factorcpt/R/package.R')

source("c:\\users\\mahrc\\dropbox\\bcf\\codes\\haeran\\tmp.R")
source("c:\\users\\mahrc\\dropbox\\bcf\\codes\\haeran\\package.R")

#####################################################

n <- 100; T <- 500
B <- 200
sim <- 100

max.q <- max(round(20, sqrt(min(n, T))))
dw <- round(min(log(T)^2, T^(6/7)/4))
sig.lev <- .05 
scales <- -(1:floor(log(log(T, 2), 2)))
rule <- round(log(T, 2)/2)

r <- 5; ar <- TRUE; alpha <- .4
cor <- TRUE; H <- min(round(n/20), 10); rho <- 0.5; beta <- .2
# vep_sigma <- diag(rep(1, n)); for(j in 1:(n/2-1)) for(k in (j+1):(n/2)) vep_sigma[j, k] <- vep_sigma[k, j] <- (beta)^abs(j-k)
phi <- sqrt(r/(1-ar*alpha^2)/((1+2*beta^2*H)/(1-rho^2)*cor))
		
chp <- round(T*c(1/3, 1/2, 3/5, 4/5))

chi.type <- c(1, 2, 0, 3)
vep.type <- c(0, 0, 1, 0)
vep.diag <- TRUE
chi.cps <- chp[chi.type > 0]
vep.cps <- chp[vep.type > 0]

init.nts.seq <- seq(1, 2.5, by=.5)
sigma.seq <- sqrt(2)*seq(1, .25, by=-.25)
prop.seq <- seq(1, .25, by=-.25)

for(ii in 4){
for(ss in 4){
for(pp in 1:length(prop.seq)){

	oracle.chi.res <- chi.res <- oracle.vep.res <- vep.res <- dcbs.res <- array(0, dim=c(sim, T))
	select.q <- est.q <- rep(0, sim)
	for(i in 1:sim){
		sd <- power.sim.data(T=T, n=n, r=r, chp=chp,
			ar=ar, alpha=alpha, chi.type=chi.type, sigma=sigma.seq[ss], chi.prop=rep(prop.seq[pp], length(chp)), 
			cor=cor, H=H, vep.type=vep.type, vep.prop=rep(prop.seq[pp], length(chp)), rho=rho, beta=beta, phi=phi*init.nts.seq[ii])
		x <- sd$x
# matplot(t(x), type="l", col=1); matlines(t(sd$chi), col=2); matlines(t(sd$vep), col=3)

		gfm <- get.factor.model(x, bn.op = 2, max.q = max.q)
		est.q[i] <- q.hat <- max(1, gfm$q.hat)
		q.seq <- sort(unique(c(q.hat, max.q, round(seq(max(1, q.hat), max(min(n - 1, 2 * q.hat), max.q), length.out = 5)))), decreasing = FALSE)

	    est.cps <- cs.list <- list()
    	for (qq in q.seq) {
        	cs <- common.seg(gfm, scales=scales, q = qq, sig.lev = sig.lev, B = B, dw = dw, do.parallel = TRUE, rule = rule)
      	  	est.cps <- c(est.cps, list(cs$est.cps))
       	 	cs.list <- c(cs.list, list(cs))
    	}
    	qq <- max(which(unlist(lapply(est.cps, length)) == max(unlist(lapply(est.cps, length)))))
	    select.q[i] <- q <- q.seq[qq]
    	cs <- cs.list[[qq]]
	    chi.res[i, est.cps[[qq]]] <- 1	    
	    
	    ocs <- oracle.common.seg(sd$chi, q=2*r+1, scales=scales, sig.lev=sig.lev, rule=rule, B=B, dw=dw, do.parallel=TRUE)
		oracle.chi.res[i, ocs$est.cps] <- 1
		
    	is <- idio.seg(gfm, q = q, scales = scales, diag = TRUE, sig.lev = sig.lev, B = B,  dw = dw, do.parallel = TRUE, rule=rule)
    	vep.res[i, is$est.cps] <- 1
    	
    	ois <- oracle.idio.seg(sd$vep, scales=scales, diag=TRUE, sig.lev=sig.lev, rule=rule, B=B, dw=dw, do.parallel=TRUE)
    	oracle.vep.res[i, ois$est.cps] <- 1
    	
    	da <- dcbs.alg(x, sig.lev=sig.lev, scales=scales, dw=dw, B=B, rule=rule)
    	dcbs.res[i, da$est.cps] <- 1
      	   	   	
	}
	ls <- list(est.q=est.q, select.q=select.q, cr=chi.res, ocr=oracle.chi.res, vr=vep.res, ovr=oracle.vep.res, dr=dcbs.res)
###	
	save(ls, file=paste("c:\\users\\mahrc\\dropbox\\bcf_sim\\multi_n", n, "T", T, "_nts", ii, "sigma", ss, "prop", pp, "_chi.chp", length(chi.cps), "_vep.cps", length(vep.cps), "_r", r, ".RData", sep=""))
}}
}

######

ii <- 4

ss <- 1

res <- c()
dis <- 2*log(T)
for(ii in 1:length(init.nts.seq)){
for(pp in 1:length(prop.seq)){

name <- paste("/users/h/dropbox/bcf_sim/multi_n", n, "T", T, "_nts", ii, "sigma", ss, "prop", pp, "_chi.chp", length(chi.cps), "_vep.cps", length(vep.cps), "_r", r, ".RData", sep=""); load(file=name)

tmp <- c(mean(apply(ls$cr, 1, sum)), sd(apply(ls$cr, 1, sum)), 
sum(apply(ls$cr[, (chi.cps[1]-dis):(chi.cps[1]+dis)], 1, max)), 
sum(apply(ls$cr[, (chi.cps[2]-dis):(chi.cps[2]+dis)], 1, max)), 
sum(apply(ls$cr[, (chi.cps[3]-dis):(chi.cps[3]+dis)], 1, max)), 
mean(apply(ls$ocr, 1, sum)), sd(apply(ls$ocr, 1, sum)), 
sum(apply(ls$ocr[, (chi.cps[1]-dis):(chi.cps[1]+dis)], 1, max)), 
sum(apply(ls$ocr[, (chi.cps[2]-dis):(chi.cps[2]+dis)], 1, max)),
sum(apply(ls$ocr[, (chi.cps[3]-dis):(chi.cps[3]+dis)], 1, max)),
mean(apply(ls$vr, 1, sum)), sd(apply(ls$vr, 1, sum)), 
sum(apply(ls$vr[, (vep.cps[1]-dis):(vep.cps[1]+dis)], 1, max)), 
mean(apply(ls$ovr, 1, sum)), sd(apply(ls$ovr, 1, sum)), 
sum(apply(ls$ovr[, (vep.cps[1]-dis):(vep.cps[1]+dis)], 1, max)),
mean(apply(ls$dr, 1, sum)), sd(apply(ls$dr, 1, sum)), 
sum(apply(ls$dr[, (chi.cps[1]-dis):(chi.cps[1]+dis)], 1, max)), 
sum(apply(ls$dr[, (chi.cps[2]-dis):(chi.cps[2]+dis)], 1, max)),
sum(apply(ls$dr[, (chi.cps[3]-dis):(chi.cps[3]+dis)], 1, max)),
sum(apply(ls$dr[, (vep.cps[1]-dis):(vep.cps[1]+dis)], 1, max)))

res <- rbind(res, tmp)

}}

write(t(round(res, 2)), ncolumns=ncol(res), sep="\t", file='/users/h/dropbox/res.txt')

ss <- 1
for(ss in 1:length(sigma.seq)){

for(ii in 1:length(init.nts.seq)){
	
pdf(paste("/users/h/dropbox/bcf/notes/sim_fig/multi_n", n, "T", T, "_nts", ii, "sigma", ss, ".pdf", sep=""), width=8, height=8)

par(mfrow=c(4, 4), mar=c(2.5, 2.5, .5, .5))

for(pp in 1:length(prop.seq)){

name <- paste("/users/h/dropbox/bcf_sim/multi_n", n, "T", T, "_nts", ii, "sigma", ss, "prop", pp, "_chi.chp", length(chi.cps), "_vep.cps", length(vep.cps), "_r", r, ".RData", sep=""); load(file=name)

rng <- range(c(apply(ls$cr, 2, sum), apply(ls$ocr, 2, sum), apply(ls$vr, 2, sum), apply(ls$ovr, 2, sum), apply(ls$dr, 2, sum)))

plot(apply(ls$cr, 2, sum), xlab="", ylab="", type="s", lwd=3, ylim=rng); abline(v=vep.cps, col=4, lwd=2, lty=3); abline(v=chi.cps, col=2, lwd=2, lty=2)
plot(apply(ls$ocr, 2, sum), xlab="", ylab="", type="s", lwd=3, ylim=rng); abline(v=vep.cps, col=4, lwd=2, lty=3); abline(v=chi.cps, col=2, lwd=2, lty=2)
plot(apply(ls$vr, 2, sum), xlab="", ylab="", type="s", lwd=3, ylim=rng); abline(v=chi.cps, col=2, lwd=2, lty=2); abline(v=vep.cps, col=4, lwd=2, lty=3)
plot(apply(ls$ovr, 2, sum), xlab="", ylab="", type="s", lwd=3, ylim=rng); abline(v=chi.cps, col=2, lwd=2, lty=2); abline(v=vep.cps, col=4, lwd=2, lty=3)
#plot(apply(ls$dr, 2, sum), xlab="", ylab="", type="s", lwd=3, ylim=rng); abline(v=chi.cps, col=2, lwd=2, lty=2); abline(v=vep.cps, col=4, lwd=2, lty=3)

}
dev.off()
}
}




