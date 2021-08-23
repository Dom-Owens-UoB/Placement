memory.limit(size=8142*10)

library(Rcpp); library(RcppArmadillo)
library(dcbsby)
# library(dfmboot)
# library(abind)
# library(fields)

library(doParallel)
cl <- makeCluster(2); registerDoParallel(cl)
# stopCluster(cl)

source("c:\\users\\mahrc\\dropbox\\bcf\\codes\\haeran\\code1.R")

source('/users/h/dropbox/bcf/codes/haeran/code1.R')

#####################################################

n <- 100; T <- 200
Bsim <- 200
sim <- 100
dw <- 5
thr <- 1-.05
M <- NULL #round(sqrt(T))

max.q <- max(20, round(sqrt(min(n, T))))
bn.op <- 2

scales <- -seq(1:floor(log(log(T, 2), 2))); ls <- length(scales)

r <- 5; ar <- TRUE; alpha <- .4
cor <- TRUE; H <- min(round(n/20), 10); rho <- 0.5; beta <- .2
# vep_sigma <- diag(rep(1, n)); for(j in 1:(n/2-1)) for(k in (j+1):(n/2)) vep_sigma[j, k] <- vep_sigma[k, j] <- (beta)^abs(j-k)
phi <- sqrt(r/(1-ar*alpha^2)/((1+2*beta^2*H)/(1-rho^2)*cor))
	
chp <- round(T*c(1/3))
chi.type <- 0
vep.type <- 2
vep.diag <- TRUE

init.nts.seq <- seq(1, 2.5, by=.5)
sigma.seq <- sqrt(2)*seq(1, .25, length.out=4)
prop.seq <- seq(1, .25, by=-.25)

if(chi.type==1) chi.res.array <- vep.res.array <- array(0, dim=c(length(init.nts.seq), length(sigma.seq), length(prop.seq), 5))
if(chi.type==2) chi.res.array <- vep.res.array <- array(0, dim=c(length(init.nts.seq), 1, 1, 5))
if(chi.type==3) chi.res.array <- vep.res.array <- array(0, dim=c(length(init.nts.seq), length(sigma.seq[1]), length(prop.seq), 5))

if(vep.type==1) chi.res.array <- vep.res.array <- array(0, dim=c(length(init.nts.seq), 1, length(prop.seq), 5))
if(vep.type==2) chi.res.array <- vep.res.array <- array(0, dim=c(length(init.nts.seq), 1, length(prop.seq), 5))

for(ii in 1:dim(chi.res.array)[1]){
	init.nts <- init.nts.seq[ii]*(chi.type!=0) + 1/init.nts.seq[ii]*(vep.type!=0)
for(ss in 1:dim(chi.res.array)[2]){
	sigma <- sigma.seq[ss]
for(pp in 1:dim(chi.res.array)[3]){
	prop <- prop.seq[pp]

	q.res <- matrix(0, nrow=sim, ncol=5); chi.res <- vep.res <- array(0, dim=c(sim, 5, 3))
	stn <- est.q <- rep(0, sim)
	for(i in 1:sim){
		sd <- power.sim.data(T=T, n=n, r=r, chp=chp,
			ar=ar, alpha=alpha, chi.type=chi.type, sigma=sigma, chi.prop=prop, 
			cor=cor, H=H, vep.type=vep.type, vep.prop=prop, rho=rho, beta=beta, phi=phi*init.nts)
		x <- sd$x
# matplot(t(x), type="l", col=1); matlines(t(sd$chi), col=2); matlines(t(sd$vep), col=3)

# par(mfrow=c(1, 2), mar=c(2.5,2.5,.5,.5)); image.plot(sd$vep[,1:chp]%*%t(sd$vep[,1:chp])/chp); image.plot(sd$vep[,-(1:chp)]%*%t(sd$vep[,-(1:chp)])/(T-chp))

		stn[i] <- sqrt(sum(diag(cov(t(sd$chi)))))/sqrt(sum(diag(cov(t(sd$vep)))))
		
		gfm <- get.factor.model(x, bn.op=bn.op, max.q = max.q, normalisation=TRUE)
		est.q[i] <- q.hat <- which.min(gfm$ic)-1
# plot(gfm$ic)
		q.seq <- round(seq(max(1, q.hat), max(2*max(1, q.hat), max.q), length.out=5))
q.seq		
		cd <- chi.detect(gfm, q.seq=q.seq, scales=scales, B=Bsim, M=M, dw=dw)
		cd <- array(unlist(cd), dim=c(3, 6, length(q.seq))); cd <- cd[, -3, , drop=FALSE]
		ind <- q.seq0 <- rep(0, 5)
		for(op in 1:5){
			tmp <- which(cd[3, op, ] >= thr)
			if(length(tmp)>0) ind[op] <- max(tmp) else ind[op] <- 1
			q.seq0[op] <- q.seq[ind[op]]
		}
		q.res[i, ] <- q.seq0
cd
		vd <- vep.detect(gfm, q.seq0, scales=scales, diag=vep.diag, B=Bsim, dw=dw)
vd
		for(op in 1:5){
			chi.res[i, op, ] <- cd[, op, ind[op]]
			vep.res[i, op, ] <- vd$res.list[[op]]
		}
	}
	res.list <- list(chi=chi.res, vep=vep.res, q.res=q.res, est.q=est.q, stn=stn)
	if(chi.type > 0) name <- paste("c:\\users\\mahrc\\dropbox\\bcf_sim\\new_chi_T",T,"n",n,"_type", chi.type, "sigma", ss, "prop", pp, "chp", chp, "stn", ii, "alpha", alpha, "beta", beta, "rho", rho, "_res.RData", sep="")
	if(vep.type > 0) name <- paste("c:\\users\\mahrc\\dropbox\\bcf_sim\\new_vep_T",T,"n",n,"_type", vep.type, "sigma", ss, "prop", pp, "chp", chp, "stn", ii, "alpha", alpha, "beta", beta, "rho", rho, "_res.RData", sep="")
	save(res.list, file=name)
	
	chi.res.array[ii, ss, pp, ] <- apply(chi.res[,,3] >= thr, 2, mean)
	vep.res.array[ii, ss, pp, ] <- apply(vep.res[,,3] >= thr, 2, mean)

}}}

save(chi.res.array, file=paste("c:\\users\\mahrc\\dropbox\\", chi.type, "_chi.RData", sep=""))
save(vep.res.array, file=paste("c:\\users\\mahrc\\dropbox\\", chi.type, "_vep.RData", sep=""))

## chi.type <- 1
for(ss in 1:dim(chi.res.array)[2]){
#	pdf(paste("/users/h/dropbox/bcf/notes/sim_fig/n", n, "T", T, "_chi", chi.type, "_sigma", ss, ".pdf", sep=""), width=8, height=6)
	x11()
	par(mfcol=c(2, dim(chi.res.array)[3]), mar=c(2.5,2.5,1,.5))
	for(pp in 1:dim(chi.res.array)[3]){
		matplot(init.nts.seq, chi.res.array[, ss, pp, -3], type="l", ylim=c(0, 1), main=prop.seq[pp], lwd=2); abline(h=1-thr, col=8, lwd=2)
		if(pp==dim(chi.res.array)[3]) legend("topright", col=1:5, bty="n", lty=1:5, legend=c("0", "1/2", "0+1/2", "max", "avg")[-3])
		matplot(init.nts.seq, vep.res.array[, ss, pp, -3], type="l", ylim=c(0, 1), lwd=2); abline(h=1-thr, col=8, lwd=2)
#		name <- paste("/users/h/dropbox/bcf_sim/chi_T",T,"n",n,"_type", 1, "sigma", ss, "prop", pp, "chp", chp, "stn", ii, "alpha", alpha, "beta", beta, "rho", rho, "_res.RData", sep=""); load(name)
#		boxplot(q.res[, -3], xlab="", ylab="", xaxt="n"); axis(1, at=1:4, labels=c("0", "1/2", "0+1/2", "max", "avg")[-3])
	}
#	dev.off()
}

## chi.type <- 2
#pdf(paste("/users/h/dropbox/bcf/notes/sim_fig/n", n, "T", T, "_chi", 2, ".pdf", sep=""), width=6, height=2)
par(mfcol=c(1, 3), mar=c(2.5,2.5,1,.5))
matplot(init.nts.seq, chi.res.array[, 1, 1, -3], type="l", ylim=c(0, 1), lwd=2); abline(h=1-thr, col=8, lwd=2)
legend("topright", col=1:5, bty="n", lty=1:5, legend=c("0", "1/2", "0+1/2", "max", "avg")[-3], cex=1)
matplot(init.nts.seq, vep.res.array[, 1, 1, -3], type="l", ylim=c(0, 1), lwd=2); abline(h=1-thr, col=8, lwd=2)
#name <- paste("/users/h/dropbox/bcf_sim/chi_T",T,"n",n,"_type", 2, "sigma", 1, "prop", 1, "chp", chp, "stn", 1, "alpha", alpha, "beta", beta, "rho", rho, "_res.RData", sep=""); load(name)
#boxplot(res.list$q.res[, -3], xlab="", ylab="", xaxt="n"); axis(1, at=1:4, labels=c("0", "1/2", "0+1/2", "max", "avg")[-3])
#dev.off()

## chi.type <- 3
#pdf(paste("/users/h/dropbox/bcf/notes/sim_fig/n", n, "T", T, "_chi", 3, ".pdf", sep=""), width=8, height=6)
par(mfcol=c(2, dim(chi.res.array)[3]), mar=c(2.5,2.5,1,.5))
for(pp in 1:dim(chi.res.array)[3]){
	matplot(init.nts.seq, chi.res.array[, 1, pp, -3], type="l", ylim=c(0, 1), lwd=2, main=prop.seq[pp]); abline(h=1-thr, col=8, lwd=2)
	if(pp==dim(chi.res.array)[3]) legend("topleft", col=1:5, bty="n", lty=1:5, legend=c("0", "1/2", "0+1/2", "max", "avg"), cex=.5)
	matplot(init.nts.seq, vep.res.array[, 1, pp, -3], type="l", ylim=c(0, 1), lwd=2, main=prop.seq[pp]); abline(h=1-thr, col=8, lwd=2)
#	name <- paste("/users/h/dropbox/bcf_sim/chi_T",T,"n",n,"_type", 3, "sigma", 1, "prop", pp, "chp", chp, "stn", 1, "alpha", alpha, "beta", beta, "rho", rho, "_res.RData", sep=""); load(name)
#	boxplot(res.list$q.res[, -3], xlab="", ylab="", xaxt="n"); axis(1, at=1:4, labels=c("0", "1/2", "0+1/2", "max", "avg")[-3])
}
dev.off()

## vep.type <- 1
#pdf(paste("/users/h/dropbox/bcf/notes/sim_fig/n", n, "T", T, "_vep", 1, ".pdf", sep=""), width=8, height=6)
par(mfcol=c(3, dim(chi.res.array)[3]), mar=c(2.5,2.5,1,.5))
for(pp in 1:dim(chi.res.array)[3]){
	matplot(1/init.nts.seq, chi.res.array[, 1, pp,], type="l", ylim=c(0, 1), main=prop.seq[pp], lwd=2); abline(h=1-thr, col=8, lwd=2)
	if(pp==dim(chi.res.array)[3]) legend("topright", col=1:5, bty="n", lty=1:5, legend=c("0", "1/2", "0+1/2", "max", "avg")[-3])
	matplot(1/init.nts.seq, vep.res.array[, 1, pp,], type="l", ylim=c(0, 1), lwd=2); abline(h=1-thr, col=8, lwd=2)
#	name <- paste("/users/h/dropbox/bcf_sim/vep_T",T,"n",n,"_type", 1, "sigma", 1, "prop", pp, "chp", chp, "stn", 1, "alpha", alpha, "beta", beta, "rho", rho, "_res.RData", sep=""); load(name)
#	boxplot(res.list$q.res[, -3], xlab="", ylab="", xaxt="n"); axis(1, at=1:4, labels=c("0", "1/2", "0+1/2", "max", "avg")[-3])	
}
#dev.off()

apply(chi.res[,,3] > thr, 2, mean)
apply(vep.res[,,3] > thr, 2, mean)
apply(vep.res[,,6] > thr, 2, mean)

plot(nts); abline(h=phi, col=2)
plot(est.q); abline(h=r, col=2)

k <-5; hist(chi.res[chi.res[, k, 3] >= thr, k, 2])
table(chi.res[chi.res[, k, 3] > thr, k, 2])

k <- 2; hist(vep.res[vep.res[, k, 3] > thr | vep.res[, k, 6] > thr, k, 2])

