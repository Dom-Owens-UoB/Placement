memory.limit(size=8142*10)

library(Rcpp); library(RcppArmadillo)
library(factorcptt)

library(doParallel)
cl <- makeCluster(4); registerDoParallel(cl)
# stopCluster(cl)

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

for(ii in 2:4){
for(ss in 1:4){
for(pp in 1:4){

if(ii==2 & ss==1) next
if(ii==2 & ss==2 & pp <=2) next

	name <- paste("c:\\users\\mahrc\\dropbox\\bcf_sim\\multi_n", n, "T", T, "_nts", ii, "sigma", ss, "prop", pp, "_chi.chp", length(chi.cps), "_vep.cps", length(vep.cps), "_r", r, ".RData", sep="")
	load(name)
	
	dcbs.res <- array(0, dim=c(sim, T))
	for(i in 1:sim){
		sd <- power.sim.data(T=T, n=n, r=r, chp=chp,
			ar=ar, alpha=alpha, chi.type=chi.type, sigma=sigma.seq[ss], chi.prop=rep(prop.seq[pp], length(chp)), 
			cor=cor, H=H, vep.type=vep.type, vep.prop=rep(prop.seq[pp], length(chp)), rho=rho, beta=beta, phi=phi*init.nts.seq[ii])
		x <- sd$x
# matplot(t(x), type="l", col=1); matlines(t(sd$chi), col=2); matlines(t(sd$vep), col=3)
    	
    	da <- dcbs.alg(x, sig.lev=sig.lev, scales=scales, dw=round(dw/2), B=B, rule=rule)
    	dcbs.res[i, da$est.cps] <- 1
      	   	   	
	}
	ls$dr <- dcbs.res	
	save(ls, file=name)
}}
}
