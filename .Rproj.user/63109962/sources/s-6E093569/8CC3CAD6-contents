
mosumvar_factor_sim_loadings <- function(N_sim=100, r =6, noise = 1, proportion = 1,  d=150, r_known = T){
  if(r_known) r_ <- r else r_ <- NULL
  #r <- max(r1,r)
  out <- matrix(NA, nrow = N_sim, ncol = 3)
  true_cps <- c(100,200)
  A_sim1_r <- A_sim1[1:r,1:r] 
  A_sim2_r <- A_sim2[1:r,1:r] 
  A_sim3_r <- A_sim3[1:r,1:r] 
  loadings1 <- matrix(rnorm(d*r, 0,1) , nrow = r, ncol = d)
 
  for (sim in 1:N_sim) {
    loadings2 <- loadings1
    factors_ <-  rbind(mosumvar::VAR.sim(100, coeffs = A_sim1_r),mosumvar::VAR.sim(100, coeffs = A_sim2_r),
                       mosumvar::VAR.sim(100, coeffs = A_sim3_r))  
    random_index <- sample(1:(d*r), floor(d*r*proportion))
    loadings2[random_index] <- loadings1[random_index] + matrix(rnorm(d*r, 0, noise) , nrow = r, ncol = d)
    data_1 <- factors_[1:100,] %*% loadings1
    data_2 <- factors_[101:200,] %*% loadings2
    data_3 <- factors_[201:300,] %*% loadings1
    data <- rbind(data_1,data_2,data_3)
    #nouniv <- mosumvar_factor(data,p=1,G=50,r=r_, method = "Score")
    univ <- mosumvar_factor(data,p=1,G=50,r=r_, univ = T, method = "Score")
    #    fcpt <- factorcpt::factor.seg.alg(t(data_), r= r_)
    #if(!is.null(nouniv$cps)) out[sim,1] <- pracma::hausdorff_dist(nouniv$cps, true_cps) / 300 else out[sim,1] <- 1
    if(!is.null(univ$cps)) out[sim,2] <- pracma::hausdorff_dist(univ$cps, true_cps) / 300 else out[sim,2] <-1
    #   out[sim,3] <- pracma::hausdorff_dist(fcpt$common.est.cps, true_cps) / 300
  }
  return(out)
}
mosumvar_factor_sim_loadings(5, 3, .1, .25, 150, F)
