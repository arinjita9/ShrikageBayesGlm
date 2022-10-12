test_that("Codes work for Negative Binomial Regression model", {
  X<-matrix(rnorm(400), ncol=4,nrow=100)
  Y<-rnbinom(100,1,0.5)
  burnin<-60
  mc_length=100
  beta_start<-matrix(rep(0.5,4), ncol=1, nrow=4)
  hsn<-horseshoe_negative_binomial_mc(Y=Y, X=X,h=4, mc_length=mc_length, starting_beta=beta_start)
  dln<-dirichlet_laplace_negative_binomial_mc(Y=Y, X=X,h=4, mc_length=mc_length,a=0.7, starting_beta=beta_start)
  dpn<- double_pareto_negative_binomial_mc(Y=Y, X=X,h=4, mc_length=mc_length, starting_beta=beta_start)
  #########Burnin Step##################
  postbeta.h<-hsn[[1]][(burnin+1):mc_length,]
  postbeta.dl<-dln[[1]][(burnin+1):mc_length,]
  postbeta.dp<-dpn[[1]][(burnin+1):mc_length,]
  ci.h<-t(apply(postbeta.h,2,quantile,probs = c(0.025,0.5,0.975)))
  ci.dl<-t(apply(postbeta.dl,2,quantile,probs = c(0.025,0.5,0.975)))
  ci.dp<-t(apply(postbeta.dp,2,quantile,probs = c(0.025,0.5,0.975)))
  ci.h
  ci.dl
  ci.dp
  })
