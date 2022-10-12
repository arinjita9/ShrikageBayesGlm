test_that("Codes work for Logistic Regression model", {
  X<-matrix(rnorm(400), ncol=4,nrow=100)
  Y<-rbinom(100,1,0.5)
  mc_length=100
  burnin<-60
  beta_start<-matrix(rep(0,4), ncol=1, nrow=4)
  hsl<-horseshoe_logistic_mc(Y=Y, X=X, mc_length=mc_length,beta_start=beta_start )
  dll<-dirichlet_laplace_logistic_mc(Y=Y, X=X, mc_length=mc_length,beta_start=beta_start, a=0.7)
  dpl<-double_pareto_logistic_mc(Y=Y, X=X, mc_length=mc_length,beta_start=beta_start)
  postbeta.h<-hsl[[1]][(burnin+1):mc_length,]
  postbeta.dl<-dll[[1]][(burnin+1):mc_length,]
  postbeta.dp<-dpl[[1]][(burnin+1):mc_length,]
  ci.h<-t(apply(postbeta.h,2,quantile,probs = c(0.025,0.5,0.975)))
  ci.dl<-t(apply(postbeta.dl,2,quantile,probs = c(0.025,0.5,0.975)))
  ci.dp<-t(apply(postbeta.dp,2,quantile,probs = c(0.025,0.5,0.975)))
   ci.h
   ci.dl
   ci.dp
  })
