test_that("Codes work for Zero-inflated model", {
  library(iZID)
  X<-matrix(rnorm(400), ncol=4,nrow=100)
  Y<-sample.zi(N=100,phi=0.2,distri='nb',r=10,p=0.6)
  hsz<-horseshoe_zinb_mc(X=X,y=Y,mc_length=100,burnin=60,thin=1)

  dlz<-dirichlet_laplace_zinb_mc(X=X,y=Y,mc_length=100,a=0.5,burnin=60,thin=1)

  dpz<-double_pareto_zinb_mc(y=Y, X=X, mc_length=100,burnin=60,thin=1)


  postbeta.h<-hsz$Beta
  postmeanbeta.h<- apply(postbeta.h,2,mean)
  cibeta.h<-t(apply(postbeta.h,2,quantile,probs = c(0.025,0.5,0.975)))
  # if zero included, not significant
  postalpha.h<-hsz$Alpha
  postmeanalpha.h<- apply(postalpha.h,2,mean)
  cialpha.h<-t(apply(postalpha.h,2,quantile,probs = c(0.025,0.5,0.975)))
  cibeta.h
  cialpha.h
  postbeta.dl<-dlz$Beta
  postmeanbeta.dl<- apply(postbeta.dl,2,mean)
  cibeta.dl<-t(apply(postbeta.dl,2,quantile,probs = c(0.025,0.5,0.975)))
  # if zero included, not significant
  postalpha.dl<-dlz$Alpha
  postmeanalpha.dl<- apply(postalpha.dl,2,mean)
  cialpha.dl<-t(apply(postalpha.dl,2,quantile,probs = c(0.025,0.5,0.975)))
   postbeta.dp<-dpz$Beta
  postmeanbeta.dp<- apply(postbeta.dp,2,mean)
  cibeta.dp<-t(apply(postbeta.dp,2,quantile,probs = c(0.025,0.5,0.975)))
  cibeta.dl
  cialpha.dl
  postalpha.dp<-dpz$Alpha
  postmeanalpha.dp<- apply(postalpha.dp,2,mean)
  cialpha.dp<-t(apply(postalpha.dp,2,quantile,probs = c(0.025,0.5,0.975)))
  cibeta.dp
  cialpha.dp

  })
