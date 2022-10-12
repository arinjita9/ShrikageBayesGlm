test_that("Codes work for Multinomial Logistic Regression model", {
  N =400; # sample size
  P = 4;
  J = 3;
  num = 1;

  beta = matrix(c(-1.3,1.2,0,0,2,-1.5,0,0), P, J-1)
  m=c(-2,-1,0,1,2)
  sd=1
  burnin=60
  X    = cbind(1, matrix(rnorm(N*(P-1),mean=m, sd=sd), # intercept included, total four covariates #each have different means
                           N, P-1));
    psi  = X %*% beta #N*p matrix
    free.p = exp(psi) / (1 + rowSums(exp(psi)))#N*p matrix
    last.p = 1 - rowSums(free.p) #N*1
    p      = cbind(free.p, last.p) #N*J
    ## Only one roll per ith sample.
    y.cat = rep(0, N)
    for (i in 1:N)
      y.cat[i] = sample.int(J, 1, prob=p[i,]) #length N # no. of J items to choose from
    y.all = model.matrix( ~ factor(y.cat) - 1 ) # N*J
    y = y.all[,-J] # removing the J th column
    n = rep(1, N)
   # a.df = data.frame("response"=y.cat, "predictors"=X[,-1])
  hsml<-horseshoe_mlogit_mc(y=y, X=X,samp=100)
  dpml<-double_pareto_mlogit_mc(y=y, X=X,samp=100)
  dlml<-dirichlet_laplace_mlogit_mc(y=y, X=X,samp=100)
  mc_length=100
  postbeta.h<- hsml$finalbeta[((burnin+1):mc_length),,]
  beta.mat.h = array(postbeta.h, dim=c(dim(postbeta.h)[1], P*(J-1)));
  beta.B.h = apply(postbeta.h, c(2,3), mean)
  ci.h<-t(matrix(apply(postbeta.h,c(2,3),quantile,probs = c(0.025,0.5,0.975)),nrow=3,ncol=dim(beta.B.h)[1]*(J-1)))
  colnames(ci.h)=c("2.5%","50%","97.5%")
  ci.h
  postbeta.dl<-  dlml$finalbeta[((burnin+1):mc_length),,]
  beta.mat.dl = array(postbeta.dl, dim=c(dim(postbeta.dl)[1], P*(J-1)));
  beta.B.dl = apply(postbeta.dl, c(2,3), mean)
  ci.dl<-t(matrix(apply(postbeta.dl,c(2,3),quantile,probs = c(0.025,0.5,0.975)),nrow=3,ncol=dim(beta.B.dl)[1]*(J-1)))
  colnames(ci.dl)=c("2.5%","50%","97.5%")
  ci.dl
  postbeta.dp<- dpml$finalbeta[((burnin+1):mc_length),,]
  beta.mat.dp = array(postbeta.dp, dim=c(dim(postbeta.dp)[1], P*(J-1)));
  beta.B.dp = apply(postbeta.dp, c(2,3), mean)
  ci.dp<-t(matrix(apply(postbeta.dp,c(2,3),quantile,probs = c(0.025,0.5,0.975)),nrow=3,ncol=dim(beta.B.dp)[1]*(J-1)))
  colnames(ci.dp)=c("2.5%","50%","97.5%")
  ci.dp
})
