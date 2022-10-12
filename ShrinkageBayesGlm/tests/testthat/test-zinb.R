#' Zero-inflated Negative Binomial Regression with Horseshoe prior
#'
#' @param X design matrix
#' @param y zero-inflated count response
#' @param mc_length no. of MCMC iterations
#' @param burnin burnin
#' @param thin thinning parameter, generally set to 1
#'
#' @export list of posterior samples of the parameters
#' @return Alpha: posterior samples of Alpha parameter
#' @return Beta: posterior samples of Beta parameter
#' @return R or R2: posterior samples of R parameter
#' @return xi: posterior samples of hyper-parameter of the prior for Beta
#' @return gamma:posterior samples of hyper-parameter of the prior for Beta
#' @return lambda:posterior samples of hyper-parameter of the prior for Beta
#' @return tausq: posterior samples of hyper-parameter of the prior for Beta
#' @return xi1:posterior samples of hyper-parameter of the prior for Alpha
#' @return gamma1:posterior samples of hyper-parameter of the prior for Alpha
#' @return lambda1:posterior samples of hyper-parameter of the prior for Alpha
#' @return tausq1:posterior samples of hyper-parameter of the prior for Alpha
#'
#' @examples
horseshoe_zinb_mc<-function(X=X,y=Y,mc_length=mc_length,burnin=burnin,thin=thin){
  # Priors -- note diffuse uniform for r
  n<-nrow(X)
  p<-ncol(X)
  beta0<-alpha0<-rep(0,p)
  # T0a<-diag(0.01,p)          # Prior precision for alpha -- tightening can help in bordeline cases
  # T0b<-diag(0.01,p)          # Prior precision for beta
  xi<-1     ## xi
  tau_sq <-1 ## tau^2
  gamma_j=as.matrix(replicate(n = p,expr = 1)) # gamma_j
  lambdaj_sq<-as.matrix( replicate(n = p,expr = 1)  )
  sigma <-10*diag(as.numeric(tau_sq%*%t(lambdaj_sq)),p,p) # sigma
  xi1<-1     ## xi
  tau_sq1 <-1 ## tau^2
  gamma_j1=as.matrix(replicate(n = p,expr = 1)) # gamma_j
  lambdaj_sq1<-as.matrix( replicate(n = p,expr = 1)  )
  sigma1 <-10*diag(as.numeric(tau_sq1%*%t(lambdaj_sq1)),p,p) # sigma
  # Inits
  beta<-alpha<-rep(0,p)

  r<-1
  Acc<-0
  y1<-rbinom(n,1,.5)        # At risk indicator (W in paper)
  y1[y>0]<-1                # If y>0, then patient is at risk w.p. 1
  n0<-length(y[y==0])       # Number of observed 0's
  q<-rep(.5,n)              # 1-p=1/(1+exp(X%*%alpha)), used for updating y1
  s<-0.025

  ############
  # Num Sims #
  ############
  # nsim<-1500			          # Number of MCMC Iterations (Used 50500 in paper)
  # thin<-1			              # Thinning interval
  # burn<-500    		          # burn
  lastit<-(mc_length-burnin)/thin	# Last stored value

  #########
  # Store #
  #########
  Beta<-Alpha<-matrix(0,lastit,p)
  R<-rep(0,lastit)
  R2<-rep(0,lastit)         # To store Gibbs updates
  xi_all= array(NA, dim=c(lastit,1))
  gamma_all=array(NA, dim=c(lastit,p))
  lambda_all=array(NA, dim=c(lastit,p))
  tau_sq_all=array(NA, dim=c(lastit,1))
  xi_all1= array(NA, dim=c(lastit,1))
  gamma_all1=array(NA, dim=c(lastit,p))
  lambda_all1=array(NA, dim=c(lastit,p))
  tau_sq_all1=array(NA, dim=c(lastit,1))
  ########
  # MCMC #
  ########
  #tmp<-proc.time()
  # k1=1
  # k=1
  for (i in 1:mc_length){

    # Update alpha
    mu<-X%*%alpha
    w<-pgdraw(1,mu)#rpg(n,1,mu)
    z<-y1-1/2                               # Or z=(y-1/2)/w with m=v%*%(T0%*%beta0+t(sqrt(w)*X)%*%sqrt(w)*z)
    sigma1 <-diag(as.numeric(tau_sq1%*%t(lambdaj_sq1)),p,p)
    v<-solve(crossprod(X*sqrt(w))+sigma1)      # Or solve(X%*%W%*%X), which is slower
    m<-v%*%(sigma1%*%alpha0+t(X)%*%z)
    alpha<-c(rmvnorm(1,m,v))
    for(k1 in 1:p){

      gamma_j1[k1] <-rigamma(1,1,(1+(1/lambdaj_sq1[k1])))
      lambdaj_sq_post_rate1=((1/gamma_j[k1])+(0.5*(alpha[k1])^2/tau_sq1))
      lambdaj_sq1[k1]<-rigamma(1,1,lambdaj_sq_post_rate1)

    }
    xi1 <-rigamma(n = 1,alpha = 1,beta = (1+1/tau_sq1))

    tau_sq_post_rate1=(1/xi)+.5*sum((alpha^2)/lambdaj_sq1)
    tau_sq1 <-rigamma(1,(p+1)/2,tau_sq_post_rate1 )
    # Update latent class indicator y1 (w in paper)
    eta<-X%*%alpha
    phi<-exp(eta)/(1+exp(eta))              # at-risk probability
    theta<-phi*(q^r)/(phi*(q^r)+1-phi)      # Conditional prob that y1=1 given y=0 -- i.e. Pr(chance zero|observed zero)
    y1[y==0]<-rbinom(n0,1,theta[y==0])      # If y=0, then draw a "chance zero" w.p. theta, if y=1, then y1=1
    n1<-sum(y1)

    # Update rho                              # Proposal variance
    rnew<-rtnorm(1,r,sqrt(s),lower=0)       # Treat r as continuous
    ratio<-sum(dnbinom(y[y1==1],rnew,q[y1==1],log=T))-sum(dnbinom(y[y1==1],r,q[y1==1],log=T))+
      dtnorm(r,rnew,sqrt(s),0,log=T) - dtnorm(rnew,r,sqrt(s),0,log=T) # Note: Uniform prior for r
    #       Proposal not symmetric
    if (log(runif(1))<ratio) {
      r<-rnew
      Acc<-Acc+1
    }

    # Update beta
    eta1<-X[y1==1,]%*%beta
    pos.int<-ifelse(y[y1==1]+r<1,1,as.integer(y[y1==1]+r))
    w<-pgdraw(pos.int,eta1)                             # Polya weights
    z<-(y[y1==1]-r)/(2*w)
    sigma <-diag(as.numeric(tau_sq%*%t(lambdaj_sq)),p,p)# Latent "response"
    v<-solve(crossprod(X[y1==1,]*sqrt(w))+sigma)
    m<-v%*%(sigma%*%beta0+t(sqrt(w)*X[y1==1,])%*%(sqrt(w)*z))
    beta<-c(rmvnorm(1,m,v))
    eta<-X%*%beta
    q<-1/(1+exp(eta))
    xi <-rigamma(n = 1,alpha = 1,beta = (1+1/tau_sq))
    tau_sq_post_rate=(1/xi)+.5*sum((beta^2)/lambdaj_sq)
    tau_sq <-rigamma(1,(p+1)/2,tau_sq_post_rate )
    for(k in 1:p){
      gamma_j[k] <-rigamma(1,1,(1+(1/lambdaj_sq[k])))
      lambdaj_sq_post_rate=((1/gamma_j[k])+(0.5*(beta[k])^2/tau_sq))
      lambdaj_sq[k]<-rigamma(1,1,lambdaj_sq_post_rate)
    }
    # Update r2 using Gibbs as in Dadaneh et al and Zhou and Carin #

    # Update latent counts, l
    l<-rep(0,n1)
    ytmp<-y[y1==1]
    for(j in 1:n1) l[j]<-sum(rbinom(ytmp[j],1,round(r/(r+1:ytmp[j]-1),6))) # Could try to avoid loop; rounding avoids numerical stability

    # Update r from conjugate gamma distribution given l and psi
    eta<-X[y1==1,]%*%beta
    psi<-exp(eta)/(1+exp(eta))
    r2<-rgamma(1,0.01+sum(l),0.01-sum(log(1-psi))) # Gamma(0.01,0.01) prior for r



    # Store
    if (i> burnin & i%%thin==0) {
      j<-(i-burnin)/thin
      Alpha[j,]<-alpha
      Beta[j,]<-beta
      R[j]<-r
      R2[j]<-r2
      xi_all[j,]= c(xi)
      gamma_all[j,]<-c(gamma_j)
      lambda_all[j,]<-c(lambdaj_sq)
      tau_sq_all[j,]<-c(tau_sq)
      xi_all1[j,]= c(xi1)
      gamma_all1[j,]<-c(gamma_j1)
      lambda_all1[j,]<-c(lambdaj_sq1)
      tau_sq_all1[j,]<-c(tau_sq1)
    }

    if (i%%100==0) print(i)

  }
  MC_Hzinb<-list(Alpha=Alpha,Beta=Beta, R=R,R2=R2,xi=xi_all,gamma=gamma_all,
                 lambda=lambda_all, tausq=tau_sq_all,xi1=xi_all1,gamma1=gamma_all1,lambda1=lambda_all1,
                 tausq1=tau_sq_all1)
}

#tot.time<-proc.time()-tmp  # 9.8 minutes to run 50K
###################Dirichlet Laplace#################
#' Zero-inflated Negative Binomial Regression with Dirichlet Laplace prior
#'
#' @param X design matrix
#' @param y zero-inflated count response
#' @param mc_length no. of MCMC iterations
#' @param burnin burnin
#' @param thin thinning parameter, generally set to 1
#'
#' @export list of posterior samples of the parameters
#' @return Alpha: posterior samples of Alpha parameter
#' @return Beta: posterior samples of Beta parameter
#' @return R or R2: posterior samples of R parameter
#' @return shi_all: posterior samples of hyper-parameter of the prior for Beta
#' @return fi_all:posterior samples of hyper-parameter of the prior for Beta
#' @return tau_sq_all:posterior samples of hyper-parameter of the prior for Beta
#' @return shi_all1: posterior samples of hyper-parameter of the prior for Alpha
#' @return fi_all1:posterior samples of hyper-parameter of the prior for Alpha
#' @return tau_sq_all1:posterior samples of hyper-parameter of the prior for Alpha
#'
#' @examples

dirichlet_laplace_zinb_mc<-function(X=X,y=Y,mc_length=mc_length,a=a,burnin=burnin,thin=thin){
  # Priors -- note diffuse uniform for r
  n<-nrow(X)
  p<-ncol(X)
  beta0<-alpha0<-rep(0,p)
  # T0a<-diag(0.01,p)          # Prior precision for alpha -- tightening can help in bordeline cases
  # T0b<-diag(0.01,p)          # Prior precision for beta
  tau_sq <-1 ## tau^2
  shi=(replicate(n = p,expr = 1)) # gamma_j
  fi<-( replicate(n = p,expr = 1)  )/p
  sigma <-10*diag(as.numeric(shi*t(fi^2)*tau_sq),p,p)# sigma ##
  tau_sq1 <-1 ## tau^2
  shi1=(replicate(n = p,expr = 1)) # gamma_j
  fi1<-( replicate(n = p,expr = 1)  )/p
  sigma1 <-10*diag(as.numeric(shi1*t(fi1^2)*tau_sq1),p,p)# sigma ##
  # Inits
  beta<-alpha<-rep(0,p)

  r<-1
  Acc<-0
  y1<-rbinom(n,1,.5)        # At risk indicator (W in paper)
  y1[y>0]<-1                # If y>0, then patient is at risk w.p. 1
  n0<-length(y[y==0])       # Number of observed 0's
  q<-rep(.5,n)              # 1-p=1/(1+exp(X%*%alpha)), used for updating y1
  s<-0.025

  ############
  # Num Sims #
  ############
  # nsim<-1500			          # Number of MCMC Iterations (Used 50500 in paper)
  # thin<-1			              # Thinning interval
  # burn<-500    		          # burn
  lastit<-(mc_length-burnin)/thin	# Last stored value

  #########
  # Store #
  #########
  Beta<-Alpha<-matrix(0,lastit,p)
  R<-rep(0,lastit)
  R2<-rep(0,lastit)         # To store Gibbs updates
  shi_all= array(NA, dim=c(lastit,p))
  fi_all=array(NA, dim=c(lastit,p))
  tau_sq_all=array(NA, dim=c(lastit,1))
  shi_all1= array(NA, dim=c(lastit,p))
  fi_all1=array(NA, dim=c(lastit,p))
  tau_sq_all1=array(NA, dim=c(lastit,1))
  ########
  # MCMC #
  ########
  #tmp<-proc.time()

  for (i in 1:mc_length){

    # Update alpha
    mu<-X%*%alpha
    w<-pgdraw(1,mu)#rpg(n,1,mu)
    z<-y1-1/2                               # Or z=(y-1/2)/w with m=v%*%(T0%*%beta0+t(sqrt(w)*X)%*%sqrt(w)*z)
    sigma1 <-diag(as.numeric(shi1*t(fi1^2)*tau_sq1),p,p)
    v<-solve(crossprod(X*sqrt(w))+sigma1)      # Or solve(X%*%W%*%X), which is slower
    m<-v%*%(sigma1%*%alpha0+t(X)%*%z)
    alpha<-c(rmvnorm(1,m,v))
    #tau--bhatta
    tau1= rgig(n=1,lambda=p*a-p, chi= 2*sum(abs(alpha)/fi1) ,psi=1  )
    tau_sq1= tau1^2

    #phi --bhatta
    T_fi1=rep(NA,p)
    for(k1 in 1:p){

      T_fi1[k1] <- rgig(n=1,lambda = a-1,chi = 2*abs(alpha[k1]),psi = 1)
    }

    fi1=T_fi1/sum(T_fi1)


    for(k1 in 1:p){

      shi1[k1] <-rgig(n=1,lambda = 1/2,chi=(alpha[k1])^2/(tau_sq1*(fi1[k1])^2),psi =1)

    }

    cand1=shi1*(fi1^2)*tau_sq1
    sigma1 <-diag(cand1,p,p)
    sigma_inv1<-diag(1/(cand1),p,p)
    # Update latent class indicator y1 (w in paper)
    eta<-X%*%alpha
    phi<-exp(eta)/(1+exp(eta))              # at-risk probability
    theta<-phi*(q^r)/(phi*(q^r)+1-phi)      # Conditional prob that y1=1 given y=0 -- i.e. Pr(chance zero|observed zero)
    y1[y==0]<-rbinom(n0,1,theta[y==0])      # If y=0, then draw a "chance zero" w.p. theta, if y=1, then y1=1
    n1<-sum(y1)

    # Update rho                              # Proposal variance
    rnew<-rtnorm(1,r,sqrt(s),lower=0)       # Treat r as continuous
    ratio<-sum(dnbinom(y[y1==1],rnew,q[y1==1],log=T))-sum(dnbinom(y[y1==1],r,q[y1==1],log=T))+
      dtnorm(r,rnew,sqrt(s),0,log=T) - dtnorm(rnew,r,sqrt(s),0,log=T) # Note: Uniform prior for r
    #       Proposal not symmetric
    if (log(runif(1))<ratio) {
      r<-rnew
      Acc<-Acc+1
    }

    # Update beta
    eta1<-X[y1==1,]%*%beta
    pos.int<-ifelse(y[y1==1]+r<1,1,as.integer(y[y1==1]+r))
    w<-pgdraw(pos.int,eta1)    #rpg(n1,y[y1==1]+r,eta1)                              # Polya weights
    z<-(y[y1==1]-r)/(2*w)                                   # Latent "response"
    v<-solve(crossprod(X[y1==1,]*sqrt(w))+sigma)
    m<-v%*%(sigma%*%beta0+t(sqrt(w)*X[y1==1,])%*%(sqrt(w)*z))
    beta<-c(rmvnorm(1,m,v))
    eta<-X%*%beta
    q<-1/(1+exp(eta))

    #tau--bhatta
    tau= rgig(n=1,lambda=p*a-p, chi= 2*sum(abs(beta)/fi) ,psi=1  )
    tau_sq= tau^2

    #fi --bhatta
    T_fi=rep(NA,p)
    for(k in 1:p){
      T_fi[k] <- rgig(n=1,lambda = a-1,chi = 2*abs(beta[k]),psi = 1)
    }

    fi=T_fi/sum(T_fi)


    for(k in 1:p){
      shi[k] <-rgig(n=1,lambda = 1/2,chi=(beta[k])^2/(tau_sq*(fi[k])^2),psi =1)

    }


    cand=shi*(fi^2)*tau_sq
    sigma <-diag(cand,p,p)
    sigma_inv<-diag(1/(cand),p,p)
    # Update r2 using Gibbs as in Dadaneh et al and Zhou and Carin #

    # Update latent counts, l
    l<-rep(0,n1)
    ytmp<-y[y1==1]
    for(j in 1:n1) l[j]<-sum(rbinom(ytmp[j],1,round(r/(r+1:ytmp[j]-1),6))) # Could try to avoid loop; rounding avoids numerical stability

    # Update r from conjugate gamma distribution given l and psi
    eta<-X[y1==1,]%*%beta
    psi<-exp(eta)/(1+exp(eta))
    r2<-rgamma(1,0.01+sum(l),0.01-sum(log(1-psi))) # Gamma(0.01,0.01) prior for r

    # Store
    if (i> burnin & i%%thin==0) {
      j<-(i-burnin)/thin
      Alpha[j,]<-alpha
      Beta[j,]<-beta
      R[j]<-r
      R2[j]<-r2
      shi_all[j,]= c(shi)
      fi_all[j,]=c(fi)
      tau_sq_all[j,]=c(tau_sq)
      shi_all1[j,]= c(shi1)
      fi_all1[j,]=c(fi1)
      tau_sq_all1[j,]=c(tau_sq1)
    }

    if (i%%100==0) print(i)

  }
  MC_DLzinb<-list(Alpha=Alpha,Beta=Beta, R=R,R2=R2,
                  shi_all= shi_all,
                  fi_all=fi_all,
                  tau_sq_all=tau_sq_all,
                  shi_all1= shi_all1,
                  fi_all1=fi_all1,
                  tau_sq_all1=tau_sq_all1)
}
#tot.time<-proc.time()-tmp  # 9.8 minutes to run 50K
###################Double Pareto#################
#' Zero-inflated Negative Binomial Regression with Double Pareto prior
#'
#' @param X design matrix
#' @param y zero-inflated count response
#' @param mc_length no. of MCMC iterations
#' @param burnin burnin
#' @param thin thinning parameter, generally set to 1
#' @export list of posterior samples of the parameters
#' @return Alpha: posterior samples of Alpha parameter
#' @return Beta: posterior samples of Beta parameter
#' @return R or R2: posterior samples of R parameter
#' @return lambda_sq_all: posterior samples of hyper-parameter of the prior for Beta
#' @return tau_sq_all:posterior samples of hyper-parameter of the prior for Beta
#' @return lambda_sq_all1:posterior samples of hyper-parameter of the prior for Alpha
#' @return tau_sq_all1: posterior samples of hyper-parameter of the prior for Alpha
#'
#' @examples
double_pareto_zinb_mc<-function(y=y, X=X, mc_length=mc_length,burnin=burnin,thin=1){
  # Priors -- note diffuse uniform for r
  n<-nrow(X)
  p<-ncol(X)
  beta0<-alpha0<-rep(0,p)
  sig_sq<-1
  neu<-1
  Eta<-1
  tauj_sq <-(replicate(n = p,expr = 1)) ## tau^2
  lambdaj_sq<-( replicate(n = p,expr = 1)  )
  sigma <-10*diag(as.numeric(sig_sq%*%t(tauj_sq)),p,p)  # sigma
  sig_sq1<-1
  neu1<-1

  Eta1<-1
  tauj_sq1 <-(replicate(n = p,expr = 1)) ## tau^2

  lambdaj_sq1<-( replicate(n = p,expr = 1)  )
  sigma1 <-10*diag(as.numeric(sig_sq1%*%t(tauj_sq1)),p,p)  # sigma
  # Inits
  beta<-alpha<-rep(0,p)

  r<-1
  Acc<-0
  y1<-rbinom(n,1,.5)        # At risk indicator (W in paper)
  y1[y>0]<-1                # If y>0, then patient is at risk w.p. 1
  n0<-length(y[y==0])       # Number of observed 0's
  q<-rep(.5,n)              # 1-p=1/(1+exp(X%*%alpha)), used for updating y1
  s<-0.025

  ############
  # Num Sims #
  ############
  lastit<-(mc_length-burnin)/thin	# Last stored value

  #########
  # Store #
  #########
  Beta<-Alpha<-matrix(0,lastit,p)
  R<-rep(0,lastit)
  R2<-rep(0,lastit)         # To store Gibbs updates
  #sig_sq_all = array(NA, dim=c(lastit,1))
  lambda_sq_all=array(NA, dim=c(lastit,p))
  tau_sq_all=array(NA, dim=c(lastit,p))
  #sig_sq_all1 = array(NA, dim=c(lastit,1))
  lambda_sq_all1=array(NA, dim=c(lastit,p))
  tau_sq_all1=array(NA, dim=c(lastit,p))
  ########
  # MCMC #
  ########
  #tmp<-proc.time()

  for (i in 1:mc_length){

    # Update alpha
    mu<-X%*%alpha
    w<-pgdraw(1,mu)#n
    z<-y1-1/2                               # Or z=(y-1/2)/w with m=v%*%(T0%*%beta0+t(sqrt(w)*X)%*%sqrt(w)*z)
    sigma1 <-diag(as.numeric(sig_sq1%*%t(tauj_sq1)),p,p)
    v<-solve(crossprod(X*sqrt(w))+sigma1)      # Or solve(X%*%W%*%X), which is slower
    m<-v%*%(sigma1%*%alpha0+t(X)%*%z)
    alpha<-c(rmvnorm(1,m,v))
    tauj1=lambdaj1=1
    for(k1 in 1:p){

      lambdaj_post_rate1=abs(alpha[k1])+Eta1
      lambdaj1[k1]<-rgamma(1,shape = neu1+1,rate=lambdaj_post_rate1)
      lambdaj_sq1[k1]=(lambdaj1[k1])^2

      tauj1[k1] <-rgig(1,0.5,chi = (alpha[k1])^2,psi = lambdaj_sq1[k1])
    }

    # Update latent class indicator y1 (w in paper)
    eta<-X%*%alpha
    phi<-exp(eta)/(1+exp(eta))              # at-risk probability
    theta<-phi*(q^r)/(phi*(q^r)+1-phi)      # Conditional prob that y1=1 given y=0 -- i.e. Pr(chance zero|observed zero)
    y1[y==0]<-rbinom(n0,1,theta[y==0])      # If y=0, then draw a "chance zero" w.p. theta, if y=1, then y1=1
    n1<-sum(y1)

    # Update rho                              # Proposal variance
    rnew<-rtnorm(1,r,sqrt(s),lower=0)       # Treat r as continuous
    ratio<-sum(dnbinom(y[y1==1],rnew,q[y1==1],log=T))-sum(dnbinom(y[y1==1],r,q[y1==1],log=T))+
      dtnorm(r,rnew,sqrt(s),0,log=T) - dtnorm(rnew,r,sqrt(s),0,log=T) # Note: Uniform prior for r
    #       Proposal not symmetric
    if (log(runif(1))<ratio) {
      r<-rnew
      Acc<-Acc+1
    }

    # Update beta
    eta1<-X[y1==1,]%*%beta
    pos.int<-ifelse(y[y1==1]+r<1,1,as.integer(y[y1==1]+r))
    w<-pgdraw(pos.int,eta1)     #n1                           # Polya weights
    z<-(y[y1==1]-r)/(2*w)
    sigma <-diag(as.numeric(sig_sq%*%t(tauj_sq)),p,p) # Latent "response"
    v<-solve(crossprod(X[y1==1,]*sqrt(w))+sigma)
    m<-v%*%(sigma%*%beta0+t(sqrt(w)*X[y1==1,])%*%(sqrt(w)*z))
    beta<-c(rmvnorm(1,m,v))
    eta<-X%*%beta
    q<-1/(1+exp(eta))

    tauj=lambdaj=1
    for(k in 1:p){
      #k=1
      lambdaj_post_rate=abs(beta[k])+Eta
      lambdaj[k]<-rgamma(1,shape = neu+1,rate=lambdaj_post_rate)
      lambdaj_sq[k]=(lambdaj[k])^2

      tauj[k] <-rgig(1,0.5,chi = (beta[k])^2,psi = lambdaj_sq[k])
    }
    # Update r2 using Gibbs as in Dadaneh et al and Zhou and Carin #

    # Update latent counts, l
    l<-rep(0,n1)
    ytmp<-y[y1==1]
    for(j in 1:n1) l[j]<-sum(rbinom(ytmp[j],1,round(r/(r+1:ytmp[j]-1),6))) # Could try to avoid loop; rounding avoids numerical stability

    # Update r from conjugate gamma distribution given l and psi
    eta<-X[y1==1,]%*%beta
    psi<-exp(eta)/(1+exp(eta))
    r2<-rgamma(1,0.01+sum(l),0.01-sum(log(1-psi))) # Gamma(0.01,0.01) prior for r



    # Store
    if (i> burnin & i%%thin==0) {
      j<-(i-burnin)/thin
      Alpha[j,]<-alpha
      Beta[j,]<-beta
      R[j]<-r
      R2[j]<-r2
      #sig_sq_all = c(sig_sq_all)
      lambda_sq_all=c(lambda_sq_all)
      tau_sq_all=c(tau_sq_all)
      #sig_sq_all1 = c(sig_sq_all1)
      lambda_sq_all1=c(lambda_sq_all1)
      tau_sq_all1=c(tau_sq_all1)
    }

    if (i%%100==0) print(i)

  }
  MC_DPzinb<-list(Alpha=Alpha,Beta=Beta, R=R,R2=R2,
                  lambda=lambda_sq_all, tausq=tau_sq_all,
                  lambda1=lambda_sq_all1,
                  tausq1=tau_sq_all1)
}
####################################################################################################################

