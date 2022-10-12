
#############Horseshoe###################################

#' Logistic regression with Horseshoe prior
#'
#' @param Y The binary response variable
#' @param X The design matrix
#' @param mc_length no. of MCMC iterations
#' @param beta_start initial parameters
#'
#' @export list of posterior samples of the parameters
#' @return sample_beta: posterior samples of the coefficients beta
#' @return W_all: posterior samples of the latent variable
#' @return xi_all: posterior samples of hyper-parameter
#' @return gamma_all: posterior samples of hyper-parameter
#' @return lambda_all: posterior samples of local parameter in the shrinkage prior
#' @return tau_sq_all: posterior samples of global parameter in the shrinkage prior

#' @examples
#' X<-matrix(rnorm(400), ncol=4,nrow=100)
#' Y<-rbinom(100,1,0.5)
#' mc_length=100
#' beta_start<-rep(0.5,4)
#' horseshoe_logistic_mc(Y=Y, X=X, mc_length=mc_length,beta_start=beta_start )

horseshoe_logistic_mc<-function(Y, X, mc_length,beta_start ){
  p=dim(X)[2]
  n=length(Y)
  Ystar<-(Y-0.5)

  #### initialization######################
  #########################################
  #########################################
  beta<-beta_start
  #########################################
  w<-pgdraw(1,X%*%beta) ###


  wd<-diag(w);dim(wd); ## diagonal matrix wd

  xi<-1     ## xi
  tau_sq <-1 ## tau^2
  gamma_j<-as.matrix(replicate(n = p,expr = 1)) # gamma_j
  lambdaj_sq<-as.matrix( replicate(n = p,expr = 1))
  sigma <-.01*diag(as.numeric(tau_sq%*%t(lambdaj_sq)),p,p)# sigma ##

  #########################################GIBBS SAMPLING#########################################
  #################################################################################################
  #################################################################################################
  #################################################################################################
  #################################################################################################
  sample_beta<-array(NA, dim=c(mc_length,p))
  W_all= array(NA, dim=c(mc_length,n))
  xi_all= array(NA, dim=c(mc_length,1))
  gamma_all=array(NA, dim=c(mc_length,p))
  lambda_all=array(NA, dim=c(mc_length,p))
  tau_sq_all=array(NA, dim=c(mc_length,1))
  #browser()
  for(i in 1:mc_length){
    w <-pgdraw(1,X%*%beta) #Augmentation
    beta_post_var=solve(t(X)%*%diag(w)%*%X + solve(sigma)) # beta_variance
    beta_post_mu= beta_post_var%*%t(X)%*%Ystar            # beta_mean
    beta <-c(mvrnorm(1,beta_post_mu ,beta_post_var))   # sample from multivariate normal for beta

    xi <-rigamma(n = 1,alpha = 1,beta = (1+1/tau_sq))

    tau_sq_post_rate=(1/xi)+.5*sum((beta^2)/lambdaj_sq)
    tau_sq <-rigamma(1,(p+1)/2,tau_sq_post_rate )

    for(j in 1:p){

      gamma_j[j] <-rigamma(1,1,(1+(1/lambdaj_sq[j])))
      lambdaj_sq_post_rate=((1/gamma_j[j])+(0.5*(beta[j])^2/tau_sq))
      lambdaj_sq[j]<-rigamma(1,1,lambdaj_sq_post_rate)

    }

    sigma <-diag(as.numeric(tau_sq%*%t(lambdaj_sq)),nrow=p,ncol=p) #

    ##storing of values:
    sample_beta[i,]<-c(beta)
    W_all[i,]<-c(w)
    xi_all[i,]= c(xi)
    gamma_all[i,]<-c(gamma_j)
    lambda_all[i,]<-c(lambdaj_sq)
    tau_sq_all[i,]<-c(tau_sq)

  }
  MC_logisticReg<-list(sample_beta,W_all,xi_all,gamma_all,lambda_all, tau_sq_all)

  return(MC_logisticReg)
}


##################DirichletLaplace##########################
#' Logistic regression with Dirichlet Laplace prior
#'
#' @param Y The binary response variable
#' @param X The design matrix
#' @param mc_length no. of MCMC iterations
#' @param beta_start initial parameters
#'
#' @export list of posterior samples of the parameters
#' @return sample_beta: posterior samples of the coefficients beta
#' @return W_all: posterior samples of the latent variable
#' @return shi_all: posterior samples of the hyper-parameters
#' @return phi_all: posterior samples of the hyper-parameters
#' @return tau_sq_all: posterior samples of the hyper-parameters

#' @examples
dirichlet_laplace_logistic_mc<-function(Y, X, mc_length,a,beta_start){

  p=dim(X)[2]
  n=length(Y)
  Ystar<-(Y-0.5)
  cutOff=.0000000000001
  #### initialization######################
  #########################################
  #########################################

  #########################################
  #beta=rep(0,p)
  beta=beta_start
  w<-pgdraw(1,X%*%beta)###??

  wd<-diag(w);dim(wd); ## diagonal matrix wd

  #a<-1/p     ## a
  tau_sq <-1 ## tau^2
  shi=(replicate(n = p,expr = 1)) # gamma_j
  #browser()
  phi<-( replicate(n = p,expr = 1)  )/p
  sigma <-.1*diag(as.numeric(shi*t(phi^2)*tau_sq),p,p)# sigma ##
  sigma_inv=solve(sigma)
  #########################################GIBBS SAMPLING#########################################
  sample_beta<-array(NA, dim=c(mc_length,p))

  W_all= array(NA, dim=c(mc_length,n))
  shi_all= array(NA, dim=c(mc_length,p))
  phi_all=array(NA, dim=c(mc_length,p))

  tau_sq_all=array(NA, dim=c(mc_length,1))

  #i=1
  for(i in 1:mc_length){
    # print(i)
    w <-pgdraw(1,abs(X%*%(beta))) #Augmentation

    beta_post_var=solve((t(X)%*%diag(w)%*%X+sigma_inv))#, error=function(e){print("error"); print((t(X)%*%diag(w)%*%X+sigma_inv));stop("a")} ) # beta_variance
    beta_post_mu= beta_post_var%*%t(X)%*%Ystar            # beta_mean
    beta <-c(mvrnorm(n = 1,mu = beta_post_mu ,Sigma = beta_post_var))   # sample from multivariate normal for beta

    #xi <-rigamma(n = 1,alpha = 1,beta = (1+1/tau_sq))

    #tau--bhatta
    tau= rgig(n=1,lambda=p*a-p, chi= 2*sum(abs(beta)/phi) ,psi=1  )
    tau_sq= tau^2

    #phi --bhatta
    T_phi=rep(NA,p)
    for(j in 1:p){
      #   #shi[j] <-rgig(1,1/2,(beta[j])^2/(tau_sq*phi[j]^2),1)
      T_phi[j] <- rgig(n=1,lambda = a-1,chi = 2*abs(beta[j]),psi = 1)
    }

    phi=T_phi/sum(T_phi)


    for(j in 1:p){
      #tau_sq_post_rate=(beta[j])^2/(shi[j]*phi[j]^2)
      #tau_sq <-rgig(1,n*a-(p/2),tau_sq_post_rate,1)
      shi[j] <-rgig(n=1,lambda = 1/2,chi=(beta[j])^2/(tau_sq*(phi[j])^2),psi =1)

    }

    #browser()
    #print(tau_sq)
    cand=shi*(phi^2)*tau_sq
    # if(min(cand)<cutOff){print(paste("here=",T_phi));cand[which(cand<cutOff)]=cutOff}
    sigma <-diag(cand,p,p)
    sigma_inv<-diag(1/(cand),p,p)

    ##storing of values:
    sample_beta[i,]<-c(beta)
    W_all[i,]<-c(w)
    shi_all[i,]= c(shi)
    phi_all[i,]<- c(phi)
    tau_sq_all[i,]<-c(tau_sq)

  }
  MC_logisticReg<-list(sample_beta,W_all,shi_all,phi_all, tau_sq_all)

  return(MC_logisticReg)
}
###########################Double Pareto#####################################
#' Logistic Regression with Double Pareto Prior
#'
#' @param Y The binary response variable
#' @param X The design matrix
#' @param mc_length no. of MCMC iterations
#' @param beta_start initial parameters
#'
#' @export list of posterior samples of the parameters
#' @return sample_beta: posterior samples of the coefficients beta
#' @return W_all: posterior samples of the latent variable
#' @return lambda_sq_all: posterior samples of the local parameter in the shrinkage prior
#' @return tau_sq_all: posterior samples of the global parameter in the shrinkage prior
#'
#' @examples
double_pareto_logistic_mc<-function(Y, X, mc_length,beta_start ){
  p=dim(X)[2]
  n=length(Y)
  Ystar<-(Y-0.5)


  ##### initialization######################
  #########################################
  #########################################
  beta<-beta_start
  #########################################
  w<-pgdraw(1,X%*%beta)###??

  wd<-diag(w);dim(wd); ## diagonal matrix wd
  sig_sq<-1
  neu<-1
  xi<-1     ## xi
  alpha<-1
  eta<-1
  tauj_sq <-(replicate(n = p,expr = 1)) ## tau^2

  lambdaj_sq<-( replicate(n = p,expr = 1)  )
  sigma <-.1*diag(as.numeric(sig_sq%*%t(tauj_sq)),p,p)

  #########################################GIBBS SAMPLING#########################################
  #################################################################################################
  #################################################################################################
  #################################################################################################
  #################################################################################################
  sample_beta <-array(NA, dim=c(mc_length,p))

  W_all = array(NA, dim=c(mc_length,n))
  sig_sq_all = array(NA, dim=c(mc_length,1))

  lambda_sq_all=array(NA, dim=c(mc_length,p))
  tau_sq_all=array(NA, dim=c(mc_length,p))
  #browser()
  for(i in 1:mc_length){
    w <-pgdraw(1,abs(X%*%(beta))) #Augmentation
    beta_post_var=solve((t(X)%*%diag(w)%*%X+solve(sigma))) # beta_variance
    #beta_post_var=solve((t(X)%*%diag(w)%*%X)) # beta_variance
    beta_post_mu= beta_post_var%*%t(X)%*%Ystar            # beta_mean
    beta <-c(mvrnorm(1,beta_post_mu ,beta_post_var))   # sample from multivariate normal for beta

    #if(tauj_sq<.000000001){      print("Tauj is small ");print(tauj_sq) }
    #sig_sq<- rigamma(n=1, (alpha+p),(xi+0.5*sum(beta^2)/tauj_sq))

    tauj=lambdaj=0*tauj_sq-1
    for(j in 1:p){

      lambdaj_post_rate=(abs(beta[j])+eta)
      #lambdaj_sq[j]<-rigamma(1,neu+1,lambdaj_sq_post_rate)
      lambdaj[j]<-rgamma(1,shape = neu+1,rate=lambdaj_post_rate)
      lambdaj_sq[j]=(lambdaj[j])^2

      tauj[j] <-rgig(1,0.5,chi = (beta[j])^2,psi = lambdaj_sq[j])
      #lambdaj_sq_post_rate=((tauj_sq/2)+eta)


    }

    sigma <-diag((tauj),p,p)

    ##storing of values:
    sample_beta[i,]<-c(beta)
    W_all[i,]<-c(w)
   # sig_sq_all[i,]= c(sig_sq)
    lambda_sq_all[i,]<-c(lambdaj_sq)
    tau_sq_all[i,]<-c(tauj_sq)

  }
  MC_logisticReg<-list(sample_beta,W_all,lambda_sq_all,tau_sq_all)

  return(MC_logisticReg)
}

