
#######Horseshoe###########
#' Multinomial logistic regression with Horseshoe prior
#'
#' @param y The categorical variable with more than two classes
#' @param X The design matrix
#' @param n sample size
#' @param m.0 initial mean parameters mostly zero, a matrix of dimension P*J-1, p=no. of variables, J=no. of categorical classes
#' @param sigma initial variance covariance matrix of dimension P*P*J-1
#' @param samp no. of MCMC samples, generally 1000
#' @param burn burnin samples, generally 500
#' @export list of posterior samples of the parameters
#' @return sample_beta: posterior samples of the coefficients beta
#' @return w: posterior samples of the latent variable
#'
#'
#' @examples
horseshoe_mlogit_mc <- function(y, X, n=rep(1,nrow(as.matrix(y))),
                     m.0=array(0, dim=c(ncol(X), ncol(y))), #dim p*J-1 mean
                     sigma=array(0, dim=c(ncol(X), ncol(X), ncol(y))), #P*P*J-1
                     samp=1000, burn=500){
  ## N - number of trials
  ## J - number of categories
  ## P - number of covariates
  ## y - N x J-1 matrix.  y_{ij} = fraction of outcomes in jth category on trial i.
  ## X - N x P design matrix
  ## n - N x 1 matrix of rolls   ##??
  ## Assume beta_J = 0 for identification.

  X = as.matrix(X);
  y = as.matrix(y);

  N = nrow(X);
  P = ncol(X);
  J = ncol(y) + 1;

  out.h = list(
    beta = array(0, dim=c(samp, P, J-1)), #> dim(out$beta)[1] 1000    4    2
    # 1000 simulations, 4 covariates, 2 class
    w    = array(0, dim=c(samp, N, J-1)),#dim(out$w)[1] 1000  100    2
    # 1000 simulations, 100 obvsns, 2 class
    xi =  rep(0,samp),
    tau_sq<-rep(0,samp),
    gamma<-vector(mode = "list", length = samp),# gamma_k
    lambda_sq<-vector(mode = "list", length = samp)#array(0,dim=c(samp,P,J-1)) # lambda_k_sq
  )

  beta = matrix(0, P, J);
  w    = matrix(0, N, J);

  ## Precompute.
  kappa = (y - 0.5)*n;

  xi<-0.1    ## xi
  tau_sq <-0.1 ## tau^2
  gamma<-matrix(0.1,P,J) # gamma_k
  lambda_sq<-matrix(0.1,P,J) # lambda_k_sq

  for (j in 1:(J-1))
    sigma[,,j] <-.01*diag(as.numeric(tau_sq%*%t(lambda_sq[,j])),P,P)


  b.0 = matrix(0, P, J-1);
  #m.0=array(0, dim=c(ncol(X), ncol(y)))
  for (j in 1:(J-1)) b.0[,j] = sigma[,,j] %*% m.0[,j]; #??

  ## A = rowSums( exp(X %*% beta[,-1]) );
  for (i in 1:(samp+burn)) {
    for (j in 1:(J-1)) {


      A = rowSums( exp(X %*% beta[,-j]) ); #N*1

      c.j   = log(A); #N*1
      eta.j = X %*% beta[,j] - c.j; #N*1

      ## omega.j
      w[,j] = pgdraw(n, eta.j)#rpg.devroye(N, n, eta.j) #N*3, last column zero

      ## beta.j
      PL.j = t(X) %*% (X * w[,j]);
      bL.j = t(X) %*% (kappa[,j] + c.j * w[,j]); ###??

      P1.j = PL.j + sigma[,,j];
      ## Can speed up using Choleksy.
      V1.j = chol2inv(chol(P1.j)); ###??
      m1.j = V1.j %*% (bL.j + b.0[,j]);

      sqrtV1.j = t(chol(V1.j));
      beta[,j] = m1.j + sqrtV1.j %*% rnorm(P);

      ## A += exp(X %*% beta[,j]) - exp(X %*% beta[,(j+1) %% (J-1)])
      xi <-rigamma(n = 1,alpha = 1,beta = (1+1/tau_sq))

      tau_sq_post_rate=(1/xi)+.5*sum((beta[,j]^2)/lambda_sq[,j])
      tau_sq <-rigamma(1,(P+1)/2,tau_sq_post_rate )

      for(k in 1:P){

        gamma[k,j] <-rigamma(1,1,(1+(1/lambda_sq[k,j])))
        lambda_sq_post_rate=(1/gamma[k,j])+(0.5*(beta[,j])^2/tau_sq)
        lambda_sq[k,j]<-rigamma(1,1,lambda_sq_post_rate[k])

      }

      sigma[,,j] <-diag(as.numeric(tau_sq%*%t(lambda_sq[,j])),nrow=P,ncol=P)
      ## Store
      if (i > burn) {
        out.h$beta[i-burn,,j] = beta[,j];
        out.h$w[i-burn,,j] = w[,j];
        out.h$tau_sq[i-burn]<-tau_sq;
        out.h$xi[i-burn]<-xi;
        #beta = array(0, dim=c(samp, P, J-1))
        out.h$gamma[[i-burn]]<-gamma;# gamma_k
        out.h$lambda_sq[[i-burn]]<-lambda_sq;
      }
    }
    # if (i %% verbose == 0) cat("Finished", i, "\n");
  }

  ab<- list(finalbeta=out.h$beta,
            w= out.h$w,
            X=X,
            y=y,
            n = n,
            xi=out.h$xi,
            lambda_sq=out.h$lambda_sq,
            gamma= out.h$gamma,
            tau_sq=out.h$tau_sq )

  return(ab)
}

###########Double Pareto###################

#' Multinomial logistic regression with Double Pareto prior
#'
#' @param y The categorical variable with more than two classes
#' @param X The design matrix
#' @param n sample size
#' @param m.0 initial mean parameters mostly zero, a matrix of dimension P*J-1, p=no. of variables, J=no. of categorical classes
#' @param sigma initial variance covariance matrix of dimension P*P*J-1
#' @param samp no. of MCMC samples, generally 1000
#' @param burn burnin samples, generally 500
#'
#' @export list of posterior samples of the parameters
#' @return sample_beta: posterior samples of the coefficients beta
#' @return w: posterior samples of the latent variable
#'
#'
#' @examples
double_pareto_mlogit_mc <- function(y, X, n=rep(1,nrow(as.matrix(y))),
                      m.0=array(0, dim=c(ncol(X), ncol(y))),
                      sigma=array(0, dim=c(ncol(X), ncol(X), ncol(y))),
                      samp=1000, burn=500)

{
  ## N - number of trials
  ## J - number of categories
  ## P - number of covariates
  ## y - N x J-1 matrix.  y_{ij} = fraction of outcomes in jth category on trial i.
  ## X - N x P design matrix
  ## n - N x 1 matrix of rolls
  ## Assume beta_J = 0 for identification.

  X = as.matrix(X);
  y = as.matrix(y);

  N = nrow(X);
  P = ncol(X);
  J = ncol(y) + 1;

  out.dp = list(
    beta = array(0, dim=c(samp, P, J-1)),
    w    = array(0, dim=c(samp, N, J-1))
    # lambdaj_sq<-array(0, dim=c(samp, P, J-1)),
    # tauj<-array(0, dim=c(samp, P, J-1))
  )

  beta = matrix(0, P, J);
  w    = matrix(0, N, J);

  ## Precompute.
  kappa = (y - 0.5)*n;
  sig_sq<-1
  neu<-1
  xi<-1     ## xi
  alpha<-1
  eta<-1
  tauj_sq <-matrix(0.1,P,J)#(replicate(n = p,expr = 1)) ## tau^2

  lambdaj<-matrix(0.1,P,J)
  lambdaj_sq<-matrix(0.1,P,J)
  #lambdaj_sq<-matrix(0.1,P,J)#( replicate(n = p,expr = 1)  )

  for (j in 1:(J-1))
    sigma[,,j] <-0.01*diag(as.numeric(sig_sq%*%t(tauj_sq[,j])),P,P)
  b.0 = matrix(0, P, J-1);
  for (j in 1:(J-1))
    b.0[,j] = sigma[,,j] %*% m.0[,j];
  ## A = rowSums( exp(X %*% beta[,-1]) );
  for (i in 1:(samp+burn)) {
    for (j in 1:(J-1)) {


      A = rowSums( exp(X %*% beta[,-j]) );

      c.j   = log(A);
      eta.j = X %*% beta[,j] - c.j;


      w[,j] = pgdraw(n, eta.j)#rpg.devroye(N, n, eta.j);

      ## beta.j
      PL.j = t(X) %*% (X * w[,j]);
      bL.j = t(X) %*% (kappa[,j] + c.j * w[,j]);

      P1.j = PL.j + sigma[,,j];
      ## Can speed up using Choleksy.
      V1.j = chol2inv(chol(P1.j));
      m1.j = V1.j %*% (bL.j + b.0[,j]);

      sqrtV1.j = t(chol(V1.j));
      beta[,j] = m1.j + sqrtV1.j %*% rnorm(P);

      ## A += exp(X %*% beta[,j]) - exp(X %*% beta[,(j+1) %% (J-1)])
      lambdaj_post_rate=(abs(beta[,j])+eta)
      for(k in 1:P){

        lambdaj[k,j]<-rgamma(1,shape = neu+1,rate=lambdaj_post_rate[k])
        lambdaj_sq[k,j]=(lambdaj[k,j])^2
        tauj_sq[k,j] <-rgig(1,0.5,chi = (beta[,j])^2,psi = lambdaj_sq[k,j])

      }

      sigma[,,j] <-diag(as.numeric(sig_sq%*%t(tauj_sq[,j])),P,P)
      ## Store
      if (i > burn) {
        out.dp$beta[i-burn,,j] = beta[,j];
        out.dp$w[i-burn,,j] = w[,j];

      }
    }
    # if (i %% verbose == 0) cat("Finished", i, "\n");
  }


  ab<- list(finalbeta=out.dp$beta,
            w= out.dp$w,
            X=X,
            y=y,
            n= n)

  return(ab)
}


###########Dirichlet Laplace################################

#' Multinomial logistic regression with Dirichlet Laplace prior
#'
#' @param y The categorical variable with more than two classes
#' @param X The design matrix
#' @param n sample size
#' @param m.0 initial mean parameters mostly zero, a matrix of dimension P*J-1, p=no. of variables, J=no. of categorical classes
#' @param sigma initial variance covariance matrix of dimension P*P*J-1
#' @param samp no. of MCMC samples, generally 1000
#' @param burn burnin samples, generally 500
#'
#' @export list of posterior samples of the parameters
#' @return sample_beta: posterior samples of the coefficients beta
#' @return w: posterior samples of the latent variable
#'
#'
#' @examples
dirichlet_laplace_mlogit_mc <- function(y, X, n=rep(1,nrow(as.matrix(y))),
                      m.0=array(0, dim=c(ncol(X), ncol(y))), #dim p*J-1 mean
                      sigma=array(0, dim=c(ncol(X), ncol(X), ncol(y))), #P*P*J-1
                      samp=1000, burn=500, a=0.8)

{
  ## N - number of trials
  ## J - number of categories
  ## P - number of covariates
  ## y - N x J-1 matrix.  y_{ij} = fraction of outcomes in jth category on trial i.
  ## X - N x P design matrix
  ## n - N x 1 matrix of rolls
  ## Assume beta_J = 0 for identification.

  X = as.matrix(X);
  y = as.matrix(y);

  N = nrow(X);
  P = ncol(X);
  J = ncol(y) + 1;

  out.dl = list(
    beta = array(0, dim=c(samp, P, J-1)), #> dim(out$beta)[1] 1000    4    2
    # 1000 simulations, 4 covariates, 2 class
    w    = array(0, dim=c(samp, N, J-1))#dim(out$w)[1] 1000  100    2
    # 1000 simulations, 100 obvsns, 2 class


  )

  beta = matrix(0, P, J);
  w    = matrix(0, N, J);

  ## Precompute.
  kappa = (y - 0.5)*n;

  tau_sq <-0.1 ## tau^2
  shi<-matrix(0.1,P,J)
  T_phi=matrix(0.1,P,J)
  phi=matrix(0.1,P,J)

  for (j in 1:(J-1))
    sigma[,,j] <-.01*diag(as.numeric(shi[,j]%*%t(phi[,j]^2)*tau_sq),P,P)

  b.0 = matrix(0, P, J-1);
  for (j in 1:(J-1)) b.0[,j] = sigma[,,j] %*% m.0[,j];

  ## A = rowSums( exp(X %*% beta[,-1]) );
  for (i in 1:(samp+burn)) {
    for (j in 1:(J-1)) {


      A = rowSums( exp(X %*% beta[,-j]) );

      c.j   = log(A);
      eta.j = X %*% beta[,j] - c.j;


      w[,j] = pgdraw(n, eta.j)#rpg.devroye(N, n, eta.j);

      ## beta.j
      PL.j = t(X) %*% (X * w[,j]);
      bL.j = t(X) %*% (kappa[,j] + c.j * w[,j]);

      P1.j = PL.j + sigma[,,j];
      ## Can speed up using Choleksy.
      V1.j = chol2inv(chol(P1.j));
      m1.j = V1.j %*% (bL.j + b.0[,j]);

      sqrtV1.j = t(chol(V1.j));
      beta[,j] = m1.j + sqrtV1.j %*% rnorm(P);

      ## A += exp(X %*% beta[,j]) - exp(X %*% beta[,(j+1) %% (J-1)])
      tau= rgig(n=1,lambda=P*a-P, chi= 2*sum(abs(beta[,j])/phi[,j]) ,psi=1  )
      tau_sq= tau^2

      #phi --bhatta

      for(k in 1:P){

        T_phi[k,j] <- rgig(n=1,lambda = a-1,chi = 2*abs(beta[,j]),psi = 1)
        #shi[k,j] <-rgig(n=1,lambda = 1/2,chi=(beta[,j])^2/(tau_sq*(phi[j])^2),psi =1)
        phi[k,j]=T_phi[k,j]/sum(T_phi[k,j])
      }




      for(k in 1:P){

        shi[k,j] <-rgig(n=1,lambda = 1/2,chi=(beta[,j])^2/(tau_sq*(phi[k])^2),psi =1)

      }

      #cand=shi*(phi^2)*tau_sq
      # if(min(cand)<cutOff){print(paste("here=",T_phi));cand[which(cand<cutOff)]=cutOff}
      sigma[,,j] <-.01*diag(as.numeric(shi[,j]%*%t(phi[,j]^2)*tau_sq),P,P)
      ## Store
      if (i > burn) {
        out.dl$beta[i-burn,,j] = beta[,j];
        out.dl$w   [i-burn,,j] = w[,j];
      }
    }
    # if (i %% verbose == 0) cat("Finished", i, "\n");
  }


  ab<- list(finalbeta=out.dl$beta,
            w= out.dl$w,
            X=X,
            y=y,
            n= n)

  return(ab)
}

