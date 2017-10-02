
library("genlasso")
library("MASS")

#function that generates observation from a normal density and a mixture of Gaussian priors centered at 0 and 
#priors with variances prior_var_1 and prior_var_2
generative_model<-function(X, prior_var_1 = 0.1, prior_var_2= 3., p= 50)
{
  u = runif(1,min= 0,max=1)
  beta = rep(0,p)
  for(k in 1:p){
    if(u <= 0.9){
      beta[k] = rnorm(1, mean = 0, sd = sqrt(prior_var_1))
    }else{
      beta[k] = rnorm(1, mean = 0, sd = sqrt(prior_var_2))
    }
  }
  y = mvrnorm(n = 1, mu = X%*%beta, Sigma = diag(n))
  return(list("beta"=beta, "y"= y))
}

##functions that performs Lasso at fixed tuning parameter and returns average coverage across selected for the
#population least squares coefficients
simulation.unadjusted.noninform<- function(prior_var = 1000., p=50, lambda.lasso)
{
  data.par = generative_model(X)
  y = data.par$y
  beta = data.par$beta
  D = diag(1,p)
  out = genlasso(y, X=X, D=D)
  
  c<-coef(out,lambda=lambda.lasso)
  coef<-as.vector(c$beta)
  E<-which(abs(coef)>0.0001)
  X.E<-X[,E]
  P.E = X.E%*%solve(t(X.E)%*%X.E)
  
  M.1 = prior_var*(X%*%t(X)) + diag(n)
  M.2 = prior_var*(X%*%t(X))%*%P.E
  M.3 = prior_var*t(P.E)%*%(X%*%t(X))%*%P.E
  mu0 = t(M.2) %*% solve(M.1)%*% y
  var = M.3 - t(M.2) %*% solve(M.1)%*% M.2
  
  UBC<-NULL
  LBC<-NULL
  coverage = 0.
  true.target = t(P.E)%*%X%*%beta
  
  for(i in 1:length(E))
  {
    LBC[i]<-mu0[i]-1.96*(sqrt(var[i,i]))
    UBC[i]<-mu0[i]+1.96*(sqrt(var[i,i]))
    coverage = coverage + ifelse((LBC[i]<true.target[i]) & (UBC[i]>true.target[i]),1,0)
  }
  
  return(coverage/length(E))
}

##functions that performs Lasso at fixed tuning parameter and returns average coverage across selected for the
#population least squares coefficients; the model for inference is the selected model
simulation.unadjusted.selected<- function(prior_var = 1000., p=50, lambda.lasso)
{
  data.par = generative_model.fixed(X)
  y = data.par$y
  beta = data.par$beta
  D = diag(1,p)
  out = genlasso(y, X=X, D=D)
  
  c<-coef(out,lambda=lambda.lasso)
  coef<-as.vector(c$beta)
  E<-which(abs(coef)>0.0001)
  X.E<-X[,E]
  P.E = X.E%*%solve(t(X.E)%*%X.E)
  
  Q<-(((prior_var)*X.E%*%t(X.E))+diag(n))
  mu0<-(prior_var)*(t(X.E)%*%solve(Q)%*%y)
  var<-(prior_var*diag(length(E)))-(prior_var^2)*(t(X.E)%*%solve(Q)%*%(X.E))
  
  
  UBC<-NULL
  LBC<-NULL
  coverage = 0.
  true.target = beta[E]
  for(i in 1:length(E))
  {
    LBC[i]<-mu0[i]-1.96*(sqrt(var[i,i]))
    UBC[i]<-mu0[i]+1.96*(sqrt(var[i,i]))
    coverage = coverage + ifelse((LBC[i]<true.target[i]) & (UBC[i]>true.target[i]),1,0)
  }
  return(coverage/length(E))
}
###test for coverage of unadjusted intervals
test.unadjusted.coverage<- function(n, p)
{ 
  #scale X
  X = matrix(rnorm(n*p), ncol=p)
  X = apply(X, 2, function(x) x/sqrt(sum(x^2)))
  #find theoretical lambda
  epsilon = mvrnorm(n=2000,mu=rep(0,n),Sigma=((sigma^2)*diag(1,n)))
  sum = 0
  for(i in 1: 2000)
  {
    sum = sum +  sigma*max(t(X)%*%epsilon[i,])
  }
  lam = 0.7*sum/2000

  tot.coverage = 0
  for(rep in 1: 100)
  {
    cov.0 = simulation.unadjusted.noninform(prior_var = 1000., p=50, lambda.lasso =lam)
    tot.coverage  = tot.coverage + cov.0
  }
  return(tot.coverage)
}