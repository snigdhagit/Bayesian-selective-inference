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
#population least squares coefficients under the selected model
simulation.adjusted.selected<- function(prior_var = 1000., p=50, lambda.lasso)
{
  data.par = generative_model.fixed(X)
  y = data.par$y
  beta = data.par$beta
  D = diag(1,p)
  out = genlasso(y, X=X, D=D)
  
  c<-coef(out,lambda=lambda.lasso)
  coef<-as.vector(c$beta)
  E<-which(abs(coef)>0.0001)
  X1<-X[,E]
  l=length(E)
  z=sign(coef[E])
  
  #computes polyhedral map
  selection<-function(lam)
  {
    #computing lasso solution at particular value of lambda=lam
    l<-length(E)
    X1<-X[,E]
    X2<-X[,-E]
    P1<-solve(t(X1)%*%X1)
    P2<-X1%*%P1
    P<-X1%*%P1%*%(t(X1))
    I<-diag(1,n)
    A0<-(1/lam)*rbind((t(X2)%*%(I-P)),-(t(X2)%*%(I-P)))
    A1<--diag(z)%*%P1%*%(t(X1))
    b0<-rbind((1-(t(X2)%*%P2%*%z)),(1+(t(X2)%*%P2%*%z)))
    b1<-(-lam)*(diag(z)%*%P1%*%z)
    A<-rbind(A0,A1)
    b<-rbind(b0,b1)
    L<-list("mat"=A,"vec"=b,"pro"=P2,"sel"=X1,"len"=l)
    return(L)
  }
  
  L<-selection(lam)
  A<-L$mat
  b<-L$vec
  P2<-L$pro
  X1<-L$sel
  l<-L$len
  
  #barrier function
  barrier<-function(z1)
  {
    if(all((b-A%*%z1)>0.00000000000000000001)==TRUE){
      return(sum(log(1+(1/(b-A%*%z1)))))} else{
        return(length(b)*(log(1+(1/0.00000000000000000001))))}
  }
  
  #print((all((b-A%*%y)>0.00000000000000000001)))
  
  #optimization obejctive for affine selection probability
  objective<-function(z1,beta1)
  {
    return(((1/2)*((norm_vec((X1%*%beta1)-z1))^2))+barrier(z1))
  }
  
  #gradient of optimization
  grad<-function(beta)
  {
    value=nlm(objective,p=y,beta1=beta,hessian=TRUE)$estimate
    return(((-beta/prior_var)+(t(X1)%*%(y-value))))
  }
  
  #run Langevin sampler to obtain samples for approximate posterior
  gamma=1/l
  g.samp=mvrnorm(n=1500,mu=rep(0,l),Sigma=diag(1,l))
  chain<-matrix(nrow=1500,ncol=l,0)
  chain[1,] = coef[E]
  
  for(i in 1:1499)
  {
    chain[i+1,]=chain[i,]+(gamma*grad(chain[i,]))+((sqrt(2*gamma))*g.samp[i,])
    print(i)
    print(chain[i,])
  }
  
  chain0= chain
  chain.bn=chain[-(1:100),]
  bci=matrix(nrow=l,ncol=2)
  for(i in 1:l)
  {
    bci[i, ] <- quantile(chain.bn[, i], probs = c(0.025, 0.975))
  }
  
  AC_langevin_barrier<- bci
  coverage = 0.
  true.target = beta[E]
  for(j in 1:length(E))
  {
    coverage = coverage + ifelse((bci[j]<true.target[j]) & (bci[j]>true.target[j]),1,0)
  }
  
  return(coverage/length(E))
}

###test to compute coverage
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
  
  tot.ad.coverage.selected =0
  for(rep in 1: 10)
  {
    cov = simulation.adjusted.selected(prior_var = 1000., p=50, lambda.lasso =lam)
    tot.ad.coverage.selected = tot.coverage + cov
    print("coverage so far")
    print(tot.ad.coverage.selected)
  }
  
  return(tot.ad.coverage.selected)
}