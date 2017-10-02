####data carved inference: computing adjusted intervals and MAP where selection is performed on a perturbed instance
#of data (data randomized by gaussian noise with 10% variance as original data)

#y.1 is the randomized response, alpha = sigma^2 + tau^2; tau^2 is the randomization variance
sampler<- function(X, X.sel, y, E, z, l, lam, alpha, sigma)
{
  selection<-function(lam)
  {
    #computing lasso solution at particular value of lambda=lam
    l<-length(E)
    X1.sel<-X.sel[,E]
    X2.sel<-X.sel[,-E]
    P1<-solve(t(X1.sel)%*%X1.sel)
    P2<-X1.sel%*%P1
    P<-X1.sel%*%P1%*%(t(X1.sel))
    I<-diag(1,n.sel)
    A0<-(1/lam)*rbind((t(X2.sel)%*%(I-P)),-(t(X2.sel)%*%(I-P)))
    A1<--diag(z)%*%P1%*%(t(X1.sel))
    b0<-rbind((1-(t(X2.sel)%*%P2%*%z)),(1+(t(X2.sel)%*%P2%*%z)))
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
  
  #data selection matrix
  indicator = rep(0,442)
  indicator[rand]=1
  S<-matrix(nrow=354,ncol=442,0)
  for(i in 1:354)
  {
    S[i,rand[i]]<-1
  }
  #selection matrices A and b
  A<-A%*%S
  
  barrier<-function(z1)
  {
    if(all((b-((A%*%z1)))>0.00000000000000000001)==TRUE){
      return(sum(log(1+(1/(b-((A%*%z1)))))))} else{
        return(length(b)*(log(1+(1/0.00000000000000000001))))}
  }
  
  objective<-function(z1,beta1)
  {
    return(((1/(2*(sigma^2)))*((norm_vec((X1%*%beta1)-z1))^2))+barrier(z1))
  }
  
  grad<-function(beta)
  {
    value=nlm(objective,p=y,beta1=beta,hessian=TRUE)$estimate
    return(((-beta/tau^2)+((t(X1)%*%(y-value))/(sigma^2))))
  }
  
  weight<-function(beta)
  {
    nlm(objective,p=y,beta1=beta,hessian=TRUE)$minimum
  }
  
  neg_adjs_posterior_barrier<-function(beta)
  {
    pr=(-(norm_vec(beta)^2))/(2*(tau^2))
    pr<-pr+weight(beta)
    logpost<-((-(norm_vec(y-(X1%*%beta)))^2)/(2*(sigma^2)))+pr
    return(-logpost)
  }
  
  fit <- nlm(neg_adjs_posterior_barrier,coef[E],hessian = TRUE)
  MAP<-fit$estimate
  #langevin random walk:
  #step size
  gamma=1/l
  g.samp=mvrnorm(n=2000,mu=rep(0,l),Sigma=diag(1,l))
  chain<-matrix(nrow=2000,ncol=l,0)
  chain[1,]<-coef[E]
  
  for(i in 1:1999)
  {
    chain[i+1,]=chain[i,]+(gamma*grad(chain[i,]))+((sqrt(2*gamma))*g.samp[i,])
    print(i)
    print(chain[i,])
  }
  
  chain0=chain
  chain.bn=chain[-(1:1000),]
  bci=matrix(nrow=l,ncol=2)
  for(i in 1:l)
  {
    bci[i, ] <- quantile(chain.bn[, i], probs = c(0.025, 0.975))
  }
  
  return(cbind(bci, MAP))
}

test.adjusted.randomized<-function(X, y)
{
  #estimating sigma^2 from full regression model
  LM = lm(y~X)
  sigma = summary(LM)$sigma
  tau = 1000
  n = 442
  p = 10
  D = diag(1,p)
  
  rand<-sort(sample(c(1:n),size=354, replace = FALSE, prob = NULL))
  y.sel<-y[rand]
  X.sel<-X[rand,]
  n.sel =354
  
  epsilon = mvrnorm(n=5000,mu=rep(0,n),Sigma=((sigma^2)*diag(1,n)))
  sum = 0
  for(i in 1: 2000)
  {
    sum = sum +  sigma*max(t(X)%*%epsilon[i,])
  }
  lam = sum/2000
  out = genlasso(y.sel, X=X.sel, D=D)
  c<-coef(out,lambda=lam)
  coef<-as.vector(c$beta)
  E<-which(abs(coef)>0.001)
  l=length(E)
  z=sign(coef[E])

  adj_inference = sampler(X, X.sel, y, E, z, l, lam, alpha, sigma)
  return(adj_inference)
}