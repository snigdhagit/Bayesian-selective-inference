#function that returns unadjusted intervals: X is predictor matrix, y is reponse,
#E is selected set (by Lasso), tau is prior variance
unadjusted_intervals<-function(X, y, E, tau)
{
  X1<-X[,E]
  l=length(E)
  Q<-(((tau^2)*X1%*%t(X1))+diag(n))
  mu0<-(tau^2)*(t(X1)%*%solve(Q)%*%y)
  var<-(tau^2*diag(l))-(tau^4)*(t(X1)%*%solve(Q)%*%(X1))
  UBC<-NULL
  LBC<-NULL
  for(i in 1:l)
  {
    LBC[i]<-mu0[i]-1.96*(sqrt(var[i,i]))
    UBC[i]<-mu0[i]+1.96*(sqrt(var[i,i]))
  }
  #storing unadjusted credible intervals in cred
  cred<-cbind(LBC,UBC)
  return(cred<-as.matrix(cred))
}

###(non-randomized): computing adjusted Bayesian intervals using a Langevin sampler and solving approximate posterior MAP
sampler<- function(X, y, E, z, l, lam, sigma)
{ 
  n = length(y)
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
  
  #estimating adjusted posterior MAP
  fit <- nlm(neg_adjs_posterior_barrier,coef[E],hessian = TRUE)
  MAP<-fit$estimate
  
  grad<-function(beta)
  {
    value=nlm(objective,p=y,beta1=beta,hessian=TRUE)$estimate
    return(((-beta/tau^2)+((t(X1)%*%(y-value))/(sigma^2))))
  }
  
  #implementing a Langevin walk
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

test.adjusted.nonrandomized<-function(X, y)
{
  #estimating sigma^2 from full regression model
  LM = lm(y~X)
  sigma = summary(LM)$sigma
  tau = 1000
  n = 442
  p = 10
  D = diag(1,p)
  
  
  library("genlasso")
  out = genlasso(y, X, D=D)
  
  epsilon = mvrnorm(n=5000,mu=rep(0,n),Sigma=((sigma^2)*diag(1,n)))
  sum = 0
  for(i in 1: 2000)
  {
    sum = sum +  sigma*max(t(X)%*%epsilon[i,])
  }
  lam = sum/2000
  
  c<-coef(out,lambda=lam)
  coef<-as.vector(c$beta)
  E<-which(abs(coef)>0.001)
  l=length(E)
  z=sign(coef[E])
  adj_inference = sampler(X, y, E, z, l, lam, sigma)
  unadj_inference = unadjusted_intervals(X, y, E, tau)
  return(list(adjusted.inf = adj_inference, unadjusted.inf = unadj_inference))
}

#running adjusted inference on diabetes data
install.packages("care")
library("care")
data(efron2004)
dim(efron2004$x)
X= efron2004$x
#normalize X
X = apply(X,2,function(x){return (x/norm_vec(x))})
y = efron2004$y
test.adjusted.nonrandomized(X, y)