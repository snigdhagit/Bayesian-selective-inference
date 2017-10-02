#computes approximate UMVUE under additive Gaussian randomization with variance tau^2 (in the univariate case)
UMVUE.compute<-function(y, sigma, tau)
{
  objective<-function(z,alpha)
  {
    (((-alpha*z))+((z^2/2)+log(1+(1/z))))
  }
  
  objective1<-function(z,alpha)
  {
    (((-alpha*z)/tau)+((z^2/2)+log(1+(1/z))))
  }
  
  approx1<-function(alpha1)
  {
    opt<-nlm(objective1,p=1,alpha=alpha1,hessian=TRUE)
    return(opt$estimate)
  }
  
  UMVUE.approx<-(y*(1+((sigma^2)/(tau^2))))-(((sigma^2)/(tau))*approx1(y))
  return(UMVUE.approx)
}
