
# 1 -----------------------------------------------------------------------

# F(x) propto x^2e^-x, 0<=x<=1. 
# Rejection sampling.
# Maximum at x = 1 of e^-1 ≈ 2.7 


# a)
samp <- function(N){
  out = c()
  K = exp(-1)
  while(length(out) < N){
    x = runif(N) # Sampling distribution
    y = runif(N,0,K) 
    index = (y < x^2*exp(-x)) # Accept?
    out=c(out,x[index])
  }
  return(out);
}

out=samp(5000)


# Now implemented through importance sampling.
N = 5000
x = runif(N) # Sample distribution
w = x^2*exp(-x) # Compute weights
w = w/sum(w) # Normalize the weights.

# b)
# Compare density estimates and estimates of the distributional means 
# between importance and rejection sampling.

# Can obtain unweighted samples by sampling repeatedly from the empirical distribution (w,theta)
x.noweight=sample(x, 5000, replace=T, prob=w)
hist(x.noweight,probability=T,xlab="x",ylab="Density",main="")
xx=seq(0,1,length=100)
lines(xx,dgamma(xx,3,1)/pgamma(1,3,1),col=2)
lines(c(0,1),c(1,1),col=3)

cat("Rejection sampling mean:",mean(out),"\n")
cat("Importance sampling mean:",sum(w*x),"\n")
cat("Importance sampling mean (unweighted):",mean(x.noweight),"\n")


# 1c ----------------------------------------------------------------------
# Sample from f(x) using a different importance sampling distribution, 
# and compare the efficiency of this to the samples obtained in (a) through the effective sample size. 
# Which sampling distribution was more efficient, and why?
# Use another sampling distribution. Compare with a) 
a=3
b=1
y = rbeta(5000,a,b) # Another sampling distribution.
w2 = y^2*exp(-y)/dbeta(y,a,b)
w2 = w2/sum(w2)
y.noweight = sample(y, 5000, replace=T, prob = w2)
hist(y.noweight,probability=T,xlab="y",ylab="Density",main="",xlim=c(0,1))
xx=seq(0,1,length=100)
lines(xx,dgamma(xx,3,1)/pgamma(1,3,1),col=2)
lines(c(0,1),c(1,1),col=3)
lines(xx,dbeta(xx,a,b),col=4,lty=2)
cat("ESS (uniform g) =",1/sum(w^2),"\n")
cat("ESS (beta g) =",1/sum(w2^2),"\n")
# clearly a more efficient estimator under the Beta(3,1) prior.

# 3 -----------------------------------------------------------------------
# Sequential importance sampling framework.

# a)
# Using importance sampling with g(x) given by U(0,1), draw 1000 weighted samples
# from f(x). What is the effective sample size?

rm(list=ls())
par(mfrow=c(1,1))

xx=seq(0,1,length=150)
plot(xx,dnorm(xx,0.5,0.01),type='l',xlab="u",ylab="density")
u=runif(1000)
w=dnorm(u,0.5,0.01)
ess=function(w) { W=w/sum(w); return(1/sum(W^2))}
cat("The ESS is",ess(w),"\n")
# 32.65/1000 = 3.2% efficiency.

# Many of the importance weights will be small, 
# but those that are drawn in the region of high density will be large, 
# and so the variance of the weights will be large. 
# Accordingly we will have a low effective sample size in this setting. 

# b) 
# How many samples do you need before your effective sample size is larger than 5000?
xx=seq(0,1,length=150)
plot(xx,dnorm(xx,0.5,0.01),type='l',xlab="u",ylab="density")
u=runif(150000)
w=dnorm(u,0.5,0.01)
ess=function(w) { W=w/sum(w); return(1/sum(W^2))}
cat("The ESS is",ess(w),"\n")
# Trial and error -> ≈ 150 000 samples.

# c)
# Try sequential importance sampling.

# Sample from f(x)^(p) using uniform sample distribution.
p=0.5
u=runif(1000)
w=dnorm(u,0.5,0.01)^p
cat("The ESS is",ess(w),"\n") # ESS ≈ 48
xx=seq(0,1,length=850)

a=sample(u,10000,replace=T,prob=w/sum(w))
hist(a,xlab="u",ylab="density",probability=T,ylim=c(0,50),
     main="Sequential importance sampling")
lines(xx,dnorm(xx,0.5,0.01),type='l')
lines(xx,dnorm(xx,0.5,0.01)^p,lty=2,col=2)

# Construct a new importance sampling distribution that is a normal approximation of the 
# weighted samples from f(x)^(1/2) and sample from  f(x) using this.
m.a=mean(a)
sd.a=sd(a)
x=rnorm(1000,m.a,sd.a)
w2=dnorm(x,0.5,0.01)/dnorm(x,m.a,sd.a)
cat("The ESS is",ess(w2),"\n") # ESS ≈ 904

# d)
# How many samples are needed (including the initial 1000) to obtain samples from
# f(x) with an effect sample size of more than 5000?
N=5000/ess(w2)*1000 * 1.05# Guess how many needed and add 5%
x=rnorm(N,m.a,sd.a)
w3=dnorm(x,0.5,0.01)/dnorm(x,m.a,sd.a)
cat("The ESS is",ess(w3),"(N =",N+1000,")\n")

# e)
# Which procedure is more efficient?
# It takes approximately 6842.88  samples to get an ESS of 5,000 under the two stage procedure, 
# compared to about 150 000 samples in a single stage procedure.

# Sequential importance sampling aims to break down potentially 
# inefficient importance sampling stages into a sequence of more 
# efficient ones in order to control the variability of the importance weights.
