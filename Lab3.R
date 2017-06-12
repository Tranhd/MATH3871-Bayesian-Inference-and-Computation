
# 1 -----------------------------------------------------------------------
# F(x) propto x^2e^-x, 0<=x<=1. Use rejection sampling.
# Maximum at x = 1 of e^-1 ≈ 2.7 
sample <- function(N){
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

out=sample(5000)
hist(out,probability=T)

xx=seq(0,1,0.001)
lines(xx,dgamma(xx,3,1)/pgamma(1,3,1),col=2) #True density (red line)
lines(c(0,1),c(1,1),col=3) #sampling distribution (green line)


# 2 -----------------------------------------------------------------------
# Approximate following integrals

#a)
u = runif(10000)
approxa <- mean(exp(exp(u)))
print(approxa)

# b)
approxb <- 2*mean(1/u^2*exp(-(1/u-1)^2))
print(approxb)

#c) 
u1 = runif(10000)
u2 = runif(10000)
approxc <- mean(exp((u1+u2)^2))
print(approxc)

#d) 


z1=1/u1-1
z2=1/u2-1
approxd <- mean(1/u1^2/u2^2*exp(-(z1+z2))*(u1<u2))
print(approxd) 
      

# 3 -----------------------------------------------------------------------
# Monte Carlo sample variance assesment.

N = seq(1000,20000, 1000) # Sample sizes
vvec = c()
estimatevec = c()
k = 1
for(i in N){
  u = runif(i)
  hu = (cos(50*u) + sin(20*u))^2
  hm = mean(hu)
  vm = (1/i^2)*sum((hu-hm)^2)
  estimatevec[k] = hm
  vvec[k] = vm
  k = k+1
}

plot(N,estimatevec, type="l", ylim=c(0.9,1.1),xlab="m",ylab="Estimate",main="")
lines(N,estimatevec+sqrt(vvec), lty=2,col=2)
lines(N,estimatevec-sqrt(vvec), lty=2,col=2)
abline(h=0.965,col=3)


# 4 -----------------------------------------------------------------------
# Same as 3) but using importance sampling.

N = seq(1000,20000, 1000) # Sample sizes
vvec = c()
estimatevec = c()
k = 1
for(i in N){
  u = rbeta(i,0.1,2)
  hu = (cos(50*u) + sin(20*u))^2
  hm = mean(hu/dbeta(u,0.1,2))
  vm = (1/i^2)*sum((hu/dbeta(u,0.1,2)-hm)^2)
  estimatevec[k] = hm
  vvec[k] = vm
  k = k+1
}

plot(N,estimatevec, type="l", ylim=c(0.5,1.5),xlab="m",ylab="Estimate",main="")
lines(N,estimatevec+sqrt(vvec), lty=2,col=2)
lines(N,estimatevec-sqrt(vvec), lty=2,col=2)
abline(h=0.965,col=3)

# This sampling distribution seems to produce a higher estimate variability than
# the uniform distribution. Why is this?

# The beta generats nearly no values higher than 0.2, almost only around 0.

