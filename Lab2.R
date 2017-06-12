# Part 1 --------

# Using the inverse transform method, write an R function to generate a random variable with the distribution function 1/2(x^2 + x)



rm(list=ls())

Inversesamp <- function(N) { 
  u = runif(N)
  return(sqrt(2*u+(1/4)) - 1/2);
  }

x = Inversesamp(10000)
hist(x, probability = T)
a = seq(0,1,0.01)
y = a+0.5
lines(a,y, pch = 7,col="red")

# Part 2 ------------
# Using the inverse transform method, write an R function to generate a random variable with the distribution function
# f(x) = exp(−2x) 0≤x<∞ or exp(2x) −∞<x<0 

Inversesamp2 <- function(N){
  u <- runif(N)
  return(c(log(2*u[u < 0.5])/2, -log(2*(1-u[u>0.5]))/2));
}

N <- 100000
hist(Inversesamp2(N),probability = T, breaks = 90)
a = seq(-5,5,0.01)
y = c(exp(2*a[a<0]),exp(-2*a[a>=0]))
lines(a,y,pch=7, col="red")

# Part 3 ------------
#Using the probability integral transform method, write an R function to generate n random samples from an 
#Exponential distribution

rm(list=ls())


Expsample <- function(lambda, N){
  u = runif(N)
  return(-log(u)/lambda)
}

N <- 10000
lambda <- 0.005
x <- Expsample(lambda,N)
hist(x, probability = T)
a = seq(0,qexp(0.99,lambda),0.1)
d = max(a)
y = lambda*exp(-lambda*a) 
lines(a,y, pch = 7,col="red")
# Do not need N to be that large, just > 100.

# Part 4 ---------------------
#If x1,...,xn are samples from an Exp(λ) distribution, 
#then  sum(xi) ∼ Gamma(n, λ). Using this result, and your answer from the previous question, 
#simulate samples from a Gamma(n, λ) distribution. 
#Produce a histogram of samples drawn from this distribution and superimpose the density function.

gammasample <- c()
lambda <- 5
N <- 1000
n <- 100
sample <- c()
for(i in 1:N){
  gammasample[i] = sum(Expsample(lambda,n))
}

hist(gammasample,probability = T)
xx = seq(qgamma(0.001,n,lambda),qgamma(0.99,n,lambda),0.01)
y = dgamma(xx, n, rate=lambda)
lines(xx,y,pch = 7,col="red")

