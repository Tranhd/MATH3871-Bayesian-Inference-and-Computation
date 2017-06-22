
# 1 -----------------------------------------------------------------------
# Data from N(theta,1)
# Model 1: theta = 0
# Model 0: theta ≠ 0, f(theta) = N(1,1)
# Priors f(0) = f(1) = 1/2

# a) ---- Estimate bayes facto with 1000 samples.
rm(list = ls())
xbar = 1
n = 10
N = 1000
theta = rnorm(N,1,1)
m1 = dnorm(xbar,0,1/sqrt(n)) #Probabilities of data given theta = 0
m0 = mean(dnorm(xbar,theta,1/sqrt(n))) # Probability of data given theta ≠ 0 (Monte carlo)
B01 = 1/(sqrt(n+1)*exp(-0.5*((n*xbar+1)^2/(n+1)-1))) # Correct bayes (B01) factor
cat(m0/m1) # Estimated B01
cat(B01) # True B01

# b) ----- 
# Determine of the estimate of B01 is affected by number of monte carlo samples, N.
# Produce 1000 replicate Bayes factor estimates when N = 500,1000,5000.
# Produce histograms and superimpose the true B01
temp1=temp2=temp3=rep(NA,1000) # To store samples
N=5000 # Largest number of N
m1=dnorm(xbar,0,1/sqrt(n)) #Probabilities of data given theta = 0
for (i in 1:1000) { #1000 replicate Bayes factor estimates
  theta=rnorm(N,1,1) # 5000 Prior sample
  temp=dnorm(xbar,theta,1/sqrt(n)) #Probabilities of data given theta ≠ 0 (Monte carlo)
  temp1[i]=mean(temp[1:500]) # Sample size 500
  temp2[i]=mean(temp[1:1000]) # Sample size 1000
  temp3[i]=mean(temp) # Sample size 5000
}

par(mfrow=c(1,3)) #1x3 subplots
hist(temp1/m1,probability=T,main="N=500",xlab="Bayes Factor") 
points(B01,0) #500 sample size histogram of B01
hist(temp2/m1,probability=T,main="N=1000",xlab="Bayes Factor", 
xlim=range(temp1/m1)) #1000 sample size histogram of B01
points(B01,0)
hist(temp3/m1,probability=T,main="N=5000",xlab="Bayes Factor", 
xlim=range(temp1/m1)) #5000 sample size histogram of B01
points(B01,0) # Superimpose true value.

# Obviously more focused around true value for N = 5000.

# c) ---------
# Estimate the Bayes Factor based on an importance sampling distribution N(m,a^2) for model 0,
# where m = 1 and a = 0.6/√11,1/√11,2/√11 and 5/√11.
# Again, obtain 1,000 replicate estimates of each of these and produce a histogram of 
# the estimate in each case, superimposing the true value of B01.

true.mean=1 # True mean
true.sd=1/sqrt(10) #true sd
N=1000 
a=c(0.6,1,2,5) # Different values of a.
B=matrix(NA,ncol=4,nrow=1000) # 1000x4 matrix
for (j in 1:4) { # For every a value
  for (i in 1:1000) { # For every replicate.
    theta=rnorm(N,true.mean,a[j]*true.sd)  # Sample from prior
    w=dnorm(xbar,theta,1/sqrt(10))*dnorm(theta,1,1)/dnorm(theta,true.mean,a[j]*true.sd) 
    #weights as true m0/sampling distribution.
    B[i,j]=mean(w)/m1 # B01 (Mean because monte carlo)
  }} 
par(mfrow=c(2,2))
r=range(B)
for (j in 1:4) {
  hist(B[,j],xlim=r,main=paste("a =",round(a[j]/sqrt(11),3)))
  points(B01,0)
}

# d) ------
# Explain the differences in the quality of the estimate of B01 for different values of a.
# When a = 1/sqrt(10) g(θ) is in fact the true posterior distribution. Hence all
# estimates of B01 will be identical as the Monte Carlo estimator will have variance
# 0 (and so the histogram plot will only have one bar, and look slightly odd). Of
# course we will never have this in practice.

# a = 2/sqrt(10) and a = 5/sqrt(10) represent more inefficient importance distributions
# that are worse approximations to the posterior. Of course, the worse the importance 
# sampling distribution, the more variable the Monte Carlo estimate of the
# marginal likelihood. Though the estimate should in be unbiased in this case.

# For a = 0.6/ 11 the importance sampling distribution does not adequately
# cover the target posterior, and so the Monte Carlo estimates are very poor. 
# The tail of these estimates goes up to 100, in this instance – 
# this does not happen for the other estimators.

# e) --------
# Using the importance distribution g(θ), which is centered at the same value as the
# prior and posterior, is of course more efficient when the variance is closer to that of
# the posterior (i.e. when 1/sqrt(10) ≤ a < 1).



# 2 -----------------------------------------------------------------------
#Repeat Question 1a) and 1b) with your choice of sample sizes, N, 
#but instead compute unbiased Monte Carlo estimates of the reciprocal Bayes Factor B01^-1 = B10
B10 = B01^-1 # True Bayes factor
xbar = 1 # data
n=10
par(mfrow=c(1,1))

N = 5000 # Number of samples for monte carlo
theta = rnorm(N,1,1/sqrt(n+1)) #
m1 = dnorm(xbar,0,1/sqrt(n))
rm0 = mean(1/dnorm(xbar,theta,1/sqrt(n))) #1/m0
B10estimate = m1*rm0
cat(B10estimate)
cat(B10)

rstore=rep(NA,1000)
N=25000
m1=dnorm(xbar,0,1/sqrt(n))
for (i in 1:1000) {
  theta=rnorm(N,1,1/sqrt(n+1))
  rstore[i]=mean(1/dnorm(xbar,theta,1/sqrt(n)))
}
hist(rstore*m1,probability=T,main=paste("N =",N),xlab="Bayes Factor")
points(B10,0)

# 3 -----------------------------------------------------------------------
# The Unknown Variance example in Lecture 7. Suppose that we have observations x1, . . . , xn
# which are iid N(θ,φ).

# model m = 0: The variance is known as φ = φ0, but θ is unknown.
# model m = 1: Both mean and variance are unknown.

# Assume that the prior for model 1 is determined by π1(θ|φ) ∼ N(m,φ/d) and 
# π1(φ) ∝ φ−t−1 exp(−bφ−1), and that the prior for model 0 by 
# π0(θ) = π1(θ|φ = φ0) ∼ N(m,φ0/d).

# a) choose (θ0,φ0) and generate data.
rm(list = ls())

theta0 = 1
phi0 = 1
n = 200
x = rnorm(n,theta0,phi0)

# b) For appropriate choices of prior parameters m,d,b and t compute the true Bayes
# Factor B01.
m = 1
d = 1
t = 2
b = 2
y = sum((x-mean(x))^2) + (n*d*(mean(x)-m)^2)/(n+d)
B01 = ((b+y/2)^(t+n/2))/(b^t*phi0^(n/2)) * gamma(t)/gamma(t + n/2) * exp(-y/(2*phi0)) # True bayes factor
m0_true = sqrt(d/(n+d))/(2*pi*phi0)^(n/2) * exp(-y/(2*phi0))
m1_true = sqrt(d/(n+d))/(2*pi)^(n/2) * gamma(t+n/2)/(gamma(t)*(b+y/2)^(t+n/2)) * b^t
## ----
# c) Estimate the Bayes Factor by Monte Carlo integration based on sampling from the
# priors.
theta0 = rnorm(n,m,phi0/d) # Prior for theta for model 0
m0samp = c()
for(i in 1:length(theta0)){
  m0samp[i] = prod(dnorm(x,theta0[i],phi0)) #Estimate the likelihood using every value for theta
}
m0 = mean(m0samp) # Take mean to monte carlo.


# Get priors from phi through importance sampling using uniform(0,10)
y = runif(10000,0,10) # Sample distribution
w = (y^(-t-1)*exp(-b*y^-1)) # Compute weights
w = w/sum(w)
phi1 = sample(y, n, replace=T, prob=w)
hist(phi1,breaks=15,density=T)
q = seq(0,10,length=100)
lines(q,(q^(-t-1)*exp(-b*q^-1)),col=2)
m1samp = c()
m1samp2 = c()

for(i in 1:length(x)){
  theta1temp = mean(rnorm(m,phi1/d))
  m1samp[i] = prod(dnorm(x,theta1temp,phi1[i]))
}

m1 = mean(m1samp)
cat(m0)
cat(m0_true)
cat(m1)
cat(m1_true)
cat(m0/m1)
cat(m0_true/m1_true)
cat(B01)



