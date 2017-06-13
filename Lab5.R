# 1 -----------------------------------------------------------------------
# Generate n = 50 samples from the N(μ,σ2) distribution with μ = 5 and σ = 2. 
# Implement a Gibbs sampler to draw samples from π(μ,σ−2|x).

mu = 5
sigma = 2
n = 50
kappa = 1
phi =5
alpha = 2
Beta =2
X = rnorm(n,mu,sigma)

mucond <- function(tau){
  return(rnorm(1,(tau*sum(X) +kappa*phi)/(tau*n+kappa),1/(tau*n*kappa)));
}

taucond <- function(mu){
  return(rgamma(1,alpha + n/2,Beta + sum((X-mu)^2))/2);
}


Gibbssample <- function(N,mu0,tau0){
  mysamples = c()
  tausamples = c()
  musamples = c()
  mui = mu0
  taui = tau0
  for(i in 1:N){
    taui = taucond(mui)
    mui = mucond(taui)
    tausamples[i] = taui
    musamples[i] = mui
  }
  mysamples[[1]] = tausamples
  mysamples[[2]] = musamples
  return(mysamples);
}

mysamples = Gibbssample(100,rnorm(phi,kappa^-1),rgamma(alpha,Beta))
tausamples = mysamples[[1]]
musamples = mysamples[[2]]
k = seq(1,length(tausamples))
par(mfrow=c(2,1))
plot(k,tausamples,'l')
plot(k,musamples,'l')
par(mfrow = c(1,1))
plot(musamples,tausamples,'l')


# 2 -----------------------------------------------------------------------
rm(list=ls())

# a)  
# Use a Gibbs sampler to compute the posterior distributions for β and φ based on 
# the three univariate conditional distributions. 
n=50
Betaprime = t(c(1,1))
phi = 1
p = 0.95
n0 = 4
S0 = 0.5
my_x = t(rep(0,n))
cov_matrix = matrix(p,n,n)
diag(cov_matrix) <- 1
x = rnorm(n,my_x,cov_matrix)
y = Betaprime[[1]] + Betaprime[[2]]*x + rnorm(n,0,1/phi)
X <- matrix(c(t(rep(1,n)),x),nrow=length(x))

xbar = mean(x)
xsquare = sum(x^2)
Sxx = sum((x-xbar)^2)
n1 = n + n0
Betahat = solve(t(X) %*% X) %*% t(X) %*% y

Beta0cond <- function(Beta0hat,Beta1hat,Beta1,xbar,phi,n,Sxx,xsquare){
  mean = Beta0hat - xbar*(Beta1-Beta1hat)
  var = (1/(phi*n*Sxx))*(1-(n*xbar^2)/(xsquare))
  return(rnorm(1,mean,var));
}

Beta1cond <- function(Beta0hat,Beta1hat,Beta0,xbar,phi,n,Sxx,xsquare){
  mean = Beta1hat - (xbar*n/xsquare)*(Beta0-Beta0hat)
  var = (1/(phi*Sxx))*(1-(n*xbar^2)/(xsquare))
  return(rnorm(1,mean,var));
}

phicond<- function(n1,Beta,Betahat,X,Y,n0,s0){
  Se = t(Y) %*% Y - t(X%*%Beta) %*% X %*% Beta
  Sb = t(Beta-Betahat) %*% t(X) %*% X %*% (Beta-Betahat) + Se + n0*S0
  return(rgamma(1,n1/2,Sb/2));
}


mysamples = c()
beta0samples = c()
beta1samples = c()
phisamples = c()
Beta0 = 1
Beta1 = 1
N=100
beta0samples[1] = Beta0
beta1samples[1] = Beta1
phisamples[1] = phi
for(i in 2:N){
  Beta0 = Beta0cond(Betahat[[1]],Betahat[[2]],Beta1,xbar,phi,n,Sxx,xsquare)
  Beta1 = Beta1cond(Betahat[[1]],Betahat[[2]],Beta0,xbar,phi,n,Sxx,xsquare)
  phi = phicond(n1,c(Beta0,Beta1),Betahat,X,y,n0,s0)
  beta0samples[i] = Beta0
  beta1samples[i] = Beta1
  phisamples[i] = phi
}
 
# Plot the time series plots of the samples obtained. 
k = seq(1,length(beta0samples))
plot(k,beta0samples,'l')
plot(k,beta1samples,'l')
plot(k,phisamples,'l')

# Plot the autocorrelation plot of the samples obtained. 
# Suggest a number k, such that keeping every k-th observation would give an 
# approximately independent sample.
k2 = seq(1,length(beta0samples),1)
plot(acf(beta0samples[k2]))
plot(acf(beta1samples[k2]))
plot(acf(phisamples[k2]))


# c -----------------------------------------------------------------------
# Using the reparameterisation discussed in Lecture 5, 
# run the Gibbs sampler again with the new parameterisation. 

alphacond <- function(phi,sigmaalpha,Ax,Betahat){
  mean = Ax %*% Betahat
  Sigma = sigmaalpha*(1/phi)
  return(mvrnorm(n=1,mean,Sigma));
}

mysamples = c()
alpha0samples = c()
alpha1samples = c()
phisamples = c()
beta0samples_para = c()
beta1samples_para = c()

Beta0 = 1
Beta1 = 1
phi = 1
N=100
Beta = matrix(data = c(Beta0,Beta1),nrow=2,ncol=1)
Ax = matrix(data=c(1,0,xbar,1),nrow=2,ncol=2)
sigmaalpha = matrix(data=c(1/n,0,0,1/Sxx),nrow=2,ncol=2)
temp = Ax %*% Beta
alpha0 = temp[1]
alpha1 = temp[2]
alpha0samples[1] = alpha0
alpha1samples[1] = alpha1
phisamples[1] = phi
beta0samples_para[1] = Beta0
beta1samples_para[1] = Beta1
for(i in 2:N){
  temp = t(t((c(alpha0 - xbar*alpha1,alpha1))))
  phi = phicond(n1,temp,Betahat,X,y,n0,s0)
  alphasamp = alphacond(phi, sigmaalpha, Ax, Betahat)
  alpha0 = alphasamp[1]
  alpha1 = alphasamp[2]
  alpha0samples[i] = alpha0
  alpha1samples[i] = alpha1
  phisamples[i] = phi
  beta0samples_para[i] = alpha0 - xbar*alpha1
  beta1samples_para[i] = alpha1
}

k = seq(1,length(beta0samples_para))
plot(acf(beta0samples))
plot(acf(beta0samples_para))