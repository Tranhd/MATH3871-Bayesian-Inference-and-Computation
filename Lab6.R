
# 1 -----------------------------------------------------------------------
## Gibs sampler

joint_density = function(x,y){x^2*exp(-x*y^2-y^2+2*y-4*x)}

# x|y ~ Gamma(3,1/(y^2 + 4))
# y|x ~ N(1/(x+1),sqrt(1/(2*(x+1))

N=1000 # Sample size
X = Y = rep(NA,N) # Initiate data arrays
X[1] = Y[1] = 0.5 # Initiate with x=y= 0.5

# Gibbs sampler:
for(i in 2:N){
  X[i] = rgamma(1,shape=3,scale=1/(Y[i-1]^2 + 4)) # Sample from conditional x|y
  Y[i] = rnorm(1,1/(X[i]+1),sqrt(1/(2*(X[i]+1)))) # Sample from conditional y|x
}

# Inspect the marginal trace plots and distributions for x and y
par(mfrow=c(2,2)) # 2x2 plot grid
hist(X, probability=T) # Distribution of X
plot(X, type="l") # Trace plot for X
hist(Y, probability=T) # Distribution of Y
plot(Y, type="l") # Trace plot for Y 

# Scatterplot of the bivariate samples with superimposed contours of the joint distribution, to check
# validity of joint samples.
xx=seq(0,3,length=100) # x between 0 and 3
yy=seq(-3,3,length=100) # y between -3 and 3
zz=matrix(NA,ncol=100,nrow=100) # "Value" matrix z[i,j] = f(x[i],y[j])
for (i in 1:100) zz[,i]=joint_density(xx,yy[i]) # Fill matrix with values. 
par(mfrow=c(1,1)) # New plot window.
contour(xx,yy,zz,xlab="x",ylab="y",col=2) # Draw contour
points(X,Y,pch=16,cex=0.3) # Superimpose joint samples.

# Looks fine, like it has converged.


# 1b ----------------------------------------------------------------------
## Construct an independence Metropolis-Hasting sampler using a Bivariate proposal distribution.
# Independance =  
# Where the proposal is the same regardless of the state of the chain: q(θ∗|θt) = q(θ∗).

ffn=function(x, y) { x^2*exp(-x*y^2-y^2+2*y-4*x) } #pi
N=3000 # Number of iterations
X=Y=rep(NA, N) # Initiate data arrays.
X[1]=Y[1]=0.5 # Initial values
for (i in c(2:N)) { #Independence Metropolis hastings sampler using bivariate proposal distribution.
  propx=rnorm(1, 0, 0.5) # Generate proposals for x,y
  propy=rnorm(1, 0, 0.5)
  # Calculate log acceptance probability for moving from X[i-1],Y[i-1] to propx,propy.
  acc=log(ffn(propx,propy))-log(ffn(X[i-1], Y[i-1])) # pi(z*) P(new state) - pi(z) P(current state)
  +dnorm(X[i-1],0,0.5,log=T)+dnorm(Y[i-1],0,0.5,log=T) 
  -dnorm(propx, 0, 0.5,log=T)-dnorm(propy,0,0.5,log=T) #q(z*,z) - q(z,z*) (P(moving from z* to z) - P(z to z*))
  if (log(runif(1))< acc && propx>0 ) { #Accept new proposal with probability acc, if x in defmängd.
    X[i]=propx; Y[i]=propy 
  } else {
    X[i]=X[i-1]; Y[i]=Y[i-1] # Else stay in previous state.
  }
}
par(mfrow=c(2,2))
hist(X, probability=T)
plot(X, type="l")
hist(Y, probability=T)
plot(Y, type="l")
# Scatterplot of the bivariate samples with superimposed contours of the joint distribution, to check
# validity of joint samples.
xx=seq(0,3,length=100) # x between 0 and 3
yy=seq(-3,3,length=100) # y between -3 and 3
zz=matrix(NA,ncol=100,nrow=100) # "Value" matrix z[i,j] = f(x[i],y[j])
for (i in 1:100) zz[,i]=joint_density(xx,yy[i]) # Fill matrix with values. 
par(mfrow=c(1,1)) # New plot window.
contour(xx,yy,zz,xlab="x",ylab="y",col=2) # Draw contour
points(X,Y,pch=16,cex=0.3) # Superimpose joint samples


# 1c ----------------------------------------------------------------------
## Construct a Random-Walk Metropolis-Hastings sampler which updates one parame- ter at a time. 

# Investigate how changing the scale of the random walk variance affects the performance of the sampler.
# Trying to scale by; 10,1,0.1
sc = 0.1
ffn=function(x, y) { x^2*exp(-x*y^2-y^2+2*y-4*x) }
N=3000 # Number of iterations
X=Y=rep(NA, N) # Initiate data arrays.
X[1]=Y[1]=0.5 # Initial values
for (i in c(2:N)) {
  propx=rnorm(1, X[i-1], sc*0.5) # Generate proposal for x using previous value (Random walk)
  acc=ffn(propx, Y[i-1])/ffn(X[i-1], Y[i-1]) # Calculate acceptance probability for propx
  # Equal to pi(proposal)/pi(previous) since proposal symmetric.
  if (runif(1)< acc && propx>0 ) { X[i]<-propx # Accept with probablity acc
  } else { X[i]<-X[i-1] }
  propy=rnorm(1, Y[i-1], sc*0.5) # Generate proposal for y using previous value (Random walk)
  acc=ffn(X[i], propy)/ffn(X[i], Y[i-1]) # Calculate acceptance probability for propy
  # Equal to pi(proposal)/pi(previous) since proposal symmetric.
  if (runif(1)< acc) { Y[i]<-propy # Accept with probablity acc
  } else { Y[i]<-Y[i-1] }
}
par(mfrow=c(1,1))
plot(X, type="l")

# sc = 10
# Chain that moves infrequently, but is able to make large jumps when it does move.
# Logical since we can move large distances but will be wrong most of the time and not accept.
# The performance of chain is undesirable.

# sc = 1
# Moves in a good way

# sc = 0.1
# Conversely the chain in the bottom panel moves frequently, but only makes small jumps. 
# Logical since were only move short distances.
# The performance of chain is undesirable.


# 2 -----------------------------------------------------------------------
## This problem considers a simple example where Gibbs sampling is not straightforward. 
# Suppose y = (y1, . . . , ym) are a random sample from the Weibull(ρ, κ) distribution.
# Suppose we place independent Gamma priors on ρ ∼ Gamma(α, β) and κ ∼ Gamma(γ, δ).
# We can write expressions for the conditional and joint distribution

# a) For your own choice of m, κ (shape) and ρ (scale, note R uses an inverse parameterisation),
# simulate m observations from the Weibull(κ, ρ) distribution.
m=100
kappa0=2
rho0=2
y=rweibull(m,kappa0,1/rho0)

# It is clear that the posterior conditionals are not standard distributions, 
# and hence the standard Gibbs sampler cannot be used here. 
# Devise an appropriate Metropolis- Hastings algorithm to approximate the posterior distribution. 

N=2500 # Nummber of iterations.
alpha=beta=gam=delta=2  # prior parameters.
logpost=function(rho,kappa,y,m,alpha,beta,gam,delta) { # The log-posteriour distribution, pi.
  m*log(kappa) +(m*kappa)*log(rho)+(kappa-1)*log(prod(y))-rho^kappa*sum(y^kappa)
  + dgamma(rho,alpha,1/beta,log=T) + dgamma(kappa,gam,1/delta,log=T)} 

kappa=rho=rep(NA, N) # Initiate kappa and rho.
kappa[1]=rho[1]=0.5
for (i in c(2:N)) { # Run metropolis hastings.
  propkappa=rnorm(1, kappa[i-1], 0.25) # Generate random wakl proposal for kappa.
  proprho=rnorm(1, rho[i-1], 0.25) # Generate random walk proposal for rho.
  # pi(proposal)-pi(current).
  acc=logpost(proprho,propkappa,y,m,alpha,beta,gam,delta) - logpost(rho[i-1],kappa[i-1],y,m,alpha,beta,gam,delta) 
  if (log(runif(1))< acc && propkappa>0 && proprho>0) { # Accept with probability acc.
    kappa[i]=propkappa; rho[i]=proprho
  } else { kappa[i]=kappa[i-1]; rho[i]=rho[i-1] }
}


# Compare the results obtained to the true values you used to simulate the data.
par(mfrow=c(2,2))
plot(kappa,type='l',xlab="Iteration",ylab="kappa")
lines(c(1,N),rep(kappa0,2),col=2)
plot(rho,type='l',xlab="Iteration",ylab="rho")
lines(c(1,N),rep(rho0,2),col=2)
hist(kappa,probability=T,xlab="kappa")
points(kappa0,0,pch=16)
hist(rho,probability=T,xlab="rho")
points(rho0,0,pch=16)