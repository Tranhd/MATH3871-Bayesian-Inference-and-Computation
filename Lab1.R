# Section 1a ------
# Compute a point estimate (i.e. the posterior mean) and a 95% central credible interval for θ.

#clear workspace
rm(list=ls())

# Gamma mean function
Gammamean <- function(aa, bb) { Gammam <- aa/bb; return(Gammam); }

# Numbers given
n <- 65
sum <- 24890
alpha <- 1
Beta <- 0.01

#Update beta distribution parameters
alpha <- alpha+sum
Beta <- Beta+n

#Plot posteriour gamma-distribution using density-gamma-function
theta = seq(0.005, 600, length = 500)
posterior = dgamma(theta, alpha, Beta)
plot(theta, posterior, type = "l")

#Calculate and print posterior mean
posterior_mean = Gammamean(alpha, Beta)
print(paste("The posterior mean is:",posterior_mean))

#Calculate and print credible-intervall of posterior mean
CR<-qgamma(c(0.025,0.975), alpha, Beta)
print(paste("The credible-intervall is: [",CR[[1]],",",CR[[2]],"]"))


# Section 1b----
# Monte Carlo estimates of the lower and upper values of the 95% credible interval for θ.

#Function to estimate credible intervall using quantiles.
CrEstimate <- function(N,alpha,Beta) { x<-rgamma(N, alpha, Beta); 
return(quantile(x, probs = c(0.025, 0.975)));} 

#Function to generate multiple credible intervall estimates
CrSample <- function(N,alpha,Beta,rep){
  x <- CrEstimate(N,alpha,Beta) #Draw sample
  mysamp = c(x[[1]],x[[2]]) 
  
  for(i in 2:rep){ #Draw additonal rep - 1 samples
    x <- CrEstimate(N,alpha,Beta)
    mysamp <- append(mysamp, x, after = length(mysamp))
  }
  return(mysamp)
}


N <- 5000 #Gamma sample size
rep <- 250 #Number of repititions 

mysamp = CrSample(N,alpha,Beta,rep)

hist(mysamp,breaks=80) #Histogram of cr-bounderies
points(CR,c(5,5),pch = 19,col="red") #Red points where correct cr-bounderies are

# For N = 500, Monte carlo estiamte seems fairly accurate
# N = 5000, spread of Monte carlo estimate decreases, much more accurate.

# Section 1b part 2 -----
# How many samples, N, are needed for the spread (i.e. min - max) of the Monte Carlo 
# estimates for each interval endpoint to be less than 0.15

#Function to find the spread for lower and higher intervall point.
CrSample_withspread <- function(N,alpha,Beta,rep){
  x <- CrEstimate(N,alpha,Beta) #Draw sample
  loweri = c(x[[1]])
  higheri = c(x[[2]])
  
  for(i in 2:rep){ #Draw additonal rep - 1 samples
    x <- CrEstimate(N,alpha,Beta)
    loweri <- append(loweri, x[[1]], after = length(loweri)) #Append lower intervallpoint
    higheri <- append(higheri, x[[2]], after = length(higheri)) #Append higher intervallpoint
  }
  return(c(max(loweri)-min(loweri),max(higheri)-min(higheri))) #Return the maximum spread for both low and high.
}

rep <- 250 #Number of repititions 
maxN <-5000 #Maxinum number samples to try
from <-500 #Minimum number samples to try
plot.new() #New plotwindow
frame()
spreads = c()
Nsamples<-seq(from,maxN,100)
for(i in Nsamples){
  spread = CrSample_withspread(i,alpha,Beta,rep) #New sample
  spreadL = spread[[1]]
  spreadH = spread[[2]]
  spreads <- append(spreads, max(spread), after = length(spreads)) #Append largest spread
}
plot(Nsamples,spreads) #Plot number of samples vs. Maximum spread.
#Just to illustrate effect of number of samples.

#Need approximatly N=50000 samples to get spread < 0.15 

# Section 1c -----
# Draw samples directly from the posterior distribution and compare to Monte Carlo.

N <- 5000
sample <- rgamma(N, alpha, Beta)
y = rpois(N,sample)
hist(y,probability = T)
x = seq(300, 500, 1)
negativebin = dnbinom(x, alpha,1-1/(Beta+1))
lines(x,negativebin,pch = 7,col="red")

# Section 1d -----

## Algebraically exact approach:
#Pros: The answers will always be exact, with no simulation error.
#Cons: It is often tedious to perform the calculations, and it is frequently not possible to perform them in closed form.
## Simulation:
#Pros: It is usually trivial to obtain estimates of whatever posterior quantities you need. 
#Cons: They are only estimates, and come with estimation/simulation error (that can be controlled).

# Section 2a ------
#  Produce marginal posterior histograms of each parameter, and scatterplots of the 6 bivariate distributions.

test = read.table("tuberculosis.txt") #Read table
names(test) <- c("alpha","c","p","mu") #Add colum names


par(mfrow=c(2,2)) #2x2 plot-grid

#Plot the histograms
hist(test[["alpha"]])
hist(test[["c"]])
hist(test[["p"]])
hist(test[["mu"]])

# bivariate -----
#New plot window.
plot.new()
frame()
par(mfrow=c(3,2)) #3x2 plot-grid

#Scatter-plot bivariate distributions
plot(test[["alpha"]],test[["c"]],pty=5,cex=0.1) 
plot(test[["alpha"]],test[["p"]],pty=5,cex=0.1) 
plot(test[["alpha"]],test[["mu"]],pty=5,cex=0.1) 
plot(test[["c"]],test[["p"]],pty=5,cex=0.1) 
plot(test[["c"]],test[["mu"]],pty=5,cex=0.1) 
plot(test[["p"]],test[["mu"]],pty=5,cex=0.1) 


# Section 2b ------ 
# Produce a histogram of the posterior distribution of Φ.
# is there any evidence that Φ < 1

par(mfrow=c(1,1)) #Standard grid

#Constant
delta <- 0.52
tau <-0.52
epss<-0.52
epsr<-0.202

#The relative fitness of the drug-resistant strains
theta <- (1-test[["c"]])*(1/tau + 1/(delta + epsr))/(1/tau + 1/(delta+epss+test[["p"]]))

# Histogram of the posterior distribution of theta
hist(theta, probability = T)

thetaless <- which(theta<1)
Plessthan1 <- length(thetaless)/length(theta) #The probability that theta is < 1.
cat(Plessthan1)
# Equal to 2%, significant evidence that the resistant strain is more evolutionarily fit than the susceptible strain

# Section 2c -----
# A central 95% credible interval
CR <- quantile(test[["p"]],c(0.025,0.975)) #95% central Credible intervall for p.
len <- CR[[2]]-CR[[1]] #Length of interval is ≈ 0.127
cat(len) #Large because posterior distribution vey skewid.

# Section 2d -------
# HDR regions

# Function to generate random p% credible intervalls
GenerateCredibleintervall <- function(p,data){
  x <- runif(1,min=-(1-p)/2,max=(1-p)/2)
  lower <- (1-p)/2 + x
  higher <- p + (1-p)/2 + x
  return(quantile((data),c(lower,higher)))
}

shortestCR = c() 
Lengthofshortest = Inf
for(i in 1:250){
  CR <- GenerateCredibleintervall(0.95,test[["p"]]) #new cr
  lens <- CR[[2]]-CR[[1]] #Length of cr
  if(lens < Lengthofshortest){ #if shorter -> update shorter.
    Lengthofshortest <- lens
    shortestCR <- CR
  }
}
print(paste("Shortest cr: [",shortestCR[[1]], ",",shortestCR[[2]],"]"))
print(paste("Length of shortest cr:",Lengthofshortest))
print(paste("Length of central cr:", len))
#Big difference.
#Likely to chose the central one as the shortest when the skewidness is minimum.
#The shortest one removes the top 5%, i.e it is between 0-95%.

# Section 3a-----
#Function to estimate pi from the buffon needle experiment.
rm(list = ls())
BuffonNeedle <- function(ll,dd,nn){ 
  l <- ll #Length of needle
  d <- dd #Distanec between lines
  n <- nn #Number of samples.
  hit <- 0 #Number of hits
  for(i in 1:n){
    dis <- runif(1,0,d/2) #Distance from the middle of the needle to nearest line
    theta <- runif(1,0,pi) #The angle from the vertical axis to the needle
    if(dis <= (l/2)*sin(theta)){ #If the half the needle projected on the horizontal axis
      # is greater than the distance to the nearest line, we hit.
      hit <- hit +1
    }
  }
  phat <- hit/n #Point estimate of probability of hitting line
  return((2*l)/(phat*d)) #Estimation of pi
}

lengthh<-0.1
dist<-1
nmax <- 2000 
nmin <- 100
nvec <- seq(nmin,nmax,10)
piestimates = c()
#Run experiment for different number of samples.
for(i in nvec){
  piestimates<-append(piestimates,BuffonNeedle(lengthh,dist,i),length(piestimates))
}
pis = array(1,dim=length(nvec))*pi #To see correct pi value
plot(nvec,piestimates,type="l") #The estimation of pi vs sample size
matplot(nvec,pis,add = T,type = "l", col = "red") #Red line representing pi

# Section 3c ------
# Buffons needle.

n <- 100
dist<-1
lengthh<-seq(0.2,1,0.1)
ratio = lengthh/dist 
N<-1000
variances= c()
for(i in 1:length(lengthh)){
  piestimates = c()
  for(j in 1:N){
    piestimates<-append(piestimates,(BuffonNeedle(lengthh[[i]],dist,n)),length(piestimates))
  }
  variances<-append(variances,var(piestimates),length(variances))
}

plot(ratio,variances, type="l")

