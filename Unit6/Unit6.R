###################### conjugate priors##############################
alpha = 5
beta = 6
theta = seq(0,100, length.out =101)
# d-use to density function
# set shape and scale
prior = dgamma(theta, shape = alpha, scale = beta)
x = 42
posterior = dgamma(theta, shape = x+alpha, scale = 1/(1+1/beta) ) 


# make plot
plot(theta, posterior, xlab = expression(theta), ylab = "density")
lines(theta, prior, lty = 4)

# posterior draws
postdraw = rgamma(2000, shape = x+alpha, scale = 1/(1+1/beta) ) 
r1= hist(postdraw, freq = F, breaks = 20, plot =F)
lines(r1, lty=3, freq = F, col = "gray90")

################## a beta-binomial model ###################

# function to compute the logarithm of the posterior density
betabinexch0=function (theta, data)
{
  eta = theta[1]
  K = theta[2]
  y = data[, 1]
  n = data[, 2]
  N = length(y)
  logf = function(y, n, K, eta) lbeta(K * eta + y, K * (1 -
                                                          eta) + n - y) - lbeta(K * eta, K * (1 - eta))
  val = sum(logf(y, n, K, eta))
  val = val - 2 * log(1 + K) - log(eta) - log(1 - eta)
  return(val)
}
# read the data
cancermortality = read.table("cancer.txt",header = T)
# contour plot after install LearnBayes package
install.packages("LearnBayes")
library(LearnBayes)
mycontour(betabinexch0,c(.0001,.003,1,20000),cancermortality, xlab=expression(eta),ylab="K")


# function to compute the logarithm of the posterior density with the tranformed parameters
betabinexch=function (theta, data)
{
  eta = exp(theta[1])/(1 + exp(theta[1]))
  K = exp(theta[2])
  y = data[, 1]
  n = data[, 2]
  N = length(y)
  logf = function(y, n, K, eta) lbeta(K * eta + y, K * (1 -
                                                          eta) + n - y) - lbeta(K * eta, K * (1 - eta))
  val = sum(logf(y, n, K, eta))
  val = val + theta[2] - 2 * log(1 + exp(theta[2]))
  return(val)
}

# contour plot
mycontour(betabinexch,c(-8,-4.5,3,16.5),cancermortality, xlab="logit eta",ylab="log K")

# calculate pi
##################### Rejection sampling
## Example 7
x = NULL
n = 100
r = 2
for(i in 1:(100*sqrt(2*pi/exp(1)))) # 100*c means 100 times accepted, becasue 100c * (1/c) equals to 100!
{
  
  u = runif(1)
  y = rcauchy(1)
  if( u <=  exp(-y^2/2)*(1+y^2)*sqrt(exp(1))/2) 
  { x = c(x,y)} #combine x and y
}

# investigate the acceptance rate
N= 1000

c = sqrt(2*pi/exp(1))

# you can identify C using the following code
X = seq(-3,3,0.01)
c = max(dnorm(X)/dcauchy(X))
# this c value is not the real c, but it is very close to the true value
u = runif(N)
y=rcauchy(N)
accept = ( u < dnorm(y)/(c*dcauchy(y)))
x = y[accept]
mean(accept)
1/c

par(mfrow=c(1,2))
hist(x)
hist(rnorm(100))


#################a stationary Markov chain
P<-rbind(c(0.5,0.5,0,0,0),c(0.25,0.5,0.25,0,0),c(0,0.25,0.5,0.25,0),
         c(0,0,0.25,0.5,0.25),c(0,0,0, 0.5,0.5))

# s be that a vector of locations that the person visit
s<-vector('numeric',50000)
s[1]<-1 # initial location
for(j in 2:50000)
  {s[j]<-sample(1:5, size=1,prob=P[s[j-1],])}

# summarize the frequencies f visits to these five states after 500, 1000, 2000, 5000, 10000, 20000, 50000
v<-c(500,1000,2000,5000,10000,20000,50000)
for(i in 1:length(v))
  print(table(s[1:v[i]])/v[i])

#verify the stationary
w<-c(0.125,0.25,0.25,0.25,0.125)
w%*%P


#########################################
#### Gibss samplers for bivariate normal
rho = 0.2
x0 = rnorm(1)
N= 10000
X = x0
Y = NULL
## for t=1
Y = rnorm(1,rho*x0,sqrt(1-rho^2))
X = rnorm(1,rho*Y[t],sqrt(1-rho^2))
## the main loop of generating samples using gibbs sampler
for(t in 2:N)
{
  Y = c(Y,rnorm(1,rho*X[t-1],sqrt(1-rho^2)))
  X = c(X,rnorm(1,rho*Y[t],sqrt(1-rho^2)))
}

## compare with the true standard normal 
par(mfrow=c(2,2))
X0 = seq(-100,100,by = 0.1)
hist(X[1000:N],breaks = 100, freq = FALSE,xlab="X")
dX0 = dnorm(X0)
lines(X0,dX0)
hist(Y[1000:N],breaks = 100, freq = FALSE,xlab="Y")
lines(X0,dX0)

## correlation 
cor(cbind(X,Y))

## compare with the samples generated using mvrnorm
library(MASS)
mu = c(0,0)
Sigma = rbind(c(1,rho),c(rho,1))
U = mvrnorm(N,mu,Sigma)
hist(U[,1],breaks = 100, freq = FALSE,xlab="X")
dX0 = dnorm(X0)
lines(X0,dX0)
hist(U[,2],breaks = 100, freq = FALSE,xlab="Y")
lines(X0,dX0)
## correlation 
cor(U)


#### hierarchical model, binomial-beta
install.packages("VGAM")
library(VGAM)
N=10000  
n= 30
a=3
b=7
X= vector("numeric",N)  
theta = vector("numeric",N)
## for t=1
theta[1]=rbeta(1,a,b)  
X[1]=rbinom(1,n,theta[1])
for (i in 2:N)
{ 
  X[i]=rbinom(1,n,theta[i-1])
  theta[i]=rbeta(1,a+X[i],n-X[i]+b)
}
## compare with the true marginal distribution 
par(mfrow=c(1,2))
X0 = seq(0,30,by = 1)
hist(X[1000:N], freq = FALSE,xlab="X")
lines(X0,dbetabinom.ab(X0,n,a,b))
theta0 = seq(0,1,by=0.01)
hist(theta[1000:N],breaks = 100, freq = FALSE,xlab=expression(theta))
lines(theta0,dbeta(theta0,a,b))


#######An example of Gibbs sampling########################
y = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)
t = c(94, 16, 63, 126, 5, 31, 1, 1, 2, 10)


#### function to draw lambda
lambda.draw = function(alpha, beta, y, t) 
{
  rgamma(length(y), y + alpha, t + beta)
}

#### function to draw beta
beta.draw = function(alpha, gamma, delta, lambda, y) 
{
  rgamma(1, length(y) * alpha + gamma, delta + sum(lambda))
}

##### Gibbs sampling 
gibbs = function(n.sims, beta.start, alpha, gamma, delta, y, t, burnin = 0, thin = 1) 
{
  beta.draws = c()
  
  lambda.draws = matrix(NA, nrow = n.sims, ncol = length(y))
  
  beta.cur = beta.start
  
  for (i in 1:n.sims) {
    # draw lambda
    lambda.cur = lambda.draw(alpha = alpha, beta = beta.cur,y = y, t = t)
    
    # draw beta
    beta.cur = beta.draw(alpha = alpha, gamma = gamma, delta = delta, lambda = lambda.cur, y = y)
    
    # burn-in and thinning
    if (i > burnin & (i - burnin)%%thin == 0) {
      lambda.draws[(i - burnin)/thin, ] = lambda.cur
      beta.draws[(i - burnin)/thin] = beta.cur
    }
  }
  
  return(list(lambda.draws = lambda.draws, beta.draws = beta.draws))
}


#### Apply to the failure data set
posterior = gibbs(n.sims = 10000, beta.start = 1, alpha = 1.5,
                  gamma = 0.01, delta = 1, y = y, t = t)
round(colMeans(posterior$lambda.draws),3)
mean(posterior$beta.draws)
round(apply(posterior$lambda.draws, 2, sd),3)
round(sd(posterior$beta.draws),3)



