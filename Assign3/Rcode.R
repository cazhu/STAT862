###################
# STAT 462/862 - Assignment 3
# Chunyang Zhu
###################
#Q1


library(MASS)
library(xtable)
# read data
admit<-read.table('admit.txt',header=T)



rank<-as.factor(rank) 



z <- matrix(, nrow=400, ncol = 4)

for(j in 1:4){

for(i in 1:400){
  if(admit$rank[i] == j){z[i,j] = 1}
  else {z[i,j] = 0}
}
  
}

# random partition
set.seed(42);
admit.new <- cbind(admit[,4], admit[,1:2], z)
colnames(admit.new) <- (cbind("admit", "gre","gpa","z1","z2", "z3", "z4"))
admit.random <- sample(admit.new)
attach(admit.random)


admit.training <- admit.random[1:200,]
admit.test <- admit.random[201:400,]

# logistic regression

model_logit = glm(admit ~ gre + gpa + z2 + z3 + z4,  family = binomial(link='logit'), data=admit.training)
summary(model_logit)
xtable(model_logit)

# prediction 
pred_logit= predict(model_logit,admit.test, type="response")
as.numeric(pred_logit>=0.5)

summary(pred_logit)

###### linear discriminant analaysis
model_lda = lda(admit ~ gre + gpa + z2 + z3 + z4, data = admit.training)
model_lda
summary(model_lda)

# prediction
pred_lda = predict(model_lda, admit.test)
summary(pred_lda)

# quadratic discriminant analysis
model_qda = qda(admit ~ gre + gpa + z2 + z3 + z4, data = admit.training)
model_qda
# prediction
pred_qda = predict(model_qda,admit.test)
summary(pred_qda)

# comparison
pred_logit_factor <- as.factor(as.numeric(pred_logit > 0.5))
data.frame(pred_lda$class, pred_qda$class, pred_logit_factor, admit.test[,7])
table(pred_lda$class, admit.test$admit[1:200])
table(pred_qda$class, admit.test$admit[1:200])
table(pred_logit_factor, admit.test$admit[1:200])
#############

# Q2(a)

# sample mean and sample variance
fsample<-function(mu,sigma,n,nsim,alpha)
{ xbar<-rep(0,nsim);xvar<-rep(0,nsim); 

# calculate confidence interval

xpred_lower <- rep(0, nsim)
xpred_upper <- rep(0, nsim)

xpred_lower_estimate <- rep(0, nsim)
xpred_upper_estimate <- rep(0, nsim)
s <- rep(0, nsim)

for(i in 1:nsim)
{set.seed(i)
  # x contains a random sample of size n of the variable X
  x<-rnorm(n,mu,sigma)
  xbar[i]<-mean(x)
  xvar[i]<-var(x)
  s[i] = sd(x)
  
  xpred_lower[i]<-xbar[i] + qnorm(alpha/2)*sigma / n^0.5
  xpred_upper[i]<-xbar[i] - qnorm(alpha/2)*sigma / n^0.5
  
  xpred_lower_estimate[i] = xbar[i] + qnorm(alpha/2) * s[i] / n^0.5;
  xpred_upper_estimate[i] = xbar[i] - qnorm(alpha/2) * s[i] / n^0.5;
  
}

cat('sample mean:',mean(xbar),"\n")
cat('sample variance:',mean(xvar),"\n")

cat('sample confidence interval: [', mean(xpred_lower), ",", mean(xpred_upper), "] \n")
cat('coverage:', mean( mu >=  xpred_lower & mu <= xpred_upper), "\n" )

cat('estimated sample confidence interval: [', mean(xpred_lower_estimate), ",", mean(xpred_upper_estimate), "] \n")
cat('esimated coverage:', mean( mu >=  xpred_lower_estimate & mu <= xpred_upper_estimate) )

# plot the samples
par(mfrow=c(2,2))  
hist(xbar);qqnorm(xbar);qqline(xbar);
hist(xvar);qqnorm(xvar);qqline(xvar);
}

# sample mean and sample variance of samples of size 10 
# from a normal distribution with mean 2 and sigma 2
fsample(2,5^0.5,10,1000,0.05)

# Q2(b)

fsample(2,5^0.5,10,1000,0.05)
fsample(2,5^0.5,10,1000,0.025)
fsample(2,5^0.5,100,1000,0.05)
fsample(2,5^0.5,100,1000,0.025)

# Q2(c)



#########################

## Q4
# (a) construct an algorithm

fbeta <- function (n, alpha, beta)
{
  x1 <- rgamma(n ,shape = alpha, scale = 1)
  x2 <- rgamma(n ,shape = beta, scale = 1)
  x<- x1/(x1+x2)
  
  x0 <- rbeta(n, alpha, beta, ncp = 0)
 # return(x)
  par(mfrow=c(2,2))  
  hist(x);qqnorm(x);qqline(x);
  hist(x0);qqnorm(x0);qqline(x0);
  
}

fbeta(10000, 2, 2)

# (b) rejection method
# based on uniform distribution
install.packages("truncnorm")
library("truncnorm")

# you can identify C using the following code

alpha = 2
beta = 2

X = seq(0,1,0.01)
c = max(dbeta(X, alpha, beta)/dtruncnorm(X))
cc = max(dbeta(X, alpha, beta)/dunif(X))

x = NULL
n = 10000
xx = NULL

for(i in 1: (10000*c))
{
  
  u = runif(1)
  yy = runif(1)
  if( u <=  dbeta(yy, alpha, beta)/(dunif(yy)*c)) 
  {xx = c(xx,yy)}
}


for(i in 1: (10000*c))
{
  
  u = runif(1)
  y = rtruncnorm(1, a= 0, b= 1, mean = 1, sd = 2)
  if( u <=  dbeta(y, alpha, beta)/(dtruncnorm(y, a = 0, b = 1, mean = 0, sd =1)*c)) 
  { x = c(x,y)}
}

# investigate the acceptance rate
N = 10000
u = runif(N)
y = rtruncnorm(N, a= 0, b= 1, mean = 0, sd = 1)
accept = (dbeta(y, alpha, beta)/(dtruncnorm(y, a = 0, b = 1, mean = 0, sd = 1)*c))
# accept2 = (dbeta(y, alpha, beta)/(dunif(y)*cc))

N = 10000
u = runif(N)
yy = runif(N)
accept2 = (dbeta(yy, alpha, beta)/(dunif(yy)*cc))


x = y[accept]
xx = yy[accept2]
x
mean(accept)
1/c

mean(accept2)
1/cc

par(mfrow=c(2,2))
hist(x);qqnorm(x);qqline(x);
x0 = rbeta(n, alpha, beta)
hist(x0);qqnorm(x0);qqline(x0)

par(mfrow=c(2,2))
hist(xx);qqnorm(xx);qqline(xx);
x0 = rbeta(n, alpha, beta)
hist(x0);qqnorm(x0);qqline(x0)

#######

# Q5

# Q5(a)


# draw random variables from f(x)


integral = function (n5){
x5 = rexp(n5, rate = 0.5)

# evaluate the mean with different n, show the figure of mean vs n
hx = 2 * (sin(x5)^2)*exp(-x5^0.5)
theta = mean(hx)
return (theta)
}

integral(10000)
integral(100000)
integral(1000000)
# Find the best n and show the final results
XX = NULL
YY = NULL
ZZ = NULL
for(i in 2: 10000) {
  YY = cbind(YY, integral(i))
  XX = cbind(XX, i)
  ZZ = cbind(ZZ, integral(i)-integral(i-1))
 }
plot(XX, YY, xlab = "n", ylab = "theta")
plot(XX, ZZ, xlab = "n", ylab = "deltatheta")

integrand <- function(x) {exp(-(x^0.5 + 0.5*x))*(sin(x))^2}
integrate(integrand, lower = 0, upper = Inf)

# Q5(b) importance sampling
install.packages("matrixStats")
library("matrixStats")
install.packages("OOmisc")
library("OOmisc")
imp<-function(M) # M is sample size
  
{
  # the range of x should be from zero to infinity
  M = 200000
  set.seed(9999)
  X1 = rlaplace(M, 1)
  X2 = rcauchy(M, location = 0, scale = 2)
  X3 = rnorm(M)
  X4 = rnorm(M, pi/2)
  
  N <- 200000
  
  components <- sample(1:3,prob=c(0.5,0.3,0.2),size=N,replace=TRUE)
  mus <- c(pi/2,3*pi/2,5*pi/2)
  sds <- sqrt(c(1,1,1))
  
  X4 <- rnorm(n=N,mean=mus[components],sd=sds[components])
  X4 = X4*(X4>0)
  X1 = X1*(X1>0); X2 = X2*(X2>0); X3 = X3*(X3>0)
  #X3 = rtruncnorm(M, a= 0, b= Inf, mean = 0, sd = 1)
  wgx1 = function (x) 
    {
    2*(sin(x))^2 * exp (-x^0.5 - 0.5*x + abs(x))
    }
  wgx2 = function (x) 
  {
    2*pi*(1+(x^2/4))* (sin(x))^2 * exp (-x^0.5 - 0.5*x)
  }
  wgx3 = function (x) 
  {
    (2*pi)^0.5 * exp(x^2/2)* (sin(x))^2 * exp (-x^0.5 - 0.5*x)
  }
  M = 100
  #0.3*exp ((x-(3*pi/2))^2/2) + 0.2* exp ((x-(5*pi/2))^2/2)
  wgx4 = function (x) 
  {
   
      (2*pi)^0.5 * 
      (
         exp ((x-(pi/2))^2/2)  
      )* (sin(x))^2 * exp (-x^0.5 - 0.5*x)
  }
  IS4 = mean((wgx4(X4)), na.rm=TRUE)
  
  IS3 = mean((wgx3(X3)), na.rm=TRUE)
  IS2 = mean((wgx2(X2)), na.rm=TRUE)
  IS1 = mean((wgx2(X1)), na.rm=TRUE)
  
  standev1=sd((wgx1(X1))/sqrt(M), na.rm=TRUE)
  c1 = (c(IS1,standev1,IS1-qnorm(0.975)*standev1,IS1+qnorm(0.975)*standev1))
  
  standev2=sd((wgx2(X2))/sqrt(M), na.rm=TRUE)
  c2 = (c(IS2,standev2,IS2-qnorm(0.975)*standev2,IS2+qnorm(0.975)*standev2))
  
  standev3=sd((wgx3(X3))/sqrt(M), na.rm=TRUE)
  c3 = (c(IS3,standev3,IS3-qnorm(0.975)*standev3,IS3+qnorm(0.975)*standev3))
  
  
  standev4=sd((wgx4(X4))/sqrt(M), na.rm=TRUE)
  c4 = (c(IS4,standev4,IS4-qnorm(0.975)*standev4,IS4+qnorm(0.975)*standev4))
  c4
  return(rbind(c1,c2,c3))
}

imp(100) 
imp(500) 

imp(1000) 
imp(2000)

imp<-function(mu=0,c=0,n=1000000,x=rnorm(n))
{
  IS = mean((mu+x>c)*((2*pi)^0.5 * exp(x^2/2)* (sin(x))^2 * exp (-x^0.5 - 0.5*x)), na.rm=TRUE)
  standev=sd((mu+x>c)*((2*pi)^0.5 * exp(x^2/2)* (sin(x))^2 * exp (-x^0.5 - 0.5*x)))/sqrt(n)
  return(c(mu,IS,standev,IS-qnorm(0.975)*standev,IS+qnorm(0.975)*standev))
}

#fix n and X in advance
n =10000000
X = rnorm(n) 
# case when c = 3
c= 0
imp(0,c,n,X)
imp(1,c,n,X)
imp(2,c,n,X)
imp(3,c,n,X)
imp(3.1,c,n,X)
imp(3.15,c,n,X) # the (near-)optimal value of mu slightly above c
imp(3.2,c,n,X)
imp(4,c,n,X)
