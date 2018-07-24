####################################################
# The following code is for section 5.1
####################################################

# simulate random variables
set.seed(999)
n<-20
rnorm(n,10,3) # normal distribution
runif(n,2,3) # uniform distribution between 2 and 3
rbeta(n,2,5) # beta distribution 

hist(rnorm(1000,10,4))
hist(runif(1000,2,3))
hist(rbeta(1000,2,5))

####################################################
# The following code is for Section 5.2
####################################################

# population and sample
set.seed(1000)
N<-10000
X<-rnorm(N,10,3)
mean(X)    # the population mean
var(X)    # the population variance
n<-500
x<-sample(X,n,replace=F)
mean(x) # the sample  mean
var(x) # the sample variance

# sample mean and sample variance
fsample<-function(mu,sigma,n,nsim)
{ xbar<-rep(0,nsim);xvar<-rep(0,nsim)
for(i in 1:nsim)
{set.seed(i)
  # x contains a random sample of size n of the variable X
  x<-rnorm(n,mu,sigma)
  xbar[i]<-mean(x)
  xvar[i]<-var(x)
}
cat('sample mean:',mean(xbar),"\n")
cat('sample variance:',mean(xvar),"\n")
# plot the samples
par(mfrow=c(2,2))  
hist(xbar);qqnorm(xbar);qqline(xbar);
hist(xvar);qqnorm(xvar);qqline(xvar);
}

# sample mean and sample variance of samples of size 10 
# from a normal distribution with mean 2 and sigma 2
fsample(2,2,10,10)

# increase simulation times to 100, 1000, 10000
fsample(2,2,10,100)
fsample(2,2,10,1000)
fsample(2,2,10,10000)

# change sample size
mu<-2;sigma<-2;n<-10;nsim<-1000
xbar<-rep(0,nsim);xbar2<-rep(0,nsim);xbar3<-rep(0,nsim)
xvar<-rep(0,nsim);xvar2<-rep(0,nsim);xvar3<-rep(0,nsim)
for(i in 1:nsim)
{set.seed(i)
  x<-rnorm(n,mu,sigma);x2<-rnorm(n*10,mu,sigma)
  x3<-rnorm(n*100,mu,sigma)
  xbar[i]<-mean(x);xbar2[i]<-mean(x2);xbar3[i]<-mean(x3)
  xvar[i]<-var(x);xvar2[i]<-var(x2);xvar3[i]<-var(x3)}
c(mean(xbar),mean(xbar2),mean(xbar3))
c(var(xbar),var(xbar2),var(xbar3))
c(mean(xvar),mean(xvar2),mean(xvar3))

# histogram and QQ-plot for xbar
par(mfcol=c(3,2))
hist(xbar);hist(xbar2);hist(xbar3)
qqnorm(xbar);qqline(xbar)
qqnorm(xbar2);qqline(xbar2)
qqnorm(xbar3);qqline(xbar3)


# histogram and QQ-plot for xvar
par(mfcol=c(3,2))
hist(xvar);hist(xvar2);hist(xvar3)
qqnorm(xvar);qqline(xvar)
qqnorm(xvar2);qqline(xvar2)
qqnorm(xvar3);qqline(xvar3)

# 95% confidence interval
alpha<-0.05
c(mean(x)-qt(1-alpha/2,n-1)*sd(x),mean(x)+qt(1-alpha/2,n-1)*sd(x))


# now repeat the above simulation for beta distribuiton 
alpha<-2;beta<-4;n<-10;nsim<-1000
xbar<-rep(NA,nsim);xbar2<-rep(NA,nsim)
xbar3<-rep(NA,nsim)
xvar<-rep(NA,nsim);xvar2<-rep(NA,nsim)
xvar3<-rep(NA,nsim)
for(i in 1:nsim)
{set.seed(i)
  x<-rbeta(n,alpha,beta);x2<-rbeta(n*10,alpha,beta)
  x3<-rbeta(n*100,alpha,beta)
  xbar[i]<-mean(x);xbar2[i]<-mean(x2);xbar3[i]<-mean(x3)
  xvar[i]<-var(x);xvar2[i]<-var(x2);xvar3[i]<-var(x3)}
c(mean(xbar),mean(xbar2),mean(xbar3))
c(var(xbar),var(xbar2),var(xbar3))
c(mean(xvar),mean(xvar2),mean(xvar3))

# histogram and QQ-plot for xbar
par(mfcol=c(3,2))
hist(xbar);hist(xbar2);hist(xbar3)
qqnorm(xbar);qqline(xbar)
qqnorm(xbar2);qqline(xbar2)
qqnorm(xbar3);qqline(xbar3)

# histogram and QQ-plot for xvar
par(mfcol=c(3,2))
hist(xvar);hist(xvar2);hist(xvar3)
qqnorm(xvar);qqline(xvar)
qqnorm(xvar2);qqline(xvar2)
qqnorm(xvar3);qqline(xvar3)



# simulated value of the mean squared error of an estimator
# try it for n=10,50, 100  
fMSE<-function(mu,sigma,n,nsim)
{xbar<-rep(NA,nsim)
for(i in 1:nsim)
{set.seed(i)
  # x contains a random sample of size n of the variable X
  x<-rnorm(n,mu,sigma);xbar[i]<-mean(x)};mse<-mean((xbar-mu)^2)
  xvar<-var(xbar);bias<-mean(xbar)-mu
  cat('case n=',n, ': mse=', mse, ', var=', xvar, ', bias=', bias,  "\n")}
fMSE(2,2,10,1000)
fMSE(2,2,50,1000)
fMSE(2,2,100,1000)


################### A case study
# function to generate data
generateData <- function(n,p,beta)
# p = number of estimators
{
  X <- matrix(runif(n*p),n,p)
  # true model does not have x0
  y <- rnorm(n,X%*%beta)
  return(list(X=X,y=y))
}

# best subset selection given data
bestsubset <- function(Data)
{
  require(leaps)
  fit = regsubsets(y~X,Data)
  b = numeric(ncol(Data$X)+1)
  names(b)  =  fit$xnames
  bb = coef(fit,which.min(summary(fit)[["cp"]]))
  b[names(bb)] = bb
  return(b)
}

# ridge regression
ridge <- function(Data,n.lam=501)
{
  require(MASS)
  lam  =  c(0,exp(seq(0,20,len=n.lam)))
  fit = lm.ridge(y~X,Data,lambda=lam)
  b = coef(fit)[which.min(fit$GCV),]
  return(b)
}

# lasso variable selection
lasso <- function(Data)
{
  require(glmnet)
  # can also use Lar
  cvfit = cv.glmnet(Data$X,Data$y)
  b = as.numeric(coef(cvfit,s=cvfit$lambda.min))
  return(b)
}

# the main body of simulation
N = 1000
n = 50
beta = c(3,-3,rep(0,8))
p = length(beta)
out = array(NA,dim=c(N,p,4),
            dimnames=list(1:N,1:p,c("Subset","Lasso","Ridge","OLS")))
for (i in 1:N)
{
  Data = generateData(n,p,beta)
  out[i,,1] = bestsubset(Data)[-1]
  # [-1]exclud beta zero
  out[i,,2] = lasso(Data)[-1]
  out[i,,3] = ridge(Data)[-1]
  out[i,,4] = coef(lm(y~X,Data))[-1]
}
bias <- apply(out,2:3,mean)-beta
variance <- apply(out,2:3,var)
MSE <- bias^2+variance
apply(MSE,2,sum)

par(mfrow=c(1,2))
ylim = range(out)
boxplot(out[,,1],col="blue",ylim=ylim,main="Subset")
boxplot(out[,,2],col="blue",ylim=ylim,main="Lasso")




