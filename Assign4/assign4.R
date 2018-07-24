

############Q1###########

###(a)#########
# independence chains #

bino = function (nsim, n, p) 
{
  vec = vector("numeric", nsim)
  # how to choose vec[1]
  vec[1] = 1
  for (i in 2:nsim) {
    y <- round (runif(1, min = 0, max = n), digits = 0)
    aprob <- min(1, (dbinom(y, n, p)/dbinom(vec[i-1], n, p))/(dunif(y, min = 0, max = n)/dunif(vec[i-1], min = 0, max = n)))
    u <- runif(1)
    if (u < aprob) 
      vec[i] = y
    else 
      vec[i] = vec[i-1]
  }
  return(vec)
  
}

vec<-bino(10000, 20, 0.7)
x0<-rbinom(10000, 20, 0.7)


c(mean(x0), mean(vec))
c(var(x0), var(vec))


par(mfrow=c(2,2))

plot(ts(vec))
hist(vec, xlim = c(0, 20), main = "Estimated")

plot(ts(x0))
hist(x0, xlim = c(0, 20), main = "Theoretical")

#######(b)#############


bino2 = function (n, sigma)
{
  vec = vector("numeric", n)
  vec[1] = 0
  for (i in 2:n) {
    y = rnorm(1, vec[i], sigma )
    y = y + vec[i-1]
    aprob = min(1, dnorm(y) / dnorm(vec[i-1]))
    u = runif(1)
    if (u < aprob) 
      vec[i] = y
    else 
      vec[i] = vec[i-1]
  }
  return(vec)
}
set.seed(998)
vec2<-bino2(10000, 0.1)
vec3<-bino2(10000, 0.5)
vec4<-bino2(10000, 10)
x0<-rnorm(10000)

rbind(cbind(mean(x0), var(x0)),
cbind(mean(vec2), var(vec2)),
cbind(mean(vec3), var(vec3)),
cbind(mean(vec4), var(vec4))
)

cbind(var(x0), var(vec2),var(vec3), var(vec4))

par(mfrow=c(4,2))

plot(ts(x0))
hist(x0, xlim = c(-4,4), main = "Theoretical")

plot(ts(vec2))
hist(vec2, xlim = c(-4,4), main = "Estimated, var = 0.01")

plot(ts(vec3))
hist(vec3, xlim = c(-4,4), main = "Theoretical, var = 0.25")

plot(ts(vec4))
hist(vec4, xlim = c(-4,4), main = "Theoretical, var = 100")

#############
## Q2
set.seed(1234)
install.packages("VGAM")
library(VGAM)


N=10000  
a=5
b=1
X= vector("numeric",N)  
theta = vector("numeric",N)
# for t = 1
theta[1]=rgamma(1,a,b)  
X[1]=rpois(1,theta[1])
# for t = N
for (i in 2:N)
{ 
  X[i]=rpois(1,theta[i-1])
  theta[i]=rgamma(1,a+X[i],1+1/b)
}
## compare with the true marginal distribution 
#par(mfrow=c(1,2))
#X0 = seq(0,6,by = 0.1)
#hist(X[1000:N], freq = FALSE,xlab="X")
#lines(X0,dbetabinom.ab(X0,n,a,b))
theta0 = seq(0,20,by=0.01)
hist(theta[1000:N],breaks = 100, freq = FALSE,xlab=expression(theta), main = "")
lines(theta0,dgamma(theta0,a,b))

