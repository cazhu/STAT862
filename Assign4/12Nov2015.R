norm = function (n, alpha) 
{
        vec = vector("numeric", n)
        vec[1] = 0
        for (i in 2:n) {
                y = runif(1, -alpha, alpha )
                y = y + vec[i-1]
                aprob = min(1, dnorm(y)/dnorm(vec[i-1]))
                u = runif(1)
                if (u < aprob) 
                    vec[i] = y
                else 
                    vec[i] = vec[i-1]
        }
        return(vec)
}
normvec<-norm(10000,1)
par(mfrow=c(2,1))
plot(ts(normvec))
hist(normvec,30)


### independence chain, generate gamma from normal with the same mean and variance 
cauchy = function (n, alpha, beta) 
{
   mu = alpha/beta
   sigma = sqrt(alpha/(beta^2))
   vec = vector("numeric", n)
   vec[1] = alpha/beta
   for (i in 2:n) {
	# need to change the proposal 
      y <- rnorm(1, , 2)
	# change dgamma to dcauchy 
      aprob <- min(1, (dgamma(y, alpha, beta)/dgamma(vec[i-1], 
        alpha, beta))/(dnorm(y, mu, sigma)/dnorm(vec[i-1], 
          mu, sigma)))
      u <- runif(1)
      if (u < aprob) 
          vec[i] = y
      else 
          vec[i] = vec[i-1]
    }
    return(vec)
}
vec<-gamm(10000,2,4)
par(mfrow=c(2,1))
plot(ts(vec))
hist(vec,20)

###q() cancel out 
# target dist and proposal distribution and parameters
### exp 13
# f = cauchy --> target dist, 
# g = y ~ N(theta^(t-1), 2) 
# burn in --> remove the first 200 samples 
# generate the fisrt item in the sequence, other values can be choosen, they evenetually converge to the target sample
# vec[i-1] is t-1 in course notes



################indepence chain#######

cauch = function (n) 
{
   vec = vector("numeric", n)
   vec[1] = 0
   for (i in 2:n) {
      y <- rnorm(1, 0, 2)
      aprob <- min(1, (dcauchy(y)/dcauchy(vec[i-1]))/(dnorm(y, 0, 2)/dnorm(vec[i-1],0, 2)))
      u <- runif(1)
      if (u < aprob) 
          vec[i] = y
      else 
          vec[i] = vec[i-1]
    }
    return(vec)
}

vec<-cauch(10000)
par(mfrow=c(2,1))
plot(ts(vec))
hist(vec,20)


gamm = function (n, alpha, beta) 
{
   mu = alpha/beta
   sigma = sqrt(alpha/(beta^2))
   vec = vector("numeric", n)
   vec[1] = alpha/beta
   for (i in 2:n) {
      y <- rnorm(1, mu, sigma)
      aprob <- min(1, (dgamma(y, alpha, beta)/dgamma(vec[i-1], 
        alpha, beta))/(dnorm(y, mu, sigma)/dnorm(vec[i-1], 
          mu, sigma)))
      u <- runif(1)
      if (u < aprob) 
          vec[i] = y
      else 
          vec[i] = vec[i-1]
    }
    return(vec)
}
vec<-gamm(10000,2,4)
par(mfrow=c(2,1))
plot(ts(vec))
hist(vec,20)



