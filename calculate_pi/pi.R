
# calculate pi by Monte Carlo simulation
pi <- function (n) {
  x = runif(n)
  y = runif(n)
  z = x^2 + y^2 - 1
  m = z[z<=0]
  return (4*length(m)/n)
}

# check the effect of sample size
x=NULL
y=NULL
z=NULL
for (i in 1: 10000){
  x[i] = i
  y[i] = pi(i)
  Z[i] = pi
}
plot(x,y, xlab = "sample size", ylab = "pi")
z =
lines(x, 3.14)