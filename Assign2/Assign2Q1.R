
### plot of the data
install.packages("lattice")
library(lattice)


## install package
install.packages("splines")
library(splines)

# adding lines in a panel
myPanel <- function(...)
{
  panel.smooth(...)
  panel.xyplot(...,alpha=.5)
}

# Generate Y

a <- 2
x <- runif(100, min = -a, max = a)
fx <- -2*(x^2)

e <- rnorm(100, mean = 0, sd = 1)

y <- fx + e

xyplot(y~x)

# Second-order polynomial model

panel.smooth <- function(x,y,...)
{
  fit <- lm(y ~ x + I(x^2))
  x0 <- seq(min(x),max(x),len=300)
  y0 = predict(fit,newdata=data.frame(x=x0),se.fit=T)
  #lpolygon(c(x0,rev(x0)),c(y0$fit-2*y0$se.fit,rev(y0$fit+2*y0$se.fit)),col="gray",border=F,...)
  llines(x0,predict(fit,data.frame(x=x0)),col="black")
}
xyplot(y~x, panel = myPanel)


  fit1 <- lm(y ~ x + I(x^2))
  knots=quantile(x,c(1/3,2/3))
  fit2 = lm(y~bs(x,knots=knots,degree=3))
  fit3 = lm(y~ns(x,knots=knots))
  fit4 <- smooth.spline(x, y,df=12)
  x0 <- seq(min(x),max(x),len=300)
  y01 = predict(fit1,newdata=data.frame(x=x0),se.fit=T)
  y02 = predict(fit2,newdata=data.frame(x=x0),se.fit=T)
  y03 = predict(fit3,newdata=data.frame(x=x0),se.fit=T)
  y04 = predict(fit4,newdata=data.frame(x=x0),se.fit=T)
  
  var1 <- (y01$se.fit)^2
  var2 <- (y02$se.fit)^2
  var3 <- (y03$se.fit)^2
  var4 <- (y04$se.fit)^2
plot(x0, var1, type ="l", xlab = "x", ylab = "variance", col = "black", ylim=c(0, 0.2))

lines(x0, var2, xlab = "x", ylab = "variance", col = "blue")

lines(x0, var3, xlab = "x", ylab = "variance", col = "red2")

lines(x0, var5, xlab = "x", ylab = "variance",col = "green")

legend(-1, 0.2, legend=c("Second-order polynomial", "Cubic spline", "Natural cubic spline", "Smoothing spline"),
       col=c("black", "blue", "red2", "green"), lwd=1, box.col = "white", cex=0.5)

y_true <- -2*(x0^2)

bias1 <- y01$fit - y_true
bias2 <- y02$fit - y_true
bias3 <- y03$fit - y_true
bias5 <- yy$fit - y_true

plot(x0, bias1, type ="l", xlab = "x", ylab = "bias", col = "black", ylim=c(-0.5, 1.2))

lines(x0, bias2, xlab = "x", ylab = "bias", col = "blue")

lines(x0, bias3, xlab = "x", ylab = "bias", col = "red2")

lines(x0, bias5, xlab = "x", ylab = "bias",col = "green")

legend(-1, 1.2, legend=c("Second-order polynomial", "Cubic spline", "Natural cubic spline", "Smoothing spline"),
       col=c("black", "blue", "red2", "green"), lwd=1, box.col = "white", cex=0.5)


plot(x0, bias1^2, type ="l", xlab = "x", ylab = "bias", col = "black", ylim=c(-0.5, 1.2))

lines(x0, bias2^2, xlab = "x", ylab = "bias", col = "blue")

lines(x0, bias3^2, xlab = "x", ylab = "bias", col = "red2")

lines(x0, bias5^2, xlab = "x", ylab = "bias",col = "green")

legend(-1, 1.2, legend=c("Second-order polynomial", "Cubic spline", "Natural cubic spline", "Smoothing spline"),
       col=c("black", "blue", "red2", "green"), lwd=1)



MSE1 <- bias1^2 + var1
MSE2 <- bias2^2 + var2
MSE3 <- bias3^2 + var3
MSE5 <- bias5^2 + var5


plot(x0, MSE1, type ="l", xlab = "x", ylab = "MSE", col = "black", ylim=c(0, 1))

lines(x0, MSE2, xlab = "x", ylab = "MSE", col = "blue")

lines(x0, MSE3, xlab = "x", ylab = "MSE", col = "red2")

lines(x0, MSE5, xlab = "x", ylab = "MSE",col = "green")

legend(-1, 1, legend=c("Second-order polynomial", "Cubic spline", "Natural cubic spline", "Smoothing spline"),
       col=c("black", "blue", "red2", "green"), lwd=1, box.col = "white", cex=0.5)


#cubic spline

panel.smooth <- function(x,y,...)
{
  knots=quantile(x,c(1/3,2/3))
  fit = lm(y~bs(x,knots=knots,degree=3))
  x0 = seq(min(x),max(x),len=300)
  y0 = predict(fit,newdata=data.frame(x=x0),se.fit=T)
  #lpolygon(c(x0,rev(x0)),c(y0$fit-2*y0$se.fit,rev(y0$fit+2*y0$se.fit)),col="gray",border=F,...)
  llines(x0,predict(fit,data.frame(x=x0)),col="black")
  lsegments(x0=knots[1],y0=min(y),x1=knots[1],y1=max(y),col="black",lty=2)
  lsegments(x0=knots[2],y0=min(y),x1=knots[2],y1=max(y),col="black",lty=2)
  
}
xyplot(y~x, panel = myPanel)


var <- (y0$se.fit)^2
plot(x0, var, xlab = "x", ylab = "variance")
y_true <- -2*(x0^2)
bias <- y0$fit - y_true
plot(x0, bias, xlab = "x", ylab = "bias")
MSE <- bias^2 + var
plot(x0, MSE, xlab = "x", ylab = "MSE")

################################################
##fitting natural spline
panel.smooth <- function(x,y,...)
{
  knots=quantile(x,c(1/3,2/3))
  fit = lm(y~ns(x,knots=knots))
  x0 = seq(min(x),max(x),len=300)
  llines(x0,predict(fit,data.frame(x=x0)),col="black")
  #lsegments(x0=min(x),y0=min(y),x1=min(x),y1=max(y),col="black",lty=2)
  #lsegments(x0=knots[1],y0=min(y),x1=knots[1],y1=max(y),col="black",lty=2)
  #lsegments(x0=knots[2],y0=min(y),x1=knots[2],y1=max(y),col="black",lty=2)
  #lsegments(x0=max(x),y0=min(y),x1=max(x),y1=max(y),col="black",lty=2)
}
xyplot(y~x, panel= myPanel)


# natural cubic spline, df = 5

panel.smooth <- function(x,y,...)
{
  x0 = seq(min(x),max(x),len=101)
  y0 = predict(fit,newdata=data.frame(x=x0),se.fit=T)
  #lpolygon(c(x0,rev(x0)),c(y0$fit-2*y0$se.fit,rev(y0$fit+2*y0$se.fit)),col="gray",border=F,...)
  llines(x0,predict(fit,data.frame(x=x0)),col="black")
  
}
xyplot(y~x, panel = myPanel)

var <- (y0$se.fit)^2
plot(x0, var, xlab = "x", ylab = "variance")
y_true <- -2*(x0^2)
bias <- y0$fit - y_true
plot(x0, bias, xlab = "x", ylab = "bias")
MSE <- bias^2 + var
plot(x0, MSE, xlab = "x", ylab = "MSE")

# choose knots for female
x = age[gender=="female"]
ns(x,df=5)

## mean and variace of the estimate
panel.smooth <- function(x,y,...)
{
  fit = lm(y~ns(x,df=5))
  x0 = seq(min(x),max(x),len=300)
  y0 = predict(fit,newdata=data.frame(x=x0),se.fit=T)
  #lpolygon(c(x0,rev(x0)),c(y0$fit-2*y0$se.fit,rev(y0$fit+2*y0$se.fit)),col="gray",border=F,...)
  llines(x0,y0$fit,col="black")
}
xyplot(y~x,panel=myPanel)

################################################
##### smoothing spline
################################################


panel.smooth <- function(x,y,...)
{
fit6<- smooth.spline(x, y, cv=T)
fit6
x0 = seq(min(x),max(x),len=300)
y0 = predict(fit,newdata=data.frame(x=x0),se.fit=T)
#lpolygon(c(x0,rev(x0)),c(y0$fit-2*y0$se.fit,rev(y0$fit+2*y0$se.fit)),col="gray",border=F,...)
llines(fit, lwd = 3, col = "black")
var <- (y0$se.fit)^2
plot(x0, var, type = "l", xlab = "x", ylab = "variance")
str(y0)
y_true <- -2*(x0^2)
bias <- y0$fit - y_true
plot(x0, bias, xlab = "x", ylab = "bias")
MSE <- bias^2 + var
}
xyplot(y ~ x, panel = myPanel)


require(mgcv)
myPanel <- function(...)
{
  panel.xyplot(...,alpha=.3)
  panel.smooth(...)
}
panel.smooth <- function(x,y,...)
{
  fit5 <- gam(y~s(x))
  xx <- seq(min(x),max(x),len=300)
  yy <- predict(fit5,newdata=data.frame(x=xx),se.fit=T)
  var5 <- yy$se.fit^2
  plot(x0, var5, type = "l", xlab = "x", ylab = "variance")
  lpolygon(c(xx,rev(xx)),c(yy$fit-1.96*yy$se.fit,rev(yy$fit+1.96*yy$se.fit)),col=rgb(.6,.6,.6,alpha=.4),border=F,...)
  llines(xx,yy$fit,col="black",lwd=2)
  llines(xx,yy$se.fit,col="black",lwd=2)
}
xyplot(y~x,panel=myPanel)
fit5

require(mgcv)
myPanel <- function(...)
{
  panel.xyplot(...,alpha=.3)
  panel.smooth(...)
}
panel.smooth <- function(x,y,...)
{
  fit <- gam(y~s(x))
  xx <- seq(min(x),max(x),len=101)
  yy <- predict(fit,newdata=data.frame(x=xx),se.fit=T)
  lpolygon(c(xx,rev(xx)),c(yy$fit-1.96*yy$se.fit,rev(yy$fit+1.96*yy$se.fit)),col=rgb(.6,.6,.6,alpha=.4),border=F,...)
  llines(xx,yy$fit,col="black",lwd=2)
}
trellis.par.set(plot.symbol=list(pch=19))
xyplot(spnbmd~age|gender,bone,panel=myPanel)


#how to add standard error into this plot?

plot(y ~ x)

f.spline <- smooth.spline(x, y,df=12)

lines(f.spline)

panel.smooth <- function(x,y,...)
{
  f.spline <- smooth.spline(x, y, df=12)
  x0 = seq(min(x),max(x),len=300)
  llines(x0,predict(f.spline, data.frame(x=x0)), col="black", lwd=3)
}
xyplot(y~x, panel=myPanel)


# confidence interval for smoothing spline
# second way
require(mgcv)
myPanel <- function(...)
{
  panel.xyplot(...,alpha=.3)
  panel.smooth(...)
}
panel.smooth <- function(x,y,...)
{
  fit <- gam(y~s(x))
  xx <- seq(min(x),max(x),len=101)
  yy <- predict(fit,newdata=data.frame(x=xx),se.fit=T)
  lpolygon(c(xx,rev(xx)),c(yy$fit-1.96*yy$se.fit,rev(yy$fit+1.96*yy$se.fit)),col=rgb(.6,.6,.6,alpha=.4),border=F,...)
  llines(xx,yy$fit,col="black",lwd=2)
}
trellis.par.set(plot.symbol=list(pch=19))
xyplot(spnbmd~age|gender,bone,panel=myPanel)

