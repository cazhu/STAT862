## read data
bone = read.table("bone.txt",header=T)
attach(bone)


### plot of the data
install.packages("lattice")
library(lattice)
xyplot(spnbmd~age|gender, bone)

## install package
install.packages("splines")
library(splines)

# adding lines in a panel
myPanel <- function(...)
{
  panel.smooth(...)
  panel.xyplot(...,alpha=.5)
}


############################################
## fitting piecewise constant model
panel.smooth <- function(x,y,...)
{
  knots=quantile(x,c(1/3,2/3))
  c1 <- mean(y[x < knots[1]])
  c2 <- mean(y[x >= knots[1] & x < knots[2]])
  c3 <- mean(y[x >= knots[2]])
  x1 <- seq(min(x),knots[1],len=33)
  x2 <- seq(knots[1],knots[2],len=33)
  x3 <- seq(knots[2],max(x),len=33)
  llines(x1,c1,col="black")
  lsegments(x0=knots[1],y0=min(y),x1=knots[1],y1=max(y),col="black",lty=2)
  llines(x2,c2,col="black")
  lsegments(x0=knots[2],y0=min(y),x1=knots[2],y1=max(y),col="black",lty=2)
  llines(x3,c3,col="black")
}
trellis.par.set(plot.symbol=list(pch=19),plot.line=list(lwd=2))
pdf("fig3_bone_piecewiseconstant.pdf")
xyplot(spnbmd~age|gender,bone,panel=myPanel)
dev.off()


# choose knots for female
x = age[gender=="female"]
knots=quantile(x,c(1/3,2/3))

############################################
## fitting piecewise linear model
panel.smooth <- function(x,y,...)
{
  knots=quantile(x,c(1/3,2/3))
  ind1 <- x < knots[1]
  ind2 <- x >= knots[1] & x < knots[2]
  ind3 <- x >= knots[2]
  x1 <- x[ind1]
  x2 <- x[ind2]
  x3 <- x[ind3]
  fit1 <- lm(y[ind1]~x1)
  fit2 <- lm(y[ind2]~x2)
  fit3 <- lm(y[ind3]~x3)
  x1 <- seq(min(x),knots[1],len=33)
  x2 <- seq(knots[1],knots[2],len=33)
  x3 <- seq(knots[2],max(x),len=33)
  llines(x1,predict(fit1,data.frame(x1=x1)),col="black")
  lsegments(x0=knots[1],y0=min(y),x1=knots[1],y1=max(y),col="black",lty=2)
  llines(x2,predict(fit2,data.frame(x2=x2)),col="black")
  lsegments(x0=knots[2],y0=min(y),x1=knots[2],y1=max(y),col="black",lty=2)
  llines(x3,predict(fit3,data.frame(x3=x3)),col="black")
}

pdf("fig3_bone_piecewiselinear.pdf")
xyplot(spnbmd~age|gender,bone,panel=myPanel)
dev.off()

############################################
## fitting continous piecewise linear model
panel.smooth <- function(x,y,...)
{
  knots=quantile(x,c(1/3,2/3))
  fit <- lm(y~bs(x,knots=knots,degree=1))
  x0 <- seq(min(x),max(x),len=101)
  llines(x0,predict(fit,data.frame(x=x0)),col="black")
  lsegments(x0=knots[1],y0=min(y),x1=knots[1],y1=max(y),col="black",lty=2)
  lsegments(x0=knots[2],y0=min(y),x1=knots[2],y1=max(y),col="black",lty=2)
}
pdf("fig3_bone_continuouspiecewiselinear.pdf")
xyplot(spnbmd~age|gender,bone,panel=myPanel)
dev.off()

############################################
## fitting qudractic spline
panel.smooth <- function(x,y,...)
{
  knots=quantile(x,c(1/3,2/3))
  fit <- lm(y~bs(x,knots=knots,degree=2))
  x0 <- seq(min(x),max(x),len=101)
  llines(x0,predict(fit,data.frame(x=x0)),col="black")
  lsegments(x0=knots[1],y0=min(y),x1=knots[1],y1=max(y),col="black",lty=2)
  lsegments(x0=knots[2],y0=min(y),x1=knots[2],y1=max(y),col="black",lty=2)
}
xyplot(spnbmd~age|gender,bone,panel=myPanel)


################################################
## fitting cubic spline
panel.smooth <- function(x,y,...)
{
  knots=quantile(x,c(1/3,2/3))
  fit = lm(y~bs(x,knots=knots,degree=3))
  x0 = seq(min(x),max(x),len=101)
  llines(x0,predict(fit,data.frame(x=x0)),col="black")
  lsegments(x0=knots[1],y0=min(y),x1=knots[1],y1=max(y),col="black",lty=2)
  lsegments(x0=knots[2],y0=min(y),x1=knots[2],y1=max(y),col="black",lty=2)
}
xyplot(spnbmd~age|gender,bone,panel=myPanel)

################################################
##fitting natural spline
panel.smooth <- function(x,y,...)
{
  knots=quantile(x,c(1/3,2/3))
  fit = lm(y~ns(x,knots=knots))
  x0 = seq(min(x),max(x),len=101)
  llines(x0,predict(fit,data.frame(x=x0)),col="black")
  lsegments(x0=min(x),y0=min(y),x1=min(x),y1=max(y),col="black",lty=2)
  lsegments(x0=knots[1],y0=min(y),x1=knots[1],y1=max(y),col="black",lty=2)
  lsegments(x0=knots[2],y0=min(y),x1=knots[2],y1=max(y),col="black",lty=2)
  lsegments(x0=max(x),y0=min(y),x1=max(x),y1=max(y),col="black",lty=2)
}
xyplot(spnbmd~age|gender,bone,panel=myPanel)

################################################
## fitting natural spline with df = 5
panel.smooth <- function(x,y,...)
{
  fit = lm(y~ns(x,df=5))
  x0 = seq(min(x),max(x),len=101)
  llines(x0,predict(fit,data.frame(x=x0)),col="black",lwd=3)
}
xyplot(spnbmd~age|gender,bone,panel=myPanel)

# choose knots for female
x = age[gender=="female"]
ns(x,df=5)

## mean and variace of the estimate
panel.smooth <- function(x,y,...)
{
  fit = lm(y~ns(x,df=5))
  x0 = seq(min(x),max(x),len=101)
  y0 = predict(fit,newdata=data.frame(x=x0),se.fit=T)
  lpolygon(c(x0,rev(x0)),c(y0$fit-2*y0$se.fit,rev(y0$fit+2*y0$se.fit)),col="gray",border=F,...)
  llines(x0,y0$fit,col="black")
}
xyplot(spnbmd~age|gender,bone,panel=myPanel)



################################################
##### smoothing spline
################################################

plot(spnbmd ~ age, data=bone, col =
       ifelse(gender=="male", "blue", "red2"),
     xlab="Age", ylab="Relative Change in Spinal BMD")
bone.spline.male <- with(subset(bone,gender=="male"),
                         smooth.spline(age, spnbmd,df =12))
bone.spline.female <- with(subset(bone, gender=="female"),
                           smooth.spline(age, spnbmd, df =12))
lines(bone.spline.male, col="blue")
lines(bone.spline.female, col="red2")
legend(20,0.20,legend=c("male", "Female"),
       col=c("blue", "red2"),lwd=2)



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