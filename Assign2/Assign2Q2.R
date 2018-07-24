library(MASS)
data(mcycle)

names(mcycle)

# (a)

mcycle.spline <- smooth.spline(mcycle$times, mcycle$accel, cv = T)
plot(mcycle)
lines(mcycleT.spline)

# (b)




mcycle5.spline <- smooth.spline(mcycle$times, mcycle$accel, df = 5)
mcycle5.spline$cv.crit
mcycle10.spline <- smooth.spline(mcycle$times, mcycle$accel, df = 10)
mcycle10.spline$cv.crit
mcycle15.spline <- smooth.spline(mcycle$times, mcycle$accel, df = 15)
mcycle15.spline$cv.crit
plot(mcycle)
lines(mcycle5.spline, col = "blue")
lines(mcycle10.spline, col = "red2")
lines(mcycle15.spline, col ="black")
legend(40, -40, legend=c("df = 5", "df = 10", "df = 15"),
       col=c("blue", "red2", "black"),lwd=1, box.col = "white")

PRESS <- 0
for(i in 1:31)
  {
mcyclecv.spline <- smooth.spline(mcycle$times, mcycle$accel, df = 5+0.5*(i-1), cv = T)
PRESS[i] <- mcyclecv.spline$cv.crit
}
plot (seq(5, 20, 0.5), PRESS, xlab = "df", ylab = "PRESS")
lines (seq(5, 20, 0.5), PRESS, xlab = "df", ylab = "PRESS")
mcycle5cv.spline <- smooth.spline(mcycle$times, mcycle$accel, df = 5, cv = T)
mcycle10cv.spline <- smooth.spline(mcycle$times, mcycle$accel, df = 10, cv = T)
mcycle15cv.spline <- smooth.spline(mcycle$times, mcycle$accel, df = 15, cv = T)
mcycle20cv.spline <- smooth.spline(mcycle$times, mcycle$accel, df = 20, cv = T)
plot(mcycle)
lines(mcycle5cv.spline, col = "blue")
lines(mcycle10cv.spline, col = "red2")
lines(mcycle15cv.spline, col ="black")
legend(40, -40, legend=c("df = 5", "df = 10", "df = 15"),
       col=c("blue", "red2", "black"),lwd=1, box.col = "white")
# (c)
GCV <- numeric(31)
a = 1
for (i in 5:20) {
  mcycle.spline <- smooth.spline(mcycle$times, mcycle$accel, df = i)
 GCV[a] <- mcycle.spline$cv.crit
 a = a + 1
  i = i+0.5
  
}
plot(seq(5,20,0.5),GCV,xlab= "df",ylab="gcv")

mcycle.spline <- smooth.spline(mcycle$times, mcycle$accel, df = seq(5,20,0.5))

str(mcycle.spline)
plot(mcycle)
lines(mcycle.spline)
mcycle.spline$cv.crit
mcycle.spline
plot(seq(5,20,0.5),mcycle.spline$cv.crit,xlab= "df",ylab="gcv")
mcycle.spline$cv.crit