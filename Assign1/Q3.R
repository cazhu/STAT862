####################################################
## Ridge Regression
####################################################

library(MASS)  

Housing <- read.table("Housing.dat", header = TRUE)

Housing = scale(Housing)
Housing = as.data.frame(Housing)
e=11
fit.ridge<-lm.ridge(MEDV~., data=Housing, lambda=seq(0,11,0.01))  
par(mfrow=c(1,2))
plot(fit.ridge)  
plot(seq(0,11,0.01),fit.ridge$GCV,xlab= expression(lambda),ylab="GCV")

select(fit.ridge)

r1 <- round(fit.ridge$coef[, which(fit.ridge$lambda == 1)], 3) 
r3 <- round(fit.ridge$coef[, which(fit.ridge$lambda == 3)], 3) 
r5 <- round(fit.ridge$coef[, which(fit.ridge$lambda == 5)], 3) 
r7 <- round(fit.ridge$coef[, which(fit.ridge$lambda == 7)], 3) 
r9 <- round(fit.ridge$coef[, which(fit.ridge$lambda == 9)], 3) 
r11 <- round(fit.ridge$coef[, which(fit.ridge$lambda == 11)], 3) 

rtotal <- rbind(r1, r3, r5, r7, r9, r11)
xtable(rtotal)

round(fit.ridge$coef[, which(fit.ridge$lambda == 4.26)], 3)

fit.ridge<-lm.ridge(MEDV~., data = Housing, lambda= 2)   
fit.ridge$coef     
