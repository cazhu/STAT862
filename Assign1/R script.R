install.packages("leaps")
library(leaps)

install.packages("xtable")
library(xtable)

install.packages("lasso2")
library(lasso2)

Grocery <- read.table("Grocery.txt", header= F)

colnames(Grocery) <- c("Y","X1","X2","X3")

# produce pairwise plots
pairs(Grocery, labels = c("Y", "X1", "X2", "X3"), main = "")

# fit the linear model
fit = lm(Grocery$Y~.,data=Grocery)


# summary of the fitted model


summary(fit)$r.squared
summary(fit)

xtable(summary(fit))

# anova table
anova(fit)
xtable(anova(fit))

### best subset selction 
leaps=regsubsets(Y~X1+X2+X3, data=Grocery,nbest=10)

#### 
pdf("bestsubset_leaps.pdf")
par(mfrow=c(2,2))
plot(leaps,scale="r2")
plot(leaps,scale="adjr2")
plot(leaps,scale="Cp")
plot(leaps,scale="bic")
dev.off()

subset<-leaps(x=Grocery[,2:4],y=Grocery[,1])
plot(x=subset$size,y=subset$Cp,xlab='size',ylab='Cp')

### 
summary(leaps)
reg.summary=summary(leaps)

### the residual sum of squares of the top 10 models
reg.summary$rss

### the R^2 of the top 10 models
reg.summary$rsq

### the adjusted R^2 of the top 10 models
reg.summary$adjr2

### the Cp of the top 10 models
reg.summary$cp



#  plot of RSS, adjusted R^2, Cp and BIC together
par(mfrow=c(2,2))
plot(reg.summary$rss, xlab="Number of Predictors", ylab="Residual Sum of Squares", type="l", xlim=c(0,4), ylim=c(min(reg.summary$rss), max(reg.summary$rss)))
points(which.min(reg.summary$rss), reg.summary$rss[which.min(reg.summary$rss)], cex=2, pch=20, col="red")

plot(reg.summary$cp, xlab="Number of Predictors", ylab="Cp", type="l", xlim=c(0,3), ylim=c(min(reg.summary$cp),max(reg.summary$cp)))
points(which.min(reg.summary$cp), reg.summary$cp[which.min(reg.summary$cp)], cex=2, pch=20, col="red")

plot(reg.summary$adjr2, xlab="Number of Predictors", ylab="Adjusted R Square", type="l", xlim=c(0,3), ylim=c(0,1))
points(which.max(reg.summary$adjr2),reg.summary$adjr2[which.max(reg.summary$adjr2)], cex=2, pch=20, col="red")

plot(reg.summary$bic, xlab="Number of Predictors", ylab="BIC", type="l", xlim=c(0,3))
points(which.min(reg.summary$bic),reg.summary$bic[which.min(reg.summary$bic)], cex=2, pch=20, col="red")


### full model
fit1<-lm(Y~ .,data=Grocery)

### null model
fit0<-lm(Y~1,data=Grocery)

### forward selection
fit.forward<-step(fit0,scope=list(lower=Y~1, upper=Y~.),direction='forward')
summary(fit.forward)

xtable(fit.forward)
xtable(anova(fit.forward))

### backward selection
fit.backward<-step(fit1,scope=list(lower=Y~1, upper=Y~.),direction='backward')
summary(fit.backward)

xtable(fit.backward)
xtable(anova(fit.backward))
       
# F-test between smaller model and full model
fit3<-lm(Y~X1+X3, data=Grocery)
summary(fit3)
anova(fit3, fit)

xtable(anova(fit3, fit))

 
### stepwise regression by AIC criterion
fit.both<-step(fit0,scope=list(lower=Y~1, upper=Y~.),direction='both')
summary(fit.both)  

# prediction
xnew = data.frame(X1 = 32000, X2 =7.5, X3= 1)

# predicted CI
# predict(fit3,xnew,se.fit = FALSE, interval = "prediction")

# confidence CI 
predict(fit3,xnew,se.fit = TRUE, interval = "confidence")

# check the assumptions
par(mfrow=c(2,2))
qqnorm(fit$res)
qqline(fit$res)
qqnorm(abs(fit$res))
qqline(abs(fit$res))
plot(fit$fitted,fit$res,xlab="Fitted values",ylab="Residuals")
plot(1:dim(Grocery)[1],fit$res,xlab="Runs",ylab="Residuals")

####################################################
## Lasso
####################################################

Z1 <- Grocery$X1*Grocery$X2
Z2 <- Grocery$X1*Grocery$X3
Z3 <- Grocery$X2*Grocery$X3

n <- nrow(Grocery)

Z4 <- rnorm(n, mean = 30, sd = 30)
Z5 <- rnorm(n, mean = 7, sd = 1)

#Add new columns into Grocery

Grocery2 <- cbind(Grocery, Z1, Z2, Z3, Z4, Z5)

#Lasso

X = as.matrix(Grocery2[,2:9])
y = Grocery2[,1]

##lars
install.packages("lars")
library(lars)
fit.lars = lars(X,y, type="lasso",trace=TRUE)
plot(fit.lars)
cv.fit.lars = cv.lars(X,y,mode="step")
bestindex = cv.fit.lars$index[which.min(cv.fit.lars$cv)]
which.min(cv.fit.lars$cv)
fit.lars$beta
bestindex
fit.lars$lambda
fit.lars$beta[6,]

###############

install.packages("glmnetcr")
library(glmnetcr)

glmnet.fit = glmnet.cr(X, y)
AIC = select.glmnet.cr(glmnet.fit, which = "AIC")
AIC
nonzero.glmnet.cr(glmnet.fit, s = AIC)

BIC = select.glmnet.cr(glmnet.fit, which = "BIC")
BIC
nonzero.glmnet.cr(glmnet.fit, s = BIC)






