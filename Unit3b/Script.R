install.packages("leaps")
install.packages("lasso2")
library(leaps)
library(lasso2)
data(Prostate)
attach(Prostate)


####################################################
## Subset Selection
####################################################

### best subset selction 
leaps=regsubsets(lpsa~lcavol+lweight+age
                 +lbph+svi+lcp+gleason+pgg45, data=Prostate,nbest=10)

leaps=regsubsets(lpsa~., data=Prostate,nvmax=10)

#### 
pdf("bestsubset_leaps.pdf")
par(mfrow=c(2,2))
plot(leaps,scale="r2")
plot(leaps,scale="adjr2")
plot(leaps,scale="Cp")
plot(leaps,scale="bic")
dev.off()

subset<-leaps(x=Prostate[,1:8],y=Prostate[,9])
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
plot(reg.summary$rss, xlab="Number of Predictors", ylab="Residual Sum of Squares", type="l", xlim=c(0,11), ylim=c(min(reg.summary$rss), max(reg.summary$rss)))
points(which.min(reg.summary$rss), reg.summary$rss[which.min(reg.summary$rss)], cex=2, pch=20, col="red")

plot(reg.summary$cp, xlab="Number of Predictors", ylab="Cp", type="l", xlim=c(0,11), ylim=c(min(reg.summary$cp),max(reg.summary$cp)))
points(which.min(reg.summary$cp), reg.summary$cp[which.min(reg.summary$cp)], cex=2, pch=20, col="red")

plot(reg.summary$adjr2, xlab="Number of Predictors", ylab="Adjusted R Square", type="l", xlim=c(0,11), ylim=c(0,1))
points(which.max(reg.summary$adjr2),reg.summary$adjr2[which.max(reg.summary$adjr2)], cex=2, pch=20, col="red")

plot(reg.summary$bic, xlab="Number of Predictors", ylab="BIC", type="l", xlim=c(0,11))
points(which.min(reg.summary$bic),reg.summary$bic[which.min(reg.summary$bic)], cex=2, pch=20, col="red")


### full model
fit1<-lm(lpsa~lcavol+lweight+age+lbph+svi+lcp+gleason+pgg45,
         data=Prostate)
fit1<-lm(lpsa~ .,data=Prostate)

### null model
fit0<-lm(lpsa~1,data=Prostate)

### forward selection
fit.forward<-step(fit0,scope=list(lower=lpsa~1, upper=fit1),direction='forward')
summary(fit.forward)

### backward selection
fit.backward<-step(fit1,scope=list(lower=lpsa~1, upper=fit1),direction='backward')
summary(fit.backward)  

### stepwise regression by AIC criterion
fit.both<-step(fit0,scope=list(lower=lpsa~1, upper=fit1),direction='both')
summary(fit.both)  


####################################################
## Ridge Regression
####################################################

library(MASS)  

# in the presence of multi-collinearity 

x1  = rnorm(30)
x2  = rnorm(30,mean=x1,sd=.01)
y   = rnorm(30,mean=5+x1+x2)
lm(y~x1+x2)$coef
lm.ridge(y~x1+x2,lambda=1)

# Prostate example
prostate = scale(Prostate)
prostate = as.data.frame(prostate)
fit.ridge<-lm.ridge(lpsa~lcavol+lweight+age
                    +lbph+svi+lcp+gleason+pgg45,
                    data=prostate, lambda=seq(0,20,0.1))  
plot(fit.ridge)  
plot(seq(0,20,0.1),fit.ridge$GCV,xlab= expression(lambda),ylab="GCV")

select(fit.ridge)
round(fit.ridge$coef[, which(fit.ridge$lambda == 6.5)], 2)    
fit.ridge<-lm.ridge(lpsa~lcavol+lweight+age  
                    +lbph+svi+lcp+gleason+pgg45,   
                    data=prostate, lambda=6.5)   
fit.ridge$coef     



####################################################
## Lasso
####################################################
install.packages("glmnet")
library(glmnet)
X = as.matrix(Prostate[,1:8])
y = Prostate$lpsa

fit = glmnet(X,y)
plot(fit)
cvfit = cv.glmnet(X,y)
plot(cvfit)
coef(fit,s=cvfit$lambda.min)
min(cvfit$cvm)

##
install.packages("lars")
library(lars)
fit.lars = lars(X,y, type="lasso",trace=TRUE)
plot(fit.lars)
cv.fit.lars = cv.lars(X,y)
bestindex = cv.fit.lars$index[which.min(cv.fit.lars$cv)]
which.min(cv.fit.lars$cv)
fit.lars$beta
bestindex
fit.lars$lambda
fit.lars$beta[6,]



install.packages("glmnetcr")
library(glmnetcr)

glmnet.fit = glmnet.cr(X, y)
AIC = select.glmnet.cr(glmnet.fit, which = "AIC")
AIC
nonzero.glmnet.cr(glmnet.fit, s = AIC)

BIC = select.glmnet.cr(glmnet.fit, which = "BIC")
BIC
nonzero.glmnet.cr(glmnet.fit, s = BIC)



detach(Prostate)


