install.packages("lasso2")
library(lasso2)
data(Prostate)
attach(Prostate)

# produce pairwise plots
pairs(cbind(lpsa,lcavol,lweight,age,lbph,svi,lcp,gleason))

# fit the linear model
fit = lm(lpsa~lcavol+lweight+age+lbph+svi+lcp+gleason,data=Prostate)


# summary of the fitted model
summary(fit)


# anova table
anova(fit)


# prediction
xnew= data.frame(lcavol = -1, lweight = 2.5, age =50,
                 lbph= 1.4, svi = 0,lcp= -1.3,gleason = 7)

# predicted CI
predict(fit,xnew,se.fit = FALSE, interval = "prediction")

# confidence CI 
predict(fit,xnew,se.fit = TRUE, interval = "confidence")


# check the assumptions
par(mfrow=c(2,2))
qqnorm(fit$res)
qqline(fit$res)
qqnorm(abs(fit$res))
qqline(abs(fit$res))
plot(fit$fitted,fit$res,xlab="Fitted values",ylab="Residuals")
plot(1:dim(Prostate)[1],fit$res,xlab="Runs",ylab="Residuals")

# compare two models
fit2 = lm(lpsa~lcavol+lweight+lbph+svi,data=Prostate)
anova(fit2,fit)



# remove the names of data Prostate
detach(Prostate)
