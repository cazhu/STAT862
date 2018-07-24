install.packages("xtable")
library(MASS)
library(xtable)
## read data
admit = read.table("admit.txt",header=T)
attach(admit)

####logistic regression
## model fitting
model_logit = glm(admit ~ . , family = binomial(link='logit'), data=admit)
require(nnet)
model_logit2 <- multinom(admit~.,data=admit)


summary(model_logit)
summary(model_logit2)
xtable(model_logit)

sum(c(1,as.numeric(admit[1,1:3]))*model_logit$coef)

# computing p hat
exp(-0.6494346)/(1+exp(-0.6494346)) 
model_logit$fit[1]