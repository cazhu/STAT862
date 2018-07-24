library(MASS)

#############  Data SAheart#######################################################
# read data
SAheart<-read.table('SAheart.txt',header=T)
attach(SAheart)



####logistic regression
## model fitting
model_logit = glm(chd ~ . , family = binomial(link='logit'), data=SAheart)
require(nnet)
model_logit2 <- multinom(chd~.,data=SAheart)

## model outputs
summary(model_logit)
str(model_logit)
model_logit$coef
model_logit$fit
as.numeric(model_logit$fit>0.5)


sum(c(1,as.numeric(SAheart[1,1:8]))*model_logit$coef)

# computing p hat
exp(0.4972209)/(1+exp(0.4972209)) 
model_logit$fit[1]




summary(model_logit2)
str(model_logit2)
model_logit2
model_logit2$fit
as.numeric(model_logit2$fit>0.5)


#### test data
Test_saheart = SAheart[1:50,1:8]

## prediction
pred_logit= predict(model_logit,Test_saheart, type="response")
as.numeric(pred_logit>=0.5)
pred_logit2= predict(model_logit2,Test_saheart)


###### linear discriminant analaysis
model_lda = lda(chd~sbp+tobacco+ldl+adiposity+
                  typea+obesity+alcohol+age,data=SAheart)
model_lda
mean(SAheart$chd)
mean(SAheart$sbp[SAheart$chd==0])
#Prior probabilities of groups 0: pie zero hat, 1: pia one hat
#Group means 0 miu zero hat, 1 miu one hat
# prediction
pred_lda = predict(model_lda,Test_saheart )

predict_lda$class

# quadratic discriminant analysis
model_qda = qda(chd~sbp+tobacco+ldl+adiposity+
                  typea+obesity+alcohol+age,data=SAheart)
model_qda

pred_qda = predict(model_qda,Test_saheart)


# comparison
cbind(pred_lda$class, pred_qda$class, pred_logit2,chd[1:50]+1)
table(pred_lda$class, chd[1:50])
table(pred_qda$class, chd[1:50])
table(pred_logit2, chd[1:50])



detach(SAheart)

##############  Data iris#######################################################

# read data
iris<-read.table('iris.txt',header=T)
attach(iris)


n= dim(iris)[1]
# randomly row-permutate data
iris2 = iris[sample(1:n),]

# linear discriminant analaysis

training = iris2[1:100,]
model_lda = lda( Species ~ SepalLength+ SepalWidth + 
                   PetalLength+ PetalWidth,data= training)

# test data
Test = iris2[101:150,1:4]

# prediction
pred_lda = predict(model_lda,Test)


# quadratic discriminant analysis
model_qda = qda(Species ~ .,data=training)
pred_qda = predict(model_qda,Test)

# comparison
data.frame(pred_lda$class, pred_qda$class,iris2[101:150,5])

table(pred_lda$class, iris2[101:150,5])
table(pred_qda$class, iris2[101:150,5])

detach(iris)




