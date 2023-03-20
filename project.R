setwd("/Users/enxilin/Desktop/csv")
bike=read.csv('Seoul Bike Sharing Demand Data Set.csv')
bike$Dayofweek=weekdays(as.Date(bike$Date,"%d/%m/%Y"))
plot(bike$Hour,bike$Count,xlab="Hour",ylab="Count",main="Bike renting demand hourly")
bike$weekstatus=ifelse(bike$Dayofweek=='Saturday'|bike$Dayofweek=='Sunday',"weekend","weekday")
bike=bike[bike$Functioning.Day=='Yes',]
bike=bike[bike$Hour=='18',]
bike=bike[ , -which(colnames(bike) %in% c('Date','Functioning.Day','Hour'))]
bike$Rain=as.factor(ifelse(bike$Rain==0, 0, 1))
bike$Snow=as.factor(ifelse(bike$Snow==0, 0, 1))
set.seed(13)
n=nrow(bike)
train=sample(n,trunc(0.7*n))
training=bike[train,]
testing=bike[-train,]
training$Seasons=as.factor(training$Seasons)
training$Holiday=as.factor(training$Holiday)
training$weekstatus=as.factor(training$weekstatus)
training$Rain=as.factor(ifelse(training$Rain==0, 0, 1))
training$Snow=as.factor(ifelse(training$Snow==0, 0, 1))

#scaterplot matrix
if(require(car)){scatterplotMatrix(training[,c(1:7)],smooth=FALSE, diagonal="histogram")}
training$Count=sqrt(training$Count)
hist(training$Count,main='square root transformation of Count',xlab="sqrt(Count)")

#interactions
bike$Count=sqrt(bike$Count)
interaction.plot(bike$Snow,bike$Seasons,bike$Count,xlab = "Snow",ylab="Count",main="Snow vs Season")
interaction.plot(bike$Snow,bike$Holiday,bike$Count,xlab = "Snow",ylab="Count",main="Snow vs Holiday")
interaction.plot(bike$Rain,bike$Seasons,bike$Count,xlab = "Rain",ylab="Count",main="Rain vs Season")
interaction.plot(bike$Snow,bike$weekstatus,bike$Count,xlab = "Snow",ylab="Count",main="Snow vs Weekstatus")
interaction.plot(bike$Rain,bike$Holiday,bike$Count,xlab = "Rain",ylab="Count",main="Rain vs Holiday")
interaction.plot(bike$Holiday,bike$Seasons,bike$Count,xlab = "Holiday",ylab="Count",main="Holiday vs Seasons")
interaction.plot(bike$Rain,bike$weekstatus,bike$Count,xlab = "Rain",ylab="Count",main="Rain vs Weekstatus")#no
interaction.plot(bike$weekstatus,bike$Holiday,bike$Count,xlab = "weekstatus",ylab="Count",main="weekstatus vs Holiday")#no
interaction.plot(bike$weekstatus,bike$Seasons,bike$Count,xlab = "weekstatus",ylab="Count",main="weekstatus vs Seasons")#no

# check multicolinearity
g=lm(Count~Temp+Hum+Wind+Visb+Dew+Solar+Rain+Snow+Seasons+Holiday+weekstatus+Rain*Seasons+Snow*weekstatus+Rain*Holiday+Holiday*Seasons,data=training)
vif(g)
library(olsrr)
ols_eigen_cindex(g)$Eigenvalue
max(ols_eigen_cindex(g)$Eigenvalue)/min(ols_eigen_cindex(g)$Eigenvalue)

##varable selection
#backward elimination:Count ~ Hum + Dew + Solar + Rain + Seasons + Holiday + weekstatus + Rain:Seasons
step(g,direction = 'backward',trace=TRUE)
back_model=lm(Count ~ Temp+ Hum + Dew + Solar + Rain + Seasons + Holiday + weekstatus + Rain:Seasons, data = training)
summary(back_model)
back_model=update(back_model,~.-Temp)
summary(back_model)
back_model=lm(Count ~ Hum + Dew + Solar + Rain + Seasons + Holiday + weekstatus + Rain:Seasons, data = training)
## conscientious approach
summary(g)
g1=update(g,~.-Seasons:Holiday)
summary(g1)
g2=update(g1,~.-Rain:Holiday)
summary(g2)
g3=update(g2,~.-Snow:weekstatus)
summary(g3)
g4=update(g3,~.-Snow)
summary(g4)
g5=update(g4,~.-Visb)
summary(g5)
g6=update(g5,~.-Wind)
summary(g6)

#ridge regression
x=cbind(training$Temp,training$Hum,training$Wind,training$Visb,training$Dew,training$Solar,training$Rain,training$Snow,
        training$Seasons,training$Holiday,training$weekstatus,training$Rain:training$Seasons,training$Snow:training$weekstatus,
        training$Rain:training$Holiday,training$Seasons:training$Holiday)
dimnames(x) = list(NULL, c("Tmep", "Hum", "Wind", "Visb", "Dew", "Solar", "Rain", "Snow", "Seasons","Holiday","weekstatus",
                           "Rain*Seasons","Snow*weekstatus","Rain:Holiday","Seasons*Holiday"))
n=length(training$Count)
m=apply(x, 2, mean)
z=x-outer(rep(1,n),m)
s=sqrt((n-1)*apply(x, 2, var))
z = z/outer(rep(1,n),s)
r= crossprod(z)
y = (training$Count - mean(training$Count))/sqrt((n-1)*var(training$Count))
v = eigen(r)$values;max(v)/min(v)
p = ncol(z)
beta = function(k) {solve(r+k*diag(p), t(z)%*%y)}
K = seq(from = 0.001, to = .3, by = .001)
THETA = matrix(nrow=p+1, ncol = length(K))
for(j in 1:(length(K))) THETA[,j] = c(K[j],beta(K[j]))
THETA = t(THETA)
matplot(THETA[,c(1,1,1)], THETA[,-1], type = "l", xlab = "k", ylab = "coef", col=1)
abline(v=0.045,col="red")
k =0.045
xnew = cbind(1,x)
px = ncol(xnew)
fit = lsfit(rbind(xnew, sqrt(k)*diag(px)), c(training$Count, rep(0, px)), int=F)
ls.print(fit)
fit_model=lm(Count~Temp+Hum+Wind+Dew+Solar+Rain+Snow+Holiday+weekstatus+Rain:Seasons+Snow:weekstatus+Rain:Holiday+Seasons:Holiday,data=training)
summary(fit_model)
g1=update(fit_model,~.-Holiday:Seasons)
summary(g1)
g2=update(g1,~.-Rain:Holiday)
summary(g2)
g3=update(g2,~.-Snow:weekstatus)
summary(g3)
g4=update(g3,~.-Snow)
summary(g4)
g5=update(g4,~.-Wind)
summary(g5)
g6=update(g5,~.-Temp)
summary(g6)

#lasso regression
#install.packages("glmnet")
library(glmnet)
cv_model=cv.glmnet(x, training$Count,alpha = 1)#nfolds=10 default
plot(cv_model)
best_lambda=cv_model$lambda.min
best_model=glmnet(x, training$Count, alpha = 1, lambda = best_lambda)
coef(best_model)
best_model=lm(Count~ Hum+Wind+Visb+Dew+Solar+Rain+Snow+Holiday+weekstatus+Seasons+Rain:Seasons, data = training)
summary(best_model)
g1=update(best_model,~.-Visb)
summary(g1)
g2=update(g1,~.-Wind)
summary(g2)
g3=update(g2,~.-Snow)
summary(g3)

#check assumption for reduced model
qqnorm(back_model$residuals,main="Normal QQ plot")
qqline(back_model$residuals)
plot(back_model$fitted.values,back_model$residuals,xlab="fitted.value",ylab="residuals", main="fitted value vs residual")
abline(h=0)
abline(h=-18)
crPlots(main,terms = ~.-Temp-Wind-Visb-Snow-Rain-Seasons)

##prediction: use R^squre, cv-mse to accurcy
#install.packages("caret")
library(caret)
testing$Seasons=as.factor(testing$Seasons)
testing$Holiday=as.factor(testing$Holiday)
testing$weekstatus=as.factor(testing$weekstatus)
testing$Rain=as.factor(ifelse(testing$Rain==0, 0, 1))
testing$Snow=as.factor(ifelse(testing$Snow==0, 0, 1))
testing$Count=sqrt(testing$Count)
newy=testing$Count

#back elimination
back_predicted=predict(back_model,newdata=testing)
sst=sum((newy - mean(newy))^2)
sse_back=sum((back_predicted - newy)^2)
rsq_back=1-sse_back/sst;rsq_back
rjust_back=1- (1 -rsq_back)* ((106 - 1)/(106-length(back_model$coefficients)-1-1));rjust_back

cv.lm=function(data, nfolds = 5) {
  n=nrow(data)
  set.seed(34)
  fold.labels=sample(rep(1:nfolds, length.out = n)) 	
  mses=matrix(NA, nrow = nfolds, ncol = 1)
  for (i in 1:nfolds) {
    test.rows=which(fold.labels == i)
    train=data[-test.rows, ]
    test=data[test.rows, ]
    current.model=lm(Count ~ Hum + Dew+  Solar+  Rain+Holiday+Seasons+weekstatus+Rain:Seasons, data = train)
    predictions=predict(current.model, newdata = test)
    test.responses=test$Count
    test.errors=test.responses - predictions
    mses[i, 1]=mean(test.errors^2)}
return(apply(mses,2,mean))}
cv.lm(testing)

#cook distance and influence plot
cutoff=4/((nrow(training)-length(back_model$coefficients)))
cutoff 
plot(back_model, which=4,cook.levels=cutoff)
influencePlot(back_model, id.method="identify", main="Influence Plot", sub="Circle size is proportial to Cook's Distance" )

# sensiticity test
row_names_df_to_remove=c('2155','6523','3091','3019','7411','4795','8227')
training1 = training[!(row.names(training) %in% row_names_df_to_remove),]
h=lm(Count~Temp+Hum+Wind+Visb+Dew+Solar+Rain+Snow+Seasons+Holiday+weekstatus+Rain*Seasons+Snow*weekstatus+Rain*Holiday+Holiday*Seasons,data=training1)
l=step(h,direction = 'backward',trace=FALSE)
summary(l)
l1=update(l,~.-Snow)
summary(l1)
l1_predicted=predict(l1,newdata=testing)
sst=sum((newy - mean(newy))^2)
sse_l1=sum((l1_predicted - newy)^2)
rsq_l1=1-sse_l1/sst;rsq_l1
rjust_l1=1- (1 -rsq_l1)* ((106 - 1)/(106-length(l1$coefficients)-1-1));rjust_l1
