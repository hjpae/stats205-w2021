### Stats 205, W2021 - Final (Hongju Pae) ###
### Heart and Estrogen/progestin Replacement Study (HERS) ###
# Data cleaning  -------------------------------------------------
rm(list = ls())
library(caret)
library(corrplot)
library(robustHD)
library(qwraps2)
library(xtable)
library(readr)

data <- read_csv("C:/Users/Kardien/Dropbox/UCI/2021 Winter Courseworks/Stats 205/final/hersdata_fixed.csv")
#min(data$LDL1[data$LDL1>=0], na.rm=TRUE) ## minus value LDL1 discarded (one subject with -20)
data$LDL1[data$LDL1<=0] <- NA
data <- na.omit(data)
attach(data)
data$ID <- NULL

plot(data)
cor(data)
corrplot(as.matrix(cor(data)), method="number")

data$age.st=standardize(age)
data$BMI.st=standardize(BMI)
data$LDL1.st=standardize(LDL1)
View(data)
attach(data)

# Exploratory DA  -------------------------------------------------
## EDA 
eda_list <-
  list("Age in years" =
         list("min"       = ~ min(age),
              "median"    = ~ median(age),
              "max"       = ~ max(age),
              "mean (sd)" = ~ qwraps2::mean_sd(age)),
       "BMI" =
         list("min"       = ~ min(BMI),
              "median"    = ~ median(BMI),
              "max"       = ~ max(BMI),
              "mean (sd)" = ~ qwraps2::mean_sd(BMI)),
       "Year 1 LDL level" =
         list("min"       = ~ min(LDL1),
              "median"    = ~ median(LDL1),
              "max"       = ~ max(LDL1),
              "mean (sd)" = ~ qwraps2::mean_sd(LDL1)),
       "Smoking status" =
         list("Yes" = ~ qwraps2::n_perc(smoking == 1),
              "No"  = ~ qwraps2::n_perc(smoking == 0)),
       "Drinking status" =
         list("Yes" = ~ qwraps2::n_perc(drinkany == 1),
              "No"  = ~ qwraps2::n_perc(drinkany == 0)),
       "Exercise status" =
         list("Yes" = ~ qwraps2::n_perc(exercise == 1),
              "No"  = ~ qwraps2::n_perc(exercise == 0)),
       "Diabetes status" =
         list("Yes" = ~ qwraps2::n_perc(diabetes == 1),
              "No"  = ~ qwraps2::n_perc(diabetes == 0)),
       "Hormone Therapy status" =
         list("Yes" = ~ qwraps2::n_perc(HT == 1),
              "No"  = ~ qwraps2::n_perc(HT == 0)),
       "Statin Therapy status" =
         list("Yes" = ~ qwraps2::n_perc(statins == 1),
              "No"  = ~ qwraps2::n_perc(statins == 0))
  )
eda_list <- summary_table(data, eda_list)
print(eda_list)

hist(LDL1, probability=TRUE, main="Year 1 LDL level", xlab="LDL1", ylab="density", ylim=c(0,0.013))
lines(density(data$LDL1), col="red")


# Simple LR  -------------------------------------------------
## models in Simple LR form 
# (Q4, 5) whole (no interactions) 
lm1 <- lm(LDL1.st ~ age.st + as.factor(smoking) + as.factor(drinkany) + as.factor(exercise) + as.factor(diabetes) + 
            BMI.st + as.factor(statins) + as.factor(HT), data = data)
summary(lm1)
car::vif(lm1) 

# (Q4, 5) no interactions except drink, exercise, diabetes 
lm2 <- lm(LDL1.st ~ age.st + as.factor(diabetes) + BMI.st + as.factor(statins) + as.factor(HT), data = data)
summary(lm2)
car::vif(lm2) 

# (Q6) HT*statin with lm2 
lm3 <- lm(LDL1.st ~ age.st + as.factor(diabetes) + BMI.st + as.factor(statins) + as.factor(HT) 
          + as.factor(statins)*as.factor(HT), data = data)
summary(lm3)
car::vif(lm3) 

# (Q7) BMI*statin with lm2 
lm4 <- lm(LDL1.st ~ age.st + as.factor(diabetes) + BMI.st + as.factor(statins) + as.factor(HT) 
          + as.factor(statins)*BMI.st, data = data)
summary(lm4)
car::vif(lm4) 

# (Q8) LDL~smoking (with lm2)
lm5 <- lm(LDL1.st ~ age.st + as.factor(diabetes) + BMI.st + as.factor(statins) + as.factor(HT)
          + as.factor(smoking), data = data)
summary(lm5)
car::vif(lm5) 

# (Q9, 10) LDL~ without diabetes (with lm2)
lm6 <- lm(LDL1.st ~ age.st + BMI.st + as.factor(statins) + as.factor(HT), data = data)
summary(lm6)
car::vif(lm6) 


# Bayesian LR - Main Effect only -------------------------------------------------
## Bayesian Linear Regression 
library(R2jags)

## lm1 (whole without interactions)
# use g-prior 
X.mat = model.matrix(LDL1.st ~ age.st + as.factor(smoking) + as.factor(drinkany) + as.factor(exercise) 
                     + as.factor(diabetes) + BMI.st + as.factor(statins) + as.factor(HT))
C0inv = t(X.mat)%*% X.mat # Precision matrix

jags.data=list(
  Y=LDL1.st,
  Xmat=X.mat,
  r=dim(X.mat)[2],
  n=dim(X.mat)[1],
  beta0=rep(0,dim(X.mat)[2]),
  C0inv=C0inv,
  gg=dim(X.mat)[1], # unit information prior
  a=0.001, b=0.001  # diffuse prior
)

reg.whole <- "model{
for(i in 1:n){
Y[i] ~ dnorm(mu[i],tau)
mu[i] <- beta[1] + beta[2]*Xmat[i,2] + beta[3]*Xmat[i,3] + beta[4]*Xmat[i,4] + beta[5]*Xmat[i,5] + beta[6]*Xmat[i,6] + beta[7]*Xmat[i,7] + beta[8]*Xmat[i,8] + beta[9]*Xmat[i,9]
# intercept, age, smoking, drink, exercise, diabetes, bmi, statin, ht
}
beta[1:r] ~ dmnorm(beta0[1:r],(tau/gg)*C0inv[1:r,1:r]) # g-prior
tau ~ dgamma(a,b)

# Estimate mean LDL for age=60, BMI=25.8, smoke=0, drink=0, diab=0, exercise=1, HT/statin varb
meanLDLht <- beta[1] + beta[2]*60 + beta[5] + beta[7]*25.8 + beta[9]  
meanLDLst <- beta[1] + beta[2]*60 + beta[5] + beta[7]*25.8 + beta[8]              
meanLDLboth <- beta[1] + beta[2]*60 + beta[5] + beta[7]*25.8 + beta[8] + beta[9]  
meanLDLnone <- beta[1] + beta[2]*60 + beta[5] + beta[7]*25.8 
meanLDLall <- beta[1] + beta[2]*60 + beta[3] + beta[5] + beta[6] + beta[7]*25.8   # with all bad health conditions and no trt

}"

jags.param=c("beta","tau","meanLDLht","meanLDLst","meanLDLboth","meanLDLnone")
jags.fit <- jags(data=jags.data, parameters.to.save = jags.param,
                 model.file=textConnection(reg.whole), n.iter=20000, n.chains=1,
                 n.burnin=5000, n.thin=1, DIC=T, digits=6)
print(jags.fit)

library(bayesplot)
library(ggplot2)
jags.mcmc=as.mcmc(jags.fit)
color_scheme_set("blue")
mcmc_intervals(jags.mcmc, 
               pars=c("meanLDLht","meanLDLst","meanLDLboth","meanLDLnone"),
               prob_outer = 0.95, # 95% - outer
               point_est = "mean")


## lm1 (whole without interactions)
# use vague independent prior 
X.mat = model.matrix(LDL1.st ~ age.st + as.factor(smoking) + as.factor(drinkany) + as.factor(exercise) 
                     + as.factor(diabetes) + BMI.st + as.factor(statins) + as.factor(HT))

jags.data=list(
  Y=LDL1.st, 
  Xmat=X.mat,
  r=dim(X.mat)[2], 
  n=dim(X.mat)[1],
  binv=1 * diag(1, dim(X.mat)[2]),
  c=0.001
)

model_reg <- "model{
for(i in 1:n){
Y[i] ~ dnorm(mu[i],tau)
mu[i] <- beta[1] + beta[2]*Xmat[i,2] + beta[3]*Xmat[i,3] + beta[4]*Xmat[i,4] + beta[5]*Xmat[i,5] + beta[6]*Xmat[i,6] + beta[7]*Xmat[i,7] + beta[8]*Xmat[i,8] + beta[9]*Xmat[i,9]
# intercept, age, smoking, drink, exercise, diabetes, bmi, statin, ht
}
beta[1:r] ~ dmnorm(rep(0,r),binv) # precision matrix
tau ~ dgamma(c,c)

# Estimate mean LDL for age=60, BMI=25.8, smoke=0, drink=0, diab=0, exercise=1, HT/statin varb
meanLDLht <- beta[1] + beta[2]*60 + beta[5] + beta[7]*25.8 + beta[9]  
meanLDLst <- beta[1] + beta[2]*60 + beta[5] + beta[7]*25.8 + beta[8]              
meanLDLboth <- beta[1] + beta[2]*60 + beta[5] + beta[7]*25.8 + beta[8] + beta[9]  
meanLDLnone <- beta[1] + beta[2]*60 + beta[5] + beta[7]*25.8 
}"

jags.param=c("beta","tau","meanLDLht","meanLDLst","meanLDLboth","meanLDLnone")
jags.fit.vip1 <- jags(data=jags.data, parameters.to.save = jags.param,
                 model.file=textConnection(model_reg), n.iter=20000, n.chains=1,
                 n.burnin=5000, n.thin=1, DIC=T, digits=6)
print(jags.fit.vip1)


## lm2 (without some predictors)
# use g-prior 
X.mat = model.matrix(LDL1.st ~ age.st + as.factor(diabetes) + BMI.st + as.factor(statins) + as.factor(HT))
C0inv = t(X.mat)%*% X.mat # Precision matrix

jags.data=list(
  Y=LDL1.st,
  Xmat=X.mat,
  r=dim(X.mat)[2],
  n=dim(X.mat)[1],
  beta0=rep(0,dim(X.mat)[2]),
  C0inv=C0inv,
  gg=dim(X.mat)[1], # unit information prior
  a=0.001, b=0.001  # diffuse prior
)

reg.whole2 <- "model{
for(i in 1:n){
Y[i] ~ dnorm(mu[i],tau)
mu[i] <- beta[1] + beta[2]*Xmat[i,2] + beta[3]*Xmat[i,3] + beta[4]*Xmat[i,4] + beta[5]*Xmat[i,5] + beta[6]*Xmat[i,6] 
# intercept, age, diabetes, bmi, statin, ht

# for LPML comparison (Q9)
like[i] <- dnorm(Y[i], mu[i],tau) 
loginvlike[i] <- -log(like[i]) # log(1/like)=log(1)-log(like)
invlike[i] <- exp(loginvlike[i]) # take the exponent of what you have above
}

beta[1:r] ~ dmnorm(beta0[1:r],(tau/gg)*C0inv[1:r,1:r]) # g-prior
tau ~ dgamma(a,b)

# Estimate mean LDL for age=60, BMI=25.8, smoke=0, drink=0, diab=0, exercise=0, HT/statin varb
meanLDLht <- beta[1] + beta[2]*60 + beta[4]*25.8 + beta[6]  
meanLDLst <- beta[1] + beta[2]*60 + beta[4]*25.8 + beta[5]              
meanLDLboth <- beta[1] + beta[2]*60 + beta[4]*25.8 + beta[5] + beta[6]  
meanLDLnone <- beta[1] + beta[2]*60 + beta[4]*25.8 
meanLDLall <- beta[1] + beta[2]*60 + beta[3] + beta[4]*25.8   # with all bad health conditions and no trt

}"

jags.param=c("beta","tau","meanLDLht","meanLDLst","meanLDLboth","meanLDLnone","meanLDLall","invlike")
jags.fit.whole2 <- jags(data=jags.data, parameters.to.save = jags.param,
                 model.file=textConnection(reg.whole2), n.iter=20000, n.chains=1,
                 n.burnin=5000, n.thin=1, DIC=T, digits=6)
print(jags.fit.whole2)

# plot the posterior - predictions and regerssion coefficient (Q5)
library(bayesplot)
library(ggplot2)
jags.mcmc=as.mcmc(jags.fit.whole2)
color_scheme_set("blue")

# predictions 
mcmc_intervals(jags.mcmc, 
               pars=c("meanLDLht","meanLDLst","meanLDLboth","meanLDLnone"),
               prob_outer = 0.95, # 95% - outer
               point_est = "mean"
               )

# regression coefficient 
mcmc_intervals(jags.mcmc, 
               pars=c("beta[2]","beta[3]","beta[4]","beta[5]","beta[6]"),
               prob_outer = 0.95, # 95% - outer
               point_est = "mean"
)

# LPML (Q9) 
CPO.diab <- 1/jags.fit.whole2$BUGSoutput$mean$invlike ## invlike is a vector of length n
LPML.diab <- sum(log(CPO.diab))
print(LPML.diab)


## lm2 (whole without interactions)
# use vague independent prior 
X.mat = model.matrix(LDL1.st ~ age.st + as.factor(diabetes) + BMI.st + as.factor(statins) + as.factor(HT))
jags.data=list(
  Y=LDL1.st, 
  Xmat=X.mat,
  r=dim(X.mat)[2], 
  n=dim(X.mat)[1],
  binv=1 * diag(1, dim(X.mat)[2]), # b value
  c=0.001
)

model_reg <- "model{
for(i in 1:n){
Y[i] ~ dnorm(mu[i],tau)
mu[i] <- beta[1] + beta[2]*Xmat[i,2] + beta[3]*Xmat[i,3] + beta[4]*Xmat[i,4] + beta[5]*Xmat[i,5] + beta[6]*Xmat[i,6] 
# intercept, age, diabetes, bmi, statin, ht
}
beta[1:r] ~ dmnorm(rep(0,r),binv) # precision matrix
tau ~ dgamma(c,c)

# Estimate mean LDL for age=60, BMI=25.8, smoke=0, drink=0, diab=0, exercise=0, HT/statin varb
meanLDLht <- beta[1] + beta[2]*60 + beta[4]*25.8 + beta[6]  
meanLDLst <- beta[1] + beta[2]*60 + beta[4]*25.8 + beta[5]              
meanLDLboth <- beta[1] + beta[2]*60 + beta[4]*25.8 + beta[5] + beta[6]  
meanLDLnone <- beta[1] + beta[2]*60 + beta[4]*25.8 
}"

jags.param=c("beta","tau","meanLDLht","meanLDLst","meanLDLboth","meanLDLnone")
jags.fit.vip2 <- jags(data=jags.data, parameters.to.save = jags.param,
                      model.file=textConnection(model_reg), n.iter=20000, n.chains=1,
                      n.burnin=5000, n.thin=1, DIC=T, digits=6)
print(jags.fit.vip2)

# draw the DIC comparison table (Q4)
comparison <- t(data.frame(c(7234.4, 7234.3, 7230.6, 7230.7)))
colnames(comparison) <- c("M1-G", "M1-VIP", "M2-G", "M2-VIP")
View(comparison)
xtable(comparison)



# Bayesian LR - Interactions -------------------------------------------------
## now will keep using g prior and lm2 as base model 
## HT*statin with lm2 (Q6) 
X.mat = model.matrix(LDL1.st ~ age.st + as.factor(diabetes) + BMI.st + as.factor(statins) + as.factor(HT)
                     + as.factor(statins)*as.factor(HT))
C0inv = t(X.mat)%*% X.mat # Precision matrix

jags.data=list(
  Y=LDL1.st,
  Xmat=X.mat,
  r=dim(X.mat)[2],
  n=dim(X.mat)[1],
  beta0=rep(0,dim(X.mat)[2]),
  C0inv=C0inv,
  gg=dim(X.mat)[1], # unit information prior
  a=0.001, b=0.001  # diffuse prior
)

jags.reg <- "model{
for(i in 1:n){
Y[i] ~ dnorm(mu[i],tau)
mu[i] <- beta[1] + beta[2]*Xmat[i,2] + beta[3]*Xmat[i,3] + beta[4]*Xmat[i,4] + beta[5]*Xmat[i,5] + beta[6]*Xmat[i,6] + beta[7]*Xmat[i,7]
# intercept, age, diabetes, bmi, statin, ht, statin*ht
}
beta[1:r] ~ dmnorm(beta0[1:r],(tau/gg)*C0inv[1:r,1:r]) # g-prior
tau ~ dgamma(a,b)

# Estimate mean LDL for age=60, BMI=25.8, smoke=0, drink=0, diab=0, exercise=0, HT/statin varb
meanLDLht <- beta[1] + beta[2]*60 + beta[4]*25.8 + beta[6] + beta[7]
meanLDLst <- beta[1] + beta[2]*60 + beta[4]*25.8 + beta[5] + beta[7]    
meanLDLboth <- beta[1] + beta[2]*60 + beta[4]*25.8 + beta[5] + beta[6] + beta[7]
meanLDLnone <- beta[1] + beta[2]*60 + beta[4]*25.8 
}"

jags.param=c("beta","tau","meanLDLht","meanLDLst","meanLDLboth","meanLDLnone")
jags.fit <- jags(data=jags.data, parameters.to.save = jags.param,
                 model.file=textConnection(jags.reg), n.iter=20000, n.chains=1,
                 n.burnin=5000, n.thin=1, DIC=T, digits=6)
print(jags.fit)

# plot the posterior - predictions and regerssion coefficient
library(bayesplot)
library(ggplot2)
jags.mcmc=as.mcmc(jags.fit)
color_scheme_set("blue")

# predictions 
mcmc_intervals(jags.mcmc, 
               pars=c("meanLDLht","meanLDLst","meanLDLboth","meanLDLnone"),
               prob_outer = 0.95, # 95% - outer
               point_est = "mean"
)

# regression coefficient 
mcmc_intervals(jags.mcmc, 
               pars=c("beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]"),
               prob_outer = 0.95, # 95% - outer
               point_est = "mean"
)


## BMI*statin with lm2 (Q7) 
X.mat = model.matrix(LDL1.st ~ age.st + as.factor(diabetes) + BMI.st + as.factor(statins) + as.factor(HT)
                     + as.factor(statins)*BMI.st)
C0inv = t(X.mat)%*% X.mat # Precision matrix

jags.data=list(
  Y=LDL1.st,
  Xmat=X.mat,
  r=dim(X.mat)[2],
  n=dim(X.mat)[1],
  beta0=rep(0,dim(X.mat)[2]),
  C0inv=C0inv,
  gg=dim(X.mat)[1], # unit information prior
  a=0.001, b=0.001  # diffuse prior
)

jags.reg <- "model{
for(i in 1:n){
Y[i] ~ dnorm(mu[i],tau)
mu[i] <- beta[1] + beta[2]*Xmat[i,2] + beta[3]*Xmat[i,3] + beta[4]*Xmat[i,4] + beta[5]*Xmat[i,5] + beta[6]*Xmat[i,6] + beta[7]*Xmat[i,7]
# intercept, age, diabetes, bmi, statin, ht, statin*BMI
}
beta[1:r] ~ dmnorm(beta0[1:r],(tau/gg)*C0inv[1:r,1:r]) # g-prior
tau ~ dgamma(a,b)

# Estimate mean LDL for age=60, BMI=25.8, smoke=0, drink=0, diab=0, exercise=0, HT/statin varb
meanLDLht <- beta[1] + beta[2]*60 + beta[4]*25.8 + beta[6]
meanLDLst <- beta[1] + beta[2]*60 + beta[4]*25.8 + beta[5] + beta[7]*25.8   
meanLDLboth <- beta[1] + beta[2]*60 + beta[4]*25.8 + beta[5] + beta[6] + beta[7]*25.8 
meanLDLnone <- beta[1] + beta[2]*60 + beta[4]*25.8 
}"

jags.param=c("beta","tau","meanLDLht","meanLDLst","meanLDLboth","meanLDLnone")
jags.fit <- jags(data=jags.data, parameters.to.save = jags.param,
                 model.file=textConnection(jags.reg), n.iter=20000, n.chains=1,
                 n.burnin=5000, n.thin=1, DIC=T, digits=6)
print(jags.fit)

# plot the posterior - predictions and regerssion coefficient
library(bayesplot)
library(ggplot2)
jags.mcmc=as.mcmc(jags.fit)
color_scheme_set("blue")

# predictions 
mcmc_intervals(jags.mcmc, 
               pars=c("meanLDLht","meanLDLst","meanLDLboth","meanLDLnone"),
               prob_outer = 0.95, # 95% - outer
               point_est = "mean"
)

# regression coefficient 
mcmc_intervals(jags.mcmc, 
               pars=c("beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]"),
               prob_outer = 0.95, # 95% - outer
               point_est = "mean"
)


# Bayesian LR - Others -------------------------------------------------
## smoking with lm2 
X.mat = model.matrix(LDL1.st ~ age.st + as.factor(diabetes) + BMI.st + as.factor(statins) + as.factor(HT)
                     + as.factor(smoking))
C0inv = t(X.mat)%*% X.mat # Precision matrix

jags.data=list(
  Y=LDL1.st,
  Xmat=X.mat,
  r=dim(X.mat)[2],
  n=dim(X.mat)[1],
  beta0=rep(0,dim(X.mat)[2]),
  C0inv=C0inv,
  gg=dim(X.mat)[1], # unit information prior
  a=0.001, b=0.001  # diffuse prior
)

jags.reg <- "model{
for(i in 1:n){
Y[i] ~ dnorm(mu[i],tau)
mu[i] <- beta[1] + beta[2]*Xmat[i,2] + beta[3]*Xmat[i,3] + beta[4]*Xmat[i,4] + beta[5]*Xmat[i,5] + beta[6]*Xmat[i,6] + beta[7]*Xmat[i,7]
# intercept, age, diabetes, bmi, statin, ht, smoking
}
beta[1:r] ~ dmnorm(beta0[1:r],(tau/gg)*C0inv[1:r,1:r]) # g-prior
tau ~ dgamma(a,b)

# Estimate mean LDL for age=60, BMI=25.8, smoke=1, drink=0, diab=0, exercise=0, HT/statin varb
meanLDLhsm <- beta[1] + beta[2]*60 + beta[4]*25.8 + beta[6] + beta[7]
meanLDLssm <- beta[1] + beta[2]*60 + beta[4]*25.8 + beta[5] + beta[7]  
meanLDLhnsm <- beta[1] + beta[2]*60 + beta[4]*25.8 + beta[6]
meanLDLsnsm <- beta[1] + beta[2]*60 + beta[4]*25.8 + beta[5]
}"

jags.param=c("beta","tau","meanLDLhsm","meanLDLssm","meanLDLhnsm","meanLDLsnsm")
jags.fit <- jags(data=jags.data, parameters.to.save = jags.param,
                 model.file=textConnection(jags.reg), n.iter=20000, n.chains=1,
                 n.burnin=5000, n.thin=1, DIC=T, digits=6)
print(jags.fit)

# plot the posterior - predictions and regerssion coefficient
library(bayesplot)
library(ggplot2)
jags.mcmc=as.mcmc(jags.fit)
color_scheme_set("blue")

# predictions 
mcmc_intervals(jags.mcmc, 
               pars=c("meanLDLhsm","meanLDLssm","meanLDLhnsm","meanLDLsnsm"),
               prob_outer = 0.95, # 95% - outer
               point_est = "mean"
)

# regression coefficient 
mcmc_intervals(jags.mcmc, 
               pars=c("beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]"),
               prob_outer = 0.95, # 95% - outer
               point_est = "mean"
)


## WITHOUT diabetes with lm2 
X.mat = model.matrix(LDL1.st ~ age.st + BMI.st + as.factor(statins) + as.factor(HT))
C0inv = t(X.mat)%*% X.mat # Precision matrix

jags.data=list(
  Y=LDL1.st,
  Xmat=X.mat,
  r=dim(X.mat)[2],
  n=dim(X.mat)[1],
  beta0=rep(0,dim(X.mat)[2]),
  C0inv=C0inv,
  gg=dim(X.mat)[1], # unit information prior
  a=0.001, b=0.001  # diffuse prior
)

jags.reg <- "model{
for(i in 1:n){
Y[i] ~ dnorm(mu[i],tau)
mu[i] <- beta[1] + beta[2]*Xmat[i,2] + beta[3]*Xmat[i,3] + beta[4]*Xmat[i,4] + beta[5]*Xmat[i,5]
# intercept, age, bmi, statin, ht

# for LPML comparison (Q9)
like[i] <- dnorm(Y[i], mu[i], tau) 
loginvlike[i] <- -log(like[i]) # log(1/like)=log(1)-log(like)
invlike[i] <- exp(loginvlike[i]) # take the exponent of what you have above
}

beta[1:r] ~ dmnorm(beta0[1:r],(tau/gg)*C0inv[1:r,1:r]) # g-prior
tau ~ dgamma(a,b)

# Estimate mean LDL for age=60, BMI=25.8, smoke=0, drink=0, diab=0, exercise=0, HT/statin varb
meanLDLht <- beta[1] + beta[2]*60 + beta[3]*25.8 + beta[5]  
meanLDLst <- beta[1] + beta[2]*60 + beta[3]*25.8 + beta[4]              
meanLDLboth <- beta[1] + beta[2]*60 + beta[3]*25.8 + beta[4] + beta[5]  
meanLDLnone <- beta[1] + beta[2]*60 + beta[3]*25.8 
}"

jags.param=c("beta","tau","meanLDLht","meanLDLst","meanLDLboth","meanLDLnone","invlike")
jags.fit <- jags(data=jags.data, parameters.to.save = jags.param,
                 model.file=textConnection(jags.reg), n.iter=20000, n.chains=1,
                 n.burnin=5000, n.thin=1, DIC=T, digits=6)
print(jags.fit)

# plot the posterior - predictions and regerssion coefficient
library(bayesplot)
library(ggplot2)
jags.mcmc=as.mcmc(jags.fit)
color_scheme_set("blue")

# predictions 
mcmc_intervals(jags.mcmc, 
               pars=c("meanLDLht","meanLDLst","meanLDLboth","meanLDLnone"),
               prob_outer = 0.95, # 95% - outer
               point_est = "mean"
)

# regression coefficient 
mcmc_intervals(jags.mcmc, 
               pars=c("beta[2]","beta[3]","beta[4]","beta[5]"),
               prob_outer = 0.95, # 95% - outer
               point_est = "mean"
)

# LPML (Q9) 
CPO <- 1/jags.fit$BUGSoutput$mean$invlike ## invlike is a vector of length n
LPML.without.diab <- sum(log(CPO))
print(LPML.without.diab)

# (p)BF (Q10)
PBF <- exp(LPML.diab - LPML.without.diab)
print(PBF)

# draw the table (Q9, 10)
assess <- t(data.frame(c(LPML.diab, LPML.without.diab, PBF, exp(LPML.without.diab - LPML.diab))))
colnames(assess) <- c("LPML.diab", "LPML.nodiab", "PBF.diab-nodiab", "PBF.nodiab-diab")
View(assess)
xtable(assess)


# Bayesian LR - Sensitivity analysis (G-prior & indp BCJ prior) (Q11) -----------------
## use lm2 - Model 2
# g-prior ... when g = n
X.mat = model.matrix(LDL1.st ~ age.st + as.factor(diabetes) + BMI.st + as.factor(statins) + as.factor(HT))
C0inv = t(X.mat)%*% X.mat # Precision matrix

jags.data=list(
  Y=LDL1.st,
  Xmat=X.mat,
  r=dim(X.mat)[2],
  n=dim(X.mat)[1],
  beta0=rep(0,dim(X.mat)[2]),
  C0inv=C0inv,
  gg=dim(X.mat)[1], # g = n
  a=0.001, b=0.001  # diffuse prior
)

model_reg <- "model{
for(i in 1:n){
Y[i] ~ dnorm(mu[i],tau)
mu[i] <- beta[1] + beta[2]*Xmat[i,2] + beta[3]*Xmat[i,3] + beta[4]*Xmat[i,4] + beta[5]*Xmat[i,5] + beta[6]*Xmat[i,6] 
# intercept, age, diabetes, bmi, statin, ht
}
beta[1:r] ~ dmnorm(beta0[1:r],(tau/gg)*C0inv[1:r,1:r]) # g-prior
tau ~ dgamma(a,b)
}"

jags.param=c("beta","tau")
jags.fit <- jags(data=jags.data, parameters.to.save = jags.param,
                        model.file=textConnection(model_reg), n.iter=20000, n.chains=1,
                        n.burnin=5000, n.thin=1, DIC=T, digits=6)
print(jags.fit)

# plot the posterior 
library(bayesplot)
library(ggplot2)
jags.mcmc=as.mcmc(jags.fit)
color_scheme_set("blue")

# regression coefficient 
mcmc_intervals(jags.mcmc, 
               pars=c("beta[2]","beta[3]","beta[4]","beta[5]","beta[6]"),
               prob_outer = 0.95, # 95% - outer
               point_est = "mean")


# g-prior ... when g = 1
X.mat = model.matrix(LDL1.st ~ age.st + as.factor(diabetes) + BMI.st + as.factor(statins) + as.factor(HT))
C0inv = t(X.mat)%*% X.mat # Precision matrix

jags.data=list(
  Y=LDL1.st,
  Xmat=X.mat,
  r=dim(X.mat)[2],
  n=dim(X.mat)[1],
  beta0=rep(0,dim(X.mat)[2]),
  C0inv=C0inv,
  gg=1, # g = 1
  a=0.001, b=0.001  # diffuse prior
)

model_reg <- "model{
for(i in 1:n){
Y[i] ~ dnorm(mu[i],tau)
mu[i] <- beta[1] + beta[2]*Xmat[i,2] + beta[3]*Xmat[i,3] + beta[4]*Xmat[i,4] + beta[5]*Xmat[i,5] + beta[6]*Xmat[i,6] 
# intercept, age, diabetes, bmi, statin, ht
}
beta[1:r] ~ dmnorm(beta0[1:r],(tau/gg)*C0inv[1:r,1:r]) # g-prior
tau ~ dgamma(a,b)
}"

jags.param=c("beta","tau")
jags.fit <- jags(data=jags.data, parameters.to.save = jags.param,
                 model.file=textConnection(model_reg), n.iter=20000, n.chains=1,
                 n.burnin=5000, n.thin=1, DIC=T, digits=6)
print(jags.fit)

# plot the posterior 
library(bayesplot)
library(ggplot2)
jags.mcmc=as.mcmc(jags.fit)
color_scheme_set("blue")

# regression coefficient 
mcmc_intervals(jags.mcmc, 
               pars=c("beta[2]","beta[3]","beta[4]","beta[5]","beta[6]"),
               prob_outer = 0.95, # 95% - outer
               point_est = "mean")


## Model 2 with VIP 
# use vague independent prior, b=1
X.mat = model.matrix(LDL1.st ~ age.st + as.factor(diabetes) + BMI.st + as.factor(statins) + as.factor(HT))
jags.data=list(
  Y=LDL1.st, 
  Xmat=X.mat,
  r=dim(X.mat)[2], 
  n=dim(X.mat)[1],
  binv=1 * diag(1, dim(X.mat)[2]), # b = 1
  c=0.001
)

model_reg <- "model{
for(i in 1:n){
Y[i] ~ dnorm(mu[i],tau)
mu[i] <- beta[1] + beta[2]*Xmat[i,2] + beta[3]*Xmat[i,3] + beta[4]*Xmat[i,4] + beta[5]*Xmat[i,5] + beta[6]*Xmat[i,6] 
# intercept, age, diabetes, bmi, statin, ht
}
beta[1:r] ~ dmnorm(rep(0,r),binv) # precision matrix
tau ~ dgamma(c,c)
}"

jags.param=c("beta","tau")
jags.fit <- jags(data=jags.data, parameters.to.save = jags.param,
                      model.file=textConnection(model_reg), n.iter=20000, n.chains=1,
                      n.burnin=5000, n.thin=1, DIC=T, digits=6)
print(jags.fit)

# plot the posterior 
library(bayesplot)
library(ggplot2)
jags.mcmc=as.mcmc(jags.fit)
color_scheme_set("blue")

# regression coefficient 
mcmc_intervals(jags.mcmc, 
               pars=c("beta[2]","beta[3]","beta[4]","beta[5]","beta[6]"),
               prob_outer = 0.95, # 95% - outer
               point_est = "mean")


# use vague independent prior 
X.mat = model.matrix(LDL1.st ~ age.st + as.factor(diabetes) + BMI.st + as.factor(statins) + as.factor(HT))
jags.data=list(
  Y=LDL1.st, 
  Xmat=X.mat,
  r=dim(X.mat)[2], 
  n=dim(X.mat)[1],
  binv=100 * diag(1, dim(X.mat)[2]), # b = 100
  c=0.001
)

model_reg <- "model{
for(i in 1:n){
Y[i] ~ dnorm(mu[i],tau)
mu[i] <- beta[1] + beta[2]*Xmat[i,2] + beta[3]*Xmat[i,3] + beta[4]*Xmat[i,4] + beta[5]*Xmat[i,5] + beta[6]*Xmat[i,6] 
# intercept, age, diabetes, bmi, statin, ht
}
beta[1:r] ~ dmnorm(rep(0,r),binv) # precision matrix
tau ~ dgamma(c,c)
}"

jags.param=c("beta","tau")
jags.fit <- jags(data=jags.data, parameters.to.save = jags.param,
                 model.file=textConnection(model_reg), n.iter=20000, n.chains=1,
                 n.burnin=5000, n.thin=1, DIC=T, digits=6)
print(jags.fit)

# plot the posterior 
library(bayesplot)
library(ggplot2)
jags.mcmc=as.mcmc(jags.fit)
color_scheme_set("blue")

# regression coefficient 
mcmc_intervals(jags.mcmc, 
               pars=c("beta[2]","beta[3]","beta[4]","beta[5]","beta[6]"),
               prob_outer = 0.95, # 95% - outer
               point_est = "mean")


# sensitivity analysis - draw table (Q11)
gprior <- read_delim("C:/Users/Kardien/Dropbox/UCI/2021 Winter Courseworks/Stats 205/final/gprior.csv", 
                     " ", escape_double = FALSE, trim_ws = TRUE)
vip <- read_delim("C:/Users/Kardien/Dropbox/UCI/2021 Winter Courseworks/Stats 205/final/vip.csv", 
                  " ", escape_double = FALSE, trim_ws = TRUE)
xtable(gprior)
xtable(vip)

