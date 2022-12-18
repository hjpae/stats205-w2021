### Stats 205, W2021 - Midterm (Hongju Pae) ###
## Problem 1: Tumor Count Data ##
Ya <- c(12, 9, 12, 14, 13, 13, 15, 8, 15, 6) # data
Yb <- c(11, 11, 10, 9, 9, 8, 7, 10, 6, 8, 8, 9, 7)
na <- length(Ya)
nb <- length(Yb)
Yasum <- sum(Ya)
Ybsum <- sum(Yb)
mean_a <- 12
mean_b <- Ybsum/nb

grid <- seq(0,20,1) 

theta_a_prior <- dgamma(grid, 120, 10)
theta_b_prior <- dgamma(grid, 12, 1)
  
theta_a_post  <- dgamma(grid, 120+Yasum, 10+na) 
theta_b_post  <- dgamma(grid, 12+Ybsum, 1+nb)


# theta a 
plot(grid, theta_a_prior, type="l", col="red", xlab="theta_a", ylab="density", ylim=c(0,1))
lines(grid, theta_a_post, col="blue")
legend("topright", c("Prior", "Posterior"), lwd=c(2,3),col=c("red", "blue"), inset=0.05)

mean_theta_a_post <- (120+Yasum) / (10+na) # mean 
var_theta_a_post <- var(theta_a_post) # variance 
q95_theta_a_post <- qgamma(c(0.025,0.975), 120+Yasum, 10+na) # 95% credible interval

theta_a <- c(mean_theta_a_post, var_theta_a_post, q95_theta_a_post)


# theta b
plot(grid, theta_b_prior, type="l", col="red", xlab="theta_b", ylab="density", ylim=c(0,1))
lines(grid, theta_b_post, col="blue")
legend("topright", c("Prior", "Posterior"), lwd=c(2,3),col=c("red", "blue"), inset=0.05)

mean_theta_b_post <- (120+Ybsum) / (10+nb) # mean
var_theta_b_post <- var(theta_b_post) # variance 
q95_theta_b_post <- qgamma(c(0.025,0.975), 12+Ybsum, 1+nb) # 95% credible interval

theta_b <- c(mean_theta_b_post, var_theta_b_post, q95_theta_b_post)


# theta a*b (joint probability distribution - need MCMC)
library(R2jags)

jags_model <- textConnection("model{
     y1 ~ dpois(theta1)
     theta1 ~ dgamma(120, 10)
     y2 ~ dpois(theta2)
     theta2 ~ dgamma(12, 1)
    joint.theta <- theta1 * theta2
}")

jags.data = list(y1 = 12, y2 = 11) # the first values from ya and yb data

jags.inits = function(){
  list("theta1" = 12, "theta2" = 8.6923) # 12 is from the question, 113/13 is from Yb/Ybsum (the mean) 
}

jags.param=c("theta1", "theta2", "joint.theta")

#jagsfit=jags(data = jags.data, n.chains = 3, inits = jags.inits,
#             parameters.to.save = jags.param, n.iter=2000, n.burnin=1000,
#             DIC = TRUE, model.file = textConnection(jags_model))

model <- jags.model(jags_model, data = jags.data, inits = jags.inits, n.chains=3, quiet=TRUE)
update(model, 10000, progress.bar="none")

samples <- coda.samples(model,
                        variable.names=jags.param,
                        n.iter=10000, progress.bar="none")

plot(samples)
summary(samples)

# conv diag
effectiveSize(samples)
gelman.diag(samples)

# data summary
mean_theta_ab_post <- 137.90 # mean
var_theta_ab_post <- 31.265^2 # variance 
q95_theta_ab_post <- as.numeric(c(84.722, 206.50)) # 95% credible interval

theta_ab <- c(mean_theta_ab_post, var_theta_ab_post, q95_theta_ab_post)

# draw the table 
theta <- t(data.frame(theta_a, theta_b, theta_ab))
colnames(theta) <- c("mean", "variance", "2.5%", "97.5%")
View(theta)

library('xtable')
xtable(theta)


## 1-3. 
# set i = 1 then run below before for loop
plot(grid, theta_prior, type="l", col="red", xlab="theta_b", ylab="density", ylim=c(0,1))
legend("topright", c("Prior", "Posterior"), lwd=c(2,3),col=c("red", "blue"), inset=0.05)

for(i in 1:50){
theta_prior <- dgamma(grid, 12*i, i)
theta_post  <- dgamma(grid, (12*i)+Ybsum, i+nb)
lines(grid, theta_prior, col="red")
lines(grid, theta_post, col="blue")
}

# set i = 50 then run below to mark i=50
lines(grid, theta_prior, col="green")
lines(grid, theta_post, col="green")


## Problem 2A: OASIS II cohort study ##
rm(list = ls())
#j = 5 #gamma distribution cannot be defined under support with 0, so dropped the trial with 0 value 

# heparin group (c)
yjc = c(2, 3, 4, 42, 4)
njc = c(122, 210, 105, 154, 70)
theta_jc = yjc / njc
lambda_jc = yjc / njc
yc = 213
nc = 5058
lambda_c = yc / nc

# placebo group (p)
yjp = c(4, 7, 9, 40, 7)
njp = c(121, 189, 109, 131, 73)
theta_jp = yjp / njp
lambda_jp = yjp / njp

# hirudin group (t) 
yt = 182
nt = 5083
lambda_t = yt / nt

# build hierarchical model using jags 
library(R2jags)
# poisson-gamma model 
jags.model = "model{
                for (i in 1:5){
                yjp[i] ~ dpois(lambda_jp[i]) # placebo group 
                lambda_jp[i] ~ dgamma(a[1], b[1])
                yjc[i] ~ dpois(lambda_jc[i]) # case-control heparin group
                lambda_jc[i] ~ dgamma(a[2], b[2])
            }
                yt ~ dpois(lambda_t) # test hirudin group 
                lambda_t ~ dgamma(a[3], b[3])
                yc ~ dpois(lambda_c) # case-control group (as whole)
                lambda_c ~ dgamma(a[4], b[4])
                

               # 2nd order hierarchy - hopefully towards gamma(12, 3)
                 for (j in 1:4){
                 a[j] ~ dnorm(12, 1)
                 b[j] ~ dnorm(3, 0.5) # distribution and hyperparams are arbitrary!!! 
            }
                 eff_cp <- lambda_jp - lambda_jc # risk difference (should be >1 if c is effective)
                 eff_tc <- lambda_c - lambda_t # risk difference (should be >1 if t is effective)
                 eff_tp <- lambda_jp - lambda_t 

                 prob_eff_cp <- step(lambda_jp - lambda_jc) # note is >0 if c is effective, not <0!!! 
                 prob_eff_tc <- step(lambda_c - lambda_t) # note is >0 if t is effective, not <0!!! 
                 prob_eff_tp <- step(lambda_jp - lambda_t)

                 rr_cp <- lambda_jp / lambda_jc # relative risk (should be >1 if c is effective)
                 rr_tc <- lambda_c / lambda_t # relative risk (should be >1 if t is effective)
                 rr_tp <- lambda_jp / lambda_t 

}"

jags.data = list(yjp = yjp, yjc = yjc, yt = yt, yc = yc, lambda_jp = lambda_jp, lambda_jc = lambda_jc, lambda_t = lambda_t, lambda_c = lambda_c) 
jags.param = c("prob_eff_tp", "prob_eff_tc", "rr_tp", "rr_tc", "eff_tp", "eff_tc", "lambda_t", "lambda_c", "lambda_jc", "lambda_jp")

jagsfit <- jags(data=jags.data, n.chains=5, inits=NULL,
                parameters=jags.param, n.iter=40000, n.burnin=10000,
                DIC=TRUE, model.file=textConnection(jags.model))
print(jagsfit)

# conv diag
jags.mcmc <- as.mcmc(jagsfit)
jags.sum <- summary(jags.mcmc)
jags.sum$statistics

library(bayesplot)
library(ggplot2)
color_scheme_set("blue")
mcmc_intervals(jags.mcmc, 
               pars=c("eff_tc", "eff_tp[1]", "eff_tp[2]", "eff_tp[3]", "eff_tp[4]", "eff_tp[5]", "lambda_c",  "lambda_jc[1]", "lambda_jc[2]", "lambda_jc[3]", "lambda_jc[4]", "lambda_jc[5]", "lambda_jp[1]", "lambda_jp[2]", "lambda_jp[3]", "lambda_jp[4]", "lambda_jp[5]", "lambda_t",  "prob_eff_tc",  "prob_eff_tp[1]", "prob_eff_tp[2]", "prob_eff_tp[3]", "prob_eff_tp[4]", "prob_eff_tp[5]", "rr_tc","rr_tp[1]", "rr_tp[2]", "rr_tp[3]", "rr_tp[4]", "rr_tp[5]"),
               prob = 0.8, # 80% intervals - inner
               prob_outer = 0.95, # 95% - outer
               point_est = "mean")

# beta-binomial model 
# jags.model = "model{
#                 for (i in 1:5){
#                 yjp[i] ~ dbin(theta_jp[i]) # placebo group 
#                 theta_jp[i] ~ dbeta(a[1], b[1])  
#                 yjc[i] ~ dbin(theta_jc[i]) # case-control heparin group
#                 theta_jc[i] ~ dbeta(a[2], b[2])
#             }

#                 # 2nd order hierarchy 
#                 for (j in 1:2){
#                 a[j] =eta[j]*mu[j]
#                 b[j] = eta[j]*(1-mu[j])
#                 eta[j] ~ dlnorm(0, 1/3); 
#                 mu[j] ~ dbeta(1, 1);   # distribution and hyperparams are arbitrary
#             }

#                 eff_cp_each <- step(theta_jp - theta_jc)
#                 eff_cp_overall <- step(mu[2] - mu[1])
# }"

#jags.data = list(yjp = yjp, yjc = yjc, theta_jp = theta_jp, theta_jc = theta_jc) 
#jags.param = c("theta_jp", "theta_jc", "eff_cp_each", "eff_cp_overall")

#jagsfit <- jags(data=jags.data, n.chains=5, inits=NULL,
#                parameters=jags.param, n.iter=2000, n.burnin=1000,
#                DIC=TRUE, model.file=textConnection(jags.model))
#print(jagsfit)


## Problem 2B: Ankylosing spondylitis II case-control study ##
rm(list = ls())
# library(epiR)
# epi.betabuster(0.332, 0.5, greaterthan, conf.level=0.35) 
# I used the toolbox http://252s-epi.vet.unimelb.edu.au:3838/epi/beta.buster/ to calculate hyperparams

# treatment group 
yt = 14
nt = 23
theta_t = yt/nt

# placebo group (historical trials) 
yp = c(23, 12, 19, 9, 39, 6, 9, 10)
np = c(107, 44, 51, 39, 139, 20, 78, 35) 
theta_p = yp/np 

library(R2jags)
# beta-binomial model
jags.model = "model{
                for (i in 1:8){
                yp[i] ~ dbin(theta_p[i], np[i]) # placebo group
                theta_p[i] ~ dbeta(a[1], b[1])
                }
                yt ~ dbin(theta_t, nt) # treatment group 
                theta_t ~ dbeta(a[2], b[2])

                # 2nd order hierarchy
                for (j in 1:2){
                a[j] =eta[j]*mu[j]
                b[j] = eta[j]*(1-mu[j])
                eta[j] ~ dlnorm(0, 1/3); # arbitrary value 
                mu[j] ~ dbeta(4.732, 8.509); # hyperparam from prior expectation (calculated by betabuster)
            }

                eff_each <- step(theta_t - theta_p) # should be >0 
                eff_overall <- step(mu[2] - mu[1])
                OR <- (theta_t / (1-theta_t)) / (theta_p / (1-theta_p))
                #eff <- theta_t - theta_p
}"

jags.data = list(yp = yp, yt = yt, np = np, nt = nt, theta_p = theta_p, theta_t = theta_t)
jags.param = c("theta_t", "theta_p", "eff_each", "eff_overall", "OR")

jagsfit <- jags(data=jags.data, n.chains=5, inits=NULL,
               parameters=jags.param, n.iter=20000, n.burnin=4000,
               DIC=TRUE, model.file=textConnection(jags.model))
print(jagsfit)

# conv diag
jags.mcmc <- as.mcmc(jagsfit)
jags.sum <- summary(jags.mcmc)
jags.sum$statistics

library(bayesplot)
library(ggplot2)
color_scheme_set("blue")
mcmc_intervals(jags.mcmc, 
               pars=c("eff_tp", "eff_cp[1]", "eff_cp[2]", "eff_cp[3]", "eff_cp[4]", "eff_cp[5]", "lambda_jc[1]", "lambda_jc[2]", "lambda_jc[3]", "lambda_jc[4]", "lambda_jc[5]", "lambda_jp[1]", "lambda_jp[2]", "lambda_jp[3]", "lambda_jp[4]", "lambda_jp[5]"),
               prob = 0.8, # 80% intervals - inner
               prob_outer = 0.95, # 95% - outer
               point_est = "mean")


## Problem 3: Pitcher's performance ##
rm(list = ls())

library(readr)
pitch <- read_csv("C:/Users/Kardien/Dropbox/UCI/2021 Winter Courseworks/Stats 205/midterm/pitch.csv")
View(pitch)

# we reduce the dataset to match the one needed for the analysis
library(dplyr)
pitch = pitch %>% select(-K9) %>%
  filter(Team=="BOS" | Team=="LAD"| Team=="LAA" |
           Team=="WSN" | Team=="SFG" | Team=="NYM"
         | Team=="NYY" | Team=="CIN" | Team=="ARI"
         | Team=="DET" | Team=="CLE" | Team=="HOU")

library(R2jags)
## 3-2. model propose: hierarchical normal model 
jags.model = "model{
              y_player ~ dnorm(mu, tau)
              y_teamyear ~ dpois(lambda)

              # 2nd order hierarchy
              mu ~ dnorm(0, 1)
              tau ~ dpois(lambda)

              # 3rd order hierarchy
              lambda ~ dgamma(1, 1)
}"


## 3-3. fitting the model - per players  
library(tidyverse)
pitch_team <- pitch %>% unite("play_team", c(Player, Team), remove = TRUE)
View(pitch_team)

pitch_2014 <- pitch_team %>% filter(Year == 2014)
pitch_2015 <- pitch_team %>% filter(Year == 2015)
pitch_2016 <- pitch_team %>% filter(Year == 2016)
pitch_2017 <- pitch_team %>% filter(Year == 2017)

library(R2jags)
# for all players 
jags.model = "model{
           for (i in 1:215){
              y_player[i] ~ dnorm(mu_p, tau_p)
           }

           for (a in 1:55){
              y_2014[a] ~ dnorm(mu_y[1], tau_y[1])
           }
           for (b in 1:54){
              y_2015[b] ~ dnorm(mu_y[2], tau_y[2])
           }
           for (c in 1:51){
              y_2016[c] ~ dnorm(mu_y[3], tau_y[3])
           }
           for (d in 1:55){
              y_2017[d] ~ dnorm(mu_y[4], tau_y[4])
           }

           for (j in 1:4){
              mu_y[j] ~ dnorm(0, 1)
              tau_y[j] ~ dpois(lambda)
           }
              mu_p ~ dnorm(0, 1)
              tau_p ~ dpois(lambda)
              lambda ~ dgamma(1, 1)
}"

jags.data = list(y_player = pitch_team$WHIP, 
                 y_2014 = pitch_2014$WHIP, 
                 y_2015 = pitch_2015$WHIP, 
                 y_2016 = pitch_2016$WHIP, 
                 y_2017 = pitch_2017$WHIP)

jags.param = c("mu_p", "tau_p", "mu_y", "tau_y", "lambda")

jagsfit <- jags(data=jags.data, n.chains=5, inits=NULL,
                parameters=jags.param, n.iter=20000, n.burnin=4000,
                DIC=TRUE, model.file=textConnection(jags.model))
print(jagsfit)

# for selective players
pitch_ck <- pitch_team %>% filter(play_team == "ClaytonKershaw_LAD")
pitch_hr <- pitch_team %>% filter(play_team == "Hyun-JinRyu_LAD")
pitch_mb <- pitch_team %>% filter(play_team == "MadisonBumgarner_SFG")
pitch_dk <- pitch_team %>% filter(play_team == "DallasKeuchel_HOU")

jags.model = "model{
for (m in 1:4){
y_ck[m] ~ dnorm(mu_ck, tau_ck)
}
for (n in 1:2){
y_hr[n] ~ dnorm(mu_hr, tau_hr)
}
for (o in 1:4){
y_mb[o] ~ dnorm(mu_mb, tau_mb)
}
for (p in 1:4){
y_dk[p] ~ dnorm(mu_dk, tau_dk)
}

for (a in 1:55){
y_2014[a] ~ dnorm(mu_y[1], tau_y[1])
}
for (b in 1:54){
y_2015[b] ~ dnorm(mu_y[2], tau_y[2])
}
for (c in 1:51){
y_2016[c] ~ dnorm(mu_y[3], tau_y[3])
}
for (d in 1:55){
y_2017[d] ~ dnorm(mu_y[4], tau_y[4])
}

for (j in 1:4){
mu_y[j] ~ dnorm(0, 1)
tau_y[j] ~ dpois(lambda)
}

mu_ck ~ dnorm(0, 1)
mu_hr ~ dnorm(0, 1)
mu_mb ~ dnorm(0, 1)
mu_dk ~ dnorm(0, 1)
tau_ck ~ dpois(lambda)
tau_hr ~ dpois(lambda)
tau_mb ~ dpois(lambda)
tau_dk ~ dpois(lambda)

lambda ~ dgamma(1, 1)
}"

jags.data = list(y_2014 = pitch_2014$WHIP, 
                 y_2015 = pitch_2015$WHIP, 
                 y_2016 = pitch_2016$WHIP, 
                 y_2017 = pitch_2017$WHIP, 
                 y_ck = pitch_ck$WHIP, 
                 y_hr = pitch_hr$WHIP, 
                 y_mb = pitch_mb$WHIP, 
                 y_dk = pitch_dk$WHIP)

jags.param = c("mu_ck", "mu_hr", "mu_mb", "mu_dk")

jagsfit <- jags(data=jags.data, n.chains=5, inits=NULL,
                parameters=jags.param, n.iter=20000, n.burnin=4000,
                DIC=TRUE, model.file=textConnection(jags.model))
print(jagsfit)

# conv diag
jags.mcmc <- as.mcmc(jagsfit)
jags.sum <- summary(jags.mcmc)
jags.sum$statistics

plot(jags.mcmc)

library(bayesplot)
library(ggplot2)
color_scheme_set("blue")
mcmc_intervals(jags.mcmc, 
               pars=c("mu_ck", "mu_hr", "mu_mb", "mu_dk"),
               prob_outer = 0.95, # 95% - outer
               point_est = "mean")


## 3-4. fitting the model - per team   

# for selected teams 
pitch_lad <- pitch %>% filter(Team == "LAD")
pitch_hou <- pitch %>% filter(Team == "HOU")
pitch_det <- pitch %>% filter(Team == "DET")

jags.model = "model{
for (i in 1:18){
y_lad[i] ~ dnorm(mu_lad, tau_lad)
}
for (j in 1:20){
y_hou[j] ~ dnorm(mu_hou, tau_hou)
}
for (k in 1:17){
y_det[k] ~ dnorm(mu_det, tau_det)
}

for (a in 1:55){
y_2014[a] ~ dnorm(mu_y[1], tau_y[1])
}
for (b in 1:54){
y_2015[b] ~ dnorm(mu_y[2], tau_y[2])
}
for (c in 1:51){
y_2016[c] ~ dnorm(mu_y[3], tau_y[3])
}
for (d in 1:55){
y_2017[d] ~ dnorm(mu_y[4], tau_y[4])
}

for (l in 1:4){
mu_y[l] ~ dnorm(0, 1)
tau_y[l] ~ dpois(lambda)
}

mu_lad ~ dnorm(0, 1)
mu_hou ~ dnorm(0, 1)
mu_det ~ dnorm(0, 1)

tau_lad ~ dpois(lambda)
tau_hou ~ dpois(lambda)
tau_det ~ dpois(lambda)

lambda ~ dgamma(1, 1)
}"

jags.data = list(y_2014 = pitch_2014$WHIP, 
                 y_2015 = pitch_2015$WHIP, 
                 y_2016 = pitch_2016$WHIP, 
                 y_2017 = pitch_2017$WHIP, 
                 y_lad = pitch_lad$WHIP, 
                 y_hou = pitch_hou$WHIP, 
                 y_det = pitch_det$WHIP)

jags.param = c("mu_lad", "mu_hou", "mu_det")

jagsfit <- jags(data=jags.data, n.chains=5, inits=NULL,
                parameters=jags.param, n.iter=20000, n.burnin=4000,
                DIC=TRUE, model.file=textConnection(jags.model))
print(jagsfit)

# conv diag
jags.mcmc <- as.mcmc(jagsfit)
jags.sum <- summary(jags.mcmc)
jags.sum$statistics

plot(jags.mcmc)

library(bayesplot)
library(ggplot2)
color_scheme_set("blue")
mcmc_intervals(jags.mcmc, 
               pars=c("mu_lad", "mu_hou", "mu_det"),
               prob_outer = 0.95, # 95% - outer
               point_est = "mean")

