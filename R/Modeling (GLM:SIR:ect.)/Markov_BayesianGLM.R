##### Question 1: Bayesian GLM #####
#1.

#a
install.packages("titanic")
library(titanic)
data(titanic_train)
library(rjags)

## Specify the data ##
y <- titanic_train$Survived
x <- titanic_train[,colnames(titanic_train) ==  "Pclass" | colnames(titanic_train) == "Fare" | colnames(titanic_train) == "Age" | colnames(titanic_train) == "Sex" | colnames(titanic_train) == "Sibsp" | colnames(titanic_train) == "Parch"]
x$Sex[x$Sex=="female"] <- 0
x$Sex[x$Sex=="male"] <- 1

#b
#Bernoulli distributed
#dbern(p)

#c
x1 <- x$Pclass
NxLev = length(levels(as.factor(x1)))

dataList = list(
  x1 = x1,
  NxLev = NxLev,
  y = y 
)
library(rjags)
library(coda)

#d
## Specify the model ##
model_string <- "model {
for (i in 1:length(y)){
## Part 1
y[i] ~ dbern(p[i])
logit(p[i]) <- beta0 + a[x1[i]]
}

beta0 ~ dnorm(0,0.001)

for(j in 1:NxLev){
a[j] ~ dnorm(0,0.001)
}

## Part 2
#Convert beta0, a[] to sum-to-zero b0, b[]
for(j in 1:NxLev){
m[j] <- beta0 + a[j]} #cell means
b0 <- mean(m[1:NxLev])
for(j in 1:NxLev){
b[j] <- m[j] - b0}
}"
writeLines(model_string,con="TEMPmodel.txt")

## Prepare JAGS for sampling from the posterior distribution ##
#Give initial values
initsList = list(
  beta0 = mean(y),
  a = rep(0,NxLev)
)

#Set MCMC parameters
parameters = c("b0", "b") # parameters to record
adaptSteps = 500  # Number of steps to "tune" the samplers. See '?adapt'
burnInSteps = 1000 # number of steps to allow the Markov chains to run before considering a reliable sample from the posterior distribution
nChains = 3 # number of parallel Markov chains to run
Iter = 5000 # number of iterations of the Markov chain (each of which provides a sample of the posterior distribution)

#Record information into JAGS
jagsModel = jags.model("TEMPmodel.txt", data=dataList, inits=initsList, n.chains=nChains, n.adapt=adaptSteps)

#Burn-in period
update(jagsModel, n.iter=burnInSteps)

#Run and record MCMC samples
codaSamples = coda.samples(jagsModel, variable.names=parameters, n.iter=Iter)

#e
# traceplot and gelman.plot are from CODA package:
traceplot(codaSamples)
traceplot(codaSamples[[1]][,"b[1]"])
traceplot(codaSamples[[2]][,"b[1]"],col="skyblue",add=T)
traceplot(codaSamples[[3]][,"b[1]"],col="gray",add=T)

#f
# gelman plot
gelman.plot(codaSamples)

#g
# cool density plot of parameter values for all chains
HDIofMCMC <- dget("HDIofMCMC.R")
DbdaDensPlot <- dget("DbdaDensPlot.R")
DBDAplColors = c("black","skyblue","gray")
DbdaDensPlot(codaSamples,"b0",plColors=DBDAplColors)

#h
# Interpret the output
logit_to_prob <- function(x){
  ret <- exp(x) / (1+exp(x))
  return(ret)
}

#i
jag_info <- data.frame(rbind(codaSamples[[1]],codaSamples[[2]],codaSamples[[3]]))
jag_info$ULb0 <- logit_to_prob(jag_info$b0)
jag_info$ULb1 <- logit_to_prob(jag_info$b0 + jag_info$b.1.)
jag_info$ULb2 <- logit_to_prob(jag_info$b0 + jag_info$b.2.)
jag_info$ULb3 <- logit_to_prob(jag_info$b0 + jag_info$b.3.)

#j
# HDI
plotPost <- dget("plotPost.R")
plotPost(jag_info$ULb0)
plotPost(jag_info$ULb1)
plotPost(jag_info$ULb2)
plotPost(jag_info$ULb3)

#k
paste("The averare Titanic passenger had a ", round(mean(jag_info$ULb0),2), "probability of surviving. There is some uncertainty about this estimate, but we can say that there is a 95% probability that this value lies between", round(quantile(jag_info$ULb0,prob=0.025),2),"and ", round(quantile(jag_info$ULb0,prob=0.975),2),".")

paste("(For all estimates of survival probability, we report the average probability as well as the 95% credible interval. The 95% CI gives a range of the uncertainty in our estimate. In the example above, our best estimate is ", round(mean(jag_info$ULb0),2), "but we are estimating this based on likelihood and statistical formulas, so there is a 95% chance that the value differs slightly but still lies between", round(quantile(jag_info$ULb0,prob=0.025),2),"and ", round(quantile(jag_info$ULb0,prob=0.025),2),".)")

paste("The probability of survival shifts depending on the class of fare the passenger took. A first class passenger had a ", round(mean(jag_info$ULb1),2), "probability of surviving (the 95% CI around that estimate has a lower bound of",round(quantile(jag_info$ULb1,prob=0.025),2),"and an upper bound of",round(quantile(jag_info$ULb1,prob=0.025),2),".")

paste("A second class passenger had a ",round(mean(jag_info$ULb2),2), "(with a 95% CI of", round(quantile(jag_info$ULb2,prob=0.025),2),",",round(quantile(jag_info$ULb2,prob=0.975),2),") probability of surviving.")

paste("A third class passenger had a ",round(mean(jag_info$ULb3),2), "(with a 95% CI of", round(quantile(jag_info$ULb3,prob=0.025),2),",",round(quantile(jag_info$ULb3,prob=0.975),2),") probability of surviving.")
 
#2
##### model using Fare #####
## Specify the data ##
x1 <- x$Fare

dataList = list(
  x1 = x1,
  y = y 
)

## Specify the model ##
model_string <- "model {
for (i in 1:length(y)){
## Part 1
y[i] ~ dnorm(mu[i] , tau)
logit(mu[i]) <- beta0 + beta1*x1[i]
}

beta0 ~ dnorm(0,0.001)
beta1 ~ dnorm(0,0.001)
tau ~ dunif(1.0E-3 , 1.0E+3)
}"
writeLines(model_string,con="TEMPmodel2.txt")

## Prepare JAGS for sampling from the posterior distribution ##
#Give initial values
initsList = list(
  beta0 = mean(y),
  beta1 = 0,
  tau = 1
)

#Set MCMC parameters
parameters = c("beta0", "beta1", "tau") # parameters to record
adaptSteps = 500  # Number of steps to "tune" the samplers. See '?adapt'
burnInSteps = 1000 # number of steps to allow the Markov chains to run before considering a reliable sample from the posterior distribution
nChains = 3 # number of parallel Markov chains to run
Iter = 5000 # number of iterations of the Markov chain (each of which provides a sample of the posterior distribution)

#Record information into JAGS
jagsModel2 = jags.model("TEMPmodel2.txt", data=dataList, inits=initsList, n.chains=nChains, n.adapt=adaptSteps)

#Burn-in period
update(jagsModel2, n.iter=burnInSteps)

#Run and record MCMC samples
codaSamples2 = coda.samples(jagsModel2, variable.names=parameters, n.iter=Iter)

# traceplot and gelman.plot are from CODA package:
traceplot(codaSamples2[[1]][,"beta1"])
traceplot(codaSamples2[[2]][,"beta1"],col="skyblue",add=T)
traceplot(codaSamples2[[3]][,"beta1"],col="gray",add=T)

# gelman plot
gelman.plot(codaSamples2)

# Interpretation
jag_info2 <- data.frame(rbind(codaSamples2[[1]],codaSamples2[[2]],codaSamples2[[3]]))

paste("The logistic transformation of the value in Survived converts the probability of 1 value to a log odds, or the log (base e) of the odds ratio. The odds ratio is the probability of Success (living) divided by the probability of Failure (dying). If the probability of Success = p, then the probability of dying = 1-p. So our y variable becomes log(p / 1-p), and we construct a linear model to predict this.")

jag_info2$LO_beta1 <- exp(jag_info2$beta1)
paste("One way to interpret the impact of Fare is to report its effects in terms of the odds ratio. So every 1 British Pound increase in Fare increases the odds of surviving by", round(mean(jag_info2$LO_beta1),2),"(95% CI: ",round(quantile(jag_info2$LO_beta1,prob=0.025),2),",",round(quantile(jag_info2$LO_beta1,prob=0.975),2),").")

paste("Another way to interpret the impact of Fare on the probability of survival, is to report the posterior distribution of having a fare of 1 (is it British pounds? I dont know), 50, and 100.")

jag_info2$p_fare1 <- logit_to_prob(jag_info2$beta0+jag_info2$beta1*1)
jag_info2$p_fare50 <- logit_to_prob(jag_info2$beta0+jag_info2$beta1*50)
jag_info2$p_fare100 <- logit_to_prob(jag_info2$beta0+jag_info2$beta1*100)

mean(jag_info2$p_fare1)
quantile(jag_info2$p_fare1,probs=c(.025,.975))

mean(jag_info2$p_fare50)
quantile(jag_info2$p_fare50,probs=c(.025,.975))

mean(jag_info2$p_fare100)
quantile(jag_info2$p_fare100,probs=c(.025,.975))

paste("A passenger who paid a Fare of 1 pound had a ",round(mean(jag_info2$p_fare1),2),"probability of surviving (95% CI: ", round(quantile(jag_info2$p_fare1,probs=.025),2),",",round(quantile(jag_info2$p_fare1,probs=.975),2),").")

paste("A passenger who paid a Fare of 50 pounds had a ",round(mean(jag_info2$p_fare50),2),"probability of surviving (95% CI: ", round(quantile(jag_info2$p_fare50,probs=.025),2),",",round(quantile(jag_info2$p_fare50,probs=.975),2),").")

paste("A passenger who paid a Fare of 100 pounds had a ",round(mean(jag_info2$p_fare100),2),"probability of surviving (95% CI: ", round(quantile(jag_info2$p_fare100,probs=.025),2),",",round(quantile(jag_info2$p_fare100,probs=.975),2),").")

# HDI
plotPost(jag_info2$LO_beta1)
plotPost(jag_info2$p_fare1)
plotPost(jag_info2$p_fare50)
plotPost(jag_info2$p_fare100)

#3 Model comparison with DIC
#a
dic_mod1 <- dic.samples(jagsModel,n.iter=3000)
mean(dic_mod1$deviance)
quantile(dic_mod1$deviance,probs=c(.025,.975))

#b
dic_mod2 <- dic.samples(jagsModel2,n.iter=3000)
mean(dic_mod2$deviance)
quantile(dic_mod2$deviance,probs=c(.025,.975))

#c
c(mean(dic_mod1$deviance),mean(dic_mod2$deviance))
#Model 1 gives a better fit. We haven't done a full model comparison with information criterion the way we did in homeworks with non-Bayesian GLMs, but you can see how we can use DIC to evaluate whether or not to include predictors in a full GLM model. A next step here would be to test Fare and Class in the same model and compare DIC with and without each predictor.