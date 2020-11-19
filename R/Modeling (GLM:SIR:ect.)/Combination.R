#Alex Fitz


#libraries
library(ggplot2) 
library(MASS)
library(pscl)
library(rjags)
library(coda)

#### Section 1 ####
#read data
df<- read.table("h2n3.txt", header = T) #Read gene bank data
head(df,5)
pop<- read.table("pop.txt", header = T) #Read populaion data
head(pop,5)

#Question 1
pc<- prcomp(df, center = T, scale. = F) #create PC object
summary(pc)

#Question 2
#Proportion of variance explained by Principal components
vars <- apply(pc$x, 2, var)  
props <- vars / sum(vars)
cumsum(props) #Cumulative sum explained by the first 150 components

qplot(1:ncol(pc$x), cumsum(props)) +xlab("PC") + #Plot cumulative sum
  ylab("Cumulative variance") +theme_minimal()

#Question 3
sum <- cumsum(props[1:150])[150] #Variance explained by the first 150 component is 99.86%
paste("The variance explained by the first 150 PC axis is", round(sum,4))

#Question 4
#Extract data frame for the 1st 150 components
pcdata<- data.frame(pc$x[,1:150], pop=pop$x) #create new data for lda

fit<- lda(pop~., data=pcdata) #perform lda
fit

#Question 5
pred<- predict(fit, pcdata) #predict values of each strain by LDA components
pred

#Question 6
head(pred$x,5) #Correct values as in the assignment

#Question 7
ldata<- data.frame(pop=factor(pop$x), pred$x) #Create a data frame for plotting

ggplot(aes(LD1, LD2, col=pop), data=ldata) +geom_point() + theme_minimal() +
  scale_color_manual(values = c("red", "orange", "yellow", "blue", "green", "purple"))

print("We can see that 2006 population clearly deviates from the other 5.")

##Part(i)
qplot(LD2, fill=pop, data=ldata, geom="histogram")+ theme_minimal()+
  scale_fill_manual(values = c("red", "orange", "yellow", "blue", "green", "purple"))
print("We can also see from plotting the histogram for LDA 2 that is clearly separates 2006 species from the remaining")

##Part(ii)
#Extract score for each allele from the PCA
allele_pc_score<- data.frame(pc$rotation[,1:150])
dim(allele_pc_score)
#Shows correct dimensions for allele_pc_score
head(allele_pc_score,5)

##Part(iii)
#Extract scaling with LD2
pc_ld_score<- fit$scaling[,2]
head(pc_ld_score)

##Part(iv)
#Extract value for each allele
allele_ld_score<- apply(allele_pc_score, 1, function(x){ #Multiply each row for each allele by the LD scaling
  sum(x*pc_ld_score)
})
head(allele_ld_score)

#Plot to identify outliers
plot(allele_ld_score, col="royal blue", ylab="Allele scaling with LDA2")

which(allele_ld_score>3 | allele_ld_score< -3) #Alleles X384.g, X384.t, X906.c, X906.t are the 4 alleles


#We can summarize data to count such mutations in 2005 and 2006

df$pop<- pop$x
mutations<- df[pop$x %in% c(2005, 2006),]
mutations<- with(mutations,rbind(tapply(X384.g, pop, mean),
                                 tapply(X384.t, pop, mean),
                                 tapply(X906.c, pop, mean),
                                 tapply(X906.t, pop, mean)))
cbind(c("X384.g", "X384.t","X906.c", "X906.t"), mutations)

#Two main loci were involved (X384 and X906)

#We can see that 2 main loci were involved which are X384 and X906
#For X384, In 2005 the G allele frequency was 0.4% compared to 29.4% in 2006 while
#the T allele frequency was 99.5% compared to 70.6% in 2006


#For X906, In 2005 the C allele frequency was 0.2% compared to 63.7% in 2006
#while In 2005 the T allele frequency was 99.8% compared to 36.3% in 2006


#### Section 2 ####
#Read in data
exp<- read.table("comm.txt", header = T)
head(exp)

#Question 1: Using AICs 
#Part I) Gaussian model
mod1<- glm(Daphnia_magna~fish*macrophytes*adapted, family=gaussian(link="identity"), data=exp)
summary(mod1) 

#Part II) Poisson model
mod2<- glm(Daphnia_magna~fish*macrophytes*adapted, family=poisson(link="log"), data=exp)
summary(mod2) 

#Part III) Negative binomial distribution
mod3<- glm.nb(Daphnia_magna~fish*macrophytes*adapted, data=exp)
summary(mod3) 

#Part IV) Zero inflated poisson model
mod4<- zeroinfl(Daphnia_magna~fish*macrophytes*adapted, data=exp, dist = "poisson")
summary(mod4) 

#Part V) Zero inflated negative binomial model
mod5<- zeroinfl(Daphnia_magna~fish*macrophytes*adapted, data=exp, dist = "negbin")
summary(mod5) 

#Part VI)AIC for each of the 5 models
AIC(mod1, mod2, mod3, mod4, mod5)
print("We can see that AIC is lowest form model 5 (Zero inflated model with negative binomial distribution")

#Here I saved the model for use in the next question 
mod5noint<- zeroinfl(Daphnia_magna ~ fish + macrophytes + adapted + fish*macrophytes + fish*adapted + 
                       macrophytes*adapted, data=exp, dist = "negbin")

#Question 2: Compare models
#Compare using AIC
AIC(mod5, mod5noint) 
print("The AIC for the model with no interation is lower than that for the model with interaction, which means that the model with no interaction is more probable")

#Compare using LRtest
1-pchisq(-2*(mod5noint$loglik- mod5$loglik), #Difference in deviance
         mod5noint$df.residual- mod5$df.residual) #Difference in df
print("P value is 0.176 which is greater than 0.05 which indicates that adding the interaction is not significant. This is also confirmed by the AIC test")

#Question 3: Compare model with interaction (our final choice) to null model
null<- zeroinfl(Daphnia_magna ~ 1, data=exp, dist = "negbin") # Create null model

#AIC
AIC(null, mod5noint) 
print("AIC is lower for the model with predictors than the Null model.")

#LR test
1-pchisq(-2*(null$loglik- mod5noint$loglik), #Difference in deviance
         null$df.residual-  mod5noint$df.residual) #Difference in df

print("The p values is less than 0.05 which indicates that our model (with predictors) is better than the null model.")

#Question 4: Plot predicted values
#Create data that contains predicted value
newdata<- expand.grid(fish=c(0,1), macrophytes=c(0,1), adapted=c(0,1)) #Add possible combinations to a new data frame
newdata$pred<- predict(mod5noint, newdata = newdata, type = "response") #Add predcited values
newdata$category<- factor(paste(newdata$fish, newdata$macrophytes),
                          levels = c("0 0", "1 0", "0 1", "1 1"),
                          labels = c("No fish/No Macrophytes", "Fish/No Macrophytes", "No Fish/Macrophytes",
                                     "Fish/Macrophytes")) #Label variables
newdata$adapted<- factor(newdata$adapted, levels = c(0,1), #Label enviroment variable
                         labels = c("Not adapted", "Adapted"))

ggplot(aes(category, pred, col=adapted, group=adapted), #Aethsetics
       data=newdata)+geom_point(size=2) + #data and point size
  theme_minimal() + #Adjust backgroud
  scale_color_manual(name= "Enviroment", values = c("orange", "royal blue")) + #colors
  geom_line() + #add connecting lines
  ylab("Predicted D.magna abundance") +xlab("Category")+ #Labels
  theme(axis.text.x = element_text(angle = 30, size=8)) #Legend


#### Section 3 ####

#Question 1: Jags Analysis
dat <- read.csv("http://www4.stat.ncsu.edu/~reich/ST590/assignments/Obama2012.csv")
head(dat)

### Specify the data ###
y <- scale(dat$PCTOBAMA)
y <- y[,1]
x <- scale(dat[,colnames(dat) ==  "DiversityIndex" | colnames(dat) == "pct_Uninsured" | colnames(dat) == "pct_unemployed"])

dataList = list(
  x1 = x[,1],
  x2 = x[,2],
  x3 = x[,3],
  y = y 
)

### Specify the model
model_string <- "model {

for (i in 1:length(y)) {
y[i] ~ dnorm(mu[i], sigma)
mu[i] <- beta0 + beta1*x1[i] + beta2*x2[i] + beta3*x3[i]}

# Priors - Beta
beta0 ~ dnorm(0,0.0001)
beta1 ~ dnorm(0,0.0001)
beta2 ~ dnorm(0,0.0001)
beta3 ~ dnorm(0,0.0001)

tau ~ dunif(1.0E-3, 1.0E+3)
sigma <- 1/sqrt(tau)
}"

#write model to external text file for JAGS
writeLines(model_string,con="model.txt")

### Prepare JAGS for sampling from the posterior distribution ##
#Give initial values
lmInfo = lm(y ~ x1 + x2 + x3, data = dataList) 
b0Init = lmInfo$coef[1]
b1Init = lmInfo$coef[2]
b2Init = lmInfo$coef[3]
b3Init = lmInfo$coef[4]
tauInit = length(y) / sum(lmInfo$res^2)

inits = list(
  beta0 = b0Init ,
  beta1 = b1Init ,
  beta2 = b2Init,
  beta3 = b3Init,
  tau = tauInit
)

#Set MCMC parameters
parameters = c("beta0", "beta1", "beta2", "beta3", "tau") # parameters to record
adaptSteps = 500  # Number of steps to "tune" the samplers. See '?adapt'
burnInSteps = 1000 # number of steps to allow the Markov chains to run before considering a reliable sample from the posterior distribution
nChains = 3 # number of parallel Markov chains to run
Iter = 5000 # number of iterations of the Markov chain (each of which provides a sample of the posterior distribution)

### Read in JAGS Model + Sample from posterior distribution
#Record information into JAGS
jagsModel = jags.model("model.txt", data=dataList, inits=inits, n.chains=nChains, n.adapt=adaptSteps)

#Burn-in period
update(jagsModel, n.iter=burnInSteps)

#Run and record MCMC samples
codaSamples = coda.samples(jagsModel, variable.names=parameters, n.iter=Iter)


#Question 2: Plots
# traceplot and gelman.plot are from CODA package:
traceplot(codaSamples)
traceplot(codaSamples[[1]][,"beta1"])
traceplot(codaSamples[[2]][,"beta1"],col="skyblue",add=T)
traceplot(codaSamples[[3]][,"beta1"],col="gray",add=T)

print("When looking at the chains produced by the traceplot() we see convergence.")

# gelman plot
gelman.plot(codaSamples)

print("When looking at the Gelman Plots we see shrink factors of slightly above 1 for all of the beta coefficients. This indicates strong convergence between the chains. A shrink factor of 1 means the variance between chains and within chains are equal to each other.")

#Density
HDIofMCMC <- dget("HDIofMCMC.R")
DbdaDensPlot <- dget("DbdaDensPlot.R")
DBDAplColors = c("black","skyblue","gray")
DbdaDensPlot(codaSamples,"beta0",plColors=DBDAplColors)
DbdaDensPlot(codaSamples,"beta1",plColors=DBDAplColors)
DbdaDensPlot(codaSamples,"beta2",plColors=DBDAplColors)
DbdaDensPlot(codaSamples,"beta3",plColors=DBDAplColors)

#Getting probabilities
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

jag_info <- data.frame(rbind(codaSamples[[1]],codaSamples[[2]],codaSamples[[3]]))
jag_info$EX_beta0 <- logit2prob(jag_info$beta0)
jag_info$EX_beta1 <- logit2prob(jag_info$beta0 + jag_info$beta1)
jag_info$EX_beta2 <- logit2prob(jag_info$beta0 + jag_info$beta2)
jag_info$EX_beta3 <- logit2prob(jag_info$beta0 + jag_info$beta3)

# HDI
plotPost <- dget("plotPost.R")
plotPost(jag_info$EX_beta0)
plotPost(jag_info$EX_beta1)
plotPost(jag_info$EX_beta2)
plotPost(jag_info$EX_beta3)

#Question 3: Verbal Summary
#Interpretation
paste("A person that is of a higher diversityindex", round(mean(jag_info$EX_beta1),2),"percent chance of voting for Obama", 
      round(mean(jag_info$EX_beta1),2),"(95% CI: ",round(quantile(jag_info$EX_beta1,prob=0.025),2),",",round(quantile(jag_info$EX_beta1,prob=0.975),2),
      "). Since this is greater than 50% it means that a higher diversity rate leads to more votes for Obama.")

paste("A person that is uninsured has a", round(mean(jag_info$EX_beta2),2),"percent chance of voting for Obama", 
      "(95% CI: ",round(quantile(jag_info$EX_beta2,prob=0.025),2),",",round(quantile(jag_info$EX_beta2,prob=0.975),2),
      "). Since this is less than 50% it means that uninsured are less likely to vote for Obama.")

paste("A person that is unemployed has a", round(mean(jag_info$EX_beta3),2),"percent chance of voting for Obama","(95% CI: ",round(quantile(jag_info$EX_beta3,prob=0.025),2),",",round(quantile(jag_info$EX_beta3,prob=0.975),2),
      "). Since this is less than 50% it means that uninsured are less likely to vote for Obama.")

paste("It can be inferred that beta 1 is positive", round(mean(jag_info$beta1),2), 
      "thus higher diversity Index leads to higher proportion of votes Obama. Similarly, a higher pct_Uninsured leads to a lower proportion of voters for Obama",round(mean(jag_info$beta2),2), 
      "and it is the same for pct_unemployed",round(mean(jag_info$beta3),2),").")

#Question 4: Plot of Posterior Distribution
# Independent variable values for Durham
DC <- x[32,]
#Since the original data is scaled we use the sd and mean to unscale it. 
y_pred <- (samps.coda[[1]][,2:4]%*%DC + samps.coda[[1]][,5])*sd(dat$PCTOBAMA) + mean(dat$PCTOBAMA)
plot(density(y_pred),main = "Density Plot of the Predicted Percent of Obama voters in Durham")
plotPost(y_pred)
print("It is clearly visible that the mean is quite high around 0.6 which further confirms the fact that Democrats are more favoured in Durham.")

jag_info2 <- data.frame(rbind(codaSamples[[1]],codaSamples[[2]],codaSamples[[3]]))
jag_info2[,2] <- jag_info2[,2]*1.5255
jag_info2[,3] <- jag_info2[,3]*-0.0683
jag_info2[,4] <- jag_info2[,4]*-1.4485
jag_info2$EX_beta0 <- logit2prob(jag_info2$beta0)
jag_info2$EX_beta1 <- logit2prob(jag_info2$beta0 + jag_info2$beta1)
jag_info2$EX_beta2 <- logit2prob(jag_info2$beta0 + jag_info2$beta2)
jag_info2$EX_beta3 <- logit2prob(jag_info2$beta0 + jag_info2$beta3)

paste("A person that is of a higher diversityindex", round(mean(jag_info2$EX_beta1),2),"percent chance of voting for Obama", 
      round(mean(jag_info2$EX_beta1),2),"(95% CI: ",round(quantile(jag_info2$EX_beta1,prob=0.025),2),",",round(quantile(jag_info2$EX_beta1,prob=0.975),2),
      "). Since this is greater than 50% it means that a higher diversity rate leads to more votes for Obama.")

paste("A person that is uninsured has a", round(mean(jag_info2$EX_beta2),2),"percent chance of voting for Obama", 
      "(95% CI: ",round(quantile(jag_info2$EX_beta2,prob=0.025),2),",",round(quantile(jag_info2$EX_beta2,prob=0.975),2),
      ").Since this is close to 50% it means that a higher uninsured rate wouldn't necessarily lead to more votes for Obama.")

paste("A person that is unemployed has a", round(mean(jag_info2$EX_beta3),2),"percent chance of voting for Obama","(95% CI: ",round(quantile(jag_info2$EX_beta3,prob=0.025),2),",",round(quantile(jag_info2$EX_beta3,prob=0.975),2),
      "). Since this is greater than 50% it means that a higher uneemployment rate leads to more votes for Obama.")

#Plot of Posterior
plotPost(jag_info2$EX_beta0)
plotPost(jag_info2$EX_beta1)
plotPost(jag_info2$EX_beta2)
plotPost(jag_info2$EX_beta3)