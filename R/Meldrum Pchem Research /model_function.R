################################################################################
##################################### INIT #####################################
################################################################################
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Meldrum Research/model")
#install.packages("mice")
#install.packages("randomForest")
library(mice)
library(randomForest)

################################################################################
################################## Function ####################################
################################################################################

prediction <- function(myData){
  # impute data using all the observations (but not using the Tg variable)
  theData <- complete(mice(myData[,colnames(myData)!="Onset.Tg"]))
  theData$Onset.Tg <- myData$Onset.Tg
  
  predictionsRF <- numeric(nrow(myData))
  predictionsLM <- numeric(nrow(myData))
  for(i in 1:nrow(myData)) {
    train <- theData[-i,]
    test  <- theData[i,]
    
    rf <- randomForest(Onset.Tg ~ ., data=train, ntree=5000)
    predictionsRF[i] <- predict(rf, test)
    
    lm <- lm(Onset.Tg ~., data=train)
    predictionsLM[i] <- predict(lm, test)
  }
  
  ############################# VARIABLE IMPORTANCE ##############################
  png(filename="rf.png")
  barplot(t(rf$importance), horiz=TRUE, las=1, main = "Importance of Variables for Random Forest Model")
  dev.off()
  
  ################################# PREDICTIONS ##################################
  
  png(filename = "rf-lr.png")
  par(mar=c(5,5,2,2), mfrow=c(1,2))
  plot(theData$Onset.Tg, predictionsRF, pch=19, xlab="real", ylab="predicted", main="Random forest")
  plot(theData$Onset.Tg, predictionsLM, pch=19, xlab="real", ylab="predicted", main="Linear regression")
  dev.off()
  
  #################################### ERRORS ####################################
  
  errRF <- theData$Onset.Tg-predictionsRF
  errLM <- theData$Onset.Tg-predictionsLM
  
  png(filename = "rf-lr-avg.png")
  par(mfrow=c(1,3))
  plot(errRF, type="h", ylim=c(-30,30), main="using random forest")
  plot(errLM, type="h", ylim=c(-30,30), main="using linear regression")
  plot(theData$Onset.Tg-mean(theData$Onset.Tg), type="h", ylim=c(-30,30), main="using average")
  dev.off()
  
  ######################### ERRORS 2 - Density ###################################
  
  png(filename = "Density.png")
  plot(density(errRF), xlab="error", main="error distributions", col="orange", lwd=3)
  points(density(errLM), type="l", col="cornflowerblue", lwd=3)
  legend("topleft", legend=paste(c("Random Forest", "Linear Regression"), round(c(mean(abs(errRF)), mean(abs(errLM))), 2)),
         col=c("orange", "cornflowerblue"), lwd=3
  )
  dev.off()
}

################################################################################
################################## LOAD DATA ###################################
################################################################################

myData <- read.csv("pca_full.csv")

#Fix dataframe
# shorten column names
colnames(myData) <- sapply(strsplit(colnames(myData), "\\.\\."), "[[", 1)
# remove missing dependant variables
myData <- myData[!is.na(myData$Onset.Tg),]
# turn binnary variables to factors
myData$Oil <- factor(myData$Oil)
myData$Acrylic <- factor(myData$Acrylic)
myData$Metal.ions <- factor(myData$Metal.ions)
myData$Water.Mixable <- factor(myData$Water.Mixable)
myData$Organic.pigment <- factor(myData$Organic.pigment)
myData$Inorganic.pigment <- factor(myData$Inorganic.pigment)


################################################################################
################################## CODE ########################################
################################################################################


run <- prediction(myData)
