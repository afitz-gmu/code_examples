################################################################################
##################################### INIT #####################################
################################################################################
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Meldrum Research/model")
#install.packages("mice")
#install.packages("randomForest")
library(mice)
library(randomForest)
library(ggplot2)


################################################################################
################################## LOAD DATA ###################################
################################################################################
set.seed(1234)

myData <- read.csv("pca_oil.csv")

# remove sample ID
myData <- myData[,-1]

# shorten column names
colnames(myData) <- sapply(strsplit(colnames(myData), "\\.\\."), "[[", 1)

# remove missing dependant variables
myData <- myData[!is.na(myData$Onset.Tg),]

# turn binnary variables to factors
#myData$Oil <- factor(myData$Oil)
#myData$Acrylic <- factor(myData$Acrylic)
#myData$Metal.ions <- factor(myData$Metal.ions)
#myData$Water.Mixable <- factor(myData$Water.Mixable)
#myData$Organic.pigment <- factor(myData$Organic.pigment)
#myData$Inorganic.pigment <- factor(myData$Inorganic.pigment)

################################################################################
################################### CLASSIFY ###################################
################################################################################

# impute data using all the observations (but not using the Tg variable)
theData <- complete(mice(myData[,colnames(myData)!="Onset.Tg"]))
theData$Onset.Tg <- myData$Onset.Tg

############################### CROSS VALIDATION ###############################

predictionsRF <- list()
predictionsLM <- list()
for(v in 0:(ncol(theData)-1)) {
  tmpData <- theData[,(1:ncol(theData))!=v]
  predictionsRF[[v+1]] <- numeric(nrow(tmpData))
  predictionsLM[[v+1]] <- numeric(nrow(tmpData))
  for(i in 1:nrow(tmpData)) {
    train <- tmpData[-i,]
    test  <- tmpData[i,]
    
    rf <- randomForest(Onset.Tg ~ ., data=train, ntree=100)
    predictionsRF[[v+1]][i] <- predict(rf, test)
    
    lm <- lm(Onset.Tg ~., data=train)
    predictionsLM[[v+1]][i] <- predict(lm, test)
  }
}

names(predictionsRF) <- c("NONE", colnames(theData)[-ncol(theData)])
names(predictionsLM) <- c("NONE", colnames(theData)[-ncol(theData)])

############################# VARIABLE IMPORTANCE ##############################

errsRF <- sapply(predictionsRF, function(pr, rl) mean(abs(pr-rl)), theData$Onset.Tg)
errsLM <- sapply(predictionsLM, function(pr, rl) mean(abs(pr-rl)), theData$Onset.Tg)

#par(mar=c(5,15,1,2))
barplot(rbind(errsRF, errsLM), beside=TRUE, col=c("orange", "cornflowerblue"),
        horiz=TRUE, las=2, main="error after excluding a variable", xlab="average error"
)

################################# PREDICTIONS ##################################

par(mar=c(5,4,2,2), mfrow=c(1,2))
plot(theData$Onset.Tg, predictionsRF[[1]], pch=19, xlab="real", ylab="predicted", main="Random forest")
plot(theData$Onset.Tg, predictionsLM[[1]], pch=19, xlab="real", ylab="predicted", main="Linear regression")

#################################### ERRORS ####################################

errRF <- theData$Onset.Tg-predictionsRF[[1]]
errLM <- theData$Onset.Tg-predictionsLM[[1]]

par(mfrow=c(1,3))
plot(errRF, type="h", ylim=c(-30,30), main="using random forest")
plot(errLM, type="h", ylim=c(-30,30), main="using linear regression")
plot(theData$Onset.Tg-mean(theData$Onset.Tg), type="h", ylim=c(-30,30), main="using average")

################################### ERRORS 2 ###################################

par(mfrow=c(1,1))
plot(density(errRF), xlab="error", main="error distributions", col="orange", lwd=3)
points(density(errLM), type="l", col="cornflowerblue", lwd=3)
legend("topleft", legend=paste(c("Random Forest", "Linear Regression"), round(c(mean(abs(errRF)), mean(abs(errLM))), 2)),
       col=c("orange", "cornflowerblue"), lwd=3
)

################################################################################
############################### EXAMPLE SCENARIO ###############################
################################################################################

# indeces for data that we will pretend we don't have Tg measurements for.
inds <- sample(nrow(myData), 10)

allData <- myData[-inds,] # the data you have for making predictions
newData <- myData[inds,]  # the new data you want to predict on

# remove the Onset.TG column from the new data - pretend it's not there
realTG <- newData$Onset.Tg # keep it just for testing our predictions
newData$Onset.Tg <- NULL

# impute both datasets separately without using the Onset.Tg
allData[,-3] <- complete(mice(allData[,-3]))
newData <- complete(mice(newData))

# classify using random forest
rf <- randomForest(Onset.Tg ~ ., data=allData)

# predict on new data
pred <- predict(rf, newData)

# compare to real (we can do it now, but this step will not be possible on a
# real new dataset since we will not know real TG values)
print(cbind(realTG, pred))


################################################################################
############################### GRAPH  #########################################
################################################################################

new <- lm(Onset.Tg~.,  data = theData)
new_p <- predict(new, theData)
print(cbind(theData$Onset.Tg, new_p))
error <- theData$Onset.Tg-new_p
mean(abs(error))

data = cbind.data.frame(predictionsLM[[1]], theData$Onset.Tg)
error <- .8
colnames(data) <- c("Predicted", "Real")

new_data <- data.frame(xcord = seq(-100, 100,1), ycord = seq(-100, 100, 1))

plot <- ggplot() + geom_point(data = data, aes(x=Real, y=Predicted), size=2) 
plot <- plot + geom_ribbon(data = new_data, aes(x=xcord, ymin = ycord-error, ymax = ycord+error), fill = "lightgrey", alpha = .4)
plot <- plot + geom_abline(slope = 1)
plot <- plot + coord_cartesian(xlim= c(-43, -25), ylim= c(-43, -25))
plot <- plot + annotate("text", x = -38, y = -25, label = "Predicted(Tg) = -17.15017 + T2*-41.4746 + Gloss*-0.01445")
plot

setEPS()
postscript("PredictedTg_Graph.eps")
plot
dev.off()

ggsave("PredictedTg_Graph.pdf", plot = plot, device = "pdf")

-17.15017031 + (-41.47469964 * T2) + (-0.01445885 * Gloss)


