setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Meldrum Research/model")

library(readr)

mdata <- read_csv("subset_tg.csv")

mdata <- as.data.frame(mdata)

# Remove all NA rows
mdata <- mdata[rowMeans(is.na(mdata))!=1,]

# Remove all NA columns
mdata <- mdata[,colMeans(is.na(mdata))==0]

colnames(mdata) <- c("ID", "Age", "T2", "Onset","Gloss", "PVC", "Oil", "Acrylic", "Metal", "Organic", "Inorganic")



model1 <- lm(Onset~ Gloss + Age + T2 + PVC + Oil + Metal + Organic , data = mdata)
summary(model1)

model2 <- step(model1)
summary(model2)
null <- lm(Gloss ~ 1, data = mdata)

AIC(null)
AIC(model1)
AIC(model2)


#mod3 <- lm(Gloss ~ T2 + PVC, data = mdata)
#AIC(mod3)


pchisq(2 * (logLik(model2) - logLik(null)), df = 6, lower.tail=FALSE)
cbind(model2$coefficients, confint(model2)[,1], confint(model2)[,2])

Pred <- cbind.data.frame(mdata$T2, mdata$PVC, mdata$Metal, mdata$Organic)
new <- 

