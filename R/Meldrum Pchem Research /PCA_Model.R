################################################################################
##################################### INIT #####################################
################################################################################
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Meldrum Research/PCA")
#install.packages("vegan")
#install.packages("readr")
#install.packages("corrplot")
library(vegan)
library(readr)
library(corrplot)

################################################################################
################################## LOAD DATA ###################################
################################################################################

theData <- read_csv("PCA data-2.csv")
theData <- as.data.frame(theData)

# Remove all NA rows
theData <- theData[rowMeans(is.na(theData))!=1,]

# Remove all NA columns
theData <- theData[,colMeans(is.na(theData))==0]

# Remove constant columns
theData <- theData[,!apply(theData, 2, function(x) length(unique(x)))==1]

# colors for variables
col <- c("#4D4D4D", "#5DA5DA", "#FAA43A", "#60BD68", "#F15854", "#B276B2",
         "#DECF3F", "#F17CB0", "#B2912F", "#4afff0", "#34bdcc", "#4f61a1",
         "#461e78", "#440a4f", "#c3fbc4", "#85f9d6", "#79c7ad", "#a6cc7a",
         "#dfff7b", "#8d7b88", "#4e414f", "#baadb5", "#2d2538", "#837a80",
         "#fff68f", "#800080", "#f8b1cc", "#c29bff", "#8d0808"
)
col <- col[1:(ncol(theData)-1)]

################################################################################
################################# CORRELATIONS #################################
################################################################################

# Correlations between variables
corrplot(cor(theData[,-1]), method="ellipse")

################################################################################
##################################### PCA ######################################
################################################################################

# run standard pca with variable scaling
pca <- prcomp(theData[,-1], scale=TRUE)

df <- theData[3:5]
df_pca <- prcomp(df, scale = TRUE)
biplot(df_pca)

####################### PROPORTION OF VARIANCE EXPLAINED #######################

# get proportion of variance explaiend by each pc
# (variance = standard deviation squared)
propVar <- pca$sdev^2/sum(pca$sdev^2)

barplot(propVar, col="cornflowerblue", main="Barplot of Variation per PC Axis", las=1, ylim=c(0, 1),
        names=paste0("PC", 1:length(propVar)),
        ylab="proportion of variance explained"
)

# First PC explains more than 40% of variance in the data
# First and second PCs together would explain around 65% of the variance
# First 3 PCs ~ 85% of variance
# According to screeplot ("elbow" method) the whole dataset can be summarized
# using 3 main principal components. We would loose around 15% of information
# in the data by doing that.

#################################### BIPLOT ####################################

# biplot
biplot(pca, xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))

# Shows which variables are correlated with 1st and 2nd principal components.
# PC1 Is mainly influenced by:
#     1) T2
#     2) Gloss
# PC2 is mainly influenced by:
#     1) Age
#     2) Metal ions
# Organic pigment, Inorganic pigment and Pigment Volume Concentration ratio
# seems to affect both PCs.

################################### LOADINGS ###################################

# We can see the same by plotting eigen vectors themselves for the first 3
# principal components
dotchart(pca$rotation[,1:3], col=col, pch=19)

# The same pattern is visible in eigenvector profiles:

#################################### SCORES ####################################

# scores of PCA represent the values of original data points on the new PCA
# axis.
# Original data had 4 groups termed "white", "black", "blue" and "red"
# We can color the points according to those groups
grCol <- tolower(sapply(strsplit(theData[,1], "_"), "[[", 1))

p <- par(bg="grey")
pairs(pca$x[,1:3], col=grCol, pch=19)
par(p)

# Seems like PC2 is most useful for distinguishing the groups. But the best
# separation comes when combining PC1 and PC3 together.

################################################################################
##################################### RDA ######################################
################################################################################

# prepare Response And Explanatory variables
datRes <- theData[,3,drop=FALSE]
colnames(datRes) <- "T2"
datExp <- theData[,-c(1,3)]

# rda function from the vegan package
RDA <- rda(datExp ~ T2, datRes, scale=TRUE)

# triplot of RDA
plot(RDA)

################################################################################
################################## REGRESSION ##################################
################################################################################

# remove the group feature and store in another variable called "dat"
dat <- theData[,-1]
colnames(dat)[2] <- "T2" # simplify the column name for T2

# scale the data
dat <- as.data.frame(scale(dat))

# run linear regression on T2 using all other variables (dot is a shortcut for that)
fit <- lm(T2 ~ ., data=dat)

# check if the explanation of T2 is significant
summary(fit) # look at the p-value on the very bottom

# which variables are the most important in explaining T2?
coef <- coefficients(fit)
round(coef, 2) # answer - pigment volume concentration, gloss and organic pigment

# the same can be seen in the summary, looking at the last column of the table -
# p-values (Pr(>|r|))
summary(fit)

# what is the proportion of variance in T2 explained by other variables?
summary(fit) # look at Multiple R-squared value
# in this case its 0.8629

p <- par(bg="white", mar=c(20,5,2,2))
barplot(abs(coef[-1]), las=3)
par(p)

# Real T2 values VS values reconstructed using other variables
plot(dat$T2, predict(fit, dat), xlab="real", ylab="predicted",
     xlim=c(-2,2), ylim=c(-2,2)
)

# Residuals of real and reconstructed values for each observation
plot(dat$T2-predict(fit, dat), type="h", ylab="real - predicted",
     main="errors"
)
