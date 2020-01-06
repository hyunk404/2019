################################
## 2019/10/29 
## Project v1
## Golub et al., data analysis
#################################



## Set environment
setRepositories()

library(multtest)
library(preprocessCore)

setwd("C:\\Users\\hyunk\\Desktop\\research assistant")


## Loading golub data
## https://science.sciencemag.org/content/286/5439/531

data(golub)
?golub


dim(golub)
golub <- t(golub)

colnames(golub) <- golub.gnames[,3]

Group <- golub.cl

Group[which(Group == 0)] <- "ALL"
Group[which(Group == 1)] <- "AML"
Group <- as.factor(Group)

## Normalization
dim(golub)
normalizedGolub <- normalize.quantiles(golub)
normalizedGolub1 <- normalize.quantiles(t(golub))
boxplot(normalizedGolub)
boxplot(normalizedGolub1)

## EDA
par(mfrow=c(1, 2))

boxplot(t(golub))
boxplot(t(normalizedGolub))


## Developing quantile normalization tool

quantileNormABC <- function(input){
  warning("This function provides quantile normalization for matrix data...")
  ncol <- ncol(input)

  tmpOutput <- matrix(NA, nrow = nrow(input), ncol = ncol(input))
  indexOutput <- matrix(NA, nrow = nrow(input), ncol = ncol(input))
  finalOutput <- matrix(NA, nrow = nrow(input), ncol = ncol(input))
  
  #Rank: Rank, Order: index, Sort: value
  for(i in 1:ncol){
    tmpOutput[,i] <- sort(input[,i])
    indexOutput[,i] <- rank(input[,i])
  }
  
  baseValue <- rowSums(tmpOutput) / ncol(tmpOutput)
  
  # Assign value to each quantile
  for(i in 1:length(baseValue)){
    finalOutput[which(indexOutput == i)] <- baseValue[i]
  }
  
  return(finalOutput)
}


## Generating test dataset
A <- c(1,5,3,5)
B <- c(2,1,6,7)
C <- c(3,2,2,6)
D <- c(4,6,1,8)

testData <- t(data.frame(A,B,C,D))

rowSums(testData)
normalize.quantiles(testData)
quantileNormABC(testData)













