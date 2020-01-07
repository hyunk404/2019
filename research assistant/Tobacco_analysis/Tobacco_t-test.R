
setRepositories()

WORK_DIR = "C:\\Users\\hyunk\\Desktop\\research assistant\\Project\\Tobacco"
DATA_DIR = "C:\\Users\\hyunk\\Desktop\\research assistant\\Project\\Tobacco\\Data"
library(data.table)
library(dplyr)

setwd(WORK_DIR)

Expression <- fread("Expression.csv", header=T, sep=",", stringsAsFactors = F)
Meta <- fread("MetaData.csv", header=T, sep=",", stringsAsFactors = F)

Expression <- data.frame(Expression)
Meta <- data.frame(Meta)

dim(Expression)
dim(Meta)

Expression <- Expression[,-c(186:206)]

probeID <- Expression$ID_REF
geneSymbol <- Expression$IDENTIFIER

Expression <- Expression[,-c(1:2)]


dim(Expression)
dim(Meta)

all.equal(colnames(Expression), Meta$SampleID)

save.image("saveData.RData")

##############################Data preprocess######################################

nrow(Meta[Meta$Group == 'non-smoker',])
nrow(Meta[Meta$Group == 'smoker',])
#length(which(Meta$Group =='smoker'))


SmokerGroup <- which(Meta$Condition =="cord blood" & Meta$Group =="smoker")
NonSmokerGroup <- which(Meta$Condition =="cord blood" & Meta$Group =="non-smoker")

subExpression <- Expression[,c(SmokerGroup, NonSmokerGroup)]
subMeta <- Meta[c(SmokerGroup, NonSmokerGroup),]

dim(subExpression)
dim(subMeta)
subMeta

Group <- factor(subMeta$Group, levels = c("non-smoker", "smoker"))

Pvalue <- c()

for(i in 1:nrow(subExpression)){
  exp <- as.numeric(subExpression[i,])
  Pvalue[i] <- t.test(exp~Group)$p.value
  
  if(i%%1000 == 0){
    print(paste0(i, "th Gene was performed..."))
  }
}

sum(Pvalue <= 0.05)

selectGene <- which(Pvalue <= 0.05)
newSet <- Expression[selectGene,]

dim(netSet)



