## 11.26 - anova
## one-way
## two-way

setRepositories()

WORK_DIR = "C:\\Users\\정현교\\Desktop\\hyunk\\research assistant\\Project\\Tobacco_anova"

setwd(WORK_DIR)
library(data.table)
library(dplyr)

Expression <- fread("Expression.csv", header=T, sep=",", stringsAsFactors = F)
Meta <- fread("MetaData.csv", header=T, sep=",", stringsAsFactors = F)

Expression <- data.frame(Expression)
Meta <- data.frame(Meta)

dim(Expression)
dim(Meta)

## preprocessing
Expression <- Expression[,-c(186:206)]

probeID <- Expression$ID_REF
geneSymbol <- Expression$IDENTIFIER

Expression <- Expression[,-c(1:2)]

dim(Expression)
dim(Meta)

all.equal(colnames(Expression), Meta$SampleID)

save.image("saveData.RData")

summarise(group_by(Meta, Condition, Group))
length(which(Meta$Group =='non-smoker'))


#### ANOVA
Group_CordBlood <- which(Meta$Condition == "cord blood")
Group_Periheral <- which(Meta$Condition == "maternal peripheral")
Group_placenta <- which(Meta$Condition == "term placenta")

subExpression <- Expression[,c(Group_CordBlood,Group_Periheral,Group_placenta)]
subMeta <- Meta[c(Group_CordBlood,Group_Periheral,Group_placenta),]

Group <- factor(subMeta$Condition, levels = c("cord blood", "maternal peripheral", "term placenta"))


#### ignore Group(smoking), one-way anova
Pvalue_aov <- c()
Pvalue_new <- c()
for(i in 1:nrow(subExpression)){
  exp <- as.numeric(subExpression[i,])
  Pvalue_aov[i] <- summary(aov(exp~Group))[[1]][1,5]
  Pvalue_new[i] <- oneway.test(exp~Group)$p.value
  if(i%%1000 == 0){
    print(paste0(i, "th Gene was performed..."))
  }
}

sum(Pvalue_aov <= 0.05)
sum(Pvalue_new <= 0.05)
sum(p.adjust(Pvalue_aov, method = "BH") <= 0.05)
sum(p.adjust(Pvalue_new, method = "BH") <= 0.05)
#visualization
plot((-log10(Pvalue_aov)),(-log10(Pvalue_new)))
abline(0,1, col = "red")

indexSig <-order(Pvalue_aov)[1:200]
geneSymbol[indexSig]

plot_data <- data.frame(Expression = as.numeric(t(as.matrix(Expression[indexSig,]))),
                        Group = rep(Group, length(indexSig)),
                        Gene =  factor(rep(geneSymbol[indexSig], each = 183)))

ggplot(plot_data, x = "Group", y = "Expression", fill = "Group", facet.by = "Gene",  scales = "free_y")

library(ggplot2)

exp <- as.numeric(Expression[indexSig[1],])
boxplot(exp~Group)


####include Group(smoking), two-way anova
smoke <- Meta$Group

SmokerGroup <- which(smoke =="smoker")
NonSmokerGroup <- which(smoke =="non-smoker")

subExpression_smoke <- Expression[,c(SmokerGroup, NonSmokerGroup)]
subMeta_smoke <- Meta[c(SmokerGroup, NonSmokerGroup),]

Group_smoke <- factor(subMeta_smoke$Group, levels = c("non-smoker", "smoker"))

##############################################################여기가 문제

onehot <- model.matrix(~Group_smoke, subExpression_smoke)
length(onehot)

Pvalue_aov2 <- c()
Pvalue_new2 <- c()

exp <- as.numeric(subExpression_smoke[1,])


for(i in 1:nrow(subExpression_smoke)){
  exp <- as.numeric(subExpression_smoke[i,])
  Pvalue_aov2[i] <- summary(aov(exp ~ Group_smoke * onehot[,2]))[[1]][1,5]
  Pvalue_new2[i] <- oneway.test(exp ~ Group_smoke + onehot[,2])$p.value
  
  if(i%%1000 == 0){
    print(paste0(i, "th Gene was performed..."))
  }
}
  ##############################################################여기까지 문제
Pvalue_aov2 <- as.numeric(Pvalue_aov2)

plot((-log10(Pvalue_aov2)),(-log10(Pvalue_new2)))
abline(0,1, col = "red")


sum(Pvalue_aov2 <= 0.05)
sum(Pvalue_new2 <= 0.05)

sum(p.adjust(Pvalue_aov2, method = "BH") <= 0.05)
sum(p.adjust(Pvalue_new2, method = "BH") <= 0.05)

indexSig <-order(Pvalue_aov2)[1:200]
geneSymbol[indexSig]

plot_data <- data.frame(Expression = as.numeric(t(as.matrix(Expression[indexSig,]))),
                        Group = rep(Group, length(indexSig)),
                        Gene =  factor(rep(geneSymbol[indexSig], each = 183)))
plot_data

library(ggpubr)
ggboxplot(head(plot_data, 4000), x = "Group", y = "Expression", fill = "Group", facet.by = "Gene",  scales = "free_y")

exp <- as.numeric(Expression[indexSig[1],])
boxplot(exp~Group)

#visualization
