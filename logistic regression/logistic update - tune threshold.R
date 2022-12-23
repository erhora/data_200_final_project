library(foreign)
library(tidyverse)
library(caret)
library(glmnet)
#library(pROC)

#2015 - 2016

best_forward_revised <- readRDS('best_forward_revised.rds') #9544 * 11
BPX_I2015_2016 <- read.xport("BPX_I2015_2016.XPT")
BPX_I2015_2016 <- BPX_I2015_2016[c('SEQN','BPXSY1','BPXDI1')]
colnames(BPX_I2015_2016) <- c('seqn','bpxsy1','bpxdi1')
complete <- left_join(best_forward_revised, BPX_I2015_2016, by = 'seqn')

#sum(is.na(complete$bpxsy1)) #2399
#sum(is.na(complete$bpxdi1)) #2399
#complete <- complete[!is.na(complete$bpxsy1), ] #7145 * 13
#sum(complete$bpxsy1 >= 140) #1024
#sum(complete$bpxdi1 >= 90) #282
#nrow(complete[complete$bpxsy1 >= 140 | complete$bpxdi1 >= 90, ]) #1108

complete <- drop_na(complete) #6361

hypertension <- c(ifelse(complete$bpxsy1 >= 140 | complete$bpxdi1 >= 90, 1, 0))
complete <- mutate(complete, hypertension = hypertension)
complete <- complete %>% 
  select(-seqn) %>% 
  select(-bpxsy1) %>% 
  select(-bpxdi1)
head(complete)

#logistic regression - variables were chosen with forward selection, 2015 - 2016 dataset

#set.seed(1)
inTrain <- sample(nrow(complete), nrow(complete)*.7)
traindata <- complete[inTrain, ]
testdata <- complete[-inTrain, ]

m1 <- glm(hypertension ~ ., family = 'binomial', data = traindata)
summary(m1)

m2 <- glm(hypertension ~ ridageyr+bmxbmi+riagendr+paq715+count_meds+section_I, 
          family = 'binomial', data = traindata)
summary(m2)
probabilities <- predict(m2, type="response", newdata = testdata)
y_pred <- ifelse(probabilities > 0.2, 1, 0) 

TP <- sum(y_pred == 1 & testdata$hypertension == 1)
TN <- sum(y_pred == 0 & testdata$hypertension == 0)
FP <- sum(y_pred == 1 & testdata$hypertension == 0)
FN <- sum(y_pred == 0 & testdata$hypertension == 1)

accuracy[idx] <- (TP+TN)/(TP+TN+FP+FN)
sensitivity[idx] <- TP/(TP+FN)
specificity[idx] <- TN/(TN+FP)

#imbalanced data - problem: predicting most as 0
#pay attention to false negatives

#use a subset of non-hypertension data

#set.seed(1)
hypertensiondata <- complete[complete$hypertension == 1,] #1108
nonhypertensiondata <- complete[complete$hypertension == 0,] #6037
sub <- sample(nrow(nonhypertensiondata), nrow(nonhypertensiondata)*.2) #1207
subset <- rbind(hypertensiondata, nonhypertensiondata[sub, ]) #2315 * 11

inTrain <- sample(nrow(subset), nrow(subset)*.7)
traindata <- subset[inTrain, ]
testdata <- subset[-inTrain, ]

m3 <- glm(hypertension ~ ridageyr+bmxbmi+riagendr+paq715+count_meds, 
          family = 'binomial', data = traindata)
summary(m3)
probabilities <- predict(m3, type="response", newdata = testdata)
y_pred <- ifelse(probabilities > 0.5, 1, 0) 

TP <- sum(y_pred == 1 & testdata$hypertension == 1)
TN <- sum(y_pred == 0 & testdata$hypertension == 0)
FP <- sum(y_pred == 1 & testdata$hypertension == 0)
FN <- sum(y_pred == 0 & testdata$hypertension == 1)

accuracy[idx] <- (TP+TN)/(TP+TN+FP+FN)
sensitivity[idx] <- TP/(TP+FN)
specificity[idx] <- TN/(TN+FP)

#Cluster the non-hypertension data first 

#hypertensiondata <- hypertensiondata %>% 
#  drop_na() #991 * 11
#nonhypertensiondata <- nonhypertensiondata %>% 
#  drop_na() #5370 * 11

#set.seed(1)
cl <- kmeans(nonhypertensiondata, 991, nstart = 1000)
centroids <- cl$centers

clustersubset <- rbind(hypertensiondata, centroids)

#set.seed(1)
inTrain <- sample(nrow(clustersubset), nrow(clustersubset)*.7)
traindata <- clustersubset[inTrain, ]
testdata <- clustersubset[-inTrain, ]

m4 <- glm(hypertension ~ ridageyr+bmxbmi+riagendr+total_caffeine+count_meds+section_I+diq010, 
          family = 'binomial', data = traindata)
summary(m4)
probabilities <- predict(m4, type="response", newdata = testdata)
y_pred <- ifelse(probabilities > 0.5, 1, 0) 

TP <- sum(y_pred == 1 & testdata$hypertension == 1)
TN <- sum(y_pred == 0 & testdata$hypertension == 0)
FP <- sum(y_pred == 1 & testdata$hypertension == 0)
FN <- sum(y_pred == 0 & testdata$hypertension == 1)

accuracy[idx] <- (TP+TN)/(TP+TN+FP+FN)
sensitivity[idx] <- TP/(TP+FN)
specificity[idx] <- TN/(TN+FP)

#add penalties

x_train <- as.matrix(traindata[c('ridageyr','bmxbmi','riagendr','paq715','total_caffeine','count_meds','section_I')])
y_train <- unlist(traindata['hypertension'])
x_test <- as.matrix(testdata[c('ridageyr','bmxbmi','riagendr','paq715','total_caffeine','count_meds','section_I')])
y_test <- unlist(testdata['hypertension'])

#lasso - larger penalties result in coefficient values that are closer to zero

cv.lasso <- cv.glmnet(x_train, y_train, alpha = 1, nfolds = 5)
fit_lasso <- glmnet(x_train, y_train, alpha = 1, lambda = cv.lasso$lambda.min)
#coef(fit_lasso)
probabilities <- predict(fit_lasso, newx = x_test, s = cv.lasso$lambda.min)
pred_lasso <- ifelse(probabilities > 0.5, 1, 0)
confusionMatrix(as.factor(pred_lasso), as.factor(y_test)) 

#ridge

cv.ridge <- cv.glmnet(x_train, y_train, alpha = 0, nfolds = 5)
fit_ridge <- glmnet(x_train, y_train, alpha = 0, lambda = cv.ridge$lambda.min)
probabilities <- predict(fit_ridge, newx = x_test, s = cv.ridge$lambda.min)
pred_ridge <- ifelse(probabilities > 0.5, 1, 0)
confusionMatrix(as.factor(pred_ridge), as.factor(y_test)) 

#k-fold cross validation

#complete
set.seed(1)
inTrain <- sample(nrow(complete), nrow(complete)*.7)
traindata <- complete[inTrain, ]
testdata <- complete[-inTrain, ]

# define training control
train_control <- trainControl(method = "cv", number = 8)

# train the model on training set
model <- train(hypertension ~ .,
               data = traindata,
               trControl = train_control,
               method = "glm",
               family=binomial())

# print cv scores
summary(model)

accuracy <- c()
sensitivity <- c()
specificity <- c()
idx <- 1
threshold <- seq(0.01, 0.50, by = 0.01)

for (i in threshold){
  probabilities <- predict(model, newdata = testdata)
  y_pred <- ifelse(probabilities > i, 1, 0) 
  
  TP <- sum(y_pred == 1 & testdata$hypertension == 1)
  TN <- sum(y_pred == 0 & testdata$hypertension == 0)
  FP <- sum(y_pred == 1 & testdata$hypertension == 0)
  FN <- sum(y_pred == 0 & testdata$hypertension == 1)
  
  accuracy[idx] <- (TP+TN)/(TP+TN+FP+FN)
  sensitivity[idx] <- TP/(TP+FN)
  specificity[idx] <- TN/(TN+FP)
  idx <- idx + 1
}

#plot

library(ggplot2)
metrics <- c(rep('accuracy',50), rep('sensitivity',50), rep('specificity',50))
results <- c(accuracy, sensitivity, specificity)
df <- data.frame(threshold = threshold, results = results, metrics = metrics)
ggplot(data = df, aes(x = threshold, y = results, group = metrics)) +
  geom_line(aes(color = metrics))
