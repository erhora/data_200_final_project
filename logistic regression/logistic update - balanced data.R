library(foreign)
library(tidyverse)
library(caret)

#2015 - 2016

best_forward_revised <- readRDS('best_forward_revised.rds') #9544 * 11
BPX_I2015_2016 <- read.xport("BPX_I2015_2016.XPT")
BPX_I2015_2016 <- BPX_I2015_2016[c('SEQN','BPXSY1','BPXDI1')]
colnames(BPX_I2015_2016) <- c('seqn','bpxsy1','bpxdi1')
complete <- left_join(best_forward_revised, BPX_I2015_2016, by = 'seqn')

sum(is.na(complete$bpxsy1)) #2399
sum(is.na(complete$bpxdi1)) #2399
complete <- complete[!is.na(complete$bpxsy1), ] #7145 * 13
sum(complete$bpxsy1 >= 140) #1024
sum(complete$bpxdi1 >= 90) #282
nrow(complete[complete$bpxsy1 >= 140 | complete$bpxdi1 >= 90, ]) #1108

hypertension <- c(ifelse(complete$bpxsy1 >= 140 | complete$bpxdi1 >= 90, 1, 0))
complete <- mutate(complete, hypertension = hypertension)
complete <- complete %>% 
  select(-seqn) %>% 
  select(-bpxsy1) %>% 
  select(-bpxdi1)
head(complete)

#logistic regression - variables were chosen with forward selection, 2015 - 2016 dataset

set.seed(1)
inTrain <- sample(nrow(complete), nrow(complete)*.7)
traindata <- complete[inTrain, ]
testdata <- complete[-inTrain, ]

m1 <- glm(hypertension ~ ., family = 'binomial', data = traindata)
summary(m1)

m2 <- glm(hypertension ~ ridageyr+bmxbmi+riagendr+paq715+total_caffeine+count_meds+section_I, 
          family = 'binomial', data = traindata)
summary(m2)
probabilities <- predict(m2, type="response", newdata = testdata)
#how to set the threshold
y_pred <- ifelse(probabilities > 0.5, 1, 0) 
confusionMatrix(as.factor(y_pred), as.factor(testdata$hypertension)) #0.8414
#    0    1
#0 1608  265
#   48   52

#imbalanced data - problem: predicting most as 0


#use a subset of non-hypertension data

hypertensiondata <- complete[complete$hypertension == 1,] #1108
nonhypertensiondata <- complete[complete$hypertension == 0,] #6037
sub <- sample(nrow(nonhypertensiondata), nrow(nonhypertensiondata)*.2) #1207
subset <- rbind(hypertensiondata, nonhypertensiondata[sub, ]) #2315 * 11

set.seed(1)
inTrain <- sample(nrow(subset), nrow(subset)*.7)
traindata <- subset[inTrain, ]
testdata <- subset[-inTrain, ]

m3 <- glm(hypertension ~ ridageyr+bmxbmi+riagendr+count_meds+section_I, 
          family = 'binomial', data = traindata)
summary(m3)
probabilities <- predict(m3, type="response", newdata = testdata)
y_pred <- ifelse(probabilities > 0.5, 1, 0) 
confusionMatrix(as.factor(y_pred), as.factor(testdata$hypertension)) #0.7591

#   0   1
#0 241  76
#1  70 219


#Cluster the non-hypertension data first 

hypertensiondata <- hypertensiondata %>% 
  drop_na() #991 * 11
nonhypertensiondata <- nonhypertensiondata %>% 
  drop_na() #5370 * 11
  
cl <- kmeans(nonhypertensiondata, 991, nstart = 1000)
centroids <- cl$centers

clustersubset <- rbind(hypertensiondata, centroids)

set.seed(1)
inTrain <- sample(nrow(clustersubset), nrow(clustersubset)*.7)
traindata <- clustersubset[inTrain, ]
testdata <- clustersubset[-inTrain, ]

m4 <- glm(hypertension ~ ridageyr+bmxbmi+riagendr+total_caffeine+count_meds+section_I, 
          family = 'binomial', data = traindata)
summary(m4)
probabilities <- predict(m3, type="response", newdata = testdata)
y_pred <- ifelse(probabilities > 0.5, 1, 0) 
confusionMatrix(as.factor(y_pred), as.factor(testdata$hypertension)) #0.7361

#   0   1
#0 221  71
#1  86 217

