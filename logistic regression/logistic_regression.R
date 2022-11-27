#install.packages("foreign")
library(foreign)
library(tidyverse)
library(caret)

#logistic regression - variables were chosen with forward selection, 2015 - 2016 dataset

#RIDAGEYR
#WTSAF2YR
#count_meds
#BMXSAD1
#section_I
#PAQ715
#BMXBMI
#WHD120
#section_E
#DBQ197

DEMO_I2015_2016 <- read.xport("DEMO_I2015_2016.XPT")
DEMO_I2015_2016 <- DEMO_I2015_2016[c('SEQN','RIDAGEYR')]

TRIGLY_I2015_2016 <- read.xport("TRIGLY_I2015_2016.XPT")
TRIGLY_I2015_2016 <- TRIGLY_I2015_2016[c('SEQN','WTSAF2YR')]

BMX_I2015_2016 <- read.xport("BMX_I2015_2016.XPT")
BMX_I2015_2016 <- BMX_I2015_2016[c('SEQN','BMXSAD1','BMXBMI')]

PAQ_I2015_2016 <- read.xport("PAQ_I2015_2016.XPT")
PAQ_I2015_2016 <- PAQ_I2015_2016[c('SEQN','PAQ715')]

WHQ_I2015_2016 <- read.xport("WHQ_I2015_2016.XPT")
WHQ_I2015_2016 <- WHQ_I2015_2016[c('SEQN','WHD120')]

DBQ_I2015_2016 <- read.xport("DBQ_I2015_2016.XPT")
DBQ_I2015_2016 <- DBQ_I2015_2016[c('SEQN','DBQ197')]

BPX_I2015_2016 <- read.xport("BPX_I2015_2016.XPT")
BPX_I2015_2016 <- BPX_I2015_2016[c('SEQN','BPXSY1','BPXDI1')]

RXQ_RX_I2015_2016 <- read.xport("RXQ_RX_I2015_2016.XPT")
RXQ_RX_I2015_2016 <- RXQ_RX_I2015_2016 %>% 
  group_by(SEQN) %>% 
  summarize(
    count_meds = n() 
  )

data2015_2016 <- list(DEMO_I2015_2016, TRIGLY_I2015_2016, BMX_I2015_2016, PAQ_I2015_2016, WHQ_I2015_2016, DBQ_I2015_2016, BPX_I2015_2016, RXQ_RX_I2015_2016) 
data2015_2016 <- data2015_2016 %>% 
  reduce(left_join, by='SEQN') %>% 
  drop_na(BPXSY1)
systolic <- c(ifelse(data2015_2016['BPXSY1'] > 140, 1, 0))
data2015_2016 <- mutate(data2015_2016, systolic = systolic)
head(data2015_2016)

set.seed(1)
inTrain <- sample(nrow(data2015_2016), nrow(data2015_2016)*.7)
traindata <- data2015_2016[inTrain, ]
testdata <- data2015_2016[-inTrain, ]

m_2015_2016_1 <- glm(systolic ~ RIDAGEYR+WTSAF2YR+count_meds+BMXSAD1+PAQ715+BMXBMI+WHD120+DBQ197,
                     family = 'binomial', data = traindata)
summary(m_2015_2016_1)

m_2015_2016_2 <- glm(systolic ~ RIDAGEYR+WTSAF2YR+BMXSAD1+PAQ715,
                     family = 'binomial', data = traindata)
summary(m_2015_2016_2)
probabilities <- predict(m_2015_2016_2, type="response", newdata = testdata)
y_pred <- ifelse(probabilities > 0.5, 1, 0)
confusionMatrix(as.factor(y_pred), as.factor(testdata$systolic)) #0.8405 


#logistic regression - cross section: 1999 - 2000, 2009 - 2010, 2017 - 2018 DEMO dataset

#SEQN: Respondent sequence number
#RIAGENDR: Gender
#RIDAGEYR: Age at Screening Adjudicated
#RIDRETH1: Race/Ethnicity
#DMQMILIT: Veteran/Military Status
#DMDMARTL: Marital Status
#DMDHHSIZ: Total number of people in the Household
#INDHHINC: Annual Household Income
#RIDEXPRG: Pregnancy Status at Exam

#1999 - 2000 DEMO dataset
DEMO1999_2000 <- read.xport("DEMO1999_2000.XPT")
DEMO1999_2000 <- DEMO1999_2000[c('SEQN','RIAGENDR','RIDAGEYR','RIDRETH1','DMQMILIT','DMDMARTL',
                                 'DMDHHSIZ','INDHHINC','RIDEXPRG')]
BPX1999_2000 <- read.xport("BPX1999_2000.XPT")
BPX1999_2000 <- BPX1999_2000[c('SEQN','BPXSAR','BPXDAR')]

data1999_2000 <- merge(DEMO1999_2000, BPX1999_2000, by = 'SEQN')

systolic <- c(ifelse(data1999_2000['BPXSAR'] > 140, 1, 0))
data1999_2000 <- mutate(data1999_2000, systolic = systolic)
head(data1999_2000)

set.seed(1)
inTrain <- sample(nrow(data1999_2000), nrow(data1999_2000)*.7)
traindata <- data1999_2000[inTrain, ]
testdata <- data1999_2000[-inTrain, ]

m_1999_2000_1 <- glm(systolic ~ RIAGENDR+RIDAGEYR+RIDRETH1+DMQMILIT+DMDMARTL+DMDHHSIZ+INDHHINC+RIDEXPRG,
         family = 'binomial', data = traindata)
summary(m_1999_2000_1)

m_1999_2000_2 <- glm(systolic ~ RIAGENDR+RIDAGEYR+RIDRETH1+DMDHHSIZ,
         family = 'binomial', data = traindata)
summary(m_1999_2000_2)
probabilities <- predict(m_1999_2000_2, type="response", newdata = testdata)
y_pred <- ifelse(probabilities > 0.5, 1, 0)
confusionMatrix(as.factor(y_pred), as.factor(testdata$systolic)) #0.8821  

m_1999_2000_3 <- glm(systolic ~ RIDAGEYR,
                     family = 'binomial', data = traindata)
summary(m_1999_2000_3)
probabilities <- predict(m_1999_2000_3, type="response", newdata = testdata)
y_pred <- ifelse(probabilities > 0.5, 1, 0)
confusionMatrix(as.factor(y_pred), as.factor(testdata$systolic)) #0.8803

#2009 - 2010 DEMO dataset
DEMO2009_2010 <- read.xport("DEMO2009_2010.XPT")
DEMO2009_2010 <- DEMO2009_2010[c('SEQN','RIAGENDR','RIDAGEYR','RIDRETH1','DMQMILIT','DMDMARTL',
                                 'DMDHHSIZ','INDHHIN2','RIDEXPRG')] #10537 * 9
BPX2009_2010 <- read.xport("BPX2009_2010.XPT")
BPX2009_2010 <- BPX2009_2010[c('SEQN','BPXSY4','BPXDI4')] #10253 * 3

data2009_2010 <- merge(DEMO2009_2010, BPX2009_2010, by = 'SEQN')
systolic <- c(ifelse(data2009_2010['BPXSY4'] > 140, 1, 0))
data2009_2010 <- mutate(data2009_2010, systolic = systolic)

set.seed(1)
inTrain <- sample(nrow(data2009_2010), nrow(data2009_2010)*.7)
traindata <- data2009_2010[inTrain, ]
testdata <- data2009_2010[-inTrain, ]

m_2009_2010_1 <- glm(systolic ~ RIAGENDR+RIDAGEYR+RIDRETH1+DMQMILIT+DMDMARTL+DMDHHSIZ+INDHHIN2+RIDEXPRG,
                     family = 'binomial', data = traindata)
summary(m_2009_2010_1)

m_2009_2010_2 <- glm(systolic ~ RIAGENDR+RIDAGEYR+RIDRETH1+DMDHHSIZ,
                     family = 'binomial', data = traindata)
summary(m_2009_2010_2)
probabilities <- predict(m_2009_2010_2, type="response", newdata = testdata)
y_pred <- ifelse(probabilities > 0.5, 1, 0)
confusionMatrix(as.factor(y_pred), as.factor(testdata$systolic)) #0.7874

m_2009_2010_3 <- glm(systolic ~ RIDAGEYR,
                     family = 'binomial', data = traindata)
summary(m_2009_2010_3)
probabilities <- predict(m_2009_2010_3, type="response", newdata = testdata)
y_pred <- ifelse(probabilities > 0.5, 1, 0)
confusionMatrix(as.factor(y_pred), as.factor(testdata$systolic)) #0.7874

#2017 - 2018 DEMO dataset
DEMO2017_2018 <- read.xport("DEMO2017_2018.XPT")
DEMO2017_2018 <- DEMO2017_2018[c('SEQN','RIAGENDR','RIDAGEYR','RIDRETH1','DMQMILIZ','DMDMARTL',
                                 'DMDHHSIZ','INDHHIN2','RIDEXPRG')] #9254 * 9
BPX2017_2018 <- read.xport("BPX2017_2018.XPT")
BPX2017_2018 <- BPX2017_2018[c('SEQN','BPXSY4','BPXDI4')] #10253 * 3

data2017_2018 <- merge(DEMO2017_2018, BPX2017_2018, by = 'SEQN')
systolic <- c(ifelse(data2017_2018['BPXSY4'] > 140, 1, 0))
data2017_2018 <- mutate(data2017_2018, systolic = systolic)

set.seed(1)
inTrain <- sample(nrow(data2017_2018), nrow(data2017_2018)*.7)
traindata <- data2017_2018[inTrain, ]
testdata <- data2017_2018[-inTrain, ]

m_2017_2018_1 <- glm(systolic ~ RIAGENDR+RIDAGEYR+RIDRETH1+DMQMILIZ+DMDMARTL+DMDHHSIZ+INDHHIN2+RIDEXPRG,
                     family = 'binomial', data = traindata)
summary(m_2017_2018_1)

m_2017_2018_2 <- glm(systolic ~ RIAGENDR+RIDAGEYR+RIDRETH1+DMDHHSIZ,
                     family = 'binomial', data = traindata)
summary(m_2017_2018_2)
probabilities <- predict(m_2017_2018_2, type="response", newdata = testdata)
y_pred <- ifelse(probabilities > 0.5, 1, 0)
confusionMatrix(as.factor(y_pred), as.factor(testdata$systolic)) #0.7458

m_2017_2018_3 <- glm(systolic ~ RIDAGEYR,
                     family = 'binomial', data = traindata)
summary(m_2017_2018_3)
probabilities <- predict(m_2017_2018_3, type="response", newdata = testdata)
y_pred <- ifelse(probabilities > 0.5, 1, 0)
confusionMatrix(as.factor(y_pred), as.factor(testdata$systolic)) #0.7458