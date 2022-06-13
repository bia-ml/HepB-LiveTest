# HepB-LiveTest
---
Title: "The Development of a Machine Learning Algorithm for Early Detection of Viral Hepatitis B Infection in Nigerian patients"
author: "Busayo Ajuwon"
date: "17 August 2021"
output:
 output:
  pdf_document:
    latex_engine: xelatex
    number_sections: yes
    toc: yes
    toc_depth: 2
---

require(ggplot2)
require(plotly)
require(psych)
require(rattle)
library(caret)


**MACHINE LEARNING EXPERIMENT**
*Name of Dataset: Allcases**

READ THE DATASET
```{r}
library(readr)
Allcases <- read_csv("D:/PhD Studies@ANU/Predicting HBsAg  in Oz/Allcases.csv")
View(Allcases)
```

DESCRIPTIVE ANALYSIS 
Any NAsin the data?
```{r}
length(which(is.na(Allcases)))
```

How many positive and negative hepatitis cases are in the data set?
```{r}
Allcases$HBSA <- as.factor(Allcases$HBSA)
summary (Allcases$HBSA)
```

Summarize data by group
```{r}
by(Allcases, Allcases$HBSA, summary)
```


DATA EXPLORATION
```{r}
ggplot(Allcases, aes(x = HBSA, fill = HBSA)) +
  geom_bar()
```

#CONTINUE WITH DATA EXPLORATION
```{r}
featurePlot(x = Allcases[, 1:20],
y = Allcases$HBSA,
plot = "box",
strip=strip.custom(par.strip.text=list(cex=.9)),
scales = list(x = list(relation="free"),
y = list(relation="free")))
```

EXAMINE THEIR DENSITY PLOTS
```{r}
featurePlot(x = Allcases[, 1:20],
y = Allcases$HBSA,
plot = "density",
strip=strip.custom(par.strip.text=list(cex=.7)),
scales = list(x = list(relation="free"),
y = list(relation="free")))
```

*Feature selection using recursive feature elimination (rfe)*
```{r warning=FALSE}
set.seed(100)
subsets <- c(1:5, 10, 15, 20)

ctrl <- rfeControl(functions = rfFuncs,
method = "repeatedcv",
repeats = 5,
verbose = FALSE)

lmProfile <- rfe(x=Allcases[, 1:20], y=Allcases$HBSA,
sizes = subsets,
rfeControl = ctrl)
lmProfile
```

**TRAINING AND TUNING THE MODEL**
Train the model and interpret the results
```{r}
set.seed(100)
# Step 1: Get row numbers for the training data
trainRowNumbers <- createDataPartition(Allcases$HBSA, p=0.7, list=FALSE)

# Step 2: Create the training dataset
trainData <- Allcases[trainRowNumbers,]

# Step 3: Create the test dataset
testData <- Allcases[-trainRowNumbers,]

```{r}
# Store X and Y for later use.
x = trainData[, 1:20]
y = trainData$HBSA
```


*Run RF algorithm on patient data and rank predictive markers according to importance*
*Define the training control*
```{r}
fitControl <- trainControl(
  method = 'repeatedcv',            # repeated cross validation
  number = 10,                      # number of folds
  repeats = 10,                    # number of repeats 
  savePredictions = 'final',       # saves predictions for optimal tuning parameter
  classProbs = T,                  # should class probabilities be returned
  summaryFunction=twoClassSummary  # results summary function
)
```

Hyperparameter tuning using tune length I take the train() function that was before, plus, additionally set the tuneLength, trControl and metric.
```{r}
#Tune hyperparameters by setting tuneLength 
set.seed(100)
model_rf = train(HBSA ~ ., data=trainData, method='rf', tuneLength=5, trControl = fitControl)
model_rf
```

```{r}
#Step 2: Predict on testData
predicted <- predict(model_rf, testData)
head(predicted)
```


```{r}
#Compute the confusion matrix
confusionMatrix(reference = testData$HBSA, data = predicted, mode='everything', positive='positive')
```

Compute variable importance
```{r}
varimp_rf <- varImp(model_rf)
plot(varimp_rf, main="Variable Importance with RF")
```

Establish decision thresholds using top-ranked predictive markers from random forest 
```{r}
set.seed(100)
model_dt = train(HBSA ~ AST + WBC + Age + ALT , data=trainData, method='rpart', tuneLength=5, trControl = fitControl)
model_dt

fancyRpartPlot(model_dt$finalModel, sub = '')
```

Random Forest Analysis of HBsAg Immunoassay Results 
```{r}
set.seed(100)
model_rf = train(HBSA ~ AST + WBC + Age + ALT , data=trainData, method='rf', tuneLength=5, trControl = fitControl)
model_rf
```

```{r}
#Step 2: Predict on testData
predicted <- predict(model_rf, testData)
head(predicted)
```

```{r}
#Compute the confusion matrix
confusionMatrix(reference = testData$HBSA, data = predicted, mode='everything', positive='positive')
```

Support Vector Machine (SVM) Analysis of HBsAg Immunoassay Results
*Train the model using Support Vector Machines*
```{r}
set.seed(100)
model_svmRadial = train(HBSA ~ AST + WBC + Age + ALT + ALB, data=trainData, method='svmRadial', tuneLength=5, trControl = fitControl)
model_svmRadial
```

```{r}
#Predict on test data
predicted <- predict(model_svmRadial, testData)
head(predicted)
```

```{r}
#Compute the confusion matrix
confusionMatrix(reference = testData$HBSA, data = predicted, mode='everything', positive='positive')
```

SVM PLOTS
ALT Kinetics Associated with the Primary Predictor Variables AST and WBC
Introduce ALT as static slice
**ALT = 20 U/L**

```{r}
m1 <- svm (HBSA~ AST + WBC + ALT, data = testData, method = 'svmRadial', gamma = 0.36, cost= 4)
plot(m1, testData, AST ~ WBC, 
     slice = list (ALT = 20))
```

**ALT = 30 U/L**
```{r}
m1 <- svm (HBSA~ AST + WBC + ALT, data = testData, method = 'svmRadial', gamma = 0.36, cost= 4)
plot(m1, testData, AST ~ WBC, 
     slice = list (ALT = 30))
```

**ALT = 50 U/L **
```{r}
m1 <- svm (HBSA~ AST + WBC + ALT, data = testData, method = 'svmRadial', gamma = 0.36, cost= 4)
plot(m1, testData, AST ~ WBC, 
     slice = list (ALT = 50))
```


**ALT = 100 U/L **
```{r}
m1 <- svm (HBSA~ AST + WBC + ALT, data = testData, method = 'svmRadial', gamma = 0.36, cost= 4)
plot(m1, testData, AST ~ WBC, 
     slice = list (ALT = 100))
```


**ALT = 1000 U/L **
```{r}
m1 <- svm (HBSA~ AST + WBC + ALT, data = testData, method = 'svmRadial', gamma = 0.36, cost= 4)
plot(m1, testData, AST ~ WBC, 
     slice = list (ALT = 1000))
```

Impact of Age on the SVM Prediction of HBsAg Immunoassay Result by AST and WBC
Introduce age as a static slice

```{r}
#Age = 15 years 
m1 <- svm (HBSA~ AST + WBC + Age, data = testData, method = 'svmRadial', gamma = 0.36, cost= 4)
plot(m1, testData, AST ~ WBC, 
     slice = list (Age = 15))
```

```{r}
#Age = 25 years 
m1 <- svm (HBSA~ AST + WBC + Age, data = testData, method = 'svmRadial', gamma = 0.36, cost= 4)
plot(m1, testData, AST ~ WBC, 
     slice = list (Age = 25))
```

```{r}
#Age = 35 years 
m1 <- svm (HBSA~ AST + WBC + Age, data = testData, method = 'svmRadial', gamma = 0.36, cost= 4)
plot(m1, testData, AST ~ WBC, 
     slice = list (Age = 35))
```

```{r}
#Age = 45 years 
m1 <- svm (HBSA~ AST + WBC + Age, data = testData, method = 'svmRadial', gamma = 0.36, cost= 4)
plot(m1, testData, AST ~ WBC, 
     slice = list (Age = 45))
```

```{r}
#Age = 55 years 
m1 <- svm (HBSA~ AST + WBC + Age, data = testData, method = 'svmRadial', gamma = 0.36, cost= 4)
plot(m1, testData, AST ~ WBC, 
     slice = list (Age = 55))
```

```{r}
#Age = 65 years 
m1 <- svm (HBSA~ AST + WBC + Age, data = testData, method = 'svmRadial', gamma = 0.36, cost= 4)
plot(m1, testData, AST ~ WBC, 
     slice = list (Age = 65))
```

#Kinetics for HBsAg positive and HBsAg negative cohorts across the age range investigated by SVM and tree-based machine learning algorithms. 
Comparison of mean WBC (A) and ALB (B) versus age at the time of testing.


WBC KINETICS 


```{r}
ggplot(aes(x = Age, y = WBC, color=HBSA), data = Allcases) +
    stat_summary(fun=mean, geom="line") +
    ylab("Mean WBC ( X 10^9 /L) ") 
```

```{r}
ggplot(aes(x = Age, y = ALB, color=HBSA), data = Allcases) +
    stat_summary(fun=mean, geom="line") +
    ylab("Mean ALB (g/L) ")
```


**RANDOM FOREST MODEL PERFORMANCE COMPARISON AFTER SUBSAMPLING** - subsampling experiments
```{r}
#Model the original unbalanced data
fitControl <- trainControl(
  method = 'repeatedcv',                   
  number = 10,                      
  repeats = 10,
  savePredictions = 'final',     
  classProbs = T,                  
  summaryFunction=twoClassSummary
) 

set.seed(100)
model_rf = train(HBSA ~ AST + WBC + Age + ALT, data=trainData, method='rf', tuneLength=5, trControl = fitControl)
model_rf

#Predict on test data
predicted <- predict(model_rf, testData)
head(predicted)

cm_original <- confusionMatrix(reference = testData$HBSA, data = predicted, mode='everything', positive='positive')
cm_original
```

```{r}
#Undersampling
fitControl <- trainControl(
  method = 'repeatedcv',                   
  number = 10,                      
  repeats = 10,
  savePredictions = 'final',     
  classProbs = T,                  
  summaryFunction=twoClassSummary,
  sampling = 'under'
) 

set.seed(100)
model_rf_under = train(HBSA ~ AST + WBC + Age + ALT, data=trainData, method='rf', tuneLength=5, trControl = fitControl)
model_rf_under

predicted <- predict(model_rf_under, testData)
head(predicted)

cm_under <- confusionMatrix(reference = testData$HBSA, data = predicted, mode='everything', positive='positive')
cm_under
```

```{r}
#Oversampling
fitControl <- trainControl(
  method = 'repeatedcv',                   
  number = 10,                      
  repeats = 10,
  savePredictions = 'final',     
  classProbs = T,                  
  summaryFunction=twoClassSummary,
  sampling = 'up'
) 

set.seed(100)
model_rf_over = train(HBSA ~ AST + WBC + Age + ALT, data=trainData, method='rf', tuneLength=5, trControl = fitControl)
model_rf_over

predicted <- predict(model_rf_over, testData)
head(predicted)

cm_over <- confusionMatrix(reference = testData$HBSA, data = predicted, mode='everything', positive='positive')
cm_over
```

```{r}
#ROSE
fitControl <- trainControl(
  method = 'repeatedcv',                   
  number = 10,  
  repeats = 10,
  savePredictions = 'final',     
  classProbs = T,                  
  summaryFunction=twoClassSummary,  
  sampling = 'rose'
) 

set.seed(100)
model_rf_rose = train(HBSA ~ AST + WBC + Age + ALT, data=trainData, method='rf', tuneLength=5, trControl = fitControl)
model_rf_rose

predicted <- predict(model_rf_rose, testData)
head(predicted)

cm_rose <- confusionMatrix(reference = testData$HBSA, data = predicted, mode='everything', positive='positive')
cm_rose
```

```{r}
#Resample models
models <- list(original = model_rf,
               under = model_rf_under,
               over = model_rf_over,
               rose = model_rf_rose)


resampling <- resamples(models)
bwplot(resampling)
```

```{r}
#Plot a figure to compare- final

library(dplyr)
comparison <- data.frame(model = names(models),
                         Sensitivity = rep(NA, length(models)),
                         Specificity = rep(NA, length(models)),
                         Precision = rep(NA, length(models)),
                         F1 = rep(NA, length(models)))



comparison <- data.frame(model = names(models))

for (name in names(models)) {
  model <- get(paste0("cm_", name))
  
  comparison[comparison$model == name, "Sensitivity"] <- model$byClass[["Sensitivity"]]
  comparison[comparison$model == name, "Specificity"] <- model$byClass[["Specificity"]]
  comparison[comparison$model == name, "Precision"] <- model$byClass[["Precision"]]
  comparison[comparison$model == name, "F1"] <- model$byClass[["F1"]]
}


library(tidyr)
comparison %>%
  gather(x, y, Sensitivity:F1) %>%
  ggplot(aes(x = x, y = y, color = model)) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 3)
