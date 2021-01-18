# Load require packages
library(e1071)
library(tidyverse)
library(RSocrata)
library(tidyr)
library(magrittr)
library(lubridate)
library(caret)
library(LiblineaR)
library(randomForest)


# load functions in package
devtools::load_all("../chigcrim/")

# Use 2019 data
crime_dat <- load_data(2019)

base_crime_dat <- crime_dat
crime_dat <- base_crime_dat
# Extract useful information from the date
crime_dat <- crime_dat %>%
  mutate(month = factor(month(ymd_hms(date, tz = "GMT"))),
         day = factor(day(ymd_hms(date, tz = "GMT"))),
         hour = factor(hour(ymd_hms(date, tz = "GMT")))) %>%
  select(-date)

# Function to convert rarely used factors into an 'other' category
otherise <- function(string_vec, n_threshold, print_summary = TRUE){
  counts <- table(string_vec)
  other_names <- names(counts[counts < n_threshold])
  string_vec[string_vec %in% other_names] <- "OTHER"
  string_vec[is.na(string_vec)] <- "OTHER"
  string_vec
  if (print_summary){
    print(paste(length(other_names), "out of", length(counts),
                "categories were converted to OTHER corresponding to",
                100*length(string_vec[string_vec == "OTHER"])/length(string_vec),
                "% of observations"))
  }
  as.factor(string_vec)
}

# Convert rarely used factors into 'other'
crime_dat$primary_type <- otherise(as.character(crime_dat$primary_type), 5000)
crime_dat$location_description <- otherise(as.character(crime_dat$location_description), 5000)

# List of features to use in SVM
keep_features <- c("primary_type", "location_description", "beat", "arrest", 
                   "domestic", "month", "hour", "community_area", "fbi_code", "longitude", "latitude")

# Prep data for use in SVM
crime_dat <- crime_dat %>% select(keep_features) %>% 
  na.omit(crime_dat) %>% mutate(fbi_code = as.factor(fbi_code))

# Make train and test data (normally want train > test but SVM won't compute with large number of rows)
sub_rows <- sample(1:nrow(crime_dat), 10000)
crime_train <- crime_dat[sub_rows,]
crime_test <- crime_dat[-sub_rows,]

# Create svm using e1071 package
svm1 <- svm(arrest ~ ., data = crime_train, type = "C-classification", kernel = "radial")
svm2 <- svm(arrest ~ latitude + longitude, data = crime_train, type = "C-classification", kernel = "radial")

# Find predictions for each moodel
preds <- predict(svm1, newdata = crime_test)
preds2 <- predict(svm2, newdata = crime_test)

# Analyse performance of SVM
confusionMatrix(preds, as.factor(crime_test$arrest))
confusionMatrix(preds2, as.factor(crime_test$arrest)) # why is this happening??

# Plot how it performs in space
plot(svm1, data = crime_train, formula = longitude ~ latitude)
plot(svm2, data = crime_train, formula = longitude ~ latitude)

# Seem to perform better without factors
crime_train2 <- crime_train %>% mutate(month = as.numeric(month), hour = as.numeric(hour))
crime_test2 <- crime_test %>% mutate(month = as.numeric(month), hour = as.numeric(hour))

svm3 <- svm(arrest ~ ., data = crime_train2, type = "C-classification", kernel = "radial")
preds3 <- predict(svm3, newdata = crime_test2)

confusionMatrix(preds3, as.factor(crime_test$arrest))

# Final model provides best normalized gini score with 0.463
MLmetrics::NormalizedGini(preds, (crime_test$arrest))
MLmetrics::NormalizedGini(preds2, (crime_test$arrest))
MLmetrics::NormalizedGini(preds3, (crime_test$arrest))

# Create a random forest to compare
rf1 <-randomForest(select(crime_train, -"arrest"), y = as.factor(crime_train$arrest))
importance(rf1)

preds4 <- predict(rf1, crime_test)

# random forest outperforms svm
confusionMatrix(preds4, as.factor(crime_test$arrest))


# cv error, trying to make generalisable
cv_error <- function(model, train_data, train_label, error, folds,...) {
  n <- nrow(train_data)
  parts <- split(sample(1:n), 1:folds)
  err = vector(length = folds)
  j = 1
  for (i in parts) {
    cv_train <- train_data[-i,]
    cv_test <- train_data[i,]
    mod1 <- model(data = cv_train, as.formula(paste(train_label, "~ .")), ...)
    preds <- predict(mod1, newdata = cv_test)
    err[j] <- error(preds, as.logical(as.data.frame(cv_test)[,train_label]))
    j <- j + 1
  }
  return(mean(err))
}

# Find CV errors of SVM and RF
CVE_svm <- cv_error(svm, crime_train, "arrest", NormalizedGini, 5, type = "C-classification")
CVE_rf1 <- cv_error(randomForest, crime_train, "arrest", NormalizedGini, 5, ntree = 300, mtry = 5)
CVE_rf2 <- cv_error(randomForest, crime_train, "arrest", NormalizedGini, 5, ntree = 1000, mtry = 5)

#less trees has better results
