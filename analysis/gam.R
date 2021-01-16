library(mgcv)
library(tidyverse)
library(lubridate)
library(spatstat)

# load functions in package
devtools::load_all("chigcrim/")

# Use 2019 data
crime_2016 <- load_data(2016)
crime_2017 <- load_data(2017)
crime_2018 <- load_data(2018)

crime_train <- rbind(crime_2016, crime_2017, crime_2018)

crime_2019 <- load_data(2019)

crime_test <- crime_2019

# Extract useful information from the date
crime_train <- crime_train %>%
  mutate(month = factor(month(ymd_hms(date, tz = "GMT"))),
         day = factor(day(ymd_hms(date, tz = "GMT"))),
         hour = factor(hour(ymd_hms(date, tz = "GMT"))),
         yday = (yday(ymd_hms(date))),
         date = date(ymd_hms(date))) 


# Convert rarely used factors into 'other'
crime_train$primary_type <- otherise(as.character(crime_train$primary_type), 5000)

crime_train$location_description <- otherise(as.character(crime_train$location_description), 5000)

# List of features to use in GAM
keep_features <- c("primary_type", "location_description", "arrest", 
                   "domestic", "month", "year", "day","date", "community_area", "fbi_code", "yday")



# Prep data for use in GAM
crime_train <- crime_train %>% select(keep_features) %>% na.omit(crime_dat) %>% 
  mutate(fbi_code = as.factor(fbi_code)) %>% mutate(community_area = as.factor(community_area)) %>%
  mutate(day = as.numeric(day))

crime_train

# Make column counting number of crimes
crime_train$n <- rep(1, nrow(crime_train))

train_dat <- aggregate(n ~ community_area + yday + year + fbi_code + date, sum, data = crime_train)

gam1 <- gam(n ~ s(as.numeric(yday), bs = "cc") + as.factor(community_area) + as.factor(fbi_code) + 
              as.factor(year), data = train_dat, family = "poisson")

predict(gam1, type = "response") - train_dat$n
plot(train_dat$n, predict(gam1, type = "response"), pch = 16, col = rgb(1,1,1, alpha = 0.5))
abline(a = 0, b = 1, col = "red")

add_n_pre <- function(dat, start_date) {
  for (ca in unique(dat$community_area)) {
    for (fbi in unique(dat$fbi_code)) {
      current_dat <-dat[dat[,"community_area"] == ca & dat[,"fbi_code"] == fbi,]
      current_dat$n_pre <- left_join(data.frame(date = current_dat$date - days(1)), current_dat, by = "date")$n
      current_dat$n_pre[is.na(current_dat$n_pre)] <- 0
      current_dat$n_pre[current_dat$date == start_date] <- NA
      dat[dat[,"community_area"] == ca & dat[,"fbi_code"] == fbi,"n_pre"] <- current_dat$n_pre
    }
  }
  return(dat)
  
}

new_data <- add_n_pre(train_dat,"2016-01-01")

gam2 <- gam(n ~ as.numeric(n_pre) + s(as.numeric(yday), bs = "cc") + 
              as.factor(year) + as.factor(fbi_code) + as.factor(community_area), data = new_data)
max(abs(predict(gam2, type = "response") - new_data[new_data$date != "2016-01-01",]$n))
plot(new_data[new_data$date != "2016-01-01",]$n, predict(gam2, type = "response"))
abline(a = 0, b = 1, col = "red")
summary(gam2)