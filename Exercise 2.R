# load dplR 
library("dendroTools")

# open TRW chronology 
data("swit272")

# Open daily temperature file
data("swit272_daily_temperatures")

# Prepare data for daily_response()
swit272_daily_temperatures <- data_transform(swit272_daily_temperatures, 
                                             date_format = "ymd",
                                             format = "daily")

# Test 1: Basic example
# lower_limit = 21, upper = 60, previous_year = TRUE, method = "cor"
t1 <- daily_response(response = swit272, 
                     env_data = swit272_daily_temperatures, 
                     row_names_subset = TRUE, 
                     method = "cor",
                     lower_limit = 21,
                     upper_limit = 60, 
                     previous_year = TRUE)

summary(t1)
plot(t1, type = 1)
plot(t1, type = 2)


# Test 2: Test the subset years option
# subset_years = c(1980, 2010), method = "lm"
t2 <- daily_response(response = swit272, 
                     env_data = swit272_daily_temperatures, 
                     row_names_subset = TRUE, 
                     method = "lm", metric = "adj.r.squared",
                     lower_limit = 21,
                     upper_limit = 60, 
                     previous_year = FALSE,
                     subset_years = c(1980, 2010))
summary(t2)
plot(t2, type = 1)
plot(t2, type = 2)


# Test 3: Use the nonlinear BRNN. Note: This will take some time
# method = "brnn"
t3 <- daily_response(response = swit272, 
                     env_data = swit272_daily_temperatures, 
                     row_names_subset = TRUE, 
                     method = "brnn",
                     lower_limit = 50,
                     upper_limit = 67, 
                     previous_year = FALSE)

summary(t3)
plot(t3, type = 1)
plot(t3, type = 2)
t3$transfer_function

# Test 4: Apply the bootstrap. Note: This will take some time
# boot = TRUE, boot_n = 100, method = "cor"
t4 <- daily_response(response = swit272, 
                     env_data = swit272_daily_temperatures, 
                     row_names_subset = TRUE, 
                     method = "cor",
                     lower_limit = 50,
                     upper_limit = 67, 
                     previous_year = FALSE, 
                     boot = TRUE, boot_n = 1000)

summary(t4)

plot(t4, type = 1)
plot(t4, type = 2)

cal_l <- t4$boot_lower
