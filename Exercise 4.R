# load dendroTools
library("dendroTools")

# open TRW chronology 
data("swit272")

# Open daily temperature and precipitation data
# (typical file with daily data from KNMI Climate explorer)
data("swit272_daily_temperatures")
swit272_monthly_temperatures <- data_transform(swit272_daily_temperatures, 
                                             date_format = "ymd",
                                             format = "monthly",
                                             monthly_aggregate_function = "auto")

data("swit272_daily_precipitation")
swit272_monthly_precipitation <- data_transform(swit272_daily_precipitation, 
                                              date_format = "ymd",
                                              format = "monthly",
                                              monthly_aggregate_function = "auto")


# 1 monthly_response()
simple_cor_monthly <- monthly_response(response = swit272, 
                                env_data = swit272_monthly_temperatures,
                                row_names_subset = TRUE, method = "cor")

summary(simple_cor_monthly)
plot(simple_cor_monthly, type = 2)


# 2 monthly_response_seascorr()
partial_cor_monthly <- monthly_response_seascorr(response = swit272,
                                       env_data_primary = swit272_monthly_temperatures,
                                       env_data_control = swit272_monthly_precipitation,
                                       row_names_subset = TRUE, remove_insignificant = FALSE)

summary(partial_cor_monthly)
plot(partial_cor_monthly, type = 2)


###############################################
# What is the benefit of applying daily data? #
###############################################

# Open daily temperature file (typical file with daily data from KNMI Climate explorer)
data("swit272_daily_temperatures")
swit272_daily_temperatures <- data_transform(swit272_daily_temperatures, 
                                             date_format = "ymd",
                                             format = "daily")

# Run basic example
simple_cor_daily <- daily_response(response = swit272, 
                             env_data = swit272_daily_temperatures, 
                             row_names_subset = TRUE, method = "cor")

summary(simple_cor_daily)
summary(simple_cor_monthly)
