# load dplR 
library("dendroTools")

# open TRW chronology 
data("swit272")

# Open daily temperature file (typical file with daily data from KNMI Climate explorer)
data("swit272_daily_temperatures")
swit272_daily_temperatures <- data_transform(swit272_daily_temperatures, 
                                             date_format = "ymd",
                                             format = "daily")

# Open daily temperature file (typical file with daily data from KNMI Climate explorer)
data("swit272_daily_precipitation")
swit272_daily_precipitation <- data_transform(swit272_daily_precipitation, 
                                             date_format = "ymd",
                                             format = "daily")
# Run basic example
simple_cor <- daily_response(response = swit272, 
                                env_data = swit272_daily_temperatures, 
                                row_names_subset = TRUE, method = "cor")

partial_cor <- daily_response_seascorr(response = swit272,
                                       env_data_primary = swit272_daily_temperatures,
                                       env_data_control = swit272_daily_precipitation,
                                       row_names_subset = TRUE)

library(ggplot2)
plot(simple_cor, type = 2) + scale_fill_gradient2()
plot(partial_cor, type = 2) + scale_fill_gradient2()

?scale_fill_gradient2

summary(simple_cor)
summary(partial_cor)
