# install.packages("dendroTools")

# load dendroTools 
library("dendroTools")

# open swit272 TRW chronology 
data("swit272")

# source ITRDB: Bigler - Sils-Maria GR Blais dal Fö - LADE - ITRDB SWIT272
# doi: https://doi.org/10.25921/37hm-5971
# URL: https://www.ncdc.noaa.gov/paleo-search/?dataTypeId=18

# See the first and last 8 rows of TRW
head(swit272, n = 8)
tail(swit272, n = 8)

# Open daily temperature file (typical file with daily data from KNMI Climate explorer)
data("swit272_daily_temperatures")
head(swit272_daily_temperatures)

# When your daily data has only two columns, i.e. 'date' and and 'climate variable', 
# you can apply data_transform() to prepare data for daily_response()
swit272_daily_temperatures <- data_transform(swit272_daily_temperatures, 
                                             date_format = "ymd",
                                             format = "daily")

dim(swit272_daily_temperatures)

# How can you do this manually?
library("lubridate")
data("swit272_daily_temperatures")
head(swit272_daily_temperatures)

class(swit272_daily_temperatures$date)
summary(swit272_daily_temperatures)

swit272_daily_temperatures$date <- ymd(swit272_daily_temperatures$date) #convert to date 
class(swit272_daily_temperatures$date)

swit272_daily_temperatures$year <- year(swit272_daily_temperatures$date) # extract year from date
head(swit272_daily_temperatures)
tail(swit272_daily_temperatures)

swit272_daily_temperatures$doy <- yday(swit272_daily_temperatures$date) #extract doy/yday from date
head(swit272_daily_temperatures)
tail(swit272_daily_temperatures)

swit272_daily_temperatures$date <- NULL  # delete date
head(swit272_daily_temperatures)
tail(swit272_daily_temperatures)

library("reshape2")
swit272_daily_temperatures <- dcast(year ~ doy, data = swit272_daily_temperatures, value.var = "t_avg")
row.names(swit272_daily_temperatures) <- swit272_daily_temperatures$year # define row.names

swit272_daily_temperatures$year <- NULL # delete year (we don't need it anymore, since it is now included in row.names)
dim(swit272_daily_temperatures) # check dimension -> there should be 366 columns (days)

# Preview your data with glimpse_daily_data() function from dendroTools
glimpse_daily_data(swit272_daily_temperatures, na.color = "black")

# You can customize all ggplots from dendroTools
library("ggplot2")
glimpse_daily_data(swit272_daily_temperatures, na.color = "white") +
  scale_fill_gradient2() + theme(legend.position = "none")
ggsave("example.png")

# Run basic example
basic_example <- daily_response(response = swit272, 
                                env_data = swit272_daily_temperatures, 
                                row_names_subset = TRUE, 
                                method = "cor")
class(basic_example)

summary(basic_example)
plot(basic_example, type = 1)
plot(basic_example, type = 2) 

# modify the legend scale

plot(basic_example, type = 2) + scale_fill_gradient2(expand = c(0, 0), breaks = c(-0.5, 0, 0.5), limits = c(-0.5,0.5)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.key.width = unit(2, "cm"))

# check all outputs 
calc <- basic_example$calculations
calc_l <- basic_example$boot_lower
calc_u <- basic_example$boot_upper

str(basic_example)

basic_example$method
basic_example$metric  
basic_example$analysed_period  
basic_example$type
basic_example$reference_window

basic_example$optimized_return  
basic_example$optimized_return_all
basic_example$transfer_function  + xlab("TRWi") + ylab("Temperature")
basic_example$temporal_stability
basic_example$cross_validation
