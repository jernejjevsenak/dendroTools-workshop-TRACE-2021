# load dendroTools
library("dendroTools")

# open TRW chronology swit272
data("swit272")

# Open daily temperature data for swit272
data("swit272_daily_temperatures")
swit272_daily_temperatures <- data_transform(swit272_daily_temperatures, 
                                             date_format = "ymd",
                                             format = "daily")

# Find the optimal climate signal
cor_run <- daily_response(response = swit272, 
                                   env_data = swit272_daily_temperatures, 
                                   row_names_subset = TRUE, method = "cor",
                                   temporal_stability_check = "running_window")

# 0 check the summary output
summary(cor_run)

# 1 check temporal stability
cor_run$temporal_stability

# 2 check the transfer function
cor_run$transfer_function

# 3 sample split test (cross-validation)
cor_run$cross_validation

# 4 Fit linear model
head(cor_run$optimized_return) # we are going to use the optimized_return output
colnames(cor_run$optimized_return)[2] <- "T_Jun18_Aug16" # change the column name

# Fit linear model
lm_model <- lm(T_Jun18_Aug16 ~ TRWi, data = cor_run$optimized_return)
summary(lm_model)

# Reconstruct climate using the swit272 chronology
reconstruction <- data.frame(linear = predict(lm_model, swit272))
reconstruction$year <- as.numeric(row.names(swit272)) # extract year

# plot the reconstructed temperatures 
library(ggplot2)
ggplot(reconstruction, aes(x = year, y = linear)) + geom_line(col = "blue") +
  scale_x_continuous(limits = c(1750, 2011)) +
  theme_bw() + xlab("Year") + ylab("Mean Temperature Jun 18 - Aug 16") 


########## Add nonlinear reconstruction using the brnn function ###########

# Fit nonlinear BRNN model
library(brnn)
brnn_model <- brnn(T_Jun18_Aug16 ~ TRWi, data = cor_run$optimized_return)

# Nonlinear reconstruction
reconstruction$brnn <- predict(brnn_model, swit272) # predict using the swit272 chronology

# Prepare data for ggplot
library(reshape2)
reconstruction <- melt(reconstruction, id.vars = "year")

# plot linear and non-linear reconstructed temperatures 
ggplot(reconstruction, aes(x = year, y = value, col = variable)) + geom_line() +
  scale_x_continuous(limits = c(1750, 2011)) +
  theme_bw() + labs(col = 'Model') +
  xlab("Year") + ylab("Mean Temperature Jun 18 - Aug 16")

