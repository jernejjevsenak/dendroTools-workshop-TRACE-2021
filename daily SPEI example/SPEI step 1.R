library("lubridate")

jernejHargreavesPET <- function(date, tavg, tdif, lat){
  
  doy <- day(date)
  phi = pi/180 * lat
  delta = 0.409*sin(2*pi/366*doy-1.39)
  dr = 1 + 0.033*cos(2*pi/365*doy)
  
  ws = try(acos(-tan(phi)*tan(delta)))
  
  if (is.nan(ws) == TRUE){
    ws = 0.1
  }
  
  Ra = (24*60/pi)*0.0820*dr*((ws*sin(phi)*sin(delta))+(cos(phi)*cos(delta)*sin(ws)))
  
  PET = 0.0023*(tavg + 17.8)*sqrt(tdif) * Ra
  
  if (PET < 0){
    PET = 0
  }
  PET
}

# 1 Define the coordinates. The coordinates are important to correctly calculate PET
coordinates <- data.frame(
  key_clim = "DUG",
  latitude = 45.60361111,
  lonngitude = 14.05583333
)

coordinates$tavg_na <- NA
coordinates$psum_na <- NA
coordinates$tmax_na <- NA
coordinates$tmin_na <- NA
coordinates$SPEI_na <- NA

for (i in 1:nrow(coordinates)){

  T_avg <- read.table("Tavg_PAD.txt", header = T)
  P_sum <- read.table("Psum_PAD.txt", header = T)
  T_max <- read.table("Tmax_PAD.txt", header = T)
  T_min <- read.table("Tmin_PAD.txt", header = T)
  
  tavg_na <- sum(is.na(T_avg))/nrow(T_avg)
  psum_na <- sum(is.na(P_sum))/nrow(P_sum)
  tmax_na <- sum(is.na(T_max))/nrow(T_max)
  tmin_na <- sum(is.na(T_min))/nrow(T_min)
  
  coordinates[i,"tavg_na"] <- tavg_na
  coordinates[i,"psum_na"] <- psum_na
  coordinates[i,"tmax_na"] <- tmax_na
  coordinates[i,"tmin_na"] <- tmin_na
  
  
  T_diff <- as.numeric((T_max[ ,4])) - as.numeric((T_min[, 4]))
  sum(T_diff < 0)
  T_diff = ifelse(T_diff < 0, 0, T_diff) # in some rare cases Tmin is greater than Tmax
  
  PET = jernejHargreavesPET(date = as.Date(T_avg[,3]), tavg = T_avg[,4], tdif = T_diff, lat = T_avg[,2])
  
  SPEI = P_sum[,4] - PET
  
  SPEI_index <- T_avg
  SPEI_index[, 4] <- SPEI
  colnames(SPEI_index)[4] <- "SPEI"
  
  SPEI_na <- sum(is.na(SPEI_index[, 4]))/nrow(SPEI_index)
  coordinates[i,"SPEI_na"] <- SPEI_na
  
  write.table(SPEI_index, paste0("SPEI_PAD.txt"))

  
}
