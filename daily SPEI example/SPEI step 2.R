# Note I am not sure if all those packages are need, but install them just in case

library("boot")
library("dplR")
library("dplyr")
library("lubridate")
library("reshape2")
library("dendroTools")
library("stringi")
library("stringr")
library("SPEI")
library("lmom")
library("zoo")
library("forcats")
library("TLMoments")
library("SPEI")
library("lmom")


# 1 An example of TRW data
chronology <- read.crn("TRW.crn")

# Extracts information from analysis and converts into a txt file 
step2_SPEI <- function(input){
  
  analysed_period_start <- list()
  analysed_period_end <- list()
  calculation <- list()
  
  calculation_lower <- list()
  calculation_upper <- list()
  
  w_length <- list()
  opti_start<- list()
  previous_year <- list()
  
  analysed_period_start[[i]] <- substr(as.character(input[4 ,2]), 1, 4)
  analysed_period_end[[i]] <- substr(as.character(input[4 ,2]), 8, 11)
  
  calculation[[i]] <- as.numeric(as.character(input[5, 2]))
  calculation_lower[[i]] <- as.numeric(as.character(input[6, 2]))
  calculation_upper[[i]]  <- as.numeric(as.character(input[7, 2]))
  
  w_length[[i]] <-  as.numeric(as.character(input[11, 2]))
  opti_start[[i]] <- as.numeric(str_extract_all(as.character(input[8, 2]),"\\(?[0-9,.]+\\)?")[[1]][1])
  previous_year[[i]] <- grepl("Previous", as.character(input[8, 2]))
  
  final_df <- data.frame(matrix(as.numeric(unlist(analysed_period_start))),
                         matrix(as.numeric(unlist(analysed_period_end))),
                         matrix(as.numeric(unlist(calculation))),
                         
                         matrix(as.numeric(unlist(calculation_lower))),
                         matrix(as.numeric(unlist(calculation_upper))),
                         
                         matrix(as.numeric(unlist(w_length))),
                         matrix(as.numeric(unlist(opti_start))),
                         matrix(unlist(previous_year))
  )
  
  colnames(final_df) <- c("start", "end", "calculation", "lower_bound", "upper_bound" ,"w_length", "opti_start", "previous_year")
  
  final_df <- dplyr::mutate(final_df, opti_end = opti_start + w_length - 1)
  
  final_df <- dplyr::select(final_df, start, end, calculation, lower_bound, upper_bound , w_length, opti_start, opti_end, previous_year)
  
  final_df
}

# This function is copy pasted from GitHub repository
# https://github.com/sbegueria/SPEI
# It is slightly modified to fit in the daily_response()
# All credits goes to authors
SPEI_daily <- function(SPEI_source_file, scale = 21){
  
  SPIEa <- SPEI_source_file[,-c(1,2, 3) ]
  SPIEa = zooreg(SPIEa, start=as.Date("1950-01-01"))
  SPIEa= as.ts(SPIEa) 
  data = SPIEa
  
  kernel=list(type='rectangular',shift=0)
  distribution='log-Logistic' 
  fit='ub-pwm'
  na.rm=FALSE
  ref.start=NULL
  ref.end=NULL
  x=FALSE
  params=NULL
  
  scale <- as.numeric(scale)
  na.rm <- as.logical(na.rm)
  x <- as.logical(x)
  #if (!exists("data",inherits=F) | !exists("scale",inherits=F)) {
  #	stop('Both data and scale must be provided')
  #}
  if (!(distribution %in% c('log-Logistic', 'Gamma', 'PearsonIII'))) {
    stop('Distrib must be one of "log-Logistic", "Gamma" or "PearsonIII"')
  }
  if (!(fit %in% c('max-lik', 'ub-pwm', 'pp-pwm'))) {
    stop('Method must be one of "ub-pwm" (default), "pp-pwm" or "max-lik"')
  }
  if ( (!is.null(ref.start) && length(ref.start)!=2) | (!is.null(ref.end) && length(ref.end)!=2) ) {
    stop('Start and end of the reference period must be a numeric vector of length two.')
  }
  
  if (!is.ts(data)) {
    data <- ts(as.matrix(data), frequency = 12, start = c(1950, 1))
  } else {
    data <- ts(as.matrix(data), frequency=frequency(data), start=start(data))
  }
  m <- ncol(data)
  fr <- frequency(data)
  
  
  coef = switch(distribution,
                "Gamma" = array(NA,c(2,m,fr),list(par=c('alpha','beta'),colnames(data),NULL)),
                "log-Logistic" = array(NA,c(3,m,fr),list(par=c('xi','alpha','kappa'),colnames(data),NULL)),
                "PearsonIII" = coef <- array(NA,c(3,m,fr),list(par=c('mu','sigma','gamma'),colnames(data),NULL))
  )
  
  dim_one = ifelse(distribution == "Gamma", 2, 3)
  
  if (!is.null(params)) {
    if (dim(params)[1]!=dim_one | dim(params)[2]!=m | dim(params)[3]!=12) {
      stop(paste('parameters array should have dimensions (', dim_one, ', ', m, ', 12)',sep=' '))
    }
  }
  
  # Loop through series (columns in data)
  if (!is.null(ref.start) && !is.null(ref.end)) {
    data.fit <- window(data,ref.start,ref.end)	
  } else {
    data.fit <- data
  }
  std <- data*NA
  for (s in 1:m) {
    # Cumulative series (acu)
    acu <- data.fit[,s]
    acu.pred <- data[,s]
    if (scale>1) {
      wgt <- kern(scale,kernel$type,kernel$shift)
      acu[scale:length(acu)] <- rowSums(embed(acu,scale)*wgt,na.rm=na.rm)
      acu[1:(scale-1)] <- NA
      acu.pred[scale:length(acu.pred)] <- rowSums(embed(acu.pred,scale)*wgt,na.rm=na.rm)
      acu.pred[1:(scale-1)] <- NA
    }
    
    # Loop through the months
    for (c in (1:fr)) {
      # Filter month m, excluding NAs
      f <- which(cycle(acu)==c)
      f <- f[!is.na(acu[f])]
      ff <- which(cycle(acu.pred)==c)
      ff <- ff[!is.na(acu.pred[ff])]
      
      # Monthly series, sorted
      month <- sort.default(acu[f], method="quick")
      
      if (length(month)==0) {
        std[f] <- NA
        next()
      }
      
      if (is.null(params)) {
        month_sd = sd(month,na.rm=TRUE)
        if (is.na(month_sd) || (month_sd == 0)) {
          std[f] <- NA
          next
        }
        
        if(distribution != "log-Logistic"){
          pze <- sum(month==0)/length(month)
          month = month[month > 0]
        }
        
        # Stop early and assign NAs if month's data is length < 4
        if(length(month) < 4){
          std[ff,s] = NA
          coef[,s,c] <- NA
          next
        }
        
        # Calculate probability weighted moments based on fit with lmomco or TLMoments
        pwm = switch(fit,
                     "pp-pwm" = pwm.pp(month,-0.35,0, nmom=3),
                     #pwm.ub(month, nmom=3)
                     TLMoments::PWM(month, order=0:2)
        )
        
        # Check L-moments validity
        lmom <- pwm2lmom(pwm)
        if ( !are.lmom.valid(lmom) || anyNA(lmom[[1]]) || any(is.nan(lmom[[1]])) ){
          next
        }
        
        # lmom fortran functions need specific inputs L1, L2, T3
        # this is handled by lmomco internally with lmorph
        fortran_vec = c(lmom$lambdas[1:2], lmom$ratios[3])
        
        # Calculate parameters based on distribution with lmom then lmomco
        f_params = switch(distribution,
                          "log-Logistic" = tryCatch(lmom::pelglo(fortran_vec), error = function(e){ parglo(lmom)$para }),
                          "Gamma" = tryCatch(lmom::pelgam(fortran_vec), error = function(e){ pargam(lmom)$para }),
                          "PearsonIII" = tryCatch(lmom::pelpe3(fortran_vec), error = function(e){ parpe3(lmom)$para })
        )
        
        # Adjust if user chose log-Logistic and max-lik
        if(distribution == 'log-Logistic' && fit=='max-lik'){
          f_params = parglo.maxlik(month, f_params)$para
        }
      } else {
        
        f_params = as.vector(params[,s,c])
        
      }
      
      # Calculate cdf based on distribution with lmom
      cdf_res = switch(distribution,
                       "log-Logistic" = lmom::cdfglo(acu.pred[ff], f_params),
                       "Gamma" = lmom::cdfgam(acu.pred[ff], f_params),
                       "PearsonIII" = lmom::cdfpe3(acu.pred[ff], f_params)				  				
      )
      
      std[ff,s] = qnorm(cdf_res)
      coef[,s,c] <- f_params
      
      # Adjust if user chose Gamma or PearsonIII
      if(distribution != 'log-Logistic'){ 
        std[ff,s] = qnorm(pze + (1-pze)*pnorm(std[ff,s]))
      }
      
    } # next c (month)
  } # next s (series)
  colnames(std) <- colnames(data)
  
  z <- list(call=match.call(expand.dots=FALSE),
            fitted=std,coefficients=coef,scale=scale,kernel=list(type=kernel$type,
                                                                 shift=kernel$shift,values=kern(scale,kernel$type,kernel$shift)),
            distribution=distribution,fit=fit,na.action=na.rm)
  if (x) z$data <- data
  if (!is.null(ref.start)) z$ref.period <- rbind(ref.start,ref.end)
  
  SPEI_source_file[,4] <- z$fitted
  
  SPEI_source_file
}

# 
SPEI_f <- read.table( "SPEI_PAD.txt", header = TRUE) # This is PET - Potential evapotranspiration for each day in a year

# You have to define the lower and upper limits, this is the range of days to be aggreagted
lower_limit = 21
upper_limit = 270

previous_year = FALSE
response = chronology
boot = FALSE
cor_method = 'pearson'

temporal_matrix_list <- list()
b2 = 1

for (i in (lower_limit:upper_limit)){

    
    SPEI_temp <- SPEI_daily(SPEI_source_file = SPEI_f, scale = i)

    vector <- SPEI_temp[, 4]
    
    nNA <- sum(is.na(vector))
    
    if (nNA > i - 1){
      
      stop("NA problem!!!!")
      
    }
    
    SPEI_temp <- SPEI_temp[c(1:(nrow(SPEI_temp)-nNA)),c(1,2,3)]
    
    SPEI_temp$values <- vector[(nNA+1):length(vector)]
    
    SPEIX_daily <- data_transform(SPEI_temp[,c(3,4)])
    
    env_data = SPEIX_daily
    
    if (previous_year == TRUE) {
      
      # FIRST, both data frames need to be arranged, the most recent year is the first one
      env_data$yearABC <- row.names(env_data)
      env_data <- dplyr::arrange(env_data, desc(yearABC))
      env_data <- years_to_rownames(env_data, "yearABC")
      env_data_previous <- env_data[-1, , F]
      env_data_current <- env_data[-nrow(env_data), ,F]
      row_names_current <- row.names(env_data_current)
      env_data <- cbind(env_data_previous, env_data_current)
      env_data <- data.frame(env_data)
      row.names(env_data) <- row_names_current
      env_data_original <- env_data
      
      #response$yearABC <- row.names(response)
      #response <- dplyr::arrange(response, desc(yearABC))
      #response <- years_to_rownames(response, "yearABC")
      #response <- data.frame(response[-nrow(response),,F ])
      #response <- data.frame(response)
      #response_original <- response
      
    }
    
  
    ncol_response <- ncol(response)
    
    colnames_response <- colnames(response)
    
    env_data$yearABC <- row.names(env_data)
    response$yearABC <- row.names(response)
    
    temporal_data <- merge(response, env_data, by = "yearABC")
    
    response <- data.frame(temporal_data[, c(2:(1 + ncol_response))],
                           row.names = temporal_data$yearABC)
    colnames(response) <- colnames_response
    
    env_data <- data.frame(temporal_data[, c((1 + ncol_response + 1):
                                               ncol(temporal_data))],
                           row.names = temporal_data$yearABC)
    
    
    temporal_matrix <- matrix(NA, nrow = 1,
                              ncol = (ncol(env_data) + lower_limit) + 1)
    # Here I create two additional temporal matrices, which will be used to store
    # lower and upper limits of bootstrap estimates
    temporal_matrix_lower <- temporal_matrix
    temporal_matrix_upper <- temporal_matrix
    
    for (j in 0: (ncol(env_data) - i)) {
      
      x <- env_data[, (j + 1)]
      
      x <- matrix(x, nrow = nrow(env_data), ncol = 1)
      
      if (boot == FALSE){
        temporal_correlation <- cor(response[, 1], x[, 1], method = cor_method)
        temporal_lower <- NA
        temporal_upper <- NA
      } else if (boot == TRUE){
        temp_df_boot <- cbind(response[, 1], x[, 1])
        
        temp_df_boot <- temp_df_boot[complete.cases(temp_df_boot), ]
        
        calc <- boot(temp_df_boot, boot_f, fun = "cor", cor.type = cor_method, R = boot_n)
        
        temporal_correlation <- colMeans(calc$t)[1]
        
        ci_int <- boot.ci(calc, conf = boot_conf_int, type = boot_ci_type)
        temporal_lower <- ci_int$norm[2]
        temporal_upper <- ci_int$norm[3]
      } else {
        print(paste0("boot should be TRUE or FALSE, instead it is ", boot))
      }
      
      temporal_matrix[1, j + 1] <- temporal_correlation
      #temporal_matrix_lower[1, j + 1] <- temporal_lower
      #temporal_matrix_upper[1, j + 1] <- temporal_upper
      
    }
    
    temporal_matrix_list[[b2]] <- temporal_matrix
    #temporal_matrix_lower_list[[b]] <- temporal_matrix_lower
    #temporal_matrix_upper_list[[b]] <- temporal_matrix_upper
    b2 = b2 +1
    
  }
  
  temporal_matrix1 <- do.call(rbind, temporal_matrix_list)
  temporal_rownames <- as.vector(seq(from = lower_limit, to = upper_limit, by = 1))
  row.names(temporal_matrix1) <- temporal_rownames

  temporal_colnames <- as.vector(seq(from = 1, to = ncol(temporal_matrix), by = 1))
  colnames(temporal_matrix1) <- temporal_colnames

  #########################################################################################
  #########################################################################################
  
  temp_result <- data.frame(temporal_matrix1)
  
heat <- temp_result
heat$temp_row_names <- seq(lower_limit, upper_limit, by = 1)
heat <- melt(heat, id.vars = c("temp_row_names"))

# Define the threshold for correlations to be 
heat$value <- ifelse(abs(heat$value) < 0.40, NA, heat$value)

ggplot(heat,aes_(x = ~as.numeric(variable), y = ~as.numeric(temp_row_names), fill = ~value)) +
  geom_tile() +
  xlab("Day of Year") +
  ylab("Season Length") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0),breaks = c(90, 180, 270)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", na.value = 'gray97', midpoint = 0, limits = c(-1, 1)) + 
  theme_minimal() +  theme(axis.text = element_text(size = 10),
                           axis.title.y = element_text(size = 18), text = element_text(size = 18),
                           axis.title.x = element_blank(),
                           plot.title = element_text(size = 16),
                           legend.title = element_blank(), 
                           legend.position = "bottom", legend.key.width = unit(3, "line"),
                           panel.background = element_rect(fill = "gray97",
                                                           colour = "gray80",
                                                           size = 0.5, linetype = "solid"))

ggsave("Example_spei.png", height = 8, width = 10)
