#################### FUNCTIONS FOR RFA : CLUSTERING ON PWM RATIO ######################
#################### authors: Philomene Le Gall Pauline Rivoire ######################



# extract seasonal positive precipitation --------------------------------------------------------

# Description
# 
# Extract from a precip time serie the non NA and >threshold precip of a given season 
# 
# Arguments
# season               string in in c("SON", ""DJF", "MAM", "JJA")
# precip_timeserie  	 precipitation time serie to process
# precip_date          dates corresponding to the precip_timeserie (same length)
#
# Value
# vector with non NA and >0 precip of the given season 


extr_seas_pos_precip <- function(season, precip_timeserie, precip_date, thshld = 0){
  stopifnot(length(precip_timeserie)==length(precip_date))
  SEASONS <- c("MAM", "JJA", "SON", "DJF")
  stopifnot(season %in% SEASONS)
  ref_seas <- which(SEASONS==season)
  months_in_season <- (((c(1,2,3)+(3*ref_seas-1))-1)%%12)+1 #gives for ex: c(12,1,2) if season == "DJF"
  
  return(precip_timeserie[(lubridate::month(precip_date)==months_in_season[1]
                           | lubridate::month(precip_date)==months_in_season[2]
                           | lubridate::month(precip_date)==months_in_season[3]) & 
                            (precip_timeserie>thshld)])
}#end for keep_seasonal_positive_precip function




# Estimation omega ratio --------------------------------------------------------------


sampling<-function(data,thres=1){ 
  #   ARGUMENTS
  # data = matrix of observation. Each row corresponds to a day and each column to a site/station
  #   VALUE
  # X, matrix of STRICTLY POSITIVE observation. Each row corresponds to a day and each column to a site/station
  nday <- nrow(data) #number of observation per station
  nb_stat<-dim(thres)[2] #number of stations
  
  data[data<thres]=NA
  X=data-thres
  return(X)}

xi.Ratio <-function(x,independence=TRUE){
  # Arguments
  #   x = temporal serie for a station
  #   n = length of temporal serie
  #
  # Values
  #   xi.R = ratio estimator
  #   
  #   a0, a1, a2, 3 = Estimator of PWM 0,1,2 and 3
  #   See : Diebolt08, Improving PWM method for the GEV distribution; and Guillou, 2009
  #
  #   Apply to precip data matrix to obtain a vector
  #   R.vect of length=nb_stat, vector of the ratio values for each site
  
  
  x<-sort(x,na.last=NA) #remove Na and sort
  
  if(independence==TRUE){# weights for PWM estimators
    weight <- (c(0:(length(x)-1))-.5)/(length(x)-1)#for PWM with cdf
    surv_weight <- 1-weight#for PWM with survival function
  }
  if(independence==FALSE){F=ecdf(x)
  weight <- F(x)
  surv_weight <- 1-weight#for PWM with survival function
  }
  
  
  a0 <- mean(x)
  a1 <- mean(weight*x)
  a2 <- mean(weight^2*x)
  a3 <- mean(weight^3*x)
  
  xi.R <-(3*a2-a0)/(2*a1-a0)
  #if(F) xi.R <- (1/xi.R)
  return(xi.R)
}
R.vec<-function(X){
  #Argument
  #   X=matrix of daily precip, a row = a day, a col = a site
  # Value 
  #   R.vect = vector of R estimates, with length = nb of sites
  
  R.vect <- apply(X, 2,  xi.Ratio)
  return(R.vect)
}