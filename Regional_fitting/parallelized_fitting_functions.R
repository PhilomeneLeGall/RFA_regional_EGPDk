######################################################################
##### Regional inference of EGPD on ERA5 precipitation data (EU) #####
#####                      All Gridpoints                        #####
#####       Philomene Le Gall, Pauline Rivoire 01.03.2021        #####
######################################################################



# required libraries ------------------------------------------------------
library(mev)
library(evd)#pgpd
library(plyr)#rbind.fill.matrix
library(zipfR)#incomplete beta function
library(foreach);library(iterators);library(parallel);library(doParallel) # parallelization


# PRELIMINARY FUNCTIONS ---------------------------------------------------
IB <- function(x,y,a,b){
  #incomplete beta as defined in appendix of Naveau et al., 2016
  z = Ibeta(y,a,b)-Ibeta(x,a,b)
  return(z)
}
H <- function(x,sigma,xi){
  z = pgpd(x,scale=sigma,shape=xi)
  return(z)
}
pG = function(kappa,u){
  fct=u^kappa# = G(u)
  return(fct)
}
Fbar <- function(x,kappa,sigma,xi){
  z = 1-pG(kappa,H(x,sigma,xi))
  return(z)
}




# Parallelized semiregional fit -------------------------------------------

fitEGPDkSemiRegCensoredIter_paral <- function(M,ncores,method="pwm",cens_thres=c(0,Inf),round=0.1,
                                              thres=0, precision=.0001,
                                              ParInit=c(0.5,0.5,0.2), loop.max = 10){
  #   DESCRIPTION
  # Fit a SEMI regional version of EGPD in Naveau et al. (2016).
  #   ARGUMENTS
  # M, matrix of positive precipitation 
  # Each row corresponds to a day, each col to a site
  # censor_thres, bounds to estimate parameters
  # ParReg, method for regionalizing shape parameter xi : 
  #         either "mean" (i.e. mean of the at-site estimates)
  #         or "pool" (pool normalized data to estimate regional xi)
  # thres = threshold of wet days (e.g. 2mm) if not null, provides conditional parameters
  
  #   VALUE
  # List of parameters (Theta,Theta_0)
  #   - the fist element of the list is the list of semi-regional flexibility parameter kappa,
  #     semi-regional scale parameters sigma and regional shape parameter xi.
  #   - the second element of the list contains initialization of parameters (i.e. at-site estimates)
  ###############################################
  #if only one temporal serie put in a matrix
  
  registerDoParallel(cores=ncores)
  
  if(is.null(dim(M))){
    M = as.matrix(M,ncol=1)}
  nday = nrow(M)
  nstat = ncol(M)
  
  Sigma_names=Kappa_names=rep(NA,nstat)
  CensoredMean = moyN  = rep(NA,nstat)
  #consider only days with data
  obs_rr = matrix(nrow=0,ncol=0)
  obs_norm = matrix(0,nrow = 0,ncol =0)
  Theta = list()
  Theta_init = list()
  kappa_init = sigma_init = xisite_init= rep(NA,nstat)
  for (station in 1:nstat){
    Sigma_names[station] = paste("sigma",station,sep="_")
    Kappa_names[station] = paste("kappa",station,sep="_")
    
    y = M[,station]
    y = y[!is.na(y)]
    y = y[y>thres]-thres#consider only positive precipitation when thres =0
    moy = mean(y,na.rm=TRUE)
    obs_rr = rbind.fill.matrix(obs_rr,t(as.matrix(y)))#and fill the matrix with NA if needed
    obs_norm = rbind.fill.matrix(obs_norm,t(as.matrix(y/moy)))#and fill the matrix with NA if needed
  }#end loop on station for data formatting
  M = t(obs_rr)
  MN=t(obs_norm)
  #definir Theta et sa forme (matrice ? liste ? dataframe ?)
  
  Theta_list <- foreach(station=1:nstat) %dopar% {
    y =na.omit(M[,station])
    fit_init = fit.extgp(y,model=1,method=method,init =ParInit, censoring = cens_thres, rounded = round,plots = FALSE)
    kappa_init = fit_init$fit$pwm[1]
    sigma_init = fit_init$fit$pwm[2]
    xisite_init = fit_init$fit$pwm[3]
    
    return(c(kappa_init, sigma_init, xisite_init))
  } #end foreach station 
  
  Theta_mat0 <- matrix(unlist(Theta_list), ncol = length(Theta_list))
  
  kappa_init <- Theta_mat0[1,]
  sigma_init <- Theta_mat0[2,]
  xisite_init <- Theta_mat0[3,]
  
  xi_init = mean(xisite_init)
  Theta_init = list("kappa"=kappa_init,"sigma"=sigma_init,"xi.reg"=rep(xi_init,nstat))
  
  Theta = Theta_init
  Theta_0 = list("kappa"=kappa_init,"sigma"=sigma_init,"xi.site"=xisite_init,"xi.reg"=rep(xi_init,nstat))
  loop <-0; increment <-precision +.1 
  log.init<-0
  Xi.old <- unique(Theta$xi.reg)
  u = cens_thres[1]
  while((loop<loop.max)&(increment>precision)){
    loop <- loop+1
    Sigma.old <- Theta$sigma
    Kappa.old <- Theta$kappa
    
    CensoredMean_list <- foreach(station=1:nstat) %dopar% {
      y = na.omit(M[,station])
      return(mean(y[y>u],na.rm=TRUE))
    } #end foreach station 
    
    CensoredMean <- unlist(CensoredMean_list)
    
    Sigma.new <- (Xi.old*CensoredMean)/(Kappa.old*IB(H(u,Sigma.old,Xi.old),
                                                     1,Kappa.old,1-Xi.old)/Fbar(u,Kappa.old,
                                                                                Sigma.old,Xi.old)-1)
    Theta$sigma <- Sigma.new
    
    moyN_list <- foreach(station=1:nstat) %dopar% {
      y = na.omit(M[,station])
      return(mean(y[y/Sigma.new[station]>u/Sigma.new[station]]/Sigma.new[station],na.rm=TRUE))
    } #end foreach station
    
    moyN <- unlist(moyN_list) # at site censored mean (data normalized by sigma.new)
    
    Kappa.new = (Xi.old*moyN)/(1*IB(H(u,Sigma.new,Xi.old),
                                    1,Kappa.old,1-Xi.old)/Fbar(u,Kappa.old,
                                                               Sigma.new,Xi.old)-1/Kappa.old)
    Theta$kappa = Kappa.new
    #update kappa and sigma in Theta
    
    increment<-max(abs(Sigma.old-Theta$sigma),abs(Kappa.old-Theta$kappa))
    
  }#end while
  
  
  return(list("Theta"=Theta,"Theta_loc"=Theta_0))
}# end function fitEGPDkSemiRegCensoredIter_paral




# Paralelized regional fit ------------------------------------------------
fitEGPDkREGCensorIter_paral <- function(M, ncores, method="pwm",cens_thres=c(0,Inf),round=0.1,
                                        thres=0, precision=.0001,
                                        ParInit=c(0.5,0.5,0.2), loop.max = 10){
  #   DESCRIPTION
  # Fit a REGIONAL  version of EGPD in Naveau et al. (2016).
  #   ARGUMENTS
  # M, matrix of positive precipitation 
  # Each row corresponds to a day, each col to a site
  # censor_thres, bounds to estimate parameters
  # ParReg, method for regionalizing shape parameter xi : 
  #         either "mean" (i.e. mean of the at-site estimates)
  #         or "pool" (pool normalized data to estimate regional xi)
  # thres = threshold of wet days (e.g. 2mm) if not null, provides conditional parameters
  
  #   VALUE
  # List of parameters: (Theta,Theta_0)
  #   - the fist element of the list is the list of regional kappa, semi-regional scale parameters sigma_i, and regional shape parameter xi
  #   - the second element of the list contains initialization parameters (at-site estimates)
  ###############################################
  #if only one temporal serie put in a matrix
  registerDoParallel(cores=ncores)
  
  if(is.null(dim(M))){
    M = as.matrix(M,ncol=1)}
  nday = nrow(M)
  nstat = ncol(M)
  
  Sigma_names=Kappa_names=rep(NA,nstat)
  CensoredMean = moyN  = rep(NA,nstat)
  #consider only days with data
  obs_rr = matrix(nrow=0,ncol=0)
  obs_norm = matrix(0,nrow = 0,ncol =0)
  Theta = list()
  Theta_init = list()
  kappasite_init = sigma_init = xisite_init= rep(NA,nstat)
  for (station in 1:nstat){
    Sigma_names[station] = paste("sigma",station,sep="_")
    
    y = M[,station]
    y = y[!is.na(y)]
    y = y[y>thres]-thres#consider only positive precipitation when thres =0
    moy = mean(y,na.rm=TRUE)
    obs_rr = rbind.fill.matrix(obs_rr,t(as.matrix(y)))#and fill the matrix with NA if needed
    obs_norm = rbind.fill.matrix(obs_norm,t(as.matrix(y/moy)))#and fill the matrix with NA if needed
  }#end loop on station for data formatting
  M = t(obs_rr)
  MN=t(obs_norm)
  #definir Theta et sa forme (matrice ? liste ? dataframe ?)
  
  Theta_list <- foreach(station=1:nstat) %dopar% {
    y =na.omit(M[,station])
    fit_init = fit.extgp(y,model=1,method=method,init =ParInit, censoring = cens_thres, rounded = round,plots = FALSE)
    kappasite_init = fit_init$fit$pwm[1]
    sigma_init = fit_init$fit$pwm[2]
    xisite_init = fit_init$fit$pwm[3]
    return(c(kappasite_init, sigma_init, xisite_init))
  }#end foreach station for initialization
  Theta_loc_mat <- matrix(unlist(Theta_list), ncol = length(Theta_list))
  Theta_loc=list("kappai" = Theta_loc_mat[1,],"sigmai"=Theta_loc_mat[2,],"xii"=Theta_loc_mat[3,])
  
  kappa_init = mean(Theta_loc$kappai)
  xi_init = mean(Theta_loc$xii)
  sigma_init = Theta_loc$sigmai
  Theta_init = list("kappa.reg"=rep(kappa_init,nstat),"sigma"=sigma_init,"xi.reg"=rep(xi_init,nstat))
  
  Theta = Theta_init
  Theta_0 = list("kappa.reg"=rep(kappa_init,nstat),"sigma"=sigma_init,"xi.reg"=rep(xi_init,nstat),
                 "kappa.site"=Theta_loc$kappai,"xi.site"=Theta_loc$xii)
  loop <-0; increment <-precision +.1 
  log.init<-0
  Xi.old <- unique(Theta$xi.reg)
  Kappa.old = unique(Theta$kappa.reg)
  u = cens_thres[1]
  while((loop<loop.max)&(increment>precision)){
    loop <- loop+1
    Sigma.old <- Theta$sigma
    
    CensoredMean_list <- foreach(station=1:nstat) %dopar% {
      y = na.omit(M[,station])
      return(mean(y[y>u],na.rm=TRUE))
    } #end foreach station 
    CensoredMean <- unlist(CensoredMean_list)#at-site censored mean 
    
    
    Sigma.new <- (Xi.old*CensoredMean)/(Kappa.old*IB(H(u,Sigma.old,Xi.old),
                                                     1,Kappa.old,1-Xi.old)/Fbar(u,Kappa.old,
                                                                                Sigma.old,Xi.old)-1)
    Theta$sigma <- Sigma.new
    
    increment<-max(abs(Sigma.old-Theta$sigma))
    
  }#end while
  
  
  return(list("Theta"=Theta,"Theta_loc"=Theta_loc))
} #end function fitEGPDkREGCensorIter_paral



# Function for parallelized bootstrap fit ---------------------------------

fitEGPDkSemiReg.boot_paral <- function(M,ncores,method="pwm",cens_thres=c(0,Inf),round=0.1,
                                       sites = "all", thres=0, precision=.0001,
                                       ParInit=c(0.5,0.5,0.2), loop.max = 10){
  #   DESCRIPTION
  # Fit a (semi)regional version of EGPD for chosen sites in a same cluster in Naveau et al. (2016).
  #   ARGUMENTS
  # M, matrix of positive precipitation 
  # Each row corresponds to a day, each col to a site
  # censor_thres, bounds to estimate parameters
  # sites, indices of GP/stations where local parameters are computed (e.g. medoid, min/max silhouette)
  # thres = threshold of wet days (e.g. 2mm) if not null, provides conditional parameters
  
  #   VALUE
  # List of parameters (Theta,Theta_0)
  #   - the fist element of the list is the list of semi-regional flexibility parameter kappa,
  #     semi-regional scale parameters sigma and regional shape parameter xi.
  #   - the second element of the list contains initialization of parameters (i.e. at-site estimates)
  ###############################################
  #if only one temporal serie put in a matrix
  
  registerDoParallel(cores=ncores)
  
  if(is.null(dim(M))){
    M = as.matrix(M,ncol=1)}
  nday = nrow(M)
  nstat = ncol(M)
  if(is.character(sites)){if(sites=="all"){sites = 1:nstat}else{stop("undefined option for choice of sites")} }
  nb_stat = length(sites)
  #parameters for all sites
  SigmaAll_names=KappaAll_names=rep(NA,nstat)
  CensoredMean = moyN  = rep(NA,nb_stat)
  #consider only days with data
  obs_rr = matrix(nrow=0,ncol=0)
  obs_norm = matrix(0,nrow = 0,ncol =0)
  Theta = list()
  Theta_init = list()
  kappaAll_init = sigmaAll_init = xisiteAll_init= rep(NA,nstat)
  for (station in 1:nstat){
    SigmaAll_names[station] = paste("sigma",station,sep="_")
    KappaAll_names[station] = paste("kappa",station,sep="_")
    
    y = M[,station]
    y = y[!is.na(y)]
    y = y[y>thres]-thres#consider only positive precipitation when thres =0
    moy = mean(y,na.rm=TRUE)
    obs_rr = rbind.fill.matrix(obs_rr,t(as.matrix(y)))#and fill the matrix with NA if needed
    obs_norm = rbind.fill.matrix(obs_norm,t(as.matrix(y/moy)))#and fill the matrix with NA if needed
  }#end loop on station for data formatting
  M = t(obs_rr)
  
  
  
  Theta_list <- foreach(station=1:nstat) %dopar% {
    y =na.omit(M[,station])
    fit_init = fit.extgp(y,model=1,method=method,init =ParInit, censoring = cens_thres, rounded = round,plots = FALSE)
    kappaAll_init = fit_init$fit$pwm[1]
    sigmaAll_init = fit_init$fit$pwm[2]
    xisiteAll_init = fit_init$fit$pwm[3]
    
    return(c(kappaAll_init, sigmaAll_init, xisiteAll_init))
  } #end foreach station for initialization
  
  Theta_mat0 <- matrix(unlist(Theta_list), ncol = length(Theta_list))
  
  kappaAll_init <- Theta_mat0[1,]
  sigmaAll_init <- Theta_mat0[2,]
  xisiteAll_init <- Theta_mat0[3,]
  
  
  xi_init = mean(xisiteAll_init)
  kappa_init = kappaAll_init[sites]; sigma_init = sigmaAll_init[sites]; xisite_init = xisiteAll_init[sites]
  Theta_init = list("kappa"=kappa_init,"sigma"=sigma_init,"xi.reg"=rep(xi_init,nb_stat))
  
  Theta = Theta_init
  Theta_0 = list("kappa"=kappa_init,"sigma"=sigma_init,"xi.site"=xisite_init,"xi.reg"=rep(xi_init,nb_stat))
  loop <-0; increment <-precision +.1 
  log.init<-0
  Xi.old <- unique(Theta$xi.reg)
  u = cens_thres[1]
  while((loop<loop.max)&(increment>precision)){
    loop <- loop+1
    Sigma.old <- Theta$sigma
    Kappa.old <- Theta$kappa
    
    CensoredMean_list <- foreach(station=1:nb_stat) %dopar% {
      site = sites[station]; y = na.omit(M[,site])
      return(mean(y[y>u],na.rm=TRUE))
    } #end foreach station 
    CensoredMean <- unlist(CensoredMean_list) #at-site censored mean
    
    Sigma.new <- (Xi.old*CensoredMean)/(Kappa.old*IB(H(u,Sigma.old,Xi.old),
                                                     1,Kappa.old,1-Xi.old)/Fbar(u,Kappa.old,
                                                                                Sigma.old,Xi.old)-1)
    Theta$sigma <- Sigma.new
    for (i in 1:nb_stat){
      site = sites[i]; y = na.omit(M[,site])
      moyN[i] = mean(y[y/Sigma.new[i]>u/Sigma.new[i]]/Sigma.new[i],na.rm=TRUE)}# at site censored mean (data normalized by sigma.new)
    
    Kappa.new = (Xi.old*moyN)/(1*IB(H(u,Sigma.new,Xi.old),
                                    1,Kappa.old,1-Xi.old)/Fbar(u,Kappa.old,
                                                               Sigma.new,Xi.old)-1/Kappa.old)
    Theta$kappa = Kappa.new
    #update kappa and sigma in Theta
    
    increment<-max(abs(Sigma.old-Theta$sigma),abs(Kappa.old-Theta$kappa))
    
  }#end while
  
  
  return(list("Theta"=Theta,"Theta_0"=Theta_0))
} #end function fitEGPDkSemiReg.boot_paral


fitEGPDkREG.boot_paral <- function(M, ncores,method="pwm",cens_thres=c(0,Inf),round=0.1,
                                   sites = "all", thres=0, precision=.0001,
                                   ParInit=c(0.5,0.5,0.2), loop.max = 10){
  #   DESCRIPTION
  # Fit a regional version of EGPD for chosen sites in a same cluster in Naveau et al. (2016).
  #   ARGUMENTS
  # M, matrix of positive precipitation 
  # Each row corresponds to a day, each col to a site
  # censor_thres, bounds to estimate parameters
  # sites, indices of GP/stations where local parameters are computed (e.g. medoid, min/max silhouette)
  # thres = threshold of wet days (e.g. 2mm) if not null, provides conditional parameters
  
  #   VALUE
  # List of parameters (Theta,Theta_0)
  #   - the fist element of the list is the regional flexibility parameter kappa,
  #     semi-regional scale parameters sigma and regional shape parameter xi.
  #   - the second element of the list contains initialization of parameters (i.e. at-site estimates)
  ###############################################
  #if only one temporal serie put in a matrix
  
  registerDoParallel(cores=ncores)
  
  if(is.null(dim(M))){
    M = as.matrix(M,ncol=1)}
  nday = nrow(M)
  nstat = ncol(M)
  if(is.character(sites)){if(sites=="all"){sites = 1:nstat}else{stop("undefined option for choice of sites")} }
  nb_stat = length(sites)
  #parameters for all sites
  SigmaAll_names=KappaAll_names=rep(NA,nstat)
  CensoredMean = moyN  = rep(NA,nb_stat)
  #consider only days with data
  obs_rr = matrix(nrow=0,ncol=0)
  obs_norm = matrix(0,nrow = 0,ncol =0)
  Theta = list()
  Theta_init = list()
  kappaAll_init = sigmaAll_init = xisiteAll_init= rep(NA,nstat)
  for (station in 1:nstat){
    SigmaAll_names[station] = paste("sigma",station,sep="_")
    KappaAll_names[station] = paste("kappa",station,sep="_")
    
    y = M[,station]
    y = y[!is.na(y)]
    y = y[y>thres]-thres#consider only positive precipitation when thres =0
    moy = mean(y,na.rm=TRUE)
    obs_rr = rbind.fill.matrix(obs_rr,t(as.matrix(y)))#and fill the matrix with NA if needed
    obs_norm = rbind.fill.matrix(obs_norm,t(as.matrix(y/moy)))#and fill the matrix with NA if needed
  }#end loop on station for data formatting
  M = t(obs_rr)
  
  
  Theta_list <- foreach(station=1:nstat) %dopar% {
    y =na.omit(M[,station])
    fit_init = fit.extgp(y,model=1,method=method,init =ParInit, censoring = cens_thres, rounded = round,plots = FALSE)
    kappaAll_init = fit_init$fit$pwm[1]
    sigmaAll_init = fit_init$fit$pwm[2]
    xisiteAll_init = fit_init$fit$pwm[3]
    
    return(c(kappaAll_init, sigmaAll_init, xisiteAll_init))
  } #end foreach station  for initialization
  
  Theta_mat0 <- matrix(unlist(Theta_list), ncol = length(Theta_list))
  
  kappaAll_init <- Theta_mat0[1,]
  sigmaAll_init <- Theta_mat0[2,]
  xisiteAll_init <- Theta_mat0[3,]
  
  
  xi_init = mean(xisiteAll_init); kappa_init = mean(kappaAll_init)
  kappasite_init = kappaAll_init[sites]; sigma_init = sigmaAll_init[sites]; xisite_init = xisiteAll_init[sites]
  Theta_init = list("kappa"=rep(kappa_init,nb_stat),"sigma"=sigma_init,"xi.reg"=rep(xi_init,nb_stat))
  
  Theta = Theta_init
  Theta_0 = list("kappa.site"=kappasite_init,"sigma"=sigma_init,"xi.site"=xisite_init)
  loop <-0; increment <-precision +.1 
  log.init<-0
  Xi.old <- unique(Theta$xi.reg)
  Kappa.old <- unique(Theta$kappa)
  u = cens_thres[1]
  while((loop<loop.max)&(increment>precision)){
    loop <- loop+1
    Sigma.old <- Theta$sigma
    #compute sigmai only for sites of interest
    CensoredMean_list <- foreach(station=1:nb_stat) %dopar% {
      site = sites[station]; y = na.omit(M[,site])
      return(mean(y[y>u],na.rm=TRUE))
    } #end foreach station at-site censored mean 
    
    CensoredMean <- unlist(CensoredMean_list)
    
    Sigma.new <- (Xi.old*CensoredMean)/(Kappa.old*IB(H(u,Sigma.old,Xi.old),
                                                     1,Kappa.old,1-Xi.old)/Fbar(u,Kappa.old,
                                                                                Sigma.old,Xi.old)-1)
    Theta$sigma <- Sigma.new
    
    
    increment<-max(abs(Sigma.old-Theta$sigma))
    
  }#end while
  
  
  return(list("Theta"=Theta,"Theta_0"=Theta_0))
}# end function fitEGPDkREG.boot_paral






# Extract positive seasonal precipitation --------------------------------------------------------

# Description
# 
# Extract from a precip time serie the non NA and >1 precip of a given season 
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




# Fit both semireg and reg at the same time -------------------------------

fitEGPDk.boot_paral <- function(M,ncores,method="pwm",cens_thres=c(0,Inf),round=0.1,
                                sites = "all", thres=0, precision=.0001,
                                ParInit=c(0.5,0.5,0.2), loop.max = 10){
  #   DESCRIPTION
  # Fit local, semi-regional and regional versions of EGPD for chosen sites in a same cluster, see Naveau et al. (2016) and Le Gall et al. (2021)
  #   ARGUMENTS
  # M, matrix of positive precipitation for sites in a same homogeneous region 
  # Each row corresponds to a day, each column to a site
  # censor_thres, bounds to estimate parameters
  # sites, indices of GP/stations where local parameters are computed (e.g. medoid, min/max silhouette)
  # thres = threshold of wet days (e.g. 2mm) if not null, provides conditional parameters
  
  #   VALUE
  # List of list of parameters (Theta_reg, Theta_semireg,Theta_0)
  #   - the fist element of the first list is the regional flexibility parameter kappa, the second element is the vector of
  #     semi-regional scale parameters sigma, the last element is the regional shape parameter xi.
  #  - the fist element of the second list is the vector of semi-regional flexibility parameters kappa, the second element is the vector of
  #     semi-regional scale parameters sigma, the last element is the regional shape parameter xi.
  #   - the second element of the list contains initialization of parameters (i.e. at-site estimates)
  ###############################################
  #if only one temporal serie put in a matrix
  
  registerDoParallel(cores=ncores)
  
  if(is.null(dim(M))){
    M = as.matrix(M,ncol=1)}
  nday = nrow(M)
  nstat = ncol(M)
  if(is.character(sites)){if(sites=="all"){sites = 1:nstat}else{stop("undefined option for choice of sites")} }
  nb_stat = length(sites)
  #parameters for all sites
  SigmaAll_names=KappaAll_names=rep(NA,nstat)
  CensoredMean = moyN  = rep(NA,nb_stat)
  #consider only days with data
  obs_rr = matrix(nrow=0,ncol=0)
  obs_norm = matrix(0,nrow = 0,ncol =0)
  ThetaR = ThetaSR = list()
  Theta_init = list()
  kappaAll_init = sigmaAll_init = xisiteAll_init= rep(NA,nstat)
  for (station in 1:nstat){
    SigmaAll_names[station] = paste("sigma",station,sep="_")
    KappaAll_names[station] = paste("kappa",station,sep="_")
    
    y = M[,station]
    y = y[!is.na(y)]
    y = y[y>thres]-thres#consider only positive precipitation when thres =0
    moy = mean(y,na.rm=TRUE)
    obs_rr = rbind.fill.matrix(obs_rr,t(as.matrix(y)))#and fill the matrix with NA if needed
    obs_norm = rbind.fill.matrix(obs_norm,t(as.matrix(y/moy)))#and fill the matrix with NA if needed
  }#end loop on station for data formatting
  M = t(obs_rr)
  
  
  
  Theta_list <- foreach(station=1:nstat) %dopar% {
    library(mev)
    y = na.omit(M[,station])
    fit_init = fit.extgp(y,model=1,method=method,init =ParInit, censoring = cens_thres, rounded = round,plots = FALSE)
    kappaAll_init = fit_init$fit$pwm[1]
    sigmaAll_init = fit_init$fit$pwm[2]
    xisiteAll_init = fit_init$fit$pwm[3]
    return(c(kappaAll_init, sigmaAll_init, xisiteAll_init))
  } #end foreach station  for initialization
  Theta_mat0 <- matrix(unlist(Theta_list), ncol = length(Theta_list))
  kappaAll_init <- Theta_mat0[1,]
  sigmaAll_init <- Theta_mat0[2,]
  xisiteAll_init <- Theta_mat0[3,]
  
  xi_init = mean(xisiteAll_init); kappa_init = mean(kappaAll_init)
  kappasite_init = kappaAll_init[sites]; sigma_init = sigmaAll_init[sites]; xisite_init = xisiteAll_init[sites]
  Theta_init = list("kappa"=rep(kappa_init,nb_stat),"sigma"=sigma_init,"xi.reg"=rep(xi_init,nb_stat))
  
  ThetaR = ThetaSR = Theta_init
  Theta_0 = list("kappa.site"=kappasite_init,"sigma"=sigma_init,"xi.site"=xisite_init)
  loop <-0; increment <-precision +.1 
  log.init<-0
  XiR.old <- unique(ThetaR$xi.reg); XiSR.old <- unique(ThetaSR$xi.reg)
  
  KappaR.old <- unique(ThetaR$kappa)
  u = cens_thres[1]
  
  CensoredMean_list <- foreach(station=1:nb_stat) %dopar% {
    site = sites[station]; y = na.omit(M[,site])
    return(mean(y[y>u],na.rm=TRUE))
  } #end foreach station at-site censored mean 
  CensoredMean <- unlist(CensoredMean_list)
  
  
  #loop for regional fit
  while((loop<loop.max)&(increment>precision)){
    loop <- loop+1
    SigmaR.old <- ThetaR$sigma
    #compute sigmai only for sites of interest
    SigmaR.new <- (XiR.old*CensoredMean)/(KappaR.old*IB(H(u,SigmaR.old,XiR.old),
                                                        1,KappaR.old,1-XiR.old)/Fbar(u,KappaR.old,
                                                                                     SigmaR.old,XiR.old)-1)
    ThetaR$sigma <- SigmaR.new
    ######################
    SigmaSR.old <- ThetaSR$sigma
    KappaSR.old <- ThetaSR$kappa
    SigmaSR.new <- (XiSR.old*CensoredMean)/(KappaSR.old*IB(H(u,SigmaSR.old,XiSR.old),
                                                           1,KappaSR.old,1-XiSR.old)/Fbar(u,KappaSR.old,
                                                                                          SigmaSR.old,XiSR.old)-1)
    ThetaSR$sigma <- SigmaSR.new
    for (i in 1:nb_stat){
      site = sites[i]; y = na.omit(M[,site])
      moyN[i] = mean(y[y/SigmaSR.new[i]>u/SigmaSR.new[i]]/SigmaSR.new[i],na.rm=TRUE)}# at site censored mean (data normalized by sigma.new)
    
    KappaSR.new = (XiSR.old*moyN)/(1*IB(H(u,SigmaSR.new,XiSR.old),
                                        1,KappaSR.old,1-XiSR.old)/Fbar(u,KappaSR.old,
                                                                       SigmaSR.new,XiSR.old)-1/KappaSR.old)
    ThetaSR$kappa = KappaSR.new
    
    
    #####################
    increment<-max(abs(SigmaR.old-ThetaR$sigma),abs(SigmaSR.old-ThetaSR$sigma),abs(KappaSR.old-ThetaSR$kappa))
    
  }#end while
  
  
  
  return(list("Theta_reg"=ThetaR,"Theta_semireg" = ThetaSR,"Theta_0"=Theta_0))
}#end fitEGPDk.boot_paral()







