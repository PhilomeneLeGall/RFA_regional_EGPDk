##### IMPROVED RFA FUNCTIONS ####
## Method in Le Gall et al., 2021
library(mev)
library(evd)
library(cluster)#pam
library(evd)#pgpd
library(plyr)#rbind.fill.matrix
library(zipfR)#incomplete beta function
library(foreach);library(iterators);library(parallel);library(doParallel) # parallelization
########################
## CLUSTERING ALGORITHM
########################
Omega.Ratio <- function(x,independence=TRUE){
  # DESCRIPTION
  # Compute the probability Weighted Moments (PWM) scale-invariant ratio described in Le Gall et al., 2021, for a single time serie.
  # ARGUMENTS
  # x               time serie (e.g. positive precipitation at a gridpoint)
  # independence    logical, if FALSE, F is estimated empirically; if TRUE, F(x) is considered as uniform 
  # VALUE
  # omega           scalar, scale invariant ratio of PWM
  # EXAMPLE
  # x = rnorm(n = 100); omega = Omega.Ratio(x)
  
  x = sort(x,na.last=NA) 
  
  if(independence == TRUE){# weights for Probability Weighted Moments (PWM) estimators
    weight = (c(0:(length(x)-1))-.5)/(length(x)-1)#estimate F(x)
  }
  if(independence == FALSE){
    f = ecdf(x)
    weight = f(x)
  }
  
  #compute the first PWMs, see : Diebolt et al., 2008 and Guillou et al., 2009
  a0 = mean(x)
  a1 = mean(weight*x)
  a2 = mean(weight^2*x)

  omega  = (3*a2-2*a1)/(2*a1-a0)
  return(omega)
}
Omega.vec<-function(X){
  # DESCRIPTION
  # Compute the scale-invariant ratio omega for a matrix of precipitation
  # ARGUMENT
  # X         matrix of positive precipitation, each row corresponds to a day, each column to a site/gridpoint
  # VALUE 
  # Omega     vector containing omega estimates for each of the ncol(X) sites/gridpoints
  # EXAMPLE
  # n = 1000
  # x1 = rgev(n, scale = 1, shape = .1); x2 = rgev(n, scale = 4, shape = .1);  x3 = rgev(n, scale = 4, shape = .4)
  # M = cbind(x1,x2,x3); omega.vect = Omega.vec(M)
  # partition = pam(x = omega.vect, k = 2, metric = "manhattan", pamonce = 5, cluster.only = TRUE)
  
  Omega = apply(X, 2,  Omega.Ratio)
  return(Omega)
}
#######################
## FIT FUNCTIONS: EGPD (Naveau et al., 2016) with flexibility function G(u) = u^kappa
######################
## PRELIMINARY FUNCTIONS
IB <- function(x,y,a,b){
  #incomplete beta as defined in appendix of Naveau et al., 2016
  z = Ibeta(y,a,b)-Ibeta(x,a,b)
  return(z)
}
H <- function(x,sigma,xi){
  # Cumulative distribution function of GPD(sigma, xi)
  z = pgpd(x,scale=sigma,shape=xi)
  return(z)
}
pG <- function(kappa,u){
  # Probability density function of the flexibility function G
  y = u^kappa
  return(y)
}
Fbar <- function(x,kappa,sigma,xi){
  # Survival EGPD(kappa, sigma, xi) function
  z = 1-pG(kappa,H(x,sigma,xi))
  return(z)
}
## (SEMI) REGIONAL FIT (NON PARALLELIZED)
fitEGPDk <- function(M,method="pwm",cens_thres = c(0,Inf),
                     sites = "all", thres = 0, precision = .0001,
                     round = 0.1, loop.max = 10){
  # DESCRIPTION
  # fit local, semi-regional and regional versions of EGPD for chosen sites in a same cluster, see Naveau et al. (2016) and Le Gall et al. (2021)
  # ARGUMENTS
  # M               matrix of positive precipitation for sites in a same homogeneous region 
  #                  Each row corresponds to a day, each column to a site
  # method          "pwm" or "mle", method used to estimate at-site parameters
  # censor_thres    bounds between which consider reliable data to estimate parameters
  # sites           indices of gridpoints/stations where local parameters are computed (e.g. medoid, min/max silhouette)
  # thres           if not null, provides conditional parameters
  # precision       the optimization scheme stops when the parameters change is lower than precision
  #   or when the number of iteration is larger than 
  # loop.max 
  
  # VALUE
  # List of list of parameters (Theta_reg, Theta_semireg,Theta_0)
  # Theta_reg       list of regional parameters of EGPD(kappa, sigma, xi). 
  #                 Regional parameter kappa and xi.reg are common between points of a same homogeneous region.
  #                 The scale parameter sigma is site-specific.
 
  # Theta_semi_reg  list of semi-regional parameters of EGPD(kappa, sigma, xi)
  #                 The shape parameter xi.reg is the only common parameter.
  #                 The flexibility and scale parameters kappa, sigma are site specific
  
  # Theta_0         list of at-site parameters: kappa, sigma and xi are all site-specific. Values of parameter used in initialization step
  ###############################################
  #if only one temporal serie put in a matrix
  
  if(is.null(dim(M))){
    M = as.matrix(M,ncol=1)}
  nday = nrow(M)
  nstat = ncol(M)
  if(is.character(sites)){if(sites=="all"){sites = 1:nstat}else{stop("undefined option for choice of sites")} }
  nb_stat = length(sites)
  #parameters for all sites
  SigmaAll_names=KappaAll_names=rep(NA,nstat)
  CensoredMean = moyN  = rep(NA,nb_stat)
  #only days with data
  obs_rr = matrix(nrow=0,ncol=0)
  ThetaR = ThetaSR = list()
  Theta_init = list()
  kappaAll_init = sigmaAll_init = xisiteAll_init= rep(NA,nstat)
  for (station in 1:nstat){
    SigmaAll_names[station] = paste("sigma",station,sep="_")
    KappaAll_names[station] = paste("kappa",station,sep="_")
    
    y = M[,station]
    y = y[!is.na(y)]
    y = y[y>thres]-thres#only positive precipitation when thres = 0
    moy = mean(y,na.rm=TRUE)
    obs_rr = rbind.fill.matrix(obs_rr,t(as.matrix(y)))#fill the matrix with NA if needed
  }
  M = t(obs_rr)
  #initialization for optimization
  ParInit = c(0.5,0.5,0.2)
  for (station in 1:nstat) {
    y = na.omit(M[,station])
    fit_init = fit.extgp(y,model=1,method=method,init =ParInit, censoring = cens_thres, rounded = round,plots = FALSE)
    kappaAll_init[station] = fit_init$fit$pwm[1]
    sigmaAll_init[station] = fit_init$fit$pwm[2]
    xisiteAll_init[station] = fit_init$fit$pwm[3]
  }#end loop on station for initialization
  xi_init = mean(xisiteAll_init); kappa_init = mean(kappaAll_init)
  kappasite_init = kappaAll_init[sites]; sigma_init = sigmaAll_init[sites]; xisite_init = xisiteAll_init[sites]
  Theta_init = list("kappa"=rep(kappa_init,nb_stat),"sigma"=sigma_init,"xi.reg"=rep(xi_init,nb_stat))
  
  ThetaR = ThetaSR = Theta_init
  Theta_0 = list("kappa.site"=kappasite_init,"sigma"=sigma_init,"xi.site"=xisite_init)
  loop <-0; increment <-precision +.1 
  log.init<-0
  XiR.old <- unique(ThetaR$xi.reg); XiSR.old <- unique(ThetaSR$xi.reg)
  
  KappaR.old <- unique(ThetaR$kappa)#regional kappa
  u = cens_thres[1]
  for (i in 1:nb_stat){
    site = sites[i]; y = na.omit(M[,site])
    CensoredMean[i] <- mean(y[y>u],na.rm=TRUE)
  }#at-site censored mean 
  #loop for (semi)regional fit
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

  return(list("Theta_reg" = ThetaR,"Theta_semireg" = ThetaSR,"Theta_0" = Theta_0))
}

## (SEMi) REGIONAL FIT (PARALLELIZED)

fitEGPDk_paral <- function(M,ncores,method="pwm",cens_thres=c(0,Inf),
                                sites = "all", thres = 0, precision = .0001,
                                round = 0.1, loop.max = 10){
  # DESCRIPTION 
  # fit local, semi-regional and regional versions of EGPD for chosen sites in a same cluster, see Naveau et al. (2016) and Le Gall et al. (2021)
  # ARGUMENTS
  # M               matrix of positive precipitation for sites in a same homogeneous region 
  #                  Each row corresponds to a day, each column to a site
  # ncores          number of cores involved in the fitting
  # method          "pwm" or "mle", method used to estimate at-site parameters
  # censor_thres    bounds between which consider reliable data to estimate parameters
  # sites           indices of gridpoints/stations where local parameters are computed (e.g. medoid, min/max silhouette)
  # thres           if not null, provides conditional parameters
  # precision       the optimization scheme stops when the parameters change is lower than precision
  #   or when the number of iteration is larger than 
  # loop.max 
  
  # VALUE
  # List of list of parameters (Theta_reg, Theta_semireg,Theta_0)
  # Theta_reg       list of regional parameters of EGPD(kappa, sigma, xi). 
  #                 Regional parameter kappa and xi.reg are common between points of a same homogeneous region.
  #                 The scale parameter sigma is site-specific.
  
  # Theta_semi_reg  list of semi-regional parameters of EGPD(kappa, sigma, xi)
  #                 The shape parameter xi.reg is the only common parameter.
  #                 The flexibility and scale parameters kappa, sigma are site specific
  
  # Theta_0         list of at-site parameters: kappa, sigma and xi are all site-specific. Values of parameter used in initialization step
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
  #only days with data
  obs_rr = matrix(nrow=0,ncol=0)
  ThetaR = ThetaSR = list()
  Theta_init = list()
  kappaAll_init = sigmaAll_init = xisiteAll_init= rep(NA,nstat)
  for (station in 1:nstat){
    SigmaAll_names[station] = paste("sigma",station,sep="_")
    KappaAll_names[station] = paste("kappa",station,sep="_")
    
    y = M[,station]
    y = y[!is.na(y)]
    y = y[y>thres]-thres#only positive precipitation when thres =0
    moy = mean(y,na.rm=TRUE)
    obs_rr = rbind.fill.matrix(obs_rr,t(as.matrix(y)))#fill the matrix with NA if needed
  }
  M = t(obs_rr)
  #initialization for optimization
  ParInit=c(0.5,0.5,0.2)
  
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
  
  
  #loop for (semi)regional fit
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
  
  
  return(list("Theta_reg" = ThetaR,"Theta_semireg" = ThetaSR,"Theta_0" = Theta_0))
}#end fitEGPDk.boot_paral()

















