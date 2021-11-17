######################################################################
##### Regional inference of EGPD on ERA5 precipitation data (EU) #####
#####                      All Gridpoints                        #####
#####       Pauline Rivoire, Philomene Le Gall 09.04.2021        #####
######################################################################
rm(list=ls(all=TRUE))
# To be adapted:
mycodes_directory="" #folder with the Rcodes for this method
mydata_directory="" #folder with precip data (ncdf files here)
output_directory="" #folder to store the outputs


library(ncdf4)


# Select the season -------------------------------------------------------
seas <- "SON"
print(seas)

if (seas=="JJA") {
  nclust <- 5; var_name = "cluster_num_5"
} else {
    nclust <- 3; var_name="cluster_num_3"
}#The optimal number of clusters get be obtained with the mean silhouette coefficient, or other metrics, or visually


# get the required functions ----------------------------------------------
source(paste0(mycodes_directory,"parallelized_fitting_functions.R"))


# get the precipitation data ----------------------------------------------
nc_precip_ERA5 = nc_open(paste0(mydata_directory,"era5_daily_precip_1979-01-02_2018-12-31_Europe_land.nc"))
lon_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lon") ; num_lon <- length(lon_ERA5)
lat_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lat") ; num_lat <- length(lat_ERA5)
times_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "time")
date_ERA5 <- format(as.POSIXlt(times_ERA5*60*60, origin="1979-01-02 12:30:00"), format = "%Y-%m-%d")
precip_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "TP") # Dimensions=lon*lat*time
nc_close(nc_precip_ERA5)


# get the PAM partition data ----------------------------------------------
nc_clusters <- nc_open(paste0(output_directory,"Omega_clusterpam_10_ERA5_EU_", seas, ".nc"))
optimal_partition <- ncvar_get(nc_clusters, varid = var_name) #select only the partition with optimal number of clusters
lon_clus <- ncvar_get(nc_clusters, varid = "longitude")
lat_clus <- ncvar_get(nc_clusters, varid = "latitude")
nc_close(nc_clusters)


List_reg_fit <- list() #store the fitting parameters
List_keep_coord <- list() #keep in memory the coordiantes of the gridpoint
for(clust_nb in 1:nclust){ #loop on the clusters
  M_pos_precip <- matrix(nrow=0,ncol=0) #matrix of seasonal positive precipititatiom 
  keep_lon <- numeric()
  keep_lat <- numeric()
  counting <- 0 #keep track of the position in the list
  for (LON in 1:num_lon) {
    for (LAT in 1:num_lat) {
        if(!is.na(optimal_partition[LON,LAT]) & optimal_partition[LON,LAT]==clust_nb){ #check we have data and the GP is indeed in the good cluster
          counting <- counting + 1
          keep_lon[counting] <- LON ; keep_lat[counting] <- LAT
          #extract seasonal positive precipitation
          y <- extr_seas_pos_precip(precip_timeserie = precip_ERA5[LON,LAT,], season = seas, precip_date = date_ERA5, thshld = 0)
          ind_1_third <- (1:length(y))[which(((1:length(y))%%3)==0)]#Only keep one third of precipitation time series to get rid of autocorelation
          M_pos_precip <- rbind.fill.matrix(M_pos_precip, t(as.matrix(y[ind_1_third])))
        }#end if in cluster
    }#end for LAT
  }#end for LON : took 1h
  rm(counting)
  
  List_keep_coord[[clust_nb]] <- list(keep_lon = keep_lon, keep_lat = keep_lat)
  List_reg_fit[[clust_nb]] <- fitEGPDkREGCensorIter_paral(M = t(M_pos_precip), ncores = floor(detectCores()/3 - 1), cens_thres = c(1,Inf))
}#end for clust_nb

save(List_reg_fit, file = paste0(output_directory,"res_REG_fit_PAM_",nclust,"clus_allGP_",seas,"_1third.Rdata"))
save(List_keep_coord, file = paste0(output_directory,"keep_coord_REG_fit_PAM_",nclust,"clus_allGP_",seas,"_1third.Rdata"))

