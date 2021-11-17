##############################################################
##### Spatial clustering on ERA5 precipitation data (EU) #####
#####               PAMonce on Omega                     #####
#####    Authors: Pauline Rivoire, Philomene Le Gall     #####
#####                 03.06.2021                         #####
##############################################################

# To be adapted:
mycodes_directory="" #folder with the Rcodes for this method
mydata_directory="" #folder with precip data (ncdf files here)
output_directory="" #folder to store the outputs


message("/!\ Next line will clean the environment (of the R session only, unfortunately :/)")
rm(list=ls(all=TRUE))

# Set seed for reproductible work
set.seed(2021)


# Required packages -------------------------------------------------------

library(ncdf4) #open nectdf files
library(lubridate) #manage date format
library(fpc) #pamk, silhouette criterion
library(cluster) # PAM
library(foreach);library(iterators);library(parallel);library(doParallel) # parallelization


# Functions for clustering ------------------------------------------------

source(paste0(mycodes_directory,"useful_functions_clustering.R"))



# Settings for parallelization --------------------------------------------

nb_cores <- floor(detectCores()/3) + 1 ; registerDoParallel(cores=nb_cores)



# Open the precip netCDF file ---------------------------------------------

## Dimensions here=lon*lat*time
nc_precip_ERA5 = nc_open(paste0(mydata_directory,"era5_daily_precip_1979-01-02_2018-12-31_Europe_land.nc"))
lon_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lon") ; num_lon <- length(lon_ERA5) #Info about longitudes
lat_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lat") ; num_lat <- length(lat_ERA5) #Info about latitudes
times_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "time")
date_ERA5 <- format(as.POSIXlt(times_ERA5*60*60, origin="1979-01-02 12:30:00"), format = "%Y-%m-%d") #Convert to date format
precip_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "TP") #/!\ can be long depending on the dataset (~2min for Era5 land precipitation)
nc_close(nc_precip_ERA5)

nmin_wetdays <- 300 #minimum number of wet days (>0mm) at the grid point to not be discarded (otherwise not enough data for fitting)


# Define the season -------------------------------------------------------

seas <- "JJA" #seas should be in c("SON", "DJF", "MAM", "JJA")
print(seas)

precip_positive_list <- list() #build a list with precipitation for each gridpoints
#keep only gridpoints without NA. We preprocessed ERA5 data to have only preicp over Europe's land
GP_kept <- !(apply(X=precip_ERA5[,,dim(precip_ERA5)[3]], MARGIN = c(1,2), FUN = is.na))
ref_in_list <- matrix(data = NA, nrow = nrow(GP_kept), ncol = ncol(GP_kept)) #keep track of the position of the gridpoint in the list

counting <- 0 #count to have the reference in the list
for (LON in 1:num_lon) {
  for (LAT in 1:num_lat) {

    if(GP_kept[LON,LAT]){ #if the gridpoint is a land point
      positive_precip <- extr_seas_pos_precip(season = seas, precip_timeserie = precip_ERA5[LON,LAT,],
                                              precip_date = date_ERA5, thshld = 0) #get the seasonal positive precipitation
      if (length(which(positive_precip>=1))>nmin_wetdays){#if the number of seasonal wet days is large enough
        counting <- counting + 1
        precip_positive_list[[counting]] <- positive_precip
        ref_in_list[LON,LAT] <- counting
      } else {
        GP_kept[LON,LAT] <- FALSE
      }#end if long enough
    }#end if kept

  }#end for LAT
}#end for LON

nb_GP_kept <- sum(GP_kept) #number of grid points kept

Omega_list <- foreach(GP=1:nb_GP_kept) %dopar% { #parallelized computation of omega for all the selected gridpoints

  RESU <- xi.Ratio(precip_positive_list[[GP]]) #function from the file "useful_functions_clustering.R"

  return(RESU)
} #end foreach GP

Omega_info <- list(Omega_vec=unlist(Omega_list), ref_in_list=ref_in_list)

save(Omega_info, file=paste0(output_directory, "Omega_vector_", seas, ".Rdata"))

load(file=paste0(output_directory, "Omega_vector_", seas, ".Rdata"))

nb_GP_kept <- length(Omega_info$Omega_vec)

dist_matrix <- dist(Omega_info$Omega_vec, method = "manhattan") #compute the distance matrix for the clustering procedure

nb_clusters = c(2:10) #Number of clusters for the partition (several options, to compare the partitions and select the best one)
list_partition <- list() #List of the output of the clustering for all numbers of clusters

#Loop that can be long, depending on the machine
list_partition <- foreach(ind=1:length(nb_clusters)) %dopar% { #paralellized loop for the clustering algorithm (one step=one number of cluster)
  NbClusters = nb_clusters[ind]
  #clustering function. pamonce = 5 optimizes the calculation and reduce the computation time
  partition_pam = pam(dist_matrix, k = NbClusters, pamonce = 5)
  return(partition_pam)
}#end foreach ind

save(list_partition, file = paste0(output_directory, "list_PAMonce_omega_10clus_", seas, ".Rdata"))



# Load partition data -----------------------------------------------------

load(file = paste0("/scratch3/pauline/Omega_Spatial_clustering/Clustering_procedure/Output_clustering/list_PAMonce_omega_10clus_", seas, ".Rdata"))
sil_crit = rep(0, length(nb_clusters))
for (ind in 1:length(nb_clusters)) {
  sil_crit[ind] <- list_partition[[ind]]$silinfo$avg.width
  print(paste("nb cluster", (ind+1), "sil", list_partition[[ind]]$silinfo$avg.width))
}#end for ind



# Store Omega and partitions in a ncdf file ---------------------------------

Omega_matrix <- matrix(data = NA, ncol = length(lat_ERA5), nrow = length(lon_ERA5))
cluster_matrix_2 <- matrix(data = NA, ncol = length(lat_ERA5), nrow = length(lon_ERA5)) #partition for 2 clusters,
cluster_matrix_3 <- matrix(data = NA, ncol = length(lat_ERA5), nrow = length(lon_ERA5)) #partition for 3 clusters,
cluster_matrix_4 <- matrix(data = NA, ncol = length(lat_ERA5), nrow = length(lon_ERA5)) #etc...
cluster_matrix_5 <- matrix(data = NA, ncol = length(lat_ERA5), nrow = length(lon_ERA5))
cluster_matrix_6 <- matrix(data = NA, ncol = length(lat_ERA5), nrow = length(lon_ERA5))
cluster_matrix_7 <- matrix(data = NA, ncol = length(lat_ERA5), nrow = length(lon_ERA5))
cluster_matrix_8 <- matrix(data = NA, ncol = length(lat_ERA5), nrow = length(lon_ERA5))
cluster_matrix_9 <- matrix(data = NA, ncol = length(lat_ERA5), nrow = length(lon_ERA5))
cluster_matrix_10 <- matrix(data = NA, ncol = length(lat_ERA5), nrow = length(lon_ERA5))
keep_coord_list <- matrix(data = NA, nrow = nb_GP_kept, ncol = 2) #keep track of the position in the list
colnames(keep_coord_list) <- c("LON ref", "LAT ref")

for (LON in 1:num_lon) {
  for (LAT in 1:num_lat) {
    ref <- Omega_info$ref_in_list[LON,LAT]
    
    Omega_matrix[LON,LAT] <- Omega_info$Omega_vec[ref]
    cluster_matrix_2[LON,LAT] <- list_partition[[1]]$clustering[ref]
    cluster_matrix_3[LON,LAT] <- list_partition[[2]]$clustering[ref] ; cluster_matrix_4[LON,LAT] <- list_partition[[3]]$clustering[ref]
    cluster_matrix_5[LON,LAT] <- list_partition[[4]]$clustering[ref] ; cluster_matrix_6[LON,LAT] <- list_partition[[5]]$clustering[ref]
    cluster_matrix_7[LON,LAT] <- list_partition[[6]]$clustering[ref] ; cluster_matrix_8[LON,LAT] <- list_partition[[7]]$clustering[ref]
    cluster_matrix_9[LON,LAT] <- list_partition[[8]]$clustering[ref] ; cluster_matrix_10[LON,LAT] <- list_partition[[9]]$clustering[ref]
    
    keep_coord_list[ref, "LON ref"] <- LON
    keep_coord_list[ref, "LAT ref"] <- LAT
    
  }#end for LAT
}#end for LON

save(keep_coord_list, file = paste0(output_directory,"coord_in_list_PAMonce_",seas,".Rdata"))

# create and write the netCDF file -- ncdf4 version
ncfname <- paste0(output_directory, "Omega_clusterpamonce_10_ERA5_EU_",seas,".nc")
# define dimensions
londim <- ncdim_def(name = "longitude", units = "degrees_east", vals = lon_ERA5)
latdim <- ncdim_def(name = "latitude", units = "degrees_north", vals = lat_ERA5)

# define variables
fillvalue <- 1e32

Omega_value <- ncvar_def(name = "Omega_value", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)

cluster_PAM_2 <- ncvar_def(name = "cluster_num_2", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)
cluster_PAM_3 <- ncvar_def(name = "cluster_num_3", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)
cluster_PAM_4 <- ncvar_def(name = "cluster_num_4", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)
cluster_PAM_5 <- ncvar_def(name = "cluster_num_5", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)
cluster_PAM_6 <- ncvar_def(name = "cluster_num_6", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)
cluster_PAM_7 <- ncvar_def(name = "cluster_num_7", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)
cluster_PAM_8 <- ncvar_def(name = "cluster_num_8", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)
cluster_PAM_9 <- ncvar_def(name = "cluster_num_9", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)
cluster_PAM_10 <- ncvar_def(name = "cluster_num_10", units = "no_unit", dim = list(londim, latdim), missval = fillvalue)

# create netCDF file and put arrays
ncout <- nc_create(ncfname,
                   list(Omega_value, cluster_PAM_2, cluster_PAM_3, cluster_PAM_4, cluster_PAM_5,
                        cluster_PAM_6, cluster_PAM_7, cluster_PAM_8, cluster_PAM_9, cluster_PAM_10),
                   force_v4=TRUE)

# put variables
ncvar_put(ncout, Omega_value, Omega_matrix)
ncvar_put(ncout, cluster_PAM_2, cluster_matrix_2)
ncvar_put(ncout, cluster_PAM_3, cluster_matrix_3)
ncvar_put(ncout, cluster_PAM_4, cluster_matrix_4)
ncvar_put(ncout, cluster_PAM_5, cluster_matrix_5)
ncvar_put(ncout, cluster_PAM_6, cluster_matrix_6)
ncvar_put(ncout, cluster_PAM_7, cluster_matrix_7)
ncvar_put(ncout, cluster_PAM_8, cluster_matrix_8)
ncvar_put(ncout, cluster_PAM_9, cluster_matrix_9)
ncvar_put(ncout, cluster_PAM_10, cluster_matrix_10)

nc_close(ncout)
