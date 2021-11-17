###############################################################
##### Convert output of the regional fittings to netcdf   #####
#####            Pauline Rivoire, 09.04.2021              #####
###############################################################

rm(list=ls())


for (seas in c("SON", "DJF", "MAM", "JJA")) {
  if(seas=="JJA"){
    nclusters = 5
  } else {
    nclusters=3
  }
  
  load(file = paste0(output_directory,"res_REG_fit_PAM_", nclusters,"clus_allGP_",seas,"_1third.Rdata"))
  load(file = paste0(output_directory,"keep_coord_REG_fit_PAM_", nclusters,"clus_allGP_",seas,"_1third.Rdata"))
  
  # get the coordinates
  nc_precip_ERA5 = nc_open("/scratch3/pauline/EGPD_on_ERA5/Data/era5_daily_precip_1979-01-02_2018-12-31_Europe_land.nc")
  lon_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lon") ; lat_ERA5 <- ncvar_get(nc = nc_precip_ERA5, varid = "lat")
  nc_close(nc_precip_ERA5)
  
  xi_mat <- matrix(data = NA, nrow = nrow(RL_reg), ncol = ncol(RL_reg))
  sigma_mat = kappa_mat = xi_mat
  
  for (clus in 1:nclusters) {
    for (GP in 1:length(List_keep_coord[[clus]]$keep_lon)) {
      #store the shape parameter
      xi_mat[List_keep_coord[[clus]]$keep_lon[GP],List_keep_coord[[clus]]$keep_lat[GP]] <- List_reg_fit[[clus]]$Theta$xi.reg[GP]
      #store the flexibility parameter
      kappa_mat[List_keep_coord[[clus]]$keep_lon[GP],List_keep_coord[[clus]]$keep_lat[GP]] <- List_reg_fit[[clus]]$Theta$kappa[GP]
      #store the scale parameter
      sigma_mat[List_keep_coord[[clus]]$keep_lon[GP],List_keep_coord[[clus]]$keep_lat[GP]] <- List_reg_fit[[clus]]$Theta$sigma[GP]
    }#end for GP
  }#end for clus
  
  
  # Create nc file ----------------------------------------------------------

  ncfname <- paste0(output_directory,"parameters_REG_fit_PAM_",nclusters,"clus_allGP_",seas,"_1third.nc")
  
  
  # define dimensions
  londim <- ncdim_def(name = "longitude", units = "degrees_east", vals = lon_ERA5)
  latdim <- ncdim_def(name = "latitude", units = "degrees_north", vals = lat_ERA5)
  
  # define variables
  fillvalue <- 1e32
  
  xi_var <- ncvar_def(name = "xi", units = "none", dim = list(londim, latdim), missval = fillvalue)
  sigma_var <- ncvar_def(name = "sigma", units = "none", dim = list(londim, latdim), missval = fillvalue)
  kappa_var <- ncvar_def(name = "kappa", units = "none", dim = list(londim, latdim), missval = fillvalue)
  
  ncout <- nc_create(ncfname, list(xi_var,sigma_var,kappa_var), force_v4=TRUE)
  
  ncvar_put(ncout, xi_var, xi_mat)
  ncvar_put(ncout, kappa_var, kappa_mat)
  ncvar_put(ncout, sigma_var, sigma_mat)
  
  nc_close(ncout)
  
  
  
}# end for seas
