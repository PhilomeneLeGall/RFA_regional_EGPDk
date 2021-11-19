This folder contains the codes to perform the regional fitting, after the clustering procedure (see Rivoire, Le Gall et al. 2021, submitted).
The fitted function are extended generalized Pareto distributions [EGPD(kappa, sigma, xi) introduced in Naveau et al. (2016)]

regional_fitting.R and semi_regional_fitting.R contain the code to perform the regional and semiregional fits. They are parallelized, and rely on functions from parallelized_fitting_functions.R.

Convert_to_ndcf.R convert the output of the fitting to a netcdf file.
