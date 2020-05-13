# Driver.R
#
# This driver is used, over a collection of tower sites, to (a) train empirical
# models used to benchmark LM simulations using flux tower data (b) predict surface
# fluxes with these trained models, and (c) write predictions to file in ALMA format
# netcdf, as though a LM.
#
# Gab Abramowitz CCRC/CLEX, UNSW 2020 (gabsun at gmail dot com)
#
# Workflow:
#
# 1. Load and partition site data
# 2. Train with training set, save empirical model parameters
# 3. Use those parameters to predict with prediction subset
# 4. Write predictions to netcdf file(s)

library(pals)
library(pryr)
library(ncdf4)
library(kohonen)
source('LoadSortData.R')
source('Train.R')
source('Predict.R')
source('WriteFiles.R')

#metfile_dir = '/Users/gab/data/flux_tower/FluxnetLSM/PLUMBER2_sample/met_inputs/'
#fluxfile_dir = '/Users/gab/data/flux_tower/FluxnetLSM/PLUMBER2_sample/flux_files/'
metfile_dir = '/Users/gab/data/flux_tower/FluxnetLSM/PLUMBER2/Nc_files/Met/'
fluxfile_dir = '/Users/gab/data/flux_tower/FluxnetLSM/PLUMBER2/Nc_files/Flux/'
data_save_dir = '/Users/gab/results/PLUMBER/empirical_models/PLUMBER2/'
logfilename = 'logEmpmodel.log'
met_varnames = c('SWdown','Tair','RH')
flux_varnames = c('Qle','Qh','NEE') #'Qle','Qh','NEE'
emodels = c('3som100') #,'2lin','3som25',3km27)
sitenumbers = c(22:22) # which sites should we build models for (in above dirs)

system(paste('rm',logfilename)) # remove any old logfile
openlog(filename=logfilename) # open log file to detail proceedings

# Load all site data:
data = LoadData(metfile_dir,fluxfile_dir,met_varnames,flux_varnames,
	data_save_dir,logfilename)

for(s in sitenumbers[1] : sitenumbers[length(sitenumbers)] ){
	# Train empirical models, parallel applied over each site:
	trainedmodel = TrainEmpiricalModels(data[[s]],emodels,met_varnames,
		flux_varnames,logfilename,data)
	# Predict using trained empirical models, parallel applied over each site:
	prediction = PredictEmpiricalFlux(trainedmodel,emodels,met_varnames,
		flux_varnames,logfilename,data)
	# create temporary list:
	trainedmodels = list()
	trainedmodels[[s]] = trainedmodel
	# Write predictions and saved models to file:
	write_to_file = write_emp_predictions(prediction,emodels,met_varnames,
		flux_varnames,logfilename,data,data_save_dir,trainedmodels)
}

closelog(filename=logfilename)
