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
library(parallel)
library(pryr)
library(ncdf4)
source('LoadSortData.R')
source('Train.R')
source('Predict.R')
source('WriteFiles.R')

runparallel = TRUE # single or multicore execution: FALSE for debugging, much faster if TRUE
#metfile_dir = '/Users/gab/data/flux_tower/FluxnetLSM/PLUMBER2_sample/met_inputs/'
#fluxfile_dir = '/Users/gab/data/flux_tower/FluxnetLSM/PLUMBER2_sample/flux_files/'
metfile_dir = '/Users/gab/data/flux_tower/FluxnetLSM/PLUMBER2/Nc_files/Met/'
fluxfile_dir = '/Users/gab/data/flux_tower/FluxnetLSM/PLUMBER2/Nc_files/Flux/'
outfile_dir = '/Users/gab/results/PLUMBER/empirical_models/model_output/'
# if dirs above are changed, this need to be deleted:
tmp_data_save_file = '/Users/gab/results/PLUMBER/empirical_models/saved_variables.Rdat'
logfilename = 'logEmpiricalmodel.log'
met_varnames = c('SWdown','Tair','RH','Wind','Precip')
flux_varnames = c('Qle','Qh','NEE')
emodels = c('1lin','2lin') #,'2lin','3km27')

system(paste('rm',logfilename)) # remove any old logfile
openlog(filename=logfilename) # open log file to detail proceedings

# Load all site data:
data = LoadData(metfile_dir,fluxfile_dir,met_varnames,flux_varnames,
	tmp_data_save_file,logfilename)

if(runparallel){
	cl = makeCluster(getOption('cl.cores', 2),type='FORK') #	Create cluster
	# Train empirical models, parallel applied over each site:
	trainedmodels = parLapply(cl=cl,data,TrainEmpiricalModels,emodels,met_varnames,
		flux_varnames,logfilename,data)
	# Predict using trained empirical models, parallel applied over each site:
	predictions = parLapply(cl=cl,trainedmodels,PredictEmpiricalFlux,emodels,met_varnames,
		flux_varnames,logfilename,data)
	write_to_file = parLapply(cl=cl,predictions,write_emp_predictions,emodels,
		met_varnames,flux_varnames,logfilename,data,outfile_dir,trainedmodels)
	stopCluster(cl)
}else{
	# Train empirical models, parallel applied over each site:
	trainedmodels = lapply(data,TrainEmpiricalModels,emodels,met_varnames,
		flux_varnames,logfilename,data)
	# Predict using trained empirical models, parallel applied over each site:
	predictions = lapply(trainedmodels,PredictEmpiricalFlux,emodels,met_varnames,
		flux_varnames,logfilename,data)
	write_to_file = lapply(predictions,write_emp_predictions,emodels,met_varnames,
		flux_varnames,logfilename,data,outfile_dir,trainedmodels)
}

closelog(filename=logfilename)
