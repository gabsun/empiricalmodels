# Predict.R
#
# Functions for using existing empirical model parameters to make predictions
# based on met data.
#
# Gab Abramowitz CCRC/CLEX, UNSW 2020 (gabsun at gmail dot com)

PredictEmpiricalFlux = function(trainedmodels,emodels,met_varnames,flux_varnames,
	logfilename,alldata){
	# For a given site (coming from lapply/parlapply), uses trained empirical models for all variables.
	writelog(paste0(' Predicting for site ',alldata[[trainedmodels$siteindex]]$name,
		':            (time ',proc.time()[3],'s)'),filename=logfilename)
	nemodels = length(emodels) # number of requested empirical model types
	emod_predictions = list() # empirical model predictions
	# Decide which empirical model functions need to be called:
	clustermods = c('som','km') # text flags that we're going to call the cluster function
	emod_simple = emodels # init. names of empirical models with all cluster models reduced to 'cluster'
	for(cl in 1:length(clustermods)){
		emod_simple = ifelse(grepl(pattern=clustermods[cl],emod_simple),'clustermod',emod_simple)
	}
	# Fetch empirical model function names:
	emod_functions = mget(paste0('predict_',emod_simple),mode='function',inherits=TRUE,envir = globalenv())
	# Call empirical model functions - each call get model for all flux variable predictions:
	for(em in 1:nemodels){
		emod_predictions[[em]] = emod_functions[[em]](trainedmodels[[em]],trainedmodels$siteindex,
			emodels[em],met_varnames,flux_varnames,logfilename,alldata)
	}
	# save so we can identify site for these model predictions later
	emod_predictions$siteindex = trainedmodels$siteindex
	return(emod_predictions)
}

predict_1lin = function(trainedmodel,siteindex,modelname,met_varnames,flux_varnames,logfilename,alldata){
	# Linear model predicting flux from SWdown only. Uses existing linear model parameters for
	# all requested fluxes to be predicted.
	writelog(paste0('  predicting with 1lin...  (',smem(),'; ',proc.time()[3],'s)'),filename=logfilename)
	nfluxes = length(flux_varnames) # number of flux predictions we'll require
	nmetvars = length(met_varnames) # number of met variables in total load request
	tsteps = alldata[[ siteindex ]]$alltsteps # number of timesteps at prediction site

	# Load data required for prediction:
	SWdownIndex = c(1:nmetvars)[met_varnames=='SWdown'] # logical mask to pull SWdown out of all met data
	SWdown = alldata[[ siteindex ]]$invars[[ SWdownIndex ]]
	predictedfluxes = matrix(NA,nfluxes,tsteps) # initialise
	for(f in 1:nfluxes){
		predictedfluxes[f,] = (trainedmodel[[f]]$gradient * SWdown) + trainedmodel[[f]]$intercept
	}
	return(predictedfluxes)
}

predict_2lin = function(trainedmodel,siteindex,modelname,met_varnames,flux_varnames,logfilename,alldata){
	# Linear model predicting flux from SWdown and Tair. Uses linear model parameters for
	# all requested fluxes to be predicted.
	writelog(paste0('  predicting with 2lin...  (',smem(),'; ',proc.time()[3],'s)'),filename=logfilename)
	nfluxes = length(flux_varnames) # number of flux predictions we'll require
	nmetvars = length(met_varnames) # number of met variables in total load request
	tsteps = alldata[[ siteindex ]]$alltsteps # number of timesteps at prediction site
	SWdownIndex = c(1:nmetvars)[met_varnames=='SWdown'] # logical mask to pull SWdown out of all met data
	TairIndex = c(1:nmetvars)[met_varnames=='Tair'] # logical mask to pull Tair out of all met data

	# Load data required for prediction:
	SWdown = alldata[[ siteindex ]]$invars[[ SWdownIndex ]]
	Tair = alldata[[ siteindex ]]$invars[[ TairIndex ]]
	predictedfluxes = matrix(NA,nfluxes,tsteps) # initialise
	for(f in 1:nfluxes){
		predictedfluxes[f,] = (trainedmodel[[f]]$gradients[1] * SWdown) +
			(trainedmodel[[f]]$gradients[2] * Tair) + trainedmodel[[f]]$intercept
	}
	return(predictedfluxes)
}

predict_clustermod = function(trainedmodel,siteindex,modelname,met_varnames,
	flux_varnames,logfilename,alldata){
	# Linear model predicting flux from SWdown and Tair. Uses linear model
	# parameters for all requested fluxes to be predicted.
	writelog(paste0('  predicting with ',modelname,
		'...  (',smem(),'; ',proc.time()[3],'s)'),filename=logfilename)
	nfluxes = length(flux_varnames) # number of flux predictions we'll require
	nmetvars = length(met_varnames) # number of met variables in total load request
	tsteps = alldata[[ siteindex ]]$alltsteps # number of timesteps at prediction site

	# Determine cluster type, number of clusters, and number of input variables:
	if(grepl('som',modelname)){
		cluster_type='som'
		loc = regexpr('som',modelname)
		nclst = as.integer(substr(modelname,(loc[1]+3),nchar(modelname)))
		nxvars = as.integer(substr(modelname,1,(loc[1]-1)))
	}else if(grepl('km',modelname)){
		cluster_type='kmeans'
		loc = regexpr('km',modelname)
		nclst = as.integer(substr(modelname,(loc[1]+2),nchar(modelname)))
		nxvars = as.integer(substr(modelname,1,(loc[1]-1)))
	}else if(grepl('cl',modelname)){
		cluster_type='clara'
		loc = regexpr('cl',modelname)
		nclst = as.integer(substr(modelname,(loc[1]+2),nchar(modelname)))
		nxvars = as.integer(substr(modelname,1,(loc[1]-1)))
	}
	# For now (this could be changed), assume that '3' input variables implies the
	# first three in the met_varnmes list of loaded variables:
	xvarnames = met_varnames[1:nxvars]

	# Load data required for prediction:
	metvars = matrix(NA,tsteps,nxvars) # initialise
	predictedfluxes = matrix(0,nfluxes,tsteps) # initialise

	for(xv in 1:nxvars){ # for each met variable:
		# First determine the index of this variable in the loaded met data:
		vIndex = c(1:nmetvars)[met_varnames==xvarnames[xv]]
		# Then load data :
		metvars[,xv] = alldata[[ siteindex ]]$invars[[ vIndex ]]
	}

	for(f in 1:nfluxes){
		nclst = trainedmodel[[f]]$nclst # SOM may change this number if not a square
		cluster_type = trainedmodel[[f]]$cluster_type
		if(cluster_type == 'kmeans'){
			xdistall = array(NA,dim=c(tsteps,nxvars,nclst)) # distance to each cluster centre for each datum
			xdist = matrix(0,tsteps,nclst) # distance to cluster centre as sum of squares over variables
			# Find distance of each time step to existing clusters:
			for(cl in 1:nclst){
				for(v in 1:nxvars){
					xdistall[,v,cl] = ( ( (metvars[,v]-trainedmodel[[f]]$metvar_means[v]) /
						trainedmodel[[f]]$metvar_sd[v] ) - trainedmodel[[f]]$clusters$centers[cl,v] )^2
					# Find cumulative squared distance over all predictor variables:
					xdist[,cl] = xdist[,cl] + xdistall[,v,cl]
				}
			}
			# For each time step, mindex is the number of the cluster it belongs to:
			mindex = apply(xdist,MARGIN=1,FUN=min_index)
		}else if(cluster_type == 'som'){
			# Scale met data first, so that mapping to clusters makes sense:
			wmetvars = matrix(NA,tsteps,nxvars) # initialise
			for(v in 1:nxvars){
				wmetvars[,v] = (metvars[,v]-trainedmodel[[f]]$metvar_means[v]) / trainedmodel[[f]]$metvar_sd[v]
			}
			# Then find which cluster each timestep belongs to:
			sortdata = map(x=trainedmodel[[f]]$clusters,newdata=wmetvars)
			mindex = sortdata$unit.classif
		}
		empflux = c()
		empflux[1:tsteps] = 0
		# Construct empirically based flux timeseries using saved regression coefficients:
		for(cl in 1:nclst){
			for(v in 1:nxvars){
				# Note that (mindex==c) is a logical mask for the timesteps belonging to cluster 'c'
				predictedfluxes[f,(mindex==cl)] = predictedfluxes[f,(mindex==cl)] +
					metvars[(mindex==cl),v] * trainedmodel[[f]]$gradients[cl,v]
			}
			predictedfluxes[f,(mindex==cl)] = predictedfluxes[f,(mindex==cl)] + trainedmodel[[f]]$intercept[cl]
		}
	}

	return(predictedfluxes)
}
