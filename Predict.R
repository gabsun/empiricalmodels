PredictEmpiricalFlux = function(trainedmodels,emodels,met_varnames,flux_varnames,
	logfilename,alldata){
	# For a given site (coming from lapply/parlapply), uses trained empirical models for all variables.
	writelog(paste0(' Predicting for site ',trainedmodels$sitename,
		':            (time ',proc.time()[3],'s)'),filename=logfilename)
	nemodels = length(emodels) # number of requested empirical model types
	emod_predictions = list() # empirical model predictions
	# Fetch empirical model function names:
	emod_functions = mget(paste0('predict_',emodels),mode='function',inherits=TRUE,envir = globalenv())
	# Call empirical model functions - each call get model for all flux variable predictions:
	for(em in 1:nemodels){
		emod_predictions[[em]] = emod_functions[[em]](trainedmodels[[em]],trainedmodels$siteindex,
			met_varnames,flux_varnames,logfilename,alldata)
	}
	# save so we can identify site for these model predictions later
	emod_predictions$siteindex = trainedmodels$siteindex
	return(emod_predictions)
}

predict_1lin = function(trainedmodel,siteindex,met_varnames,flux_varnames,logfilename,alldata){
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

predict_2lin = function(trainedmodel,siteindex,met_varnames,flux_varnames,logfilename,alldata){
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


PredictEmpFlux = function(infile,xvarnames,yvarname,emod){
	# Predicts empirically based flux for a single site.
	library(ncdf) # load netcdf library
	nxvars = length(xvarnames) # number of independent variables
	# Check consistency of inputs to function:
	if(nxvars != length(emod$grad[1,])){
		CheckError('Number of dependent vars in call to EmpFlux inconsistent')
	}
	# Determine number of time steps in testing site data:
	fid = open.ncdf(infile)
	timing = GetTimingNcfile(fid) # in PALS package
	close.ncdf(fid)
	ntsteps = timing$tsteps
	yvar=c() # dependent variable
	xvar = matrix(NA,ntsteps,nxvars) # declare x var data matrix
	# Get test site dependent variable data:
	for(v in 1:nxvars){
		if(xvarnames[v] != 'Qair'){
			tmpx = GetFluxnetVariable(xvarnames[v],infile,'blah')
			xvar[,v] = tmpx$data
		}else{
			tmpx = GetFluxnetVariable(xvarnames[v],infile,'blah')
			tmpTair = GetFluxnetVariable('Tair',infile,'blah')
			tmpPSurf = GetFluxnetVariable('PSurf',infile,'blah')
			xvar[,v] = Spec2RelHum(tmpx$data,tmpTair$data,tmpPSurf$data)
		}
	}
	if(emod$type=='kmeans'){
		cat('Sorting testing data into existing clusters... \n')
		nclst = length(emod$xclst$withinss)
		xdistall = array(NA,dim=c(ntsteps,nxvars,nclst))
		xdist = matrix(0,ntsteps,nclst)
		# Find distance of each time step to existing clusters:
		for(c in 1:nclst){
			for(v in 1:nxvars){
				xdistall[,v,c] = ((xvar[,v]-emod$xmn[v])/emod$xsd[v] - emod$xclst$centers[c,v])^2
				xdist[,c] = xdist[,c] + xdistall[,v,c]
			}
		}
		# For each time step, minvals is its distance to the cluster it belongs to:
		minvals = apply(xdist,1,function(x) min(x))
		cat('Constructing empirical flux estimate... \n')
		empflux = c()
		empflux[1:ntsteps] = 0
		# Construct empirically based flux timeseries using saved regression coefficients:
		for(c in 1:nclst){
			for(v in 1:nxvars){
				 empflux[(xdist[,c] == minvals)] = empflux[(xdist[,c] == minvals)] +
				 	xvar[(xdist[,c] == minvals),v]*emod$grad[c,v]
			}
			empflux[(xdist[,c] == minvals)] = empflux[(xdist[,c] == minvals)] + emod$int[c]
		}
	}else if(emod$type=='mlr'){
		empflux = c()
		empflux = 0 # initialise
		for(v in 1:nxvars){
			empflux = empflux + emod$grad[1,v]*xvar[,v]
		}
		empflux = empflux + emod$int[1]
	}
	return(empflux)
}
