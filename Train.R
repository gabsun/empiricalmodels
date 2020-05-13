# Train.R
#
# Functions to train empirical models
#
# Gab Abramowitz CCRC/CLEX, UNSW 2020 (gabsun at gmail dot com)

TrainEmpiricalModels = function(sitedata,emodels,met_varnames,flux_varnames,logfilename,
	alldata,removeflagged=TRUE){
	# For a given site (coming from lapply/parlapply), trains all empirical models for all flux variables.
	# First, define all shorthand representations of cluster+regression emp models:
	writelog(paste0(' Training model for site ',sitedata$name,
		':            (time ',proc.time()[3],'s)'),filename=logfilename)
	nemodels = length(emodels) # number of requested empirical model types
	emod_pars = list() # empirical model parameters
	# Decide which empirical model functions need to be called:
	clustermods = c('som','km') # text flags that we're going to call the cluster function
	emod_simple = emodels # init. names of empirical models with all cluster models reduced to 'cluster'
	for(cl in 1:length(clustermods)){
		emod_simple = ifelse(grepl(pattern=clustermods[cl],emod_simple),'clustermod',emod_simple)
	}
	# Fetch empirical model function names. emod_functions will functions named, e.g. train_lin, or
	# train_2lin, or train_cluster for models involving 'km' or 'som':
	emod_functions = mget(paste0('train_',emod_simple),mode='function',inherits=TRUE,envir = globalenv())
	# Call empirical model functions - each call gets model for all flux variable predictions:
	for(em in 1:nemodels){
		emod_pars[[em]] = emod_functions[[em]](sitedata,emodels[em],met_varnames,
			flux_varnames,logfilename,alldata,removeflagged)
	}
	# save so we can identify site for these model parameters later
	emod_pars$siteindex = sitedata$index
	return(emod_pars)
}

train_1lin = function(sitedata,modelname,met_varnames,flux_varnames,logfilename,alldata,removeflagged){
	# Linear model predicting flux from SWdown only. Gets linear model parameters for
	# all requested fluxes to be predicted.
	nfluxes = length(flux_varnames) # number of flux predictions we'll require
	nmetvars = length(met_varnames) # number of met variables in total load request
	nsites = length(alldata) - 1 # number of sites to collate training data over (less current site)
	SWdownIndex = c(1:nmetvars)[met_varnames=='SWdown'] # logical mask to pull SWdown out of all met data

	# Collate data for training (this'll suck up some memory):
	SWdown=c(1:sitedata$trainset_tsteps) # initialise to size of number of training time steps
	fluxes = matrix(NA,nfluxes,sitedata$trainset_tsteps)
	if(removeflagged) fluxqc = matrix(NA,nfluxes,sitedata$trainset_tsteps)
	ctr = 1
	for(s in 1:nsites){
		tsteps = alldata[[ sitedata$train_index[s] ]]$alltsteps # timesteps in this training site
		# First get SWdown:
		SWdown[ctr:(ctr+tsteps-1)] = alldata[[ sitedata$train_index[s] ]]$invars[[ SWdownIndex ]]
		# Then get requested fluxes:
		for(f in 1:nfluxes){
			fluxes[f,ctr:(ctr+tsteps-1)] = alldata[[ sitedata$train_index[s] ]]$outvars[[ f ]]
			if(removeflagged) fluxqc[f,ctr:(ctr+tsteps-1)] = alldata[[ sitedata$train_index[s] ]]$outvarsqc[[ f ]]
		}
		ctr = ctr + tsteps
	}
	# Report to gauge memory usage of data replication:
	writelog(paste0('  1lin data written, now calculating...  (',smem(),'; ',proc.time()[3],'s)'),
		filename=logfilename)
	# Now find linear model parameters for each flux:
	fluxpars = list()
	for(f in 1:nfluxes){
		fluxpars[[f]] = list()
		# Fit linear model:
		if(removeflagged){
			tmp_pars = lm(fluxes[f,]~SWdown,subset=(fluxqc[f,]==0),na.action=na.omit)
		}else{
			tmp_pars = lm(fluxes[f,]~SWdown,na.action=na.omit)
		}
		# Save regression paramters:
		fluxpars[[f]]$intercept = tmp_pars$coefficients[1]
		fluxpars[[f]]$gradient = tmp_pars$coefficients[2]
	}
	# Add string to indicate which independent variables were used (for writing prediction file):
	fluxpars$eminputs = 'SWdown'

	return(fluxpars)
}

train_2lin = function(sitedata,modelname,met_varnames,flux_varnames,logfilename,alldata,removeflagged){
	# Linear model predicting flux from SWdown and Tair. Gets linear model parameters for
	# all requested fluxes to be predicted.
	nfluxes = length(flux_varnames) # number of flux predictions we'll require
	nmetvars = length(met_varnames) # number of met variables in total load request
	nsites = length(alldata) - 1 # number of sites to collate training data over (less current site)
	SWdownIndex = c(1:nmetvars)[met_varnames=='SWdown'] # logical mask to pull SWdown out of all met data
	TairIndex = c(1:nmetvars)[met_varnames=='Tair'] # logical mask to pull Tair out of all met data

	# Collate data for training (this'll suck up some memory):
	metvars = matrix(NA,sitedata$trainset_tsteps,2) # initialise to size of number of training time steps
	fluxes = matrix(NA,nfluxes,sitedata$trainset_tsteps)
	if(removeflagged) fluxqc = matrix(NA,nfluxes,sitedata$trainset_tsteps)
	ctr = 1
	for(s in 1:nsites){
		tsteps = alldata[[ sitedata$train_index[s] ]]$alltsteps # timesteps in this training site
		# First get SWdown, then Tair:
		metvars[ctr:(ctr+tsteps-1),1] = alldata[[ sitedata$train_index[s] ]]$invars[[ SWdownIndex ]]
		metvars[ctr:(ctr+tsteps-1),2] = alldata[[ sitedata$train_index[s] ]]$invars[[ TairIndex ]]
		# Then get requested fluxes:
		for(f in 1:nfluxes){
			fluxes[f,ctr:(ctr+tsteps-1)] = alldata[[ sitedata$train_index[s] ]]$outvars[[ f ]]
			if(removeflagged) fluxqc[f,ctr:(ctr+tsteps-1)] = alldata[[ sitedata$train_index[s] ]]$outvarsqc[[ f ]]
		}
		ctr = ctr + tsteps
	}
	# Report to gauge memory usage of data replication:
	writelog(paste0('  2lin data written, now calculating...  (',smem(),'; ',proc.time()[3],'s)'),
		filename=logfilename)
	# Now find linear model parameters for each flux:
	fluxpars = list()
	for(f in 1:nfluxes){
		fluxpars[[f]] = list()
		# Fit linear model:
		if(removeflagged){
			tmp_pars = lm(fluxes[f,]~metvars,subset=(fluxqc[f,]==0),na.action=na.omit)
		}else{
			tmp_pars = lm(fluxes[f,]~metvars,na.action=na.omit)
		}
		# Save regression paramters:
		fluxpars[[f]]$intercept = tmp_pars$coefficients[1]
		fluxpars[[f]]$gradients = tmp_pars$coefficients[2:3]
	}
	# Add string to indicate which independent variables were used (for writing prediction file):
	fluxpars$eminputs = 'SWdown, Tair'
	return(fluxpars)
}

train_clustermod = function(sitedata,modelname,met_varnames,flux_varnames,logfilename,alldata,removeflagged){
	# Nonlinear model predicting flux from SWdown, Tair and RH. Gets linear model parameters for
	# each cluster after clustering met data into 27 clusters, with all requested fluxes to be predicted.
	nfluxes = length(flux_varnames) # number of flux predictions we'll require
	nmetvars = length(met_varnames) # number of met variables in total load request
	nsites = length(alldata) - 1 # number of sites to collate training data over (less current site)

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

	# Collate data for training (this'll suck up some memory):
	metvars = matrix(NA,sitedata$trainset_tsteps,nxvars) # initialise to size of number of training time steps
	fluxes = matrix(NA,nfluxes,sitedata$trainset_tsteps)
	if(removeflagged) fluxqc = matrix(NA,nfluxes,sitedata$trainset_tsteps)
	ctr = 1
	for(s in 1:nsites){
		tsteps = alldata[[ sitedata$train_index[s] ]]$alltsteps # timesteps in this training site
		for(xv in 1:nxvars){ # for each met variable:
			# First determine the index of this variable in the loaded met data:
			vIndex = c(1:nmetvars)[met_varnames==xvarnames[xv]]
			# Then load data to training matrix:
			metvars[ctr:(ctr+tsteps-1),xv] = alldata[[ sitedata$train_index[s] ]]$invars[[ vIndex ]]
		}
		# Then get requested flux data:
		for(f in 1:nfluxes){
			fluxes[f,ctr:(ctr+tsteps-1)] = alldata[[ sitedata$train_index[s] ]]$outvars[[ f ]]
			if(removeflagged) fluxqc[f,ctr:(ctr+tsteps-1)] = alldata[[ sitedata$train_index[s] ]]$outvarsqc[[ f ]]
		}
		ctr = ctr + tsteps
	}
	ntsteps = length(metvars[,1])

	# Report to gauge memory usage of data replication:
	writelog(paste0('  ',modelname,' data written, now clustering ',ntsteps,' timesteps using ',cluster_type,
		'...  (',smem(),'; ',proc.time()[3],'s)'),filename=logfilename)

	# To begin, cluster met data based on ALL data (including time steps that have gap-filled fluxes).
	# This saves reclustering for each flux.

	# First scale data:
	xsd=c();xmn=c()
	newmetvars = matrix(NA,ntsteps,nxvars)
	for(v in 1:nxvars){
		xsd[v] = sd(metvars[,v])
		xmn[v] = mean(metvars[,v])
		newmetvars[,v] = (metvars[,v]-xmn[v]) / xsd[v]
	}
	# Then cluster met variables:
	if(cluster_type=='kmeans'){
		xclst = kmeans(newmetvars,nclst,iter.max = 30,nstart=10,algorithm='Lloyd') # nclst is assumed # clusters
		dataclusters = xclst$cluster
	}else if(cluster_type=='clara'){
		xclst = clara(x=newmetvars,k=nclst,metric = "euclidean",stand = FALSE,
			samples = 50, sampsize = 500, mediods.x=FALSE, keep.data =FALSE, pamLike = TRUE)
	}else if(cluster_type=='som'){
		xclst = som(newmetvars,grid=somgrid(floor(sqrt(nclst)),floor(sqrt(nclst)),"rectangular"),rlen=50)
		nclst = (floor(sqrt(nclst)))^2
		dataclusters = xclst$unit.classif
	}

	# Now find fit between met vars and each flux:
	fluxpars = list()
	for(f in 1:nfluxes){
		fluxpars[[f]] = list()
		writelog(paste0('  ...now regressing ',flux_varnames[f],'...'),filename=logfilename)
		intcpt = c()
		grad = matrix(NA,nclst,nxvars) # regression coefficients
		# Reshape data for Qc info:
		if(removeflagged){
			fluxqc_mask = (fluxqc[f,]==0)
		}else{
			fluxqc_mask = ! is.na(fluxes[f,]) # should all be TRUE
		}
		# Find linear parameters for each cluster:
		for(cl in 1:nclst){
			# Find timesteps that belong to this cluster, after gap-filled data removed:
			tsteps_in_cluster = (dataclusters[fluxqc_mask]==cl)
			# Get linear fit data:
			tmp = lm(fluxes[f,tsteps_in_cluster]~metvars[tsteps_in_cluster,],na.action=na.omit)
			intcpt[cl] = tmp$coefficients[1]
			grad[cl,] = tmp$coefficients[2:(nxvars+1)]
		}

		# Save regression paramters:
		fluxpars[[f]]$clusters = xclst
		fluxpars[[f]]$intercept = intcpt
		fluxpars[[f]]$gradients = grad
		fluxpars[[f]]$metvar_means = xmn
		fluxpars[[f]]$metvar_sd = xsd
		fluxpars[[f]]$nclst = nclst
		fluxpars[[f]]$cluster_type = cluster_type
		fluxpars[[f]]$train_tsteps = ntsteps
	}
	# Add string to indicate which independent variables were used (for writing prediction file):
	fluxpars$eminputs = 'SWdown, Tair, RH'
	return(fluxpars)
}
