# Train.R
#
# Functions to train empirical models
#
# Gab Abramowitz CCRC/CLEX, UNSW 2020 (gabsun at gmail dot com)

TrainEmpiricalModels = function(sitedata,emodels,met_varnames,flux_varnames,logfilename,alldata,removeflagged=TRUE){
	# For a given site (coming from lapply/parlapply), trains all empirical models for all flux variables.
	library(pals)
	source('Train.R')
	writelog(paste0(' Training model for site ',sitedata$name,
		':            (time ',proc.time()[3],'s)'),filename=logfilename)
	nemodels = length(emodels) # number of requested empirical model types
	emod_pars = list() # empirical model parameters
	# Fetch empirical model function names:
	emod_functions = mget(paste0('train_',emodels),mode='function',inherits=TRUE,envir = globalenv())
	# Call empirical model functions - each call get model for all flux variable predictions:
	for(em in 1:nemodels){
		emod_pars[[em]] = emod_functions[[em]](sitedata,met_varnames,flux_varnames,logfilename,alldata,removeflagged)
	}
	# save so we can identify site for these model parameters later
	emod_pars$siteindex = sitedata$index
	return(emod_pars)
}

train_1lin = function(sitedata,met_varnames,flux_varnames,logfilename,alldata,removeflagged){
	# Linear model predicting flux from SWdown only. Gets linea model parameters for
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

train_2lin = function(sitedata,met_varnames,flux_varnames,logfilename,alldata,removeflagged){
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





#
#
# 	# Find out how many time steps in total training set:
# 	ntsteps = 0 # initialise number of time steps
# 	dstart = c() # index in all data vector of data set start
# 	dend = c() # index in all data vector of data set end
# 	for(s in 1:nsites){
# 		fid = open.ncdf(metfiles[s])
# 		timing = GetTimingNcfile(fid) # in PALS package
# 		close.ncdf(fid)
# 		dstart[s] = ntsteps + 1
# 		ntsteps = ntsteps + timing$tsteps
# 		dend[s] = ntsteps
# 	}
# 	yvar=c() # dependent variable
# 	xvar = matrix(NA,ntsteps,nxvars) # declare x variable data matrix
# 	if(removeflagged){
# 		xvar_qc = matrix(NA,ntsteps,nxvars)
# 		yvar_qc=c()
# 	}
# 	# Now poulate x-data matrix and y-data vector:
# 	for(s in 1:nsites){
# 		cat('Fetching data set ',s,'\n')
# 		# Get independent variables:
# 		for(v in 1:nxvars){
# 			# If we're using humidity as a benchmark input, change to relative humidity:
# 			if(xvarnames[v] != 'Qair'){
# 				tmpx = GetFluxnetVariable(xvarnames[v],metfiles[s],'blah')
# 				xvar[(dstart[s]:dend[s]),v] = tmpx$data
#
# 			}else{
# 				tmpx = GetFluxnetVariable(xvarnames[v],metfiles[s],'blah')
# 				tmpTair = GetFluxnetVariable('Tair',metfiles[s],'blah')
# 				tmpPSurf = GetFluxnetVariable('PSurf',metfiles[s],'blah')
# 				xvar[(dstart[s]:dend[s]),v] = Spec2RelHum(tmpx$data,tmpTair$data,tmpPSurf$data)
# 			}
# 			if(removeflagged & tmpx$qcexists){
# 				xvar_qc[(dstart[s]:dend[s]),v] = tmpx$qc
# 			}else if(removeflagged & ! tmpx$qcexists){
# 				# If no QC flag exists, assume data are original:
# 				xvar_qc[(dstart[s]:dend[s]),v] = 1
# 			}
# 		}
# 		tmpy = GetFluxnetVariable(yvarname,fluxfiles[s],'blah')
# 		yvar = c(yvar,tmpy$data)
# 		if(removeflagged & tmpy$qcexists){
# 			yvar_qc = c(yvar_qc,tmpy$qc)
# 		}else if(removeflagged & ! tmpy$qcexists){
# 			# If no QC flag exists, assume data are original:
# 			tmp = c(1:length(tmpy$data))*0 + 1
# 			yvar_qc = c(yvar_qc,tmp)
# 		}
# 	}
# 	# Construct training data set only from observed, and not
# 	# gap-filled data, if requested:
# 	if(removeflagged){
# 		cat('Using only non-gap-filled data to train with... \n')
# 		# Convert to logical mask using default 1=>TRUE, 0=>FALSE
# 		flagmask = (as.logical(xvar_qc[,1]) & as.logical(yvar_qc))
# 		if(nxvars>1){
# 			for(v in 1:(nxvars-1)){
# 				flagmask = (flagmask & as.logical(xvar_qc[,(v+1)]))
# 			}
# 		}
# 		xvar_new = xvar[flagmask,]
# 		yvar_new = yvar[flagmask]
# 		xvar = xvar_new
# 		yvar = yvar_new
# 		cat('...using',length(yvar_new),'of original',ntsteps,'time steps \n')
# 		ntsteps = length(yvar_new)
# 	}
# 	# Choose empirical model type:
# 	if(emod=='kmeans'){
# 		cat('Clustering',ntsteps,'time steps from',nsites,'sites... \n')
# 		# First sd-weight variables:
# 		wxvar = matrix(NA,ntsteps,nxvars) # declare x var data matrix
# 		xsd = c()
# 		xmn = c()
# 		for(v in 1:nxvars){
# 			xsd[v] = sd(xvar[,v])
# 			xmn[v] = mean(xvar[,v])
# 			wxvar[,v] = (xvar[,v]-xmn[v]) / xsd[v]
# 		}
# 		# Cluster dependent variables:
# 		xclst = kmeans(wxvar,eparam,iter.max = 50,nstart=3) # eparam assumed # clusters
# 		cat('At least',min(xclst$size),'data in each cluster \n')
# 		cat('Regressing cluster number: ')
# 		intcpt = c()
# 		grad = matrix(NA,eparam,nxvars) # regression coefficients
# 		for(c in 1:eparam){
# 			cat(c,' ')
# 			tmp = lm(yvar[(xclst$cluster==c)]~
# 				xvar[(xclst$cluster==c),],na.action=na.omit)
# 			intcpt[c] = tmp$coefficients[1]
# 			grad[c,] = tmp$coefficients[2:(nxvars+1)]
# 		}
# 		cat('\n')
# 	}else if(emod=='mlr'){
# 		# Perform linear regression
# 		intcpt = c()
# 		grad = matrix(NA,1,nxvars) # regression coefficients
# 		xclst = NA
# 		xsd = NA
# 		xmn=NA
# 		tmp = lm(yvar~xvar,na.action=na.omit)
# 		intcpt[1] = tmp$coefficients[1]
# 		grad[1,] = tmp$coefficients[2:(nxvars+1)]
# 	}else{
# 		CheckError('Unknown empirical model choice.')
# 	}
# 	empmod = list(xclst=xclst,grad=grad,int=intcpt,
# 		xsd=xsd,xmn=xmn,type=emod)
# 	return(empmod)
# }
