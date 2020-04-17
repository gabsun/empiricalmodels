# WriteFiles.R
#
# Writes empirical model predictions to file
#
# Gab Abramowitz CCRC/CLEX, UNSW 2020 (palshelp at gmail dot com)

write_emp_predictions = function(predictions,emodels,met_varnames,flux_varnames,
	logfilename,alldata,outfile_dir,trainedmodels,removeflagged=TRUE){
	# Write empirical model predictions to file for a single site (coming from lapply/parlapply)
	writelog(paste0(' Writing model predictions for site ',alldata[[predictions$siteindex]]$name,
		':            (time ',proc.time()[3],'s)'),filename=logfilename)
	nemodels = length(emodels) # number of requested empirical model types
	# Write to file separately for each empirical model type:
	for(em in 1:nemodels){
		writemodel(predictions[[em]],predictions$siteindex,met_varnames,
			flux_varnames,logfilename,alldata,removeflagged,modelname=emodels[em],outfile_dir,
			trainedmodels[[predictions$siteindex]][[em]]$eminputs)
	}
}
return()

writemodel = function(prediction,siteindex,met_varnames,flux_varnames,logfilename,
	alldata,removeflagged,modelname,outfile_dir,eminputs){
	# Creates an (empirical) model output file for a single site and all requested fluxes

	# Define netcdf file for storing benchmark simulations:
	missing_value=NcMissingVal # default missing value for all variables
	# Define x, y dimensions
	xd = ncdim_def('x',vals=c(1),units='')
	yd = ncdim_def('y',vals=c(1),units='')
	# Define time dimension:
	td = ncdim_def('time',unlim=TRUE,units=alldata[[siteindex]]$details$timeunits$value,
		vals=alldata[[siteindex]]$details$timedata)
	# Begin defining variables (time independent first):
	vars = list() # initialise
	# Define latitude variable:
	vars[[1]] = ncvar_def('latitude',units='degrees_north',
		dim=list(xd,yd),missval=missing_value,longname='Latitude')
	# Define longitude variable:
	vars[[2]] = ncvar_def('longitude',units='degrees_east',
		dim=list(xd,yd),missval=missing_value,longname='Longitude')

	# Define time varying flux variables for netcdf output file:
	for(fv in 1:length(flux_varnames)){
		# for this flux variable, first get units and longname
		# (just using current site info, but could be any, since they're all the same)
		attnames = names(alldata[[siteindex]]$outvarsatts[[fv]])
		unitsindex = c(1:length(attnames))[attnames=='units']
		longnameindex = c(1:length(attnames))[attnames=='long_name']
		# The define this flux variable:
		vars[[(fv+2)]]=ncvar_def(name=flux_varnames[fv],
			units=alldata[[siteindex]]$outvarsatts[[fv]][[unitsindex]],dim=list(xd,yd,td),
			missval=missing_value,longname=alldata[[siteindex]]$outvarsatts[[fv]][[unitsindex]])
	}

	# Create benchmark netcdf file:
	filename = paste0(outfile_dir,modelname,'/Emp',modelname,'_',alldata[[siteindex]]$name,'.nc')
	ncid = nc_create(filename=filename,vars=vars)

	# Add empirical model attributes to each variable:
	for(fv in 1:length(flux_varnames)){
		ncatt_put(ncid,flux_varnames[fv],'empirical_model_inputs',eminputs)
		ncatt_put(ncid,flux_varnames[fv],'empirical_model_type',modelname)
	}

	# Write global attributes:
	ncatt_put(ncid,varid=0,attname='Production_time',attval=as.character(Sys.time()))
	ncatt_put(ncid,varid=0,attname='Training_sites',attval=alldata[[siteindex]]$trainset_names)
	ncatt_put(ncid,varid=0,attname='Production_source',
		attval=paste('Out of sample empirical model simulation:',alldata[[siteindex]]$name,
		'data was NOT used in training regression parameters for this prediction.'))
	if(removeflagged){ # Note if gap-filled data has been excluded
		ncatt_put(ncid,varid=0,attname='Exclusions',
 			attval='Empirical model trained only on non-gap-filled data')
	}
	ncatt_put(ncid,varid=0,attname='Contact',attval='Gab Abramowitz: gabsun@gmail.com')

	# Add time-independent variable data to file:
	ncvar_put(ncid,'latitude',vals=alldata[[siteindex]]$details$lat)
	ncvar_put(ncid,'longitude',vals=alldata[[siteindex]]$details$lon)

	# Add each variable to file:
	for(fv in 1:length(flux_varnames)){
		ncvar_put(ncid,flux_varnames[fv],vals=prediction[fv,])
	}
	nc_close(ncid)
	return()

	# add version of this code!

}
