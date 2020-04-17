# LoadSortData.R
#
# Functions to load flux tower data for empirical model construction
#
# Gab Abramowitz CCRC/CLEX, UNSW 2020 (gabsun at gmail dot com)

LoadData = function(metfile_dir,fluxfile_dir,met_var,flux_var,
  tmp_data_save_file,logfilename){
  writelog(paste('Memory used @ start LoadData:',smem()),filename=logfilename)

  # Check if local save file for this exists, load if it does:
  if(file.exists(tmp_data_save_file)){
    load(file=tmp_data_save_file)
    writelog(paste0('Loaded site data from ',tmp_data_save_file,
      ' (',smem(),')'),filename=logfilename)
  }else{
    # Otherwise read flux tower data files to collect data:
    # List files in met and flux directories:
    metfiles=list.files(path=metfile_dir)
    fluxfiles=list.files(path=fluxfile_dir)
    # Check number of flux and met files match:
    if(length(metfiles)!=length(fluxfiles)){
      message='Number of met files differs from number of flux files.'
      writelog(message,filename=logfilename)
      stop(message)
    }
    # Make sure first 6 chars of names in both met and flux files match:
    if(!all(substr(metfiles,start=1,stop=6)==substr(fluxfiles,start=1,stop=6))){
      message='Sites or order of sites in met files differs from flux files.'
      writelog(message,filename=logfilename)
      stop(message)
    }
    nsites = length(metfiles) # number of sites
    tsteps_per_site = c()
    # Now start to load data:
    data = list()
    namesvector = c() # will be used to report training site names in output file
    for(s in 1:nsites){ # over each site
      data[[s]] = list()
      data[[s]]$name = substr(metfiles[s],start=1,stop=6) # site name
      namesvector[s] = data[[s]]$name # will be used to report training site names in output file
      data[[s]]$index = s # site index
      data[[s]]$invars = list() # will contain 'input' or met forcing variables
      data[[s]]$outvars = list() # will contain 'output' or flux variables
      #data[[s]]$invarsqc = list() # will contain qc flag for input variables
      data[[s]]$outvarsqc = list() # will contain qc flag for output variables
      data[[s]]$outvarsatts= list() # will contain units and longname for each flux var
      # Get met variables:
      mid = nc_open(paste0(metfile_dir,metfiles[s]),write=FALSE,readunlim=FALSE)
      timing = GetTimingNcfile(mid) # in pals package
      data[[s]]$alltsteps = timing$tsteps # number of timesteps for site
      tsteps_per_site[s] = timing$tsteps # number of timesteps for site, local variable
      for(mv in 1:length(met_var)){
        data[[s]]$invars[[mv]] = ncvar_get(mid,met_var[mv],collapse_degen=FALSE)
        #data[[s]]$invarsqc[[mv]] = ncvar_get(mid,paste0(met_var[mv],'_qc'),collapse_degen=FALSE)
      }
      # Save site ancillary data (lat/lon etc) to write to emp model prediction files:
      data[[s]]$details = get_sitedetails(mid)
      nc_close(mid)
      # Now do the same for flux variables:
      fid = nc_open(paste0(fluxfile_dir,fluxfiles[s]),write=FALSE,readunlim=FALSE)
      timing = GetTimingNcfile(fid) # in pals package
      # Make sure number of timesteps match with met file:
      if(timing$tsteps != data[[s]]$alltsteps){
        message=paste('Number of time steps differs between met and flux files in',data[[s]]$name)
        writelog(message,filename=logfilename)
        stop(message)
      }
      for(fv in 1:length(flux_var)){
        data[[s]]$outvars[[fv]] = ncvar_get(fid,flux_var[fv],collapse_degen=FALSE)
        data[[s]]$outvarsqc[[fv]] = ncvar_get(fid,paste0(flux_var[fv],'_qc'),collapse_degen=FALSE)
        data[[s]]$outvarsatts[[fv]] = ncatt_get(fid,flux_var[fv])
      }
      nc_close(fid)
      # Report to log file:
      writelog(paste0('Read ',paste(met_var,collapse=', '),', ',paste(flux_var,collapse=', '),
        ' at site ',data[[s]]$name,'; memory used: ',smem(),'; time ',proc.time()[3]),
        filename=logfilename)
    }
    # Now collect indexing infromation for each site to enable training:
    for(s in 1:nsites){ # over each site
      # Create index vector for all sites but this one (training mask):
      data[[s]]$train_index = CreateTrainIndex(sitenumber=s,numberofsites=nsites)
      data[[s]]$trainset_tsteps = sum(tsteps_per_site[data[[s]]$train_index])
      data[[s]]$trainset_names = paste0(namesvector[data[[s]]$train_index],collapse=', ')
    }
    writelog(paste('Writing site data to',tmp_data_save_file,'to speed loading next
      time...'),filename=logfilename)
    save(data,file=tmp_data_save_file)
  }
  return(data)
}

get_sitedetails = function(mid){
  detail = list()
  detail$timeunits=ncatt_get(mid,'time','units')
	detail$timedata=ncvar_get(mid,'time')
	detail$lat=ncvar_get(mid,'latitude')
	detail$lon=ncvar_get(mid,'longitude')
  return(detail)
}

CreateTrainIndex = function(sitenumber,numberofsites){
  # Create index vector that excludes prediction site from index vector of all sites
  full_index = c(1:numberofsites)
  index_with_site_as_zero = ifelse(full_index==sitenumber,yes=0,no=full_index)
  logical_full_index = index_with_site_as_zero>0
  shortened_index_vector = full_index[logical_full_index]
  return(shortened_index_vector)
}

smem = function(){
  # Just to make memory reporting to lof file a little more intuitive,
  # Reports in Mb instead of b (ish), with units
  short_mem = ceiling(mem_used()/1000000)
  report = paste0(short_mem,'Mb')
  return(report)
}
