# -*- coding: utf-8 -*-
'''
:author Sven Gastauer
:licence MIT
'''
import pandas as pd
import time


def raw2meta_extract(fn):
	"""
	Reasds raw2 files including GPS and enginerring information

	Parameters
	----------
	fn : string
		Path and filenmae of *.raw2 file

	Returns
	-------
	data : pandas DataFrame
		CTD (Salinity, Temperature, Fluorescence, Pressure), Pitch and Roll, Compass information
	gps : pandas DataFrame
		GPS position information
	zoog : pandas DataFrame
		Zoocam grayscale values

	"""
	pgain = 0.04
	poff = -10
	tgain = 0.001
	toff = -5
	sgain = 0.001
	soff = -1
	delta_t = 8
			
	
	#get file index
	print(time.ctime() + ": Processing "+fn)		
	print(time.ctime() + ": Generating file index...")		
	with open(fn) as f:
		list2 = [row.split()[0] for row in f]
	
	##########################################
	#read files
	##########################################
	
	f = open(fn)
	raw2 = f.readlines()
	f.close()
	
	print(time.ctime() + ": Loading CF_DIVE")		
	
	##########################################
	# CF_DIVE 0F
	##########################################
	
		
	'''
	This packet marks the present:
		Nsurf = Dive-Set Number 
		Ncyc = Cycle Number
		Npro = the profile number
		uxti0 = the UNIX time that the Dive-Set
		uxti1 = The Unix time this specific cycle began
		
		For the 0901 code, the Dive-Set Number is only incremented after 
		surface communications (GPS and SBD) are attempted (multiple cycles 
		between surface communications will not increment the Dive-Set 
		Number, but will increment the Cycle Number).  
		This packet should be used to set Nsurf, Ncyc, Npro for all 
		proceeding packets, until the next CF_DIVE packet is encountered.  
	'''
	
	cf_dive_idx = [i for i, j in enumerate(list2) if j == '0f']
	cf_dive_raw = [raw2[i].split() for i in cf_dive_idx]
	cf_dive = pd.DataFrame(cf_dive_raw)
	cf_dive = cf_dive.iloc[:,1:]
	
	cf_dive.columns = ['Nsurf','Ncyc','Npro','uxti0','uxti1','Dow','Month',
					   'day','Time','Year']
	cf_dive = cf_dive.astype(dtype = {'Nsurf':'int64','Ncyc':'int64',
									  'Npro':'int64','uxti0':'int64',
									  'uxti1':'int64'})
	
	##########################################			
	# CF_PDAT  11
	##########################################
	
	print(time.ctime() + ": Loading CF_PDAT")
	edat_idx = [i for i, j in enumerate(list2) if j == '11']
	edat_raw = [raw2[i].split() for i in edat_idx]
	edat = pd.DataFrame(edat_raw)
	edat = edat.iloc[:,1:9]
	edat.columns = ['Nsurf','Ncyc','Npro','time','pressure','temperature',
					'salinity','fluorescence']
	edat = edat.astype(dtype = {'Nsurf':'int64','Ncyc': 'int64','Npro': 'int64',
								'time':'float','pressure':'float',
								'temperature':'float','salinity':'float',
								'fluorescence':'float'} )
	edat['pressure']=edat['pressure'] * pgain + poff #pressure as a double; step 1 of conversion
	#still need to find pmin and do p=p-pmin to convert to dBar	
	sal_cond = edat['salinity'] > 0
	edat.loc[sal_cond, 'salinity'] = edat.loc[sal_cond,'salinity']  * sgain + soff
	sal_cond = edat['temperature'] > 0
	edat.loc[sal_cond, 'temperature'] = edat.loc[sal_cond,'temperature']  * tgain + toff
	
	for var in ['salinity','temperature','fluorescence']:
		cond = edat[var] <= 0
		edat.loc[cond, var] = float('nan')
	
	edat = pd.merge(edat,cf_dive)
	edat['Dive_start_time'] = pd.to_datetime(edat.uxti0, unit='s')
	edat['Dive_start_time'] = edat['Dive_start_time'].dt.tz_localize('UTC')
	#add time_of_measure
	edat['time_of_measure'] = edat['Dive_start_time'] + pd.to_timedelta(edat['time'].astype('str') + 'seconds')
	#edat.time_of_measure = edat.time_of_measure.dt.tz_localize('UTC')
	edat['time_of_measure_PDT'] = edat.time_of_measure - pd.to_timedelta(delta_t, unit='hours') #transform to local time as defined -8 hours not ST
	#correct pressure
	edat['pressure'] = edat.pressure - edat.pressure.min() #Correct pressure
	
	##########################################				
	#CF_EDAT 21
	##########################################
	
	pr_idx = [i for i, j in enumerate(list2) if j == '21']
	pr_raw = [raw2[i].split() for i in pr_idx]
	pr = pd.DataFrame(pr_raw)
	pr = pr.iloc[:,1:7]
	pr.columns = ['Nsurf','Ncyc','Npro','compass','pitch','roll']
	pr = pr.astype(dtype = {'Nsurf':'int64','Ncyc': 'int64',
							'Npro': 'int64','compass':'float',
							'pitch':'float','roll':'float'})
	pr.loc[:,['compass','pitch', 'roll']] /= 10
	
	print(time.ctime() + "Loading CF_GPS1")
	
	##########################################
	#CF_GPS1--start of dive-set 01
	##########################################
	
	gps1_idx = [i for i, j in enumerate(list2) if j == '01']
	gps1_raw = [raw2[i].split() for i in gps1_idx]
	gps1 = pd.DataFrame(gps1_raw)
	gps1 = gps1.iloc[:,[1,3,4,5,6,13]]
	gps1.columns = ['Nsurf_start','year','yr_day_start','lat_start', 'lon_start',
					'UTC_time_fix_start']
	gps1 = gps1.astype(dtype = {'Nsurf_start':'int64', 'year':'int64',
								'yr_day_start':'float','lat_start': 'float', 
								'lon_start': 'float'})
	
	base_date = pd.to_datetime(gps1['year'].astype('str') + '/01/01 00:00:00')
	
	gps1['UTC_time_fix_start'] = base_date + pd.to_timedelta((gps1['yr_day_start']-1).astype('str') + ' days')
	
	print(time.ctime() + ": Loading CF_GPS2")
	
	##########################################
	#CF_GPS2--end of dive-set 02
	##########################################
	
	gps2_idx = [i for i, j in enumerate(list2) if j == '02']
	gps2_raw = [raw2[i].split() for i in gps2_idx]
	gps2 = pd.DataFrame(gps2_raw)
	gps2 = gps2.iloc[:,[1,3,4,5,6,13]]
	gps2.columns = ['Nsurf_end', 'year','yr_day_end',
					  'lat_end', 'lon_end','UTC_time_fix_end']
	gps2 = gps2.astype(dtype = {'Nsurf_end':'int64', 'year':'int64',
								'yr_day_end':'float','lat_end': 'float', 
								'lon_end': 'float'})
	
	base_date = pd.to_datetime(gps2['year'].astype('str') + '/01/01 00:00:00')
	
	gps2['UTC_time_fix_end'] = base_date + pd.to_timedelta((gps2['yr_day_end']-1).astype('str') + ' days')
	print(time.ctime() + "Loading CF_ZOOG")
	
	##########################################
	#CF_ZOOG this is the zooglider grayscale value
	##########################################
	
	zoog_idx = [i for i, j in enumerate(list2) if j == 'b4']
	zoog_raw = [raw2[i].split() for i in zoog_idx]
	zoog = pd.DataFrame(zoog_raw)
	#dt = pd.to_datetime(zoog.iloc[:,7] +' '+ zoog.iloc[:,8] + ' ' + zoog.iloc[:,9] + 
	#					' ' + zoog.iloc[:,10] + ' ' + zoog.iloc[:,11])
	zoog = zoog.iloc[:,[1,2,3,4,5,6]]
	#zoog['date'] = dt
	zoog.columns = ['zstart', 'zstop','n_img','n_err', 'avg','unix_tstamp']
	zoog = zoog.astype(dtype = {'zstart':'int64', 'zstop':'int64','n_img':'int64','n_err': 'int64', 'avg': 'float', 'unix_tstamp':'float'})
	zoog['UTC_time'] = pd.to_datetime(zoog.unix_tstamp[0], unit='s')
	zoog.UTC_time = zoog.UTC_time.dt.tz_localize('UTC')
	zoog['PDT_time'] = zoog.UTC_time - pd.to_timedelta(delta_t, unit='hours')

	##########################################
	#Export
	##########################################
	
	print(time.ctime() + ": Preparing data for export")
	
	##GPS
	gps = pd.merge(gps1, gps2, left_on = 'Nsurf_start', right_on = 'Nsurf_end')
	gps = gps[['Nsurf_start','Nsurf_end', 'UTC_time_fix_start', 'UTC_time_fix_end','lon_start', 'lon_end', 'lat_start', 'lat_end']]
	
	## Data
	data = pd.concat([edat, pr.iloc[:,3:]], sort=False,axis=1)
	#only keep important info
	data = data[['Nsurf','Ncyc','Npro', 'pitch', 'roll', 'compass','pressure','temperature','salinity','fluorescence','uxti0','Dive_start_time','time_of_measure','time_of_measure_PDT']]
	print(time.ctime() + ": Completed")
	return data, gps, zoog