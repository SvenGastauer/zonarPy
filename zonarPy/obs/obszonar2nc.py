# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 09:15:35 2019

@author: sveng
"""

from netCDF4 import Dataset
import xarray as xr
import numpy as np
import pandas as pd
import time
import glob
import os
import os.path

from zonar_reader import Zonar
from meta_reader import raw2meta_extract
from read_sat import read_sat
from help_fun import compute_c, absorption

def read_all_zonar(zdir,start=0,end=0):
	"""
	reads acoustic data from a given folder containing the raw data

	Parameters
	----------
	zdir : string
		Path to acoustic raw folder.
	start : integer, optional
		Dive number at which to start reading the files, 0 starts at the first available dive. The default is 0.
	end : integer, optional
		Dive number at which to stop reading the files, 0 stops at the last available file. The default is 0.

	Returns
	-------
	start : pandas DataFrame
		starting values for the acoustic raw data. Includes non calibration settings.
	all_raws : pandas DataFrame
		Acoustic data in raw counts, by beam, burst, and Ping. Includes timestamp and available meta information.

	"""
	acfn = glob.glob(zdir + '\\B*')
	if len(acfn)>0:
		if end == 0 : end = len(acfn)
		if end > start:
			all_raws = pd.DataFrame()
			z=Zonar()
			#get starting values for frequencies
			start = z.read_one_dive(z.read_raw(acfn[1]))[2]
			#z.add_depth()
		for a in range(0, end):
			all_raws = all_raws.append(z.read_one_dive(z.read_raw(acfn[a]))[0])
			all_raws['gn']= np.array(start['gn'][all_raws['beam']-1])
	else:
		start=[];all_raws=[]
	return start, all_raws

def get_cal(filename, cal0):
	"""
	Converts calibraiton information into xarray for nc inclusion

	Parameters
	----------
	filename : string
		Output nc filename.
	all_raws : pandas DataFrame
		all_raws acoustic data and meta data
	cal0 : defaultdict
		defaultdict containing the acoustic calibraiton information.

	Returns
	-------
	None.
	Adds cal to nc file (filename)

	"""
	print(time.ctime() + ': Gathering Calibration data...')
	cal=pd.DataFrame()
	cal['Frequency'] = cal0['Frequency']
	cal['Gain'] = cal0['Gain']
	cal['gn'] = cal0['gn']
	cal['beam_deg'] = cal0['beam_deg']
	cal['beam_rad'] = cal0['beam_rad']
	cal['tau'] = cal0['tau']
	cal['blank'] = cal0['blank']
	cal['dt'] = cal0['dt']
	cal['tPing'] = cal0['tPing']
	cal['tScan'] = cal0['tScan']
	cal['nScan'] = cal0['nScan']
	cal['tWait'] = cal0['tWait']
	cal['nBin'] = cal0['nBin'] 
	cal['TSGain'] = cal0['TS_Gain']
	cal['CalNoise'] = cal0['CalNoise']
	cal['SoureLevel'] = cal0['sl']
	cal['alpha'] = cal0['alpha']
	cal['c'] = cal0['cspeed']
	cal['TS_cal'] = cal0['TS_cal']
	cal['TS0'] = cal0['TS0']
	cal['Gain_TS'] = cal0['TS_Gain']
	
	cal_arr = xr.Dataset(cal)

	cal_arr.Frequency.attrs['units'] = 'kHz'
	cal_arr.Frequency.attrs['standard_name'] = 'acoustic_frequency'
	
	cal_arr.Gain.attrs['Gain'] = 'dB re 1m-1'
	cal_arr.Gain.attrs['standard_name'] = 'System_Gain'
	
	cal_arr.Gain.attrs['Gain_TS'] = 'dB re 1m2'
	cal_arr.Gain.attrs['standard_name'] = 'Calibration_Gain'
	
	cal_arr.tau.attrs['units'] = 's'
	cal_arr.tau.attrs['standard_name'] = 'Pulse_Duration'
	
	cal_arr.beam_deg.attrs['units'] = 'degrees'
	cal_arr.beam_deg.attrs['standard_name'] = '3dB_beamwidth_in_degrees'
	
	cal_arr.beam_rad.attrs['units'] = 'radians'
	cal_arr.beam_rad.attrs['standard_name'] = '3dB_beamwidth_in_radians'
	
	cal_arr.alpha.attrs['standard_name'] = 'attenuation_coefficient_during_calibration'
	
	cal_arr.dt.attrs['units'] = 'mirco-seconds'
	cal_arr.dt.attrs['standard_name'] = 'Period_between_scans'
	cal_arr.dt.attrs['description'] = 'Normally dt is 200 us which equals sampling rate of 5 kHz'
	
	cal_arr.blank.attrs['units'] = 'mirco-seconds'
	cal_arr.blank.attrs['standard_name'] = 'Period_between_end-of-transmit_and_the_first_scan'
	
	cal_arr.tPing.attrs['units'] = 'ms'
	cal_arr.tPing.attrs['standard_name'] = 'Ping_interval'
	
	cal_arr.tScan.attrs['units'] = 'ms'
	cal_arr.tScan.attrs['standard_name'] = 'Scan_duration'
	
	cal_arr.to_netcdf(filename, 'a', group='Zonar/Calibration')

	
#environment
def get_env(filename,env):
	"""
	Adds CTD data to xarray for nc conversion

	Parameters
	----------
	filename : string
		nc output filename.
	env : pandas DataFrame
		Environmental data

	Returns
	-------
	None.
	Adds environment to nc file (filename)

	"""
	print(time.ctime() + ': Gathering Environment data...')
	var = ['temperature','salinity','fluorescence']
	xy = ['pressure','Nsurf','Dive_start_time']#,'lon_start','lat_start','lon_end','lat_end']
	env['numtime'] = env.time_of_measure.astype(int)
	env = env.dropna(subset=['fluorescence'])
	print('Processing: time')
	measuretime = pd.pivot_table(env,values='numtime', index=xy[0], columns=xy[1])
	print('Processing: ' + var[0])
	tp = pd.pivot_table(env,values=var[0],index=xy[0], columns=xy[1:])
	print('Processing: ' + var[1])
	sp = pd.pivot_table(env,values=var[1],index=xy[0], columns=xy[1])
	print('Processing: ' + var[2])
	fp = pd.pivot_table(env,values=var[2],index=xy[0], columns=xy[1])
	e_arr = xr.Dataset(data_vars = {var[0]:(('Depth','Dive'),np.array(tp)),
								 var[1]:(('Depth','Dive'),np.array(sp)),
								 var[2]:(('Depth','Dive'),np.array(fp)),
								 'time':(('Depth','Dive'),np.array(measuretime))},
								 coords={'Dive':np.array(tp.columns.get_level_values('Nsurf')),
				 'Depth': np.array(tp.index)})
	#e_arr['StartTime'] = ('Dive',np.array(tp.columns.get_level_values('Dive_start_time')))
	e_arr.temperature.attrs['standard_name'] = 'Ambient_Water_Temperature'
	e_arr.temperature.attrs['units'] = 'Degrees_Celsius'
	e_arr.time.attrs['standard_name'] = 'time_of_measure'
	e_arr.salinity.attrs['standard_name'] = 'Ambient_Water_Salinity'
	e_arr.salinity.attrs['units'] = 'PSU'
	e_arr.fluorescence.attrs['standard_name'] = 'Ambient_Water_Fluorescence'
	e_arr.fluorescence.attrs['units'] = 'RFU'
	e_arr.Depth.attrs['standard_name'] = 'Depth_of_measurement'
	e_arr.Depth.attrs['units'] = 'dBar'
	e_arr.Dive.attrs['standard_name'] = 'Dive_number'
	e_arr.to_netcdf(filename, 'a', group = 'Environment' )
#gps
def get_gps(filename, gps):
	"""
	Adds gps data to xarray for nc conversion

	Parameters
	----------
	filename : string
		nc output filename.
	gps : pandas DataFrame
		GPS data

	Returns
	-------
	None.
	Adds gps to nc file (filename)

	"""
	print(time.ctime() + ': Gathering GPS data...')
	gps = gps.reset_index()
	gps = gps.set_index('Nsurf_start')
	gps_arr = xr.Dataset(data_vars ={'Lon_start':('Dive',gps.lon_start),
								  'Lon_end':('Dive',gps.lon_end),
								  'Lat_start':('Dive',gps.lat_start),
								  'Lat_end':('Dive',gps.lat_end),
								  'Time_start':('Dive',gps.UTC_time_fix_start),
								  'Time_end':('Dive',gps.UTC_time_fix_end)})
	gps_arr.Lon_start.attrs['standard_name'] = 'Longitude_at_the_start_of_the_dive'
	gps_arr.Lon_end.attrs['standard_name'] = 'Longitude_at_the_end_of_the_dive'
	gps_arr.Lat_start.attrs['standard_name'] = 'Latitude_at_the_start_of_the_dive'
	gps_arr.Lat_end.attrs['standard_name'] = 'Latitude_at_the_end_of_the_dive'
	gps_arr.Time_start.attrs['standard_name'] = 'Time_at_the_start_of_the_dive - UTC'
	gps_arr.Time_end.attrs['standard_name'] = 'Time_at_the_end_of_the_dive - UTC'
	gps_arr.to_netcdf(filename, 'a', group = 'GPS' )
		

#zooglider grey values
def get_zoog(filename, zoog):
	"""
	Adds zooglider data to xarray for nc conversion

	Parameters
	----------
	filename : string
		nc output filename.
	zoog : pandas DataFrame
		zoog data

	Returns
	-------
	None.
	Adds zoog to nc file (filename)

	"""
	zoog = zoog.reset_index()
	zoog = zoog.set_index('UTC_time','PDT_time')
	zoog_arr = xr.Dataset(zoog)
	zoog_arr.to_netcdf(filename, 'a', group = 'zoog' )
	
#acoustics
def get_ac(filename, all_raws, beam):
	"""
	Adds CTD data to xarray for nc conversion

	Parameters
	----------
	filename : string
		nc output filename.
	all_raws : pandas DataFrame
		acoustic data
	beam : integer
		slect beam for which to get the acoustic data

	Returns
	-------
	None.
	Adds acoustic data to nc file (filename)

	"""
	print(time.ctime() + ': Gathering Acoustic data for beam: ' + str(beam) )
	
	sub_raws0 = all_raws[(all_raws.beam == beam) & (all_raws.Ping == 0)]
	sub_raws1 = all_raws[(all_raws.beam == beam) & (all_raws.Ping == 1)]
	sub_raws2 = all_raws[(all_raws.beam == beam) & (all_raws.Ping == 2)]
	sub_raws3 = all_raws[(all_raws.beam == beam) & (all_raws.Ping == 3)]
	
	ac_temp0 = pd.pivot_table(sub_raws0, values='Raw',index = ['nscan'], columns = ['ping_time','nBurst','dive'])
	ac_temp1 = pd.pivot_table(sub_raws1, values='Raw',index = ['nscan'], columns = ['ping_time','nBurst','dive'])
	ac_temp2 = pd.pivot_table(sub_raws2, values='Raw',index = ['nscan'], columns = ['ping_time','nBurst','dive'])
	ac_temp3 = pd.pivot_table(sub_raws3, values='Raw',index = ['nscan'], columns = ['ping_time','nBurst','dive'])
	
	dz = pd.pivot_table(sub_raws0, values='dz',index = ['nscan'], columns = ['ping_time','nBurst','dive'])

	
	ac_arr = xr.Dataset(data_vars={'Ping1':(('nScan','Burst'),np.array(ac_temp0)),
								'Ping2':(('nScan','Burst'),np.array(ac_temp1)),
								'Ping3':(('nScan','Burst'),np.array(ac_temp2)),
								'Ping4':(('nScan','Burst'),np.array(ac_temp3)),
								'Range':(('nScan','Burst'),np.array(dz))},
					 coords = {'nScan':np.array(ac_temp1.index),
				  'Burst': np.array(ac_temp0.columns.get_level_values('nBurst'))})
	ac_arr['Time'] = ('Burst',np.array(ac_temp0.columns.get_level_values('ping_time')))
	ac_arr['Dive'] = ('Burst' ,np.array(ac_temp0.columns.get_level_values('dive')))
	
	
	#set meta data
	ac_arr.Ping1.attrs['units'] = '40 counts per dB re V'
	ac_arr.Ping1.attrs['standard_name'] = 'Raw_acoustic_data'
	ac_arr.Ping2.attrs['units'] = '40 counts per dB re V'
	ac_arr.Ping2.attrs['standard_name'] = 'Raw_acoustic_data'
	ac_arr.Ping3.attrs['units'] = '40 counts per dB re V'
	ac_arr.Ping3.attrs['standard_name'] = 'Raw_acoustic_data'
	ac_arr.Ping4.attrs['units'] = '40 counts per dB re V'
	ac_arr.Ping4.attrs['standard_name'] = 'Raw_acoustic_data'
	ac_arr.Time.attrs['standard_name'] = 'Datetime_of_the_burst'
	ac_arr.Burst.attrs['standard_name'] = 'Burst_number'
	ac_arr.Range.attrs['standard_name'] = 'Depth'
	ac_arr.Range.attrs['units'] = 'dBar'
	ac_arr.nScan.attrs['standard_name'] = 'Sample_number'
	ac_arr.Burst.attrs['description'] = 'Each burst consists of 4 consequetive pings'
	
	ac_arr.to_netcdf(filename, 'a',group='Zonar/' + 'Beam_' + str(beam))

	return(ac_arr)
	#return ac_arr	
#get sat file
def  get_sat(filename,header,gps,engineering,profile, zoocam, zonar, misc):
	"""
	add satllite transmitted data to xarray for nc conversion

	Parameters
	----------
	filename : string
		nc filename.
	header : TYPE
		DESCRIPTION.
	gps : TYPE
		DESCRIPTION.
	engineering : TYPE
		DESCRIPTION.
	profile : TYPE
		DESCRIPTION.
	zoocam : TYPE
		DESCRIPTION.
	zonar : TYPE
		DESCRIPTION.
	misc : TYPE
		DESCRIPTION.

	Returns
	-------
	None.

	"""
	print(time.ctime() + ': Processing sat data...' )
	
	def create_var(dic, var, nam=None,vnam='Sat/misc', filename=filename):
		if nam is None: nam = var
		for i in range(len(var)):
			print(time.ctime() + ': Writing ' + vnam + '/' + nam[i] )
			if len(dic[var[i]]) > 0: 
				xr.Dataset(dic[var[i]]).to_netcdf(filename, 'a',group = vnam + '/' + nam[i])
	
	#gps
	if len(gps)>0: xr.Dataset(gps).to_netcdf(filename, 'a',group = 'Sat/gps')
	
	
	#profile
	if len(profile)>0: xr.Dataset(profile['Profile_data']).to_netcdf(filename, 'a',group = 'Sat/profile/Profile_data')
	if len(profile)>0: xr.Dataset(pd.DataFrame({'Start':profile['Start_Profile_data']})).to_netcdf(filename, 'a',group = 'Sat/profile/Start')

	#zoocam
	#var = [k for k in zoocam.keys()]
	create_var(zoocam, ['Zoocam'], vnam='Sat/zoocam')
	#zonar
	var = [k for k in zonar.keys()]
	create_var(zonar, var, vnam='Sat/zonar')
	
	#engineering
	for k in engineering.keys():
		print(k)
		if len(engineering[k])>0: 
			if k == 'Engineering_TS' :
				eng = engineering[k]
				engdf = pd.DataFrame({'sec_since_start' : eng.iloc[:,1:13][eng['Code'] == 'ET00'].values.flatten().astype(int),
						  'pressure_counts' : eng.iloc[:,1:13][eng['Code'] == 'ET01'].values.flatten().astype(int),
						  'heading10deg' : eng.iloc[:,1:13][eng['Code'] == 'ET02'].values.flatten().astype(int),
						  'pitch10deg' : eng.iloc[:,1:13][eng['Code'] == 'ET03'].values.flatten().astype(int),
						  'roll10deg' : eng.iloc[:,1:13][eng['Code'] == 'ET04'].values.flatten().astype(int),
						  'pitchPotCounts' : eng.iloc[:,1:13][eng['Code'] == 'ET05'].values.flatten().astype(int),
						  'rollPotCounts' : eng.iloc[:,1:13][eng['Code'] == 'ET05'].values.flatten().astype(int),
						  'dive#' : np.tile(eng['dive#'][eng['Code'] == 'ET00'],12).flatten().astype(int)})
				dd = xr.Dataset(engdf)
				dd.set_coords(['dive#', 'sec_since_start'])  
				dd.to_netcdf(filename, 'a',group = str('Sat/engineering/'+k))
			else:
				xr.Dataset(engineering[k]).to_netcdf(filename, 'a',group = str('Sat/engineering/'+k))
	
	#header
	var = [k for k in header.keys()]
	var.remove('argos')
	var.remove('optical_sensor')
	create_var(header, var, vnam='Sat/header')
	
	#misc
	print(time.ctime() + ': Writing misc...' )
	var = [k for k in misc.keys()]
	var.remove('Shore_comm')
	nam = ['Dive Info', 'Mission ID','Email message','Engineering',\
		'Parameters list','Parameter values' ]
	create_var(misc, var,nam, vnam='Sat/misc')
	
	
#create nc file
def generate_nc(filename, env, gps, zoog, all_raws, miss, header,gps_sat,
				engineering,profile, zoocam, zonar, misc, cal):
	print(time.ctime() + ': Preparing Zonar data...')
	### ZONAR Preparation
	#get number of pings
	all_raws['sample_time'] = pd.to_datetime(all_raws['dive_time']).values.astype(np.int64) +  \
		( all_raws['dt']/1000 * (all_raws['nscan'] + 1) + all_raws['blank'] + \
   all_raws['tau'] * 1000 + all_raws['Ping'] * all_raws['tPing'])/1000 * 10**9
	#all_raws['sample_time'] = pd.to_datetime(all_raws['real_time'])
	all_raws['ping_time'] = pd.to_datetime(all_raws['dive_time']).values.astype(np.int64) +  \
		(all_raws['Ping'] * all_raws['tPing'])/1000 * 10**9
	all_raws['ping_time'] = pd.to_datetime(all_raws['ping_time'])
	        
	#pings = all_raws.Ping.unique()
	beams = all_raws.beam.unique()
	
	print(time.ctime() + ': Create netcdf file...')
	
	dataset = Dataset(filename, 'w')
	#add meta data
	dataset.description = header['mission'].Description[0]
	dataset.history = 'Created ' + time.ctime(time.time())
	dataset.serial = header['vehicle'].SN[0]
	dataset.Spray_serial = header['mission'].Spray_SN[0]
	dataset.EEPROM_Ver = header['vehicle'].EEPROM__Ver[0]
	dataset.Deployment_ID = header['mission'].Deployment_ID[0]
	dataset.Mission_Name = header['mission'].E_Name[0]
	dataset.argos_ID = header['argos'].Argos_ID[0]
	
	dataset.mission = miss
	
	dataset.close()
	
	#add environment
	get_env(filename, env)
	#add GPS
	get_gps(filename, gps)
	#add zoog
	get_zoog(filename, zoog)
	#add cal
	get_cal(filename, cal)
	#add sat
	get_sat(filename, header,gps,engineering,profile, zoocam, zonar, misc)
	
	#add raw acoustics
	for beam in beams:
		_ = get_ac(filename, all_raws,beam)
		
		

def get_cal_val(fn, date):
	caldf = pd.read_table(fn, sep='\t{1,}',comment='#',header=None, \
					   engine='python')
	caldf.columns = pd.read_table(fn, sep='\s{2,}',skiprows=7, nrows=1,\
							   header=None, engine='python').values[0]
	caldf = caldf.set_index(pd.to_datetime(caldf['#monYr'], format='%b%y'))
	return caldf.iloc[caldf.index.get_loc(date,method='nearest')]


def get_missions(mdir):
	return [i for i in next(os.walk(mdir))[1] if '20' in i]

def missions_to_nc(mdir, missions, outdir='', acdir='zonar_flash',
				   raw2dir = 'raw2',satdir = 'in-situ_sat_file',
				   d_start=0,d_end=0, TS0=None, Scal=0, Tcal=20, dcal=5, 
				   Gain = [54,54], calfn=None, force=False):
	for miss in missions:
		print(time.ctime() + ': Processing ' + miss )
		outfn = outdir + miss + '.nc'
		if os.path.isfile(outfn) and force==False:
			print('NetCDF for ', miss,' already exists...skipping...')
			continue
		else:
			#get raw_dirs
			rdir = [i for i in next(os.walk(mdir+miss))[1]if 'raw_data' in i]
			if len(rdir) > 0:
				raw_dir = mdir + miss + '\\' + rdir [0] 
			else:
				print(time.ctime() + ': No RAW DATA FOLDER found for ' + \
		  miss + ' ...skipping mission...')
				continue
			#get satellite data
			sat_dir = mdir + miss + '\\' + satdir
			satfn = glob.glob(sat_dir +'\\*.sat')
			if len(satfn)>0:
				header, gps_sat,engineering, profile, zoocam, zonar,\
					 misc = read_sat(satfn[0])
			else:
				print(time.ctime() + ': No SAT file found for ' + miss + \
		  ' ...skipping mission...')
				continue
			
			#get acoustics zonar data
			ac_dir = raw_dir + '\\' + acdir + '\\'
			start, all_raws = read_all_zonar(ac_dir,d_start, d_end)
			if len(all_raws) == 0:
				print(time.ctime() + ': No ACOUSTIC DIVE DATA found for ' + \
		  miss + ' ...skipping mission...')
				continue
			#get acoustic calibration data
			# get the cal information from file if available
			if calfn is not None:
				print(time.ctime() + ': Reading acoustic cal values from file...')
				cv = get_cal_val(calfn, pd.to_datetime(all_raws.start_time.values[0]))
				NL = cv[cv.index.isin(['N[200]','N[1000]'])].values
				TS = cv[cv.index.isin(['T[200]','T[1000]'])].values
				print('NL = ' + str(NL) + ' TS  = ' + str(TS))
				Spray = cv.Spray
				znr = cv.Zonar
				monYr = cv[0]
			else:
				ctemp = Zonar().init_cal()
				NL = ctemp['CalNoise']
				TS = [87, 82]
				Spray = header['vehicle']['SN'][0]
				znr = 'NA'
				monYr = pd.to_datetime(all_raws.start_time.values[0]).strftime('%b%y')
				
			#update calibration values
			print(time.ctime() + ': Updating calibration data')
			#TS0 - the theoretical TS of the used sphere
			if TS0 is None: TS0 = [-53.47, -50.92]
			if Gain is None: Gain = [54, 54]
			# get the source level approximation
			alpha = np.array([absorption(f=f, S=np.array(Scal), 
								T=np.array(Tcal), D=np.array([dcal]))/1000 \
					 for f in start['freq']]).flatten()
			c = compute_c(dcal,Scal,Tcal)
			
			
			cal = Zonar().init_cal(CalNoise = NL,
						  Gain = Gain,
						  TS_cal = TS,
						  TS0 = TS0,
						  alpha = alpha,
						  cspeed = c,
						  Spray = Spray,
						  Zonar = znr,
						  monYr = monYr,
						  tau = start.pulse.values,
						  blank = start.blank.values,
						  dt = start.dt.values,
						  tScan = start.tScan.values,
						  tPing = start.tPing.values,
						  nScan = start.nScan.values,
						  tWait = start.tWait.values,
						  nBin = start.nBin.values,
						  gn = start.gn.values,
						  Frequency = start.freq.values)
			del cal['Noise']
			# SL = Source level, 40log10(z) + 2 alpha z = 2 x Transmission 
			#Loss, G0 = system gain
			SL = np.array(NL)/cal['gn'] + np.array(cal['Gain']) +\
				 40 * np.log10(dcal)+ 2 * alpha * dcal 
			cal['sl'] = SL
			#get expected TS in counts
			TS_c = TS0 + SL - 40*np.log10(dcal) - 2 * alpha * dcal
			#get the TS cal Gain value
			Gcal = TS_c - np.array(TS) #dB re V
			cal['TS_Gain'] = Gcal
				
			#get physical data
			raw2_dir = raw_dir + '\\' + raw2dir
			rawfn = glob.glob(raw2_dir + '\\*.raw2')
			if len(rawfn) > 0:
				env, gps, zoog = raw2meta_extract(rawfn[0])
			else:
				print(time.ctime() + ': No RAW2 file found for ' + 
		  miss + ' ...skipping mission...')
				continue
				
			#generate netCDF file
			generate_nc(filename=outfn,env=env, gps=gps, zoog=zoog, 
			   all_raws=all_raws, miss=miss,header=header,gps_sat=gps_sat,
			   engineering=engineering,profile=profile, zoocam=zoocam, 
			   zonar=zonar, misc=misc,cal=cal)
			
def get_AllRoiCounts(rdir='Z:\\Zooglider_end_to_end_results', pat = '/20*/',
					 outdir=''):
	"""
	Get all the ROI counts from the Zoocam data

	Parameters
	----------
	rdir : string, optional
		Path to the folder containing the Zooglider data. This can be a parent\
			 folder as the pattern matching will be recursive. The default is \
				 'Z:\\Zooglider_end_to_end_results'.
	pat : string, optional
		Regex Pattern to be matched. The default is '/20*/'.
	outdir : string, optional
		Folder to which the csv files should be written. The default is ''.

	Returns
	-------
	cdat : TYPE
		DESCRIPTION.

	"""
	ldirs = glob.glob(rdir + pat)
	for ldir in ldirs:
		miss = os.path.basename(os.path.normpath(ldir))
		print(time.ctime() + ': Processing ' + miss )
		cdir = glob.glob(ldir + '/*_converted_to_csv')[0] + '\\Physical_data_at_Interpolated_Depths\\'
		int_fn = sorted(glob.glob(cdir + '*_interpolated_data.csv'))
		meas_fn = sorted(glob.glob(cdir + '*_roi_counts.csv'))
		if len(meas_fn)>0:
			cdat = pd.DataFrame()
			for i in range(len(int_fn)):
				print(time.ctime() + ': ' + str(i+1) + ' of ' + str(len(int_fn)) + ' - ' + str(np.round((i+1)/len(int_fn) * 100)) + '% done')
				tmp = pd.read_csv(int_fn[i])
				tmp2 = pd.read_csv(meas_fn[i])
				tmp2['filename'] = tmp2['Unnamed: 0'].str.split('.png', n=1,expand=True)[0]
				tmp2.drop(columns=['Unnamed: 0'], inplace=True)
				cdat = cdat.append(tmp.merge(tmp2))
				cdat.loc[:,'mission'] = miss
				outfn = outdir + miss + '.csv.gz'
			print(time.ctime() + ': Writing ' + outfn)
			cdat.to_csv(outfn, compression = 'gzip')
		else:
			print(time.ctime() + ': No converted_to_csv folder found for ' + ldir) 
			continue
	return cdat