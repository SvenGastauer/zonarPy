# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 15:55:13 2020

@author: sveng
"""
import os
import glob
import time
import numpy as np
import pandas as pd
import xarray as xr
import datetime
from astral import Astral

from help_fun import haversine, compute_c, absorption
from zonar_reader import Zonar
from meta_reader import raw2meta_extract
from read_sat import read_sat
from zonar2nc import generate_nc


def get_paths(mdir, acdir='zonar_flash',raw2dir = 'raw2',satdir = 'in-situ_sat_file'):
	"""
	Parameters
	----------
	mdir : TYPE
		DESCRIPTION.
	acdir : TYPE, optional
		DESCRIPTION. The default is 'zonar_flash'.
	raw2dir : TYPE, optional
		DESCRIPTION. The default is 'raw2'.
	satdir : TYPE, optional
		DESCRIPTION. The default is 'in-situ_sat_file'.

	Returns
	-------
	ac_dir : TYPE
		DESCRIPTION.
	rawfn : TYPE
		DESCRIPTION.
	satfn : TYPE
		DESCRIPTION.
	miss : TYPE
		DESCRIPTION.

	"""
	rdir = [i for i in next(os.walk(mdir))[1]if 'raw_data' in i]
	miss = os.path.basename(os.path.normpath(mdir))
	if len(rdir) > 0:
		raw_dir = mdir + '\\' + rdir [0] 
	else:
		print(time.ctime() + ': No RAW DATA FOLDER found for ' + miss + ' ...skipping mission...')
		return [],[],[],[]
	#get acoustics zonar data
	ac_dir = raw_dir + '\\' + acdir + '\\'
	#get physical data
	raw2_dir = raw_dir + '\\' + raw2dir
	rawfn = glob.glob(raw2_dir + '\\*.raw2')
	#get satellite data
	sat_dir = mdir + '\\' + satdir
	satfn = glob.glob(sat_dir +'\\*.sat')
	return ac_dir, rawfn, satfn, miss


def get_cal_val(fn, date):
	caldf = pd.read_table(fn, sep='\t{1,}',comment='#',header=None, engine='python')
	caldf.columns = pd.read_table(fn, sep='\s{2,}',skiprows=7, nrows=1,header=None, engine='python').values[0]
	caldf = caldf.set_index(pd.to_datetime(caldf['#monYr'], format='%b%y'))
	return caldf.iloc[caldf.index.get_loc(date,method='nearest')]

		
def read_all_zonar(zdir,start=0,end=0):
	acfn = glob.glob(zdir + '\\B*')
	s=start
	cal = []
	if len(acfn)>0:
		if end == 0 : end = len(acfn)
		if end > s:
			all_raws = pd.DataFrame()
			z=Zonar()
			#get starting values for frequencies
			_,_,start,cal = z.read_one_dive(z.read_raw(acfn[1]))
			#z.add_depth()
		for a in range(s, end):
			all_raws = all_raws.append(z.read_one_dive(z.read_raw(acfn[a]))[0])
			all_raws['gn']= np.array(start['gn'][all_raws['beam']-1])
	else:
		start=[];all_raws=[]
	return start, all_raws, cal

def update_cal(cal, raws, **kwargs):
	#TS=None, Ncal=None, N=None, alpha=None,cspeed=None,TS0=[-53.47, -50.92],beam=[9.8,4], gn=[40,40],G0=[54,54]
	#update the cal values provided
	cal.update(kwargs)
	#update raw data
	if 'Gain' in kwargs:
		raws['Gain'] = np.array(cal['Gain'])[raws['beam']-1]
	if 'Noise'  in kwargs:
		raws['Noise'] = np.array(cal['Noise'])[raws['beam']-1]
	if 'CalNoise' in kwargs:
		raws['CalNoise'] = np.array(cal['CalNoise'])[raws['beam']-1]
	if 'sl' in kwargs:
		raws['sl'] = np.array(cal['sl'])[raws['beam']-1]
	if 'tau' in kwargs:
		raws['tau'] = np.array(cal['tau'])[raws['beam']-1]
	if 'beam_deg' in kwargs:
		raws['beam_deg'] = np.array(cal['beam_deg'])[raws['beam']-1]
	if 'beam_rad' in kwargs:
		raws['beam_rad'] = np.array(cal['beam_rad'])[raws['beam']-1]
	return cal, raws	
	
def get_data(ac_dir, rawfn, satfn, miss, TS0=None, NL=[1540,1530], TS=[87,82], \
			 Scal=0, Tcal=20, dcal=5, Gain = [54,54], calfn=None):
	
	##########################################################################
	# RAW2 DATA (GPS, ENV, ZOOG)
	##########################################################################
	
	#get physical data
	if len(rawfn) > 0:
		env, gps, zoog = raw2meta_extract(rawfn[0])
		gps_c = pd.DataFrame({'Time':np.append(gps.UTC_time_fix_start.values,gps.UTC_time_fix_end.values),
									   'Lon':np.append(gps.lon_start.values,gps.lon_end.values),
									   'Lat':np.append(gps.lat_start.values,gps.lat_end.values)}).sort_values(by=['Time'])
		gps_c['dist'] = haversine(gps_c.Lat.shift(), gps_c.Lon.shift(), gps_c.loc[0:, 'Lat'], gps_c.loc[0:, 'Lon'])
		gps_c.loc[0,'dist'] = 0
		gps_c.dist = gps_c.dist.cumsum()
		
		print(time.ctime() + ': Computing soundspeed and alpha for each CTD point')
		#compute alpha for env
		env['alpha200'] = env.apply(lambda x: absorption(f=200,T=x['temperature'],S=x['salinity'],D=x['pressure']),axis=1)/1000
		env['alpha1000'] = env.apply(lambda x: absorption(f=1000,T=x['temperature'],S=x['salinity'],D=x['pressure']),axis=1)/1000
		
		#compute c for env
		env['c'] = env.apply(lambda x: compute_c(x['pressure'],x['salinity'], x['temperature']), axis=1)
		
	else:
		print(time.ctime() + ': No RAW2 file found for ' + miss + ' ...skipping mission...')
		return
	
	###########################################################################
	# SAT DATA (header, gps_sat, engineering, profile, zoocam, zonar, misc)
	###########################################################################
	
	#get satellite data
	if len(satfn)>0:
		header, gps_sat,engineering, profile, zoocam, zonar, misc = read_sat(satfn[0])
	else:
		print(time.ctime() + ': No SAT file found for ' + miss + ' ...skipping mission...')
		return
	
	##########################################################################
	# GET ACOUSTICS DATA
	##########################################################################
	
	#get acoustics zonar data
	print(time.ctime() + ': Start reading Zonar data')
	start, all_raws, cal = read_all_zonar(ac_dir)
	
	#get acoustic sample and ping time
	all_raws['sample_time'] = pd.to_datetime(all_raws['dive_time']).values.astype(np.int64) +  ( all_raws['dt']/1000 * (all_raws['nscan'] + 1) + all_raws['blank'] + all_raws['tau'] * 1000 + all_raws['Ping'] * all_raws['tPing'])/1000 * 10**9
	all_raws['ping_time'] = pd.to_datetime(all_raws['dive_time']).values.astype(np.int64) +  (all_raws['Ping'] * all_raws['tPing'])/1000 * 10**9
	all_raws['ping_time'] = pd.to_datetime(all_raws['ping_time'])
	
	#get salinity, temperature, fluorescence and longitude / latitude for each acoustic ping to compute c
	print(time.ctime() + ': Merging environmental data to Zonar data')
	env = env.dropna(subset=['fluorescence'])
	all_raws['lon'] = np.interp(all_raws.ping_time.astype('int64'),pd.to_datetime(gps_c.Time).astype('int64'), gps_c.Lon)
	all_raws['lat'] = np.interp(all_raws.ping_time.astype('int64'),pd.to_datetime(gps_c.Time).astype('int64'), gps_c.Lat)
	all_raws['temperature'] = np.interp(all_raws.ping_time.astype('int64'),pd.to_datetime(env.time_of_measure).astype('int64'), env.temperature)
	all_raws['salinity'] = np.interp(all_raws.ping_time.astype('int64'),pd.to_datetime(env.time_of_measure).astype('int64'), env.salinity)
	all_raws['fluorescence'] = np.interp(all_raws.ping_time.astype('int64'),pd.to_datetime(env.time_of_measure).astype('int64'), env.fluorescence)
	#compute c for each ping
	all_raws['c'] = np.interp(all_raws.ping_time.astype('int64'),pd.to_datetime(env.time_of_measure).astype('int64'), env.c)
	
	#update depth of acoustic samples
	print(time.ctime() + ': Updating Zonar depth information based on new environmental info')
	z0 = ( all_raws['blank']  + all_raws['tau']*1000 / 2 ) * all_raws['c']/2/1000 #center of first scan
	dz = z0 + all_raws['nscan'] *  all_raws['c']/2/1000 * all_raws['dt'] * 0.001
	zb = all_raws['press'] / 10
	z = zb + dz * np.cos(17*np.pi/180) #constant is conversion to depth from distance from transducer due to angle of assent at 17deg
	all_raws['dz'] = dz
	all_raws['z'] = z
	
	
	#compute attenuation for each ping
	print(time.ctime() + ': Getting attenuation at each position')
	all_raws['alpha'] = 0
	all_raws.loc[all_raws['Frequency']==200,'alpha'] = np.interp(all_raws[all_raws['Frequency']==200].ping_time.astype('int64'),
														   pd.to_datetime(env.time_of_measure).astype('int64'), env.alpha200)
	all_raws.loc[all_raws['Frequency']==1000,'alpha'] = np.interp(all_raws[all_raws['Frequency']==1000].ping_time.astype('int64'),
														   pd.to_datetime(env.time_of_measure).astype('int64'), env.alpha1000)
	
	#update nomwl, k, a, psiD	 - norminal wavel length, wave number, active radius, equivalent beam angle
	all_raws['nomwl'] = all_raws['c'] / all_raws['Frequency']
	all_raws['k'] = 2 * np.pi / all_raws['nomwl']
	all_raws['a'] = 1.6 / ( all_raws['k'] * np.sin(all_raws['beam_rad']/2))#active area
	all_raws['psiD'] = 10 * np.log10( 5.78 / ( ( all_raws['k'] * all_raws['a'] ) ** 2))#equivalent beam angle in steradians
	
	########################################
	# CAL PREPARATION UPDATE
	#########################################
	# get the cal information from file if available
	if calfn is not None:
		print(time.ctime() + ': Reading acoustic cal values from file...')
		cv = get_cal_val(calfn, gps.UTC_time_fix_start.min())
		NL = cv[cv.index.isin(['N[200]','N[1000]'])].values
		TS = cv[cv.index.isin(['T[200]','T[1000]'])].values
		print('NL = ' + str(NL) + ' TS  = ' + str(TS))
	
	#update calibration values
	print(time.ctime() + ': Updating calibration data')
	#TS0 - the theoretical TS of the used sphere
	if TS0 is None: TS0 = [-53.47, -50.92]
	
	# get the source level approximation
	alpha = np.array([absorption(f=f, S=np.array(Scal), T=np.array(Tcal), D=np.array([dcal]))/1000 for f in start['freq']]).flatten()
	
	# SL = Source level, 40log10(z) + 2 alpha z = 2 x Transmission Loss, G0 = system gain
	SL = np.array(NL)/start['gn'] + np.array(cal['Gain']) + 40 * np.log10(dcal)\
		 + 2 * alpha * dcal 
	
	#get expected TS in counts
	TS_c = TS0 + SL - 40*np.log10(dcal) + 2 * alpha * dcal
	
	#get the TS cal Gain value
	Gcal = TS_c - np.array(TS) #dB re V
	
	#########################################
	# NOISE / PASSIVE DIVE 
	#########################################
	
	#find the listening dive and use mean as noise measurement
	def get_lNoise(all_raws, freq, start):
		db = all_raws[all_raws['Frequency']==freq]['Raw'] / start['gn'][start['freq']==freq].values
		db = pd.DataFrame(db)
		db['Dive'] = all_raws[all_raws['Frequency']==freq]['dive']
		db = db.reset_index()
		
		#find the listening dive (i.e. the dive with minimum mean)
		db['lin'] = 10**(db['Raw']/10)
		ldive = db[['Dive','lin']].groupby('Dive').mean().idxmin()
		#get noise level from listen dive
		N  = db[db['Dive'] == ldive[0]]['Raw'].min()#27#zCal['Noise'][b]
		return N
			
			
	N = np.array([get_lNoise(all_raws,f, start) for f in start['freq']])
	#linear forms of the calibration and noise values
	n = 10**(np.array(N) / 10)
	
	
	c = all_raws['c'].mean()
	
	tau = np.array(start['pulse'])
	alpha_mean = all_raws.groupby('Frequency')['alpha'].mean().values	
	
	#update cal values
	cal, all_raws = update_cal(cal, all_raws, 
			Gain = Gain,Noise=N,TScal = Gcal.values, tau = tau, cspeed=c,
			TS_Gain = Gcal.values, CalNoise = np.array(NL)/40, sl = SL.values,
			alpha_cal=alpha,alpha=alpha_mean)
	#d0 = db[psel] 
	
	#create selection criteria for frequencies
	sel = all_raws['Frequency']
	sel = sel.replace(1000,1)
	sel = sel.replace(200,0)
	
	#create dB of raw value
	d0 = all_raws.Raw / start['gn'][sel].values
	d0 = 10**(d0/10)
	
	#compute SNR
	SNR = (( d0 - n[sel])/n[sel])
	SNR[SNR<0.1] = 0.1
	SNR = 10 * np.log10(SNR)
	
	#add sl, TScal and tau to all_raws	
	all_raws.loc[:,'sl'] = SL.values[sel]
	all_raws.loc[:,'TScal'] = Gcal.values[sel]
	all_raws.loc[:,'tau'] = tau[sel]
	
	#remove values with SNR < 3
	all_raws.loc[:,'SNR']=SNR
	
	dz = all_raws['nscan'].values *( all_raws['c'] /2/1000 )* 200/1000
	z0 = (all_raws['blank'] + (all_raws['tau'])/2)*(c/2/1000)
	d = z0 + dz
	
	#linearize raw
	all_raws.loc[:,'dB'] = 10*np.log10(d0)# - n[sel])

	#all_raws = all_raws[all_raws['SNR']>3]
	
	print(time.ctime() + ': Generating Sv')
	all_raws['Sv'] = all_raws['dB'].values - all_raws['sl'].values - ( 10 *np.log10( all_raws['c'].values * all_raws['tau'].values/1000  / 2  ) ) - all_raws['psiD'].values + ( 20 * np.log10(d)) + (2 *all_raws['alpha'].values * d)+all_raws['TScal'].values
	#Sv = d1.flatten() - SL - ( 10 *np.log10( c[SNR>3] * tau  / 2  ) ) - PSI + ( 20 * np.log10(d[SNR > 3])) + (2 *alpha[SNR>3] * d[SNR > 3])+G
			
	#get TS
	print(time.ctime() + ': Generating TS')		
	all_raws['TS'] = all_raws['dB'] - all_raws['sl'] + ( 40 * np.log10(d)) + (2 * all_raws['alpha'] * d)+all_raws['TScal']
	
	return  env, gps, zoog, all_raws, cal, miss, header, gps_sat,engineering, profile, zoocam, zonar, misc


def get_mission(mdir, TS0=None, NL=[1540,1530], TS=[76,81], Scal=0, Tcal=0, dcal=5, Gain=[54,54], **kwargs):
	ac_dir, rawfn,satfn,miss = get_paths(mdir)
	env, gps, zoog, all_raws, cal, miss, header, gps_sat,engineering, profile, zoocam, zonar, misc = get_data(ac_dir, rawfn, satfn, miss, TS0 = TS0, NL = NL, TS = TS, Scal = Scal, Tcal = Tcal, dcal = dcal, Gain = Gain)
	return env, gps, zoog, all_raws, cal, miss, header, gps_sat,engineering, profile, zoocam, zonar, misc


def mission2csv(mdir, outfn=None, TS0=None, NL=[1540,1530], TS=[76,81], Scal=0, Tcal=0, dcal=5, Gain=[54,54], **kwargs):
	ac_dir, rawfn,satfn,miss = get_paths(mdir)
	raw = get_data(ac_dir, rawfn, satfn, miss, TS0 = TS0, NL = NL, TS = TS, Scal = Scal, Tcal = Tcal, dcal = dcal, Gain = Gain)[0]
	out = raw[['dive','beam','Frequency','Ping','nBurst', 'nscan','z','dz','ping_time','lon','lat','c','alpha','temperature','salinity','fluorescence','SNR','Sv','TS', 'dB']]
	out.columns = ['Dive','Beam','Frequency','Ping','Burst','nSample','Depth', 'Range', 'Time','Longitude','Latitude','c_speed','alpha', 'temperature', 'salinity','fluorescence','SNR','Sv','TS', 'dB']
	
	if outfn is None:
		outfn = miss + '.csv.gz'
	
	print(time.ctime() + ': Writing to ' + str(outfn))
		
	out.to_csv(outfn, index = False, compression='gzip')
	

def get_AllRoiCounts(rdir='Z:\\Zooglider_end_to_end_results'):
	ldirs = glob.glob(rdir + '/20*/')
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
				outfn = miss + '.csv.gz'
			print(time.ctime() + ': Writing ' + outfn)
			cdat.to_csv(outfn, compression = 'gzip')
		else:
			print(time.ctime() + ': No converted_to_csv folder found for ' + ldir) 
			continue
	return cdat

def get_missions(mdir):
	return [i for i in next(os.walk(mdir))[1] if '20' in i]

def missions_to_csv(mdir):
	missions = get_missions(mdir)
	missions =[mdir + mi for mi in missions[2:]]
	for miss in missions:
		mission2csv(miss)
		
def miss2nc(filename, mdir, TS0=None, NL=[1540,1530], TS=[76,81], Scal=0, Tcal=0, \
			dcal=5, Gain=[54,54]):
	env, gps, zoog, all_raws, cal, miss, header, gps_sat,engineering, profile, \
		zoocam, zonar, misc = get_mission(mdir, TS0=TS0, NL=NL, TS=TS, Scal=Scal, \
									Tcal=Tcal, dcal=dcal, Gain=Gain)
	generate_nc(filename, env, gps, zoog, all_raws, miss, header,gps_sat,\
			 engineering,profile, zoocam, zonar, misc, cal)
	
#update for new nc funcion and add option for csv
def missions_to_nc(mdir, missions=0, outdir='', force=False, calfn = None, csv = True):
	if missions == 0: 
		missions = get_missions(mdir)
		missions = missions[5:]
	for miss in missions:
		print(time.ctime() + ': Processing ' + miss )
		outfn = outdir + miss + '.nc'
		if os.path.isfile(outfn) and force==False:
			print('NetCDF for ', miss,' already exists...skipping...')
			continue
		else:
			#get raw_dirs
			ac_dir, rawfn,satfn,miss = get_paths(mdir+miss)

			if len(ac_dir) == 0 or len(rawfn) == 0 or len(satfn) == 0 or len(miss) == 0:
				#print(time.ctime() + ': No RAW DATA FOLDER found for ' + miss+ ' ...skipping mission...')
				continue
			
			#get all data
			env, gps, zoog, all_raws, cal, miss, header, gps_sat,engineering, \
				profile, zoocam, zonar, misc = get_data(ac_dir, rawfn, satfn, \
											miss, TS0=None, Scal=0, Tcal=20, \
												dcal=5, calfn=calfn)
			#generate netCDF file
			generate_nc(filename=outfn,env=env, gps=gps, zoog=zoog, \
			   all_raws=all_raws, miss=miss,header=header,gps_sat=gps_sat,\
				   engineering=engineering,profile=profile, zoocam=zoocam, \
					   zonar=zonar, misc=misc, cal=cal)
			
			#write csv
			if csv == True:
				out = all_raws[['dive','beam','Frequency','Ping','nBurst', 'nscan','z','dz','ping_time','lon','lat','c','alpha','temperature','salinity','fluorescence','SNR','Sv','TS', 'dB']]
				out.columns = ['Dive','Beam','Frequency','Ping','Burst','nSample','Depth', 'Range', 'Time','Longitude','Latitude','c_speed','alpha', 'temperature', 'salinity','fluorescence','SNR','Sv','TS', 'dB']
				outfn = outdir + miss + '.csv.gz'
				print(time.ctime() + ': Writing to ' + str(outfn))
				out.to_csv(outfn, index = False, compression='gzip')
		