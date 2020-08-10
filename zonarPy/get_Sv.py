# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 08:49:47 2020

@author: sveng
"""

import xarray as xr
import numpy as np
import pandas as pd
from .help_fun import haversine, compute_c, absorption
import matplotlib.pyplot as plt
import time
import os
from astral import Astral


def get_Sv(miss,beam=1,ping=1,
		   add_daytime=True,
		   rem_bot=True,
		   bthresh=-50,
		   pitchfix=False,
		   pitch_plot=False, 
		   fix_c = False, 
		   fix_alpha=False,
		   fix_sl=False):
	"""
	Transform netCDF to Sv and TS with environmental data

	Parameters
	----------
	miss : path
		Path to netCDF file of the mission.
	beam : integer, optional
		For which beam the data is extracted, 1 for 200 kHz, 2 for 1000 kHz, \
			0 for 200 & 1000 kHz. The default is 1.
	ping : integer, optional
		For which ping the data will be extracted with 4 pings available \
			(1-4), 0 extracts data fro all pings. The default is 1.
	add_daytime: boolean
		Add daytime information to the output if True. The default is True.
	rem_bot: boolean
		Removes bottom below zonar detected bottom if True. 
		Effectively removes values above bthresh, the bottom threshold 
		within a 5 m range of the detected bottom. The default is True.
	bthresh: float
		Threshold value for the bottom detection, values above this threshold, 
		within a 5 m range of the detected bottom will be removed if 
		rem_bot is True. The default is -50, ignored if rem_bot is False.
	pitch_plot : boolean, optional
		True of False if a plot with the recorded pitch data should be plotted.\
			 The default is False.

	Returns
	-------
	Sv_out : pandas dataFrame
		A dataframe containing the Acoustic, Environmental and calibration \
			information.

	"""
	Sv_out = pd.DataFrame()
	beams = [1,2] if beam == 0 else [beam]
	print(time.ctime() + ': Starting...')
	for beam in beams:
		print(time.ctime() + ': Collecting environment for beam ' + str(beam))
		
		#TScal = TScal if type(TScal) is int else TScal[beam-1]
		
		
		#get environmental data - depth/pressure dBar, Dive, fluorescence, salinity, StartTime, temperature, time
		env = xr.open_dataset(miss, group='Environment').to_dataframe()
		env = env.dropna()
		env=env.reset_index()
		env = env.sort_values(by='time')
		
		#get info on active or passive
		aop = xr.open_dataset(miss, group='Sat/zonar/Zonar_beam').to_dataframe()
		aop['active'] = True
		aop.loc[aop['pulse'] == '0','active'] = False
		aop = aop.set_index(['dive#','f1'])
		aop = aop['active']
	
		#get calibraiton coefficients for CTD/engineering data
		calco = xr.open_dataset(miss, group='Sat/header/calibration').to_dataframe()
		
		#get GPS
		print(time.ctime() + ' Collecting GPS for beam ' + str(beam))
		gps = xr.open_dataset(miss, group='GPS').to_dataframe()
		gps = gps.reset_index()
		gps['Dive'] += 1
		
		gps_c = pd.DataFrame({'Time':np.append(gps.Time_start.values,gps.Time_end.values),
									   'Lon':np.append(gps.Lon_start.values,gps.Lon_end.values),
									   'Lat':np.append(gps.Lat_start.values,gps.Lat_end.values)}).sort_values(by=['Time'])
		gps_c['dist'] = haversine(gps_c.Lat.shift(), gps_c.Lon.shift(), gps_c.loc[0:, 'Lat'], gps_c.loc[0:, 'Lon'])
		gps_c.loc[0,'dist'] = 0
		gps_c.dist = gps_c.dist.cumsum()
		
				
		#get acoustic data of the selected beam
		zB1 = xr.open_dataset(miss, group='Zonar/Beam_' + str(beam)).to_dataframe()
		zB1 = zB1.dropna()
		zCal = xr.open_dataset(miss, group='Zonar/Calibration').to_dataframe()
		f = zCal['Frequency'][beam-1]
		
		#get lat lon and dist for env
		print(time.ctime() + ' Interpolating environment for beam ' + str(beam))
		env['lon'] = np.interp(pd.to_datetime(env['time']).astype('int64'),pd.to_datetime(gps_c.Time).astype('int64'), gps_c.Lon)
		env['lat'] = np.interp(pd.to_datetime(env['time']).astype('int64'),pd.to_datetime(gps_c.Time).astype('int64'), gps_c.Lat)
		env['dist'] = np.interp(pd.to_datetime(env['time']).astype('int64'),pd.to_datetime(gps_c.Time).astype('int64'), gps_c.dist)
		env['alpha'] = env.apply(lambda x: absorption(f=f,T=x['temperature'],S=x['salinity'],D=x['Depth']),axis=1)
		
		#engineering
		eng = xr.open_dataset(miss, group='Sat/engineering/Engineering_TS').to_dataframe()
		eng['pressure'] = eng['pressure_counts'] * float(calco[calco['Code'] == 'CP12']['gain']) + float(calco[calco['Code'] == 'CP12']['offset'])
		
		#env values for acoustics
		zB1['lon'] = np.interp(zB1.Time.astype('int64'),pd.to_datetime(gps_c.Time).astype('int64'), gps_c.Lon)
		zB1['lat'] = np.interp(zB1.Time.astype('int64'),pd.to_datetime(gps_c.Time).astype('int64'), gps_c.Lat)
		zB1['dist2D'] = np.interp(zB1.Time.astype('int64'),pd.to_datetime(gps_c.Time).astype('int64'), gps_c.Lon)
		zB1['temp'] = np.interp(zB1.Time.astype('int64'), pd.to_datetime(env['time']).astype('int64'), env.temperature)
		zB1['sal'] = np.interp(zB1.Time.astype('int64'),pd.to_datetime(env['time']).astype('int64'), env.salinity)
		zB1['fluo'] = np.interp(zB1.Time.astype('int64'),pd.to_datetime(env['time']).astype('int64'), env.fluorescence)
		zB1['fluo'] = np.interp(zB1.Time.astype('int64'),pd.to_datetime(env['time']).astype('int64'), env.fluorescence)
		zB1['depth_glider'] = np.interp(zB1.Time.astype('int64'),pd.to_datetime(env['time']).astype('int64'), env.Depth)
		
		if fix_alpha == False:
			zB1['alpha'] = np.interp(zB1.Time.astype('int64'),pd.to_datetime(env['time']).astype('int64'), env.alpha)/1000
		else:
			zB1['alpha'] = fix_alpha[beam-1]
			
		if pitchfix == False:
			for d in zB1['Dive'].unique():
				if eng[(eng['dive#'] == d)].shape[0] != 0:
					zB1.loc[zB1['Dive']==d,'Pitch'] = \
						pd.DataFrame(np.interp(zB1[zB1['Dive']==d][['depth_glider']],
							 eng[(eng['dive#'] == d) & (eng['pitch10deg']>0)]['pressure'],
							 eng[(eng['dive#']==d) & (eng['pitch10deg']>0)]['pitch10deg']/10)).values
				else:
					zB1.loc[zB1['Dive']==d,'Pitch'] = 17
		else:
			zB1['Pitch']  = pitchfix
			
		#get c for all acoustics values
		if fix_c == False:
			print(time.ctime() + ': Computing sound speed')
			zB1['c'] = compute_c(zB1.depth_glider, zB1.sal, zB1.temp, zB1.lat)
		else:
			print('Setting sound speed')
			zB1['c'] = fix_c
		#get attenuation in dB per m
		
		# ZOnar header info 
		#zHead = xr.open_dataset(miss, group='Sat/zonar/Zonar_beam').to_dataframe()

		z0 = (zCal['blank'][beam-1] + zCal['tau'][beam-1] / 2)  * zB1['c'].values/2/1000 #center of first scan
		dz = z0 + zB1.index.get_level_values('nScan').values *  zB1['c'].values/2/1000 * zCal['dt'][beam-1] * 0.001
		z = zB1["depth_glider"]+dz * np.cos(zB1['Pitch'].values*np.pi/180) #constant is conversion to depth from distance from transducer due to angle of assent at 17deg
		
		#select the selected ping
		#transform into dB re V, divide by counts per dB
		
		if ping == 0:
			pings = [1,2,3,4]
		else:
			pings = [ping]
		psels = ['Ping' +str(i) for i in pings]
	
		
		for psel in psels:
			print(time.ctime() + ': Processing Beam ' + str(beam) + ' for ' + psel)
			db = zB1[psel] / 40
			db = pd.DataFrame(db)
			db['Dive'] = zB1.Dive
			db = db.reset_index()
			
			#find the listening dive (i.e. the dive with minimum mean)
			db['lin'] = 10**(db[psel]/10)
			ldive = db[['Dive','lin']].groupby('Dive').mean().idxmin()
			
			db = db.set_index(['Dive','Burst','nScan'])
			
			
			#get noise level from listen dive
			N  = db.iloc[db.index.get_level_values('Dive') == ldive[0]][psel].min()#27#zCal['Noise'][b]
			#G0 = zCal['Gain'][beam-1] #system gain from calibration info
			#SL approximated as SL = EL_noise + 2TL + G0
			#SL = [120.54,122.409][beam-1] #####
			if fix_sl == False:
				SL = zCal['SoureLevel'][beam-1] #get source level from calibraiton info
				#SL = [120.54,120.289][beam-1] #####zCal['SoureLevel'][beam-1] #get source level from calibraiton info
			else:
				SL = fix_sl[beam-1]
			zB1['nomwl'] = zB1['c'] / zCal['Frequency'][beam-1]
			zB1['k'] = 2 * np.pi / zB1['nomwl']
			zB1['a'] = 1.6 / ( zB1['k'] * np.sin(zCal['beam_rad'][beam-1]/2))#active area
	
			zB1['PSI'] = 10 * np.log10( 5.78 / ( ( zB1['k'] * zB1['a'] ) ** 2))#equivalent beam angle in steradians
			PSI = zB1['PSI'].mean()
			
			Gcal = zCal.Gain_TS[beam-1]
			
			#linear forms of the calibration and noise values
			n = 10**(N/10)
			#g = 10**(G0/10)
			c = zB1['c'].values
			tau = zCal['tau'][beam-1]/1e3
			#bdms = 1
			#dz = zB1.reset_index()['nScan'].values *( c /2/1000 )* 200/1000
			#z0 = (bdms + (tau*1000)/2)*(c/2/1000)
			d = dz
			
			#alpha = zB1['alpha200'] #/ 1000
			alpha = zB1['alpha'] 
			G= Gcal
			
			d0 = db[psel] 
			d0 = 10**(d0/10)
			SNR = (( d0.values - n)/n)
			SNR[SNR<0.1] = 0.1
			SNR = 10 * np.log10(SNR)
			db['SNR'] = SNR
			
			#remove values with SNR < 3
			d0 = d0.values[SNR > 3]
			d1 = 10*np.log10(d0 - n) 
			
			Sv = d1.flatten() - SL - ( 10 *np.log10( c[SNR>3] * tau  / 2  ) ) - PSI + ( 20 * np.log10(d[SNR > 3])) + (2 *alpha[SNR>3] * d[SNR > 3])+G
			#Sv = d1.flatten() - SL - ( 10 *np.log10( c[SNR>3] * tau  / 2  ) ) - PSI + ( 20 * np.log10(d[SNR > 3])) + (2 *alpha[SNR>3] * d[SNR > 3])+G
			
			Sv = pd.DataFrame({'Sv': Sv})
			Sv['Depth'] = z[SNR>3]
			Sv['frequency'] = f
			Sv['beam'] = beam
			Sv['ping'] = psel
			Sv['Dive'] = zB1['Dive'][SNR>3]
			Sv['TS'] = d1.flatten() - SL + ( 40 * np.log10(d[SNR > 3])) + (2 *alpha[SNR>3] * d[SNR > 3])+G
			Sv['SNR'] = SNR[SNR>3]
			Sv['alpha'] = alpha[SNR>3]
			Sv['c'] = c[SNR>3]
			Sv['dz'] = dz[SNR>3]
			Sv['Sal'] = zB1['sal'][SNR>3]
			Sv['Fluo'] = zB1['fluo'][SNR>3]
			Sv['Temp'] = zB1['temp'][SNR>3]
			Sv['GDepth'] = zB1['depth_glider'][SNR>3]
			Sv['lon'] = zB1['lon'][SNR>3]
			Sv['lat'] = zB1['lat'][SNR>3]
			Sv['Time'] = zB1['Time'][SNR>3]
			Sv['dist2D'] = zB1['dist2D'][SNR>3]
			
			#get daytime 
			print(time.ctime() + ': Getting times of the day')
			a = Astral()
			a.solar_depression=12
			divesub = env.groupby('Dive')[['time','lon','lat']].mean()
			dtime=pd.to_datetime(divesub.time.values)
			dtime = dtime.tz_localize('UTC')
			sunpos = pd.DataFrame([a.sun_utc(dtime.mean().normalize(), \
									divesub.lat.values.mean(), 
									divesub.lon.values.mean())]).transpose()
			sunpos = sunpos[0].dt.tz_convert('America/Los_Angeles')[[0,1,3,4]]
			sunpos['prenight'] = sunpos['dawn'].normalize()+pd.Timedelta(value=1,unit='ms')
			sunpos['night']= sunpos['dawn'].normalize()+pd.Timedelta(value=23.99999,unit='hours')
			sunpos = sunpos.sort_values()
			cuts = sunpos.dt.hour.values*60+sunpos.dt.minute.values
			
			divesub['time'] = pd.to_datetime(divesub.time).dt.tz_localize('UTC')
			divesub['local'] = divesub.time.dt.tz_convert('America/Los_Angeles')
			
			divesub['daytime'] = pd.cut(divesub.local.dt.hour * 60 + \
							 divesub.local.dt.minute,\
								 cuts,\
									 labels=['night','dawn','day','dusk','nnight']).astype(str)
			divesub.loc[divesub['daytime']=='nnight','daytime'] = 'night'
			
			env = env.merge(divesub['daytime'].reset_index())
			
			#add bottom depth detected
			print(time.ctime() + ': Getting Bottom Depth')
			bd = xr.open_dataset(miss, group='Sat/zonar/z').to_dataframe()
			bd = bd.loc[bd.Code == 'ZT01']
			bd['z'] = bd['z']/10
			bd = bd.rename(columns={"Dive#": "Dive"})
			bd = bd.groupby('Dive').max().z+20
			Sv['BDepth'] = Sv.Dive.map(bd)
			
			#remove below bottom data
			Sv = Sv.loc[(Sv.BDepth > (Sv.Depth-5)) & (Sv.Sv < -50)]			
			#add daytime information
			Sv['daytime'] = Sv.Dive.map(divesub['daytime'])
			
			print(time.ctime() + ': Combining Data')
			Sv_out = pd.concat([Sv_out,Sv])
			
			

	return Sv_out,env,aop

'''
def ncToPiv2D(miss, outdir = None,
			  f = 0, p = 0, 
			  v = ['sv'], 
			  r = 1, 
			  rem_passive=True,
			  dzmin=1,dzmax=20,
			  fix_c=False, 
			  fix_alpha=False,
			  fix_sl=False,
			  pitchfix = False):
	"""
	Transforms netcdf data into 2D pivot table

	Parameters
	----------
	miss : path
		Path and filename of netcdf file. The default is 0.
	outdir : path
		Path to which the csv files will be stored. If outdir is None, the \
			path of the input file will be selected. The default is None.
	p : integer, optional
		Ping selection, between 1 and 4 for one ping or 0 for all pings. \
			The default is 0.
	v : list, optional
		List of variables to be extracted sv for Sv, ts for TS, Temp for \
			Temperature, Sal for Salinity, Fluo for FLuorescence, c for \
				sound speed, alpha for acoustic absorption coefficient, \
					lon for longitude, lat for latitude, time for time. 
					The default is ['sv'].
	r : numeric, optional
		Vertical resolution to which the data will be summarised. The default is 1.
	dzmin : numeric
		Minimum distance from transducer to use the data [m]. The default is 1.
	dzmax : numeric
		Maximum distance from the transducer to use [m]. The default is 20.

	Returns
	-------
	None.

	"""
	#transform frequency into beam value
	if f == 200:
		b = 1
	elif f == 1000:
		b = 2
	else:
		b = 0
		
	#check outdir
	if outdir is None:
		outdir = os.path.dirname(miss)
	
	#Translate data from netcdf into Sv, TS, etc.
	Sv,env,aop = get_Sv(miss, beam=b, ping=p,
			 fix_c=fix_c,fix_alpha=fix_alpha,fix_sl=fix_sl,pitchfix=pitchfix)
	
	if f == 0:
		s200 = Sv.loc[Sv.frequency==200][Sv.loc[Sv.frequency==200,'Dive'].isin(Sv.loc[Sv.frequency==1000,'Dive'].unique())]
		Sv = Sv.loc[Sv.frequency==1000][Sv.loc[Sv.frequency==1000,'Dive'].isin(Sv.loc[Sv.frequency==200,'Dive'].unique())]
		Sv = Sv.append(s200)
	Sv.loc[Sv.Dive.isin(aop.reset_index().loc[aop.reset_index().active==False,'dive#'].unique()),'Sv'] = np.nan
	
	
	
	#Sv = Sv.reset_index()
	Svdz = Sv[Sv['dz'] >= dzmin]
	Svdz = Svdz[Svdz['dz'] <= dzmax]
	
	#linearize Sv
	Svdz['sv'] = 10**(Svdz['Sv']/10)
	Svdz['ts'] = 10**(Svdz['TS']/10)
	
	Svdz = Svdz.reset_index().set_index([Svdz['Dive'],Svdz.reset_index()['Burst'],Svdz.reset_index()['Depth']])
	Svdz['Depth_r'] = np.round(Svdz['Depth']/r)*r
	
	Svdz['time'] = pd.to_datetime(Svdz.Time).to_numpy().astype(float)
	
	sub = Svdz[['sv','ts','frequency','Depth_r','lon','lat','time','Sal','Fluo','Temp','c']].groupby(['Dive','Depth_r','frequency']).mean()
	
	if f == 0:
		freqs = [200,1000]
	else:
		freqs = [f]
	for svar in v:
		if svar in ['sv', 'ts']:
			for f in freqs:
				print(time.ctime() + ': Processing ' + svar + ' @ ' + str(f) + 'kHz')
				fnout = os.path.basename(miss)[:-3] + '_' + svar + '_' + str(f) + 'kHz_vres_' + str(r) + 'm_per_dive_rb_' + str(dzmin) + '_' + str(dzmax) + '.csv'
				sfpiv = 10 * np.log10(sub.query('frequency=='+str(f)).pivot_table(values=svar,index='Depth_r',columns='Dive'))
				print(time.ctime() + ': Writing ' + fnout)
				sfpiv.to_csv(outdir + fnout)

		else:
			print(time.ctime() + ': Processing ' + svar )
			fnout = os.path.basename(miss)[:-3] + '_' + svar + '_vres_' + str(r) + 'm_per_dive_rb_' + str(dzmin) + '_' + str(dzmax) + '.csv'
			sfpiv = (sub.query('frequency=='+str(freqs[0])).pivot_table(values=svar,index='Depth_r',columns='Dive'))
			print(time.ctime() + ': Writing ' + fnout)
			sfpiv.to_csv(outdir + fnout)
	if len(freqs)>1:
		print(time.ctime() + ': Generating Sv Delta')
				
		Sv200 = pd.read_csv(outdir + os.path.basename(miss)[:-3] + \
					  '_' + v[0] + '_' + str(200) + 'kHz_vres_' + str(r) + \
						  'm_per_dive_rb_' + str(dzmin) + '_' + str(dzmax) + '.csv')
			
		Sv1000 = pd.read_csv(outdir + os.path.basename(miss)[:-3] + \
					   '_' + v[0] + '_' + str(1000) + 'kHz_vres_' + str(r) + \
						   'm_per_dive_rb_' + str(dzmin) + '_' + str(dzmax) + '.csv')
		SvDelta = Sv200 - Sv1000
		SvDelta['Depth_r'] = Sv200.Depth_r
		SvDelta = SvDelta.set_index('Depth_r')
		fnout = os.path.basename(miss)[:-3] + '_' + svar + '_' + \
			'Delta_vres_' + str(r) + 'm_per_dive_rb_' + str(dzmin) + '_' + str(dzmax) + '.csv'
		SvDelta.to_csv(outdir + fnout)
	
'''
def get_active_dives(miss):
	#get info on active or passive
	aop = xr.open_dataset(miss, group='Sat/zonar/Zonar_beam').to_dataframe()
	aop['active'] = True
	aop.loc[aop['pulse'] == '0','active'] = False
	aop = aop.set_index(['dive#','f1'])
	aop = aop['active']
	return aop
