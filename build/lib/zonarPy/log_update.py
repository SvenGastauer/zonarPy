# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 19:53:14 2020

@author: sveng
"""

import requests
import pandas as pd
import time
from read_sat import read_sat
import numpy as np
from astral import Astral
from datetime import datetime

def getSatFile(sn, outfn):
	link = "https://cce_lter:man0hman@spray.ucsd.edu/data/active/00"+ str(sn) + ".sat"
	print(time.ctime() + ': Getting sat data...' )

	f = requests.get(link)
	#print(f.text)
	f.raise_for_status() # ensure we notice bad responses
	print(time.ctime() + ': Writing sat data...' )

	file = open(outfn + ".sat", "w")
	file.write(f.text)
	file.close()


def sat2log(fn, outfn=None, xlsx=True):
	header, gps, engineering, profile, zoocam, zonar, misc = read_sat(fn)
	
	print(time.ctime() + ': Gathering info from sat file...' )
	
	gps['GPS_info'] = gps['GPS_info'].drop_duplicates()
	gps['GPS_info']  = gps['GPS_info'][gps['GPS_info']['mission_status'].astype(int).between(1,2, inclusive=True)]
	zoocam['Zoocam'] = zoocam['Zoocam'].drop_duplicates()
	engineering['Flight_Line'] = engineering['Flight_Line'].drop_duplicates()
	
	NDives = gps['GPS_info']['dive#'].to_numpy().astype(int).max()
	cols = gps['GPS_info'].columns.drop('Code')
	gps0 = gps['GPS_info'][cols].apply(pd.to_numeric, errors='coerce').groupby('dive#').mean()
	Latitude = gps0['Latitude'].to_numpy().astype(float)	   
	Longitude = gps0['Longitude'].to_numpy().astype(float)
	sel = gps0.index.values#np.arange(0,2*len(gps['GPS_info']['dive#'][gps['GPS_info']['dive#'].to_numpy().astype(int)>0].unique()),2)
	#Latitude = gps['GPS_info']['Latitude'][gps['GPS_info']['dive#'].to_numpy().astype(int) > 0 ].to_numpy().astype(float)[sel]
	#Longitude = gps['GPS_info']['Longitude'][gps['GPS_info']['dive#'].to_numpy().astype(int) > 0 ].to_numpy().astype(float)[sel]
	#Zcam
	
	Dive_no = gps['GPS_info']['dive#'][gps['GPS_info']['dive#'].to_numpy().astype(int)>0].unique().astype(int)
	Dir_no = np.zeros(len(Dive_no))
	Dir_no[0:len(zoocam['Zoocam']['outDir'])] = zoocam['Zoocam']['outDir'].to_numpy().astype(int)[1:]
	
	Zmax= engineering['Flight_Line']['Zmax'].to_numpy().astype(float)
	
	zc_active = np.zeros(len(Dive_no))
	zc_active[0:len(zoocam['Zoocam']['UseZcam'])] = zoocam['Zoocam']['UseZcam'].to_numpy().astype(int)[1:]
	
	Gb_rem_USB=zoocam['Zoocam']['Gb_rem_USB'].to_numpy().astype(int)
	Gb_rem_SD=zoocam['Zoocam']['Gb_rem_SD'].to_numpy().astype(int)
	
	ImageFile_Size=np.round(zoocam['Zoocam']['Mb_trans'].to_numpy().astype(int)/1000,1)
	
	zonar['Zonar_beam']['dive#'] = zonar['Zonar_beam']['dive#'].astype(int)
	zs_active = np.zeros(len(Dive_no))
	zs_active[zonar['Zonar_beam']['pulse'].groupby(zonar['Zonar_beam']['dive#']).sum().astype(int)[1:].values > 0] = 1
	
	#Start_YD = yd.groupby('Dive').min()
	#End_YD = yd.groupby('Dive').max()
	#Start_YD = misc['Email_mess']['Year-day'][misc['Email_mess']['dive#'].to_numpy().astype(int) > 0 ].to_numpy().astype(float)[np.arange(0,2*NDives,2)]
	#End_YD =  misc['Email_mess']['Year-day'][misc['Email_mess']['dive#'].to_numpy().astype(int) > 0 ].to_numpy().astype(float)[np.arange(1,2*NDives,2)]
	#Duration = np.array((End_YD - Start_YD ) * 60).astype(float)
	
	dtime = pd.to_datetime( gps['GPS_info']['month'] + gps['GPS_info']['day'] + gps['GPS_info']['year'] + gps['GPS_info']['time'],
						   format = '%b%d%Y%H:%M')
	dtime = dtime.dt.tz_localize('UTC')
	
	
	dtime_PST = dtime.dt.tz_convert('America/Los_Angeles')
	
	hours_UTC = dtime.dt.hour + dtime.dt.minute/60
	hours_PST = dtime_PST.dt.hour + dtime_PST.dt.minute/60
	
	yd = gps['GPS_info']['day'].to_numpy().astype(float) + hours_UTC/24
	Start_YD = np.array(yd[gps['GPS_info']['mission_status']=='1'])
	End_YD = np.array(yd[gps['GPS_info']['mission_status']=='2'])
	Duration = np.array((End_YD - Start_YD ) * 24*60).astype(float)
	
	Start_time_UTC = np.array(hours_UTC[gps['GPS_info']['mission_status']=='1'])
	End_time_UTC = np.array(hours_UTC[gps['GPS_info']['mission_status']=='2'])
	Start_time_PST = np.array(hours_PST[gps['GPS_info']['mission_status']=='1'])
	End_time_PST = np.array(hours_PST[gps['GPS_info']['mission_status']=='2'])
	Start_Date_PST = dtime_PST.dt.strftime('%d-%b-%Y')[gps['GPS_info']['mission_status']=='1']
	
	interval = np.insert(np.array([np.array(Start_time_PST)[i] - np.array(End_time_PST)[i-1] for i in range(1,len(End_time_PST))]),0,0)*60
	elapsed_time = dtime_PST[gps['GPS_info']['mission_status']=='2'] - dtime_PST[gps['GPS_info']['mission_status']=='1'].to_numpy()[0] 
	elapsed_time = np.array(elapsed_time.dt.total_seconds())/3600
	
	a = Astral()
	
	sunpos = np.array([a.solar_elevation(dtime.to_numpy()[sel[i]], Latitude[i], Longitude[i]) for i in range(len(Latitude))])
	twilight_index = np.where(sunpos > 0,0,np.NaN)
	twilight_index = np.where(sunpos < -12,2,twilight_index)
	twilight_index[np.where(np.isnan(twilight_index))] = twilight_index[np.where(np.isnan(twilight_index))[0]-1]+1
	while len(np.where(np.isnan(twilight_index))[0]) > 0:
		twilight_index[np.where(np.isnan(twilight_index))] = twilight_index[np.where(np.isnan(twilight_index))[0]-1]
	twilight_index = twilight_index.astype(int)
	
	zsr = 1/zoocam['Zoocam']['zCamMod'].to_numpy().astype(float) * 2 #every xth image sampled at 2 Hz 
	
	micol = pd.MultiIndex.from_arrays([['Dive Start Location','Dive Start Location',
										'Dive No','Dive No','Output','ZooCam',
						'ZooCam','USB storage','SD storage','Image File','Zonar',
						'UTC','UTC','Duration','Interval between dive','Start Time',
						'End Time','Start Time','Start Date','End Time','Elapsed Time',
						'Twilight Index','Zoocam'],
						['Latitude','Longitude','No. dives','Dive no.',
		  'Directory No.','Zmax','active','remaining','remaining','size (Gb)',
		  'active','Start_YD','End_YD','min','min','UTC','UTC','PST','PST','PST','hr',
		  '0=day,1=sunset,2=night,3=sunrise','Zoocam Sample Rate']])
										
	ll=pd.DataFrame(columns=micol,index=Dive_no)
	ll.loc[:,('Dive Start Location','Latitude')]=Latitude
	ll.loc[:,('Dive Start Location','Longitude')]=Longitude
	ll.loc[:,('Dive No','No. dives')] = NDives
	ll.loc[:,('Dive No','Dive no.')] = Dive_no
	ll.loc[:,('Output','Directory No.')] = Dir_no
	ll.loc[:,('ZooCam','Zmax')] = Zmax
	ll.loc[:,('ZooCam','active')] = zc_active[Dive_no-1]
	ll.loc[:,('USB storage','remaining')] = Gb_rem_USB[Dive_no-1]
	ll.loc[:,('SD storage','remaining')] = Gb_rem_SD[Dive_no-1]
	ll.loc[:,('Image File', 'size (Gb)')] = ImageFile_Size[Dive_no-1]
	ll.loc[:,('Zonar','active')] = zs_active
	ll.loc[:,('UTC','Start_YD')] = Start_YD
	ll.loc[:,('UTC','End_YD')] = End_YD
	ll.loc[:,('Duration','min')] = Duration
	ll.loc[:,('Interval between dive','min')] = interval
	ll.loc[:,('Start Time','UTC')] = Start_time_UTC
	ll.loc[:,('End Time','UTC')] = End_time_UTC
	ll.loc[:,('Start Date','PST')] = Start_Date_PST.values
	ll.loc[:,('Start Time','PST')] = Start_time_PST
	ll.loc[:,('End Time','PST')] = End_time_PST
	ll.loc[:,('Elapsed Time','hr')] = elapsed_time
	ll.loc[:,('Twilight Index','0=day,1=sunset,2=night,3=sunrise')] = twilight_index
	ll.loc[:,('Zoocam','Zoocam Sample Rate')] = zsr[Dive_no-1]
	
	ll.loc[ll.index[1]:,('Dive No','No. dives')] = ''
	
	if xlsx is True:
		print(time.ctime() + ': Writing log file...' + outfn + '.xlsx')
		ll.to_excel(outfn + '.xlsx')
	
	return ll
