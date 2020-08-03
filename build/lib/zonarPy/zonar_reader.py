# -*- coding: utf-8 -*-
'''
:author Sven Gastauer
:licence MIT
'''
#load packages
import numpy as np
import pandas as pd
import os
from datetime import datetime
from collections import defaultdict

#for plotting
import plotly
import plotly.graph_objs as go

#pickle for saving - obsolete
#import pickle as pkl

class Zonar():
	def read_raw(self, fname):
		"""
		reads in the raw acoustic data from the Zonar dive data

		Parameters
		----------
		fname : string
			Path to Zonar raw file

		Returns
		-------
		bfile - numpy array of the raw data

		"""
		try:
			print('\nReading file %s...'%os.path.split(fname)[1])
			bfile = np.fromfile(fname, dtype='uint8')
		except FileExistsError:
			print("No such file")
		return bfile
	
	def read_header(self, bfile, k=1):
		"""
		reads in the header data of the raw file.
				

		Parameters
		----------
		bfile : array
			Numpy array from raw zonar data as created by Zonar.read_raw()
		k : int, optional
			Starting position of the header in the file. The default is 1.

		Returns
		-------
		int
			PID - Packet ID for the FOLLOWING (!) packet
		k : int
			Ending position of the header in the file.
		head : pandas DataFrame
			Datatframe containing the header information.
			Header contains:
				Name  	Description 	Format
				---- 	----------- 	-------
				ver 	software version 	uint8 (1 byte)
				size 	header size 	uint8 (1 byte)
				pid 	Packet ID for the following packet 	uint8 (1 byte)
				nbytes 	Size of the following packet 	i4 (4 bytes)
				chksum 	Checksum for the following packet 	i4 (4 bytes)
		"""
		# define the header datatypes to get the information and the 
		# correct size
		dthead = np.dtype([('ver','uint8'),
					   ('size', 'uint8'),
					   ('pid', 'uint8'),
					   ('nbytes','i4'),
					   ('chksum','i4')])
		# reade the header information into pandas dataframe from the buffered
		# np array
		head = pd.DataFrame(np.frombuffer(bfile[k:k+dthead.itemsize],dthead))
		# define the positon in the file after the header
		k= k + dthead.itemsize
		return head['pid'], k, head
	
	def read_start(self, bfile, k):
		"""
		Reads the Dive Start packet (56 bytes of info).
		Contains the Configuration of the Zonar for this dive.

		Parameters
		----------
		bfile : array
			Numpy array from raw zonar data as created by Zonar.read_raw()
		k : int, optional
			Starting position of Dive info in the file. 
			This position is provided by the previously read in part.


		Returns
		-------
		k : int
			Ending position of the Dive configuration information in the file.
		start : pandas DataFrame
			Dive Starting information, includes information shared by both frequencies:
				Name 	Format 	Description
				---- 	------ 	-----------
				ver 	uint8 	packet version
				diveNo 	i2  	Dive number (from glider)
				sprayTm i4 	 	Glider Time [s] (from glider)
				zooTm 	i4 	 	Zonar Time [s]
				nPings 	uint8 	Number of pings in each Burst
				dtBurst uint8 	Number of seconds between bursts, if in continuous mode
				warmup 	i2 	 	time to warmup the electronics before transmit [ms]
				tbin 	i2 	 	time to average data into each bin [ms]
				barker 	uint8 	TRUE if use barker-code, else straight freq
				nBit 	uint8 	Number of bits in barker code
				
		freqstart : pandas DataFrame
			Frequency specific Dive Starting information:
				Name 	Format 	Description
				---- 	------ 	-----------
				freq 	i2  	nominal frequency [kHz]
				pulse 	i2  	transmit pulse duration [ms]
				blank 	i2  	time between end of transmit and the first scan [ms]
				dt 	i2  	period between scans [ms] e.g. dt=200 ms -> 5 kHz sample rate
				tScan 	i2  	duration to take scans [ms], i.e. tWait = tPing
				tPing 	i2  	time interval between pings [ms]
				nScan 	i2  	number of a/d scans taken, i.e. nScan = floor(1000 * tScan / dt)
				tWait 	i2  	time before next scan [ms], i.e. tWait = tPing - nScan * dt - blank - pulse
				nBin 	i2  	number of scans to average per bin, i.e. nBin = round(1000 * tBin / dt)
				gn 	 	i2  	number of a/d counts per 1 dB re V		

		"""
		csync   = int('aa', 16) # 0xAA = expected sync character
		sync	= bfile[k] #read the sync position
		k += 1
		if sync == csync:		
			#start
			dtstart = np.dtype([('ver', 'uint8'),
								('diveNo','i2'),
								('sprayTm', 'i4'),
								('zooTm','i4'),
								('nPings', 'uint8'),
								('dtBurst', 'uint8'),
								('warmup', 'i2'),
								('tbin', 'i2'),
								('barker', 'uint8'),
								('nBit', 'uint8')])
			start = pd.DataFrame(np.frombuffer(bfile[k:k+dtstart.itemsize],dtstart))
			k = k + dtstart.itemsize
				   
			dtfreq = np.dtype([('freq', 'i2'),
							   ('pulse', 'i2'),
							   ('blank', 'i2'),
							   ('dt', 'i2'),
							   ('tScan', 'i2'),
							   ('tPing', 'i2'),
							   ('nScan', 'i2'),
							   ('tWait', 'i2'),
							   ('nBin', 'i2'),
							   ('gn', 'i2')])
			freqstart = pd.DataFrame()
			for i in range(0,2):
				freqstart = freqstart.append(pd.DataFrame(np.frombuffer(bfile[k:k+dtfreq.itemsize],dtfreq)), ignore_index=True)
				k=k+dtfreq.itemsize
		else:
			k = -1
			start = 0
			freqstart = 0
		return k, start, freqstart
	
	def read_burst(self, bfile, k):
		"""
		reads one selected burst header and data

		Parameters
		----------
		bfile : array
			Numpy array from raw zonar data as created by Zonar.read_raw()
		k : int, optional
			Starting position of Dive info in the file. 
			This position is provided by the previously read in part.


		Returns
		-------
		k : int
			Ending position of the Dive configuration information in the file.
		
		burst : pandas DataFrame
			Burst header data (12 bytes of header data):
				
				Name 	Format 	Description
				---- 	------ 	-----------
				ver 	uint8 	Allows for future versions of the data packet, set to 0
				beam 	uint8 	Frequency ID with 1 = 200 kHz, 2= 1000 kHz, refer to dive start package for exact frequency
				nP 	 	uint8 	Number of pings in the burst
				nS 	 	i2 	 	Number of scans per Ping
				press 	i2 	 	pressure from glider [LSB, 0.1 dBar] with 1 LSB = 0.1 dBar
				zooTm 	i4 	 	time stamp: can be synched with glider time using dive_start info [s] (unix epoch time)
				
		rs : pandas DataFrame
			Burst data:
				2 bytes for each nP * nS. Sequential per ping (nS for ping number 1, nS for ping number 2, ....)
				
		"""
		csync = int('aa', 16)  # 0xAA = expected sync character
		sync  = bfile[k] #sync char found
		k = k + 1
		if csync == sync:
			#read burst header data
			dtburst = np.dtype([('ver', 'uint8'),
								('beam', 'uint8'),
								('nP', 'uint8'),
								('nS', 'i2'),
								('press', 'i2'),
								('zooTm', 'i4'),])  
			#get burst data
			burst = pd.DataFrame(np.frombuffer(bfile[k : k + dtburst.itemsize], dtburst))
			k = k + dtburst.itemsize
			nB = burst['nP'] * burst['nS'] #total scans over all pings
			rs = pd.DataFrame(np.frombuffer(bfile[k : int(k + nB[0] * 2)], dtype='i2', count=nB[0]))
			ping_lst = range(0,burst['nP'][0])
			rs['Ping'] = np.repeat(ping_lst, burst['nS'][0])
			k = k +  nB[0] * 2
		else:
			k = -1
			burst = 0
			rs = 0
		return k, burst, rs		
	
	def read_avg(self, bfile,k):
		"""
		reads the average packet data

		Parameters
		----------
		bfile : array
			Numpy array from raw zonar data as created by Zonar.read_raw()
		k : int, optional
			Starting position of average pacet data in the file. 
			This position is provided by the previously read in part.


		Returns
		-------
		k : int
			Ending position of the Dive configuration information in the file.
		
		burst : pandas DataFrame
			Burst header data (12 bytes of header data):
				
				Name 	Format 	Description
				---- 	------ 	-----------
				ver 	uint8 	Allows for future versions of the data packet, set to 0
				beam 	uint8 	Frequency ID with 1 = 200 kHz, 2= 1000 kHz, refer to dive start package for exact frequency
				nBin 	uint8 	Number of averaged bins = floor(total number of Scans / scans per bin)
				press 	i2 	 	pressure from glider [LSB, 0.1 dBar] with 1 LSB = 0.1 dBar
				zooTm 	i4 	 	time stamp: can be synched with glider time using dive_start info [s] (unix epoch time)
				
		rs : pandas DataFrame
			Burst data:
				2 bytes for each averaged return data
				
		
		"""
		csync = int('aa', 16)  # 0xAA = expected sync character
		sync  = bfile[k]; k += 1  #sync char found
		if sync == csync: # match!!! continue
					
			dtavg = np.dtype([('ver', 'uint8'),
					  ('beam', 'uint8'),
					  ('nBin', 'uint8'),
					  ('press','i2'),
					  ('zooTm','i4')])
			avg = pd.DataFrame(np.frombuffer(bfile[k : k + dtavg.itemsize], dtavg))
			k = k + dtavg.itemsize
			
			rs = np.frombuffer(bfile[k : int(k + avg['nBin'][0] * 2)], dtype='i2', count=avg['nBin'][0])
			k = k +  avg['nBin'][0] * 2 + 2
		else:
			k = -1
			avg = 0
			rs = 0
		return k, avg, rs
	
	def init_cal(self, **kwargs):
		"""
		intializes original calibration values
		!!!These need to be corrected for the latest calibration values!!!
		
		Parameters
		----------
		**kwargs : np array of size 2
		Any calibration values can be overwritten:
			Gain 	float 	calibration gain
			Noise 	float 	Noise
			CalNoise 	float 	calibration noise
			sl 	 	float 	Source level
			tau 	float  	transmit pulse duration [ms]
			beam_deg 	float 	beam angle in degrees
			alpha 	float 	anntenuation coefficient
			cspeed 	 	float 	sound speed in surrounding medium [m/s]
			gn 	 	int 	Number of a/d counts per dB re V

		Returns
		-------
		cal : defaultdict
			Gain 	float 	calibration gain
			Noise 	float 	Noise
			CalNoise 	float 	calibration noise
			sl 	 	float 	Source level
			tau 	float  	transmit pulse duration [ms]
			beam_deg 	float 	beam angle in degrees
			alpha 	float 	anntenuation coefficient
			cspeed 	 	float 	sound speed in surrounding medium [m/s]
			gn 	 	int 	Number of a/d counts per dB re V
			beam_rad 	float 	beam angle in radians

		"""
		cal = defaultdict(list)
		cal["Gain"] = [54, 54] #system gain
		cal["TS_Gain"] = [0,0] #calibraiton gain
		cal["Noise"] =  [27, 39] #noise
		cal["CalNoise"] = [31, 37] #calibration noise
		cal["sl"] = [103, 113] #source elvel
		cal["tau"] = [0.006, 0.006] #pulse duration [s]
		cal["beam_deg"] = [9.8, 4] #beam angle
		cal["alpha"] = [0.054, 0.38] #attenuation coefficient
		cal["cspeed"] = 1500 #sound speed in surrounding medium
		cal['gn'] = [40,40]
		cal.update(kwargs) #update with provided values
		
		cal["beam_rad"] = list(np.array(cal["beam_deg"]) * np.pi / 180) #convert beam angle from degrees to radians
		
		return cal
	
	def update_cal(self, raws, start, cal, **kwargs):
		"""
		updates calibration data in the raw output

		Parameters
		----------
		raws : pandas DataFrame
			raw output see read_one_dive() for details
		start : pandas DataFrame
			read_start output()
		cal : defaultdict
			init_cal output
		**kwargs : np array of size 2
		Any calibration values can be overwritten:
			Gain 	float 	calibration gain
			Noise 	float 	Noise
			CalNoise 	float 	calibration noise
			sl 	 	float 	Source level
			tau 	float  	transmit pulse duration [ms]
			beam_deg 	float 	beam angle in degrees
			alpha 	float 	anntenuation coefficient
			cspeed 	 	float 	sound speed in surrounding medium [m/s]
			gn 	 	int 	Number of a/d counts per dB re V

		Returns
		-------
		cal : defaultdict
			Calibration output, see cal_init() for details
		raws : pandas DataFrame
			raw data output see read_one_dive() for details

		"""
		
		#update the cal values provided
		cal.update(kwargs)
		
		#update raw data
		raws['Frequency'] = np.array(start['freq'])[raws['beam']-1]
		raws['Gain'] = np.array(cal['Gain'])[raws['beam']-1]
		raws['Noise'] = np.array(cal['Noise'])[raws['beam']-1]
		raws['CalNoise'] = np.array(cal['CalNoise'])[raws['beam']-1]
		raws['sl'] = np.array(cal['sl'])[raws['beam']-1]
		raws['tau'] = np.array(cal['tau'])[raws['beam']-1]
		raws['beam_deg'] = np.array(cal['beam_deg'])[raws['beam']-1]
		raws['beam_rad'] = np.array(cal['beam_rad'])[raws['beam']-1]
		raws['alpha'] = np.array(cal['alpha'])[raws['beam']-1]
		raws['c'] = cal['cspeed']
		
		raws['nomwl'] = raws['c'] / raws['Frequency']
		raws['k'] = 2 * np.pi / raws['nomwl']
		raws['a'] = 1.6 / ( raws['k'] * np.sin(raws['beam_rad']/2))#active area
		raws['psiD'] = 10 * np.log10( 5.78 / ( ( raws['k'] * raws['a'] ) ** 2))#equivalent beam angle in steradians
		
		#update depth with new cal data
		dz,raws = self.add_depth(start, raws)
		#update Sv with the new claibration data
		raws = self.get_Sv(start, raws)
		
		return cal, raws
		
		
	
	def read_one_dive(self, bfile, **kwargs):
		"""
		reads one full dive

		Parameters
		----------
		bfile : array
			Numpy array from raw zonar data as created by Zonar.read_raw()
		**kwargs : TYPE
			kwargs can be calibration values:
				Gain - calibration gain, defaults 54, 43
				Noise, noise defaults to 27, 39
				CalNoise, calibration noise defaults to 31, 37
				sl, source level, defaults to 103, 113
				tau, pulse duration, defaults to 0.006, 0.006
				beam_deg, beam angle defaults to 9.8, 4
				alpha, attenuation coefficient, defaults to 0.054, 0.38
				cspeed, sound speed in surrounding medium, defaults to 1500
				gn, number of a/d counts per dB re V, defaults to 40, 40

		Returns
		-------
		raws : pandas DataFrame
			Raw data including header and calibration information:
				
				Name 	Format 	Description
				---- 	------ 	-----------
				Raw 	uint8 	Raw count data
				Ping 	int 	Ping number
				dive 	i2  	Dive number (from glider)
				beam 	uint8 	Frequency ID with 1 = 200 kHz, 2= 1000 kHz, refer to dive start package for exact frequency
				nBurst 	uint8 	Number of burst
				press 	i2 	 	pressure from glider [LSB, 0.1 dBar] with 1 LSB = 0.1 dBar
				zooTm 	i4 	 	time stamp: can be synched with glider time using dive_start info [s] (unix epoch time)
				sprayTm i4 	 	Glider Time [s] (from glider)
				nP 	 	uint8 	Number of pings in the burst
				nS 		i2  	number of a/d scans taken, i.e. nScan = floor(1000 * tScan / dt)
				Frequency 	 	i2  	nominal frequency [kHz]
				Gain 	float 	calibration gain
				Noise 	float 	Noise
				CalNoise 	float 	calibration noise
				sl 	 	float 	Source level
				tau 	float  	transmit pulse duration [ms]
				beam_deg 	float 	beam angle in degrees
				beam_rad 	float 	beam angle in radians
				alpha 	float 	anntenuation coefficient
				c 	 	float 	sound speed in surrounding medium [m/s]
				nomwl 	float 	nominal wave length, lambda = c / Frequency [m]
				k 	 	float 	wave number, 2pi / lambda
				a 	 	float 	active area 1.6 * (k / sin(beam_rad/2))
				psiD 	float 	equivalent beam area 10 * log10 (5.78 / ((ka)**2))
				nScan 	i2  	number of a/d scans taken, i.e. nScan = floor(1000 * tScan / dt)
				dt 	 	i2  	period between scans [ms] e.g. dt=200 ms -> 5 kHz sample rate
				blank 	i2  	time between end of transmit and the first scan [ms]
				tScan 	i2  	duration to take scans [ms], i.e. tWait = tPing
				tPing 	i2  	time interval between pings [ms]
				tWait 	i2  	time before next scan [ms], i.e. tWait = tPing - nScan * dt - blank - pulse
				dz 	 	float 	distance from transducer with the center of first scan being z0 = blank + tau * 1000 / 2 * c / 2 / 1000 , dz becomes z0 + nscan * c / 2 / 1000 * dt * 0.001
				z 	 	float 	pressure [dBar] or depth [m] with the pressure in dBar zb = press/10, z becomes z = zb + dz * cos(17 * pi / 180), assuming a glider tilt angle of 17 degrees
				start_time 	date 	%Y-%m-%d %H:%M:%S of Glider time, sprayTm
				dive_time 	date 	%Y-%m-%d %H:%M:%S of zooTm, Zonar Time
				
		avgRaws : pandas DataFrame
			Averaged Data:
				0 	 	int 	chksum
				dive 	int 	Dive number
				beam 	uint8 	Frequency ID with 1 = 200 kHz, 2= 1000 kHz, refer to dive start package for exact frequency
				nAvg 	int 	Number of averaged bin
				nBin 	uint8 	Number of averaged bins = floor(total number of Scans / scans per bin)
				zooTm 	i4 	 	time stamp: can be synched with glider time using dive_start info [s] (unix epoch time)
				sprayTm i4 	 	Glider Time [s] (from glider)
				nP 	 	uint8 	Number of pings in the burst
				
		freqstarts : pandas DataFrame
			Frequency specific information:
				freq 	int 	Nominal Frequency [kHz]
				pulse 	int 	Pulse duration [ms]
				blank 	int 	blanking time [ms]
				dt 	 	int 	period between scans [ms] e.g. dt=200 ms -> 5 kHz sample rate
				tPing 	i2  	time interval between pings [ms]
				tScan 	i2  	duration to take scans [ms], i.e. tWait = tPing
				tWait 	i2  	time before next scan [ms], i.e. tWait = tPing - nScan * dt - blank - pulse
				nBin 	uint8 	Number of averaged bins = floor(total number of Scans / scans per bin)
				gn 	 	int 	Number of a/d counts per dB re V

		cal : defaultdict
			Calibration information:
				Gain 	float 	calibration gain
				Noise 	float 	Noise
				CalNoise 	float 	calibration noise
				sl 	 	float 	Source level
				tau 	float  	transmit pulse duration [ms]
				beam_deg 	float 	beam angle in degrees
				alpha 	float 	anntenuation coefficient
				cspeed 	 	float 	sound speed in surrounding medium [m/s]
				gn 	 	int 	Number of a/d counts per dB re V
				beam_rad 	float 	beam angle in radians
				
		"""
		print('Processing dive...')
		
		#initialise cal and update with input if any is provided
		cal = self.init_cal(**kwargs)
		
		#start timer
		t0 = datetime.now()
		#-- init values ----------------------
		id_miss  = int('c0', 16) #header id for start-mission (future use).
		id_start = int('c1', 16) #start-dive info (settings)
		id_end   = int('c2', 16) #end-dive info
		id_burst = int('c3', 16) #burst raw data
		id_avg   = int('c4', 16) #burst avg data
		id_EOF   = int('ff', 16) #end-of-file 
		
		#initialise empty dataframes
		#starts = pd.DataFrame()
		freqstarts = pd.DataFrame()
		bursts = pd.DataFrame()
		raws = pd.DataFrame()
		avgs = pd.DataFrame()
		avgRaws = pd.DataFrame()
		
		#set starting values
		k = 0 
		#k = 1
		id = 0
		nBurst = 0
		nAvg = 0
		nb = len(bfile)
		
		#loop through file
		while id > -1  and id < 255 and k < nb:  
			# id = -1 at EOF
			# added k<nb sep17
			# find next header in file
			id, k,head = self.read_header(bfile, k + 1 ) #= -1 if eof
			id = id[0]
			# updates k = index to the next packet
			if id == id_start:
				#print( '\nstart dive packet ')
				k, start, freqstart = self.read_start(bfile, k)
				#starts = start #starts.append(start)
				freqstarts = freqstart #freqstarts.append(freqstart)
				
			elif id == id_end:
				print('\nend dive packet ')
				#k = read_end( k )
				
			elif id ==  id_burst:
				#print('burst packet');
				nBurst = nBurst + 1
				k, burst, rs = self.read_burst(bfile, k)
				#k=k[0]
				bursts = bursts.append(burst)
				rs['dive'] = start['diveNo'][0]
				rs['beam'] = burst['beam'][0]
				rs['nBurst'] = nBurst
				rs['press'] = burst['press'][0]
				rs['zooTm'] = burst['zooTm'][0]
				rs['sprayTm'] = start['sprayTm'][0]
				rs['nP'] = burst['nP'][0]
				rs['nS'] = burst['nS'][0]
				raws = raws.append(rs)
				
			elif id ==  id_avg:
				#print('avg packet');
				nAvg += 1
				k, avg, rs = self.read_avg( bfile , k)
				avgs = avgs.append(avg)
				rs = pd.DataFrame(rs)
				rs['dive'] = start['diveNo'][0]
				rs['beam'] = avg['beam'][0]
				rs['nAvg'] = nAvg
				rs['press'] = avg['press'][0]
				rs['zooTm'] = avg['zooTm'][0]
				rs['sprayTm'] = start['sprayTm']
				rs['nP'] = avg['nBin'][0]
				
				avgRaws = avgRaws.append(rs)
				#kavg = kavg + 1
			elif id ==  id_EOF:
				print('EOF')
				#truncate_zonar(nBurst) #tidy up global
				#return
			else:
				string = 'Unknown PID %i'%int(id)
				print(string)
				k = k + head['nbytes']
				k = k[0]
		
		raws.rename(columns={0:'Raw'}, inplace=True)	   
		if len(raws)>0:
			raws['Frequency'] = np.array(freqstarts['freq'])[raws['beam']-1]
			raws['Gain'] = np.array(cal['Gain'])[raws['beam']-1]
			raws['Noise'] = np.array(cal['Noise'])[raws['beam']-1]
			raws['CalNoise'] = np.array(cal['CalNoise'])[raws['beam']-1]
			raws['sl'] = np.array(cal['sl'])[raws['beam']-1]
			raws['tau'] = np.array(cal['tau'])[raws['beam']-1]
			raws['beam_deg'] = np.array(cal['beam_deg'])[raws['beam']-1]
			raws['beam_rad'] = np.array(cal['beam_rad'])[raws['beam']-1]
			raws['alpha'] = np.array(cal['alpha'])[raws['beam']-1]
			raws['c'] = cal['cspeed']
			
			raws['nomwl'] = raws['c'] / raws['Frequency']
			raws['k'] = 2 * np.pi / raws['nomwl']
			raws['a'] = 1.6 / ( raws['k'] * np.sin(raws['beam_rad']/2))#active area
			raws['psiD'] = 10 * np.log10( 5.78 / ( ( raws['k'] * raws['a'] ) ** 2))#equivalent beam angle in steradians
			
			
			dz,raws = self.add_depth(freqstarts, raws)
			#raws = self.get_Sv(freqstarts, raws)
			
			#start time
			t1 = np.array(raws['sprayTm'])
			raws['start_time'] = [pd.Timestamp.utcfromtimestamp(x).strftime('%Y-%m-%d %H:%M:%S') for x in t1]
			
			#dive time
			t2 = t1 + np.array(raws['zooTm']) - np.array(raws['zooTm'])[0]
			raws['dive_time'] = [datetime.utcfromtimestamp(x).strftime('%Y-%m-%d %H:%M:%S') for x in t2]
				
		t1 = datetime.now()
		t2 = t1 - t0
	
		#remove depths greater than 20 from 1000 kHz
	
	
		print("Finished after " + str(np.round(divmod(t2.total_seconds(),3600)[1],2))+ " seconds.")
		return raws, avgRaws, freqstarts, cal 
	
	def add_depth(self, start, raws):
		"""
		add depth or pressure information to dataset

		Parameters
		----------
		start : pandas DataFrame
			read_start output()
		raws : pandas DataFrame
			raw output see read_one_dive() for details
		
		Returns
		-------
		dz : np array
			dz 	 	float 	distance from transducer with the center of first scan being z0 = blank + tau * 1000 / 2 * c / 2 / 1000 , dz becomes z0 + nscan * c / 2 / 1000 * dt * 0.001
		raws : pandas DataFrame
			raw output from read_one_dive() with updated nscan, dt, blan, tScan, tPing,tWiat, dz and z information

		"""
		#number of pings and samples to get the depth variable into the right shape
		nP = (raws['nP'].unique())
		nScan = np.array(start['nScan'])
		nP = np.repeat(nP, int(len(nScan) / len(nP)))
		 
		nscan = []
		for i in range(0,len(start)): 
			nscan = np.append(nscan, np.tile(np.arange(nScan[i]),int(nP[i])))
		   
		raws['nscan'] = np.tile(nscan, (int(len(raws) / len(nscan))))
		raws['dt'] = np.array(start['dt'])[raws['beam']-1]
		dz = raws['nscan']  * raws['c']/2/1000 * raws['dt'] * 0.001
		raws['blank'] = np.array(start['blank'])[raws['beam']-1]
		raws['tScan'] = np.array(start['tScan'])[raws['beam']-1]
		raws['tPing'] = np.array(start['tPing'])[raws['beam']-1]
		raws['tWait'] = np.array(start['tWait'])[raws['beam']-1]
		
		z0 = ( raws['blank']  + raws['tau'] * 1000 / 2 ) * raws['c']/2/1000 #center of first scan
		dz = z0 + raws['nscan'] *  raws['c']/2/1000 * raws['dt'] * 0.001
		zb = raws['press'] / 10
		z = zb + dz * np.cos(17*np.pi/180) #constant is conversion to depth from distance from transducer due to angle of assent at 17deg
		raws['dz'] = dz
		raws['z'] = z
		return dz, raws
	
	
	def echogram(self, raws, beamsel, cmap = 'zonar', zmin=-999, zmax=-999, **kwargs):
		
		#set color map
		if cmap == "zonar":
			cmap = [[0.0, 'rgb(49,54,149)'], 
					 [0.1111111111111111, 'rgb(69,117,180)'], 
					 [0.2222222222222222, 'rgb(116,173,209)'], 
					 [0.3333333333333333, 'rgb(171,217,233)'], 
					 [0.4444444444444444, 'rgb(224,243,248)'], 
					 [0.5555555555555556, 'rgb(254,224,144)'], 
					 [0.6666666666666666, 'rgb(253,174,97)'], 
					 [0.7777777777777778, 'rgb(244,109,67)'], 
					 [0.8888888888888888, 'rgb(215,48,39)'], 
					 [1.0, 'rgb(165,0,38)']]
			
		color_scale = cmap
			
		#create the heatmap
		sub = raws[raws['beam'] == beamsel]
		
		#set color limits
		if zmin == -999:
			zmin = np.min(sub['Sv'])
		if zmax == -999:
			zmax = np.max(sub['Sv'])
		
		fig = plotly.tools.make_subplots(rows=2,
										 cols = int((sub['Ping'].max()+1)/2),
										 subplot_titles=["Ping "+str(i) for i in range(1,int(sub['Ping'].max())+2)],
										 shared_xaxes=True,
										 shared_yaxes=True,
										 vertical_spacing=0.1,
										 horizontal_spacing=0.1)
		#define number of columns and rows
		cols = list(np.repeat(np.arange(1, int(len(raws['Ping'].unique())/2)+1),2))
		rows = list(np.array(np.tile(np.arange(1, int(len(raws['Ping'].unique())/2)+1),2)))
		
		#loop through the different pings
		for b in range(0,int(sub['Ping'].max() + 1)):
		
			sub1 = sub[sub['Ping'] == b]
			trace = go.Heatmap(z = sub1['Sv'],
							   x = sub1['dive_time'],#np.floor((sub1['nBurst']-1)/2),
							   y = sub1['dz'],
							   zmin = zmin,
							   zmax = zmax,
							   colorscale = color_scale
							   )
			fig.append_trace(trace, row = int(rows[b]), col = int(cols[b]))
			
			
		#fig['layout'].update(yaxis = dict(autorange='reversed', title="Depth [m]"))
		#fig['layout'].update(xaxis = dict( title = "Burst Number"))
		fig['layout']['xaxis1'].update(title='Burst Number')
		fig['layout']['xaxis2'].update(title='Burst Number')
		
		fig['layout']['yaxis1'].update(title='Distance from transducer [m]', autorange='reversed')
		fig['layout']['yaxis2'].update(title='Distance from transducer [m]', autorange='reversed')
		return fig
	
	'''
	# OBSOLETE
	def save_object(self, obj, filename):
		with open(filename, 'wb') as output:  # Overwrites any existing file.
			pkl.dump(obj, output, pkl.HIGHEST_PROTOCOL)
	
	
	fname = 'Z:/Zonar/2018_01_31-02_12_SD_Trough_4/zonar_flash\\B0103'
	#bfile = np.fromfile(fname,dtype='uint8')
	z=Zonar()
	bfile = z.read_raw(fname)
	raws, avg, start,cal = z.read_one_dive(bfile)		
	'''