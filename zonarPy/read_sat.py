"""
Created on Mon Nov 18 09:30:05 2019

@author: sveng
"""

#SAT file reader
import pandas as pd
import time
import numpy as np
from astral import Astral

def read_sat(fn):
	"""
	read sat file, the information send via satellite during the voyage
	
	Parameters
	----------
	fn : str
		path to sat file inclduing filename

	Returns
	-------
	header : dict
		Mission header information, including non-acoustic calibration values.
	gps : pandas DataFrame
		GPS position information of the individual dives
	engineering : pandas DataFrame
		Engineering informaiton, including drift and heading info, pitch and roll, compass
	profile : pandas DataFrame
		Zoocam profile information
	zoocam : pandas DataFrame
		Zoocam summary profile information
	zonar : pandas DataFrame
		Zonar summary information
	misc : pandas DataFrame
		Misc Data

	"""
	print(time.ctime() + ': Processing sat data...' )

	f = open(fn)
	sat = f.readlines()
	f.close()
	
	#generate index
	
	with open(fn) as f:
		list2 = [row.split()[0] for row in f]
		##############
		# HEADER LINES
		##############
		
		h_mod = [i for i, j in enumerate(list2) if 'MOD' in j]
		h_md = [i for i, j in enumerate(list2) if 'MD' in j]
		h_vn = [i for i, j in enumerate(list2) if 'VN' in j]
		h_vo = [i for i, j in enumerate(list2) if 'VO' in j]
		h_va = [i for i, j in enumerate(list2) if 'VA' in j]
		h_cp = [i for i, j in enumerate(list2) if 'CP' in j]
		h_ct = [i for i, j in enumerate(list2) if 'CT' in j]
		h_cs = [i for i, j in enumerate(list2) if 'CS' in j]
		h_co = [i for i, j in enumerate(list2) if 'CO' in j]
		h_cd = [i for i, j in enumerate(list2) if 'CD' in j]
		h_ch = [i for i, j in enumerate(list2) if 'CH' in j]
		h_cx = [i for i, j in enumerate(list2) if 'CX' in j]
		
		#############
		# GPS DATA
		#############
		gps_G = [i for i, j in enumerate(list2) if j == 'G']
		gps_R = [i for i, j in enumerate(list2) if j == 'R']
		gps_r = [i for i, j in enumerate(list2) if j == 'r']
		gps_W = [i for i, j in enumerate(list2) if j == 'W']
		gps_w = [i for i, j in enumerate(list2) if j == 'w']
		
		###########
		# ENGINEERING DATA
		###########
		en_EC = [i for i, j in enumerate(list2) if 'EC' in j]
		en_ED = [i for i, j in enumerate(list2) if 'ED' in j]
		en_EF01 = [i for i, j in enumerate(list2) if 'EF01' in j]
		en_EF02 = [i for i, j in enumerate(list2) if 'EF02' in j]
		en_EN = [i for i, j in enumerate(list2) if 'EN' in j]
		en_EP01 = [i for i, j in enumerate(list2) if 'EP01' in j]
		en_EP02 = [i for i, j in enumerate(list2) if 'EP02' in j]
		en_ET = [i for i, j in enumerate(list2) if 'ET' in j]
		#en_e = [i for i, j in enumerate(list2) if 'e' in j]
		en_EY = [i for i, j in enumerate(list2) if 'EY' in j]
		en_ES = [i for i, j in enumerate(list2) if 'ES' in j]
		en_EW = [i for i, j in enumerate(list2) if 'EW' in j]
		en_EX = [i for i, j in enumerate(list2) if 'EX' in j]
		
		#########
		# PROFILE DATA
		#########
		p_D = [i for i, j in enumerate(list2) if j == 'D']
		p_p = [i for i, j in enumerate(list2) if j == 'p']
		
		########
		# ZOOCAM
		########
		zc_ZC = [i for i, j in enumerate(list2) if 'ZC' in j]
		zc_ZS01 = [i for i, j in enumerate(list2) if 'ZS01' in j]
		zc_ZS02 = [i for i, j in enumerate(list2) if 'ZS02' in j]
		
		#######
		# ZONAR
		########
		zc_ZN00 = [i for i, j in enumerate(list2) if 'ZN00' in j]
		zc_ZN01 = [i for i, j in enumerate(list2) if 'ZN01' in j]
		zc_ZN02 = [i for i, j in enumerate(list2) if 'ZN02' in j]
		zc_ZT01 = [i for i, j in enumerate(list2) if 'ZT01' in j]
		zc_ZT02 = [i for i, j in enumerate(list2) if 'ZT02' in j]
		z = [i for i, j in enumerate(list2) if 'z' in j]
		
		########
		# DOPPLER DATA
		########
		#dp_B1 = [i for i, j in enumerate(list2) if j == 'B1']
		#dp_AS = [i for i, j in enumerate(list2) if j == 'AS']
		#dp_AT = [i for i, j in enumerate(list2) if j == 'AT']
		#dp_a = [i for i, j in enumerate(list2) if j == 'a']
		
		#######
		# MISC
		#######
		misc_dive = [i for i, j in enumerate(list2) if j == '!dive']
		misc_E = [i for i, j in enumerate(list2) if j == 'E']
		misc_M = [i for i, j in enumerate(list2) if j == 'M']
		misc_S = [i for i, j in enumerate(list2) if j == 'S']
		misc_SBD = [i for i, j in enumerate(list2) if j == 'SBD']
		misc_X = [i for i, j in enumerate(list2) if j == 'X']
		misc_x = [i for i, j in enumerate(list2) if j == 'x']
	f.close()
	
	
	#######################
	# HEADER INFORMATION
	#######################
	#MODIFY version # comment
	if len(h_mod) > 0:
		df_MOD = pd.DataFrame([sat[i] for i in h_mod])
		df_MOD.columns = ['Comment']
	else:
		df_MOD = []
	
	#mission description
	if len(h_md) > 0:
		MD = [sat[i].split() for i in h_md]
		df_MD = pd.DataFrame(MD[0][0:3] + [' '.join(MD[0][3:])]).transpose() 
		df_MD.columns = ['Spray_SN','Deployment_ID','E_Name','Description']
	else:
		df_MD = []
	
	#vehicle information
	if len(h_vn) > 0:
		df_VN = pd.DataFrame([sat[i].split() for i in h_vn]) 
		df_VN.columns = ['code','SN','nsensor','CTD_Type','EEPROM__Ver'] + ['s'+ str(i+1) for i in range(int(df_VN[2][0]))]
	else:
		df_VN = []
	
	#VO string:  string = ‘FLUOR’ for fluorometer, ‘OBS’ for optical backscatter
	#is used to describe the type of optical sensor mounted on Spray.
	if len(h_vo) > 0:
		df_VO = pd.DataFrame([sat[i].split() for i in h_vo]) 
	else:
		df_VO = []
	
	#VA ptt_id = Argos PTT ID.
	#Argos ID = Argos PTT ID installed on this Spray.
	if len(h_va) > 0:
		df_VA = pd.DataFrame([sat[i].split()[0:2] for i in h_va]) 
		df_VA.columns = ['Code','Argos_ID']
	else:
		df_VA = []
	
	## Calibration coefficients
	# CP - Pressure CT - Temperature CS - Salinity CO - Optics
	if len(h_cp + h_ct + h_cs) > 0:
		df_cal = [sat[i].split() for i in h_cp + h_ct + h_cs]
		df_cal = [df_cal[i]  + (5-len(df_cal[i]))*[0] + [df_cal[i][0][1]] + [df_cal[i][0][2]] + [df_cal[i][0][3]] for i in range(len(df_cal))]
		df_cal = pd.DataFrame(df_cal)
		df_cal.columns = ['Code','offset','gain','off2','gn2','stype','stype_c','cal_format']
	else:
		df_cal = []
	
	#Optical Coefficients
	if len(h_co)>0:
		df_CO = pd.DataFrame([sat[i].split() for i in h_co]) 
		df_CO.columns = ['Code','offset','gain','hg','ogain','off2','gn2']
	else:
		df_CO=[]
	
	#Doppler Backscatter (ABS) Coefficients
	if len(h_cd)>0:
		df_CD = pd.DataFrame([sat[i].split() for i in h_cd]) 
		df_CD.columns = ['Code','bin_ave_start','bin_ave_end','off2','gn2','bd','alpha','name']
	else:
		df_CD = []
	
	#compass calibration
	if len(h_ch)>0:
		df_CH = pd.DataFrame([sat[i].split() for i in h_ch]) 
		df_CH.columns = ['Code','MagVar','UseMC','file_name']
	else:
		df_CH = []
	
	#DO calibration
	if len(h_cx)>0:
		df_CX = pd.DataFrame([sat[i].split() for i in h_cx]) 
		df_CX.columns = ['Code','MagVar','UseMC','file_name']
	else:
		df_CX = []
	
	header = {'modify':df_MOD,
			  'mission':df_MD,
			  'vehicle':df_VN,
			  'optical_sensor':df_VO,
			  'argos':df_VA,
			  'calibration':df_cal,
			  'cal_optical':df_CO,
			  'cal_doppler':df_CD,
			  'cal_compass':df_CH,
			  'cal_DO':df_CX}
	
	##########################
	# GPS INFORMATION
	##########################
	
	if len(gps_G)>0:
		df_G = pd.DataFrame([sat[i].split() for i in gps_G]) 
		df_G.columns = ['Code','dive#','mission_status','day','month','year',
					 'time','valid_flag','Lat_d','Lat_m','Lon_d','Lon_m', 'Tfix',
					 'Nsat','minSNR','meanSNR', 'maxSNR','HDOP','GPS_HEALTH',
					 'Wing_ind_Roll_stat', 'Latitude','Longitude']
	else:
		df_G = []
		
	
	if len(gps_W)>0:
		df_W = pd.DataFrame([sat[i].split() for i in gps_W]) 
		df_W.columns = ['Code','dive#','n_waypoints']
	else:
		df_W = []
	
	
	if len(gps_w)>0:
		df_w = pd.DataFrame([sat[i].split() for i in gps_w]) 
		df_w.columns = ['Code','waypoint', 'Latitude','Longitude']
	else:
		df_w = []
	
	
	if len(gps_R)>0:
		df_R = pd.DataFrame([sat[i].split() for i in gps_R]) 
		df_R.columns = ['Code','dive#','n_waypoints','route_entry','end-of-route-action',
					 'direction','current-bucking','current-crossing-angle',
					 'max_dive_ccm','steering_manual','max_dive_sm','steering_point',
					 'max_dive_sp', 'cor_angle_min','cor_angle_max']
	else:
		df_R = []
		
	if len(gps_r)>0:
		df_r = pd.DataFrame([sat[i].split() for i in gps_r]) 
		df_r.columns = ['Code','route_entry_index','detect_mode','watch-circle_radius','angle_approach']
	else:
		df_r = []
	
	gps = {'GPS_info':df_G,
		   'Waypoints':df_W,
		   'Waypoints_values':df_w,
		   'Route_info':df_R,
		   'Route_data':df_r}
	
	
	#######################
	# ENGINEERING Data
	#######################
	
	if len(en_EC)>0:
		df_EC = pd.DataFrame([sat[i].split() for i in en_EC]) 
		df_EC.columns = ['Code','dive#','Ntries','Nsent','SBDI_STAT',
					 'SBD_SHOREW_STAT','T_SBD','Wing_used']
	else:
		df_EC = []
	
	if len(en_ED)>0:
		df_ED = pd.DataFrame([sat[i].split() for i in en_ED]) 
		df_ED.columns = ['Code','dive#','Pstart','Pavg','Pend',
					 'Drift_Tm','Num_Samp','Pmp_Tm']
	else:
		df_ED = []
	
	
	if len(en_EF01)>0:
		df_EF01 = pd.DataFrame([sat[i].split() for i in en_EF01]) 
		df_EF01.columns = ['Code','dive#','Navg','Psurf','Zmax',
					 'Pitch', 'Altimeter','ADP_alt','Roll_ERR',
					 'EXEC_STATUS']
	else:
		df_EF01 = []
		
	if len(en_EF02)>0:
		df_EF02 = pd.DataFrame([sat[i].split() for i in en_EF02]) 
		df_EF02.columns = ['Code','dive#','Navg','Psurf','Zmax',
					 'Pitch', 'Altimeter','ADP_alt','Z_AT_ALT','Roll_ERR',
					 'EXEC_STATUS']
		df_EF02['Bottom_depth'] = df_EF02['Z_AT_ALT'].to_numpy().astype(float) + df_EF02['Altimeter'].to_numpy().astype(float)
	else:
		df_EF02 = []
	
	
	if len(en_EN)>0:
		df_EN = pd.DataFrame([sat[i].split() for i in en_EN]) 
		df_EN.columns = ['Code','dive#','DRx','Dry','WLAT','WLON', 
					  'Tleave','Tend','Tslow','Heading_desired','TF']
	else:
		df_EN = []
	
	if len(en_EP01)>0:
		df_EP01 = pd.DataFrame([sat[i].split() for i in en_EP01]) 
		df_EP01.columns = ['Code','dive#','Depth_MAX','Volts','Amps','NBAD_AMP', 
					  'T_max_amp','Max_Amps','Vacuum','T_PUMP_1','T_PUMP_2']
	else:
		df_EP01 = []
		
	if len(en_EP02)>0:
		df_EP02 = pd.DataFrame([sat[i].split() for i in en_EP02]) 
		df_EP02.columns = ['Code','dive#','Depth_MAX','Volts','Amps','NBAD_AMP', 
					  'T_max_amp','Max_Amps','Vacuum','T_PUMP_1','T_PUMP_2','VENT_TM','AIR_FLAG']
	else:
		df_EP02 = []
		
	if len(en_ET)>0:
		df_ET = pd.DataFrame([sat[i].split() for i in en_ET]) 
		df_ET.columns = ['Code','dive#','N_datum_time_series','N_datum_output',
					  'Packet_ID']
		#n of time series lines
		ne = np.ceil(df_ET.iloc[:,2].astype(int) / df_ET.iloc[:,3].astype(int)).astype(int)
		eng = pd.DataFrame()
		code = ['time [s] since start-of-dive (wing down and valve opened).',
		  'Pressure [counts]  (process the same as the profile data).',
		  'Heading [degrees*10] from the compass.',
		  'Pitch Angle [degrees*10] from the compass tilt sensor.',
		  'Roll Angle [degrees*10] from the compass tilt sensor.',
		  'Pitch potentiometer [counts].',
		  'Roll potentiometer [counts].']

		for i in range(len(en_ET)):
			tmp = pd.DataFrame(sat[en_ET[i]+1 : en_ET[i]+ne[i]]).iloc[:,0].str.split(expand=True)
			tmp['Code'] = df_ET.iloc[i,0]
			tmp['dive#'] = df_ET.iloc[i,1]
			tmp['Packet_ID'] = df_ET.iloc[i,4]
			tmp['Sensor'] = code[int(tmp['Code'][0][-1])]
			eng = eng.append(tmp)
		df_ET = eng
		
	else:
		df_ET = []
	
	if len(en_EY)>0:
		df_EY = pd.DataFrame([sat[i].split() for i in en_EY]) 
		df_EY.columns = ['Code','dive#', 'Voltage_Bat','Pressure_counts','Exception_stat','Stack_reg']
	else:
		df_EY = []
		  
	if len(en_ES)>0:
		df_ES = pd.DataFrame([sat[i].split() for i in en_ES]) 
		df_ES.columns = ['Code','dive#', 'Valid','Surf_drift_t','Surf_drift_dx','Surf_drift_dy',
					  'Surf_speed','Surf_dir']
	else:
		df_ES = []
	
	
	if len(en_EW)>0:
		df_EW = pd.DataFrame([sat[i].split() for i in en_EW])
		df_EW.columns = ['Code','dive#','Watt_hyd_motor','Watt_pitch','Watt_roll',
					  'Watt_dummy','Time_motor_ran_hyd', 'Time_motor_ran_pitch',
					  'Time_motor_ran_roll','Time_motor_ran_dummy','Volt_Sat',
					  'Volt_pump_bat','SPI_error_GPIO','SPI_error_HYD', 
					  'SPI_error_PITCH','SPI_error_roll','watchdog_error']
	else:
		df_EW = []
		
	if len(en_EX)>0:
		df_EX = pd.DataFrame([sat[i].split() for i in en_EX]) 
		df_EX.columns = ['Code','dive#', 'Profile#','Surfacing#','Cycle#','Unix_time']
	else:
		df_EX = []
	
	fl = df_EF01 if len(df_EF01)>0 else df_EF02
	pl = df_EP01 if len(df_EP01)>0 else df_EP02
	engineering = {'Comm_line':df_EC,'Drift_mode_line':df_ED,'Flight_Line':fl,
				   'Nav_Line':df_EN,'Pump_Line':pl,'Engineering_TS':df_ET,
				   'Stack_Reg_Stat':df_EY,'Surf_drift':df_ES,'Motor_Energy':df_EW,
				   'Profile_Dive_count':df_EX}
	
	###########################
	# PROFILE DATA
	###########################
	if len(p_D)>0:
		df_D = pd.DataFrame([sat[i].split()[0:3] for i in p_D]) 
		df_D.columns = ['Code','dive#', 'npts']
	else:
		df_D = []
	
	if len(p_p)>0:
		df_p = pd.DataFrame([sat[i].split() + [0]*(8-len(sat[i].split())) for i in p_p]) 
		
		df_p.columns = ['Code','dive#', 'packet_ind','pressure_counts','temperature_counts','salinity_counts','optical_counts','DO_counts']
	else:
		df_p = []
	
	
	
	profile = {'Start_Profile_data':p_D,'Profile_data':df_p}
	
	######################
	# ZOOCAM
	######################
	
	if len(zc_ZC)>0:
		df_ZC = pd.DataFrame([sat[i].split() for i in zc_ZC]) 
		if df_ZC.shape[1] ==15:
			df_ZC.columns = ['Code','dive#', 'start-status','stop-status','NumImage',
						  'NumError','AvgCenter','Unix_time','zState', 
						  'outDir','Mb_trans','Trans_speed','Gb_rem_SD','Gb_rem_USB',
						  'UseZcam']
		elif df_ZC.shape[1] ==10:
			df_ZC.columns = ['Code','dive#', 'start-status','stop-status','NumImage',
						  'NumError','AvgCenter','Unix_time','zState', 
						  'outDir']
		else:
			cnams = ['Code','dive#', 'start-status','stop-status','NumImage',
						  'NumError','AvgCenter','Unix_time','zState', 'startNtry',
						  'outDir','Mb_trans','Trans_speed','Gb_rem_SD','Gb_rem_USB',
						  'UseZcam','UvLED','zCamMod','Unix_time2']
			df_ZC.columns = cnams[0:df_ZC.shape[1]]
	
	else:
		df_ZC = []
		
	if len(zc_ZS01)>0:
		df_ZS01 = pd.DataFrame([sat[i].split() for i in zc_ZS01])
		df_ZS01.columns = ['Code','dive#', 'version', 'type','#dataPoints', 'nCol', 'pID'] 
		
		zmin = np.hstack([zc_ZS01])+1
		zmax = zmin + np.ceil(df_ZS01['#dataPoints'].to_numpy().astype(int)/df_ZS01['nCol'].to_numpy().astype(int))
		z = np.hstack([np.arange(zmin[i],zmax[i]) for i in range(len(zmin))]).astype(int)
		
		zz = pd.DataFrame([sat[i].split() for i in z]) 
		zz =zz.drop(zz.columns[0], axis=1).to_numpy()
		zz=zz[zz!=None].astype(int)
		df_z = pd.DataFrame({'Code':np.repeat(df_ZS01['Code'].to_numpy(),df_ZS01['#dataPoints'].to_numpy().astype(int)),
					   'Dive#':np.repeat(df_ZS01['dive#'].to_numpy().astype(int),df_ZS01['#dataPoints'].to_numpy().astype(int)),
					   '#dataPoints': np.hstack([np.flip(np.arange(df_ZS01['#dataPoints'].to_numpy().astype(int)[i])) for i in range(len(df_ZS01['#dataPoints']))]),
					   'z':zz})
		
	else:
		zc_ZS01 = []
	
	if len(zc_ZS02)>0:
		df_ZS02 = pd.DataFrame([sat[i].split() for i in zc_ZS01 + zc_ZS02])
		df_ZS02.columns = ['Code','dive#', 'version', 'type','#dataPoints', 'nCol', 'pID'] 
		
		zmin = np.hstack([zc_ZS01 + zc_ZS02])+1
		zmax = zmin + np.ceil(df_ZS02['#dataPoints'].to_numpy().astype(int)/df_ZS02['nCol'].to_numpy().astype(int))
		z = np.hstack([np.arange(zmin[i],zmax[i]) for i in range(len(zmin))]).astype(int)
		
		zz = pd.DataFrame([sat[i].split() for i in z]) 
		zz =zz.drop(zz.columns[0], axis=1).to_numpy()
		zz=zz[zz!=None].astype(int)
		df_z = pd.DataFrame({'Code':np.repeat(df_ZS02['Code'].to_numpy(),df_ZS02['#dataPoints'].to_numpy().astype(int)),
					   'Dive#':np.repeat(df_ZS02['dive#'].to_numpy().astype(int),df_ZS02['#dataPoints'].to_numpy().astype(int)),
					   '#dataPoints': np.hstack([np.flip(np.arange(df_ZS02['#dataPoints'].to_numpy().astype(int)[i])) for i in range(len(df_ZS02['#dataPoints']))]),
					   'z':zz})
	else:
		zc_ZS02 = []
		
		
	zoocam = {'Zoocam':df_ZC,'ROI_TS':zc_ZS01, 'ROI_count':zc_ZS02}
	
	##########################
	# ZONAR
	##########################
	
	if len(zc_ZN00)>0:
		df_ZN00 = pd.DataFrame([sat[i].split() for i in zc_ZN00]) 
		df_ZN00.columns = ['Code','dive#', 'Unix_time_last','zState','status',
					  'nPings','dtBurst','warmup','tBin', 'barker','nBit']
	else:
		df_ZN00 = []
	
	if len(zc_ZN01 + zc_ZN02)>0:
		df_ZN = pd.DataFrame([sat[i].split() for i in zc_ZN01 + zc_ZN02]) 
		df_ZN.columns = ['Code','dive#', 'f1','pulse','blank',
					  'dt','tScan','tPing','nScan', 'tWait','nBin']
	else:
		df_ZN = []
	 #+ zc_ZT02 + zc_ZS01 + zc_ZS02
	if len(zc_ZT01+ zc_ZT02 + zc_ZS01 + zc_ZS02)>0:
		df_ZT = pd.DataFrame([sat[i].split() for i in zc_ZT01+ zc_ZT02 + zc_ZS01 + zc_ZS02]) 
		df_ZT.columns = ['Code','dive#', 'beam#','bin','nScan','nCol','Xpid']
		zmin = np.hstack([zc_ZT01 + zc_ZT02+ zc_ZS01 + zc_ZS02])+1
		zmax = zmin + np.ceil(df_ZT['nScan'].to_numpy().astype(int)/df_ZT['nCol'].to_numpy().astype(int))
		z = np.hstack([np.arange(zmin[i],zmax[i]) for i in range(len(zmin))]).astype(int)
		
	else:
		df_ZT = []
		z=[]	
	if len(z)>0:
		zz = pd.DataFrame([sat[i].split() for i in z]) 
		zz =zz.drop(zz.columns[0], axis=1).to_numpy()
		zz=zz[zz!=None].astype(int)
		df_z = pd.DataFrame({'Code':np.repeat(df_ZT['Code'].to_numpy(),df_ZT['nScan'].to_numpy().astype(int)),
					   'Dive#':np.repeat(df_ZT['dive#'].to_numpy().astype(int),df_ZT['nScan'].to_numpy().astype(int)),
					   'beam#':np.repeat(df_ZT['beam#'].to_numpy().astype(int),df_ZT['nScan'].to_numpy().astype(int)),
					   'bin':np.repeat(df_ZT['bin'].to_numpy().astype(int),df_ZT['nScan'].to_numpy().astype(int)),
					   'index': np.hstack([np.flip(np.arange(df_ZT['nScan'].to_numpy().astype(int)[i])) for i in range(len(df_ZT['nScan']))]),
					   'z':zz})
		#df_z.columns = ['z']
	else:
		df_z = []
	
	zonar = {'Zonar_gen':df_ZN00,'Zonar_beam':df_ZN,'Zonar_depth-time':df_ZT,'z':df_z}
	
	#####################
	# MISC
	#####################
	
	if len(misc_dive)>0:
		df_dive = pd.DataFrame([sat[i].split() for i in misc_dive]) 
		df_dive.columns = ['Code','dive#', 'OBSOLETE','email_rec_date','email_rec_time']
	else:
		df_dive = []
	
	
	if len(misc_E)>0:
		df_E = pd.DataFrame([sat[i].split() for i in misc_E]) 
		df_E.columns = ['Code','dive#', 'max_depth','Battery_Volts','Pump_current',
					 'altimeter','heading_comp', 'pitch','east','north','lat_wp',
					 'lon_wp', 'Ntries','Nsent','Navg','Wing','SBD','Stat_shore_com',
					 'Surf_pressure','Tpump1','Tpump2','Vacuum', 'time_max_amp','amp_max',
					 'error_int','time_last_mes_SBD','time_leave_surf','time_from_surf',
					 'exception_status']
	else:
		df_E = []
	
	
	if len(misc_M)>0:
		df_M = pd.DataFrame([sat[i].split() for i in misc_M]) 
		df_M.columns = ['Code','dive#', 'year','month','mission_id']
	else:
		df_M = []
		
	if len(misc_S)>0:
		df_S = pd.DataFrame([sat[i].split() for i in misc_S]) 
	else:
		df_S = []
	
	if len(misc_SBD)>0:
		df_SBD = pd.DataFrame([sat[i].split()[0:11] for i in misc_SBD]) 
		df_SBD.columns = ['Code','dive#', 'Spray_SN','SBD_imei','Year-day','valid_flag',
					   'MOMSN_1','bytes_mess','lat','lon','stdev_pos']
	else:
		df_SBD = []
	
	if len(misc_X)>0:
		df_X = pd.DataFrame([sat[i].split() for i in misc_X]) 
		df_X.columns = ['Code','dive#', 'ver','N_param']
	else:
		df_X = []
		
	if len(misc_x)>0:
		df_x = pd.DataFrame([sat[i].split() for i in misc_x]) 
		df_x.columns = ['Code','address', 'value']
	else:
		df_x = []
	
	misc = {'Dive_info':df_dive,'Mission_ID':df_M,'Shore_comm':df_S,'Email_mess':df_SBD,'Engineering':df_E,'Param_list': df_X,'Param_val':df_x}
	

	return header,gps,engineering, profile,zoocam, zonar, misc

def get_twilight_index(dtime, Longitude,Latitude,cf=2):
	a = Astral()
	sunpos = np.array([a.solar_elevation(dtime[cf*i], Latitude[i], Longitude[i]) for i in range(len(Latitude))])
	twilight_index = np.where(sunpos > 0,0,np.NaN)
	twilight_index = np.where(sunpos < -12,2,twilight_index)
	twilight_index[np.where(np.isnan(twilight_index))] = twilight_index[np.where(np.isnan(twilight_index))[0]-1]+1
	while len(np.where(np.isnan(twilight_index))[0]) > 0:
		twilight_index[np.where(np.isnan(twilight_index))] = twilight_index[np.where(np.isnan(twilight_index))[0]-1]
	twilight_index = twilight_index.astype(int)
	return twilight_index
