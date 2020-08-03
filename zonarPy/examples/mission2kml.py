# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 22:46:59 2020

@author: sven
"""
from kml2 import get_kml_track
import glob
import time
outdir = r'C:\\Users\\sven\\Documents\\Zonar\\data\\kml\\'

#define folder containing the netCDF files
mdir = r'C:\\Users\\sven\\Documents\\Zonar\\data\\nc_zonar\\'
#get a list of all available mission files
missions = glob.glob(mdir+'*.nc')
	
for miss in missions:
	print(time.ctime() + ': Processing ' + miss)
	get_kml_track(miss, outdir)