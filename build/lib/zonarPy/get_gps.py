# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 22:43:17 2020

@author: sven
"""

import xarray as xr

def get_gps(miss):
	"""
	Get the detected bottom depth from a Zooglider netCDF file

	Parameters
	----------
	miss : path and filename of Zooglider netCDF file
		String containing the path and filename of a Zooglider netCDF file 
		containing the Satelite information, where the detected bottom 
		depth is stored.

	Returns
	-------
	dataframe lon and lat start and end for each dive with time

	"""
	#get gps
	gps = xr.open_dataset(miss, group='GPS').to_dataframe()
	
	return(gps)