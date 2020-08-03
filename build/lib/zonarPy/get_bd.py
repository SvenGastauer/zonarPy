# -*- coding: utf-8 -*-
"""
Created on Fri May 29 12:12:23 2020

@author: sven
"""
import xarray as xr

def get_bd(miss):
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
	dataframe bottom depth by dive in meters.

	"""
	#add bottom depth detected
	bd = xr.open_dataset(miss, group='Sat/zonar/z').to_dataframe()
	bd = bd.loc[bd.Code == 'ZT01']
	bd['z'] = bd['z']/10
	bd = bd.rename(columns={"Dive#": "Dive"})
	bd = bd.groupby('Dive').max().z
	
	return(bd)