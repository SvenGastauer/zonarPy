# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 14:29:07 2020

@author: sven
"""

# This is a recreation of the example found in the KML Reference:
# http://code.google.com/apis/kml/documentation/kmlreference.html#gxtrack
import os
import time
import simplekml
import xarray as xr
import pandas as pd
from get_gps import get_gps
import numpy as np

# Data for the track
def get_kml_track(miss,outdir):
	mn = os.path.basename(miss)[:-3]
	print(time.ctime() + ': Processing ' + mn)
	
	print(time.ctime() + ': Getting environment' )
	
	env = xr.open_dataset(miss, group='Environment').to_dataframe()
	env = env.reset_index()
	tmelt = pd.melt(env, id_vars =['Depth','Dive'], value_vars =['time']) 
	tmelt['temperature'] = pd.melt(env, id_vars =['Depth','Dive'], value_vars =['temperature'])['value'] 
	tmelt['salinity'] = pd.melt(env, id_vars =['Depth','Dive'], value_vars =['salinity'])['value'] 
	tmelt['fluorescence'] = pd.melt(env, id_vars =['Depth','Dive'], value_vars =['fluorescence'])['value'] 
	tmelt = tmelt[tmelt.value > 0]
	
	print(time.ctime() + ': Getting GPS ')
	
	gps = get_gps(miss)
	gmelt = pd.melt(gps.reset_index(),id_vars = ['Dive'], value_vars = ['Time_start','Time_end'])
	gmelt['Lon'] = pd.melt(gps.reset_index(),id_vars = ['Dive'], value_vars = ['Lon_start','Lon_end'])['value']
	gmelt['Lat'] = pd.melt(gps.reset_index(),id_vars = ['Dive'], value_vars = ['Lat_start','Lat_end'])['value']
	gmelt = gmelt.sort_values('value')
	gmelt = gmelt.rename(columns={"value": "Time"})
	
	print(time.ctime() + ': Converting for kml')
	
	tmelt['Lon'] = np.interp(tmelt.value.astype('int64'),gmelt['Time'].astype('int64'), gmelt.Lon)
	tmelt['Lat'] = np.interp(tmelt.value.astype('int64'),gmelt['Time'].astype('int64'), gmelt.Lat)
	
	
	tmelt['Pos'] = (tmelt.Dive.astype(str))
	gm = tmelt.sort_values(by='Dive')
	
	g = gm.iloc[np.arange(0,gm.shape[0],2)]
	when = pd.to_datetime(g.value).values.astype(str).tolist()
	coord = g.apply(lambda X: (X["Lon"],X["Lat"],-X["Depth"]) ,axis=1).values.tolist()
	
	temperature = g.temperature.values.tolist()
	salinity = g.salinity.values.tolist()
	fluo = g.fluorescence.values.tolist()
	dive = g.Dive.values.tolist()
	depth = g.Depth.values.tolist()
	
	# Create the KML document
	kml =  simplekml.Kml(name="Zooglider - " + mn, open=1)
	doc =  kml.newdocument(name='Zooglider Mission', snippet= simplekml.Snippet(mn))
	doc.lookat.gxtimespan.begin = str(pd.to_datetime(g.value.min()))
	doc.lookat.gxtimespan.end = str(pd.to_datetime(g.value.max()))
	doc.lookat.longitude = g.Lon.mean()
	doc.lookat.latitude = g.Lat.mean()
	#doc.lookat.range = 2300.000000
	# Create a folder
	fol = doc.newfolder(name='Dives')
	# Create a schema for extended data: heart rate, cadence and power
	schema =  kml.newschema()
	schema.newgxsimplearrayfield(name='dive', type=simplekml.Types.int, displayname='Dive #')
	schema.newgxsimplearrayfield(name='depth', type=simplekml.Types.int, displayname='Depth')
	schema.newgxsimplearrayfield(name='temp', type=simplekml.Types.int, displayname='Temperature')
	schema.newgxsimplearrayfield(name='sal', type=simplekml.Types.int, displayname='Salinity')
	schema.newgxsimplearrayfield(name='fluo', type=simplekml.Types.float, displayname='Fluorescence')
	# Create a new track in the folder
	trk = fol.newgxtrack(name=os.path.basename(miss)[:-3])
	# Apply the above schema to this track
	trk.extendeddata.schemadata.schemaurl = schema.id
	# Add all the information to the track
	trk.newwhen(when) # Each item in the give nlist will become a new <when> tag
	trk.newgxcoord(coord) # Ditto
	trk.extendeddata.schemadata.newgxsimplearraydata('dive', dive) # Ditto
	trk.extendeddata.schemadata.newgxsimplearraydata('depth', depth) # Ditto
	trk.extendeddata.schemadata.newgxsimplearraydata('temp', temperature) # Ditto
	trk.extendeddata.schemadata.newgxsimplearraydata('sal', salinity) # Ditto
	trk.extendeddata.schemadata.newgxsimplearraydata('fluo', fluo) # Ditto
	# Styling
	trk.stylemap.normalstyle.iconstyle.icon.href = 'zg.ico'#'http://earth.google.com/images/kml-icons/track-directional/track-0.png'
	#trk.stylemap.normalstyle.linestyle.color = 'ff6600'
	trk.stylemap.normalstyle.linestyle.width = 6
	trk.stylemap.highlightstyle.iconstyle.icon.href = 'zg.ico'#'http://earth.google.com/images/kml-icons/track-directional/track-0.png'
	#trk.stylemap.highlightstyle.iconstyle.scale = 1.2
	#trk.stylemap.highlightstyle.linestyle.color = 'ff6600'
	trk.stylemap.highlightstyle.linestyle.width = 8
	trk.altitudemode = simplekml.AltitudeMode.relativetoground
	# Save the kml to file
	print(time.ctime() + ': Saving - ' + outdir + mn + ".kml")
	
	kml.save(outdir + mn + ".kml")
	print(time.ctime() + ': Done...')
	