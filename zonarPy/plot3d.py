00# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 12:19:28 2019

@author: sveng
"""

import pandas as pd
import numpy as np
from mayavi import mlab

def plot3d(data, xvar,yvar,zvar,fvar,ival='Dive',day=None, filament=None, cmap='RdYlBu', rev_cmap=True, alpha=1,cmin=None,cmax=None,title='', temp = None, chl=None, intdat=None):
    """
    Create 3D curtain plot

    Parameters
    ----------
    data : dataframe
        Pandas dataframe containing the needed data.
    xvar : character
        Name of the x axis variable, in the data dataframe (normally lon - or however longitude is name in the dataframe) .
    yvar : character
        Name of the y axis variable, in the data dataframe  (normally lat - or however latitude is name in the dataframe) .
    zvar : character
        Name of the z axis variable, in the data dataframe (normally a rounded Depth or Pressure value).
    fvar : character
        Name of the variable used for the filling (for example temperature, salinity, Sv1000 or whtaever you want it to be).
    ival : character, optional
        This is a grouping value, if not using 'Dive', make sure to define it. The default is 'Dive'.
    day : dataframe, optional
        Reference to a dataframe containing day/night information, if defined, this will put a bar on top of the plot. The default is None.
    filament : dataframe, optional
        This will put a bar below the plot, for example to delimit certain regions. The default is None.
    cmap : colormap, optional
        Character describing a valid colormap. The default is 'RdYlBu'.
    rev_cmap : boolean, optional
        Define if the colormap should be reversed. The default is True.
    alpha : float, optional
        Sets the transparenct. The default is 1.
    cmin : float, optional
        Define the color scale minimum. The default is None.
    cmax : float, optional
        Define the color scale maximum. The default is None.
    title : character, optional
        Add a title to the plot. The default is ''.
    temp : TYPE, optional
        Add temperature background map. The default is None.
    chl : TYPE, optional
        Add chlorophyll background map. The default is None.
    intdat : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.

    """
    ne=0
    #bring x,y,z data into shape, create 3 pivot tables
    x = pd.pivot(data=data.reset_index(), index=ival,columns=zvar,values = xvar)
    y = pd.pivot(data=data.reset_index(), index=ival,columns=zvar,values = yvar)
    z = -pd.pivot(data=data.reset_index(), index=ival,columns=zvar,values = zvar)
    
    #create pivot table with fill value, where ival become the index
    fill = pd.pivot(data=data.reset_index(), index=ival,columns=zvar,values = fvar)
    
    #define cmax and cmin if they are not given
    if cmax is None:
        cmax = np.ceil(np.nanmax(fill))
    if cmin is None:
        cmin = np.floor(np.nanmin(fill))
        
    #create a mayavi figure
    fig = mlab.figure(size = (2048,1526),\
                bgcolor = (1,1,1), fgcolor = (0.5, 0.5, 0.5), figure=title)
    
    #create a mesh from the x,y,z and fill info
    svp = mlab.mesh(x,y,z,scalars=fill,colormap = cmap, opacity=alpha, extent=[0.2,1,0.1,1,0.22,.69], figure=fig, representation='surface')
    
    #reverse color map if needed
    if rev_cmap == True:
        svp.module_manager.scalar_lut_manager.reverse_lut = True
    #set color range
    svp.module_manager.scalar_lut_manager.data_range = (cmin, cmax)
    #svp.module_manager.scalar_lut_manager.use_below_range_color = True
    #svp.module_manager.scalar_lut_manager.below_range_color = (0,0,0,0)
    #svp.module_manager.scalar_lut_manager.lut.number_of_colors = 64

    #add colorbar
    cb = mlab.colorbar(object=svp, title=title, orientation='vertical')
    cb.scalar_bar.unconstrained_font_size = True
    cb.label_text_property.font_size=20
    cb.title_text_property.font_size=46
    
    #make plot window
    mlab.gcf()
    
    #set plot limits
    xmin=x.values[np.isfinite(x.values)].min();xmax=x.values[np.isfinite(x.values)].max()
    ymin=y.values[np.isfinite(y.values)].min();ymax=y.values[np.isfinite(y.values)].max()
    zmin=z.values[np.isfinite(z.values)].min();zmax=z.values[np.isfinite(z.values)].max()
    
    #change the aesthetics
    ax = mlab.axes(nb_labels=5, 
           xlabel='Longitude',
           ylabel='Latitude',
           zlabel='Depth [m]',
           ranges=[xmin,xmax,ymin,ymax,zmin,zmax])
    ax.axes.font_factor=0.5
    ax.axes.label_format = '    %6.2f'
    
    #add day night at top 
    if day is not None:
        lut = np.zeros((2, 4))
        lut[1,:] = [0,0,0,100]
        lut[0,:] = [100,100,0,100]
        p3d=mlab.plot3d(day.lon.values, day.lat.values, day['depth'].values,
                 day['numSun'].values, tube_radius=0.001,
                 figure=fig, opacity=1,extent=[0.2,1,0.1,1,0.695,0.71])
        p3d.module_manager.scalar_lut_manager.lut.number_of_colors = 2
        p3d.module_manager.scalar_lut_manager.lut.table = lut
        
    #add line at bottom of plot
    if filament is not None:
        mlab.plot3d(filament.lon.values, filament.lat.values, filament['fildep'].values,
                 filament['numfil'].values, tube_radius=0.01,
                 figure=fig, opacity=1,extent=[0.2,1,0.1,1,0.19,0.21],colormap='Set1')
    
    #add chl background map
    if chl is not None:
        if ne == 0:
            zl = 0.2
        else:
            zl = 0.72
        ne+=1
        x1 = chl.longitude.values.min()
        x2 = chl.longitude.values.max()
        y1 = chl.latitude.values.min()
        y2 = chl.latitude.values.max()
        dx = xmax - xmin
        dy = ymax - ymin
        xl1 = -(xmin-x1)/dx
        xl2 = -(xmax-x2)/dx
        yl1 = -(ymin-y1)/dy
        yl2 = -(ymax-y2)/dy

        zs = pd.pivot_table(data=chl, index='latitude', columns='longitude', values='chlorophyll').transpose()
        #cmap2 =  matplotlib.cm.get_cmap()
        #cmap2.set_bad('white',1.)
        sat1 = mlab.imshow(zs, figure=fig,extent=[0.2+xl1,1+xl2,0.1+yl1,1+yl2,zl,zl],opacity=0.8,colormap=cmap, interpolate=False)
        sat1.module_manager.scalar_lut_manager.reverse_lut = True
        sat1.module_manager.scalar_lut_manager.data_range = (0, 10)
        sat1.module_manager.scalar_lut_manager.lut.nan_color = (0, 0, 0, 0)
        cb2 = mlab.colorbar(object=sat1, title='CHla', orientation='horizontal')
        cb2.scalar_bar.unconstrained_font_size = True
        cb2.label_text_property.font_size=20
        cb2.title_text_property.font_size=46
        
    #add temp backgorund map
    if temp is not None:
        if ne == 0:
            zl = 0.2
        else:
            zl = 0.72
        ne+=1
        
        x1 = temp.longitude.values.min()
        x2 = temp.longitude.values.max()
        y1 = temp.latitude.values.min()
        y2 = temp.latitude.values.max()
        dx = xmax - xmin
        dy = ymax - ymin
        xl1 = -(xmin-x1)/dx
        xl2 = -(xmax-x2)/dx
        yl1 = -(ymin-y1)/dy
        yl2 = -(ymax-y2)/dy

        zs2 = pd.pivot_table(data=temp, index='latitude', columns='longitude', values='sst').transpose()
        sat2 = mlab.imshow(zs2, figure=fig,extent=[0.2+xl1,1+xl2,0.1+yl1,1+yl2,zl,zl], colormap=cmap,opacity=0.8, interpolate=False)
        sat2.module_manager.scalar_lut_manager.reverse_lut = True
        sat2.module_manager.scalar_lut_manager.lut.nan_color = (0, 0, 0, 0)
        #sat2.module_manager.scalar_lut_manager.data_range = (10, 25)
    
    #ignore this, as this was specific for the 38 kHz map for calibration
    if intdat is not None:
        x1 = intdat.Lon_S.values.min()
        x2 = intdat.Lon_S.values.max()
        y1 = intdat.Lat_S.values.min()
        y2 = intdat.Lat_S.values.max()
        dx = xmax - xmin
        dy = ymax - ymin
        xi1 = -(xmin-x1)/dx
        xi2 = -(xmin-x2)/dx
        yi1 = -(ymin-y1)/dy
        yi2 = -(ymin-y2)/dy
        
        xi = pd.pivot(data=intdat, index='Ping_S',columns='Depth_mean',values = 'Lon_S')
        yi = pd.pivot(data=intdat, index='Ping_S',columns='Depth_mean',values = 'Lat_S')
        zi = pd.pivot(data=intdat, index='Ping_S',columns='Depth_mean',values = 'Depth_mean')
        
        filli = pd.pivot(data=intdat, index='Ping_S',columns='Depth_mean',values = 'Sv38')
        

        ivp = mlab.mesh(xi,yi,zi,scalars=filli,colormap = cmap, opacity=alpha, extent=[xi1,xi2,yi1,yi2,0.2,.7], figure=fig, representation='surface')
        if rev_cmap == True:
            ivp.module_manager.scalar_lut_manager.reverse_lut = True
        ivp.module_manager.scalar_lut_manager.data_range = (cmin, cmax)
    
    #set initial view    
    mlab.view(90,110,-2)




'''
env = xr.open_dataset(fname, group = 'Environment')
meta = env.to_dataframe()
env.close()
meta = meta.reset_index()
pmin = meta.Depth.min()
meta.Depth = meta.Depth - pmin
#for some reason the meta data and the acoustic dates are one day apart
meta['deltatime'] = meta['deltatime'].fillna(0)
meta['time_of_measure'] = meta['StartTime'] + pd.to_timedelta(meta['deltatime'].astype('str') + 'seconds')

meta['time_of_measure'] = pd.to_datetime(meta['time_of_measure']) -  pd.Timedelta(days=1)
meta['StartTime'] = pd.to_datetime(meta['StartTime']) -  pd.Timedelta(days=1)
meta['EndTime'] = pd.to_datetime(meta['EndTime']) -  pd.Timedelta(days=1)

#remove meta where no temperature is available
meta = meta.dropna()
#get time passed since dive start
meta['t_passed'] = meta['time_of_measure'] - meta['StartTime']
#get total dive duration
meta['dt_tot'] = meta['StartTime'] - meta['EndTime']
#get proportion  of dive passed
meta['prop_time_el'] = meta['t_passed'] / meta['dt_tot']
meta['lat'] = meta['Lat_start'] + meta['prop_time_el'] * (meta['Lat_end'] - meta['Lat_start'])
meta['lon'] = meta['Lon_start'] + meta['prop_time_el'] * (meta['Lon_end'] - meta['Lon_start'])
meta = meta.sort_values(by='time_of_measure')
mlab.figure(size = (1024,768),\
                bgcolor = (1,1,1), fgcolor = (0.5, 0.5, 0.5))
mlab.plot3d(meta.lon.values, meta.lat.values, -meta.Depth.values, meta.salinity, extent=[0,1,0,1,0,1], tube_radius=0.003)
'''

'''
# plot with matplotlib
fig = plt.figure(figsize=plt.figaspect(0.5))

#sv plot
ax = fig.add_subplot(1, 2, 1, projection='3d')

color_dimension = sv.values # change to desired fourth dimension
minn, maxx = color_dimension[np.isfinite(color_dimension)].min(), color_dimension[np.isfinite(color_dimension)].max()
minn=-80; maxx=-30
norm = matplotlib.colors.Normalize(minn, maxx)
m = plt.cm.ScalarMappable(norm=norm, cmap='RdYlBu_r')
m.set_array([])
svcolors = m.to_rgba(color_dimension)

ax.plot_surface(x, y, z, rstride=1, cstride=1, linewidth=0,
                antialiased=False, facecolors=svcolors)
#ax.plot_surface(x,y,z, facecolors=svcolors, vmin=minn, vmax=maxx, shade=False)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

#temp plot
ax = fig.add_subplot(1, 2, 2, projection='3d')

color_dimension = temp.values # change to desired fourth dimension
minn, maxx = color_dimension[np.isfinite(color_dimension)].min(), color_dimension[np.isfinite(color_dimension)].max()
norm = matplotlib.colors.Normalize(minn, maxx)
m = plt.cm.ScalarMappable(norm=norm, cmap='RdYlBu_r')
m.set_array([])
tempcolors = m.to_rgba(color_dimension)

ax.plot_surface(x,y,z, rstride=1, cstride=1, linewidth=0,
                antialiased=False, facecolors=tempcolors, vmin=minn, vmax=maxx, shade=False)


'''