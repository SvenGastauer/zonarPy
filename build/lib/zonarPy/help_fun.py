# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 16:19:19 2020

@author: sveng
"""


import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from datetime import datetime

# vectorized haversine function
def haversine(lat1, lon1, lat2, lon2, to_radians=True, earth_radius=6371):
	"""
	slightly modified version: of http://stackoverflow.com/a/29546836/2901002
	
	Calculate the great circle distance between two points
	on the earth (specified in decimal degrees or in radians)
	
	All (lat, lon) coordinates must have numeric dtypes and be of equal length.
	
	
	"""
	if to_radians:
		lat1 = np.radians(lat1)
		lat2 = np.radians(lat2)
		lon1 = np.radians(lon1)
		lon2 = np.radians(lon2)
	a = np.sin((lat2-lat1)/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2-lon1)/2.0)**2
	
	return earth_radius * 2 * np.arcsin(np.sqrt(a))


#compute c for each point, we need 
def compute_c(D, S, T, lat = None, method = 'mackenzie'):
	'''
	
	Description
	-----------
	
	Estiamtes the sound speed in water given the provided conditions using either methods descrubed by Mackenzie (1981), Coppens (1981) or Leroy (2008)
	
	If method == 'mackenzie':
	    
	    Calculate speed of sound in seawater based on MacKenzie (1981)
	    The empirical equation generally holds validity for a temperature range between 2 and 30 degrees Celsius, Salinities between 25 and 40 parts per thousand and a depth range between 0 and 8000 m
	    
	    source: Mackenzie, K. V. (1981). Nine term equation for sound speed in the oceans. The Journal of the Acoustical Society of America, 70(3), 807-812.
	    http://asa.scitation.org/doi/abs/10.1121/1.386920
	
	If method == 'coppens':
	    
	    Calculates speed of sound in seawater based on Coppens (1981)
	    The empirical equation generally holds validity for a temperature range between 0 and 35 degrees Celsius, Salinities between 0 and 45 parts per thousand and a depth range between 0 and 4000 m
	    
	    source: Coppens, A. B. (1981). Simple equations for the speed of sound in Neptunian waters. The Journal of the Acoustical Society of America, 69(3), 862-863.
	    http://asa.scitation.org/doi/abs/10.1121/1.385486
	    
	If method == 'leroy':
	    
	    A new equation for the accurate calculation of sound speed in all oceans. The Journal of the Acoustical Society of America, 124(5), 2774-2782.
	    Returns the sound speed according to Leroy et al (2008). This "newer" equation should solve the sound speed within 0.2 m/s for all seas, including the Baltic and Black sea, based on Temperature, Salinity and Latitude. Exceptions are some seas with anomalities close to the bottom. The equation was specifically designed to be used in marine acoustics.
	    
	    source: Leroy, C. C., Robinson, S. P., & Goldsmith, M. J. (2008). A new equation for the accurate calculation of sound speed in all oceans. The Journal of the Acoustical Society of America, 124(5), 2774-2782.
	    http://asa.scitation.org/doi/abs/10.1121/1.2988296
	
	Parameters
	----------
	
	D: float
	    Depth in meters
	
	S: float 
	    Salinity in parts per thousands
	
	T: float
	    Temperature in degrees Celsius
	
	lat: float
	    Latitude needed for Leroy method
	
	Examples
	---------
	
	c_Mackenzie1981(100,35,10)
	
	Returns
	-------
	
	c: float
	    Estimated sound velocity in m/s
	
	'''
	
	if method == 'mackenzie':
	    c = 1448.96 + 4.591 * T - 5.304 * 10**(-2) * (T**2) + 2.374 * \
	    (10**(-4)) * (T**3) + 1.340 * (S-35) + 1.630 * (10**(-2)) * D + \
	    1.675 * (10**(-7)) * (D**2) - 1.025 * (10**(-2)) * T * (S - 35) - \
	    7.139 * (10**(-13)) * T * (D**3)
	
	if method == 'coppens':
	    t = T/10
	    D = D/1000
	    c0 = 1449.05 + 45.7 * t - 5.21 * (t**2)  + 0.23*(t**3)  + \
	    (1.333 - 0.126*t + 0.009*(t**2)) * (S - 35)
	    c = c0 + (16.23 + 0.253*t)*D + (0.213-0.1*t)*(D**2)  + \
	    (0.016 + 0.0002*(S-35))*(S- 35)*t*D
		   
	if method == 'leroy':
	    c = 1402.5 + 5 * T - 5.44 * 10**(-2) * T**2 + 2.1 * 10**(-4) * T**3 + \
		1.33 * S - 1.23 * (10**(-2)) * S * T + 8.7 * (10**(-5)) * \
		S * T** 2 + 1.56 * (10**(-2)) * D + 2.55 * (10**(-7)) * D**2 - 7.3 * \
		(10**(-12)) * D**3 + 1.2 *(10**(-6)) * D * (lat - 45) - 9.5 *\
		(10**(-13)) * T * D**3 + 3 * (10**(-7)) * T**2 * D + 1.43 * \
		(10**(-5)) * S * D
		
	return c

def absorption(f,S,T,D,pH=None, method='Doonan'):
	return alpha_FrancoisGarrison(f=f,T=T,S=S,D=D,pH=8)
'''
         Returns the Absorption coefficient [dB/km] for
         the given acoustic frequency (f [kHz]), salinity
         (S, [ppt]), temperature (T, [degC]), and depth
         (D, [m]).
        
         Note that the salinity units are ppt, not the more
         modern psu. See the comment on page 2 of Doonan et al.
         for a discussion of this. For most purposes, using psu
         values in place of ppt will make no difference to the
         resulting alpha value.
        
         
         By default function implements the formula given in
         Doonan, I.J.; Coombs, R.F.; McClatchie, S. (2003).
         The Absorption of sound in seawater in relation to
         estimation of deep-water fish biomass. ICES
         Journal of Marine Science 60: 1-9.
        
         Note that the paper has two errors in the formulae given
         in the conclusions.
        
         Optional argument 'method', allows the user to specify the
         Absorption formula to use. Default is Doonan etal (2003)
         other possibilities are:
         'doonan'  - Doonan et al (2003)
         'fandg'   - Francois & Garrison (1982)
         
         Based on and modified from Matlab code written  by 
         Gavin Macaulay, August 2003 and Yoann Ladroit 2015
         
        def get_method():
            if method == 'fandg':
                m = 'Francois and Garrison (1982)'
            elif method == 'Doonan':
                m = 'Doonan et al (2003)'
            return m
        
        print("Computing absorption, using %s" % get_method())
        if method == 'Doonan':
           #Check if temperature is higher then 20 or the frequency is above 120 kHz,
           #outside of Doonan's validity...IF so change method to Francois and Garrison
            if any(T) > 20 or f < 10 or f>120:
                method = 'fandg'
                print("Switching method to %s" % get_method())
        
            c = 1412 + 3.21 * T + 1.19 * S + 0.0167 * D
            A2 = 22.19 * S * (1.0 + 0.017 * T)
            f2 = 1.8 * 10**(7.0 - 1518/(T+273))
            #f2 = 1.8 * 10**7.0 * np.exp( -1818/(T+273.1))
            P2 = np.exp(-1.76e-4*D)
            A3 = 4.937e-4 - 2.59e-5*T + 9.11e-7*T*T - 1.5e-8*T*T*T
            P3 = 1.0 - 3.83e-5 * D + 4.9e-10 * D * D
            alpha = A2 * P2 * f2 * f * f / ( f2 * f2 + f * f ) / c + A3 * P3 * f * f
        
        elif method == 'fandg':
            alpha = alpha_FrancoisGarrison(f,T,S,D,pH=8)
        #print("Absorption is %.2f dB/km for %i kHz" %(alpha,f))
        return(alpha)
'''
		
    

def alpha_FrancoisGarrison(f,T,S,D,pH):
	'''calculation of absorption according to:
	Francois & Garrison, J. Acoust. Soc. Am., Vol. 72, No. 6, December 1982
	
	Parameters:
		
		f [int]: frequency (kHz)
		T [float]: Temperature (degC)
		S [float]: Salinity (ppt)
		D [float]: Depth (m)
		pH [float]: Acidity
	
	Returns:
		alpha [dB/km]
		
	'''

	# Total absorption = Boric Acid Contrib. + Magnesium Sulphate Contrib. + Pure Water Contrib.
	Kelvin = 273	# for converting to Kelvin (273.15)

	# Measured ambient temp
	T_kel = Kelvin + T

	# Calculate speed of sound (according to Francois & Garrison, JASA 72 (6) p1886)
	c = 1412 + 3.21 * T + 1.19 * S + 0.0167 * D

	# Boric acid contribution
	A1 = (8.86 / c ) * 10**(0.78 * pH - 5)
	P1 = 1
	f1 = 2.8 * np.sqrt(S / 35) * 10**(4 - 1245 / T_kel)
	Boric = (A1 * P1 * f1 * f**2)/(f**2 + f1**2)

	# MgSO4 contribution
	A2 = 21.44 * (S / c) * (1 + 0.025 * T)
	P2 = 1 - 1.37 * 10**(-4) * D + 6.2 * 10**(-9) * D**2
	f2 = (8.17 * 10**(8 - 1990/T_kel)) / (1 + 0.0018 * (S - 35))
	MgSO4 = (A2 * P2 * f2 * f**2)/(f**2 + f2**2)

	# Pure water contribution
	if T <= 20:
		A3 = 4.937 * 10**(-4) - 2.59 * 10**(-5) * T + 9.11 * 10**(-7) * T**2 - 1.50 * 10**(-8) * T**3
	else:
		A3 = 3.964 * 10**(-4) - 1.146 * 10**(-5) * T + 1.45 * 10**(-7) * T**2 - 6.50 * 10**(-10) * T**3
	P3 = 1 - 3.83 * 10**(-5) * D + 4.9 * 10**(-10) * D**2
	H2O = A3 * P3 * f**2

	# Total absorption
	alpha = Boric + MgSO4 + H2O

	return alpha
def tvg(f,s,t,d,R):
	alpha=absorption(f,s,t,d)/1000
	tvg = 20 * np.log10(R) + 2 * alpha * R
	return tvg
	