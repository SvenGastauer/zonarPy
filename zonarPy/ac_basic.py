# -*- coding: utf-8 -*-
"""
Created on Mon May  4 12:48:05 2020

@author: sven
"""
import numpy as np

def nearfield(f,c,theta):
	"""
	Compute the nearfield

	Parameters
	----------
	f : numeric
		Transducer Frequency in kHz [kHz].
	c : numeric
		Ambient sound speed [m/s].
	theta : numeric
		3dB angle or beam width in degrees [degrees].

	Returns
	-------
	Rnf : numeric
		Range of the nearfield for the given conditions in meters [m].

	"""
	lmbd = c/ ( f * 1000)
	k = 2* np.pi / lmbd
	a = 1.6 / (k * np.sin((theta * np.pi/180) / 2))
	Rnf = (2*a)**2 / lmbd
	return Rnf

def eba(f,c,theta):
	"""
	Compute the equivalent beam angle for a circular transducer.

	Parameters
	----------
	f : numeric
		Transducer Frequency in kHz [kHz].
	c : numeric
		Ambient sound speed [m/s].
	theta : numeric
		3dB angle or beam width in degrees [degrees].

	Returns
	-------
	EBA : numeric
		equivalent beam angle in dB [dB].

	"""
	lmbd = c/ ( f * 1000)
	k = 2* np.pi / lmbd
	a = 1.6  / (k * np.sin((theta * np.pi/180) / 2))
	EBA = 10 * np.log10( 5.78 / ( ( k * a ) ** 2))#equivalent beam angle in steradians
	return EBA

def vol_samp(f,c,theta,tau,R,start=0):	
	f = f*1000
	Rtot = R+start
	Vtot = 10**(eba(f,c,theta)/10) * Rtot**2 * c * tau / 2
	V0 = 10**(eba(f,c,theta)/10) * start**2 * c * tau / 2
	V = Vtot - V0
	return V
def footprint_radius(theta,R):
	return R * np.tan(theta * np.pi / 180 / 2)

def footprint_area(theta, R):
	return np.pi * footprint_radius(theta,R)**2


'''
vol_samp(f=200,c=1450,theta=9.8,tau=6/1000,R=10)
vol_samp(f=1000,c=1450,theta=4,tau=6/1000,R=10)
#Zonar
nearfield(200,1480,9.8)

nearfield(1000,1480,4)

c=1450;f=200000
0.045**2/(c/f)

c=1450;f=1000000
0.022**2/(c/f)
'''