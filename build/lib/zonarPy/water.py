# -*- coding: utf-8 -*-
'''
:author Sven Gastauer
:licence MIT
'''
class water():

    def compute_c(self,D, S, T, lat = None, method = 'mackenzie'):
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
    
    