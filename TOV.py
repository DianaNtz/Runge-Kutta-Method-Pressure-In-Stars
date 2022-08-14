"""
@author: https://github.com/DianaNtz
"""
import numpy as np
import matplotlib.pyplot as plt
dr=100 #in m
R=30000 #star radius in m
Rfinal=R
r0=dr
steps=int(-(r0-Rfinal)/dr)
c=299792458 #speed of light in m/s
rho0=6*10**16 #central density of the star in kg/m^3
G=6.67259*10**-11 #gravitational constant in m^3/(kg s^2)
rs=2*G/c**2*(4*np.pi/3*rho0*R**3) #schwarzschild radius
P0=rho0*c**2*(1-np.sqrt(1-rs/R))/(3*np.sqrt(1-rs/R)-1) #central pressure (TOV)
