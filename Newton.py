"""
@author: https://github.com/DianaNtz
"""
import numpy as np
import matplotlib.pyplot as plt
dr=100 #in m
R=30000 #stars radius in m
Rfinal=R
r0=dr
steps=int(-(r0-Rfinal)/dr)
rho0=6*10**16 #central density of the star in kg/m^3
G=6.67259*10**-11 #gravitational constant in m^3/(kg s^2)
P0=4*np.pi/3*rho0*R**3*G/(2*R)*rho0 #central pressure (Newton)