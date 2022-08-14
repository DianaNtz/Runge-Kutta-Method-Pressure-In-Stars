"""
@author: https://github.com/DianaNtz
"""
import numpy as np
import matplotlib.pyplot as plt
dr=100 #in m
R=210000 #in m
Rfinal=R
r0=dr
steps=int(-(r0-Rfinal)/dr)
rho0=2*10**16 #central density of the star in kg/m^3
G=6.67259*10**-11 #gravitational constant in m^3/(kg s^2)
n=1 #dimensionless
k=2 #m^5/kg/s^2 for n=1
P0=k*rho0**(1+1/n) #central pressure 