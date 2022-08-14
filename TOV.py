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
def rho(r,rho0,P):
    if(r<=R):
        rho=rho0
    else:
        rho=0
    return rho    
def Mass(r,rho0,Parr):
    n=int((r-r0)/dr)
    rint=np.linspace(r0,r,n) 
    M=0
    for i in range(0,n):
        M=4*np.pi*rint[i]**2*rho(rint[i],rho0,Parr[i])*dr+M 
    return M
def f(r,P,Parr):
    TOV=-G*(Mass(r,rho0,Parr)+4*np.pi*r**3*P/c**2)/(r**2*(1-2*G*Mass(r,rho0,Parr)/(r*c**2)))*(rho(r,rho0,P)+P/c**2)
    return TOV