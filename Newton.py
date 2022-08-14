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
    Newton=-G*Mass(r,rho0,Parr)/r**2*rho(r,rho0,P)
    return Newton