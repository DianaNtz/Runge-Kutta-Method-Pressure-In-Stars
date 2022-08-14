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
def rho(r,rho0,P):
    return ((1/k)*P)**(n/(n+1))   
def Mass(r,rho0,Parr):
    n=int((r-r0)/dr)
    rint=np.linspace(r0,r,n) 
    M=0
    for i in range(0,n):
        M=4*np.pi*rint[i]**2*rho(rint[i],rho0,Parr[i])*dr+M 
    return M
def f(r,P,Parr):   
    Lane=-G/r**2*Mass(r,rho0,Parr)*rho(r, rho0, P)
    return Lane