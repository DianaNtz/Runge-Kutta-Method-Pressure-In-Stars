"""The code below was written by @author: https://github.com/DianaNtz and is a fourth order Runge Kutta implementation. 
It solves in particular the equation of hydrostatic equilibrium for a spherical symmetric star with polytropic equation of state 
and compares it with the analytical solution for polytropic index n=1."""
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
P=np.empty(steps+1, dtype='double')
r=np.empty(steps+1, dtype='double')
Pn=P0
rn=r0
#Runge Kutta fourth
for i in range(0,steps+1):
    r[i]=rn
    P[i]=Pn
    k1=dr*f(rn,Pn,P)
    k2=dr*f(rn+0.5*dr,Pn+0.5*k1,P)
    k3=dr*f(rn+0.5*dr,Pn+0.5*k2,P)
    k4=dr*f(rn+dr,Pn+k3,P)
    rn=rn+dr
    Pn=Pn+k2*(2/6)+k1*(1/6)+k3*(2/6)+k4*(1/6)
#analytical solution of the Lane Emden equation for n=1  
gamma=1/n+1
zeta=np.sqrt(4*np.pi*G*(gamma-1)/(k*gamma))*r*rho0**((2-gamma)/2)
theta=np.sin(zeta)/zeta
rhoa=theta**n*rho0
Pa=k*rhoa**gamma
ax1 = plt.subplots(1, sharex=True, figsize=(10,5))
plt.plot(r/1000,Pa,color='black',linestyle='-',linewidth=2.5,label = "$P_a(r)$ ")          
plt.plot(r/1000,P,color='deepskyblue',linestyle='-.',linewidth=2.5,label="$P(r)$ ")
plt.xlabel("distance $r$ in [km]",fontsize=19) 
plt.ylabel(r'pressure $P$ in [N/$m^2$]',fontsize=19,labelpad=10)
plt.xlim([r0/1000,Rfinal/1000])
plt.ylim([0,P0+P0*0.13])  
plt.xticks(fontsize= 17)
plt.yticks(fontsize= 17)
plt.legend(loc=1,fontsize=19,handlelength=3) 
plt.savefig("Lane.png",dpi=120)
plt.show()