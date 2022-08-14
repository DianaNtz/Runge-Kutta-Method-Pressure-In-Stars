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
P0=rho0*c**2*(1-np.sqrt(1-rs/R))/(3*np.sqrt(1-rs/R)-1) #central pressure 
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
#analytical solution for constant density star
Pa=rho0*c**2*(np.sqrt(1-r**2*rs/R**3)-np.sqrt(1-rs/R))/(3*np.sqrt(1-rs/R)-np.sqrt(1-r**2*rs/R**3))
ax1 = plt.subplots(1, sharex=True, figsize=(10,5))
plt.plot(r,Pa,color='black',linestyle='-',linewidth=2.5,label = "$P_a(r)$ ")          
plt.plot(r,P,color='deepskyblue',linestyle='-.',linewidth=2.5,label="$P(r)$ ")
plt.xlabel("distance $r$ in [m]",fontsize=19) 
plt.ylabel(r'pressure $P$ in [N/$m^2$]',fontsize=19,labelpad=10)
plt.xlim([r0,Rfinal])
plt.ylim([0,P0+P0*0.13])  
plt.xticks(fontsize= 17)
plt.yticks(fontsize= 17)
plt.legend(loc=1,fontsize=19,handlelength=3) 
plt.savefig("TOV.png",dpi=250)
plt.show()