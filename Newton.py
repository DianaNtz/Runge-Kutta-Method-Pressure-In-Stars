"""The code below was written by @author: https://github.com/DianaNtz and is a fourth order Runge Kutta implementation. 
It calculates in particular the pressure radius relation of a spherical symmetric Newtonian star with constant density 
and compares it with its analytical solution."""
import numpy as np
import matplotlib.pyplot as plt
dr=100 #in m
R=30000 #stars radius in m
Rfinal=R
r0=dr
steps=int(-(r0-Rfinal)/dr)
rho0=6*10**16 #central density of the star in kg/m^3
G=6.67259*10**-11 #gravitational constant in m^3/(kg s^2)
P0=4*np.pi/3*rho0*R**3*G/(2*R)*rho0 #central pressure 
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
Pa= 4*np.pi/3*rho0*R**3*G/(2*R)*rho0*(1-(r/R)**2)
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
plt.savefig("Newton.png",dpi=120)
plt.show()
