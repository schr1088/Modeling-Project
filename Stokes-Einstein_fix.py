# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 14:19:54 2017

@author: Peter
"""

#Modeling of Stokes-Einstein equation,diffusion of spherical particles through a medium w low reynolds number
#D=kb*T/6*pi*eta*r


from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
from pylab import linspace

pi = 3.14 
K = 1.38**-23  #boltzmann constant m^2kgs^-2K^-1
T = 278.15 #K absolute temperature
nu = 1.5182 #mPa*s eta dynamic viscosity of water @ 5C
R = np.arange(0.1,1,0.05) #radius of spherical particle after 1 m Reynolds num too high
D = (K*T/(6*pi*nu*R)) #diffusion constant
dr = np.diff(R) #change in radius
t_final = 10.
nx = 41
dx = 2 / (nx -1)
sigma = .2
dt = sigma * dx**2 / nu
nt = 25

#values from change in radius
D_06 = 0.01
D_03 = 0.02
D_02 = 0.03
D_01 = 0.06

#Df=np.ones((D,R),dtype=object)

m = np.ones(nx) #array
m[int(0.5 / dx):int(1 / dx +1 )] = 2

um = np.ones(nx)

for n in range(nt):
    um = m.copy()
    for i in range(1, nx -1):
        m[i] = um[i] + nu *dt/dx**2 *(um[i+1] - 2*um[i]+um[i-1])


m_06= np.ones(nx)
m_06[int(0.5/ dx):int(1 / dx +1)] = 2

for n in range(nt):
    um = m_06.copy()
    for i in range(1, nx -1):
        m_06[i] = um[i] + D_06*nu *dt/dx**2 *(um[i+1] - 2*um[i]+um[i-1])

m_03=np.ones(nx)
m_03[int(0.5/ dx):int(1 / dx +1)] = 2

for n in range(nt):
    um = m_03.copy()
    for i in range(1, nx -1):
        m_03[i] = um[i] + D_03*nu *dt/dx**2 *(um[i+1] - 2*um[i]+um[i-1])

m_02=np.ones(nx)
m_02[int(0.5/ dx):int(1 / dx +1)] = 2

for n in range(nt):
    um = m_02.copy()
    for i in range(1, nx -1):
        m_02[i] = um[i] + D_02*nu *dt/dx**2 *(um[i+1] - 2*um[i]+um[i-1])

m_01=np.ones(nx)
m_01[int(0.5/ dx):int(1 / dx +1)] = 2

for n in range(nt):
    um = m_01.copy()
    for i in range(1, nx -1):
        m_01[i] = um[i] + D_01*nu *dt/dx**2 *(um[i+1] - 2*um[i]+um[i-1])


#Plotting    

    
fig = plt.figure(1)
plt.title('Stokes-Einstein Diffusion Constant')
plt.xlabel('Diffusion Constant')
plt.ylabel('Radius of Particle (m)')
plt.grid(True)
plt.yscale('linear')
plt.xscale('linear')
plt.plot(D, R, 'b-', linewidth=4)

fig2 = plt.figure()
plt.title('Stokes-Einstein Diffusion Constant Logarithmic Scale')
plt.xlabel('Diffusion Constant')
plt.ylabel('Radius of Particle (m)')
plt.grid(True)
plt.yscale('log')
plt.xscale('log')
plt.plot(D, R, 'r-', linewidth=4)

fig3 = plt.figure(3)
plt.title('1-D Diffusion without Stokes-Einstein Coefficient')
plt.ylabel('Concentration')
plt.xlabel('Time')
plt.plot(linspace(0,2,nx),m, 'r')

fig4 = plt.figure(4)
plt.title('1-D Diffusion with Stokes-Einsteins Coefficient')
plt.ylabel('Concentration')
plt.xlabel('Time')
plt.plot(linspace(0,2,nx),m, 'r',m_06, 'b',m_03,'k',m_02, 'm',m_01,'g');

plt.show()