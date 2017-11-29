#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 16:13:08 2017

@author: hannahbanks
"""
import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
from functions import *
import scipy as sp


def gradient_Kepler(x,y,phi_0,e):
    A = np.cos(phi_0)
    B = np.sin(phi_0)
    dy = -(x+e*A*np.sqrt(x**2+y**2))
    dx = y + B*e*np.sqrt(x**2+y**2)
    
    """
    dy = (2*x*np.sqrt(x**2+y**2)*(1+e*e*A*A) + 2*e*A*(2*x*x + y*y) + 2*e*B*y*x + 2*A*B*y*e*e*np.sqrt(x*x + y*y))
    dx = -(2*y*np.sqrt(x**2+y**2)*(1+e*e*B*B) + 2*e*B*(x*x + 2*y*y) + 2*e*A*y*x + 2*A*B*x*e*e)
    """
    return dy/dx
    
    
    

def new_point_func(position, velocity,DM):
    #position and velocity at the point of exit
    speed = modulus(velocity)
    radial_dist = modulus(position)
    E = 0.5*DM.m*speed*speed - (G*mass_sun*DM.m*(radial_dist**-1)*(r_chi**-3))
    

   
    if E > 0:
        print("ESCAPED")
        return None
    
    r_hat = position/modulus(position)
    
    x_0 = np.dot(position,r_hat)
    
    perp = velocity - (np.dot(r_hat,velocity))*r_hat
    
    perp_hat = perp/modulus(perp)
    
    y_0 = np.dot(position,perp_hat)
    
    
    vel_phi = np.dot(velocity,perp_hat)
    vy_0 = vel_phi
    vx_0 = np.dot(velocity,r_hat)
    
    L = DM.m*radial_dist*vel_phi
    
    r_0 = L**2/(DM.m*DM.m*mass_sun*G*r_chi**-3)
    
   
    
    
    eps_square = 1 + (2*E*r_0)/(G*DM.m*mass_sun*r_chi**-3)
    
    eps = np.sqrt(eps_square)
    
    if vx_0 > 0:
         phi_0 = 2*np.pi - np.arccos((r_0 - radial_dist)/(eps*radial_dist))
    else: 
        
        phi_0 = np.arccos((r_0 - radial_dist)/(eps*radial_dist))
  
    
    
    
    
    
    
    """
    cos_psi = abs(np.dot(velocity,r_hat))/speed
    psi = np.arccos(cos_psi)
    sin_angle = (1/np.tan(psi))*r_0/(radial_dist*eps)
    phi_0 = -1*np.arcsin(sin_angle) 
    """
    
    
    
    
    
       
    phiplot = np.arange(0,2*np.pi,0.001)
    x = []
    y = []
    
    for i in phiplot:
       rplot = r_0/(1 + (eps*np.cos(i - phi_0)))
       pointgraph = rplot*np.cos(i)*r_hat + rplot*np.sin(i)*perp_hat
       x.append(pointgraph[0])
       y.append(pointgraph[1])
       
    plt.plot(x,y)
    
    circle1 =plt.Circle((0,0),Solar_Radius/r_chi,color='r')
    fig=plt.gcf()
    """
    cos_phi_phi_0 = (1-(r_0/radial_dist) )/eps
    phi = -np.arccos(cos_phi_phi_0)+ phi_0
    """
   
    
    phi_1 = np.arccos((r_0 - radial_dist)/(radial_dist*eps)) + phi_0
    
    phi_2 = 2*np.pi - np.arccos((r_0 - radial_dist)/(radial_dist*eps)) + phi_0
    
    if phi_1 % (2*np.pi) == 0:
        phi = phi_2
    else:
        phi = phi_1
    
    fig.gca().add_artist(circle1)
    #these are the components in the 2D perpendicular basis
    x_new = (Solar_Radius/r_chi)*np.cos(phi)
    y_new = (Solar_Radius/r_chi)*np.sin(phi)
    
    new_position =   x_new*r_hat + y_new*perp_hat
    
    plt.plot(Solar_Radius/r_chi,0, 'mx', markersize=20,label='Initial Position',markeredgewidth='3')
    plt.plot(new_position[0],new_position[1], 'mx', markersize=20,label='Initial Position',markeredgewidth='3')
    
    #new velocity
    grad_0 = gradient_Kepler(x_0,y_0,phi_0,eps)
    grad = gradient_Kepler(x_new,y_new,phi_0,eps)
    
    
    vx_new_square = (speed**2)/(1 + grad*grad)
    vy_new_square = (speed**2) - vx_new_square
   
    
    if grad_0*grad >0:
        if vx_0>0:
            vx_new = -np.sqrt(vx_new_square)
        else:
            vx_new = np.sqrt(vx_new_square)
        if vy_0 >0:
            vy_new =  - np.sqrt(vy_new_square) 
        else:
            vy_new =  np.sqrt(vy_new_square) 
    else:
         
        vx_new = - np.sqrt(vx_new_square)
        if vy_0 >0:
            vy_new =  np.sqrt(vy_new_square) 
        else:
            vy_new =  - np.sqrt(vy_new_square)
  
        
    
    """
    if vx_0 > 0:
        vx_new = -np.sqrt(vx_new_square)
    else:
        vx_new = np.sqrt(vx_new_square)
    if vy_0 > 0:
        vy_new =  -np.sqrt(vy_new_square) 
    else: 
        vy_new = +np.sqrt(vy_new_square) 
    """
    new_velocity = (vx_new*r_hat) + (vy_new*perp_hat)

    
    
    return [new_position,new_velocity]

def f(z,t):
    xx = z[0]
    vx = z[1]
    yy = z[2]
    vy = z[3]
    zz = z[4]
    vz = z[5]
    
    r = ((xx**2)+(yy**2) +(zz**2))**0.5
    ax= -G*(r_chi**-3)*mass_sun*xx*r**(-3.)
    ay= -G*(r_chi**-3)*mass_sun*yy*r**(-3.)
    az = -G*(r_chi**-3)*mass_sun*zz*r**(-3)
    
    return [vx,ax,vy,ay,vz,az]




t = sp.linspace(0,60000,10000)

def trajectory(x,y,vx,vy):
    initial = [x,vx,y,vy,0,0]
    soln = spi.odeint(f,initial,t)
    x = soln[:,0]
    vx = soln[:,1]
    y = soln[:,2]
    vy = soln[:,3]
    plt.figure()
    plt.plot(x,y)
    plt.xlabel("Horizontal Displacement ")
    plt.ylabel("Vertical Displacement ")
   
    
    plt.plot(boundary,0, 'mx', markersize=20,label='Initial Position',markeredgewidth='3')
    
    
    circle1 =plt.Circle((0,0),boundary,color='r') # adding mars to the plot
    fig=plt.gcf()
    fig.gca().add_artist(circle1)
    plt.legend(numpoints=1,loc='upper left')
