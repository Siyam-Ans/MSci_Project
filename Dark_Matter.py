#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 18:07:33 2017

@author: hannahbanks
"""
import numpy as np
import matplotlib.pyplot as plt
class DM:
    def __init__(self,r0,v0,m,cs):

        """create dark matter particle at initial position r0, with initial velocity v0
        mass m and mean free path l (in terms of rchi), distances in terms of DM scales"""

        self.r0=np.array(r0/r_chi)
        self.v0 = np.array(v0/r_chi)
        self.position=[(self.r0)]
        self.velocity=[(self.v0)]
        self.plot_points = [(self.r0)]
        self.m = m
        self.time = [0]
        self.cs = cs/(r_chi*r_chi)
        self.leftsun = False
        self.escaped = False
        self.position_leave = []
        self.position_reenter = []
        self.velocity_leave = []
        self.velocity_reenter = []
        self.leavetime = []
       

    def append(self,r,v,t):

        """append a new position and velocity to list """

        if isinstance(r,np.ndarray)== True:
            self.position.append(r)
        else:
            self.position.append(np.array(r))
            
        if isinstance(v,np.ndarray)==True:
            self.velocity.append(v)
        else:
            self.velocity.append(np.array(v))
        
        self.time.append(t)
    
    def append_plot(self,r):
        if isinstance(r,np.ndarray)==True:
            self.plot_points.append(r)
        else:
            self.plot_points.append(np.array(r))
    

    def r(self):

        """returns current position of DM particle"""
        
        return (self.position[-1])    


    def v(self):

        """returns current velocity of DM particle"""

        return (self.velocity[-1])
   
   
    def t(self):
        
        """returns current velocity of DM particle"""
        
        return (self.time[-1])
    
    """
    def velocities(self):
        
        return list of modulus of velocities of DM i.e. actually speeds (for the histogram), have inner and outer velocities are lists when DM is at two different temperatures
        
        velocities = []
        
        for i in range(0,len(self.velocity)):
            v = modulus(r_chi * self.velocity[i])
            velocities.append(v)
            if  modulus(self.position[i]) > 12:
                outervelocities.append(v)
            else:
                inner_velocities.append(v)
                
        return velocities

"""
    def points(self):
        
        """Return the list of positions of the DM particle"""
        
        return self.plot_points


    def plotDM(self):
        
        """2D plot of the projected DM path in the x,z plane"""
        
        points = self.points()
        xco_ords=[]
        zco_ords=[]        
        for i in points:
            xco_ords.append(i[0])
            zco_ords.append(i[2])
        print(xco_ords,zco_ords)
        plt.plot(xco_ords,zco_ords)