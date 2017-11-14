#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 11:23:06 2017

@author: hannahbanks
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

class DM:
    def __init__(self,r0,v0,m,tcs):

        """create dark matter particle"""

        self.r0=np.array(r0)
        self.v0 = np.array(v0)
        self.position=[(self.r0)]
        self.velocity=[(self.v0)] 
       

    def append(self,r,v):

        """append a new position direction to list """

        if isinstance(r,np.ndarray)== True:
            self.position.append(r)
        else:
            self.position.append(np.array(r))

        if isinstance(v,np.ndarray)==True:
            self.direction.append(v)
        else:
            self.position.append(np.array(v))

    def r(self):

        """get current position of ray"""
        return (self.position[-1])    

    def v(self):

        """get current velocity of ray"""

        return (self.direction[-1])

    def points(self):
        return self.position


    def plotDM(self):
        #just a 2D plot for now
        points = self.points()
        xco_ords=[]
        
        zco_ords=[]        
        for i in points:
            xco_ords.append(i[0])
            zco_ords.append(i[2])
            plt.plot(zco_ords,xco_ords)
            
def collisionpoint(r,v,d):
     mod_v = np.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
     unit_v = v/mod_v
     return r + (d*unit_v)
 
    
        
    