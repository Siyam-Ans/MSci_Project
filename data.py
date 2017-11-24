#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 06:41:56 2017

@author: hannahbanks
"""
import numpy as np

k_B = 1.38064852e-23
T_c = 1.544e7  #Temperature of the core 
G = 6.67e-11
rho_c = 1.489e5 #Mass denisty of the core 
mass_DM = 8.9e-27 #DM mass= 5GeV
r_chi = np.sqrt((3*k_B*T_c)/(2*np.pi*G*rho_c*mass_DM)) #The DM cloud Scale
av_den = 1410
mass_sun = 1.989e30

w_square = G*av_den*np.pi*(4/3)
w = np.sqrt(w_square)


Solar_Radius = 695700000
boundary = Solar_Radius/r_chi



radius,temp,rho = np.loadtxt("data.txt", skiprows = 9, usecols=(1,2,3), unpack = 1)

radius = (radius*Solar_Radius)/r_chi
rho = rho*1000*(r_chi)**3

def get_parameters(distance):
    i = 0
    while distance > radius[i]:
        i +=1
    if i != 0:
        i -=1
    return [temp[i],rho[i]]



