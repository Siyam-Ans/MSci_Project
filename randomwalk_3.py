#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 20:56:33 2017

@author: hannahbanks
"""
from functions import *
from RK45_3 import *
from data import *
from ellipse import *
import numpy as np

def get_interpolated_time(DM,time):
    if time > DM.leavetime[-1]:
        print("YES")
        i = len(DM.leavetime)  - 1
        print(i)
    else:
        i = 0
        while time > DM.leavetime[i]:
            i +=1
        
        if i != 0:
            i -=1
    
    return i

def random_walk3(DMparticle,iterations):
    l_chi = mfp(DMparticle, rho_c) 
    
    
    
    """DM particle position and velocity after each collision to trace out the random walk"""
    
    for i in range(iterations):
        
        
        extinction = ran.uniform(0,1)
        tau = -1*np.log(extinction) 
        print("tau is", tau)
        time = RK45_time(chi,integrand,0.00000005,tau,0.01,0.01,0.2,10,1000)
        print(time)
        
        if DMparticle.leftsun == False:
            collision_point = get_point_speed(DMparticle.r(),DMparticle.v(),time)[0]
            distance_from_origin = modulus(collision_point)
            
        else:
           
            index = get_interpolated_time(DMparticle,time)
            print(index)
            collision_point = get_point_speed(DMparticle.position_reenter[index],DMparticle.velocity_reenter[index],time - DMparticle.leavetime[index])[0]
            distance_from_origin = modulus(collision_point)
            print(distance_from_origin)
        
       
      
       
        T = get_parameters(distance_from_origin)[0]
        
        
        
        v_nx = (transform_method_MB(T))/r_chi
        v_ny = (transform_method_MB(T))/r_chi
        v_nz = (transform_method_MB(T))/r_chi
        #The velocity of the nuclei that DM is going to collide with
        v_n = np.array([v_nx,v_ny,v_nz])
        
        cos_theta = ran.uniform(-1,1)
        theta = np.arccos(cos_theta)
        
        phi = (ran.uniform(0,2*np.pi))
        new_velocity = collision(DMparticle.v(),v_n,theta,phi)
        
        
        
        DMparticle.append(collision_point,new_velocity,time)
        DMparticle.append_plot(collision_point)
        DMparticle.leftsun = False
        DMparticle.leavetime = []
        DMparticle.position_leave = []
        DMparticle.position_reenter = []
        DMparticle.velocity_leave = []
        DMparticle.velocity_reenter = []
        print(i)
    print('K is',l_chi/r_chi)
        
    return