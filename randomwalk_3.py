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

def random_walk3(DMparticle,iterations):
    l_chi = mfp(DMparticle, rho_c) 
    
    
    
    """DM particle position and velocity after each collision to trace out the random walk"""
    
    for i in range(iterations):
        DM.leftsun = False
        
        extinction = ran.uniform(0,1)
        tau = -1*np.log(extinction) 
        
        time = RK45_time(chi,integrand,0.00000005,tau,0.01,0.01,0.2,10)
        
        
        if DMparticle.leftsun == False:
            collision_point = get_point_speed(DMparticle.r(),DMparticle.v(),time)[0]
            print("yes")
        else:
            collision_point = get_point_speed(DMparticle.position_reenter,DMparticle.velocity_reenter,time - DMparticle.leavetime)[0]
            print("no")
        
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
        
        
        
        DMparticle.append(collision_point,new_velocity,10)
        DMparticle.append_plot(collision_point)
        print(i)
    print('K is',l_chi/r_chi)
        
    return