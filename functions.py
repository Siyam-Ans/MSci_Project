#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 18:08:52 2017

@author: hannahbanks
"""
import numpy as np
import scipy as sp
import random as ran
from data import *

#git test

def modulus(array):
    """Modulus of a vector"""
    if isinstance(array,np.ndarray)== True:
        return np.sqrt(np.dot(array,array))
    else: 
        array = np.array(array)
        return np.sqrt(np.dot(array,array))

def transform_method_MB(T,m_N = 1.67e-27):
    
    """transformation method for the M-B distribution"""
    
    x = ran.uniform(0,1)
    y = (sp.special.erfinv(2*x - 1))*np.sqrt((k_B*T*2)/m_N)
    return y

def cart2polar(vector):
    r = np.sqrt(np.dot(vector,vector))
    phi = np.arctan(vector[1]/vector[0])
    theta = np.arctan(vector[2]/(np.sqrt((vector[0])**2+(vector[1])**2)))

    return np.array([r, theta, phi])
    

def polar2cart(vector):
    r = vector[0]
    theta = vector[1]
    phi = vector[2]
    x = r*np.cos(theta)*np.sin(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    
    return np.array(x,y,z)


def scale_variables(r_chi, vector):
    return vector/r_chi


def unscale_variables(r_chi, vector):
    return vector*r_chi

def MB(v,T=10000):
    
    """points on a Maxwell-Boltzmann curve for a given temperature and a velocity"""
    
    return ((mass_DM/(2*k_B*np.pi*T))**(3./2))*(4*np.pi*v*v)*np.exp((-1*mass_DM*v*v)/(2*k_B*T))  

 
def collisionpoint(r,v,d):
    
    
        
    """Determining the point of collision inside the box, collision points outside the box will be reflected inside"""
        
    #r is initial position, v is iniial velocity, d is distance to collision
    
    unit_v =v/modulus(v)
    
    collision_point = r + (d*unit_v)
    """    
    if collision[0]>box_size:
        collision[0] = 2*box_size - collision[0]
        self.v()[0] = - self.v()[0]
    if collision[1]>box_size:
        collision[1] = 2*box_size - collision[1]
        self.v()[1] = - self.v()[1]
    if collision[2]>box_size:
        collision[2] = 2*box_size - collision[2]
        self.v()[2] = - self.v()[2]
        """
    return collision_point

def collision(w, v_n, theta,phi, m_N = 1.67e-27, m_DM = 8.9e-27 ):
    
    """Dynamics n to find the velocity of the DM particle after the collision in the lab frame"""
    #w is DM speed before collision, v_n is nucleus speed, m_N and m_DM are the masses of the nucleus and DM respectively 
    #theta and phi are picked out of a hat
    
    s = (m_DM*w + m_N*v_n)/(m_DM+m_N) #COM velocity
    t = w - s  #Velocity of incoming DM in centre of mass fram
    ty = t[1]
    tz = t[2]
    ratio = -tz/ty
    vect = [0,ratio,1]
    
    p_unit = vect/modulus(vect)
                
    t_unit = t/modulus(t)
    
    q_unit = np.cross(p_unit , t_unit)
    
    t_prime = modulus(t)*np.sin(theta)*np.cos(phi)*p_unit + \
    modulus(t)*np.sin(theta)*np.sin(phi)*q_unit - modulus(t)*np.cos(theta)*t_unit
    
    v_DM = s + t_prime
    
    return v_DM

def mfp(DM, density):
    """takes in an unscaled density and works out the mfp in m"""
    number_density = density/(1.67e-27)
    mfp = 1/(DM.cs*r_chi*r_chi*number_density)
    return mfp
        
def transform_method_distance(l_chi):
    
    """transformation method for mfp pdf"""
    
    x = ran.uniform(0,1)
    y = -1*l_chi*np.log(1-x)  
    return y/r_chi
 
def intercept(DM):
    pos = DM.r()
    
    unit_direction = DM.v()/modulus(DM.v())
    
    dot = np.dot(unit_direction,pos)
    
    
    
    
    discriminent = (dot*dot)-np.dot(pos,pos)+(boundary*boundary)
    length_1 = -dot + np.sqrt(discriminent)
    length_2 = -dot-np.sqrt(discriminent)
    point_1 = DM.r()+length_1*unit_direction
    point_2 = DM.r()+length_2*unit_direction
    if modulus(point_1)==boundary:
        return point_1,length_1
    else:
        return point_2,length_2





def reflection(DM,travel):
    
    boundary = Solar_Radius/r_chi
    direction = DM.v()/modulus(DM.v())
    interception,length = intercept(DM)
        
    normal = (-1*interception)/modulus(interception)
    speed = modulus(DM.v())
    new_direction = direction - 2*normal*np.dot(normal,direction)
    excess = travel - length
    new_direction_n = new_direction/modulus(new_direction)
    velocity = speed*new_direction_n
    final_point = collisionpoint(interception,velocity,excess)
    
     #need to put correct time in
    return final_point, velocity, interception


    
