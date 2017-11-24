#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 12:01:09 2017

@author: hannahbanks
"""

import numpy as np
import random as ran
from functions import *
from ellipse import *
from data import *




def func(DM,t):
    evaluation_point = get_point_speed(DM.r(),DM.v(),t)[0]
    distance = modulus(evaluation_point)
    result = distance - (Solar_Radius/r_chi)
    return result


def intersection_ellipse(DM,t0,t1):
    # t0 is a point you know is in the sun
    #t1 is outsoide the sun
    print("HELP")
    i=0
    if func(DM,t0)*func(DM,t1) > 0:
        print("No ROOT")
        return
    else:
        while abs((t1-t0)/2) > 0.00001:
            midpoint = (t1+t0)/2
            if func(DM,midpoint)==0:
                return midpoint
            elif func(DM,t0)*func(DM,midpoint) <0:
                t1 = midpoint
            else:
                t0 = midpoint
        print(func(DM,midpoint))
        return midpoint

"""             
    evaluation_point = get_point_speed(DM.r(),DM.v(),t1)[0]
    
    distance = modulus(evaluation_point)
    test = distance - (Solar_Radius/r_chi)
   
    t2 = (t1+t0)/2
    
    while  abs(test) > 0.001*Solar_Radius/r_chi and i<60:
        
        
        evaluation_point = get_point_speed(DM.r(),DM.v(),t2)[0]
        i = i+1
        distance =(modulus(evaluation_point))
        print(distance)
        
        test = distance - (Solar_Radius/r_chi)
        
        if test > 0:
            
            t2 = (t2+t0)/2
        else: 
            t2 = (t2+t1)/2
        
    time = t2
    print("THe number of iteractions was: ",i)
    return time
"""
        
    
    
    
def get_point_speed(r,v,t):
    a = np.arctan(-v[0]/(r[0]*w))

    b = np.arctan(-v[1]/(r[1]*w))
    c = np.arctan(-v[2]/(r[2]*w))
    A = r[0]/np.cos(a)
    B = r[1]/np.cos(b)
    C = r[2]/np.cos(c)
    x = A*np.cos(w*t+a)
    y = B*np.cos(w*t+b)
    z = C*np.cos(w*t+c)
    vx = -A*w*np.sin((w*t)+a)
    vy = -B*w*np.sin((w*t)+b)
    vz = -C*w*np.sin((w*t)+c)
    return np.array([x,y,z]), np.array([vx,vy,vz])

def integrand(DM,t_0,t):
    cross_section = DM.cs
    
    evaluation_point = get_point_speed(DM.r(),DM.v(),t)[0]
    distance = modulus(evaluation_point)
    print("start",distance)
    
    if DM.leftsun == True:
        
        evaluation_point = get_point_speed(DM.position_reenter,DM.velocity_reenter,t - DM.leavetime)[0]
        distance = modulus(evaluation_point)
        print(distance)
    
    elif distance > (boundary):
        DM.leftsun = True
        time = intersection_ellipse(DM,0,t)
        DM.leavetime = time
        leave_velocity = get_point_speed(DM.r(),DM.v(),time)[1]
        
        if modulus(leave_velocity)> np.sqrt(2*G*(r_chi**-3)*mass_sun/(Solar_Radius/r_chi)):
            
            print(ESCAPED)
            return
        #point is the point where it leaves the sun
        leave_point = get_point_speed(DM.r(),DM.v(),time)[0]
        DM.position_leave = leave_point
        #work out where it comes back in using the KEPLER ORBIT
        new_point,new_velocity = new_point_func(leave_point,leave_velocity,DM)
        DM.position_reenter = new_point
        DM.velocity_reenter = new_velocity
        
        evaluation_point = get_point_speed(DM.position_reenter, DM.velocity_reenter, t - DM.leavetime)[0]
        distance = modulus(evaluation_point)
        print(distance)
    else:
        print(distance)
        None
        
        
        
    print(DM.leftsun,"end",distance)
    density = get_parameters(distance)[1]
    speed = modulus(get_point_speed(DM.r(),DM.v(),t)[1])
    number_density = density/1.67e-27
    
    f = number_density*cross_section*speed
    
    return f

  


def RK45_time(DM, integrand, h_init, tau, atol,rtol, minscale,maxscale):
    h = h_init
    
    
    safety_factor = 0.9
    alpha = 0.14
    beta =  0.08
    errold = 1e-4
    t = [0]
    tau_d = [0]
    
    
    a = np.arctan(-DM.v()[0]/(DM.r()[0]*w))
    b = np.arctan(-DM.v()[1]/(DM.r()[1]*w))
    c = np.arctan(-DM.v()[2]/(DM.r()[2]*w))
    A = DM.r()[0]/np.cos(a)
    B = DM.r()[1]/np.cos(b)
    C = DM.r()[2]/np.cos(c)
    
    
    
    
    while tau_d[-1] < tau:
        """
        if r[-1]+h >b:
            h = b-r[-1] 
            
       """
        
       
        k1 = (h*integrand(DM,0,t[-1]))
        k2 = (h*integrand(DM,0,t[-1]+0.25*h))
        k3 = (h*integrand(DM,0,t[-1]+h*(3/8)))
        k6 = (h*integrand(DM,0,t[-1]+0.5*h))
        k4 = (h*integrand(DM,0,t[-1]+(12/13)*h))
        
        k5 = (h*integrand(DM,0,t[-1]+h))
        
        
        
        
        
        
        
        fourth_order = tau_d[-1]+ (k1*(5179/57600)) + (k3*(7571/16695)) + (k4*(393/640)) + (k5*(-92097/339200)) + (k6*(187/2100))
        new_tau = tau_d[-1] + (35./384)*k1 + (500./1113)*k3 + (125./192)*k4 + (-2187./6784)*k5 + (11./84)*k6
        
        delta = abs(new_tau - fourth_order)
        sk = atol + rtol*max(abs(new_tau), abs(tau_d[-1]))
        
        
        err = (delta)/sk
        if err == 0:
            scale = maxscale
        else:
            scale = safety_factor*(err**(-1*alpha))*(errold**beta)
            if scale > maxscale:
                scale = maxscale
            elif scale < minscale:
                scale = minscale
        
        
        if err <= 1:
    
            t.append(t[-1]+h)
            tau_d.append(new_tau)
            h =h*scale
            errold = max(err,1e-4)
           
                
        else:
            print("rk45error")
            h = h*scale
          
    tau_difference = tau_d[-1]-tau_d[-2]
    fraction = (tau_d[-1]- tau)/tau_difference
    t_difference = t[-1]-t[-2]
    required = t[-1] - (fraction*t_difference)
    
      
    return required
 
           
            