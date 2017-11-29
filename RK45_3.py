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
import scipy as sp




def func(t,r,v,DM):
    evaluation_point = get_point_speed(r,v,t)[0]
    
    distance = modulus(evaluation_point)
    
    result = distance - (Solar_Radius/r_chi)
    return result

def intersection_ellipse(t0,t1,r,v,DM):
    print(func(t0,r,v,DM))
    print(func(t1,r,v,DM))
    #using Brents method
    root = sp.optimize.brentq(func,t0,t1,args=(r,v,DM))
    
    return root


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

def plot_SHO(r,v):
    t = sp.linspace(0,60000,10000)
    
    a = np.arctan(-v[0]/(r[0]*w))

    b = np.arctan(-v[1]/(r[1]*w))
    c = np.arctan(-v[2]/(r[2]*w))
    A = r[0]/np.cos(a)
    B = r[1]/np.cos(b)
    C = r[2]/np.cos(c)
    x = []
    y=[]
    for i in t:
        x_i = A*np.cos(w*i+a)
        y_i = B*np.cos(w*i+b)
        x.append(x_i)
        y.append(y_i)
        
    
  
    plt.plot(x,y,'b')
    return 
    
def integrand(DM,t_0,t):
    #t is the time since the LAST COLLISION i.e. local time along this trajectory
    cross_section = DM.cs
    
    evaluation_point = get_point_speed(DM.r(),DM.v(),t)[0]
    
    distance = modulus(evaluation_point)
    
    
    
    if DM.leftsun == True:
        
        evaluation_point = get_point_speed(DM.position_reenter[-1],DM.velocity_reenter[-1],t - DM.leavetime[-1])[0]
        distance = modulus(evaluation_point)
        
        
        while distance > boundary:
            
            time = intersection_ellipse(0.01,t - DM.leavetime[-1],DM.position_reenter[-1],DM.velocity_reenter[-1],DM)
          
            leave_point,leave_velocity = get_point_speed(DM.position_reenter[-1],DM.velocity_reenter[-1],time)
            DM.position_leave.append(leave_point)
            DM.velocity_leave.append(leave_velocity)
            new_point,new_velocity = new_point_func(leave_point,leave_velocity,DM)
            DM.position_reenter.append(new_point)
            DM.velocity_reenter.append(new_velocity)
            
            new_leave_time = (DM.leavetime[-1]) + time
            DM.leavetime.append(new_leave_time)
           
            evaluation_point = get_point_speed(DM.position_reenter[-1], DM.velocity_reenter[-1], t - DM.leavetime[-1])[0]
            distance = modulus(evaluation_point)
            print(distance)
        
      
        
    
    elif distance > (boundary):
        print("LEFT SUN")
        DM.leftsun = True
        time = intersection_ellipse(t_0,t,DM.r(),DM.v(),DM)
        DM.leavetime.append(time)
        leave_point,leave_velocity = get_point_speed(DM.r(),DM.v(),time)
        
        
        if modulus(leave_velocity)> np.sqrt(2*G*(r_chi**-3)*mass_sun/(Solar_Radius/r_chi)):
            
            print(ESCAPED)
            return
        #point is the point where it leaves the sun
        
        DM.position_leave.append(leave_point)
        DM.velocity_leave.append(leave_velocity)
        #work out where it comes back in using the KEPLER ORBIT
        new_point,new_velocity = new_point_func(leave_point,leave_velocity,DM)
     
        DM.position_reenter.append(new_point)
        DM.velocity_reenter.append(new_velocity)
        
        plot_SHO(DM.position_reenter[-1],DM.velocity_reenter[-1])
        
        evaluation_point = get_point_speed(DM.position_reenter[-1], DM.velocity_reenter[-1], t - DM.leavetime[-1])[0]
        distance = modulus(evaluation_point)
        
    else:
        
        None
        
        
        
    print(distance)
    density = get_parameters(distance)[1]
    speed = modulus(get_point_speed(DM.r(),DM.v(),t)[1])
    number_density = density/1.67e-27
    
    f = number_density*cross_section*speed
    
    return f

  


def RK45_time(DM, integrand, h_init, tau, atol,rtol, minscale,maxscale,absmax):
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
        index_0 = len(DM.leavetime)
       
        k1 = (h*integrand(DM,0,t[-1]))
        k2 = (h*integrand(DM,t[-1],t[-1]+0.25*h))
        k3 = (h*integrand(DM,t[-1]+0.25*h,t[-1]+h*(3/8)))
        k6 = (h*integrand(DM,t[-1]+h*(3/8),t[-1]+0.5*h))
        k4 = (h*integrand(DM,t[-1]+(0.5*h), t[-1]+(12/13)*h))
        
        k5 = (h*integrand(DM,t[-1]+(12/13*h),t[-1]+h))
        
        
        
        
        
        
        
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
       
            
        print("scale",scale)
        if err <= 1:
    
            t.append(t[-1]+h)
            tau_d.append(new_tau)
            if h*scale<absmax:
                h =h*scale
            else:
                h = absmax
            errold = max(err,1e-4)
            
            print("tau",tau_d[-1],h)
           
                
        else:
            DM.leavetime = DM.leavetime[:index_0]
            DM.position_leave = DM.position_leave[:index_0]
            DM.position_reenter = DM.position_reenter[:index_0]
            DM.velocity_leave = DM.velocity_leave[:index_0]
            DM.velocity_reenter = DM.velocity_reenter[:index_0]
            print("rk45error")
            h = h*scale
          
    tau_difference = tau_d[-1]- tau_d[-2]
    fraction = (tau_d[-1]- tau)/tau_difference
    t_difference = t[-1]-t[-2]
    required = t[-1] - (fraction*t_difference)
    print(required,t[-1],t[-2])
      
    return required
 
           
            