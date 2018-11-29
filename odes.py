#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###
# Name: Jacob Anabi,
# Student ID: 2294644,
# Email: anabi@chapman.edu,
# Course: PHYS220/MATH220/CPSC220 Fall 2018
# Assignment: CW11
###

import numpy as np
from matplotlib import pyplot as plt

def gen_method(r0, N):
    t = np.linspace(0, 5*(2*np.pi), 5*N, endpoint=True) # domain points
    delta_t = (2*np.pi)/N # delta t
    r = np.zeros((5*N, 2)) # 5Nx2 array
    r[0] = r0 # initial r-value
    J = np.array(([0,1],[-1,0])) # antisymmetric matrix (meaning its transpose equals its negative)
    I = np.array(([1,0],[0,1])) # identity matrix
    for i in range(1, 5*N):
        r[i] = (I*np.cos(t[i]) + J*np.sin(t[i]))@(r0)
    return r

def euler_method(r0, N):
    delta_t = (2*np.pi)/N # delta t
    r = np.zeros((5*N, 2)) # 5Nx2 array
    r[0] = r0 # initial r-value
    J = np.array(([0,1],[-1,0])) # antisymmetric matrix (meaning its transpose equals its negative)
    for i in range(1, 5*N):
        r[i] = r[i-1] + delta_t*(J@(r[i-1])) # euler's method
    return r # r vector

def heun_method(r0, N):
    delta_t = (2*np.pi)/N # delta t
    r = np.zeros((5*N, 2)) #creates the 5N by 2 array
    r[0] = r0
    J = np.array(([0,1],[-1,0]))
    for i in range(1, 5*N):
        slopeL = (J@(r[i-1])) #take the slope at current point
        r[i] = r[i-1] + delta_t*(J@(r[i-1])) #use Euler's method to approximate next point's location
        slopeR = (J@(r[i])) #take the slope at next point
        r[i] = 1/2(slopeL + slopeR)#then take the average to get the true location
    return r

def gen_plot(x, y, labels=[""], linestyles=[""], xlabel="", ylabel="", title=""):
    # Plotting
    fig = plt.figure(figsize=(8,6)) # Setting funciton figure size (width, height)
    axes = plt.axes() # Creating function plot axes
    for i in range(len(y)): # plot each s
        axes.plot(x, y[i], label=labels[i], linestyle=linestyles[i]) # plotting graph    axes.legend() # axes legend
    axes.legend() # add legend
    plt.xlabel(xlabel) # x-axis label for graph
    plt.ylabel(ylabel) # y-axis label for graph
    plt.title(title) # the title of the graph
    plt.legend() #plot legend

    plt.show() # show the graph