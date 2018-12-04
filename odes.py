#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###
# Name: Jacob Anabi, Gage Kizzar
# Student ID: 2294644, 2291700
# Email: anabi@chapman.edu, kizzar@chapman.edu
# Course: PHYS220/MATH220/CPSC220 Fall 2018
# Assignment: CW11
###

"""
odes.py Module Decription:
    This module allows computation and plotting of some vector r(t)=[x(t),v(t)]
    using Euler's method, Heun's method, 2nd Order Runge-Kutta method, and 4th Order Runge-Kutta method.
    The actual definition and computation is also provided.
"""

import numpy as np
from matplotlib import pyplot as plt

def gen_method(r0, N):
    """
    gen_method function description:
        This method computes the vector r(t)'s true value.

        Args:
            r0 - the initial r-value
            N - the number of steps in each period
    """
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
    """
    euler_method function description:
        This method computes the vector r(t)'s using Euler's method.

        Args:
            r0 - the initial r-value
            N - the number of steps in each period
    """
    delta_t = (2*np.pi)/N # delta t
    r = np.zeros((5*N, 2)) # 5Nx2 array
    r[0] = r0 # initial r-value
    J = np.array(([0,1],[-1,0])) # antisymmetric matrix (meaning its transpose equals its negative)
    for i in range(1, 5*N):
        r[i] = r[i-1] + delta_t*(J@(r[i-1])) # euler's method
    return r # r vector

def heun_method(r0, N):
    """
    heun_method function description:
        This method computes the vector r(t)'s using Heun's method.

        Args:
            r0 - the initial r-value
            N - the number of steps in each period
    """
    delta_t = (2*np.pi)/N # delta t
    r = np.zeros((5*N, 2)) # creates the 5N by 2 array
    r[0] = r0
    J = np.array(([0,1],[-1,0]))
    for i in range(1, 5*N):
        r[i] = r[i-1]+delta_t*(J@(r[i-1])) # r_intermediate
        r[i] = r[i-1]+(delta_t/2)*(J@(r[i-1])+J@(r[i])) # actual r
    return r

def runge_kutta_2ndOrd(r0, N):
    """
    runge_kutta_2ndOrd function description:
        This method computes the vector r(t)'s using Runge Kutta 2nd Order's method.

        Args:
            r0 - the initial r-value
            N - the number of steps in each period
    """
    delta_t = (2*np.pi)/N # delta t
    r = np.zeros((5*N, 2)) # creates the 5N by 2 array
    r[0] = r0
    J = np.array(([0,1],[-1,0]))
    for i in range(1, 5*N):
        k1 = delta_t*(J@(r[i-1]))
        k2 = delta_t*(J@(r[i-1] + k1))
        r[i] = r[i-1] + 1/2*(k1+k2)
    return r

def runge_kutta_4thOrd(r0, N):
    """
    runge_kutta_4thOrd function description:
        This method computes the vector r(t)'s using Runge Kutta 4th Order's method.

        Args:
            r0 - the initial r-value
            N - the number of steps in each period
    """
    delta_t = (2*np.pi)/N # delta t
    r = np.zeros((5*N, 2)) # creates the 5N by 2 array
    r[0] = r0
    J = np.array(([0,1],[-1,0]))
    for i in range(1, 5*N):
        k1 = delta_t*(J@(r[i-1]))
        k2 = delta_t*(J@(r[i-1] + k1/2))
        k3 = delta_t*(J@(r[i-1] + k2/2))
        k4 = delta_t*(J@(r[i-1] + k3))
        r[i] = r[i-1] + (k1+2*k2+2*k3+k4)/6
    return r

def gen_plot(x, y, labels=[""], linestyles=[""], xlabel="", ylabel="", title=""):
    """
    gen_plot function description:
        This function plots some generic x and y values

        Args:
            x - the domain
            y - an 2D array of multiple range values based on the domain
            labels - list of legend labels
            linestyles - list of linestyles
            xlabel - label for the domain
            ylabel - label for the range
            title - title of the graph
    """
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