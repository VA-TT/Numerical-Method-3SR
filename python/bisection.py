# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 14:19:38 2025

@author: User
"""
import numpy as np
import matplotlib.pyplot as plt

errors = []
steps = []
step = 0
max_iteration = 1000
x1 = -100
x2 = 100
tolerance = 1e-8

def myFunction(x):
    return x**3 - 2

def bisection(a ,b, eps = 1e-8, step=0):
    while (step < max_iteration):
        fa = myFunction(a)
        fb = myFunction(b)
    
        if (fa * fb >0):
            print("No root could be found in this range!")
            return np.nan
    
        error = abs(a - b)/2
        errors.append(error)
        steps.append(step)
        x_average = (a + b) * 0.5
        f_average = myFunction(x_average)
        if error < eps:
            return x_average
        if (fa * fb <0):
            if (fa * f_average <0):
                return bisection(a, x_average, eps, step + 1)
            elif (fb * f_average < 0) :
                return bisection(x_average, b, eps, step + 1)
            

root = bisection(x1, x2, tolerance)
print ("Approximative root is: ", root)


plt.semilogy(steps, errors)
plt.xlabel(r"Iteration $n^{th}$")
plt.ylabel("Error")
plt.title("Convergence of Bisection Method")
plt.grid(True)
plt.show()