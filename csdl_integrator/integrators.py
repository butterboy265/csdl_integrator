import numpy as np
import csdl_alpha as csdl

"""
csdl_integrator.py

A Python library providing numerical integration methods using CSDL.

Functions:
    trapezoidIntegrator - Implements the trapezoidal integration method.
    rk4Integrator - Implements the 4th-order Runge-Kutta integration method.
    csdlIntegrator - General-purpose integrator that supports multiple methods.
"""

# TRAPEZOID INTEGRATOR
def trapezoidIntegrator(func, timeInt, initCond, *args, numSteps=100):
    t_0 = timeInt[0]
    t_F = timeInt[1]  # Assigning initial and final timestamps

    h = (t_F - t_0) / numSteps  # Determine step size
 
    t = csdl.linear_combination(t_0, t_F, numSteps + 1) # Time vector

    # Initialize the solution array (2D array for multiple ODEs)
    y = np.zeros((numSteps + 1, initCond.shape[0]))  # Shape: (numSteps + 1, number of equations)
    y = csdl.Variable(value=y)
     
    y = y.set(csdl.slice[0, :], initCond) # Set initial condition
    
    for i in csdl.frange(numSteps):
        t_i = t[i]
        y_i = y[i]

        dydt_i = func(t_i, y_i, *args)  # Get the derivatives at the current step
        y_i1_guess = y_i + h * dydt_i  # Initial guess for the next step
        
        dydt_i1 = func(t_i + h, y_i1_guess, *args)  # Calculate the derivatives at the guessed next step
        y_i1 = y_i + (h / 2) * (dydt_i + dydt_i1)  # Update using the trapezoidal rule

        y = y.set(csdl.slice[i + 1, :], y_i1)

    return t, y

# RK4 INTEGRATOR
def rk4Integrator(func, timeInt, initCond, *args, numSteps=100):
    t_0 = timeInt[0]
    t_F = timeInt[1]  # Assigning initial and final timestamps

    h = (t_F - t_0) / numSteps  # Determine step size

    t = csdl.linear_combination(t_0, t_F, numSteps + 1)  # Time vector

    # Initialize the solution array (2D array for multiple ODEs)
    y = np.zeros((numSteps + 1, initCond.shape[0]))  # Shape: (numSteps + 1, number of equations)
    y = csdl.Variable(value=y)

    y = y.set(csdl.slice[0, :], initCond) # Set initial condition

    for i in csdl.frange(numSteps):
        t_i = t[i]
        y_i = y[i]

        k_1 = func(t_i, y_i, *args)
        k_2 = func(t_i + h/2, y_i + h * (k_1/2), *args)
        k_3 = func(t_i + h/2, y_i + h * (k_2/2), *args)
        k_4 = func(t_i + h, y_i + h * k_3, *args)

        y_i1 = y_i + (h/6) * (k_1 + 2*k_2 + 2*k_3 + k_4)
        
        y = y.set(csdl.slice[i + 1, :], y_i1)  # Store the result in the solution array

    return t, y


# Main integrator function that takes a method as input
def integrate(func, timeInt, initCond, *args, method='trapezoid', numSteps=100):
    # Dictionary to map method names to functions
    methods = {
        'trapezoid': trapezoidIntegrator,
        'rk4': rk4Integrator
    }

    if method not in methods:
        raise ValueError(f"Unknown method '{method}'. Available methods: {list(methods.keys())}")

    # Call the selected method
    integrator = methods[method]
    return integrator(func, timeInt, initCond, *args, numSteps=numSteps)