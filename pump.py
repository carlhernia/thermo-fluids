# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 22:09:58 2022

@author: Carlos Hernandez, William Wilches
"""

import pint
import numpy as np
import matplotlib.pyplot as pl
import CoolProp.CoolProp as cp
import scipy.optimize as sc

##### Unit Conversion w/ Pint ################################################

u = pint.UnitRegistry()
Q = u.Quantity

def convert(val, old_unit, new_unit):
    original = Q(val, old_unit)
    new = original.to(new_unit)
    return new.magnitude

def swameeJain(Re, D, e):
    return 0.25 * (np.log10((e/D)/3.7 + 5.4/(Re**0.9)))**(-2)

def colebrookFunc(f, *args):
    Re, D, e = args
    return 1 / np.sqrt(f) + 2.0 * np.log10((e/D)/3.7 + 2.51/(Re*np.sqrt(f)))
    
def colebrookEq(Re, D, e):
    f0 = swameeJain(Re, D, e)
    return sc.fsolve(colebrookFunc, f0, args=(Re, D, e))

# fluid properties
fluid = 'Water'

# rho = cp.PropsSI('D', 'Q', 0, 'P', P, fluid) # density, kg/m^3
# rho = convert(rho, 'kg/m^3', 'lb/ft^3')
rho = 62.3 # density, kg/m^3

# mu = cp.PropsSI('V', 'Q', 0, 'P', P, fluid)
# mu = convert(mu, 'Pa*s', 'lb/ft/s')
mu = 6.556E-4 # viscosity, lb/ft/s

# diameter selection criteria
e = 0.000005 # pipe surface relative roughness, ft
vmax = 10 # limit erosional velocity, ft/s
Hmax = 3 # limit friciotnal pipe losses per 100ft of pipe, ft

g = 32.2 # gravity, ft/s^2

h_equip = 15 # equipment loss, ft
# H_L = h_major + h_minor + h_equip


# schedule 40 standard pipe sizes
diameters = np.array([0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 3.5,
                      4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 24])

Vdot = convert(600, 'gal/min', 'ft**3/s') # volumetric flow rate

for D in diameters:
    D = convert(D, 'in', 'ft') # pipe diameter, ft
    A = np.pi * D**2 / 4 # pipe area, ft^2
    v = Vdot / A # average flow velocity, ft/s
    Re = rho * v * D / mu # Reynolds number
    
    if Re < 2300:
        f = 64 / Re
    else:
        f = colebrookEq(Re, D, e)
    
    H = f * 100 / D * v**2 / (2 * g) # friction loss per 100ft
    
    if v <= vmax and H <= Hmax:
        break
    
D = convert(D, 'ft', 'in')
print('Economical pipe size:', D, 'inches')


















