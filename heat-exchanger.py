# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 20:22:29 2022

@author: Carlos Hernandez, William Wilches
"""

import pint
import numpy as np
import matplotlib.pyplot as pl
import CoolProp.CoolProp as cp

##### Unit Conversion w/ Pint ################################################

u = pint.UnitRegistry()
Q = u.Quantity

def convert(val, old_unit, new_unit):
    original = Q(val, old_unit)
    new = original.to(new_unit)
    return new.magnitude

def fricfac(Re_D):
    if Re_D > 2300:
        return (1.58 * np.log(Re_D) - 3.28)**(-2)
    else:
        return 16 / Re_D

def tube_nusselt(D_h, L_t, Re_D, Pr):
    if Re_D > 2300:
        return (fricfac(Re_D) / 2 * (Re_D - 1000) * Pr) / (1 + 12.7 * (fricfac(Re_D) / 2)**(1/2) * (Pr**(2/3) - 1))
    else:
        return 1.86 * (D_h * Re_D * Pr / L_t)**(1/3)

def shell_nusselt(Re_D, Pr):
    return 0.36 * Re_D**(0.55) * Pr**(1/3)

L_t = convert(254, 'mm', 'm') # tube length, m
D_s = convert(50.8, 'mm', 'm') # shell inside diamter, mm
N_p = 1 # nubmer of passes

# B = L_t / (N_b + 1) # baffle spacing
B = convert(25.4, 'mm', 'm') # baffle spacing, assumed, mm
N_b = L_t / B - 1 # baffle number, derived

d_r = 1.3 # diameter ratio
d_o = convert(3.175, 'mm', 'm') # tube outside diameter, mm
d_i = d_o / d_r # tube inside diamter

P_r = 1.25 # tube pitch ratio
P_t = P_r * d_o # tube pitch

C_t = P_t - d_o # tube clearance
A_c = D_s * C_t * B / P_t # cross-flow area of shell

layout = 'triangular' # change to square or triangular based on tube layout

if layout == 'square':
    D_e = 4 * (P_t**2 - np.pi * d_o**2 / 4) / (np.pi * d_o) # D_e for square pitch layout
else:
    D_e = 4 * ((3**(0.5) * P_t**2) / 4 - (np.pi * d_o**2) / 8) / (np.pi * d_o / 2) # D_e for triangular pitch layout

# tube count constant
if N_p == 1:
    CTP = 0.93
elif N_p == 2:
    CTP = 0.90
else:
    CTP = 0.85

# tube layout constant
# CL = 1 # square-pitch layout
CL = 0.866 # triangular-pitch layout

A_shaded = CL * P_t**2 # shaded area
N_t = int((CTP) * (np.pi * D_s**2 / 4) / A_shaded)

##############################################################################

# engine oil given values
fluid_h = 'engine oil'
m_dot_h = 0.23 # flow rate, kg/s
T_h_i = 120 # initial oil temp, C
T_h_o = 115 # max final oil temp, C
R_fi = 0.176E-3 # oil fouling factor, m**2*K/W

# engine coolant given values
fluid_c = 'engine coolant (50% ethylene glycol)'
m_dot_c = 0.47 # flow rate, kg/s
T_c_i = 90 # initial coolant temp, C
T_c_o = 100 # final coolant temp, C
R_fo = 0.353E-3 # coolant fouling factor, m**2*K/W

# other given values
k_w = 42.7 # thermal conductivity, W/m/K
P_L_max = 10E3 # max pressure drop, kPa

T_h_av = (T_h_i + T_h_o) / 2 # avg oil temp, C
T_c_av = (T_c_i + T_c_o) / 2 # avg coolant temp, C

##### use coolprop here, maybe #####
# engine oil table properties
rho_1 = 828 # kg/m**3
c_p_1 = 2307 # J/kg/K
k_1 = 0.135 # W/m/K
mu_1 = 1.027E-2 # N*s/m
Pr_1 = 175

# engine coolant table properties
rho_2 = 1020 # kg/m**3
c_p_2 = 3650 # J/kg/K
k_2 = 0.442 # W/m/K
mu_2 = 0.08E-2 # N*s/m
Pr_2 = 6.6

##### tube side (engine oil) #################################################

A_c1 = (np.pi * d_i**2 / 4) * (N_t / N_p) # cross-flow area, m**2

v_1 = m_dot_h / (rho_1 * A_c1) # velocity, m/s
Re_1 = rho_1 * v_1 * d_i / mu_1 # Reynolds number, engine oil

Nu_1 = tube_nusselt(d_i, L_t, Re_1, Pr_1) # oil nusselt number
h_1 = Nu_1 * k_1 / d_i # oil conv heat transfer coefficient, W/m**2/K

##### shell side (50% ethylene glycol) #######################################

A_c2 = D_s * C_t * B / P_t # cross-flow area, m**2

v_2 = m_dot_c / (rho_2 * A_c2) # velocity, m/s
Re_2 = rho_2 * v_2 * D_e / mu_2 # Reynolds number, engine coolant

Nu_2 = shell_nusselt(Re_2, Pr_2) # coolant nusselt number
h_2 = Nu_2 * k_2 / D_e # coolant conv heat transfer coefficient, W/m**2/K

A_i = np.pi * d_i * L_t * N_t # total heat transfer area for oil, m**2
A_o = np.pi * d_o * L_t * N_t # total heat transfer area for coolant, m**2

# overall heat transfer coefficient, W/K
UA_o = 1 / (1 / (h_1 * A_i) + R_fi / A_i + np.log(d_o/d_i) / (2 * np.pi * k_w * L_t * N_t) + R_fo / A_o + 1/ (h_2 * A_o))

##### epsilon-NTU method #####################################################

C_1 = m_dot_h * c_p_1 # heat capacity for oil, W/K
C_2 = m_dot_c * c_p_2 # heat capactiy for coolant, W/K

C_min = min([C_1, C_2]) # min heat capacity, W/K
C_max = max([C_1, C_2]) # max heat capacity, W/K

C_r = C_min / C_max # heat capacity ratio

NTU = UA_o / C_min # number of transfer units
NTU_1 = NTU / N_p

# effectiveness for shell-and-tube hx
eff = 2 * (1 + C_r + (1 + C_r**2)**(0.5) * (1 + np.exp(-NTU_1 * (1 + C_r**2)**(0.5))) / (1 - np.exp(-NTU_1 * (1 + C_r**2)**(0.5))))**(-1)

T_h_o = T_h_i - eff * C_min / C_1 * (T_h_i - T_c_i) # oil outlet temp
T_c_o = T_c_i + eff * C_min / C_2 * (T_h_i - T_c_i) # coolant outlet temp

q = eff * C_min * (T_h_i - T_c_i) # heat transfer rate, W

delP_1 = 4 * (fricfac(Re_1) * L_t / d_i + 1) * N_p * 1/2 * rho_1 * v_1**2 # oil pressure loss, kPa
delP_2 = fricfac(Re_2) * D_s / D_e * (N_b + 1) * 1/2 * rho_2 * v_2**2 # coolant pressure loss, kPa
print('Engine oil pressure loss:', delP_1, 'kPa')
print('Engine coolant pressure loss:', delP_2, 'kPa')
print('')

beta = (A_o + A_i) / (np.pi * D_s**2 / 4 * L_t) # surface density, m**2/m**3
print('Surface density:', beta, 'm^2/m^3')
print('Heat transfer rate', q, 'W')



