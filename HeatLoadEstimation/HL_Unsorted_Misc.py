"""
Heat load estimation for the RADES-BabyIAXO cryostat.

Multiple unsorted approaches and calculations.



09.06.2021
"""

import math
import numpy as np
import scipy as sp
from scipy import constants
from scipy import optimize
from scipy import integrate
import matplotlib.pyplot as plt

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

################################################################################
## 1.1 Definition of the functions used in the later calculation ##

# Function: calculating the surface of a cylinder depending on the length and the diameter
def surface_cylinder(L, D):
    A = sp.constants.pi * D * L + 2 * sp.constants.pi * D**2 / 4
    return A #m²

#Function: calculating the heat load of conduction according to Fourier's law
def Q_cond(A_surface, d_distance, T_hot, T_cold, lambda_cond):
    return A_surface * lambda_cond * 1/d_distance * (T_hot - T_cold) #W


# Function: calculating the heat load of thermal radiation for a cold surface which is enclosed by a warm surface
# According to "Tieftemperaturtechnologie" by Frey and Ekin for diffuse reflection
def Q_therm_rad(A_hot, A_cold, T_hot, T_cold, epsilon_hot, epsilon_cold):
    # Emission coefficent for an enclosed vessel
    epsilon_effective = 1/ (1/epsilon_cold + A_cold/A_hot * (1/epsilon_hot - 1)) #-
    # Simple equation with an effective emission coefficent
    Q_rad = epsilon_effective * sp.constants.Stefan_Boltzmann * A_cold * (T_hot**4 - T_cold**4) #W
    return Q_rad #W

# Function: calculating the heat load of gas conduction between the hot vessel and the colder MLI
# Assumption: mean free path >> distance between the walls
# According to Grohmann lecture and Frey and Ekin (Ekin already simplified a lot)
def Q_gas_cond(A_cold, T_hot, T_cold, pressure_vac):
    # Equation to calculate acommodation coefficent alpha_0
    # alpha_hot = 0.3
    # alpha_cold = 1
    # alpha_0 = alpha_cold * alpha_hot / (alpha_hot + alpha_cold * (1-alpha_hot) * A_cold/A_hot)
    # print(alpha_0)

    # Good assumption for helium according to Grohmann lecture and Ekin
    alpha_0 = 0.5 #-   ->conservative value
    kappa = 5/3 #-

    # Arithmetic mean for K according to Grohmann lecture
    T_mean = (T_hot + T_cold)/2.0 #K
    # "Constant" factor K
    K = (sp.constants.R / (8 * sp.constants.pi * 4*10**-3 * T_mean)) * (kappa+1)/(kappa-1) #W/(m²PaK)
    # Final equation
    Q_gas = A_cold * alpha_0 * K * pressure_vac * (T_hot - T_cold) #W
    return Q_gas #W

# Function: calculating the heat transport in MLI due to thermal radiation and conduction
# The model's values are based on experience with the LHC cryostat
# According to Parma
def q_MLI_Parma(N_MLI, T_hot, T_cold):
    # Empirical model parameters
    alpha = 1.401*10**-4
    beta = 3.741*10**-9

    return (beta/(N_MLI+1) * (T_hot**4-T_cold**4)) + (alpha/(N_MLI+1) * (T_hot+T_cold)/2 * (T_hot-T_cold)) #W/m²

# Function: calculating the heat transport in MLI due to thermal radiation and conduction
# The model's values are based on experience with the LHC cryostat
# According to Riddone
def q_MLI_Riddone(N_MLI, T_hot, T_cold):
    # Empirical model parameters
    alpha = 1.401*10**-4
    beta = 3.741*10**-9

    return (beta/(N_MLI) * (T_hot**4-T_cold**4)) + (alpha/(N_MLI) * (T_hot+T_cold)/2 * (T_hot-T_cold)) #W/m²


# Integral thermal conductivity of PTFE according to NIST
# Integral form, just put in the upper and the lower boundaries
# Data and equation range: 4 - 300 K
# Curve fit % error relative to data: 5
def lambda_PTFE(T_cold, T_hot):
    result_int = sp.integrate.quad(lambda x: 10**(2.7380 + -30.677*math.log(x,10) + 89.430*math.log(x,10)**2 + -136.99*math.log(x,10)**3
                                              + 124.69*math.log(x,10)**4 + -69.556*math.log(x,10)**5 + 23.320*math.log(x,10)**6
                                              + -4.3135*math.log(x,10)**7 + 0.33829*math.log(x,10)**8), T_cold, T_hot)
    return result_int[0] #W/m

# Integral thermal conductivity of Cu RRR = 50 according to NIST
# Integral form, just put in the upper and the lower boundaries
# Data and equation range: 4 - 300 K
# Curve fit % error relative to data: 2
def lambda_Cu_50(T_cold, T_hot):
    result_int = sp.integrate.quad(lambda x: 10**((1.8743 + -0.6018*x**0.5 + 0.26426*x + -0.051276*x**1.5 + 0.003723*x**2)
                                                 /(1 + -0.41538*x**0.5 + 0.13294*x + -0.0219*x**1.5 + 0.0014871*x**2)), T_cold, T_hot)
    return result_int[0] #W/m

# Integral thermal conductivity of Cu OFHC RRR = 100 according to NIST
# Data and equation range: 4 - 300 K
# Curve fit % error relative to data: 1
def lambda_Cu_100(T_cold, T_hot):
    result_int = sp.integrate.quad(lambda x: 10**((2.2154 + -0.88068*x**0.5 + 0.29505*x + -0.04831*x**1.5 + 0.003207*x**2)
                                                 /(1 + -0.47461*x**0.5 + 0.13871*x + -0.02043*x**1.5 + 0.001281*x**2)), T_cold, T_hot)
    return result_int[0] #W/m

# Integral thermal conductivity of Cu RRR = 300 according to NIST
# Integral form, just put in the upper and the lower boundaries
# Data and equation range: 4 - 300 K
# Curve fit % error relative to data: 1
def lambda_Cu_300(T_cold, T_hot):
    result_int = sp.integrate.quad(lambda x: 10**((1.357 + 2.669*x**0.5 + -0.6683*x + 0.05773*x**1.5 + 0*x**2)
                                                 /(1 + 0.3981*x**0.5 + -0.1346*x + 0.01342*x**1.5 + 0.0002147*x**2)), T_cold, T_hot)
    return result_int[0] #W/m

# Specific Heat capacity  of Cu according to NIST
# Data and equation range: 4 - 300 K
# Curve fit % error relative to data: 10 (T < 15K); 5 (T ≥ 15K)
def cp_Cu(T):
    return 10**(-1.91844 + -0.15973*math.log(T,10) + 8.61013*math.log(T,10)**2 + -18.996*math.log(T,10)**3
                + 21.9661*math.log(T,10)**4 + -12.7328*math.log(T,10)**5 + 3.54322*math.log(T,10)**6
                + -0.3797*math.log(T,10)**7 + 0*math.log(T,10)**8) #J/kgK



# Thermal conductivity of SST304 according to NIST
# Data and equation range: 1 - 300 K
# Curve fit % error relative to data: 2
def lambda_SST304(T):
    return 10**(-1.4087 + 1.3982*math.log(T,10) + 0.2543*math.log(T,10)**2 + -0.6260*math.log(T,10)**3
                + 0.2334*math.log(T,10)**4 + 0.4256*math.log(T,10)**5 + -0.4658*math.log(T,10)**6
                + 0.1650*math.log(T,10)**7 + -0.0199*math.log(T,10)**8) #W/(mK)

# Integral thermal conductivity of SST304 according to NIST
# Integral form, just put in the lower and upper boundaries
# Data and equation range: 1 - 300 K
# Curve fit % error relative to data: 2
def lambda_int_SST304(T_cold, T_hot):
    result_int = sp.integrate.quad(lambda x: 10**(-1.4087 + 1.3982*math.log(x,10) + 0.2543*math.log(x,10)**2 + -0.6260*math.log(x,10)**3
                                                  + 0.2334*math.log(x,10)**4 + 0.4256*math.log(x,10)**5 + -0.4658*math.log(x,10)**6
                                                  + 0.1650*math.log(x,10)**7 + -0.0199*math.log(x,10)**8), T_cold, T_hot)
    return result_int[0] #W/m

################################################################################
## 1.2 Geometry data and boundary conditions ##

# Length and diameter of the cold mass, thermal shield and vacuum chamber, respectively
L_CM = 9.90 #m
D_CM = 0.60 #m
L_TS = 9.95 #m
D_TS = 0.65 #m
L_VC = 10.0 #m
D_VC = 0.70 #m

# Temperature of the respective vessel/layer
T_CM = 4.0 #K
T_TS = 35.0 #K
T_VC = 293.0 #K

# Pressure in the vacuum system
pressure_vac = 10**-4 #Pa

# Emission coefficent of the respective vessel/layer
epsilon_TS = 0.03 #-
epsilon_VC = 0.08 #-
epsilon_MLI = 0.023 #-


# Surface of the ideal cylinders of the cold mass, thermal shield and vacuum chamber, respectively
A_CM = surface_cylinder(L_CM, D_CM) #m²
A_TS = surface_cylinder(L_TS, D_TS) #m²
A_VC = surface_cylinder(L_VC, D_VC) #m²


################################################################################
## 2.1.1 Calculation of the heat load with effective heat conduction coeff. and equilibrium on the first MLI layer ##
## Variation of T_TS and T_CM and calculation of the respective heat load ##

# Effective heat conduction coefficent of the respective MLI
lambda_eff_MLI_30 = 4*10**-5 #W/(mK)
lambda_eff_MLI_20 = 4*10**-5 #W/(mK)

# MLI thickness (value from Frey, Tieftemperaturtechnologie for 25 layers)
d_MLI_eff = 0.0125 #m

# print("q_MLI_Tieftemp: ", lambda_eff_MLI_30/d_MLI_eff)
"""
# Definition of an function to be able to use the newton method to find the MLI temperaure
def func_newton(x, A_hot, A_cold, T_hot, epsilon_hot, epsilon_MLI, pressure_vac, lambda_eff, d_MLI_eff, T_cold):
    Q_rad = Q_therm_rad(A_hot, A_cold, T_hot, x, epsilon_hot, epsilon_MLI)
    Q_gas = Q_gas_cond(A_cold, T_hot, x, pressure_vac)
    Q_cond_MLI = Q_cond(A_cold, d_MLI_eff, x, T_cold, lambda_eff)
    return Q_rad + Q_gas - Q_cond_MLI

# Create lists for the different temperatures to be calculated
T_TS_lst = [35.0, 45.0, 50.0, 55.0, 65.0, 75.0, 80.0] #K
T_CM_lst = [i for i in np.arange(4,15.5,0.5)] #K

# Allocate a matrix for the heat load depending on the anchoring temperatures
Q_TS_matrix = np.zeros([len(T_TS_lst), len(T_CM_lst)], dtype=float)
Q_CM_matrix = np.zeros([len(T_TS_lst), len(T_CM_lst)], dtype=float)

# Two for loops to calculate every point of the matrix
for i in range(0,len(T_TS_lst)):
    for j in range(0,len(T_CM_lst)):
        # Find fitting T_MLI_30 temperature for the respective T_TS and T_CM pair using the newton method
        T_MLI_30 = sp.optimize.newton(func_newton, T_VC, args=(A_VC, A_TS, T_VC, epsilon_VC, epsilon_MLI, pressure_vac, lambda_eff_MLI_30, d_MLI_eff, T_TS_lst[i]),
                                      rtol=10**-4, disp=True)
        # Find fitting T_MLI_30 temperature for the respective T_TS and T_CM pair using the newton method
        T_MLI_20 = sp.optimize.newton(func_newton, T_TS, args=(A_TS, A_CM, T_TS_lst[i], epsilon_TS, epsilon_MLI, pressure_vac, lambda_eff_MLI_20, d_MLI_eff, T_CM_lst[j]),
                                      rtol=10**-4, disp=True)

        # Save the sum of both heat loads on a matrix for postprocessing
        Q_TS_matrix[i][j] = Q_cond(A_TS, d_MLI_eff, T_MLI_30, T_TS_lst[i], lambda_eff_MLI_30)
        Q_CM_matrix[i][j] = Q_cond(A_CM, d_MLI_eff, T_MLI_20, T_CM_lst[j], lambda_eff_MLI_20)
"""
################################################################################

################################################################################
## 2.1.2 Calculation of the heat load with effective heat conduction coeff. ##
## Variation of T_TS and T_CM and calculation of the respective heat load ##
"""
# Create lists for the different temperatures to be calculated
T_TS_lst = [35.0, 45.0, 50.0, 55.0, 65.0, 75.0, 80.0] #K
T_CM_lst = [i for i in np.arange(4,15.5,0.5)] #K

# Allocate a matrix for the heat load depending on the anchoring temperatures
Q_TS_matrix = np.zeros([len(T_TS_lst), len(T_CM_lst)], dtype=float)
Q_CM_matrix = np.zeros([len(T_TS_lst), len(T_CM_lst)], dtype=float)

# Two for loops to calculate every point of the matrix
for i in range(0,len(T_TS_lst)):
    for j in range(0,len(T_CM_lst)):
        # Save the sum of both heat loads on a matrix for postprocessing
        Q_TS_matrix[i][j] = Q_cond(A_TS, d_MLI_eff, T_VC, T_TS_lst[i], lambda_eff_MLI_30)
        Q_CM_matrix[i][j] = Q_cond(A_CM, d_MLI_eff, T_TS_lst[i], T_CM_lst[j], lambda_eff_MLI_20)

print(Q_TS_matrix)
print(Q_CM_matrix)
"""
################################################################################


################################################################################
## 2.2 Calculation of the heat load for an MLI system with the model from Parma ##
## Variation of T_TS and T_CM and calculation of the respective heat load ##
"""
# Create lists for the different temperatures to be calculated
# T_TS_lst = [35.0, 45.0, 50.0, 55.0, 65.0, 75.0, 80.0] #K
# T_CM_lst = [i for i in np.arange(4,15.5,0.5)] #K

T_TS_lst = [35.0, 45.0, 50.0, 60.0, 80.0] #K
T_CM_lst = [4.0, 5.0, 10.0, 15.0] #K

# Allocate a matrix for the heat load depending on the anchoring temperatures
Q_TS_matrix = np.zeros([len(T_TS_lst), len(T_CM_lst)], dtype=float)
Q_CM_matrix = np.zeros([len(T_TS_lst), len(T_CM_lst)], dtype=float)

# Two for loops to calculate every point of the matrix
for i in range(0,len(T_TS_lst)):
    for j in range(0,len(T_CM_lst)):
        # Save the sum of both heat loads on a matrix for postprocessing
        Q_TS_matrix[i][j] = A_TS * q_MLI_Parma(30, T_VC, T_TS_lst[i]) + Q_gas_cond(A_TS, T_VC, T_TS_lst[i], pressure_vac)
        Q_CM_matrix[i][j] = A_CM * q_MLI_Parma(10, T_TS_lst[i], T_CM_lst[j]) + Q_gas_cond(A_CM, T_TS_lst[i], T_CM_lst[j], pressure_vac)

print(Q_TS_matrix)
print(Q_CM_matrix)
"""
################################################################################


################################################################################
## 2.3 Calculation of the heat load with the model from Parma and no MLI layer on the cold mass ##
## Variation of T_TS and T_CM and calculation of the respective heat load ##
"""
# Emission coefficent for the cold mass assuming a polished SST surface
epsilon_CM = 0.07

# Create lists for the different temperatures to be calculated
# T_TS_lst = [35.0, 45.0, 50.0, 55.0, 65.0, 75.0, 80.0] #K
# T_CM_lst = [i for i in np.arange(4,15.5,0.5)] #K
T_TS_lst = [35.0, 45.0, 50.0, 80.0] #K
T_CM_lst = [4.0, 5.0, 10.0, 15.0] #K

# Allocate a matrix for the heat load depending on the anchoring temperatures
Q_TS_matrix = np.zeros([len(T_TS_lst), len(T_CM_lst)], dtype=float)
Q_CM_matrix = np.zeros([len(T_TS_lst), len(T_CM_lst)], dtype=float)

# Two for loops to calculate every point of the matrix
for i in range(0,len(T_TS_lst)):
    for j in range(0,len(T_CM_lst)):
        # Save the sum of both heat loads on a matrix for postprocessing
        Q_TS_matrix[i][j] = A_TS * q_MLI_Parma(30, T_VC, T_TS_lst[i])
        # Q_CM_matrix[i][j] = A_CM * q_MLI_Parma(0, T_TS_lst[i], T_CM_lst[j])
        Q_CM_matrix[i][j] = (Q_therm_rad(A_TS, A_CM, T_TS_lst[i], T_CM_lst[j], epsilon_TS, epsilon_CM)
                           + Q_gas_cond(A_CM, T_TS_lst[i], T_CM_lst[j], 10**-2))

# print(Q_TS_matrix)
# print(Q_CM_matrix)
"""
################################################################################


################################################################################
## 2.3.2 Calculation of the heat load for an MLI system with the model from Riddone ##
## Variation of T_TS and T_CM and calculation of the respective heat load ##
"""
# Create lists for the different temperatures to be calculated
# T_TS_lst = [35.0, 45.0, 50.0, 55.0, 65.0, 75.0, 80.0] #K
# T_CM_lst = [i for i in np.arange(4,15.5,0.5)] #K

T_TS_lst = [35.0, 45.0, 50.0, 60.0, 80.0] #K
T_CM_lst = [4.0, 5.0, 10.0, 15.0] #K

# Allocate a matrix for the heat load depending on the anchoring temperatures
Q_TS_matrix = np.zeros([len(T_TS_lst), len(T_CM_lst)], dtype=float)
Q_CM_matrix = np.zeros([len(T_TS_lst), len(T_CM_lst)], dtype=float)

# Two for loops to calculate every point of the matrix
for i in range(0,len(T_TS_lst)):
    for j in range(0,len(T_CM_lst)):
        # Save the sum of both heat loads on a matrix for postprocessing
        Q_TS_matrix[i][j] = A_TS * q_MLI_Riddone(30, T_VC, T_TS_lst[i]) + Q_gas_cond(A_TS, T_VC, T_TS_lst[i], pressure_vac)
        Q_CM_matrix[i][j] = A_CM * q_MLI_Riddone(10, T_TS_lst[i], T_CM_lst[j]) + Q_gas_cond(A_CM, T_TS_lst[i], T_CM_lst[j], pressure_vac)

# print(Q_TS_matrix)
# print(Q_CM_matrix)
"""
################################################################################






################################################################################
## 2.7 Thickness of the thermal shield to yield an acceptable dT of the shield ##
# Assuming radial uniformity to get a 1D problem
"""
# Assuming a thermal shield at 50K = T_0
T_0 = 50.0
# Specific heat load on the thermal shield covered with 30 layers of MLI
q_TS = q_MLI_Parma(30, T_VC, T_0)
# Thermal conductivity coefficent @ 50K for Cu-RRR=30
lambda_Cu_30_50K = 732   #W/(mK)
# Maximal prescribed dT of the thermal shield
dT_max = 10.0 #K

# Minimal thickness as a function of the maximum dT
thickness_TS = (q_TS * L_TS*L_TS) / (2 * lambda_Cu_30_50K * dT_max)
# print("Thickness TS: ", thickness_TS)

# Function for calculation of the local temperature with Neumann BC
# Note: x is the dimensionless length
def temp_distribution_shield_Neumann(x, L, thickness, T_0, q, heat_cond_coeff):
    # return T_0 - (q * L*L * x*x)/(2 * heat_cond_coeff * thickness) + (q * L*L * x)/(heat_cond_coeff * thickness)
    return T_0 - (q * L*L)/(heat_cond_coeff * thickness) * (0.5*x*x - x)

def temp_distribution_shield_Dirichlet(x, L, thickness, T_0, q, heat_cond_coeff):
    return T_0 - (q * L*L)/(heat_cond_coeff * thickness) * (0.5*x*x - 0.5*x)

# print(temp_distribution_shield(0.5, L_TS, thickness_TS, T_0, q_TS, lambda_Cu_30_50K))


# # Creating a plot of the temperature distribution along the length of the
# x_lst = [i for i in np.arange(0,1.05,0.05)]
# T_TS_local_lst = [0] * len(x_lst)
# for i in range(0,len(x_lst)):
#     T_TS_local_lst[i] = temp_distribution_shield_Neumann(x_lst[i], L_TS, thickness_TS, T_0, q_TS, lambda_Cu_30_50K)
# plt.plot(x_lst, T_TS_local_lst)
# plt.xlabel("Dimensionless length x / -")
# plt.ylabel("Temperature TS / K")
# plt.show()
"""
################################################################################

################################################################################
## 2.8 Thermal shield service box temperautre gradient ##
# Assuming radial uniformity to get a 1D problem
"""
# Geometry shield: 60 cm diameter and 1 m length
# Thickness shield = 1.5 mm
# Material: Al -> 711 W/mK @ 35 K
# Specific Heat load: 5.3 W per 2.45 m² of the SB TS
thickness = 1.5e-3 #m
cond_coeff = 711 #W/mK
specific_heat_SB = 5.3/2.45 #W/m²
length = 1.0 #m
temperature_TS_SB = 33 #K

# Function for calculation of the local temperature with Neumann BC
# Note: x is the dimensionless length
# x = 0 is the beginning of the shield, x = 1 is the end of the shield
def temp_distribution_shield_Neumann(x, L, thickness, T_0, q, heat_cond_coeff):
    # return T_0 - (q * L*L * x*x)/(2 * heat_cond_coeff * thickness) + (q * L*L * x)/(heat_cond_coeff * thickness)
    return T_0 - (q * L*L)/(heat_cond_coeff * thickness) * (0.5*x*x - x)

def temp_distribution_shield_Dirichlet(x, L, thickness, T_0, q, heat_cond_coeff):
    return T_0 - (q * L*L)/(heat_cond_coeff * thickness) * (0.5*x*x - 0.5*x)

# print(temp_distribution_shield_Neumann(0.5, length, thickness, temperature_TS_SB, specific_heat_SB, cond_coeff))
"""
################################################################################





################################################################################
## 3. Mechanical stability of the vessels ##

## 3.1 Using Parma's rule of thumb
# t/r >3.7%
# t=13mm @ 35cm radius
#
# print(0.037*0.3)

# # Thickness of the different vessels using rule of thumb from Parma
# t_VC = 13 mm
# t_TS = 12 mm
# t_CM = 11.1 mm

# # Weight of the vessels - ballpark
# density_SST304 = 8030 #kg/m³
# density_Cu = 8900 #kg/m³
# print("weight_TS: ", 0.002 * density_Cu * A_TS)
# print("weight_CM: ", 0.005 * density_SST304 * A_CM)

#
# weight_VC = density_SST304 * A_VC * 0.037*0.35 #kg
# weight_TS = density_Cu * A_TS * 0.037*0.325 #kg
# weight_CM = density_SST304 * A_CM * 0.037*0.3 + density_SST304 * A_CM * 0.037*0.3 #kg
#
# print("weight_VC: ", weight_VC)
# print("weight_TS: ", weight_TS)
# print("weight_CM: ", weight_CM)




# ## 3.2 Using Tresca failure criterion ##
# # sigma_admissible > p*r/t
# # sigma_admissible = yield strength / safety factor
# # Values for the yield strength according to Ekin
# yield_SST304 = 240 *10**6
# yield_Cu = 70 *10**6
# # Chosing an adequate safety factor
# safety_factor = 2
# # Calculating sigma admissible
# sigma_a_SST304 = yield_SST304 / safety_factor
# sigma_a_Cu = yield_Cu / safety_factor
# # Print the result
# print("Thickness Tresca failure criterion: ", 2* 1.013*10**5 * 0.30 / sigma_a_SST304)


# ## 3.3 Using bending structural analysis ##
# t_TS = 0.008 #m
# r_TS = 0.325 #m
# l_TS = 9.9 #m
# # Youngs modulus of Copper
# YM_Cu = 117*10**9
#
# # Calculation of bending stress
# w_TS = density_Cu * A_TS * t_TS * 9.81 #N
# moment_TS = (w_TS * l_TS) / 8
# I_TS = sp.pi*0.25*(r_TS**4 - (r_TS-t_TS)**4)
# sigma_B_TS = (moment_TS*r_TS)/I_TS
# print("sigma_B_TS: ", sigma_B_TS*10**-6)
#
# #Calculation of the bending height assuming ideal phytagoras law
# bending_TS = l_TS/2 * ((1+sigma_B_TS/YM_Cu)**2 -1)**0.5
# print("bending_TS: ", bending_TS)

################################################################################



## 4. Mechanical analysis of el. equipment insert
"""
radius_insert = 0.2 #m
length_insert = 0.3 #m
press_max = 1.013*10**5 #Pa

## Tresca failure criterion
# sigma_admissible > p*r/t
# sigma_admissible = yield strength / safety factor
# Values for the yield strength according to https://www.curbellplastics.com/Research-Solutions/Plastic-Material-Properties/G10-FR-4-Glass-Epoxy-Properties
yield_G10 = 260 *10**6 #Pa
# Chosing an adequate safety factor
safety_factor = 2
# Calculating sigma admissible
sigma_a_G10 = yield_G10 / safety_factor
# Print the result
print("Thickness Tresca failure criterion: ", 1.013*10**5 * radius_insert / sigma_a_G10)

## Buckling analysis
# Infinite thin tube
# Values for the Young's modulus according to https://www.curbellplastics.com/Research-Solutions/Plastic-Material-Properties/G10-FR-4-Glass-Epoxy-Properties
# E_G10 = 16547417504 #Pa
# Values for the Young's modulus according to G:\Departments\TE\Groups\CRG\Sections\CI\Projects_and_Services\Cryolab\He II pressurized G10 insert\1 - Project documents\G10 data sheets (Von Roll)
E_G10 = 20e9 #Pa

# Poission's ratio: just assumption
poisson_G10 = 0.25

# thickness of insert
thickness_insert = 0.002 #m
print("Critical pressure from buckling analysis: ", 10**-5 * E_G10/(4*(1-poisson_G10**2)) * (thickness_insert/radius_insert)**3)



A_insert = sp.pi*0.25* ((2*radius_insert)**2-(2*radius_insert-2*thickness_insert)**2) # m²
print("Heat load insert SST: ", A_insert*1/1*lambda_int_SST304(50,300))
print("Heat load insert G10 TS: ", A_insert*1/(1*length_insert)*(133.2-8.48))
print("Heat load insert G10 CM: ", A_insert*1/(0.4)*(8.48-0.0901))

#Additional radiation:
A_thermPlate = sp.pi*0.25* 0.4**2 # m²
print("Radiation on therm plate: ", 67.5e-3*A_thermPlate)
"""









################################################################################
## 5. Postprocessing and visualization ##


# # Restructuring because the values need to be in a 2D array (not lists) and in equal size for 3D plotting
# T_TS_plot_array = np.transpose(np.array([np.array(T_TS_lst)])) * np.full((1,len(T_CM_lst)), 1)
# T_CM_plot_array = np.full((len(T_TS_lst),1), 1) * np.array(T_CM_lst)


# # Plotting the heat loads at the respective temperaure pairs in a 3D plot
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# ax.plot_surface(T_TS_plot_array, T_CM_plot_array, Q_TS_matrix)
#                 #cmap=cm.coolwarm, linewidth=0, antialiased=False)
# ax.plot_surface(T_TS_plot_array, T_CM_plot_array, Q_CM_matrix)

# ax.set_xlabel('T_TS / K')
# ax.set_ylabel('T_CM / K')
# ax.set_zlabel('Q_TS/Q_CM / W')
# # ax.set_zlabel('Q_CM / W')

# # plt.show()



# ## Simple plot to show Parma model ##
# T_warm = [i for i in np.arange(35,81,1.0)]
# HL = [0] * len(T_warm)
# for i in range(0,len(T_warm)):
#     HL[i] = A_CM * q_MLI_Parma(20, T_warm[i], 4.5)
# plt.plot(T_warm, HL)
# plt.xlabel("Temperature TS / K")
# plt.ylabel("Heat load CM @4.5K / W")
# plt.show()











"""
################################################################################
## Parma Slides Exercise ##

A_VV = sp.pi
A_TS = sp.pi * 0.8
A_CM = sp.pi * 0.6

HL_CM = Q_therm_rad(2.51, 1.88, 80.0, 2.0, 0.1, 0.06)
# HL_CM = Q_therm_rad(3.14, 2.51, 293, 80.0, 0.2, 0.1)
# Q_therm_rad(A_hot, A_cold, T_hot, T_cold, epsilon_hot, epsilon_cold)

q_MLI = q_MLI_Parma(30, 290, 80.0) * 2.51
# q_MLI_Parma(N_MLI, T_hot, T_cold)

print("HL_CM: ", HL_CM)
print("q_MLI: ", q_MLI)

# a) Bare cold mass:                                HL_CM = 73 W
# b) Cold mass wrapped with 1 layer of Al foil:     HL_CM = 41 W
# c) Cold mass wrapped with 30 layers of MLI:       HL_CM = 1.96 W
#       calculated the heat load with the formula and 290K and not 293K because there is vacuum in between
# d) Addition of thermal shield actively cooled:    HL_CM = 0.29 W,     HL_TS = 79 W
# e) Wrapping of MLI around thermal shield:         HL_CM = 0.29 W,     HL_TS = 2.38 W
# f) Adding 1 Al foil around cold mass:             HL_CM = 0.187 W,    HL_TS = 2.38 W

################################################################################
"""


"""
################################################################################
## Heat Load Cooper Shield ##
# Using the equations from Boyan's master thesis

# Heat conduction coeff lambda from the literature for a vacuum pressure under 10^-4 mbar and 20 layers of MLI
lambda_apparent = 2.4 * 10**-3; #W/mK

# Surface of the cryostats assuming cylindrical shape
# Main cryostat assuming an value from Boyan's master thesis
A_mc = 1.89 + 0.56; #m^2
# Remote cryostat assuming the cooper shield sitting in the middle between the bore and the cavity - d = 65cm
A_rc = 10 * 3.14159 * 0.65 + 3.14159 * 0.5 * 0.65**2; #m^2

# Temperatures
T_warm = 295; #K
T_shield = 50; #K

# Heat loads Q
Q_mc = lambda_apparent * A_mc * (T_warm - T_shield); #W
Q_rc = lambda_apparent * A_rc * (T_warm - T_shield); #W

print("Heat load main cryostat: ", Q_mc, " W")
print("Heat load remote cryostat: ", Q_rc, " W")
################################################################################
"""
