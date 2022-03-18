# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 22:39:46 2021

@author: msiodlac
"""
import scipy as sp


################################################################################
## Heat load of the sensor wires ##

# Assumptions:
# -> 4 wires per sensor
# -> AWG 36, 0.127mm
# -> Manganin
# -> TS @ 50K and CM @ 4K

# # Cross section of one sensor wire
# A_Sensor = sp.pi * 0.25 * 0.000127**2 #m²
# # Sum of the heat load of 4 wires on the specific shield
# Q_Sensor_TS = 4 * 1/1 * A_Sensor * (4465.87-262.69) #W
# Q_Sensor_CM = 4 * 1/1 * A_Sensor * (262.69-0.81) #W
# # Print results: per sensor per meter
# print("Q_Sensor_TS: ", Q_Sensor_TS)
# print("Q_Sensor_CM: ", Q_Sensor_CM)

################################################################################

################################################################################
## Heat load and cooldown time of the calibration system (CS) ##

# Assumptions:
# -> Copper RRR = 100
# -> CM @ 4.5 K, Calibration system dT = 20 K

# Temperature boundary conditions
T_CM = 4.5 #K
T_0 = T_CM + 20 #K
T_end = 4.51 #K

# Length to the CM
L_CS = 0.1 #m
# Mass of calibration system CS
m_CS = 0.03 #kg

# Data extracted from Ekin “Experimental Techniques for Low Temperature Measurements”
# Fitted with Excel
# CAREFUL: Equation only acceptable between 4 and 50 K !
def cp_Brass(T):
    return 0.0025676366 * T**2.7623630373 #J/(kgK)


# ## Coax cable Qaxial RG402-CU -> see datasheet
# # Equivalent cross section of the cooper in the coax cable
# A_coax_Cu = sp.constants.pi * 0.25 * (0.00092**2 + 0.00358**2 - 0.00299**2) #m²
# # Cross section of PTFE in the coax cable
# A_coax_PTFE = sp.constants.pi * 0.25 * (0.00299**2 - 0.00092**2) #m²
# # Heat load of the calibration system during operation
# Q_CS_CM = 1/L_CS * (A_coax_Cu*lambda_Cu_100(T_CM,T_0) + A_coax_PTFE*lambda_PTFE(T_CM,T_0)) #W
# # Integration of the transient heat transfer equation
# result_int = sp.integrate.quad(lambda x: ( m_CS * cp_Brass(x)
#                                          / (1/L_CS * (A_coax_Cu*lambda_Cu_100(T_CM,x) + A_coax_PTFE*lambda_PTFE(T_CM,x))) )
#                                          , T_end, T_0)
# # Transfering the target value
# t_end = result_int[0]



## Coax cable Keycom Semi-Rigid ULT-05
# https://www.keycom.co.jp/eproducts/upj/upj2/page.htm
# Function for the low conductivity caox cable from the data sheet
# Conservative conductivity and then linear dT -> see excel for explanation
def Q_ULT05(T_hot, T_cold, L):
    return 0.000065 * 1/L * (T_hot-T_cold)

# Heat load of the calibration system during operation
Q_CS_CM = Q_ULT05(T_0, T_CM, L_CS) #W
# Integration of the transient heat transfer equation
result_int = sp.integrate.quad(lambda x: ( m_CS * cp_Brass(x)/Q_ULT05(x, T_CM, L_CS) )
                                         , T_end, T_0)
# Transfering the target value
t_end = result_int[0]


# Print results:
print("Calibration System cooldown time: ", t_end/60 , " min")
print("Q_CS_CM: ", Q_CS_CM)

################################################################################
