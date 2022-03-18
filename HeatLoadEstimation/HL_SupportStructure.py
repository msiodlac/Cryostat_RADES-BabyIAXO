# -*- coding: utf-8 -*-
"""
Function to estimate the heat load through the support structure for the BabyIAXO RADES cryostat.


Created on Mon Dec  6 20:43:36 2021
@author: msiodlac
"""
import math
import scipy as sp

# Function: calculating the surface of a cylinder depending on the length and the diameter
def surface_cylinder(L, D):
    A = sp.constants.pi * D * L + 2 * sp.constants.pi * D**2 / 4
    return A #m²

# Integral thermal conductivity of G10 normal direction according to NIST
# https://trc.nist.gov/cryogenics/materials/G-10%20CR%20Fiberglass%20Epoxy/G10CRFiberglassEpoxy_rev.htm
# Data / equation range: 4 - 300 K / 10 - 300 K
# Curve fit % error relative to data: 5
def lambda_G10(T_cold, T_hot):
    result_int = sp.integrate.quad(lambda x: 10**(-4.1236 + 13.788*math.log(x,10) + -26.068*math.log(x,10)**2 + 26.272*math.log(x,10)**3
                                              + -14.663*math.log(x,10)**4 + 4.4954*math.log(x,10)**5 + -0.6905*math.log(x,10)**6
                                              + 0.0397*math.log(x,10)**7 + 0*math.log(x,10)**8), T_cold, T_hot)
    return result_int[0] #W/m


################################################################################
## Heat load of the support structure ##

## Assumptions:
    # Thickness of the TS and CM
    # TS is made out of Cu
    # CM is primarily made out of SST304 with Cu coating, but it could also be made completely out of Cu
        # -> be conservative and use Cu density!
    # CM @ 4 K, TS @ 50 K
    # Thermal conductivity integrals from: Duthil, Material Properties at Low temperature, CAS 2013
        # = Cryolab seminar 08
    # Length between thermalization
    # Perfect thermalization of the support structure with the TS

# Length between thermalization
length_pillar = 0.05 #m

## Weight of the TS and CM
density_Cu = 8900 #kg/m³
density_SST304 = 8030 #kg/m³

thickness_CM = 5e-3 #m
thickness_TS = 2e-3 #m

# Assumption of ideal cylinder shape of CM and TS
A_CM = surface_cylinder(10., 0.6) #m²
A_TS = surface_cylinder(10., 0.65) #m²

weight_CM = density_Cu * thickness_CM * A_CM #kg
weight_TS = density_Cu * thickness_TS * A_TS #kg


# ## SST304 support structure
# # Equivalent cross section of one SST304 pillar
# # Load due to the weight of the vessels
# load_CM = weight_CM*9.81 #N
# load_TS = weight_TS*9.81 #N
# # Cross section using the load, the yield strength and a safety factor of 2
# # Yield strength of SST304 according to Ekin “Experimental Techniques for Low Temperature Measurements”
# yield_strength = 240e+6 #Pa
# A_SS_CM = load_CM * 2 / yield_strength #m²
# A_SS_TS = (load_CM+load_TS) * 2 / yield_strength #m²
# # Print the equivalent diameter of the support structure
# print("d_SS_TS: " ,(4*A_SS_TS/sp.pi)**0.5)
# print("d_SS_CM: " ,(4*A_SS_CM/sp.pi)**0.5)
# # Heat load on the different stages
# Q_SS_TS = A_SS_TS * 1/length_pillar * (3077 - 136) #W
# Q_SS_CM = A_SS_CM * 1/length_pillar * (136 - 0.4) #W


## G10 support structure
# Equivalent cross section of one G10 pillar
# Load due to the weight of the vessels
load_CM = weight_CM*9.81 #N
load_TS = weight_TS*9.81 #N
# Cross section using the load, the yield strength and a safety factor of 2
# Tensile strength of G10 according to https://www.curbellplastics.com/Research-Solutions/Plastic-Material-Properties/G10-FR-4-Glass-Epoxy-Properties
# There is no yield strength, because G10 has no plastic deformation, so I just use this
yield_strength = 262e+6 #Pa
A_SS_CM = load_CM * 3 / yield_strength #m²
A_SS_TS = (load_CM+load_TS) * 3 / yield_strength #m²
# Print the equivalent diameter of the support structure
print("d_SS_TS: " ,(4*A_SS_TS/sp.pi)**0.5)
print("d_SS_CM: " ,(4*A_SS_CM/sp.pi)**0.5)

# Heat load on the different stages
# Using thermal cond. integrals from Cryolab seminar 08
# Q_SS_TS = A_SS_TS * 1/length_pillar * (133.2 - 8.48) #W
# Q_SS_CM = A_SS_CM * 1/length_pillar * (8.48 - 0.0901) #W
# Using thermal cond. from NIST
Q_SS_TS = A_SS_TS * 1/length_pillar * lambda_G10(50,300) #W
Q_SS_CM = A_SS_CM * 1/length_pillar * lambda_G10(4.5,50) #W



## Print the results
print("Q_SS_TS: ", Q_SS_TS)
print("Q_SS_CM: ", Q_SS_CM)

# Three of these rods for stability
# Each of these structures is supposed to be able to carry all the weight
# I am aware that this is not clear at all but it is only supposed to be an estimation for the HL
# A more accurate description of the support structure is not possible now
print("Q_SS_TS Total: ", 3*Q_SS_TS)
print("Q_SS_CM Total: ", 3*Q_SS_CM)


## Joanna Liberadzka zirconium oxide ball transfer units
# See report from Torsten Koettig
# Each unit can support 500 N
# print("Q_SS_Ball_TS: ", load_TS/500 * 1.152*10**-6 * 500**0.1661 * (300**2.195 - 50**2.195))
# print("Q_SS_Ball_CM: ", load_CM/500 * 6.837*10**-6 * 500**0.418 * (50**1.594 - 4**1.594))

# -> They are not the limiting factor and have a higher effective thermal conductivity !
# -> They can be ignored in the heat load estimation

################################################################################
