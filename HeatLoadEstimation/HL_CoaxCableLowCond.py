"""
Heat load estimation for a coax cable.

Numerical approach to simulate the complicated behavior between conductor and dielectric.


06.07.2021
"""

import math
import numpy as np
import scipy as sp
# from scipy import constants
# from scipy import optimize
from scipy import integrate
import matplotlib.pyplot as plt

# from matplotlib import cm

import csv

################################################################################
## Definition of help functions ##

# Thermal conductivity of Cu OFHC RRR=100 according to NIST
# Data and equation range: 4 - 300 K
# Curve fit % error relative to data: 1
def lambda_Cu_100(T):
    return 10**((2.2154 + -0.88068*T**0.5 + 0.29505*T + -0.04831*T**1.5 + 0.003207*T**2)
               /(1 + -0.47461*T**0.5 + 0.13871*T + -0.02043*T**1.5 + 0.001281*T**2)) #W/(mK)

# Integral thermal conductivity of Cu OFHC RRR = 100 according to NIST
# Data and equation range: 4 - 300 K
# Curve fit % error relative to data: 1
def lambda_int_Cu_100(T_cold, T_hot):
    result_int = sp.integrate.quad(lambda x: 10**((2.2154 + -0.88068*x**0.5 + 0.29505*x + -0.04831*x**1.5 + 0.003207*x**2)
                                                 /(1 + -0.47461*x**0.5 + 0.13871*x + -0.02043*x**1.5 + 0.001281*x**2)), T_cold, T_hot)
    return result_int[0] #W/m

# Thermal conductivity of PTFE according to NIST
# Data and equation range: 4 - 300 K
# Curve fit % error relative to data: 5
def lambda_PTFE(T):
    return 10**(2.7380 + -30.677*math.log(T,10) + 89.430*math.log(T,10)**2 + -136.99*math.log(T,10)**3
                + 124.69*math.log(T,10)**4 + -69.556*math.log(T,10)**5 + 23.320*math.log(T,10)**6
                + -4.3135*math.log(T,10)**7 + 0.33829*math.log(T,10)**8) #W/(mK)

# Thermal conductivity of PTFE according to lakeshore
# https://www.lakeshore.com/docs/default-source/product-downloads/literature/lstc_appendixi_l.pdf?sfvrsn=5f2ab85b_4
# Fit with excel: R=0.9973
# def lambda_PTFE(T):
#     return -2.66883E-14*T**6 + 2.66000E-11*T**5 - 1.04573E-08*T**4 + 2.06223E-06*T**3 - 2.14191E-04*T**2 + 1.12517E-02*T - 8.03724E-03

# Integral thermal conductivity of PTFE according to NIST
# Integral form, just put in the upper and the lower boundaries
# Data and equation range: 4 - 300 K
# Curve fit % error relative to data: 5
def lambda_int_PTFE(T_cold, T_hot):
    result_int = sp.integrate.quad(lambda x: 10**(2.7380 + -30.677*math.log(x,10) + 89.430*math.log(x,10)**2 + -136.99*math.log(x,10)**3
                                              + 124.69*math.log(x,10)**4 + -69.556*math.log(x,10)**5 + 23.320*math.log(x,10)**6
                                              + -4.3135*math.log(x,10)**7 + 0.33829*math.log(x,10)**8), T_cold, T_hot)
    return result_int[0] #W/m

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

# Thermal conductivity of Brass according to lakeshore
# https://www.lakeshore.com/docs/default-source/product-downloads/literature/lstc_appendixi_l.pdf?sfvrsn=5f2ab85b_4
# Fit with excel: R=0.9986
def lambda_Brass(T):
    return -6.60416E-14*T**6 + 1.25445E-10*T**5 - 1.02859E-07*T**4 + 4.11791E-05*T**3 - 9.19342E-03*T**2 + 1.32025*T - 1.16560 #W/(mK)

# Integral thermal conductivity of Brass according to lakeshore
# https://www.lakeshore.com/docs/default-source/product-downloads/literature/lstc_appendixi_l.pdf?sfvrsn=5f2ab85b_4
# Integral form, just put in the lower and upper boundaries
# Fit with excel: R=0.9986
def lambda_int_Brass(T_cold, T_hot):
    result_int = sp.integrate.quad(lambda x: -6.60416E-14*x**6 + 1.25445E-10*x**5 - 1.02859E-07*x**4 + 4.11791E-05*x**3 - 9.19342E-03*x**2 + 1.32025*x - 1.16560, T_cold, T_hot)
    return result_int[0] #W/m

# Thermal conductivity of Silver RRR=30 according to cryocomp

# Fit with excel: R=0.99265
def lambda_Ag(T):
    if 3 <= T < 146:
        return -9.68287E-09*T**6 + 5.08813E-06*T**5 - 1.05617E-03*T**4 + 1.08844E-01*T**3 - 5.66662E+00*T**2 + 1.28638E+02*T - 1.90982E+02 #W/(mK)
    elif 146 <= T:
        return 421 #W/(mK)
    else:
        raise Exception("Temperature out of correlation range!")


# Integral thermal conductivity of Silver RRR=30 according to cryocomp

# Integral form, just put in the lower and upper boundaries
def lambda_int_Ag(T_cold, T_hot):
    # Both temperatures inside the correlation range
    if 3 <= T_cold < 146 and 3 <= T_hot < 146:
        result_int = sp.integrate.quad(lambda x: -9.68287E-09*x**6 + 5.08813E-06*x**5 - 1.05617E-03*x**4 + 1.08844E-01*x**3 - 5.66662E+00*x**2 + 1.28638E+02*x - 1.90982E+02, T_cold, T_hot)
        return result_int[0] #W/m
    # Only the colder temperature in the correlation range and the hoter temperature in the constant range
    elif 3 <= T_cold < 146 and 146 <= T_hot:
        T_border = 146 #K
        result_int = sp.integrate.quad(lambda x: -9.68287E-09*x**6 + 5.08813E-06*x**5 - 1.05617E-03*x**4 + 1.08844E-01*x**3 - 5.66662E+00*x**2 + 1.28638E+02*x - 1.90982E+02, T_cold, T_border)
        return result_int[0] + 421 * (T_hot - T_border) #W/m
    # Both temperatures in the constant range
    elif 3 <= T_cold < 146 and 146 <= T_hot:
        return 421 * (T_hot - T_cold) #W/m
    else:
        raise Exception("Temperatures out of correlation range!")


# Mean thermal conductivity of Silver plated Brass (SBS)
# Weighed with the cross section and the correlations from above
def lambda_SBS(T):
    # Fraction of cross section
    d_Brass = 0.00051 - 2*3.5e-6 #m
    d_Ag    = 0.00051 #m
    A_Brass = sp.pi * 0.25 * d_Brass**2 #m²
    A_Ag    = sp.pi * 0.25 * (d_Ag**2 - d_Brass**2) #m²
    fraction_Ag = A_Ag / (A_Brass + A_Ag)

    return (1-fraction_Ag) * lambda_Brass(T) + fraction_Ag * lambda_Ag(T) #W/(mK)

# Integral mean thermal conductivity of Silver plated Brass (SBS)
# Weighed with the cross section and the correlations from above
def lambda_int_SBS(T_cold, T_hot):
    # Fraction of cross section
    d_Brass = 0.00051 - 2*3.5e-6 #m
    d_Ag    = 0.00051 #m
    A_Brass = sp.pi * 0.25 * d_Brass**2 #m²
    A_Ag    = sp.pi * 0.25 * (d_Ag**2 - d_Brass**2) #m²
    fraction_Ag = A_Ag / (A_Brass + A_Ag)

    return (1-fraction_Ag) * lambda_int_Brass(T_cold, T_hot) + fraction_Ag * lambda_int_Ag(T_cold, T_hot) #W/m



# Conductance of Cu/Cu interface according to Gmelin, Thermal boundary resistance of mechanical contacts between solids at sub-ambient temperatures
# Curve fit of the data with excel
# Unit: W/(m²K)
def conductance_CuCu(T):
    return 1 * ( 5.307104*10**-7*T**4 - 2.669007*10**-4*T**3 + 3.895959*10**-2*T**2 - 2.848293*10**-2*T + 21.20419 )

# Conductance of Cu/PTFE/Cu and Cu/Nylon according to Gmelin, Thermal boundary resistance of mechanical contacts between solids at sub-ambient temperatures
# Logarithmic fit of the data with excel
# Unit: W/(m²K)
def conductance_CuPolymer(T):
    return 38.39937511 * T**0.352276319

# Conductance of SST/IN/SST and SST/Cu according to Gmelin, Thermal boundary resistance of mechanical contacts between solids at sub-ambient temperatures
# Logarithmic fit of the data with excel
# Unit: W/(m²K)
def conductance_SSTCu(T):
    return 5.1329802128714 * T**0.944152524326247

################################################################################

################################################################################
## Function of the model ##
def coax_cable_HL_model(T_TS, pos_TS, length_TS, T_RT, T_CM, length_x, N):

    ## Initialization ##

    # Heat intercept with the TS
    # T_TS = 50 #K
    # pos_TS = 2 #m
    # length_TS = 0.1 #m
    # Boundary conditions
    # T_RT = 300 #K
    # T_CM = 4 #K
    # Length of the coax cable
    # length_x = 3 #m
    # Iteration points in x direction
    # N = 301

    dx = length_x/(N-1) #m

    # Outer diameters of the different layers
    d_SBS = 0.00051 #m
    d_PTFE = 0.00162 #m
    d_Cu = 0.00162 + 2*5e-6 #m
    d_SST = 0.0022 #m

    # Equivalent thickness of the different layers
    dy_SBS = d_SBS/2 #m
    dy_PTFE = d_PTFE/2 - d_SBS/2 #m
    dy_Cu = d_Cu/2 - d_PTFE/2 #m
    dy_SST = d_SST/2 - d_Cu/2 #m

    # Cross section area of the respective layer
    A_SBS = sp.pi * 0.25 * d_SBS**2 #m²
    A_PTFE = sp.pi * 0.25 * (d_PTFE**2 - d_SBS**2) #m²
    A_Cu = sp.pi * 0.25 * (d_Cu**2 - d_PTFE**2) #m²
    A_SST = sp.pi * 0.25 * (d_SST**2 - d_Cu**2) #m²

    # Radial effective area for the heat transfer in rings to be multiplied with the respective conduction coefficent
    A_SBS2PTFE = 2*sp.pi * dx * 1/math.log((dy_SBS)/(0.5*dy_SBS)) #m
    A_PTFE2SBS = 2*sp.pi * dx * 1/math.log((0.5*dy_PTFE+dy_SBS)/(dy_SBS)) #m
    A_PTFE2Cu  = 2*sp.pi * dx * 1/math.log((dy_PTFE+dy_SBS)/(0.5*dy_PTFE+dy_SBS)) #m
    A_Cu2PTFE  = 2*sp.pi * dx * 1/math.log((0.5*dy_Cu+dy_PTFE+dy_SBS)/(dy_PTFE+dy_SBS)) #m
    A_Cu2SST   = 2*sp.pi * dx * 1/math.log((dy_Cu+dy_PTFE+dy_SBS)/(0.5*dy_Cu+dy_PTFE+dy_SBS)) #m
    A_SST2Cu   = 2*sp.pi * dx * 1/math.log((0.5*dy_SST+dy_Cu+dy_PTFE+dy_SBS)/(dy_Cu+dy_PTFE+dy_SBS)) #m
    A_SST2TS   = 2*sp.pi * dx * 1/math.log((dy_SST+dy_Cu+dy_PTFE+dy_SBS)/(0.5*dy_SST+dy_Cu+dy_PTFE+dy_SBS)) #m



    # Initialization of the specific temperatur vector of the different materials
    T_SBS = np.zeros(N, dtype = float) #K
    T_PTFE = np.zeros(N, dtype = float) #K
    T_Cu = np.zeros(N, dtype = float) #K
    T_SST = np.zeros(N, dtype = float) #K
    # Initialization of length direction vector = Calculation points
    # x_lst = [i for i in np.arange(0,length_x+dx,dx)] #m
    x_lst = np.zeros(N, dtype = float)
    for i in range(len(x_lst)):
        x_lst[i] = i * dx

    # Boundary conditions at warm end
    T_SBS[0] = T_RT #K
    T_PTFE[0] = T_RT #K
    T_Cu[0] = T_RT #K
    T_SST[0] = T_RT #K

    # Boundary conditions at cold end
    # Dirichlet boundary condition
    T_SBS[N-1] = T_CM #K
    T_PTFE[N-1] = T_CM #K
    T_Cu[N-1] = T_CM #K
    T_SST[N-1] = T_CM #K

    # Neumann boundary condition
    # q_SBS_0 = 10**4 #W/m²
    # q_PTFE_0 = 10**1 #W/m²
    # q_Cu_0 = 10**4 #W/m²
    # q_SST_0 = 10**4 #W/m²
    # T_SBS[N-1] = T_SBS[N-2] - 1/lambda_SBS(T_SBS[N-2]) * dx * q_SBS_0 #K
    # T_PTFE[N-1] = T_PTFE[N-2] - 1/lambda_PTFE(T_PTFE[N-2]) * dx * q_PTFE_0 #K
    # T_Cu[N-1] = T_Cu[N-2] - 1/lambda_Cu_100(T_Cu[N-2]) * dx * q_Cu_0 #K
    # T_SST[N-1] = T_SST[N-2] - 1/lambda_SST304(T_SST[N-2]) * dx * q_SST_0 #K


    # # Preparing the plot that updates every while loop
    # plt.ion()
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # line1, = ax.plot(x_lst, T_SBS, 'b-')
    # line2, = ax.plot(x_lst, T_PTFE, 'r--')
    # line3, = ax.plot(x_lst, T_Cu, 'y-')
    # line4, = ax.plot(x_lst, T_SST, 'g--')
    # plt.xlabel("length / m")
    # plt.ylabel("temperature / K")


    # Starting conditions for the temperature vectors
    for i in range(1,N-1):

        # Linear Equation between the two end points
        T_SBS[i] = T_SBS[0] + i * (T_SBS[N-1] - T_SBS[0]) / (N-1)
        T_PTFE[i] = T_PTFE[0] + i * (T_PTFE[N-1] - T_PTFE[0]) / (N-1)
        T_Cu[i] = T_Cu[0] + i * (T_Cu[N-1] - T_Cu[0]) / (N-1)
        T_SST[i] = T_SST[0] + i * (T_SST[N-1] - T_SST[0]) / (N-1)


    # Initialization of temperature vectors and factor matrix
    T_vector = np.zeros((4*N,1), dtype = float)
    T_vector_new = np.zeros((4*N,1), dtype = float)
    lambda_vector = np.zeros((4*N), dtype = float)
    conductance_vector = np.zeros((4*N), dtype = float)
    B_vector = np.zeros((4*N,1), dtype = float)
    A_matrix = np.zeros((len(T_vector),len(T_vector)), dtype = float)

    # Initialization for T_vector
    for i in range(N):
        T_vector[i,0] = T_SST[i]
        T_vector[i+N,0] = T_Cu[i]
        T_vector[i+2*N,0] = T_PTFE[i]
        T_vector[i+3*N,0] = T_SBS[i]


    ################################################################################

    ################################################################################
    ## Iteration ##

    epsilon = 1

    while epsilon > 10**-5:

        # lambda_vector as a function of T_vector in W/mK
        for i in range(4*N):
            # lambda for the SST layer
            if  0 <= i <= N-1:
                lambda_vector[i] = lambda_SST304(T_vector[i,0]) #W/mK
            # lambda for the Cu layer
            if  N <= i <= 2*N-1:
                lambda_vector[i] = lambda_Cu_100(T_vector[i,0]) #W/mK
            # lambda for the PTFE layer
            elif 2*N <= i <= 3*N-1:
                lambda_vector[i] = lambda_PTFE(T_vector[i,0]) #W/mK
            # lambda for the SBS layer
            else:
                lambda_vector[i] = lambda_SBS(T_vector[i,0]) #W/mK

        # conductance_vector as a function of T_vector in W/K
        for i in range(4*N):
            # Conductance for the interface between the SST layer and the heat intercept (TS)
            if  0 <= i <= N-1:
                conductance_vector[i] = conductance_SSTCu(0.5 * (T_vector[i,0] + T_TS)) * sp.pi * d_Cu * dx #W/K
            # Conductance for the interface between the SST and the Cu layer
            elif N <= i <= 2*N-1:
                conductance_vector[i] = conductance_SSTCu(0.5 * (T_vector[i,0] + T_vector[i-N,0])) * sp.pi * d_SST * dx #W/K
            # Conductance for the interface between the Cu and the PTFE layer
            elif 2*N <= i <= 3*N-1:
                conductance_vector[i] = conductance_CuPolymer(0.5 * (T_vector[i,0] + T_vector[i-N,0])) * sp.pi * d_PTFE * dx #W/K
            # Conductance for the interface between the PTFE and the SBS layer
            else:
                conductance_vector[i] = conductance_CuPolymer(0.5 * (T_vector[i,0] + T_vector[i-N,0])) * sp.pi * d_SBS * dx #W/K


        # B_vector as a function of T_vector and lambda_vector
        for i in range(4*N):
            # Dirichlet boundary conditions
            if i==0 or i==N-1 or i==N or i==2*N-1 or i==2*N or i==3*N-1 or i==3*N or i==4*N-1:
                B_vector[i,0] = T_vector[i,0]
            # Heat intercept at the TS
            if pos_TS - length_TS/2 <= i*dx <= pos_TS + length_TS/2:
                B_vector[i,0] = T_TS * ( 1/(lambda_vector[i] * A_SST2TS) + 1/conductance_vector[i] )**-1


        # A_matrix as a function of the lambda_vector
        for i in range(4*N):

            # SST layer
            if  0 <= i <= N-1:
                # Boundary conditions
                if i == 0 or i == N-1:
                    A_matrix[i,i] = 1
                # Cells not at the boundary
                else:
                    # Heat along the layer
                    A_matrix[i,i-1] = - 0.5*(lambda_vector[i-1] + lambda_vector[i]) * A_SST * 1/dx
                    A_matrix[i,i+1] = - 0.5*(lambda_vector[i+1] + lambda_vector[i]) * A_SST * 1/dx
                    # Heat to Cu
                    A_matrix[i+N,i] = - ( 1/(lambda_vector[i] * A_SST2Cu) + 1/conductance_vector[i+N] + 1/(lambda_vector[i+N] * A_Cu2SST) )**-1
                    # Coefficients of T_i,i
                    A_matrix[i,i] = ( 0.5*(lambda_vector[i-1] + lambda_vector[i]) * A_SST * 1/dx
                                    + 0.5*(lambda_vector[i+1] + lambda_vector[i]) * A_SST * 1/dx
                                    + ( 1/(lambda_vector[i] * A_SST2Cu) + 1/conductance_vector[i+N] + 1/(lambda_vector[i+N] * A_Cu2SST) )**-1
                                    )
                    # Heat to the TS
                    if pos_TS - length_TS/2 <= i*dx <= pos_TS + length_TS/2:
                        A_matrix[i,i] += ( 1/(lambda_vector[i] * A_SST2TS) + 1/conductance_vector[i] )**-1

            # Cu layer
            elif N <= i <= 2*N-1:
                # Boundary conditions
                if i == N or i == 2*N-1:
                    A_matrix[i,i] = 1
                # Cells not at the boundary
                else:
                    # Heat along the layer
                    A_matrix[i,i-1] = - 0.5*(lambda_vector[i-1] + lambda_vector[i]) * A_Cu * 1/dx
                    A_matrix[i,i+1] = - 0.5*(lambda_vector[i+1] + lambda_vector[i]) * A_Cu * 1/dx
                    # Heat to SST
                    A_matrix[i-N,i] = - ( 1/(lambda_vector[i-N] * A_SST2Cu) + 1/conductance_vector[i] + 1/(lambda_vector[i] * A_Cu2SST) )**-1
                    # Heat to PTFE
                    A_matrix[i+N,i] = - ( 1/(lambda_vector[i] * A_Cu2PTFE) + 1/conductance_vector[i+N] + 1/(lambda_vector[i+N] * A_PTFE2Cu) )**-1
                    # Coefficients of T_i,i
                    A_matrix[i,i] = ( 0.5*(lambda_vector[i-1] + lambda_vector[i]) * A_Cu * 1/dx
                                    + 0.5*(lambda_vector[i+1] + lambda_vector[i]) * A_Cu * 1/dx
                                    + ( 1/(lambda_vector[i-N] * A_SST2Cu) + 1/conductance_vector[i] + 1/(lambda_vector[i] * A_Cu2SST) )**-1
                                    + ( 1/(lambda_vector[i] * A_Cu2PTFE) + 1/conductance_vector[i+N] + 1/(lambda_vector[i+N] * A_PTFE2Cu) )**-1
                                    )

            # PTFE layer
            elif 2*N <= i <= 3*N-1:
                # Boundary conditions
                if i == 2*N or i == 3*N-1:
                    A_matrix[i,i] = 1
                # Cells not at the boundary
                else:
                    # Heat along the layer
                    A_matrix[i,i-1] = - 0.5*(lambda_vector[i-1] + lambda_vector[i]) * A_PTFE * 1/dx
                    A_matrix[i,i+1] = - 0.5*(lambda_vector[i+1] + lambda_vector[i]) * A_PTFE * 1/dx
                    # Heat to SST
                    A_matrix[i-N,i] = - ( 1/(lambda_vector[i-N] * A_Cu2PTFE) + 1/conductance_vector[i] + 1/(lambda_vector[i] * A_PTFE2Cu) )**-1
                    # Heat to SBS
                    A_matrix[i+N,i] = - ( 1/(lambda_vector[i] * A_PTFE2SBS) + 1/conductance_vector[i+N] + 1/(lambda_vector[i+N] * A_SBS2PTFE) )**-1
                    # Coefficients of T_i,i
                    A_matrix[i,i] = ( 0.5*(lambda_vector[i-1] + lambda_vector[i]) * A_PTFE * 1/dx
                                    + 0.5*(lambda_vector[i+1] + lambda_vector[i]) * A_PTFE * 1/dx
                                    + ( 1/(lambda_vector[i-N] * A_Cu2PTFE) + 1/conductance_vector[i] + 1/(lambda_vector[i] * A_PTFE2Cu) )**-1
                                    + ( 1/(lambda_vector[i] * A_PTFE2SBS) + 1/conductance_vector[i+N] + 1/(lambda_vector[i+N] * A_SBS2PTFE) )**-1
                                    )

            # SBS layer
            else:
                # Boundary conditions
                if i == 3*N or i == 4*N-1:
                    A_matrix[i,i] = 1
                # Cells not at the boundary
                else:
                    # Heat along the layer
                    A_matrix[i,i-1] = - 0.5*(lambda_vector[i-1] + lambda_vector[i]) * A_SBS * 1/dx
                    A_matrix[i,i+1] = - 0.5*(lambda_vector[i+1] + lambda_vector[i]) * A_SBS * 1/dx
                    # Heat to PTFE
                    A_matrix[i-N,i] = - ( 1/(lambda_vector[i-N] * A_PTFE2SBS) + 1/conductance_vector[i] + 1/(lambda_vector[i] * A_SBS2PTFE) )**-1
                    # Coefficients of T_i,i
                    A_matrix[i,i] = ( 0.5*(lambda_vector[i-1] + lambda_vector[i]) * A_SBS * 1/dx
                                    + 0.5*(lambda_vector[i+1] + lambda_vector[i]) * A_SBS * 1/dx
                                    + ( 1/(lambda_vector[i-N] * A_PTFE2SBS) + 1/conductance_vector[i] + 1/(lambda_vector[i] * A_SBS2PTFE) )**-1
                                    )


        # Calculation of the T_vector that fullfiles the A*T=B equation
        T_vector_new = np.matmul(np.linalg.inv(A_matrix),B_vector)


        # Check for convergence
        epsilon = 0

        for i in range(len(T_vector_new)):
            epsilon += abs(T_vector[i] - T_vector_new[i])

        epsilon *= 1/sum(T_vector_new)
        print("epsilon: ", epsilon)

        # Relaxation factor
        relax_factor = 0.5
        # Relaxation to the calculated values
        for i in range(len(T_vector_new)):
            T_vector[i] = relax_factor*T_vector_new[i] + (1-relax_factor)*T_vector[i]


        # # Updating the plot every while loop
        # for i in range(N):
        #     T_SST[i] = T_vector_new[i,0]
        #     T_Cu[i] = T_vector_new[i+N,0]
        #     T_PTFE[i] = T_vector_new[i+2*N,0]
        #     T_SBS[i] = T_vector_new[i+3*N,0]

        # line1.set_ydata(T_SBS)
        # line2.set_ydata(T_PTFE)
        # line3.set_ydata(T_Cu)
        # line4.set_ydata(T_SST)
        # fig.canvas.draw()
        # fig.canvas.flush_events()

    ################################################################################

    ################################################################################
    ## Postprocessing ##

    # Retransfer to the temperature vector of the respective layer
    for i in range(N):
        T_SST[i] = T_vector_new[i,0]
        T_Cu[i] = T_vector_new[i+N,0]
        T_PTFE[i] = T_vector_new[i+2*N,0]
        T_SBS[i] = T_vector_new[i+3*N,0]

    # print("T_vector: ", T_PTFE)
    # print("lambda_vector: ", lambda_vector)
    # print(B_vector)
    # print(A_matrix)
    # print(T_vector_new)

    # lambda effective to compare to measured findings of an RF cable
    # effective radial conductivity coefficient for the dielectric, by summing up all the conductances of the boundaries and the conductivty coefficients
    lambda_eff_dielectric = np.zeros((N), dtype = float) #W/(mK)
    for i in range(N):
        lambda_eff_dielectric[i] = ( ( 1/conductance_vector[i+N] + 1/(lambda_vector[i+2*N] * A_Cu2SST) + 1/(lambda_vector[i+2*N] * A_Cu2PTFE)
                                     + 1/conductance_vector[i+2*N] + 1/(lambda_vector[i+3*N] * A_PTFE2Cu) + 1/(lambda_vector[i+3*N] * A_PTFE2SBS)
                                     + 1/conductance_vector[i+3*N])**-1
                                    * (2*sp.pi * dx * 1/math.log((dy_Cu+dy_PTFE+dy_SBS)/(dy_SBS)))**-1 )

    print("lambda_eff_dielectric: ", lambda_eff_dielectric)


    # Actual heat load on the TS
    Q_TS_Coax = 0
    for i in range(N):
        if pos_TS - length_TS/2 <= i*dx <= pos_TS + length_TS/2:
            Q_TS_Coax += (T_SST[i]-T_TS) * ( 1/(lambda_vector[i] * A_SST2TS) + 1/conductance_vector[i] )**-1


    # Actual heat loads on the CM
    # Simple linear Fourier heat transport
    # Use lambda with the higher T to be conseravtive
    q_SBS = 0.5*(lambda_SBS(T_SBS[N-2]) + lambda_SBS(T_SBS[N-1])) * 1/dx * (T_SBS[N-2] - T_SBS[N-1])
    q_PTFE = 0.5*(lambda_PTFE(T_PTFE[N-2]) + lambda_PTFE(T_PTFE[N-1])) * 1/dx * (T_PTFE[N-2] - T_PTFE[N-1])
    q_Cu = 0.5*(lambda_Cu_100(T_Cu[N-2]) + lambda_Cu_100(T_Cu[N-1])) * 1/dx * (T_Cu[N-2] - T_Cu[N-1])
    q_SST = 0.5*(lambda_SST304(T_SST[N-2]) + lambda_SST304(T_SST[N-1])) * 1/dx * (T_SST[N-2] - T_SST[N-1])

    Q_SBS = A_SBS * q_SBS
    Q_PTFE = A_PTFE * q_PTFE
    Q_Cu = A_Cu * q_Cu
    Q_SST = A_SST * q_SST

    Q_CM_Coax = Q_SBS + Q_PTFE + Q_Cu + Q_SST #W

    print("q_SBS: ", q_SBS, "W/m²")
    print("q_PTFE: ", q_PTFE, "W/m²")
    print("q_Cu: ", q_Cu, "W/m²")
    print("q_SST: ", q_SST, "W/m²")

    print("Q_SBS: ", Q_SBS, "W")
    print("Q_PTFE: ", Q_PTFE, "W")
    print("Q_Cu: ", Q_Cu, "W")
    print("Q_SST: ", Q_SST, "W")

    print("Q_TS_Coax: ", Q_TS_Coax, "W")
    print("Q_CM_Coax: ", Q_CM_Coax, "W")

    # Heat load scaled with the inverse of the COP of a Carnot refrigerator
    # Used to weigh the total heat load and optimize the position of the heat intercept
    Q_Carnot_total = (T_RT - T_TS)/T_TS * Q_TS_Coax + (T_RT - T_CM)/T_CM * Q_CM_Coax


    ## Plots

    # fig = plt.subplots(1, 1)
    # plt.plot(x_lst, T_SBS, 'b-', label='T_SBS')
    # plt.plot(x_lst, T_PTFE, 'r--', label='T_PTFE')
    # plt.plot(x_lst, T_Cu, 'y-', label='T_Cu')
    # plt.plot(x_lst, T_SST, 'g--', label='T_SST')
    # plt.xlabel("length / m")
    # plt.ylabel("temperature / K")

    # plt.legend()
    # plt.tight_layout()


    # ## Showing the temperature difference along the length of the cable
    # T_diff_SST = abs(T_SBS - T_SST) #
    # T_diff_Cu = abs(T_SBS - T_Cu) #K
    # T_diff_PTFE = abs(T_SBS - T_PTFE) #K

    # fig = plt.subplots(1, 1)
    # plt.plot(x_lst, T_diff_SST, 'b-', label='T_diff_SBSSST')
    # plt.plot(x_lst, T_diff_Cu, 'r-', label='T_diff_SBSCu')
    # plt.plot(x_lst, T_diff_PTFE, 'g-', label='T_diff_SBSPTFE')
    # plt.xlabel('length / m')
    # plt.ylabel('|dT| SBS / K')

    # plt.legend()
    # plt.tight_layout()


    ################################################################################

    ################################################################################
    ## Saving of the data in a csv file ##

    # Name of the file
    filename = "data_CoaxCable_SST.csv"
    # Data of the file
    fields = ['x', 'T_SST', 'T_Cu', 'T_PTFE', 'T_SBS', 'dT', 'T_average']
    # Restructuring into rows
    rows = np.zeros((N,7), dtype = float)
    for i in range(N):
        rows[i,0] = x_lst[i]
        rows[i,1] = T_SST[i]
        rows[i,2] = T_Cu[i]
        rows[i,3] = T_PTFE[i]
        rows[i,4] = T_SBS[i]
        rows[i,5] = (T_SBS[i] - T_SST[i])
        rows[i,6] = 1/4 * (T_SST[i] + T_Cu[i] + T_PTFE[i] + T_SBS[i])
    # Writing to csv file
    with open('Data/' + filename, 'w', newline='') as csvfile:
        # Creating a csv writer object
        csvwriter = csv.writer(csvfile)
        # Writing the data
        csvwriter.writerow(fields)
        csvwriter.writerows(rows)

    return Q_Carnot_total



################################################################################

################################################################################
## Function calling ##

# Boundary conditions
T_RT = 293 #K
T_CM = 4.5 #K
# Length of the coax cable
length_x = 1.0 #m
# Heat intercept with the TS
T_TS = 50 #K
pos_TS = 0.5*length_x#m
length_TS = 0.1 #m

# Resolution in x direction
N = 301


Q_Carnot_total = coax_cable_HL_model(T_TS, pos_TS, length_TS, T_RT, T_CM, length_x, N)




################################################################################
## Optimization of the heat intercept position ##

# pos_TS_opti = np.array([i for i in np.arange(0.1, length_x, 0.1)])
# Q_Carnot_opti = np.zeros(len(pos_TS_opti), dtype = float)

# for i in range(len(pos_TS_opti)):
#     Q_Carnot_opti[i] = coax_cable_HL_model(T_TS, pos_TS_opti[i], length_TS, T_RT, T_CM, length_x, N)

# pos_TS_opti /= length_x
# print("pos_TS_opti: ", pos_TS_opti)
# print("Q_Carnot_opti: ", Q_Carnot_opti)
#
# plt.plot(pos_TS_opti, Q_Carnot_opti)
# plt.xlabel("position heat intercept / -")
# plt.ylabel("COP scaled heat load / W")
# plt.show()

################################################################################

################################################################################
## Simple Equation: heat load of the coax cables ##
"""
# Assumptions:
# -> setup from above
# -> coax cable Qaxial RG402-CU -> see datasheet
# -> Copper RRR = 100
# -> RT @ 293 K, TS @ 50K and CM @ 4.5K

# Outer diameters of the different layers
d_Brass = 0.00051 - 2*3.5e-6 #m
d_Ag    = 0.00051 #m
d_PTFE  = 0.00162 #m
d_Cu    = 0.00162 + 2*5e-6 #m
d_SST   = 0.0022 #m

# Cross section area of the respective layer
A_Brass = sp.pi * 0.25 * d_Brass**2 #m²
A_Ag    = sp.pi * 0.25 * (d_Ag**2 - d_Brass**2) #m²
A_PTFE  = sp.pi * 0.25 * (d_PTFE**2 - d_Ag**2) #m²
A_Cu    = sp.pi * 0.25 * (d_Cu**2 - d_PTFE**2) #m²
A_SST   = sp.pi * 0.25 * (d_SST**2 - d_Cu**2) #m²

# Heat load according to Fourier conduction
Q_Coax_TS = 1/pos_TS * ( A_Brass*lambda_int_Brass(T_TS,T_RT) + A_Ag*lambda_int_Ag(T_TS,T_RT)
                       + A_PTFE*lambda_int_PTFE(T_TS,T_RT)
                       + A_Cu*lambda_int_Cu_100(T_TS,T_RT) + A_SST*lambda_int_SST304(T_TS,T_RT) ) #W
Q_Coax_CM = 1/(length_x-pos_TS) * ( A_Brass*lambda_int_Brass(T_CM,T_TS) + A_Ag*lambda_int_Ag(T_CM,T_TS)
                                  + A_PTFE*lambda_int_PTFE(T_CM,T_TS)
                                  + A_Cu*lambda_int_Cu_100(T_CM,T_TS) + A_SST*lambda_int_SST304(T_CM,T_TS) ) #W

# Print results: per cable in watt
print("Q_Coax_TS: ", Q_Coax_TS, "W")
print("Q_Coax_CM: ", Q_Coax_CM, "W")

# Print fraction of PTFE heat load
print("Fraction Brass TS: ", 1/pos_TS * A_Brass*lambda_int_Brass(T_TS,T_RT)/Q_Coax_TS *100, "%" )
print("Fraction Brass CM: ", 1/(length_x-pos_TS) * A_Brass*lambda_int_Brass(T_CM,T_TS)/Q_Coax_CM *100, "%" )
print("Fraction Ag TS: ", 1/pos_TS * A_Ag*lambda_int_Ag(T_TS,T_RT)/Q_Coax_TS *100, "%" )
print("Fraction Ag CM: ", 1/(length_x-pos_TS) * A_Ag*lambda_int_Ag(T_CM,T_TS)/Q_Coax_CM *100, "%" )
print("Fraction PTFE TS: ", 1/pos_TS * A_PTFE*lambda_int_PTFE(T_TS,T_RT)/Q_Coax_TS *100, "%" )
print("Fraction PTFE CM: ", 1/(length_x-pos_TS) * A_PTFE*lambda_int_PTFE(T_CM,T_TS)/Q_Coax_CM *100, "%" )
print("Fraction Cu TS: ", 1/pos_TS * A_Cu*lambda_int_Cu_100(T_TS,T_RT)/Q_Coax_TS *100, "%" )
print("Fraction Cu CM: ", 1/(length_x-pos_TS) * A_Cu*lambda_int_Cu_100(T_CM,T_TS)/Q_Coax_CM *100, "%" )
print("Fraction SST TS: ", 1/pos_TS * A_SST*lambda_int_SST304(T_TS,T_RT)/Q_Coax_TS *100, "%" )
print("Fraction SST CM: ", 1/(length_x-pos_TS) * A_SST*lambda_int_SST304(T_CM,T_TS)/Q_Coax_CM *100, "%" )

"""
################################################################################



