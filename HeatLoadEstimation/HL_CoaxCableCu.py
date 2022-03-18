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

# import xlsxwriter
import csv

################################################################################
## Definition of help functions ##

# Thermal conductivity of Cu OFHC RRR=100 according to NIST
# Data and equation range: 4 - 300 K
# Curve fit % error relative to data: 1
def lambda_Cu_100(T):
    return 10**((2.2154 + -0.88068*T**0.5 + 0.29505*T + -0.04831*T**1.5 + 0.003207*T**2)
               /(1 + -0.47461*T**0.5 + 0.13871*T + -0.02043*T**1.5 + 0.001281*T**2))

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
                + -4.3135*math.log(T,10)**7 + 0.33829*math.log(T,10)**8)

# Integral thermal conductivity of PTFE according to NIST
# Integral form, just put in the upper and the lower boundaries
# Data and equation range: 4 - 300 K
# Curve fit % error relative to data: 5
def lambda_int_PTFE(T_cold, T_hot):
    result_int = sp.integrate.quad(lambda x: 10**(2.7380 + -30.677*math.log(x,10) + 89.430*math.log(x,10)**2 + -136.99*math.log(x,10)**3
                                              + 124.69*math.log(x,10)**4 + -69.556*math.log(x,10)**5 + 23.320*math.log(x,10)**6
                                              + -4.3135*math.log(x,10)**7 + 0.33829*math.log(x,10)**8), T_cold, T_hot)
    return result_int[0] #W/m

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
    d_Cu_inner = 0.00092 #m
    d_PTFE = 0.00299 #m
    d_Cu_outer = 0.00358 #m

    # Equivalent thickness of the different layers
    dy_Cu_inner = d_Cu_inner/2 #m
    dy_PTFE = d_PTFE/2 - d_Cu_inner/2 #m
    dy_Cu_outer = d_Cu_outer/2 - d_PTFE/2 #m

    # Cross section area of the respective layer
    A_Cu_inner = sp.pi * 0.25 * d_Cu_inner**2 #m²
    A_PTFE = sp.pi * 0.25 * (d_PTFE**2 - d_Cu_inner**2) #m²
    A_Cu_outer = sp.pi * 0.25 * (d_Cu_outer**2 - d_PTFE**2) #m²

    # Radial effective area for the heat transfer in rings to be multiplied with the respective conduction coefficent
    A_Cuo2TS = 2*sp.pi * dx * 1/math.log((dy_Cu_outer+dy_PTFE+dy_Cu_inner)/(0.5*dy_Cu_outer+dy_PTFE+dy_Cu_inner)) #m
    A_Cuo2PTFE = 2*sp.pi * dx * 1/math.log((0.5*dy_Cu_outer+dy_PTFE+dy_Cu_inner)/(dy_PTFE+dy_Cu_inner)) #m
    A_PTFE2Cuo = 2*sp.pi * dx * 1/math.log((dy_PTFE+dy_Cu_inner)/(0.5*dy_PTFE+dy_Cu_inner)) #m
    A_PTFE2Cui = 2*sp.pi * dx * 1/math.log((0.5*dy_PTFE+dy_Cu_inner)/dy_Cu_inner) #m
    A_Cui2PTFE = 2*sp.pi * dx * 1/math.log(dy_Cu_inner/(0.5*dy_Cu_inner)) #m


    # Initialization of the specific temperatur vector of the different materials
    T_Cu_inner = np.zeros(N, dtype = float) #K
    T_PTFE = np.zeros(N, dtype = float) #K
    T_Cu_outer = np.zeros(N, dtype = float) #K
    # Initialization of length direction vector = Calculation points
    # x_lst = [i for i in np.arange(0,length_x+dx,dx)] #m
    x_lst = np.zeros(N, dtype = float)
    for i in range(len(x_lst)):
        x_lst[i] = i * dx

    # Boundary conditions at warm end
    T_Cu_inner[0] = T_RT #K
    T_Cu_outer[0] = T_RT #K
    T_PTFE[0] = T_RT #K

    # Boundary conditions at cold end
    # Dirichlet boundary condition
    T_Cu_inner[N-1] = T_CM #K
    T_PTFE[N-1] = T_CM #K
    T_Cu_outer[N-1] = T_CM #K

    # Neumann boundary condition
    # q_Cu_inner_0 = 10**4 #W/m²
    # q_PTFE_0 = 10**1 #W/m²
    # q_Cu_outer_0 = 10**4 #W/m²
    # T_Cu_inner[N-1] = T_Cu_inner[N-2] - 1/lambda_Cu_100(T_Cu_inner[N-2]) * dx * q_Cu_inner_0 #K
    # T_PTFE[N-1] = T_PTFE[N-2] - 1/lambda_PTFE(T_PTFE[N-2]) * dx * q_PTFE_0 #K
    # T_Cu_outer[N-1] = T_Cu_outer[N-2] - 1/lambda_Cu_100(T_Cu_outer[N-2]) * dx * q_Cu_outer_0 #K


    # # Preparing the plot
    # plt.ion()
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # line1, = ax.plot(x_lst, T_Cu_inner, 'b-')
    # line2, = ax.plot(x_lst, T_PTFE, 'r-')
    # line3, = ax.plot(x_lst, T_Cu_outer, 'b--')
    # plt.xlabel("length / m")
    # plt.ylabel("temperature / K")


    # Starting conditions for the temperature vectors
    for i in range(1,N-1):

        # # Linear Equation between the two end points
        T_Cu_inner[i] = T_Cu_inner[0] + i * (T_Cu_inner[N-1] - T_Cu_inner[0]) / (N-1)
        T_PTFE[i] = T_PTFE[0] + i * (T_PTFE[N-1] - T_PTFE[0]) / (N-1)
        T_Cu_outer[i] = T_Cu_outer[0] + i * (T_Cu_outer[N-1] - T_Cu_outer[0]) / (N-1)

        # Polynomic fit for the temperature distribution without radial heat transfer
        # x = i * dx
        # T_Cu_inner[i] = -7.6063*x**5 + 49.05*x**4 - 94.206*x**3 + 62.473*x**2 - 147.2*x + 300.54
        # T_PTFE[i] = -1.847*x**5 + 11.599*x**4 - 26.292*x**3 + 25.469*x**2 - 100.8*x + 300.69
        # T_Cu_outer[i] = -7.6063*x**5 + 49.05*x**4 - 94.206*x**3 + 62.473*x**2 - 147.2*x + 300.54


    # Initialization of temperature vectors and factor matrix
    T_vector = np.zeros((3*N,1), dtype = float)
    T_vector_new = np.zeros((3*N,1), dtype = float)
    lambda_vector = np.zeros((3*N), dtype = float)
    conductance_vector = np.zeros((3*N), dtype = float)
    B_vector = np.zeros((3*N,1), dtype = float)
    A_matrix = np.zeros((len(T_vector),len(T_vector)), dtype = float)

    # Initialization for T_vector
    for i in range(N):
        T_vector[i,0] = T_Cu_outer[i]
        T_vector[i+N,0] = T_PTFE[i]
        T_vector[i+2*N,0] = T_Cu_inner[i]


    ################################################################################

    ################################################################################
    ## Iteration ##

    epsilon = 1

    while epsilon > 10**-5:

        # lambda_vector as a function of T_vector
        for i in range(3*N):
            # lambda for the outer Copper layer
            if  0 <= i <= N-1:
                lambda_vector[i] = lambda_Cu_100(T_vector[i,0])
            # lambda for the PTFE layer
            elif N <= i <= 2*N-1:
                lambda_vector[i] = lambda_PTFE(T_vector[i,0])
            # lambda for the inner Copper layer
            else:
                lambda_vector[i] = lambda_Cu_100(T_vector[i,0])

        # conductance_vector as a function of T_vector in W/K
        for i in range(3*N):
            # Conductance for the interface between the outer Copper layer and the heat intercept
            if  0 <= i <= N-1:
                conductance_vector[i] = conductance_CuCu(0.5 * (T_vector[i,0] + T_TS)) * sp.pi * d_Cu_outer * dx
            # Conductance for the interface between the outer Copper and the PTFE layer
            elif N <= i <= 2*N-1:
                conductance_vector[i] = conductance_CuPolymer(0.5 * (T_vector[i,0] + T_vector[i-N,0])) * sp.pi * d_PTFE * dx
            # Conductance for the interface between the inner Copper and the PTFE layer
            else:
                conductance_vector[i] = conductance_CuPolymer(0.5 * (T_vector[i,0] + T_vector[i-N,0])) * sp.pi * d_Cu_inner * dx


        # B_vector as a function of T_vector and lambda_vector
        for i in range(3*N):
            # Dirichlet boundary conditions
            if i==0 or i==N-1 or i==N or i==2*N-1 or i==2*N or i==3*N-1:
                B_vector[i,0] = T_vector[i,0]
            # Heat intercept at the TS
            if pos_TS - length_TS/2 <= i*dx <= pos_TS + length_TS/2:
                B_vector[i,0] = T_TS * ( 1/(lambda_vector[i] * A_Cuo2TS) + 1/conductance_vector[i] )**-1


        # A_matrix as a function of the lambda_vector
        for i in range(3*N):

            # Cu_outer layer
            if  0 <= i <= N-1:
                # Boundary conditions
                if i == 0:
                    A_matrix[i,i] = 1
                elif i == N-1:
                    A_matrix[i,i] = 1
                else:
                    A_matrix[i,i-1] = - 0.5*(lambda_vector[i-1] + lambda_vector[i]) * A_Cu_outer * 1/dx
                    A_matrix[i,i+1] = - 0.5*(lambda_vector[i+1] + lambda_vector[i]) * A_Cu_outer * 1/dx
                    # Heat to PTFE
                    A_matrix[i+N,i] = - ( 1/(lambda_vector[i] * A_Cuo2PTFE) + 1/conductance_vector[i+N] + 1/(lambda_vector[i+N] * A_PTFE2Cuo) )**-1
                    A_matrix[i,i] = ( 0.5*(lambda_vector[i-1] + lambda_vector[i]) * A_Cu_outer * 1/dx
                                    + 0.5*(lambda_vector[i+1] + lambda_vector[i]) * A_Cu_outer * 1/dx
                                    + ( 1/(lambda_vector[i] * A_Cuo2PTFE) + 1/conductance_vector[i+N] + 1/(lambda_vector[i+N] * A_PTFE2Cuo) )**-1
                                    )
                    # Heat to the TS
                    if pos_TS - length_TS/2 <= i*dx <= pos_TS + length_TS/2:
                        A_matrix[i,i] += ( 1/(lambda_vector[i] * A_Cuo2TS) + 1/conductance_vector[i] )**-1


            # PTFE layer
            elif N <= i <= 2*N-1:

                # Boundary conditions
                if i == N:
                    A_matrix[i,i] = 1
                elif i == 2*N-1:
                    A_matrix[i,i] = 1
                # Cells not at the boundary
                else:
                    A_matrix[i,i-1] = - 0.5*(lambda_vector[i-1] + lambda_vector[i]) * A_PTFE * 1/dx
                    A_matrix[i,i+1] = - 0.5*(lambda_vector[i+1] + lambda_vector[i]) * A_PTFE * 1/dx
                    # Heat to outer Copper
                    A_matrix[i-N,i] = - ( 1/(lambda_vector[i-N] * A_Cuo2PTFE) + 1/conductance_vector[i] + 1/(lambda_vector[i] * A_PTFE2Cuo) )**-1
                    # Heat to inner Copper
                    A_matrix[i+N,i] = - ( 1/(lambda_vector[i] * A_PTFE2Cui) + 1/conductance_vector[i+N] + 1/(lambda_vector[i+N] * A_Cui2PTFE) )**-1
                    A_matrix[i,i] = ( 0.5*(lambda_vector[i-1] + lambda_vector[i]) * A_PTFE * 1/dx
                                    + 0.5*(lambda_vector[i+1] + lambda_vector[i]) * A_PTFE * 1/dx
                                    + ( 1/(lambda_vector[i-N] * A_Cuo2PTFE) + 1/conductance_vector[i] + 1/(lambda_vector[i] * A_PTFE2Cuo) )**-1
                                    + ( 1/(lambda_vector[i] * A_PTFE2Cui) + 1/conductance_vector[i+N] + 1/(lambda_vector[i+N] * A_Cui2PTFE) )**-1
                                    )


            # Cu_inner layer
            else:
                # Boundary conditions
                if i == 2*N:
                    A_matrix[i,i] = 1
                elif i == 3*N-1:
                    A_matrix[i,i] = 1
                else:
                    A_matrix[i,i-1] = - 0.5*(lambda_vector[i-1] + lambda_vector[i]) * A_Cu_inner * 1/dx
                    A_matrix[i,i+1] = - 0.5*(lambda_vector[i+1] + lambda_vector[i]) * A_Cu_inner * 1/dx
                    # Heat to PTFE
                    A_matrix[i-N,i] = - ( 1/(lambda_vector[i-N] * A_PTFE2Cui) + 1/conductance_vector[i] + 1/(lambda_vector[i] * A_Cui2PTFE) )**-1
                    A_matrix[i,i] = ( 0.5*(lambda_vector[i-1] + lambda_vector[i]) * A_Cu_inner * 1/dx
                                    + 0.5*(lambda_vector[i+1] + lambda_vector[i]) * A_Cu_inner * 1/dx
                                    + ( 1/(lambda_vector[i-N] * A_PTFE2Cui) + 1/conductance_vector[i] + 1/(lambda_vector[i] * A_Cui2PTFE) )**-1
                                    )


        # Calculation of the T_vector that fullfiles the A*T=B equation
        T_vector_new = np.matmul(np.linalg.inv(A_matrix),B_vector)


        # Retransfer to the temperature vector of the respective layer
        for i in range(N):
            T_Cu_outer[i] = T_vector_new[i,0]
            T_PTFE[i] = T_vector_new[i+N,0]
            T_Cu_inner[i] = T_vector_new[i+2*N,0]




        # Check for convergence
        epsilon = 0

        for i in range(len(T_vector_new)):
            epsilon += abs(T_vector[i] - T_vector_new[i])

        epsilon *= 1/sum(T_vector_new)
        print("epsilon: ", epsilon)

        # Relaxation factor
        relax_factor = 0.2
        # Relaxation to the calculated values
        for i in range(len(T_vector_new)):
            T_vector[i] = relax_factor*T_vector_new[i] + (1-relax_factor)*T_vector[i]



        # # Updating the plot every while loop
        # line1.set_ydata(T_Cu_inner)
        # line2.set_ydata(T_PTFE)
        # line3.set_ydata(T_Cu_outer)
        # fig.canvas.draw()
        # fig.canvas.flush_events()

    ################################################################################

    ################################################################################
    ## Postprocessing ##

    print("T_vector: ", T_PTFE)
    # print("lambda_vector: ", lambda_vector)
    # print(B_vector)
    # print(A_matrix)
    # print(T_vector_new)

    # lambda effective to compare to measured findings of an RF cable
    lambda_eff_outer = np.zeros((N), dtype = float)
    lambda_eff_inner = np.zeros((N), dtype = float)
    lambda_eff_CuCu = np.zeros((N), dtype = float)
    for i in range(N):
        lambda_eff_outer[i] = ( ( 1/(lambda_vector[i] * A_Cuo2PTFE) + 1/conductance_vector[i+N] + 1/(lambda_vector[i+N] * A_PTFE2Cuo) )**-1
                              * (2*sp.pi * dx * 1/math.log((0.5*dy_Cu_outer+dy_PTFE+dy_Cu_inner)/(0.5*dy_PTFE+dy_Cu_inner)))**-1 )
        lambda_eff_inner[i] = ( ( 1/(lambda_vector[i+N] * A_PTFE2Cui) + 1/conductance_vector[i+2*N] + 1/(lambda_vector[i+2*N] * A_Cui2PTFE) )**-1
                              * (2*sp.pi * dx * 1/math.log((0.5*dy_PTFE+dy_Cu_inner)/(0.5*dy_Cu_inner)))**-1 )
        lambda_eff_CuCu[i] = ( ( 1/conductance_vector[i+N] + 1/(lambda_vector[i+N] * A_PTFE2Cuo) + 1/(lambda_vector[i+N] * A_PTFE2Cui) + 1/conductance_vector[i+2*N] )**-1
                              * (2*sp.pi * dx * 1/math.log((dy_PTFE+dy_Cu_inner)/(dy_Cu_inner)))**-1 )

    # print("lambda_eff_outer: ", lambda_eff_outer)
    # print("lambda_eff_inner: ", lambda_eff_inner)
    print("lambda_eff_CuCu: ", lambda_eff_CuCu)


    # Actual heat load on the TS
    Q_TS_Coax = 0
    for i in range(N):
        if pos_TS - length_TS/2 <= i*dx <= pos_TS + length_TS/2:
            Q_TS_Coax += (T_Cu_outer[i]-T_TS) * ( 1/(lambda_vector[i] * A_Cuo2TS) + 1/conductance_vector[i] )**-1


    # Actual heat loads on the CM
    # Simple linear Fourier heat transport
    # Use lambda with the higher T to get be conseravtive
    q_Cu_inner = 0.5*(lambda_Cu_100(T_Cu_inner[N-2]) + lambda_Cu_100(T_Cu_inner[N-1])) * 1/dx * (T_Cu_inner[N-2] - T_Cu_inner[N-1])
    q_PTFE = 0.5*(lambda_PTFE(T_PTFE[N-2]) + lambda_PTFE(T_PTFE[N-1])) * 1/dx * (T_PTFE[N-2] - T_PTFE[N-1])
    q_Cu_outer = 0.5*(lambda_Cu_100(T_Cu_outer[N-2]) + lambda_Cu_100(T_Cu_outer[N-1])) * 1/dx * (T_Cu_outer[N-2] - T_Cu_outer[N-1])

    Q_Cu_inner = A_Cu_inner * q_Cu_inner
    Q_PTFE = A_PTFE * q_PTFE
    Q_Cu_outer = A_Cu_outer * q_Cu_outer

    print("q_Cu_inner: ", q_Cu_inner)
    print("q_PTFE: ", q_PTFE)
    print("q_Cu_outer: ", q_Cu_outer)

    print("Q_Cu_inner: ", Q_Cu_inner)
    print("Q_PTFE: ", Q_PTFE)
    print("Q_Cu_outer: ", Q_Cu_outer)

    print("Q_TS_Coax: ", Q_TS_Coax)
    print("Q_CM_Coax: ", Q_Cu_inner + Q_PTFE + Q_Cu_outer)

    # Heat load scaled with the inverse of the COP of a Carnot refrigerator
    # Used to weigh the total heat load and optimize the position of the heat intercept
    Q_Carnot_total = (T_RT - T_TS)/T_TS * Q_TS_Coax + (T_RT - T_CM)/T_CM * (Q_Cu_inner + Q_PTFE + Q_Cu_outer)




    ## Plot

    # Preparing the plot
    # fig1 = plt.subplots(1, 1)
    # ax1 = fig1.add_subplot(111)
    # line1, = ax1.plot(x_lst, T_Cu_inner, 'b-')
    # line2, = ax1.plot(x_lst, T_PTFE, 'r-')
    # line3, = ax1.plot(x_lst, T_Cu_outer, 'b--')
    # plt.xlabel("length / m")
    # plt.ylabel("temperature / K")

    # fig = plt.subplots(1, 1)
    # plt.plot(x_lst, T_Cu_inner, 'b:', label='T_Cu_inner')
    # plt.plot(x_lst, T_PTFE, 'r-', label='T_PTFE')
    # plt.plot(x_lst, T_Cu_outer, 'b--', label='T_Cu_outer')
    # plt.xlabel("length / m")
    # plt.ylabel("temperature / K")

    # plt.legend()
    # plt.tight_layout()


    # ## Showing the temperature difference along the length of the cable
    # T_diff_outer = abs(T_Cu_inner - T_Cu_outer) #
    # T_diff_inner = abs(T_Cu_inner - T_PTFE) #K

    # fig = plt.subplots(1, 1)
    # plt.plot(x_lst, T_diff_outer, 'b-', label='T_diff_CuCu')
    # plt.plot(x_lst, T_diff_inner, 'r-', label='T_diff_CuPTFE')
    # plt.xlabel('length / m')
    # plt.ylabel('|dT| Inner Cu / K')

    # plt.legend()
    # plt.tight_layout()


    ################################################################################

    ################################################################################
    ## Saving of the data ##

    # # Filter for the data
    # #N_Filter = 101
    # N_Filter = N
    # numpy_data = np.zeros((4,N_Filter), dtype = float)
    # for i in range(0, N_Filter):
    #     numpy_data[0,i] = x_lst[int((N-1)/(N_Filter-1) * i)]
    #     numpy_data[1,i] = T_Cu_inner[int((N-1)/(N_Filter-1) * i)]
    #     numpy_data[2,i] = T_PTFE[int((N-1)/(N_Filter-1) * i)]
    #     numpy_data[3,i] = T_Cu_outer[int((N-1)/(N_Filter-1) * i)]
    #
    #
    # # Create an xlsx sheet
    # outWorkbook = xlsxwriter.Workbook("N301_Dirichlet_at2m.xlsx")
    # outSheet = outWorkbook.add_worksheet()
    #
    # # Writing of the headers
    # outSheet.write("A1", "x")
    # outSheet.write("B1", "T_Cu_inner")
    # outSheet.write("C1", "T_PTFE")
    # outSheet.write("D1", "T_Cu_outer")
    # outSheet.write("A2", "m")
    # outSheet.write("B2", "K")
    # outSheet.write("C2", "K")
    # outSheet.write("D2", "K")
    #
    # # Writing of the data
    # for i in range(0, N_Filter):
    #     outSheet.write(i+2, 0, numpy_data[0,i])
    #     outSheet.write(i+2, 1, numpy_data[1,i])
    #     outSheet.write(i+2, 2, numpy_data[2,i])
    #     outSheet.write(i+2, 3, numpy_data[3,i])
    #
    # outWorkbook.close()

    ## Save the data in a csv file

    # Name of the file
    filename = "data_CoaxCable_Cu.csv"
    # Data of the file
    fields = ['x', 'T_Cu_outer', 'T_PTFE', 'T_Cu_inner', 'dT', 'T_average']
    # Restructuring into rows
    rows = np.zeros((N,6), dtype = float)
    for i in range(N):
        rows[i,0] = x_lst[i]
        rows[i,1] = T_Cu_outer[i]
        rows[i,2] = T_PTFE[i]
        rows[i,3] = T_Cu_inner[i]
        rows[i,4] = (T_Cu_inner[i] - T_Cu_outer[i])
        rows[i,5] = 1/3 * (T_Cu_outer[i] + T_PTFE[i] + T_Cu_inner[i])
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
length_x = 3.41 #m
# Heat intercept with the TS
T_TS = 50 #K
pos_TS = 0.25*length_x #m
length_TS = 0.1 #m

# Resolution in x direction
N = 301


Q_Carnot_total = coax_cable_HL_model(T_TS, pos_TS, length_TS, T_RT, T_CM, length_x, N)




################################################################################
## Optimization of the heat intercept position ##

# pos_TS_opti = np.array([i for i in np.arange(0.1, length_x, 0.1)])
# Q_Carnot_opti = np.zeros(len(pos_TS_opti), dtype = float)
#
# for i in range(len(pos_TS_opti)):
#     Q_Carnot_opti[i] = coax_cable_HL_model(T_TS, pos_TS_opti[i], length_TS, T_RT, T_CM, length_x, N)
#
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
d_Cu_inner = 0.00092 #m
d_PTFE = 0.00299 #m
d_Cu_outer = 0.00358 #m

# Cross section area of the respective layer
A_Cu_inner = sp.pi * 0.25 * d_Cu_inner**2 #m²
A_PTFE = sp.pi * 0.25 * (d_PTFE**2 - d_Cu_inner**2) #m²
A_Cu_outer = sp.pi * 0.25 * (d_Cu_outer**2 - d_PTFE**2) #m²

# Heat load according to Fourier conduction
Q_Coax_TS = 1/pos_TS * ((A_Cu_inner+A_Cu_outer)*lambda_int_Cu_100(T_TS,T_RT) + A_PTFE*lambda_int_PTFE(T_TS,T_RT)) #W
Q_Coax_CM = 1/(length_x-pos_TS) * ((A_Cu_inner+A_Cu_outer)*lambda_int_Cu_100(T_CM,T_TS) + A_PTFE*lambda_int_PTFE(T_CM,T_TS)) #W

# Print results: per cable in watt
print("Q_Coax_TS: ", Q_Coax_TS, "W")
print("Q_Coax_CM: ", Q_Coax_CM, "W")

# Print fraction of PTFE heat load
print("Fraction PTFE TS: ", 1/pos_TS * A_PTFE*lambda_int_PTFE(T_TS,T_RT)/Q_Coax_TS *100, "%" )
print("Fraction PTFE CM: ", 1/(length_x-pos_TS) * A_PTFE*lambda_int_PTFE(T_CM,T_TS)/Q_Coax_CM *100, "%" )
"""
################################################################################



