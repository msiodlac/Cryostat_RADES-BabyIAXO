"""
Heat load estimation for the tuning rod.

Numerical approach to simulate the complicated behavior between solid conduction from room temperature,
thermal radiation from the TS and solid conduction to the cavities.


19.07.2021
"""

import math
import numpy as np
import scipy as sp
from scipy import constants
import matplotlib.pyplot as plt

# import sys
# PATH = r'C:\Users\msiodlac\cernbox\Documents\0_cern_msiodlac\BoyanModels\py-plot';
# sys.path.insert(0, PATH)
# from plotter import plotter

from matplotlib import cm

import xlsxwriter
import csv

################################################################################
## Definition of functions ##

# Function: calculating the heat load of thermal radiation for a cold surface which is enclosed by a warm surface
# According to "Tieftemperaturtechnologie" by Frey and Ekin  for diffuse reflection
def Q_therm_rad(A_hot, A_cold, T_hot, T_cold, epsilon_hot, epsilon_cold):
    # Emission coefficent for an enclosed vessel
    epsilon_effective = 1/ (1/epsilon_cold + A_cold/A_hot * (1/epsilon_hot - 1)) #-
    # Simple equation with an effective emission coefficent
    Q_rad = epsilon_effective * sp.constants.Stefan_Boltzmann * A_cold * (T_hot**4 - T_cold**4) #W
    return Q_rad #W

# Thermal conductivity of G10 normal direction according to NIST
# https://trc.nist.gov/cryogenics/materials/materialproperties.htm
# Data range: 4 - 300 K
# Equation range: 10 - 300 K
# Curve fit % error relative to data: 5
def lambda_G10_normal(T):
    return 10**(-4.1236 + 13.788*math.log(T,10) + -26.068*math.log(T,10)**2 + 26.272*math.log(T,10)**3
                + -14.663*math.log(T,10)**4 + 4.4954*math.log(T,10)**5 + -0.6905*math.log(T,10)**6
                + 0.0397*math.log(T,10)**7 + 0*math.log(T,10)**8) #W/mK

# # Thermal conductivity of SST304 according to NIST
# # https://trc.nist.gov/cryogenics/materials/materialproperties.htm
# # Data range: 4 - 300 K
# # Equatin range: 1 - 300 K
# # Curve fit % error relative to data: 2
# def lambda_SST304(T):
#     return 10**(-1.4087 + 1.3982*math.log(T,10) + 0.2543*math.log(T,10)**2 + -0.6260*math.log(T,10)**3
#                 + 0.2334*math.log(T,10)**4 + 0.4256*math.log(T,10)**5 + -0.4658*math.log(T,10)**6
#                 + 0.1650*math.log(T,10)**7 + -0.0199*math.log(T,10)**8) #W/mK


################################################################################

################################################################################
## Function of the model ##
def tuning_rod_HL_model(N, T_RT, T_TS, pos_TS, length_TS, T_CM, pos_cav1_1, pos_cav1_2, pos_cav2_1, pos_cav2_2, pos_cav3_1, pos_cav3_2, length_cav, epsilon_TS, epsilon_Rod, relaxation_factor):

    ## Initialization ##

    # # Iteration points in x direction
    # N = 500
    # # Boundary conditions
    # T_RT = 300 #K
    # # Heat intercept with the TS
    # T_TS = 50 #K
    # pos_TS = 0.3 #m
    # length_TS = 0.1 #m
    # # Heat intercept at the cavities
    # T_CM = 10 #K
    # pos_cav1 = 2 #m
    # pos_cav2 = 5 #m
    # pos_cav3 = 8 #m
    # length_cav = 0.05 #m
    # # Relaxation factor for stability of the simulation
    # relaxation_factor = 0.01


    # Geometry data of the rod
    diameter_rod = 0.012 #m
    thickness_rod = 0.003 #m
    length_x = 13 #m
    dx = length_x/(N-1) #m

    # Cross section and circumference of the hollow rod
    A_rod_cross = sp.pi * 0.25 * (diameter_rod**2-(diameter_rod-2*thickness_rod)**2) #m²
    circumference_rod = sp.pi * diameter_rod #m

    # Geometry data of the TS
    diameter_TS_bore = 0.65 #m
    circumference_TS = sp.pi * diameter_TS_bore #m
    position_bore = length_x - 10.0 #m
    diameter_TS_transfer = 0.15 #m
    circumference_TS_transfer = sp.pi * diameter_TS_transfer #m



    # Initialization of the specific temperature and lambda vector
    T_rod = np.zeros((N,1), dtype = float) #K
    T_rod_new = np.zeros((N,1), dtype = float) #K
    lambda_vector = np.zeros((N), dtype = float) #W/mK
    Q_rad_vector = np.zeros((N), dtype = float) #W
    B_vector = np.zeros((N,1), dtype = float)
    A_matrix = np.zeros((N,N), dtype = float)


    # Initialization of external heat source vectors
    # Q_rad = np.zeros(N, dtype = float) #W
    # Q_intercept = np.zeros(N, dtype = float) #W
    # Q_cav1 = np.zeros(N, dtype = float) #W
    # Q_cav2 = np.zeros(N, dtype = float) #W
    # Q_cav3 = np.zeros(N, dtype = float) #W


    # Boundary conditions at warm end
    T_rod[0,0] = T_RT #K

    # Boundary conditions at cold end
    # Neumann boundary condition
    # Adiabatic end temperature
    q_rod_N = 0 #W/m²

    # Preparing the plot
    # Initialization of length direction vector = Calculation points
    x_lst = np.zeros(N, dtype = float)
    for i in range(len(x_lst)):
        x_lst[i] = i * dx
    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1, = ax.plot(x_lst, T_rod, 'b-')
    plt.xlabel("length / m")
    plt.ylabel("temperature / K")


    # Starting conditions for the temperature vectors
    for i in range(1,N):
        T_rod[i,0] = T_TS


    # Convergence criterium initialization
    epsilon = 1

    ################################################################################

    ################################################################################
    ## Iteration ##

    while epsilon > 10**-4:

        # lambda_vector as a function of T_rod
        for i in range(N):
            lambda_vector[i] = lambda_G10_normal(T_rod[i,0])

        # Q_rad_vector as a function of T_rod
        for i in range(N):
            # rod enclosed by the service box vacuum chamber
            if i*dx < pos_TS - length_TS/2:
                Q_rad_vector[i] = Q_therm_rad(circumference_TS_transfer*dx, circumference_rod*dx, T_RT, T_rod[i,0], epsilon_TS, epsilon_Rod)
            # rod enclosed by the service box or the transfer line TS
            elif i*dx < position_bore:
                Q_rad_vector[i] = Q_therm_rad(circumference_TS_transfer*dx, circumference_rod*dx, T_TS, T_rod[i,0], epsilon_TS, epsilon_Rod)
            # rod enclosed by the bore TS
            else:
                Q_rad_vector[i] = Q_therm_rad(circumference_TS*dx, circumference_rod*dx, T_TS, T_rod[i,0], epsilon_TS, epsilon_Rod)

        # print("Q_rad_vector: ", Q_rad_vector)
        # print("T_rod: ", T_rod)

        # B_vector as a function of T_rod, lambda_vector and Q_rad_vector
        for i in range(N):
            # Dirichlet boundary conditions
            if i == 0 or i == N-1:
                B_vector[i,0] = T_rod[i,0]
            else:
                # Heat load through thermal radiation
                B_vector[i,0] = Q_rad_vector[i]
                # if i*dx < pos_TS + length_TS/2:
                #     B_vector[i,0] = 0
                # if pos_TS + length_TS/2 <= i*dx:
                #     B_vector[i,0] = Q_rad_vector[i]

                # Heat intercept at the TS
                if pos_TS - length_TS/2 <= i*dx <= pos_TS + length_TS/2:
                    B_vector[i,0] += T_TS * lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))

                # Heat intercept at the first cavity
                if pos_cav1_1 - length_cav/2 <= i*dx <= pos_cav1_1 + length_cav/2:
                    B_vector[i,0] += T_CM * lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))
                if pos_cav1_2 - length_cav/2 <= i*dx <= pos_cav1_2 + length_cav/2:
                    B_vector[i,0] += T_CM * lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))
                # Heat intercept at the second cavity
                if pos_cav2_1 - length_cav/2 <= i*dx <= pos_cav2_1 + length_cav/2:
                    B_vector[i,0] += T_CM * lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))
                if pos_cav2_2 - length_cav/2 <= i*dx <= pos_cav2_2 + length_cav/2:
                    B_vector[i,0] += T_CM * lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))
                # Heat intercept at the third cavity
                if pos_cav3_1 - length_cav/2 <= i*dx <= pos_cav3_1 + length_cav/2:
                    B_vector[i,0] += T_CM * lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))
                if pos_cav3_2 - length_cav/2 <= i*dx <= pos_cav3_2 + length_cav/2:
                    B_vector[i,0] += T_CM * lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))


        # A_matrix as a function of the lambda_vector
        for i in range(N):

            # Boundary conditions
            if i == 0:
                A_matrix[i,i] = 1
            elif i == N-1:
                A_matrix[i,i] = 1
            else:
                A_matrix[i,i-1] = - 0.5*(lambda_vector[i-1] + lambda_vector[i]) * A_rod_cross * 1/dx
                A_matrix[i,i+1] = - 0.5*(lambda_vector[i+1] + lambda_vector[i]) * A_rod_cross * 1/dx

                A_matrix[i,i] = ( 0.5*(lambda_vector[i-1] + lambda_vector[i]) * A_rod_cross * 1/dx
                                + 0.5*(lambda_vector[i+1] + lambda_vector[i]) * A_rod_cross * 1/dx
                                )

                # Heat intercept at the TS
                if pos_TS - length_TS/2 <= i*dx <= pos_TS + length_TS/2:
                    A_matrix[i,i] += lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))

                # Heat intercept at the first cavity
                if pos_cav1_1 - length_cav/2 <= i*dx <= pos_cav1_1 + length_cav/2:
                    A_matrix[i,i] += lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))
                if pos_cav1_2 - length_cav/2 <= i*dx <= pos_cav1_2 + length_cav/2:
                    A_matrix[i,i] += lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))
                # Heat intercept at the second cavity
                if pos_cav2_1 - length_cav/2 <= i*dx <= pos_cav2_1 + length_cav/2:
                    A_matrix[i,i] += lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))
                if pos_cav2_2 - length_cav/2 <= i*dx <= pos_cav2_2 + length_cav/2:
                    A_matrix[i,i] += lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))
                # Heat intercept at the third cavity
                if pos_cav3_1 - length_cav/2 <= i*dx <= pos_cav3_1 + length_cav/2:
                    A_matrix[i,i] += lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))
                if pos_cav3_2 - length_cav/2 <= i*dx <= pos_cav3_2 + length_cav/2:
                    A_matrix[i,i] += lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))


        # Calculation of the T_vector that fullfiles the A*T=B equation
        T_rod_new = np.matmul(np.linalg.inv(A_matrix),B_vector)

        # Neumann Boundary Condition
        T_rod_new[N-1,0] = T_rod_new[N-2,0] - 1/lambda_G10_normal(T_rod[N-2]) * dx * q_rod_N #K


        # Check for convergence
        epsilon = 0
        for i in range(N):
            epsilon += abs(T_rod[i,0] - T_rod_new[i,0])

        epsilon *= 1/sum(T_rod_new)
        print("epsilon: ", epsilon)


        # Creating a new temperature vector with a relaxation factor to avoid overshoot
        for i in range(N):
            T_rod[i,0] = (1-relaxation_factor) * T_rod[i,0] + relaxation_factor * T_rod_new[i,0]


        # Updating the plot every while loop
        line1.set_ydata(T_rod)
        fig.canvas.draw()
        fig.canvas.flush_events()

    ################################################################################

    ################################################################################
    ## Postprocessing ##

    # Actual heat load on the TS
    Q_TS_Rod = 0
    for i in range(N):
        if pos_TS - length_TS/2 <= i*dx <= pos_TS + length_TS/2:
            Q_TS_Rod += (T_rod_new[i,0]-T_TS) * lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))

    # Actual heat loads on the CM
    # Simple linear Fourier heat transport
    Q_CM1_Rod = 0
    Q_CM2_Rod = 0
    Q_CM3_Rod = 0
    for i in range(N):
        if pos_cav1_1 - length_cav/2 <= i*dx <= pos_cav1_1 + length_cav/2:
            Q_CM1_Rod += (T_rod_new[i,0]-T_CM) * lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))
        if pos_cav1_2 - length_cav/2 <= i*dx <= pos_cav1_2 + length_cav/2:
            Q_CM1_Rod += (T_rod_new[i,0]-T_CM) * lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))

        if pos_cav2_1 - length_cav/2 <= i*dx <= pos_cav2_1 + length_cav/2:
            Q_CM2_Rod += (T_rod_new[i,0]-T_CM) * lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))
        if pos_cav2_2 - length_cav/2 <= i*dx <= pos_cav2_2 + length_cav/2:
            Q_CM2_Rod += (T_rod_new[i,0]-T_CM) * lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))

        if pos_cav3_1 - length_cav/2 <= i*dx <= pos_cav3_1 + length_cav/2:
            Q_CM3_Rod += (T_rod_new[i,0]-T_CM) * lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))
        if pos_cav3_2 - length_cav/2 <= i*dx <= pos_cav3_2 + length_cav/2:
            Q_CM3_Rod += (T_rod_new[i,0]-T_CM) * lambda_vector[i] * 2*sp.pi * dx * 1/math.log(diameter_rod/(diameter_rod-thickness_rod))


    print("Q_TS_Rod: ", Q_TS_Rod)
    print("Q_CM1_Rod: ", Q_CM1_Rod)
    print("Q_CM2_Rod: ", Q_CM2_Rod)
    print("Q_CM3_Rod: ", Q_CM3_Rod)
    print("Q_CM_Rod: ", Q_CM1_Rod + Q_CM2_Rod + Q_CM3_Rod)

    ################################################################################

    ################################################################################
    ## Saving of the data ##

    # # Filter for the data
    # N_Filter = 101
    # # N_Filter = N
    # numpy_data = np.zeros((4,N_Filter), dtype = float)
    # for i in range(0, N_Filter):
    #     numpy_data[0,i] = x_lst[int((N-1)/(N_Filter-1) * i)]
    #     numpy_data[1,i] = T_rod[int((N-1)/(N_Filter-1) * i)]
    #     numpy_data[2,i] = T_PTFE[int((N-1)/(N_Filter-1) * i)]
    #     numpy_data[3,i] = T_Cu_outer[int((N-1)/(N_Filter-1) * i)]
    #
    #
    # # Create an xlsx sheet
    # outWorkbook = xlsxwriter.Workbook("N10001_Dirichlet_at2m.xlsx")
    # outSheet = outWorkbook.add_worksheet()
    #
    # # Writing of the headers
    # outSheet.write("A1", "x")
    # outSheet.write("B1", "T_rod")
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
    filename = "data_T_rod_MLI_Test.csv"
    # Data of the file
    fields = ['x', 'T_rod']
    # Restructuring into rows
    rows = np.zeros((N,2), dtype = float)
    for i in range(N):
        rows[i,0] = x_lst[i]
        rows[i,1] = T_rod[i]
    # Writing to csv file
    with open('Data/' + filename, 'w', newline='') as csvfile:
        # Creating a csv writer object
        csvwriter = csv.writer(csvfile)
        # Writing the data
        csvwriter.writerow(fields)
        csvwriter.writerows(rows)



    return Q_CM1_Rod + Q_CM2_Rod + Q_CM3_Rod

################################################################################

################################################################################
## Function calling ##

# Iteration points in x direction
N = 500
# Boundary conditions
T_RT = 293 #K
# Heat intercept with the TS
T_TS = 50 #K
pos_TS = 0.3 #m
length_TS = 0.1 #m
# Heat intercept at the cavities
T_CM = 4.5 #K
pos_cav1_1 = 4.0 #m
pos_cav1_2 = 5.0 #m
pos_cav2_1 = 6.0 #m
pos_cav2_2 = 8.0 #m
pos_cav3_1 = 9.0 #m
pos_cav3_2 = 12.0 #m
length_cav = 0.05 #m

# Emissivities
# Emissivity for the TS made out of Copper
epsilon_TS = 0.05   #0.03; this was the first value but being conservative here doesn't hurt us!
# Emissivity for the rod made out of Al because of the 1-single-layer MLI
# For 1-single-layer MLI: 0.06
# For uncovered G10: 1.0
epsilon_Rod = 0.06

# Relaxation factor for stability of the simulation
relaxation_factor = 0.02


Q_CM_total = tuning_rod_HL_model(N, T_RT, T_TS, pos_TS, length_TS, T_CM, pos_cav1_1, pos_cav1_2, pos_cav2_1, pos_cav2_2, pos_cav3_1, pos_cav3_2, length_cav, epsilon_TS, epsilon_Rod, relaxation_factor)
