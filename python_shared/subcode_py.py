# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 18:29:16 2022

@author: msiodlac
"""

import sys
PATH = r'C:\Users\msiodlac\cernbox\Documents\0_cern_msiodlac\python_shared';
sys.path.insert(0, PATH)

import hepak_wrapper as hp
import numpy as np
import os
from nptdms import TdmsFile
from datetime import datetime as dt

# import sympy
# from sympy.solvers import nsolve
# from sympy.abc import z
from scipy.optimize import fsolve
from scipy.interpolate import interp1d

def tdms_importer(data_folder_path, prefix, file_dates, suffix):
    "Function that imports data from the LabView TDMS files"

    time = np.array([])

    for index, file_date in enumerate(file_dates):
        file_path = os.path.join(data_folder_path, prefix + file_date + suffix)
        print(file_path)
        tdms_file = TdmsFile(file_path)

        # read and format time from the tdms_file
        # raw_time = tdms_file.channel_data(tdms_file.groups()[0], "Time")
        raw_time = tdms_file["Time"]["Time Stamp"]
        # formatted_time = np.array([dt.strptime(date_str, '%H:%M:%S, %d-%b-%Y ') for date_str in raw_time])
        # formatted_time = np.array([dt.strptime(date_str, '%H:%M:%S.%f,%d-%b-%Y') for date_str in raw_time])
        formatted_time = np.array([dt64.astype(dt) for dt64 in raw_time])
        time = np.concatenate([time, formatted_time])


        # construct matrix with temperatures from the tdms_file
        # http://akuederle.com/create-numpy-array-with-for-loop
        volt_data_file = np.empty((0, len(raw_time)))
        channels_data_file = np.empty((0, len(raw_time)))
        for channel in tdms_file["Raw Data"].channels():
            volt_data_file = np.append(volt_data_file, [channel[:]], axis=0)
        for channel in tdms_file["Converted Data"].channels():
            channels_data_file = np.append(channels_data_file, [channel[:]], axis=0)

        # https://simtk.org/svn/openknee/utl/SimVitro/tdms_plottingPF.py
        if index == 0:
            nr_volt_channels = np.shape(volt_data_file)[0]
            nr_channels = np.shape(channels_data_file)[0]
            channel_names = np.array([channel.name for channel in tdms_file["Converted Data"].channels()])[-nr_channels:]
            volt_ch_names = np.array([channel.name for channel in tdms_file["Raw Data"].channels()])[0:nr_volt_channels]
            # channel_names = tdms_file["Converted Data"].channels()[-nr_channels:]
            # volt_ch_names = tdms_file["Raw Data"].channels()[0:nr_volt_channels]

            channels_data = np.empty((np.shape(channels_data_file)[0],0))
            volt_ch_data = np.empty((np.shape(volt_data_file)[0],0))

        channels_data = np.concatenate((channels_data, channels_data_file), axis=1)
        volt_ch_data = np.concatenate((volt_ch_data, volt_data_file), axis=1)

    return time, channels_data, channel_names, volt_ch_data, volt_ch_names




def calcSma(data, smaPeriod):
    "This function is calculating the sliding mean avg across a window of points."
    "    smaPeriod - size of the window (i.e. nr points)"
    "    data - data to be used"

    j = next(i for i, x in enumerate(data) if x is not None)
    our_range = range(len(data))[j + smaPeriod - 1:]
    empty_list = [np.nan] * (j + smaPeriod - 1)
    sub_result = [np.mean(data[i - smaPeriod + 1: i + 1]) for i in our_range]

    return np.array(empty_list + sub_result)

def moving_average(data) :
    "Calculation of the mean average as a mean of 7 values (3 neigbhoring points each and the point itself)."
    ret = np.cumsum(data, dtype=float)
    n = len(ret)
    avg = np.zeros(n)

    avg[0] = ret[0]
    avg[1] = ret[2]/3
    avg[2] = ret[4]/5
    avg[3] = ret[6]/7
    avg[n-3] = (ret[n-1]-ret[n-6])/5
    avg[n-2] = (ret[n-1]-ret[n-4])/3
    avg[n-1] = (ret[n-1]-ret[n-2])
    for i in range(4,n-3):
        avg[i] = (ret[i+3]-ret[i-4])/7
    return avg

def ftt_signal(signal_1D, sampling_rate):
    'Get magnitudes of frequencies for single-sided spectrum plot'

    # signal - signal data
    # sampling_rate - sampling rate in [s]

    N = len(signal_1D)                  # number of points
    f_sampling = 1/sampling_rate        # frequency
    # f_nyquist = f_sampling/2            # Nyquest frequency = max freq to be detected by FT
    magnitudes = abs(np.fft.fft(signal_1D))   # Fourier transform to get magnitudes at each frequency
    bin_vals = np.linspace(0,N-1,N)
    fax_Hz = bin_vals*f_sampling/N      # double-sided frequencies
    N_2 = np.ceil(N/2)                  # single sided spectrum
    frequencies = fax_Hz[0:int(N_2)-1]
    magnitudes = magnitudes[0:int(N_2)-1]

    return frequencies, magnitudes



def V_to_T(T_sensor_name, V_data):

    T_data_corrected = []
    # counter = 0
    counter2 = 0

    for y in V_data:

        if T_sensor_name == "TT855 (2nd stage)":
            # N929
            K1 = 0.436683194689976517
            K2 = 8.373732242267578840
            K3 = -22.9241785244084895
            K4 = 200.8169402340427040
            K5 = -505.582973906304687
            K6 = 664.4906277218833570
            K7 = -172.808795905206352
            # Calculate the resistance
            x = y/1e-5 #Ohm
            T = K1 + K2*1000/x + K3*(1000/x)**2 + K4*(1000/x)**3 + K5*(1000/x)**4 + K6*(1000/x)**5 + K7*(1000/x)**6

        elif T_sensor_name == "TT856" or T_sensor_name == "TT856 (HP out)":
            # N930
            K1 = -0.481209481076803058
            K2 = 26.075956588378176100
            K3 = -183.3263356694951650
            K4 = 902.93139760382473500
            K5 = -2076.907537329941990
            K6 = 2398.7866686843335600
            K7 = -736.7690399140119550
            # Calculate the resistance
            x = y/1e-5 #Ohm
            T = K1 + K2*1000/x + K3*(1000/x)**2 + K4*(1000/x)**3 + K5*(1000/x)**4 + K6*(1000/x)**5 + K7*(1000/x)**6

        elif T_sensor_name == "TT857":
            # N1076
            K1 = -0.348001181751897093
            K2 = 19.024468487245030700
            K3 = -134.5337326875887810
            K4 = 715.84556631185114400
            K5 = -1687.982530876994130
            K6 = 1964.1275345049798500
            K7 = -488.3361798129044470
            # Calculate the resistance
            x = y/1e-5 #Ohm
            T = K1 + K2*1000/x + K3*(1000/x)**2 + K4*(1000/x)**3 + K5*(1000/x)**4 + K6*(1000/x)**5 + K7*(1000/x)**6

        elif T_sensor_name == "TT858" or T_sensor_name == "TT858 (LP in)":
            # N928
            K1 = -1.66467124978225911
            K2 = 52.84870405925903470
            K3 = -376.868188262917101
            K4 = 1596.085047425702210
            K5 = -3316.84907789155841
            K6 = 3364.702600967139010
            K7 = -962.895990639925003
            # Calculate the resistance
            x = y/1e-5 #Ohm
            T = K1 + K2*1000/x + K3*(1000/x)**2 + K4*(1000/x)**3 + K5*(1000/x)**4 + K6*(1000/x)**5 + K7*(1000/x)**6

        elif T_sensor_name == "TT381":
            # Calculate the resistance
            x = y/1e-3 #Ohm
            # standard curve from Labview
            T = 9.35055+39.4229*np.log(x)-35.1221*np.log(x)**2+21.4523*np.log(x)**3-6.1763*np.log(x)**4+0.849805*np.log(x)**5-0.0262748*np.log(x)**6
            # standard curve from Aleks
            # T = 2.8738332e+01 + 2.3135617e+00*x+1.4637766e-03*x**2+-1.5835948e-06*x**3

        elif T_sensor_name == "TT835":
            # Calculate the resistance
            x = y/1e-3 #Ohm
            # standard curve from Labview
            T = 9.35055+39.4229*np.log(x)-35.1221*np.log(x)**2+21.4523*np.log(x)**3-6.1763*np.log(x)**4+0.849805*np.log(x)**5-0.0262748*np.log(x)**6
            # standard curve from Aleks
            # T = 2.8738332e+01 + 2.3135617e+00*x+1.4637766e-03*x**2+-1.5835948e-06*x**3

        elif T_sensor_name == "TT850":
            # Calculate the resistance
            x = y/1e-3 #Ohm
            # standard curve from Labview
            T = 9.35055+39.4229*np.log(x)-35.1221*np.log(x)**2+21.4523*np.log(x)**3-6.1763*np.log(x)**4+0.849805*np.log(x)**5-0.0262748*np.log(x)**6
            # standard curve from Aleks
            # T = 2.8738332e+01 + 2.3135617e+00*x+1.4637766e-03*x**2+-1.5835948e-06*x**3

        elif T_sensor_name == "TT851":
            # Calculate the resistance
            x = y/1e-3 #Ohm
            # standard curve from Labview
            T = 9.35055+39.4229*np.log(x)-35.1221*np.log(x)**2+21.4523*np.log(x)**3-6.1763*np.log(x)**4+0.849805*np.log(x)**5-0.0262748*np.log(x)**6
            # standard curve from Aleks
            # T = 2.8738332e+01 + 2.3135617e+00*x+1.4637766e-03*x**2+-1.5835948e-06*x**3

        elif T_sensor_name == "TT851":
            # Calculate the resistance
            x = y/1e-3 #Ohm
            # standard curve from Labview
            T = 9.35055+39.4229*np.log(x)-35.1221*np.log(x)**2+21.4523*np.log(x)**3-6.1763*np.log(x)**4+0.849805*np.log(x)**5-0.0262748*np.log(x)**6
            # standard curve from Aleks
            # T = 2.8738332e+01 + 2.3135617e+00*x+1.4637766e-03*x**2+-1.5835948e-06*x**3

        elif T_sensor_name == "TT852" or T_sensor_name == "TT852 (HP in)":
            # Calculate the resistance
            x = y/1e-3 #Ohm
            # standard curve from Labview
            T = 9.35055+39.4229*np.log(x)-35.1221*np.log(x)**2+21.4523*np.log(x)**3-6.1763*np.log(x)**4+0.849805*np.log(x)**5-0.0262748*np.log(x)**6
            # standard curve from Aleks
            # T = 2.8738332e+01 + 2.3135617e+00*x+1.4637766e-03*x**2+-1.5835948e-06*x**3

        elif T_sensor_name == "TT853" or T_sensor_name == "TT853 (LP out)":
            # Calculate the resistance
            x = y/1e-3 #Ohm
            # standard curve from Labview
            T = 9.35055+39.4229*np.log(x)-35.1221*np.log(x)**2+21.4523*np.log(x)**3-6.1763*np.log(x)**4+0.849805*np.log(x)**5-0.0262748*np.log(x)**6
            # standard curve from Aleks
            # T = 2.8738332e+01 + 2.3135617e+00*x+1.4637766e-03*x**2+-1.5835948e-06*x**3

        elif T_sensor_name == "TT840":
            # Calculate the resistance
            x = y/1e-3 #Ohm
            # standard curve from Labview
            T = 9.35055+39.4229*np.log(x)-35.1221*np.log(x)**2+21.4523*np.log(x)**3-6.1763*np.log(x)**4+0.849805*np.log(x)**5-0.0262748*np.log(x)**6
            # standard curve from Aleks
            # T = 2.8738332e+01 + 2.3135617e+00*x+1.4637766e-03*x**2+-1.5835948e-06*x**3



        # elif T_sensor_name == "TT800":
        #     # same as TT810 because it stays at warm anyways
        #     if x < 29.5699:
        #         T = 1.1697845e+01 + 3.2348376e+01*np.log(x)+-3.0710738e+01*np.log(x)**2+2.3909298e+01*np.log(x)**3+-9.3619015e+00*np.log(x)**4+1.8733248e+00*np.log(x)**5+-1.3349663e-01*np.log(x)**6
        #     elif x >= 29.5699 and x < 70.73:
        #         T = -9.1629179e+01 + -1.5571503e+03*np.log(x)+2.1374694e+03*np.log(x)**2+-1.1220743e+03*np.log(x)**3+2.9210700e+02*np.log(x)**4+-3.7892951e+01*np.log(x)**5+1.9812156e+00*np.log(x)**6
        #     elif x >= 70.73:
        #         # T = 2.8738332e+01 + 2.3135617e+00*x+1.4637766e-03*x**2+-1.5835948e-06*x**3  # standard Pt100 curve: new, mistake corrected
        #         T = 2.8881500e+01 + 2.3135617e+00*x+1.4637766e-03*x**2-1.5835948e-06*x**3      # adjusted to pass through 273.16 K point
        #     else:
        #         print("Invalid sensor resistance")
        #         return



        else:
            counter2 = counter2 + 1
            T = 0
            if counter2 == 1:
                print(f"V_to_T function message: Invalid T sensor name for calibration. Sensor name: {T_sensor_name}")


        # print(x)
        # print(offset)
        T_data_corrected = np.append(T_data_corrected, T)

    return T_data_corrected


def V_to_P(P_sensor_name, V_data):

    P_data_corrected = []
    # counter = 0
    counter2 = 0

    for x in V_data:

        if P_sensor_name == "PT850 (bara)":
            # Range corresponds to 0 and 25 bar
            # Signal is between 4 and 20 mA or 1 and 5 V
            P = (-6.25 + 6.25*x)*1e5 #Pa

        elif P_sensor_name == "PDT801 (mbar)":
            # Range corresponds to 0 and 50 mbar
            # Signal is between 4 and 20 mA or 1 and 5 V
            # This sensor also has a positive offset of 0.4 mbar that we want to substract here
            P = (-12.5 + 12.5*x)*1e2 - 0.4e2 #Pa

        elif P_sensor_name == "PDT802 (mbar)":
            # Range corresponds to 0 and 620 mbar
            # Signal is between 4 and 20 mA or 1 and 5 V
            P = (-155. + 155.*x)*1e2 #Pa

        else:
            counter2 = counter2 + 1
            P = 0
            if counter2 == 1:
                print(f"V_to_P function message: Invalid P sensor name. Sensor name: {P_sensor_name}")

        # print(x)
        # print(offset)
        P_data_corrected = np.append(P_data_corrected, P)

    return P_data_corrected

def V_to_RPM(sensor_name, V_data):
    "Caclulate the rpm of the Cryofan from the voltage signal: Böhmwind CryoFan"
    RPM_data_corrected = []
    for x in V_data:
        # Signal between 0 and 4 V
        # Corresponds to 0 and 22000 rpm
        RPM = 5500*x

        RPM_data_corrected = np.append(RPM_data_corrected, RPM)

    return RPM_data_corrected



def MFT_calc(T_TT851, T_TT852, Q_EH732, correction_situ_dT, p_PT850):
    "MFT calculation as a function of the temperature sensors TT851 and TT852 and the voltage drop over EH732"

    # TT851 in K
    T_in = T_TT851 #K
    # TT852 in K
    T_out = T_TT852 #K
    # EH732 in W
    Q_EH = Q_EH732 #V
    # nominal pressure from PT850
    p_nominal = p_PT850 #Pa

    # Corrective factors for the offset of the two sensors to each other
    # Corrective situ factor measured in situ during the measurement for no heat applied
    # This is the corrective factor for parasitic heat loads in the setup and is only valid for this specific measuring point
    # correction_situ_dT

    # Heat capactiy calculated with HEPAK
    # Here the two absolute values of the sensors are used even when they are very inaccurate
    cp_in = hp.HeCalc(14, 0, 1, p_nominal, 2, T_in, 1) #J/(kgK)
    cp_out = hp.HeCalc(14, 0, 1, p_nominal, 2, T_out, 1) #J/(kgK)

    # Corrected temperature difference
    dT = T_out - T_in - correction_situ_dT #K
    return Q_EH/dT/(0.5*(cp_in + cp_out)) #kg/s

def Q_dT(T_in, T_out, mass_flow, p_nominal):
# Heat load as a function of the inlet temperature, outlet temperature, mass flow, and pressure

    # Heat capactiy calculated with HEPAK
    # Here the two absolute values of the sensors are used even when they are very inaccurate
    cp_in = hp.HeCalc(14, 0, 1, p_nominal, 2, T_in, 1) #J/(kgK)
    cp_out = hp.HeCalc(14, 0, 1, p_nominal, 2, T_out, 1) #J/(kgK)

    # Heat load
    Q_HL = mass_flow * (0.5*(cp_in + cp_out)) * (T_out-T_in) #W

    return Q_HL

def Uncertainty_calc_sensor(sensor_name, avg_value, std_dev):

    counter2 = 0
    # Coverage factor describes the probalitiy in a normal distribution
    # First every uncertainty needs to be brought to the same coverage factor and then the combined uncertainty can be multiplied with the chosen coverage factor
    coverage_factor = 2

    if sensor_name == "TT855 (2nd stage)":
        # N929
        K1 = 0.436683194689976517
        K2 = 8.373732242267578840
        K3 = -22.9241785244084895
        K4 = 200.8169402340427040
        K5 = -505.582973906304687
        K6 = 664.4906277218833570
        K7 = -172.808795905206352

        # Calculate back the resistance corresponding to the averaged temperature
        def f(z):
            return -avg_value + -avg_value + K1 + K2*1000*z**(-1) + K3*1000**2*z**(-2) + K4*1000**3*z**(-3) + K5*1000**4*z**(-4) + K6*1000**5*z**(-5) + K7*1000**6*z**(-6)
        solution = fsolve(f,1900.0)
        R = solution[0] #Ohm
        # Standard current for TVO
        I = 10e-6 #A
        U = R*I #V
        # Sensistivity coefficients
        dT_dU = -K2*(1000*I)*U**(-2) - 2*K3*(1000*I)**2*U**(-3) - 3*K4*(1000*I)**3*U**(-4) - 4*K5*(1000*I)**4*U**(-5) - 5*K6*(1000*I)**5*U**(-6) - 6*K7*(1000*I)**6*U**(-7)
        dT_dI =  K2*1000/U + 2*K3*(1000/U)**2*I + 3*K4*(1000/U)**3*I**2 + 4*K5*(1000/U)**4*I**3 + 5*K6*(1000/U)**5*I**4 + 6*K7*(1000/U)**6*I**5

        # Uncertainties according to GUM and the respective datasheets
        u_A_noise = std_dev #K
        # Calibration uncertainty for now, data for CX SD 1050 assumed (discussed with TK)
        # Because I could not find a TMi TVO A2 calibration uncertainty
        T_points_acc = [1.4, 4.2, 10, 20, 30, 50, 100, 300, 400]
        err_acc_points = [4e-3, 4e-3, 4e-3, 8e-3, 9e-3, 12e-3, 17e-3, 46e-3, 74e-3]
        interp = interp1d(T_points_acc, err_acc_points)
        u_B_cal = interp(avg_value) #K
        # Fit uncertainty from the corresponding data sheet
        T_points_fit = [1.490, 4.233, 10, 12, 15, 18, 20, 25, 238]
        err_fit_points = [4.69e-3, 3.971e-3, 0.197e-3, 3.435e-3, 9.084e-3, 5.4e-3, 1.543e-3, 7.432e-3, 17.626e-3]
        interp_fit = interp1d(T_points_fit, err_fit_points)
        u_B_fit = interp_fit(avg_value)

        u_B_volt = (80e-6*U + 0.2*57e-6)/3 #V #coverage factor = 3 see documentation
        u_B_current = (5.2e-8)/3**0.5 #A # rectangular distribution assumed

        # Squared combined uncertainty
        u_c_sqrd = u_A_noise**2 + u_B_cal**2 + u_B_fit**2 + (dT_dU*u_B_volt)**2 + (dT_dI*u_B_current)**2
        # Final expanded uncertainty
        u_e = coverage_factor*u_c_sqrd**0.5


    elif sensor_name == "TT856" or sensor_name == "TT856 (HP out)":
        # N930
        K1 = -0.481209481076803058
        K2 = 26.075956588378176100
        K3 = -183.3263356694951650
        K4 = 902.93139760382473500
        K5 = -2076.907537329941990
        K6 = 2398.7866686843335600
        K7 = -736.7690399140119550

        # Calculate back the resistance corresponding to the averaged temperature
        # solution = sympy.solve(-avg_value + K1 + K2*1000*z**(-1) + K3*1000**2*z**(-2) + K4*1000**3*z**(-3) + K5*1000**4*z**(-4) + K6*1000**5*z**(-5) + K7*1000**6*z**(-6), "z")
        def f(z):
            return -avg_value + K1 + K2*1000*z**(-1) + K3*1000**2*z**(-2) + K4*1000**3*z**(-3) + K5*1000**4*z**(-4) + K6*1000**5*z**(-5) + K7*1000**6*z**(-6)
        solution = fsolve(f,1900.0)
        R = solution[0] #Ohm
        # Standard current for TVO
        I = 10e-6 #A
        U = R*I #V
        # Sensistivity coefficients
        dT_dU = -K2*(1000*I)*U**(-2) - 2*K3*(1000*I)**2*U**(-3) - 3*K4*(1000*I)**3*U**(-4) - 4*K5*(1000*I)**4*U**(-5) - 5*K6*(1000*I)**5*U**(-6) - 6*K7*(1000*I)**6*U**(-7)
        dT_dI =  K2*1000/U + 2*K3*(1000/U)**2*I + 3*K4*(1000/U)**3*I**2 + 4*K5*(1000/U)**4*I**3 + 5*K6*(1000/U)**5*I**4 + 6*K7*(1000/U)**6*I**5

        # Uncertainties according to GUM and the respective datasheets
        u_A_noise = std_dev #K
        # Calibration uncertainty for now, data for CX SD 1050 assumed (discussed with TK)
        # Because I could not find a TMi TVO A2 calibration uncertainty
        T_points_acc = [1.4, 4.2, 10, 20, 30, 50, 100, 300, 400]
        err_acc_points = [4e-3, 4e-3, 4e-3, 8e-3, 9e-3, 12e-3, 17e-3, 46e-3, 74e-3]
        interp = interp1d(T_points_acc, err_acc_points)
        u_B_cal = interp(avg_value) #K
        # Fit uncertainty from the corresponding data sheet
        T_points_fit = [1.490, 4.233, 10, 12, 15, 18, 20, 25, 238]
        err_fit_points = [4.779e-3, 5.557e-3, 2.534e-3, 1.346e-3, 4.624e-3, 8.708e-3, 1.936e-3, 0.142e-3, 11.939e-3]
        interp_fit = interp1d(T_points_fit, err_fit_points)
        u_B_fit = interp_fit(avg_value)

        u_B_volt = (80e-6*U + 0.2*57e-6)/3 #V #coverage factor = 3 see documentation
        u_B_current = (5.2e-8)/3**0.5 #A # rectangular distribution assumed

        # Squared combined uncertainty
        u_c_sqrd = u_A_noise**2 + u_B_cal**2 + u_B_fit**2 + (dT_dU*u_B_volt)**2 + (dT_dI*u_B_current)**2
        # Final expanded uncertainty
        u_e = coverage_factor*u_c_sqrd**0.5

    elif sensor_name == "TT857":
        # N1076
        K1 = -0.348001181751897093
        K2 = 19.024468487245030700
        K3 = -134.5337326875887810
        K4 = 715.84556631185114400
        K5 = -1687.982530876994130
        K6 = 1964.1275345049798500
        K7 = -488.3361798129044470

        # Calculate back the resistance corresponding to the averaged temperature
        # solution = sympy.solve(-avg_value + K1 + K2*1000*z**(-1) + K3*1000**2*z**(-2) + K4*1000**3*z**(-3) + K5*1000**4*z**(-4) + K6*1000**5*z**(-5) + K7*1000**6*z**(-6), "z")
        def f(z):
            return -avg_value + K1 + K2*1000*z**(-1) + K3*1000**2*z**(-2) + K4*1000**3*z**(-3) + K5*1000**4*z**(-4) + K6*1000**5*z**(-5) + K7*1000**6*z**(-6)
        solution = fsolve(f,1900.0)
        R = solution[0] #Ohm
        # Standard current for TVO
        I = 10e-6 #A
        U = R*I #V
        # Sensistivity coefficients
        dT_dU = -K2*(1000*I)*U**(-2) - 2*K3*(1000*I)**2*U**(-3) - 3*K4*(1000*I)**3*U**(-4) - 4*K5*(1000*I)**4*U**(-5) - 5*K6*(1000*I)**5*U**(-6) - 6*K7*(1000*I)**6*U**(-7)
        dT_dI =  K2*1000/U + 2*K3*(1000/U)**2*I + 3*K4*(1000/U)**3*I**2 + 4*K5*(1000/U)**4*I**3 + 5*K6*(1000/U)**5*I**4 + 6*K7*(1000/U)**6*I**5

        # Uncertainties according to GUM and the respective datasheets
        u_A_noise = std_dev #K
        # Calibration uncertainty for now, data for CX SD 1050 assumed (discussed with TK)
        # Because I could not find a TMi TVO A2 calibration uncertainty
        T_points_acc = [1.4, 4.2, 10, 20, 30, 50, 100, 300, 400]
        err_acc_points = [0.1e-3, 4e-3, 4e-3, 8e-3, 9e-3, 12e-3, 17e-3, 46e-3, 74e-3]
        interp = interp1d(T_points_acc, err_acc_points)
        u_B_cal = interp(avg_value) #K
        # Fit uncertainty from the corresponding data sheet
        T_points_fit = [1.490, 4.233, 10, 12, 15, 18, 20, 25, 238]
        err_fit_points = [0.959e-3, 3.427e-3, 8.11e-3, 2.466e-3, 4.04e-3, 1.98e-3, 2.235e-3, 4.392e-3, 16.003e-3]
        interp_fit = interp1d(T_points_fit, err_fit_points)
        u_B_fit = interp_fit(avg_value)

        u_B_volt = (80e-6*U + 0.2*57e-6)/3 #V #coverage factor = 3 see documentation
        u_B_current = (5.2e-8)/3**0.5 #A # rectangular distribution assumed

        # Squared combined uncertainty
        u_c_sqrd = u_A_noise**2 + u_B_cal**2 + u_B_fit**2 + (dT_dU*u_B_volt)**2 + (dT_dI*u_B_current)**2
        # Final expanded uncertainty
        u_e = coverage_factor*u_c_sqrd**0.5

    elif sensor_name == "TT858" or sensor_name == "TT858 (LP in)":
        # N928
        K1 = -1.66467124978225911
        K2 = 52.84870405925903470
        K3 = -376.868188262917101
        K4 = 1596.085047425702210
        K5 = -3316.84907789155841
        K6 = 3364.702600967139010
        K7 = -962.895990639925003

        # Calculate back the resistance corresponding to the averaged temperature
        # solution = sympy.solve(-avg_value + K1 + K2*1000*z**(-1) + K3*1000**2*z**(-2) + K4*1000**3*z**(-3) + K5*1000**4*z**(-4) + K6*1000**5*z**(-5) + K7*1000**6*z**(-6), "z")
        def f(z):
            return -avg_value + K1 + K2*1000*z**(-1) + K3*1000**2*z**(-2) + K4*1000**3*z**(-3) + K5*1000**4*z**(-4) + K6*1000**5*z**(-5) + K7*1000**6*z**(-6)
        solution = fsolve(f,1900.0)
        R = solution[0] #Ohm
        # Standard current for TVO
        I = 10e-6 #A
        U = R*I #V
        # Sensistivity coefficients
        dT_dU = -K2*(1000*I)*U**(-2) - 2*K3*(1000*I)**2*U**(-3) - 3*K4*(1000*I)**3*U**(-4) - 4*K5*(1000*I)**4*U**(-5) - 5*K6*(1000*I)**5*U**(-6) - 6*K7*(1000*I)**6*U**(-7)
        dT_dI =  K2*1000/U + 2*K3*(1000/U)**2*I + 3*K4*(1000/U)**3*I**2 + 4*K5*(1000/U)**4*I**3 + 5*K6*(1000/U)**5*I**4 + 6*K7*(1000/U)**6*I**5

        # Uncertainties according to GUM and the respective datasheets
        u_A_noise = std_dev #K
        # Calibration uncertainty for now, data for CX SD 1050 assumed (discussed with TK)
        # Because I could not find a TMi TVO A2 calibration uncertainty
        T_points_acc = [1.4, 4.2, 10, 20, 30, 50, 100, 300, 400]
        err_acc_points = [4e-3, 4e-3, 4e-3, 8e-3, 9e-3, 12e-3, 17e-3, 46e-3, 74e-3]
        interp = interp1d(T_points_acc, err_acc_points)
        u_B_cal = interp(avg_value) #K
        # Fit uncertainty from the corresponding data sheet
        T_points_fit = [1.490, 4.233, 10, 12, 15, 18, 20, 25, 238]
        err_fit_points = [7.018e-3, 1.289e-3, 1.177e-3, 2.354e-3, 3.498e-3, 1.876e-3, 0.581e-3, 7.277e-3, 12.451e-3]
        interp_fit = interp1d(T_points_fit, err_fit_points)
        u_B_fit = interp_fit(avg_value)

        u_B_volt = (80e-6*U + 0.2*57e-6)/3 #V #coverage factor = 3 see documentation
        u_B_current = (5.2e-8)/3**0.5 #A # rectangular distribution assumed

        # Squared combined uncertainty
        u_c_sqrd = u_A_noise**2 + u_B_cal**2 + u_B_fit**2 + (dT_dU*u_B_volt)**2 + (dT_dI*u_B_current)**2
        # Final expanded uncertainty
        u_e = coverage_factor*u_c_sqrd**0.5

    elif sensor_name == "TT381":
        # standard curve from Labview
        # T = 9.35055+39.4229*np.log(x)-35.1221*np.log(x)**2+21.4523*np.log(x)**3-6.1763*np.log(x)**4+0.849805*np.log(x)**5-0.0262748*np.log(x)**6
        # standard curve from Aleks
        # T = 2.8738332e+01 + 2.3135617e+00*x+1.4637766e-03*x**2+-1.5835948e-06*x**3

        # Calculate back the resistance corresponding to the averaged temperature
        def f(z):
            return -avg_value + 2.8738332e+01 + 2.3135617*z+1.4637766e-03*z**2+-1.5835948e-06*z**3
        solution = fsolve(f,41.0)
        R = solution[0] #Ohm
        # Standard current for Pt100
        I = 1e-3 #A
        U = R*I #V
        # Sensistivity coefficients
        dT_dU = 2.3135617/I + 2*1.4637766e-03*I**(-2)*U - 3*1.5835948e-06*I**(-3)*U**2
        dT_dI = -2.3135617*U - 2*1.4637766e-03*I**(-3)*U**2 - 3*1.5835948e-06*I**(-4)*U**3

        # Uncertainties according to GUM and the respective datasheets
        u_A_noise = std_dev #K
        u_B_fit = 1/3 *(0.3 + 0.005*abs(avg_value-273.15)) #K
        u_B_volt = (80e-6*U + 0.2*57e-6)/3 #V #coverage factor = 3 see documentation
        u_B_current = (5.2e-8)/3**0.5 #A # rectangular distribution assumed

        # Squared combined uncertainty
        u_c_sqrd = u_A_noise**2 + u_B_fit**2 + (dT_dU*u_B_volt)**2 + (dT_dI*u_B_current)**2
        # Final expanded uncertainty
        u_e = coverage_factor*u_c_sqrd**0.5

    elif sensor_name == "TT835" or sensor_name == "TT850" or sensor_name == "TT851" or sensor_name == "TT852" or sensor_name == "TT852 (HP in)" or sensor_name == "TT853" or sensor_name == "TT853 (LP out)" or sensor_name == "TT840":
        # standard curve from Labview
        # T = 9.35055+39.4229*np.log(x)-35.1221*np.log(x)**2+21.4523*np.log(x)**3-6.1763*np.log(x)**4+0.849805*np.log(x)**5-0.0262748*np.log(x)**6
        # standard curve from Aleks
        # T = 2.8738332e+01 + 2.3135617e+00*x+1.4637766e-03*x**2+-1.5835948e-06*x**3

        # Calculate back the resistance corresponding to the averaged temperature
        def f(z):
            return -avg_value + 2.8738332e+01 + 2.3135617*z+1.4637766e-03*z**2+-1.5835948e-06*z**3
        solution = fsolve(f,12.0)
        R = solution[0] #Ohm
        # Curve of IEC 60751 - I got a a bit different values with this correlation but I do not have time to assess why!
        # R = 100*(1 + 3.9083e-3*(avg_value-273.15) - 5.775e-7*(avg_value-273.15)**2 - 4.183e-12*(avg_value-273.15-100)*(avg_value-273.15)**3)

        # Standard current for Pt100
        I = 1e-3 #A
        U = R*I #V
        # Sensistivity coefficients
        dT_dU = 2.3135617/I + 2*1.4637766e-03*I**(-2)*U - 3*1.5835948e-06*I**(-3)*U**2
        dT_dI = -2.3135617*U - 2*1.4637766e-03*I**(-3)*U**2 - 3*1.5835948e-06*I**(-4)*U**3

        # Uncertainties according to GUM and the respective datasheets
        u_A_noise = std_dev #K
        u_B_fit = 1/3 *(0.3 + 0.005*abs(avg_value-273.15)) #K
        u_B_volt = (80e-6*U + 0.2*57e-6)/3 #V #coverage factor = 3 see documentation
        u_B_current = (5.2e-8)/3**0.5 #A # rectangular distribution assumed

        # Squared combined uncertainty
        u_c_sqrd = u_A_noise**2 + u_B_fit**2 + (dT_dU*u_B_volt)**2 + (dT_dI*u_B_current)**2
        # Final expanded uncertainty
        u_e = coverage_factor*u_c_sqrd**0.5

    elif sensor_name == "PT850" or sensor_name == "PT850 (bara)":
        # Signal as current between 4 and 20 mA or 1 and 5 V
        # This signal is transformed into voltage via high precision resistors of 250 Ohm
        # Back calculation of the voltage signal
        # Range corresponds to 0 and 25 bar
        U = (avg_value*1e-5 + 6.25)/6.25 #V

        # Sensitivity coefficient
        dP_dU = 6.25*1e5

        # Uncertainties according to GUM and the respective datasheets
        u_A_noise = std_dev #Pa
        u_B_acc = (0.25e-2*25e5)/3**0.5 #Pa # rectangular distribution assumed
        u_B_volt = (85e-6*U + 5*96e-6)/3  #V #coverage factor = 3 see documentation

        # Squared combined uncertainty
        u_c_sqrd = u_A_noise**2 + u_B_acc**2 + (dP_dU*u_B_volt)**2
        # Final expanded uncertainty
        u_e = coverage_factor*u_c_sqrd**0.5


    elif sensor_name == "PDT801" or sensor_name == "PDT801 (mbar)":
        # Signal as current between 4 and 20 mA or 1 and 5 V
        # This signal is transformed into voltage via high precision resistors of 250 Ohm
        # Back calculation of the voltage signal
        # Range corresponds to 0 and 50 mbar
        U = (avg_value*1e-2 + 12.5)/12.5 #V
        # Sensitivity coefficient
        dP_dU = 12.5*1e2

        # Uncertainties according to GUM and the respective datasheets
        u_A_noise = std_dev #Pa
        u_B_volt = 85e-6*U + 5*96e-6  #V #coverage factor = 3 see documentation

        # Datasheet Rosemount 3051C page 53
        # This absolute accuracy corresponds to 3sigma!!
        # Range 2 (tag on the sensor) -> defines Upper Range Limit (URL)
        URL = 621e2 #Pa
        u_B_acc = (0.14/100*50e2 + 0.2/100*URL)/3 #Pa

        # Squared combined uncertainty
        u_c_sqrd = u_A_noise**2 + u_B_acc**2 + (dP_dU*u_B_volt)**2
        # Final expanded uncertainty
        u_e = coverage_factor*u_c_sqrd**0.5

    elif sensor_name == "PDT802" or sensor_name == "PDT802 (mbar)":
        # Signal as current between 4 and 20 mA or 1 and 5 V
        # This signal is transformed into voltage via high precision resistors of 250 Ohm
        # Back calculation of the voltage signal
        # Range corresponds to 0 and 620 mbar
        U = (avg_value*1e-2 + 155.)/155. #V
        # Sensitivity coefficient
        dP_dU = 155.*1e2

        # Uncertainties according to GUM and the respective datasheets
        u_A_noise = std_dev #Pa
        u_B_volt = 85e-6*U + 5*96e-6  #V #coverage factor = 3 see documentation

        # Datasheet Rosemount 3051C page 53
        # This absolute accuracy corresponds to 3sigma!!
        # Range 2 (tag on the sensor) -> defines Upper Range Limit (URL)
        URL = 621e2 #Pa
        u_B_acc = (0.14/100*620e2 + 0.2/100*URL)/3 #Pa

        # Squared combined uncertainty
        u_c_sqrd = u_A_noise**2 + u_B_acc**2 + (dP_dU*u_B_volt)**2
        # Final expanded uncertainty
        u_e = coverage_factor*u_c_sqrd**0.5

    elif sensor_name == "CryoFan (rpm)":
        # The CryoFan is perfect!
        # Joke.
        # I just do not care
        u_e = 0

    else:
        counter2 = counter2 + 1
        u_e = 0
        if counter2 == 1:
            print(f"Uncertainty_calc_sensor function message: Invalid sensor name for calculation. Sensor name: {sensor_name}")


    # print(x)
    # print(offset)
    # T_data_corrected = np.append(T_data_corrected, T)

    return u_e

def Uncertainty_calc_EH(sensor_name, voltage, current):
    "Uncertainty calculation according to GUM for the electric heaters in W"

    # Sensitivity coefficients
    dQ_dU = current #W/V
    dQ_dI = voltage #W/A

    # Uncertainties according to GUM and the respective datasheets
    # 2nd stage heater
    if sensor_name == "EH731":
        if voltage < 10:
            u_B_volt = 30e-6*voltage + 5e-6*10 #V
        elif voltage < 100:
            u_B_volt = 45e-6*voltage + 9e-6*100 #V
        else:
            print("Error: Too much voltage for EH731!")
    # MFT heater
    elif sensor_name == "EH732":
        u_B_volt = voltage*0.05e-2 + 5e-3 #V
    # Cavity heater
    elif sensor_name == "EH835":
        u_B_volt = voltage*0.05e-2 + 5e-3 #V
    else:
        print(f"Uncertainty_calc_EH function message: Invalid T sensor name for calibration. Sensor name: {sensor_name}")
    u_B_current = current*0.1e-2 + 2e-3 #A

    # Squared combined uncertainty
    u_c_sqrd = (dQ_dU*u_B_volt)**2 + (dQ_dI*u_B_current)**2 #W*W
    # Final combined uncertainty
    u_c = 2*u_c_sqrd**0.5 #W

    return u_c

def Uncertainty_calc_massflow(T_in, T_out, Q_EH, correction_situ_dT, p_nominal, std_dev_Tin, std_dev_Tout):
    "Uncertainty calculation according to GUM for the mass flow in kg/s"
    # Assuming that the uncertainty of the pressure does not influence the uncertainty of the mass flow because the dependency of cp(p) is small

    # Heat capactiy calculated with HEPAK
    # Here the two absolute values of the sensors are used even when they are very inaccurate
    cp_in = hp.HeCalc(14, 0, 1, p_nominal, 2, T_in, 1) #J/(kgK)
    cp_out = hp.HeCalc(14, 0, 1, p_nominal, 2, T_out, 1) #J/(kgK)

    # Estimating the derivative of the cp at the respective temperatures dy/dx = (y_2-y_1)/(x_2-x_1)
    d_cp_in = (hp.HeCalc(14, 0, 1, p_nominal, 2, T_in+0.1, 1)-hp.HeCalc(14, 0, 1, p_nominal, 2, T_in-0.1, 1))/0.2 #J/(kgK)/K
    d_cp_out = (hp.HeCalc(14, 0, 1, p_nominal, 2, T_out+0.1, 1)-hp.HeCalc(14, 0, 1, p_nominal, 2, T_out-0.1, 1))/0.2 #J/(kgK)/K

    # Corrected temperature difference
    dT = T_out - T_in - correction_situ_dT #K
    #Q_EH/dT/(0.5*(cp_in + cp_out)) #kg/s

    # Sensitivity coefficients
    dM_dQ = 1/dT/(0.5*(cp_in + cp_out)) #kg/s/W
    dM_dTin = -2*Q_EH * (d_cp_in*dT - (cp_in+cp_out))/((cp_in+cp_out)*dT)**2 #kg/s/K
    dM_dTout = -2*Q_EH * (d_cp_out*dT + (cp_in+cp_out))/((cp_in+cp_out)*dT)**2 #kg/s/K

    # Uncertainties according to GUM and the respective datasheets
    u_c_EH = Uncertainty_calc_EH("EH732", 15, 152.2e-3) #W
    # We care only about the differential temperature thus all the Type B uncertainties cancel each other out and only the noise is relevant
    u_A_noise_Tin = std_dev_Tin #K
    u_A_noise_Tout = std_dev_Tout #K

    # Squared combined uncertainty
    u_c_sqrd = (dM_dQ*u_c_EH)**2 + (dM_dTin*u_A_noise_Tin)**2 + (dM_dTout*u_A_noise_Tout)**2 #(kg/s)**2
    # Final expanded uncertainty
    u_e = 2*u_c_sqrd**0.5 #kg/s

    return u_e #kg/s

def Uncertainty_calc_QCryofan(T_in, T_out, mass_flow, p_nominal, std_dev_Tin, std_dev_Tout, uncertainty_mass_flow):
    "Uncertainty calculation according to GUM for the Cryofan Heat load in W"
    # Assuming that the uncertainty of the pressure does not influence the uncertainty of the mass flow because the dependency of cp(p) is small

    # Heat capactiy calculated with HEPAK
    # Here the two absolute values of the sensors are used even when they are very inaccurate
    cp_in = hp.HeCalc(14, 0, 1, p_nominal, 2, T_in, 1) #J/(kgK)
    cp_out = hp.HeCalc(14, 0, 1, p_nominal, 2, T_out, 1) #J/(kgK)

    # Estimating the derivative of the cp at the respective temperatures dy/dx = (y_2-y_1)/(x_2-x_1)
    d_cp_in = (hp.HeCalc(14, 0, 1, p_nominal, 2, T_in+0.1, 1)-hp.HeCalc(14, 0, 1, p_nominal, 2, T_in-0.1, 1))/0.2 #J/(kgK)/K
    d_cp_out = (hp.HeCalc(14, 0, 1, p_nominal, 2, T_out+0.1, 1)-hp.HeCalc(14, 0, 1, p_nominal, 2, T_out-0.1, 1))/0.2 #J/(kgK)/K

    # Corrected temperature difference
    dT = T_out - T_in #K
    #Q_EH/dT/(0.5*(cp_in + cp_out)) #kg/s

    # Sensitivity coefficients
    dQ_dM = dT*(0.5*(cp_in + cp_out)) #Ws/kg
    dQ_dTin = 0.5*mass_flow * (d_cp_in*dT - (cp_in+cp_out)) #W/K
    dQ_dTout = 0.5*mass_flow * (d_cp_out*dT + (cp_in+cp_out)) #W/K

    # Uncertainties according to GUM and the respective datasheets
    u_c_mass_flow = uncertainty_mass_flow/2 #kg/s
    # We care only about the differential temperature thus all the Type B uncertainties cancel each other out and only the noise is relevant
    u_A_noise_Tin = std_dev_Tin #K
    u_A_noise_Tout = std_dev_Tout #K

    # Squared combined uncertainty
    u_c_sqrd = (dQ_dM*u_c_mass_flow)**2 + (dQ_dTin*u_A_noise_Tin)**2 + (dQ_dTout*u_A_noise_Tout)**2 #W**2
    # Final expanded uncertainty
    u_e = 2*u_c_sqrd**0.5 #kg/s

    return u_e #kg/s

