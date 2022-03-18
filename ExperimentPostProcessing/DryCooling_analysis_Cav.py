# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 18:52:00 2022

@author: msiodlac
"""

import numpy as np
import os
import sys
PATH = r'C:\Users\msiodlac\cernbox\Documents\0_cern_msiodlac\python_shared';
sys.path.insert(0, PATH)

import hepak_wrapper as hp

from shutil import rmtree
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from nptdms import TdmsFile
from datetime import datetime as dt
from scipy.optimize import curve_fit
from subcode_py import *
from scipy.ndimage.interpolation import shift
from prettytable import PrettyTable, ALL
import datetime

#
plt.close('all')

# import data
data_folder_path = 'C:\\Users\\msiodlac\\cernbox\\Documents\\4_Experiment\\Data\\FirstCooldown'
prefix = "RemoteCoolingCavity_MarcTest_"
file_dates_MFT = ['2022_02_17_08_52_05', '2022_02_24_13_05_51']
file_dates_Cav = ['2022_02_21_13_13_47', '2022_02_22_13_13_48', '2022_02_23_13_13_49']
suffix = ".tdms"


# results_folder_path = r"Results"

# for future deletion of uneven regions
start_Cav_cst_time = [69086, 74200, 76836, 77802, 84233, 88099, 92026, 94564, 99000, 101194, 104199, 107398, 110196, 112479,
                      166991, 170801, 176311, 180862, 188680, 191115, 193711, 195840, 199722, 202002, 204495]
end_Cav_cst_time = [71215, 75384, 77412, 78172, 85583, 89162, 93295, 95397, 100086, 102297, 105230, 108353, 111201, 113203,
                    168788, 171850, 177971, 182336, 189472, 192000, 194407, 196588, 200281, 202891, 205504]
start_Cav_cst_time = [x-300 for x in end_Cav_cst_time]

# import data
time, channels_data, channel_names, volt_ch_data, volt_ch_names = tdms_importer(data_folder_path, prefix, file_dates_Cav, suffix)


## Calculate the values from the voltage data
# This allows to change the calibration curve or add an offset or something

channels_data_offcrr = np.copy(channels_data)

# Calculate sensor data directly from the voltages
for index1, sensor_name in enumerate(channel_names):
    # Cryofan rpm calculate from volt signal
    if sensor_name == "CryoFan (rpm)":
        channels_data_offcrr[index1] = 0
        channels_data_offcrr[index1] = V_to_RPM(sensor_name, volt_ch_data[index1])
    # Pressure sensors
    elif sensor_name == "PT850 (bara)" or sensor_name == "PDT801 (mbar)" or sensor_name == "PDT802 (mbar)":
        print("sensor name ", sensor_name)
        channels_data_offcrr[index1] = 0
        channels_data_offcrr[index1] = V_to_P(sensor_name, volt_ch_data[index1])
    # Temperature sensors
    else:
        print("sensor name ", sensor_name)
        channels_data_offcrr[index1] = 0
        channels_data_offcrr[index1] = V_to_T(sensor_name, volt_ch_data[index1])


# PLOT RAW DATA
plt.rcParams['axes.grid'] = True        # make grid on all figures by default
plt.legend(loc='upper right')

# Start of test timestamp
timestamp_0 = time[0].timestamp()
time_zeroed = np.array([time1.timestamp()-timestamp_0 for time1 in time[:]])

fig1, (T_plot, T_offcrr_plot) = plt.subplots(1, 2)
fig1.set_figheight(10)
fig1.set_figwidth(20)
# fig2, (p_plot, p_offcrr_plot) = plt.subplots(1, 2)
# fig2.set_figheight(10)
# fig2.set_figwidth(20)

# # Plot raw & corrected temperatures
# for p_data, p_data_offcrr, data_label in zip(channels_data, channels_data_offcrr, channel_names):
#     if data_label == 'PT850 (bara)':
#         p_plot.plot(time_zeroed, p_data, label='PT850')
#         p_offcrr_plot.plot(time_zeroed, p_data_offcrr, label='PT850')
#     if data_label == 'PDT801 (mbar)':
#         p_plot.plot(time_zeroed, p_data, label='PDT801')
#         p_offcrr_plot.plot(time_zeroed, p_data_offcrr, label='PDT801')
#     if data_label == 'PDT802 (mbar)':
#         p_plot.plot(time_zeroed, p_data, label='PDT802')
#         p_offcrr_plot.plot(time_zeroed, p_data_offcrr, label='PDT802')

# p_plot.set_title('Measured p')
# p_plot.set_xlabel('Time [s]')
# p_plot.set_ylabel('p [bar/mbar]')
# p_plot.legend()

# p_offcrr_plot.set_title('Offset corrected P')
# p_offcrr_plot.set_xlabel('Time [s]')
# p_offcrr_plot.set_ylabel('p [bar/mbar]')
# p_offcrr_plot.legend()


# Plot raw & corrected temperatures
for T_data, T_data_offcrr, data_label in zip(channels_data, channels_data_offcrr, channel_names):
    if data_label == 'TT855 (2nd stage)' or data_label == 'TT858 (LP in)':
        T_plot.plot(time_zeroed, T_data, label=data_label)
        T_offcrr_plot.plot(time_zeroed, T_data_offcrr, label=data_label)

T_plot.set_title('Measured T')
T_plot.set_xlabel('Time [s]')
T_plot.set_ylabel('T [K]')
T_plot.legend()

T_offcrr_plot.set_title('Offset corrected T')
T_offcrr_plot.set_xlabel('Time [s]')
T_offcrr_plot.set_ylabel('T [K]')
T_offcrr_plot.legend()




# # Manual selection of even regions
# points = plt.ginput(n=22, timeout=0)
# print(points)
# ind_temp = [int(round(i[0])) for i in points]
# ind_temp_start = ind_temp[0::2]
# ind_temp_end = ind_temp[1::2]
# print(f'start_m_cst_time = {ind_temp_start}')
# print(f'end_m_cst_time = {ind_temp_end}')
# sys.exit()
# np.savetxt('yourfilename', points)



## PostProcessing

## Moving average to act as filter
channels_data_offcrr_filter = np.copy(channels_data_offcrr)
for index, data in enumerate(channels_data_offcrr):
    channels_data_offcrr_filter[index] = moving_average(channels_data_offcrr[index])


# Calculate all average parameters per point and store in avg_par matrix
# Calculate all the standard deviation and store in st_devs matrix
avg_par = np.zeros((len(channel_names), len(start_Cav_cst_time)))
st_devs = np.zeros((len(channel_names), len(start_Cav_cst_time)))
# Calculate all the expandend uncertainty and store in uncert_e matrix
uncert_e = np.zeros((len(channel_names), len(start_Cav_cst_time)))


for index, (data_label, data) in enumerate(zip(channel_names, channels_data_offcrr_filter)):
    for index1, (time_start, time_end) in enumerate(zip(start_Cav_cst_time, end_Cav_cst_time)):
        # indices of the points where T stage is constant for m_dot for this iteration in T data file
        m_cst_range_inx = np.where((time_zeroed>time_start) & (time_zeroed<time_end))[0]
        # mT_cst_inx = cst_data_ind[np.where((cst_data_ind>m_cst_range_inx[0]) & (cst_data_ind<m_cst_range_inx[-1]))]
        # determine avg parameters
        avg_par[index][index1] = np.average(data[m_cst_range_inx])
        st_devs[index][index1] = np.std(data[m_cst_range_inx])
        uncert_e[index][index1] = Uncertainty_calc_sensor(data_label, avg_par[index][index1], st_devs[index][index1])


## Effectiveness calculation

eps_HP = np.zeros(len(start_Cav_cst_time))
eps_LP = np.zeros(len(start_Cav_cst_time))
eps_avg = np.zeros(len(start_Cav_cst_time))

for i in range(len(start_Cav_cst_time)):
    T_HP_in = avg_par[7][i] #K
    T_HP_out = avg_par[0][i] #K
    T_LP_in = avg_par[2][i] #K
    T_LP_out = avg_par[8][i] #K
    p_nominal = avg_par[13][i] #Pa
    # See diagram for the location of the pressures
    # Assumption here, that the pressure sensor PT850 is located at the lowest pressures directly before the cryofan
    # And only the pressure drop across the CFHEX is substantial in the system
    # Attention:
    # For Run1 we had to close the LP in isolate valve, thus the pressure measurement was not possible
    # The effectiveness can only be calculated with the nominal pressure
    p_HP_in = p_nominal #Pa
    p_HP_out = p_nominal #Pa
    p_LP_in = p_nominal #Pa
    p_LP_out = p_nominal #Pa

    dh_HP_max = hp.HeCalc(9, 0, 1, p_HP_in, 2, T_HP_in, 1) - hp.HeCalc(9, 0, 1, p_HP_out, 2, T_LP_in, 1)
    dh_LP_max = hp.HeCalc(9, 0, 1, p_LP_out, 2, T_HP_in, 1) - hp.HeCalc(9, 0, 1, p_LP_in, 2, T_LP_in, 1)
    dh_max = min(dh_HP_max, dh_LP_max)

    eps_HP[i] = (hp.HeCalc(9, 0, 1, p_HP_in, 2, T_HP_in, 1) - hp.HeCalc(9, 0, 1, p_HP_out, 2, T_HP_out, 1)) / dh_max
    eps_LP[i] = (hp.HeCalc(9, 0, 1, p_LP_out, 2, T_LP_out, 1) - hp.HeCalc(9, 0, 1, p_LP_in, 2, T_LP_in, 1)) / dh_max
    eps_avg[i] = 0.5 * (eps_HP[i] + eps_LP[i])



