# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 11:06:16 2022

@author: msiodlac
"""

# import sys
# PATH = r'C:\Users\msiodlac\cernbox\Documents\0_cern_msiodlac\BoyanModels\py-plot';
# sys.path.insert(0, PATH)
# from plotter import plotter

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.ticker import MultipleLocator as ML

# from pandas import *
import pandas as pd
import numpy as np

################################################################################

# colors = ['r','g', 'b', 'k']
# markers = ['.', 'v', 's', 'd']
# linestyles = ['solid', 'dotted', 'dashed', 'dashdot']



plt.style.use('default')
# Activation of LaTeX font in matplotlib
plt.rcParams['text.usetex'] = True
# plt.rcParams('font', family='serif', serif='Charter')
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Charter'] + plt.rcParams['font.serif']
# Adjust the font size
plt.rcParams['font.size'] = 11
# Adjust the font size of just one element: title, x-axis, y-axis, legend, text, annotation
plt.rcParams['legend.fontsize'] = 11
# Direction of the tick: in, out, inout
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.size'] = 5.0
plt.rcParams['xtick.minor.size'] = 3.0
plt.rcParams['ytick.major.size'] = 5.0
plt.rcParams['ytick.minor.size'] = 3.0
# plt.rcParams['axes.linewidth'] = 1.0




# plt.annotate('Plotting style = ' + style,  xy = (1, 0.05), ha = 'left', va = 'center')


# plt.plot(err, z, s = 7, label = r'$\Sigma(x) = \gamma x^2 \sin(\theta)$')


# ax = plt.gca()

# ax.xaxis.set_minor_locator(MultipleLocator(.5))
# ax.yaxis.set_minor_locator(MultipleLocator(.005))
# plt.legend()
# plt.xlabel('$x$', labelpad = 10)
# plt.ylabel('$\phi$', labelpad = 10);
# plt.savefig('professional_plot.png', dpi = 300, pad_inches = .1, bbox_inches = 'tight')

"""
################################################################################
## Single vs Series vs Parallel T CM##

# Reading the data
data = pd.read_csv('Results/data_OptiLayoutCompare.csv')

# converting column data to list
mass_flow = data['mass_flow'].tolist()
T_single = data['T_single'].tolist()
T_series = data['T_series'].tolist()
T_parallel = data['T_parallel'].tolist()

# Size of the plot
# textwidth = 418 pt (pixels)
# pagewidth A4 = 595 pt = 21.0 cm = 8.268 in
width = 0.6 * 418 * 8.268/595 #in
# golden ratio
# height = width / 1.618
# quadratic shape
height = width


fig, ax = plt.subplots(figsize=(width, height))
# fig.tight_layout()
fig.subplots_adjust(left=.125, bottom=.12, right=.97, top=.97)

# Data lines, ticks, description
plt.plot(mass_flow, T_single, color='darkorange', ls='--', marker='o', label = 'Single')
plt.plot(mass_flow, T_parallel, color='steelblue', ls='--', marker='o', label = 'Parallel')
plt.plot(mass_flow, T_series, color='seagreen', ls='--', marker='o', label = 'Series')
ax.set_ylabel('temperature CM outlet (K)')
ax.set_xlabel('mass flow (g/s)')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_minor_locator(ML(0.1))
ax.yaxis.set_minor_locator(ML(1.0))
ax.set_xlim(0.1, 1.0)
ax.set_ylim(4, 10)

# Final commands
# ax.grid()
plt.legend()
plt.show()
# plt.savefig("Results/Plot_SingleSeriesParallel.pdf", dpi=600)

################################################################################
"""
"""
################################################################################
## Mass flow measurements ##

# Reading the data
data_Run1 = pd.read_csv('Results/data_MFT_Run1.csv')
data_Run2 = pd.read_csv('Results/data_MFT_Run2.csv')

# converting column data to list
RPM_R1 = data_Run1['RPM'].tolist()
mass_flow_R1 = data_Run1['MassFlow'].tolist()
uncert_R1 = data_Run1['Uncertainty'].tolist()
RPM_R2 = data_Run2['RPM'].tolist()
mass_flow_R2 = data_Run2['MassFlow'].tolist()
uncert_R2 = data_Run2['Uncertainty'].tolist()


# Size of the plot
# textwidth = 418 pt (pixels)
# pagewidth A4 = 595 pt = 21.0 cm = 8.268 in
width = 0.6 * 418 * 8.268/595 #in
# golden ratio
# height = width / 1.618
# quadratic shape
height = width


fig, ax = plt.subplots(figsize=(width, height))
# fig.tight_layout()
fig.subplots_adjust(left=.155, bottom=.12, right=.97, top=.97)

# Data lines, ticks, description
plt.errorbar(RPM_R1, mass_flow_R1, yerr=uncert_R1, capsize=3, elinewidth=1, color='darkorange', ls='', marker='o', label = '12 bar')
plt.errorbar(RPM_R2, mass_flow_R2, yerr=uncert_R2, capsize=3, elinewidth=1, color='steelblue', ls='', marker='o', label = '18 bar')
ax.set_ylabel('mass flow rate (g/s)')
ax.set_xlabel('circualtor impeller speed (rpm)')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_minor_locator(ML(1000.0))
ax.yaxis.set_minor_locator(ML(0.1))
ax.set_xlim(500, 7500)
ax.set_ylim(0.0, 1.3)

# Final commands
# ax.grid()
plt.legend()
plt.show()
plt.savefig("Results/Plot_Exp_MFT.pdf", dpi=600)

################################################################################


################################################################################
## Cryofan heat load measurements ##

# Reading the data
data_Run1 = pd.read_csv('Results/data_QFan_Run1.csv')
data_Run2 = pd.read_csv('Results/data_QFan_Run2.csv')

# converting column data to list
mass_flow_R1 = data_Run1['MassFlow'].tolist()
m_uncert_R1 = data_Run1['M_uncertainty'].tolist()
Q_Cryofan_R1 = data_Run1['Q_Cryofan'].tolist()
uncert_R1 = data_Run1['Uncertainty'].tolist()
Model_R1 = data_Run1['Model'].tolist()
mass_flow_R2 = data_Run2['MassFlow'].tolist()
m_uncert_R2 = data_Run2['M_uncertainty'].tolist()
Q_Cryofan_R2 = data_Run2['Q_Cryofan'].tolist()
uncert_R2 = data_Run2['Uncertainty'].tolist()
Model_R2 = data_Run2['Model'].tolist()


# Size of the plot
# textwidth = 418 pt (pixels)
# pagewidth A4 = 595 pt = 21.0 cm = 8.268 in
width = 0.6 * 418 * 8.268/595 #in
# golden ratio
# height = width / 1.618
# quadratic shape
height = width


fig, ax = plt.subplots(figsize=(width, height))
# fig.tight_layout()
fig.subplots_adjust(left=.155, bottom=.12, right=.97, top=.97)

# Data lines, ticks, description
plt.errorbar(mass_flow_R1, Q_Cryofan_R1, yerr=uncert_R1, capsize=3, elinewidth=1, color='darkorange', ls='', marker='o', label = '12 bar experiment')
plt.plot(mass_flow_R1, Model_R1, color='darkorange', ls='', marker='x', label = '12 bar model')
plt.errorbar(mass_flow_R2, Q_Cryofan_R2, yerr=uncert_R2, capsize=3, elinewidth=1, color='steelblue', ls='', marker='o', label = '18 bar experiment')
plt.plot(mass_flow_R2, Model_R2, color='steelblue', ls='', marker='x', label = '18 bar model')
ax.set_xlabel('mass flow rate (g/s)')
ax.set_ylabel('circulator heat load (W)')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_minor_locator(ML(0.1))
ax.yaxis.set_minor_locator(ML(5.0))
ax.set_xlim(0.0, 1.0)
ax.set_ylim(0., 35)

# Final commands
# ax.grid()
plt.legend()
plt.show()
plt.savefig("Results/Plot_Exp_QCryofan.pdf", dpi=600)

################################################################################


################################################################################
## Model comparison ##

# Reading the data
data = pd.read_csv('Results/data_Model_Comparison.csv')


# converting column data to list
T_Run1 = data['Run1'].tolist()
uncert_R1 = data['Uncert1'].tolist()
T_Run2 = data['Run2'].tolist()
uncert_R2 = data['Uncert2'].tolist()
Model_R1 = data['Model1'].tolist()
Model_R2 = data['Model2'].tolist()

# Size of the plot
# textwidth = 418 pt (pixels)
# pagewidth A4 = 595 pt = 21.0 cm = 8.268 in
width = 1 * 418 * 8.268/595 #in
# golden ratio
# height = width / 1.618
# quadratic shape
height = width/1.618


fig, ax = plt.subplots(figsize=(width, height))
# fig.tight_layout()
fig.subplots_adjust(left=.08, bottom=.13, right=.965, top=.97)

# Labels
x = np.arange(len(T_Run1))
text = ['TT853','TT850','TT854','TT852','TT856', 'TT855', 'TT858', 'TT853']


# Data lines, ticks, description
plt.errorbar(x, T_Run1, yerr=uncert_R1, capsize=3, elinewidth=1, color='darkorange', ls='--', marker='o', label = '12 bar, 0.5 g/s experiment')
plt.plot(x, Model_R1, color='darkorange', ls='', marker='x', label = '12 bar, 0.5 g/s model')
plt.errorbar(x, T_Run2, yerr=uncert_R2, capsize=3, elinewidth=1, color='steelblue', ls='--', marker='o', label = '18 bar, 0.6 g/s experiment')
plt.plot(x, Model_R2, color='steelblue', ls='', marker='x', label = '18 bar, 0.6 g/s model')
plt.xticks(np.arange(len(text)), text, rotation=45)
# ax.set_xlabel('mass flow (g/s)')
ax.set_ylabel('temperature (K)')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
# ax.xaxis.set_minor_locator(ML(0.1))
ax.yaxis.set_minor_locator(ML(5.0))
ax.set_xlim(0.0, 7.0)
ax.set_ylim(5., 65)

# Final commands
# ax.grid()
plt.legend()
plt.show()
plt.savefig("Results/Plot_Exp_Model_Compare.pdf", dpi=600)

################################################################################
"""








































