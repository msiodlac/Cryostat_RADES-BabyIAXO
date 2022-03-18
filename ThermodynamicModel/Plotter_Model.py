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
ax.set_xlabel('mass flow rate (g/s)')
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
plt.savefig("Results/Plot_SingleSeriesParallel.pdf", dpi=600)

################################################################################

################################################################################
## 15 bar vs 20 bar vs 25 bar: T CM out##

# Reading the data
data = pd.read_csv('Results/data_OptiPressure.csv')

# converting column data to list
mass_flow = data['mass_flow'].tolist()
bar10 = data['10bar'].tolist()
bar15 = data['15bar'].tolist()
bar20 = data['20bar'].tolist()
bar25 = data['25bar'].tolist()

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
plt.plot(mass_flow, bar10, color='darkorange', ls='--', marker='o', label = '10 bar')
plt.plot(mass_flow, bar15, color='steelblue', ls='--', marker='o', label = '15 bar')
plt.plot(mass_flow, bar20, color='seagreen', ls='--', marker='o', label = '20 bar')
plt.plot(mass_flow, bar25, color='darkred', ls='--', marker='o', label = '25 bar')
ax.set_ylabel('temperature CM outlet (K)')
ax.set_xlabel('mass flow rate (g/s)')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_minor_locator(ML(0.1))
ax.yaxis.set_minor_locator(ML(0.5))
ax.set_xlim(0.35, 1.0)
ax.set_ylim(4, 6)

# Final commands
# ax.grid()
plt.legend()
plt.show()
plt.savefig("Results/Plot_PressureOpti.pdf", dpi=600)

################################################################################

################################################################################
## Model result ##

# Reading the data
data = pd.read_csv('Results/data_Model_Proposed.csv')


# converting column data to list
x = data['Number'].tolist()
T_Model = data['Model'].tolist()

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
text = ['1','2','3','4','5', '6', '7', '8', '1']


# Data lines, ticks, description
plt.plot(x, T_Model, color='darkorange', ls='--', marker='o', label = 'model results')
plt.xticks(np.arange(len(text)), text, rotation=45)
# ax.set_xlabel('mass flow (g/s)')
ax.set_ylabel('temperature (K)')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
# ax.xaxis.set_minor_locator(ML(0.1))
ax.yaxis.set_minor_locator(ML(5.0))
ax.set_xlim(0.0, 8.0)
ax.set_ylim(0., 50)

# Final commands
# ax.grid()
plt.legend()
plt.show()
plt.savefig("Results/Plot_Model_Proposed.pdf", dpi=600)

################################################################################









































