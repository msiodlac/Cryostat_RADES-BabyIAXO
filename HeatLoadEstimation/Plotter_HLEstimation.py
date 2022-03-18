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
## Single vs Series vs Parallel ##

# Reading the data
data = pd.read_csv('Data/data_OptiLayoutCompare.csv')

# converting column data to list
mass_flow = data['mass_flow'].tolist()
T_single = data['T_single'].tolist()
T_series = data['T_series'].tolist()
T_parallel = data['T_series'].tolist()

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
fig.subplots_adjust(left=.14, bottom=.16, right=.99, top=.97)

# Data lines, ticks, description
plt.plot(mass_flow, T_single, color='steelblue', ls='-', marker='', label = 'Single')
plt.plot(mass_flow, T_series, color='darkorange', ls=':', marker='', label = 'T CM out Series')
plt.plot(mass_flow, T_parallel, color='darkorange', ls=':', marker='', label = 'T CM out Parallel')
ax.set_ylabel('outgoing CM temperature (K)')
ax.set_xlabel('mass flow (g/s)')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_minor_locator(ML(1))
ax.yaxis.set_minor_locator(ML(25))
ax.set_xlim(0, 13)
ax.set_ylim(4.5,300)

# Final commands
# ax.grid()
plt.legend()
plt.show()
# plt.savefig("Data/Plot_SingleSeriesParallel.pdf", dpi=600)

################################################################################



################################################################################
## Tuning rod temperature distribution ##
"""
# Reading the data
data_MLI = pd.read_csv('Data/data_T_rod_MLI.csv')
data_NoMLI = pd.read_csv('Data/data_T_rod_NoMLI.csv')
# converting column data to list
x = data_MLI['x'].tolist()
T_MLI = data_MLI['T_rod'].tolist()
T_NoMLI = data_NoMLI['T_rod'].tolist()

# Size of the plot
# textwidth = 418 pt (pixels)
# pagewidth A4 = 595 pt = 21.0 cm = 8.268 in
width = 1.1 * 418 * 8.268/595 #in
# golden ratio
# height = width / 1.618
# quadratic shape
# height = width
height = width / (1.5*1.618)

fig, ax = plt.subplots(figsize=(width, height))
# fig.tight_layout()
fig.subplots_adjust(left=.14, bottom=.16, right=.99, top=.97)

# Data lines, ticks, description
plt.plot(x, T_MLI, color='steelblue', ls='-', marker='', label = 'MLI ($\epsilon = 0.06$)')
plt.plot(x, T_NoMLI, color='darkorange', ls=':', marker='', label = 'No MLI ($\epsilon = 1$)')
ax.set_ylabel('temperature tuning rod (K)')
ax.set_xlabel('length tuning rod (m)')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_minor_locator(ML(1))
ax.yaxis.set_minor_locator(ML(25))
ax.set_xlim(0, 13)
ax.set_ylim(4.5,300)

# Annotation
plt.plot([3.0,3.0], [4.5,300], color='dimgrey', ls='--', marker='')
plt.annotate('service box\nand\ntransfer line', (0.95,144))
plt.annotate('BabyIAXO bore', (3.1,175))

# Final commands
# ax.grid()
plt.legend()
plt.show()
plt.savefig("Data/Plot_TuningRod.pdf", dpi=600)
"""
################################################################################

################################################################################
## Cu coaxial cable temperature distribution ##
"""
# Reading the data
data_CoaxCu = pd.read_csv('Data/data_CoaxCable_Cu.csv')
data_CoaxCu_NoTS = pd.read_csv('Data/data_CoaxCable_Cu_NoTS.csv')
# converting column data to list
x = data_CoaxCu['x'].tolist()
T_Cu_outer = data_CoaxCu['T_Cu_outer'].tolist()
T_PTFE = data_CoaxCu['T_PTFE'].tolist()
T_Cu_inner = data_CoaxCu['T_Cu_inner'].tolist()
dT = data_CoaxCu['dT'].tolist()
T_NoTS = data_CoaxCu_NoTS['T_average'].tolist()

# Size of the plot
# textwidth = 418 pt (pixels)
# pagewidth A4 = 595 pt = 21.0 cm = 8.268 in
width = 0.55 * 418 * 8.268/595 #in
# golden ratio
# height = width / 1.618
# quadratic shape
height = width

# Structuring of the two plots inside a grid and their ticks, labels, ...
fig = plt.figure(figsize=(width,(1+1/3)*height))
fig.subplots_adjust(left=.157, bottom=.10, right=.965, top=.98)
gs = fig.add_gridspec(2, hspace=0.1, height_ratios=[3, 1])
axs = gs.subplots(sharex=True)

# Annotation
rect1 = matplotlib.patches.Rectangle((0.2,0), 0.1, 300, color='lightgray')
axs[0].add_patch(rect1)
axs[0].annotate('heat\nintercept', (0.21,225))

axs[0].plot(x, T_Cu_outer, color='steelblue', ls='-', marker='', label = 'Cu outer')
axs[0].plot(x, T_PTFE, color='darkorange', ls='-', marker='', label = 'PTFE')
axs[0].plot(x, T_Cu_inner, color='steelblue', ls='--', marker='', label = 'Cu inner')
axs[0].plot(x, T_NoTS, color='black', ls='-', marker='', label = 'No interc.')
axs[0].set_ylabel('temperature Cu coaxial cable (K)')
axs[0].xaxis.set_ticks_position('both')
axs[0].yaxis.set_ticks_position('both')
axs[0].tick_params('x', labelbottom=False)
axs[0].yaxis.set_minor_locator(ML(25))
axs[0].set_ylim(0, 300)
axs[0].legend()

rect2 = matplotlib.patches.Rectangle((0.2,0), 0.1, 10, color='lightgray')
axs[1].add_patch(rect2)
axs[1].plot(x, dT, color='black', ls='-', marker='',)
axs[1].set_ylabel('$\Delta T_{rad}$ (K)')
axs[1].set_xlabel('length coaxial cable (m)')
axs[1].xaxis.set_ticks_position('both')
axs[1].yaxis.set_ticks_position('both')
axs[1].xaxis.set_minor_locator(ML(0.1))
axs[1].yaxis.set_minor_locator(ML(2.5))
axs[1].set_xlim(0, 1)
axs[1].set_ylim(0, 10)


# Final commands
# axs[0].grid()
plt.show()
# plt.savefig("Data/Plot_CoaxCableCu.pdf", dpi=600)
"""
################################################################################

################################################################################
## SST coaxial cable temperature distribution ##
"""
# Reading the data
data_CoaxSST = pd.read_csv('Data/data_CoaxCable_SST.csv')
data_CoaxSST_NoTS = pd.read_csv('Data/data_CoaxCable_SST_NoTS.csv')
# converting column data to list
x = data_CoaxSST['x'].tolist()
T_SST = data_CoaxSST['T_SST'].tolist()
T_Cu = data_CoaxSST['T_Cu'].tolist()
T_PTFE = data_CoaxSST['T_PTFE'].tolist()
T_SBS = data_CoaxSST['T_SBS'].tolist()
dT = data_CoaxSST['dT'].tolist()
T_NoTS = data_CoaxSST_NoTS['T_average'].tolist()

# Size of the plot
# textwidth = 418 pt (pixels)
# pagewidth A4 = 595 pt = 21.0 cm = 8.268 in
width = 0.55 * 418 * 8.268/595 #in
# golden ratio
# height = width / 1.618
# quadratic shape
height = width

# Structuring of the two plots inside a grid and their ticks, labels, ...
fig = plt.figure(figsize=(width,(1+1/3)*height))
fig.subplots_adjust(left=.157, bottom=.10, right=.965, top=.98)
gs = fig.add_gridspec(2, hspace=0.1, height_ratios=[3, 1])
axs = gs.subplots(sharex=True)

# Annotation
rect1 = matplotlib.patches.Rectangle((0.45,0), 0.1, 300, color='lightgray')
axs[0].add_patch(rect1)
axs[0].annotate('heat\nintercept', (0.46,88))

axs[0].plot(x, T_SST, color='darkgreen', ls='-', marker='', label = 'SST')
axs[0].plot(x, T_Cu, color='steelblue', ls='--', marker='', label = 'Cu')
axs[0].plot(x, T_PTFE, color='darkorange', ls='-.', marker='', label = 'PTFE')
axs[0].plot(x, T_SBS, color='indigo', ls=':', marker='', label = 'SPB')
axs[0].plot(x, T_NoTS, color='black', ls='-', marker='', label = 'No interc.')
axs[0].set_ylabel('temperature SST coaxial cable (K)')
axs[0].xaxis.set_ticks_position('both')
axs[0].yaxis.set_ticks_position('both')
axs[0].tick_params('x', labelbottom=False)
axs[0].yaxis.set_minor_locator(ML(25))
axs[0].set_ylim(0, 300)
axs[0].legend()

rect2 = matplotlib.patches.Rectangle((0.45,0), 0.1, 10, color='lightgray')
axs[1].add_patch(rect2)
axs[1].plot(x, dT, color='black', ls='-', marker='',)
axs[1].set_ylabel('$\Delta T_{rad}$ (K)')
axs[1].set_xlabel('length coaxial cable (m)')
axs[1].xaxis.set_ticks_position('both')
axs[1].yaxis.set_ticks_position('both')
axs[1].xaxis.set_minor_locator(ML(0.1))
axs[1].yaxis.set_minor_locator(ML(2.5))
axs[1].set_xlim(0, 1)
axs[1].set_ylim(0, 10)


# Final commands
# axs[0].grid()
plt.show()
# plt.savefig("Data/Plot_CoaxCableSST.pdf", dpi=600)
"""

################################################################################













































################################################################################
## Boyan Plotter ##

# # saving_path = './'+str(Q_dot)+'W_'+str(d)+'mm_'+str(L)+'m_'+str(int(m_dot*1e6))+'mgs\\'
# saving_path = './'+'Plots'+'\\'

# font = {
#     'family': 'sans-serif',
#     'color':  'black',
#     'weight': 'normal',
#     'size': 11,
#     }

# font_legend = {
#     'family': 'sans-serif',
#     'weight': 'normal',
#     'size': 11,
#     }

# data_Rod =  {
#     'globals:': {
#         'font': font,
#         'font_legend': font_legend
#     },
#     'plots' : [
#             {
#                 'lines' : [
#                     {
#                         'x' : x,
#                         'y' : T_MLI,
#                         'label' : 'MLI ($\epsilon = 0.06$)',
#                         'color' : 'blue',
#                         'marker' : '',
#                         'linestyle' : 'solid'
#                     },
#                     {
#                         'x' : x,
#                         'y' : T_NoMLI,
#                         'label' : 'No MLI ($\epsilon = 1$)',
#                         'color' : 'blue',
#                         'marker' : '',
#                         'linestyle' : 'dashed'
#                     }
#                 ],
#                 # 'title' : 'Capillary: '+str(Q_dot)+' W, '+str(L)+' m, '+str(m_dot*10**6)+' mg/s',
#                 # 'description' : 'Capillary: '+str(Q_dot)+' W, '+str(L)+' m, '+str(m_dot*10**6)+' mg/s',
#                 # 'description' : ' TEST',
#                 'xlabel' : 'length tuning rod (m)',
#                 'ylabel' : 'temperature tuning rod (K)',
#                 'labelsize' : 11,
#                 'legend' : 'best',
#                 'grid' : True,
#                 'ij' : [0, 0],
#                 'xlim': (0, 13),
#                 'ylim': (4, 300)
#             }
#         ],
#    'nrows' : 1,
#    'ncols' : 1,
#    'sharex' : False,
#    'sharey' : False,
#    'figsize' : (5,3),#(12,6),
#    'save' : False,
#    'name' : 'TuningRod_TempDistribution',
#    'path' : saving_path
# }

# plotter(data_Rod)

