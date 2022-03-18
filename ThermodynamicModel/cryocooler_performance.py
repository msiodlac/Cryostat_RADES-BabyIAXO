"""
TWO-STAGE CRYOCOOLER PERFORMANCE FROM MEASUREMENT DATA

This module implements a performance map of a two-stage cryocooler from measurement data.
Interpolation of the data allows for calculation of:
    - Cooling powers of both stages (if both temperatures are given)
    - Temperatures at both stages (if both cooling powers are given)
    - Temperature and cooling power of second stage (if both values for the first stage are given)
at arbitrary input values within the convex hull of the measured performance map.

Furthermore, the performance map can be plotted. 

The input file should have the following columns: First stage power, Second stage power,
first stage temperature, second stage temperature. Spaces or tabs can be the seperators. 
The first two lines are skipped and can be used for comments.

[AO] 10.03.20
- corrected code to display all the cryocooler performance lines
"""

import numpy as np
#from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import CloughTocher2DInterpolator
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker


# If the temperature values are not available for a specific combination of
# P1 and P2, then add the line anyways in the file and set the temperatures to 0.
# The code will delete that point and create the plot.
def remove_null_points(x, y):
    # if both x and y are 0, don't plot them since it is an invalid value
    x_null_values = np.where(x == 0)[0]
       
    remove_points = [];
    
    for val in x_null_values:
        if x[val] == 0 and y[val] == 0:
            remove_points.append(val)
    x = np.delete(x, remove_points)
    y = np.delete(y, remove_points)
    
    return x, y

def plot_performance_map_2(performance_data_filename, points_to_plot, path, ax):
    """Given the name of the input file this function plots the heat exchanger performance map.
    """
    
    font_legend = {
    'family': 'sans-serif',
    'weight': 'normal',
    'size': 10,
    }
    font = {
    'family': 'sans-serif',
    'color':  'black',
    'weight': 'normal',
    'size': 11,
    }    
    # Read columns from file to arrays
    P1_data, P2_data, T1_data, T2_data = np.loadtxt(performance_data_filename, skiprows=2, unpack=True)
    
    
    # Allocate arrays with numbers of different values used for the two cooling powers as dimensions
    P1 = np.zeros((len(np.unique(P1_data)), len(np.unique(P2_data))))
    P2 = np.zeros((len(np.unique(P1_data)), len(np.unique(P2_data))))
    T1 = np.zeros((len(np.unique(P1_data)), len(np.unique(P2_data))))
    T2 = np.zeros((len(np.unique(P1_data)), len(np.unique(P2_data))))

    # Fill arrays with the data such that corresponding values have the same position in the arrays
    for row, P1_value in enumerate(np.unique(P1_data)):
        for column, P2_value in enumerate(np.unique(P2_data)):
            
            select = np.logical_and(P1_data == P1_value, P2_data == P2_value)
            P1[row][column] = P1_data[select]
            P2[row][column] = P2_data[select]
            T1[row][column] = T1_data[select]
            T2[row][column] = T2_data[select]
        
    
    
    # P1
    P1_u = P1[:,0] 
    for index, P in enumerate(P1_u):
        x = T1[index, :]
        y = T2[index, :]
        
        x, y = remove_null_points(x,y)
       
        if index == 0:
            ax.annotate(str(int(P)), # this is the text
                         (x[-1],y[-1]), # this is the point to label
                         textcoords="offset points", # how to position the text
                         xytext=(0,20), # distance from text to points (x,y)
                         ha='center',
                         fontsize=font['size']) # horizontal alignment can be left, right or center
        elif index == 1:
            ax.annotate(str(int(P)), # this is the text
                         (x[-1],y[-1]), # this is the point to label
                         textcoords="offset points", # how to position the text
                         xytext=(0,10), # distance from text to points (x,y)
                         ha='center',
                         fontsize=font['size']) # horizontal alignment can be left, right or center
        elif index == 2:
            ax.annotate(str(int(P)), # this is the text
                         (x[-1],y[-1]), # this is the point to label
                         textcoords="offset points", # how to position the text
                         xytext=(0,3), # distance from text to points (x,y)
                         ha='center',
                         fontsize=font['size']) # horizontal alignment can be left, right or center
        elif index == 3:
            ax.annotate(str(int(P)), # this is the text
                         (x[-1],y[-1]), # this is the point to label
                         textcoords="offset points", # how to position the text
                         xytext=(0,10), # distance from text to points (x,y)
                         ha='center',
                         fontsize=font['size']) # horizontal alignment can be left, right or center 
        elif index == 4:
            ax.annotate(str(int(P)), # this is the text
                         (x[-1],y[-1]), # this is the point to label
                         textcoords="offset points", # how to position the text
                         xytext=(0,20), # distance from text to points (x,y)
                         ha='center',
                         fontsize=font['size']) # horizontal alignment can be left, right or center
        elif index == 5:
            ax.annotate(str(int(P)), # this is the text
                         (x[-1],y[-1]), # this is the point to label
                         textcoords="offset points", # how to position the text
                         xytext=(0,10), # distance from text to points (x,y)
                         ha='center',
                         fontsize=font['size']) # horizontal alignment can be left, right or center
        elif index == 6:
            ax.annotate(str(int(P)), # this is the text
                         (x[-1],y[-1]), # this is the point to label
                         textcoords="offset points", # how to position the text
                         xytext=(0,10), # distance from text to points (x,y)
                         ha='center',
                         fontsize=font['size']) # horizontal alignment can be left, right or center
                  

        
        lines = ax.plot(x,y, color='r', linewidth=1)        
        
    # P2
    P2_u = P2[0,:]  
    for index, P in enumerate(P2_u):
       
        x = T1[:, index]
        y = T2[:, index]
        
        x, y = remove_null_points(x,y)
        if index is not 3:
            if index is 4:                
                ax.annotate(str(int(P)), # this is the text
                     (x[0],y[0]), # this is the point to label
                     textcoords="offset points", # how to position the text
                     xytext=(-15,0), # distance from text to points (x,y)
                     ha='left',
                     fontsize=font['size']) # horizontal alignment can be left, right or center
            else:
                ax.annotate(str(int(P)), # this is the text
                     (x[0],y[0]), # this is the point to label
                     textcoords="offset points", # how to position the text
                     xytext=(-10,0), # distance from text to points (x,y)
                     ha='left',
                     fontsize=font['size']) # horizontal alignment can be left, right or center             
        else:
               ax.annotate(str(P), # this is the text
                 (x[0],y[0]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(-20,0), # distance from text to points (x,y)
                 ha='left',
                 fontsize=font['size']) # horizontal alignment can be left, right or center         

        
        lines = ax.plot(x,y, color='b', linewidth=1)
        
    for point in points_to_plot:
        ax.plot(point['T_stage_1'], point['T_stage_2'], 'go')
        
        

    
    red_patch = mpatches.Patch(color='red', label='First stage')
    blue_patch = mpatches.Patch(color='blue', label='Second stage')
    ax.legend(handles=[red_patch, blue_patch], loc='best', title='Cooling power [W]', title_fontsize=font_legend['size'], prop=font_legend)
    ax.set_xlim(sorted(T1_data)[1]*0.6, np.max(T1)*1.1)
    ax.set_ylim(sorted(T2_data)[1]*0.6, np.max(T2)*1.3)
    ax.set_xlabel('First stage temperature [K]', fontdict=font)
    ax.set_ylabel('Second stage temperature [K]', fontdict=font)
    ax.set_title('')
    ax.tick_params(labelsize = font['size'])
    
   
    start, end = ax.get_ylim()
    ax.yaxis.set_ticks(np.arange(2, 22, 2))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))    
    
    ax.margins(x=0)
    ax.grid(linewidth=0.5)

    return ax
 
    


def plot_performance_map(performance_data_filename):
    """Given the name of the input file this function plots the heat exchanger performance map.
    """
    # Read columns from file to arrays
    P1_data, P2_data, T1_data, T2_data = np.loadtxt(performance_data_filename, skiprows=2, unpack=True)

    # Allocate arrays with numbers of different values used for the two cooling powers as dimensions
    P1 = np.zeros((len(np.unique(P1_data)), len(np.unique(P2_data))))
    P2 = np.zeros((len(np.unique(P1_data)), len(np.unique(P2_data))))
    T1 = np.zeros((len(np.unique(P1_data)), len(np.unique(P2_data))))
    T2 = np.zeros((len(np.unique(P1_data)), len(np.unique(P2_data))))

    # Fill arrays with the data such that corresponding values have the same position in the arrays
    for row, P1_value in enumerate(np.unique(P1_data)):
        for column, P2_value in enumerate(np.unique(P2_data)):
            select = np.logical_and(P1_data == P1_value, P2_data == P2_value)
            P1[row][column] = P1_data[select]
            P2[row][column] = P2_data[select]
            T1[row][column] = T1_data[select]
            T2[row][column] = T2_data[select]
    
    # Calculate contour levels; upper boundary is only shown if highest level is below the boundary
    levels_P1 = np.linspace(np.amin(P1), np.amax(P1) - 1e-5, len(np.unique(P1)))    
    levels_P2 = np.linspace(np.amin(P2), np.amax(P2) - 1e-5, len(np.unique(P2)))
    
    # [AO] Correction: all boundaries shown
    levels_P1 = np.linspace(np.amin(P1) + 1e-5, np.amax(P1) - 1e-5, len(np.unique(P1)))    
    levels_P2 = np.linspace(np.amin(P2) + 1e-5, np.amax(P2) - 1e-5, len(np.unique(P2)))
    
    # Create contour plots which are iso-lines of power in a coordinate system of temperatures
    contour_P1 = plt.contour(T1, T2, P1, levels_P1, colors='r')
    contour_P2 = plt.contour(T1, T2, P2, levels_P2, colors='b')

    # Create the legend to distinguish first and second stage iso-lines
    contour_P1.collections[0].set_label('First stage')
    contour_P2.collections[0].set_label('Second stage')
    plt.legend(title='Cooling power [W]', loc='lower right')

    for i in range(len(levels_P1)):
        plt.annotate('{:g}'.format(round(levels_P1[i])), xy=(T1[i][-1] - 0.3, T2[i][-1] + 0.05), color='r')                    
    for i in range(1, len(levels_P2)):
        plt.annotate('{:.1f}'.format(levels_P2[i]), xy=(T1[0][i] - 1.2, T2[0][i] - 0.05), color='b')
    # Label individual iso-lines with their respective value

    margin_1 = (max(T1_data) - min(T1_data))/10
    margin_2 = (max(T2_data) - min(T2_data))/10
    plt.axis([min(T1_data) - margin_1, max(T1_data) + margin_1, min(T2_data) - margin_2, max(T2_data) + margin_2])
    
    # Configure the plot layout and save it to the specified file 
    plt.xlabel('First stage temperature [K]')
    plt.ylabel('Second stage temperature [K]')
#    plt.grid()
    

def P1P2_interp(performance_data_filename):
    """Given the name of the data file, this function returns two functions which allow calculation 
    of the two cooling powers given the two stage temperatures. Inputs to the returned functions
    have to be named explicitly in the function call.
    
    Sample usage:
        P1, P2 = P1P2_interp('performance_data.txt')
        P1(T1=50, T2=10)
    """
    # Read columns from file to arrays
    P1_data, P2_data, T1_data, T2_data = np.loadtxt(performance_data_filename, skiprows=2, unpack=True)   

    # Create array of points at which the values are given in the data file
    points = np.stack((T1_data, T2_data), axis=1)

    # Interpolate linearly between given data points 
    P1_interp = CloughTocher2DInterpolator(points, P1_data)
    P2_interp = CloughTocher2DInterpolator(points, P2_data)
   
    # Define functions in order to enforce named arguments
    def P1(*,T1 , T2):
        result = P1_interp(T1, T2)
        if not hasattr(T1, "__len__"):
            result = float(result)
        return result

    def P2(*,T1 , T2):
        result = P2_interp(T1, T2)
        if not hasattr(P1, "__len__"):
            result = float(result)
        return result

    return P1, P2

def T1T2_interp(performance_data_filename):
    """Given the name of the data file, this function returns two functions which allow calculation 
    of the two stage temperatures given the two cooling powers. Inputs to the returned functions
    have to be named explicitly in the function call.
    
    Sample usage:
        T1, T2 = T1T2_interp('performance_data.txt')
        T1(P1=10, P2=1)
    """
    # Read columns from file to arrays
    P1_data, P2_data, T1_data, T2_data = np.loadtxt(performance_data_filename, skiprows=2, unpack=True)   

    # Create array of points at which the values are given in the data file
    points = np.stack((P1_data, P2_data), axis=1)
    # Interpolate linearly between given data points 
    T1_interp = CloughTocher2DInterpolator(points, T1_data)
    T2_interp = CloughTocher2DInterpolator(points, T2_data)

    # Define functions in order to enforce named arguments
    def T1(*,P1 , P2):
        result = T1_interp(P1, P2)
        if not hasattr(P1, "__len__"):
            result = float(result)
        return result

    def T2(*,P1 , P2):
        result = T2_interp(P1, P2)
        if not hasattr(P1, "__len__"):
            result = float(result)
        return result

    return T1, T2

def T2P2_interp(performance_data_filename):
    """Given the name of the data file, this function returns two functions which allow calculation
    of temperature and cooling power of the second stage given the both values for the first stage.
    Inputs to the returned functions have to be named explicitly in the function call.
    
    Sample usage:
        T2, P2 = T2P2_interp('performance_data.txt')
        T2(P1=10, T1=50)
    """
    # Read columns from file to arrays
    P1_data, P2_data, T1_data, T2_data = np.loadtxt(performance_data_filename, skiprows=2, unpack=True)   

    # Create array of points at which the values are given in the data file
    points = np.stack((T1_data, P1_data), axis=1)

    # Interpolate linearly between given data points 
    T2_interp = CloughTocher2DInterpolator(points, T2_data)
    P2_interp = CloughTocher2DInterpolator(points, P2_data)

    # Define functions in order to enforce named arguments
    def T2(*,T1 , P1):
        result = T2_interp(T1, P1)
        if not hasattr(P1, "__len__"):
            result = float(result)
        return result

    def P2(*,T1 , P1):
        result = P2_interp(T1, P1)
        if not hasattr(P1, "__len__"):
            result = float(result)
        return result

    return T2, P2
