# -*- coding: utf-8 -*-
"""
Thermodynamic model to calculate the dry cavity cooling system in tzhe CERN Cryolab.
All the components are strucutred into functions in the present file and then connected accordingly.
The piping between the components is included.
Model describes the experimental setup at the status of 02/2022.

Created on Thu Jan 27 17:17:05 2022

@author: msiodlac
"""

import numpy as np
import scipy as sp
# import matplotlib.pyplot as plt
from cryocooler_performance import T1T2_interp
from cryocooler_performance import P1P2_interp

import sys
PATH = r'C:\Users\msiodlac\cernbox\Documents\0_cern_msiodlac\python_shared\hepak-py-wrap';
sys.path.insert(0, PATH)
import hepak_wrapper as hp

##############################################################################

##############################################################################
## Definition of Functions ##

# Function: calculating the heat load of thermal radiation for a cold surface which is enclosed by a warm surface
# According to "Tieftemperaturtechnologie" by Frey and Ekin for diffuse reflection
def Q_therm_rad(A_hot, A_cold, T_hot, T_cold, epsilon_hot, epsilon_cold):
    # Emission coefficent for an enclosed vessel
    epsilon_effective = 1/ (1/epsilon_cold + A_cold/A_hot * (1/epsilon_hot - 1)) #-
    # Simple equation with an effective emission coefficent
    Q_rad = epsilon_effective * sp.constants.Stefan_Boltzmann * A_cold * (T_hot**4 - T_cold**4) #W
    return Q_rad #W

def HEX_CryoCooler(T_in, p_in, m_dot, d_inner, T_w, d_coldhead, l_HEX, Q_ex, N):
    """
    Function to calculate the temperature and pressure distribution for a capillary HEX around a coldhead of the CryoCooler.
    Heat and flow correlations are based on Torsten's old book: "Studienbrief: Wärme und Stoffübertragung in Strömungen" from Erwin Richter.


    Parameters
    ----------
    Flow conditions
        T_in : K
            Inlet temperature in K, e.g. T_in = 10 #K.
        p_in : Pa
            Inlet pressure in Pa, e.g. p_in = 3e5 #Pa.
        m_dot : kg/s
            Mass flux in kg/s, e.g. m_dot = 0.5e-3 #kg/s
            High massflow gives longest HEX.
        d_inner : m
            Inner diameter of the capillary in m, e.g. d_inner = 0.003 #m
            Largest inner diameter gives longest HEX (but lowest dP).

    Cryocooler specs
        T_w : K
            Wall temperature of the cryocooler stage in K (assumed to be constant), e.g. T_w = 5 #K.
        d_coldhead : m
            Diameter of the coldhead in m, e.g. d_coldhead = 0.1 # m.
        l_HEX : m
            Length of the HEX capillary in m, e.g. l_HEX = 1.1394 #m.
        Q_ex : W
            External heat input to the CryoCooler stage, e.g. TS which is thermalized directly onto the stage.

    Simulation specs
        N :
            Resolution of the calculation, i.e. number of increments, e.g. N = 501.

    Returns
    -------
    state_out : dict (h, T, p)
        Enthalpy, temperature, and pressure of the following state structured as a python dictionary.

    """

    # Increment length
    dx = l_HEX/(N-1) #m

    # ## Heat correlation specs according to "Studienbrief: Wärme und Stoffübertragung in Strömungen" from Erwin Richter
    # # Bending radius around the coldhead
    # R_k = 0.5 * (d_coldhead + d_inner) #Largest inner radius of bending gives largest diameter
    # Re_prime = 16.5*((R_k/d_inner)**(1/2))
    # Re_primeprime = 15200*((d_inner/R_k)**0.28)

    ## Heat correlation specs according to VDI Heat Atlas 2013
    # Bending diameter around the coldhead
    # See VDI Heatatlas: D = D_w for h < D_w, which assumed here
    D = d_coldhead + d_inner
    # Critical Reynolds number between the flow regimes
    Re_crit = 2300 * (1 + 8.6*(d_inner/D)**0.45)

    # Heat transfer areas
    A_cross = (np.pi/4)*(d_inner**2)
    A_heat = np.pi*d_inner*dx

    # Extracted heat of the cryocooler
    Q_extracted = 0

    # Allocation of temperature and pressure vector
    T = np.zeros(N, dtype = float)
    p = np.zeros(N, dtype = float)

    # Initialization with the inlet conditions
    T[0] = T_in
    p[0] = p_in

    # Start of for loop
    ##########################################################################
    for i in range(0,N-1):

        # Calculate flow properties as a function of the present T and p
        Cp     = hp.HeCalc(14, 0, 1, p[i], 2, T[i], 1) #J/(kgK)
        Rho    = hp.HeCalc(3,  0, 1, p[i], 2, T[i], 1) #kg/m³
        Mu     = hp.HeCalc(25, 0, 1, p[i], 2, T[i], 1) #Pas
        Lambda = hp.HeCalc(26, 0, 1, p[i], 2, T[i], 1) #W/(m*K)

        # Calculate dimensionless numbers
        u = m_dot / (A_cross * Rho)
        Re = (Rho * u * d_inner)/Mu
        Pr = (Cp * Mu)/Lambda

        # Prandtl number at the wall
        Cp_w     = hp.HeCalc(14, 0, 1, p[i], 2, T_w, 1) #J/(kgK)
        Mu_w     = hp.HeCalc(25, 0, 1, p[i], 2, T_w, 1) #Pas
        Lambda_w = hp.HeCalc(26, 0, 1, p[i], 2, T_w, 1) #W/(m*K)
        Pr_w = (Cp_w * Mu_w)/Lambda_w


        ## Friction coefficient f (Darcy coefficient)
        # !! Correlations according to VDI Heat Atlas 2013 !!
        # Assumption: technical smooth pipe and correlation is valid even for not isothermal flow

        # Laminar regime (equation of Mishra and Gupta)
        if 1 < Re * (d_inner/D)**0.5 and Re < Re_crit:
            f = 64/Re * (1 + 0.033 * (np.log10(Re * (d_inner/D)**0.5))**4)

        # Turbulent regime (equation of Mishra and Gupta)
        elif Re_crit < Re and Re < 10e5:
            f = 0.3164/(Re**0.25) * (1 + 0.095 * (d_inner/D)**0.5 * Re**0.25)

        else:
            raise Exception('No flow correlation for this Re regime.')


        ## Effective heat transfer coefficient alpha
        # !! Correlations according to VDI Heat Atlas 2013 !!

        # Laminar regime (equation of Schmidt and Gnielinski)
        if Re < Re_crit:
            m = 0.5 + 0.2903 * (d_inner/D)**0.194
            Nu = (3.66 + 0.08 * (1 + 0.8 * (d_inner/D)**0.9) * Re**m * Pr**(1/3)) * (Pr/Pr_w)**0.14

        # Transitional regime (linear interpolation)
        elif Re_crit < Re and Re < 2.2*10e4:

            # Nusselt number at Re = Re_crit
            m = 0.5 + 0.2903 * (d_inner/D)**0.194
            Nu_lam = (3.66 + 0.08 * (1 + 0.8 * (d_inner/D)**0.9) * Re_crit**m * Pr**(1/3)) * (Pr/Pr_w)**0.14

            # Nusselt number at Re = 2.2*10e4
            Zeta = 0.3164/(2.2*10e4**0.25) + 0.03 * (d_inner/D)**0.5
            Nu_turb = ( (Zeta/8 * 2.2*10e4 * Pr) / (1 + 12.7 * (Zeta/8)**0.5 * (Pr**(2/3) - 1)) ) * (Pr/Pr_w)**0.14

            # Mixed Nusselt number
            Gamma = (2.2*10e4 - Re)/(2.2*10e4 - Re_crit)
            Nu = Gamma * Nu_lam + (1-Gamma) * Nu_turb

        # Turbulent regime (equation of Gnielinski)
        elif 2.2*10e4 < Re:
            Zeta = 0.3164/(Re**0.25) + 0.03 * (d_inner/D)**0.5
            Nu = ( (Zeta/8 * Re * Pr) / (1 + 12.7 * (Zeta/8)**0.5 * (Pr**(2/3) - 1)) ) * (Pr/Pr_w)**0.14

        else:
            raise Exception('No heat correlation for this Re regime.')


        # Heat transfer to the wall
        alpha = Nu*Lambda/d_inner
        dQ = alpha * A_heat * (T[i] - T_w)

        ## Compressible Calculation Grohmann equations
        JT = hp.HeCalc(21, 0, 1, p[i], 2, T[i], 1) #K/Pa
        Chi = hp.HeCalc(19, 0, 1, p[i], 2, T[i], 1) #1/Pa
        Beta = 1/T[i] * hp.HeCalc(17, 0, 1, p[i], 2, T[i], 1) #1/K
        Phi = (JT*Cp + u**2*Chi) / (Cp - u**2*Beta)

        A_cross = (np.pi/4)*(d_inner**2)
        m_flux = m_dot / A_cross
        q_flux = -dQ/A_heat
        dpdx = ( ( (-f*m_flux**2)/(2*Rho*d_inner)
                  + (4*q_flux*m_flux*Beta)/(Rho*d_inner*(Cp-u**2*Beta)) )
                / ( 1 - (m_flux**2*(Chi+Beta*Phi))/Rho ) )

        dTdx = ( ( (-f*m_flux**2*Phi)/(2*Rho*d_inner)
                  + (4*q_flux)/(m_flux*d_inner*(Cp-u**2*Beta)) )
                / ( 1 - (m_flux**2*(Chi+Beta*Phi))/Rho ) )

        # Update temperature and pressure for the next increment
        p[i+1] = p[i] + dpdx * dx
        T[i+1] = T[i] + dTdx * dx


        # Extracted heat of the cryocooler
        Q_extracted += dQ

    # End of for loop
    ##########################################################################

    # Add the additional external heat which is extracted directly at the stage
    # This influences the stage temperature and therefore the fluid temperature after the stage
    Q_extracted += Q_ex

    # print("Heat extracted: ", Q_extracted, " W")

    T_out = T[N-1]
    p_out = p[N-1]
    h_out = hp.HeCalc(9, 0, 1, p_out, 2, T_out, 1) #W/kg

    state_out = {"h": h_out, "T": T_out, "p": p_out, "Q_cooling": Q_extracted}
    return state_out


def CFHX(m_dot, p_nominal, epsilon, p_HP_in, T_HP_in, p_LP_in, T_LP_in):
    """
    Function to calculate the temperature and pressure change inside a CFHX from the Cryolab.
    The results consist of the outlet pressure and temperature for the high and low pressure side.
    The calculation is based on experimental findings for the CFHX from Aleksandra Onufrena <aleksandra.onufrena@cern.ch>.

    For this particular setup the inlet parameters depends on the results therefore an iteratively calculation is compulsory.

    Parameters
    ----------
    m_dot : kg/s
        Mass flow through both the HighPressure and the LowPressure side of the CFHX.
    epsilon : -
        Effectiveness of the CFHX.
        Using a more conservative value, because we are in a different pressure regime than measured.
    epsilon = 0.965
    p_HP_in : Pa
        Inlet pressure on the high pressure side in Pa.
    T_HP_in : K
        Inlet temperature on the high pressure side in K.
    p_LP_in : Pa
        Inlet pressure on the low pressure side in Pa.
    T_LP_in : K
        Inlet temperature on the low pressure side in K.

    Returns
    -------
    state_out : dict (h_HP, T_HP, p_HP, h_LP, T_LP, p_LP)
        Enthalpy, temperature, and pressure of the outlet states on both the low and the high pressure side structured as a python dictionary.

    """


    # Data from Aleks:
        # Length CFHX = 22 cm
        # D_in = 23.6 mm, D_out 40.5 mm
        # T range: 40 K - 5 K
        # m_dot = 0.5 g/s
        # p = 1 bar
        # Effectiveness: 97.4 %
        # dp_HP = 4.8 mbar (= dp23)
        # dp_LP = 5 mbar (= dp78)

    # Geometry of the CFHX
    A_HP = 0.25 * np.pi * 0.0236**2 #m²
    A_LP = 0.25 * np.pi * (0.0405**2 - 0.0236**2) #m²


    ## Calculation of the outgoing pressure with the scaled pressure drops

    # Scaling of the pressure drop with the Darcy--Weisbach equation
    # dp = f *L/D_i * 0.5 * Rho * u**2
    dp_HP_Aleks = 4.8e2 #Pa
    dp_LP_Aleks = 5.0e2 #Pa
    # Mean density with the arithmetic mean of the temperature range values
    Rho_Aleks = hp.HeCalc(3,  0, 1, 1e5, 2, 0.5*(40+5), 1) #kg/m³
    u_HP_Aleks = 0.5e-3/(A_HP*Rho_Aleks) #m/s
    u_LP_Aleks = 0.5e-3/(A_LP*Rho_Aleks) #m/s

    # Mean density of the two inlet temperatures and the nominal pressure to be able to compare the dp
    Rho = hp.HeCalc(3,  0, 1, p_nominal, 2, 0.5*(T_HP_in + T_LP_in), 1) #kg/m³
    u_HP = m_dot/(A_HP*Rho) #m/s
    u_LP = m_dot/(A_LP*Rho) #m/s

    # Actual scaling
    dp_HP = Rho/Rho_Aleks * u_HP**2/u_HP_Aleks**2 * dp_HP_Aleks #Pa
    dp_LP = Rho/Rho_Aleks * u_LP**2/u_LP_Aleks**2 * dp_LP_Aleks #Pa

    # Calculation of the outgoing pressure with the scaled pressure drops
    p_HP_out = p_HP_in - dp_HP #Pa
    p_LP_out = p_LP_in - dp_LP #Pa


    ## Calculation of the outgoing temperatures using the effectiveness
    # Asumming that the effectiveness is the same for both the HP and the LP side!

    # Check which stream restricts the heat exchange -> Pinch point
    # See "Compact heat exchangers" by Kays, London : Chapter 7
    dh_HP_max = hp.HeCalc(9, 0, 1, p_HP_in, 2, T_HP_in, 1) - hp.HeCalc(9, 0, 1, p_HP_out, 2, T_LP_in, 1)
    dh_LP_max = hp.HeCalc(9, 0, 1, p_LP_out, 2, T_HP_in, 1) - hp.HeCalc(9, 0, 1, p_LP_in, 2, T_LP_in, 1)

    # The maximum possible heat transfer corresponds to the restricting one
    dh_max = min(dh_HP_max, dh_LP_max)

    # Calculating the specific enthalpy with all known pressures and temperatures
    h_HP_in = hp.HeCalc(9, 0, 1, p_HP_in, 2, T_HP_in, 1) #J/kg
    h_LP_in = hp.HeCalc(9, 0, 1, p_LP_in, 2, T_LP_in, 1) #J/kg

    # Calculating the outgoing enthalpies
    h_HP_out = h_HP_in - epsilon * dh_max #J/kg
    h_LP_out = h_LP_in + epsilon * dh_max #J/kg

    # Calculation of the temperatures dependend on the specific enthalpy and the pressure
    T_HP_out = hp.HeCalc(2, 0, 1, p_HP_out, 9, h_HP_out, 1) #K
    T_LP_out = hp.HeCalc(2, 0, 1, p_LP_out, 9, h_LP_out, 1) #K

    # Cross check the dp scaling
    # print("u_HP_Aleks", u_HP_Aleks)
    # print("u_HP", u_HP)
    # print("Rho_Aleks", Rho_Aleks)
    # print("Rho", Rho)
    # print("dp_HP", dp_HP)
    # print("dp_HP/dp_HP_Aleks ", dp_HP/dp_HP_Aleks)
    # print("dp_LP/dp_LP_Aleks ", dp_LP/dp_LP_Aleks)

    # Output of the results
    state_out = {"h_HP": h_HP_out, "T_HP": T_HP_out, "p_HP": p_HP_out,
                 "h_LP": h_LP_out, "T_LP": T_LP_out, "p_LP": p_LP_out}
    return state_out


def CryoFan(m_dot, p_in, T_in, p_out, T_out):
    """
    Function to calculate the temperature and pressure change inside a CryoFan.
    Noordenwind CryoFan is modelled.
    The outlet pressure is chosen beforehand and tuned with the RPM of the CryoFan.
    The outlet temperature is calculated with the inlet temperature and the heat dissipated by the CryoFan.
    The heat load consists of a friction, a dynamic, and a static fraction and is etimated with the help of the datasheet.

    The calculated outlet temperature depends already on the outlet temperature therefore an iteratively calculation is compulsory.

    Parameters
    ----------
    m_dot : kg/s
        Mass flow through the CryoFan.
    p_in : Pa
        Inlet pressure of the low pressure side in Pa.
    T_in : K
        Inlet temperature of the low pressure side in K.
    p_out : Pa
        Outlet pressure of the high pressure side in Pa.
        This value is chosen and not changed
    T_out : K
        Outlet temperature of the high pressure side in K.
        This is the value from the preciding calculation loop.

    Returns
    -------
    state_out : dict (h, T, p)
        Enthalpy, temperature, and pressure of the following state structured as a python dictionary.

    """

    def Q_Boehmwind(vol_flow, Rho):
        # Efficiency of the Boehmwind CryoFan in -.
        # Fit function of measured data for rpm = 22000.
        # Parameter: Volume flow needs to be in m³/h.
        efficiency_Boehmwind = 0.01 *(1.5962e-3*vol_flow**4 - 1.0414e-2*vol_flow**3 - 2.8084*vol_flow**2 + 2.3715e1*vol_flow + 9.1550) #-
        # Dynamic loss of the Boehmwind CryoFan in W/rho.
        # Fit function of measured data for rpm = 22000.
        # Parameter: Volume flow needs to be in m³/h.
        dynLoss_Boehmwind = -3.1011e-4*vol_flow**4 - 3.0597e-3*vol_flow**3 + 1.6961e-2*vol_flow**2 + 2.9853e-1*vol_flow + 4.6333e-2 #W/rho

        # Friction loss
        Q_friction = dynLoss_Boehmwind * Rho #W
        # Dynamic heat load
        Q_dynamic = Q_friction/efficiency_Boehmwind - Q_friction #W
        # Static heat load
        # Using the given value for operation at 30 K and 20 bara
        Q_static = 7.0 #W

        return Q_friction + Q_dynamic + Q_static

    # Calculation of a mean rho
    Rho_in = hp.HeCalc(3,  0, 1, p_in, 2, T_in, 1) #kg/m³
    Rho_out = hp.HeCalc(3,  0, 1, p_out, 2, T_out, 1) #kg/m³
    Rho = 0.5 * (Rho_in + Rho_out) #kg/m³
    # Calculation of a mean cp
    Cp_in = hp.HeCalc(14, 0, 1, p_in, 2, T_out, 1) #J/(kgK)
    Cp_out = hp.HeCalc(14, 0, 1, p_out, 2, T_out, 1) #J/(kgK)
    Cp = 0.5 * (Cp_in + Cp_out) #J/(kgK)
    # Mean volume flow
    vol_flow = m_dot / Rho * 3600 #m³/h

    ## Heat loads
    # Estimating the different heat loads that are applied on the system by the cryofan
    # Static heat load and the heat load from the fan efficiency will be dissipated across the fan.
    # Friction losses will occur in the piping of system and application.
    # Since the friction losses are small in the respective application it is assumed that all friction loss occurs at the CryoFan aswell!
    # ->Tested the friction loss in a remote cooling application and it was negligible
    # Boehmwind CryoFan
    # Call of the function for the Boehmwind CryoFan
    Q_CryoFan = Q_Boehmwind(vol_flow, Rho)

    # New temperature due to the heat load of the Cryofan
    T_out = T_in + Q_CryoFan/(Cp * m_dot)

    # Prepare the output of the results
    h_out = hp.HeCalc(9, 0, 1, p_out, 2, T_out, 1) #J/kg

    print("Cryofan heat load: ", Q_CryoFan)

    state_out = {"h": h_out, "T": T_out, "p": p_out}
    return state_out


def SimplePipe(T_in, p_in, m_dot, d_inner, l_pipe, N, BC_flag, BC_value):
    """
    Function to calculate the temperature and pressure distribution for a horizontal capillary.
    You can chose between a fixed heat load or a fixed wall temperature as boundary condition.
    Heat and flow correlations are based on VDI Heatatlas 2013 and Grohmann's equation.

    Parameters
    ----------
    Flow conditions
        T_in : K
            Inlet temperature in K, e.g. T_in = 10 #K.
        p_in : Pa
            Inlet pressure in Pa, e.g. p_in = 3e5 #Pa.
        m_dot : kg/s
            Mass flux in kg/s, e.g. m_dot = 0.5e-3 #kg/s
        d_inner : m
            Inner diameter of the capillary in m, e.g. d_inner = 0.003 #m
        l_pipe : m
            Length of the pipe in m, e.g. l_pipe = 3 #m

    Simulation specs
        N :
            Resolution of the calculation, i.e. number of increments, e.g. N = 500.

    Boundary condition
        BC_flag :
            0 : Choose 0 for a fixed pipe-area specific heat flux in W/m², specified with BC_value, e.g. BC_flag = 0. \n
            1 : Choose 1 for a fixed wall tempeature in K, specified with BC_value, e.g. BC_flag = 1.
        BC_value :
            Value which is either used as the specific heat flux or the wall temperature depending on the BC_flag value, e.g. BC_value = 5 # W/m² / K. \n
            Choose BC_flag = 0 and BC_value = 0 for adiabatic conditions!

    Returns
    -------
    state_out : dict (h, T, p)
        Enthalpy, temperature, and pressure of the following state structured as a python dictionary.

    """

    # Increment
    dx = l_pipe/(N-1) #m
    # Cross section
    A_cross = (np.pi/4)*(d_inner**2)
    # Mass flux
    m_flux = m_dot / A_cross #kg/sm²

    # Length vector for position of the increment
    x = np.zeros(N, dtype = float)
    for i in range(N):
        x[i] = i * dx

    # Allocation of temperature and pressure vector
    T = np.zeros(N, dtype = float)
    p = np.zeros(N, dtype = float)

    # Initialization with the inlet conditions
    T[0] = T_in
    p[0] = p_in

    # # Try heat load function for being more realistic
    # # Creating a linear dependency for the heat load to take the radial conduction of the cooling surface into account
    # q_flux = np.zeros(N, dtype = float)
    # for i in range(0,N):
    #     q_flux[i] = BC_value * 1.5 - BC_value/(N-1) * i

    # Average specific heat flux of the pipe: Initialization
    q_pipe = 0 #W/m²

    ###########################################################################
    ## Compressible Calculation ##

    for i in range(0,N-1):

        ## Properties at the present temperature and pressure

        #Specific heat
        Cp = hp.HeCalc(14, 0, 1, p[i], 2, T[i], 1) #J/(kgK)
        # Density
        Rho = hp.HeCalc(3,  0, 1, p[i], 2, T[i], 1) #kg
        # Dynamic viscosity
        Mu = hp.HeCalc(25, 0, 1, p[i], 2, T[i], 1) #Pas
        # Thermal conductivity coefficient
        Lambda = hp.HeCalc(26, 0, 1, p[i], 2, T[i], 1) #W/(m*K)

        # Joule-Thomson coefficient
        JT = hp.HeCalc(21, 0, 1, p[i], 2, T[i], 1) #K/Pa
        # Isothermal compressibility
        Chi = hp.HeCalc(19, 0, 1, p[i], 2, T[i], 1) #1/Pa
        # Expansivity (isobaric)
        Beta = 1/T[i] * hp.HeCalc(17, 0, 1, p[i], 2, T[i], 1) #1/K


        # Velocity
        u = m_dot / (A_cross * Rho)
        # Reynolds number
        Re = (Rho * u * d_inner)/Mu
        # Prandtl number
        Pr = (Cp * Mu)/Lambda


        ## Friction coefficient f (Darcy coefficient)
        # !! Correlations according to VDI Heat Atlas 2013 !!
        # Assumption: technical smooth pipe

        # Laminar regime
        if Re < 3000:
            f = 64/Re

        # Transitional regime (equation of Blasius)
        elif Re > 3000 and Re < 10000:
            f = 0.3164/(Re)**(1/4)

        # Turbulent regime (equation of Konakov)
        elif Re > 10000:
            f = (1.8*(np.log10(Re))-1.5)**-2

        else:
            raise Exception('No flow correlation for this Re regime.')


        ## Effective heat transfer coefficient alpha
        # !! Correlations according to VDI Heat Atlas 2013 !!

        # Laminar regime
        if Re < 2300:
            # Local Nusselt number
            Nu_local_1 = 3.66
            if x[i] == 0:
                Nu_local_2 = 1.615 * (Re * Pr * d_inner/x[i+1])**(1/3)
            else:
                Nu_local_2 = 1.615 * (Re * Pr * d_inner/x[i])**(1/3)
            Nu = (Nu_local_1**3 + 0.7**3 + (Nu_local_2-0.7)**3)**(1/3)

        # Transitional regime
        elif Re > 2300 and Re < 10000:

            # Nusselt number at Re = 2300
            Nu_lam_1 = 3.66
            if x[i] == 0:
                Nu_lam_2 = 1.615 * (2300 * Pr * d_inner/x[i+1])**(1/3)
            else:
                Nu_lam_2 = 1.615 * (2300 * Pr * d_inner/x[i])**(1/3)
            Nu_lam = (Nu_lam_1**3 + 0.7**3 + (Nu_lam_2-0.7)**3)**(1/3)

            # Nusselt number at Re = 10000
            Zeta = (1.8*(np.log10(10000))-1.5)**-2
            if x[i] == 0:
                Nu_turb = ( ( Zeta/8 * 10000 * Pr * (1+ 1/3 *(d_inner/x[i+1])**(2/3)) )
                          / ( 1 + (12.7 * (Zeta/8)**0.5 * (Pr**(2/3) - 1)) ) )
            else:
                Nu_turb = ( ( Zeta/8 * 10000 * Pr * (1+ 1/3 *(d_inner/x[i])**(2/3)) )
                          / ( 1 + (12.7 * (Zeta/8)**0.5 * (Pr**(2/3) - 1)) ) )

            # Mixed Nusselt number
            Gamma = (Re - 2300)/(10000 - 2300)
            Nu = (1-Gamma) * Nu_lam + Gamma * Nu_turb

        # Turbulent regime
        elif Re > 10000:
            Zeta = (1.8*(np.log10(Re))-1.5)**-2
            if x[i] == 0:
                Nu = ( ( Zeta/8 * Re * Pr * (1+ 1/3 *(d_inner/x[i+1])**(2/3)) )
                     / ( 1 + (12.7 * (Zeta/8)**0.5 * (Pr**(2/3) - 1)) ) )
            else:
                Nu = ( ( Zeta/8 * Re * Pr * (1+ 1/3 *(d_inner/x[i])**(2/3)) )
                     / ( 1 + (12.7 * (Zeta/8)**0.5 * (Pr**(2/3) - 1)) ) )

        else:
            raise Exception('No heat correlation for this Re regime.')

        # Heat transfer to the wall
        alpha = Nu * Lambda/d_inner

        # # Crosscheck the heat transfer coefficient
        # if i == 0 or i == N-2:
        #     print("alpha: ", alpha, " (", i, ")")

        # Calculation specific to the chosen boundary condition
        if BC_flag == 0:
            # q_flux is fixed
            q_flux = BC_value #W/m²
        else:
            # Wall temperature is fixed
            q_flux = alpha * (BC_value - T[i]) #W/m²


        ## Grohmann / Arp flow equations (see my thesis)

        # Support variable
        Phi = (JT*Cp + u**2*Chi) / (Cp - u**2*Beta)
        # Presure gradient
        dpdx = ( ( (-f*m_flux**2)/(2*Rho*d_inner)
                 + (4*q_flux*m_flux*Beta)/(Rho*d_inner*(Cp-u**2*Beta)) )
               / ( 1 - (m_flux**2*(Chi+Beta*Phi))/Rho ) )
        # Temperature gradient
        dTdx = ( ( (-f*m_flux**2*Phi)/(2*Rho*d_inner)
                 + (4*q_flux)/(m_flux*d_inner*(Cp-u**2*Beta)) )
               / ( 1 - (m_flux**2*(Chi+Beta*Phi))/Rho ) )

        ## Transfering the results to the next increment
        p[i+1] = p[i] + dpdx * dx
        T[i+1] = T[i] + dTdx * dx

        # Average specific heat flux of the pipe: Summing up
        q_pipe += q_flux

        # End of for loop
        ######################################################################

    T_out = T[N-1]
    p_out = p[N-1]
    h_out = hp.HeCalc(9, 0, 1, p_out, 2, T_out, 1) #J/kg

    # Average specific heat flux of the pipe: averaging
    q_pipe *= 1/(N-1)

    state_out = {"h": h_out, "T": T_out, "p": p_out, "x_plot": x, "T_plot": T, "q_pipe": q_pipe}
    return state_out

def FlowRestriction(T_in, p_in, m_dot_out, d_inner, f):
    """
    Function to calculate the pressure drop over a single flow restriction, e.g. arc, merger.
    Isothermal flow is always assumed.
    Calculation is based on VDI Heatatlas 2013.

    Parameters
    ----------
    Flow conditions
        T_in : K
            Inlet temperature in K, e.g. T_in = 10 #K.
        p_in : Pa
            Inlet pressure in Pa, e.g. p_in = 3e5 #Pa.
        m_dot_out : kg/s
            Mass flow of the outgoing stream in kg/s, e.g. m_dot = 0.5e-3 #kg/s
        d_inner : m
            Inner diameter of the capillary in m, e.g. d_inner = 0.003 #m
        f : m
            Darcy friction coefficient.

    Returns
    -------
    state_out : dict (h, T, p)
        Enthalpy, temperature, and pressure of the following state structured as a python dictionary.

    """

    # Cross section
    A_cross = (np.pi/4)*(d_inner**2)

    # Assumption isenthalpic flow!
    h_in = hp.HeCalc(9, 0, 1, p_in, 2, T_in, 1) #J/kg

    # Iteration for the calculation of p_out even though the influence is probably negligible
    # I checked it and for 20 bar it really is negligible
    dp = 0.0
    p_out = 0.0
    for i in range(5):
        p_out = p_in - dp
        T_out = hp.HeCalc(2, 0, 1, p_out, 9, h_in, 1)
        Rho_out = hp.HeCalc(3, 0, 1, p_out, 2, T_out, 1) #kg/m³
        # Velocity of the outgoing flow
        u_out = m_dot_out/(A_cross*Rho_out) #m/s

        # Calculation of the dp with Bernoulli equation and resistance coefficient (see VDI Heatatlas 2013)
        dp = f * Rho_out * 0.5 * u_out**2


    h_out = hp.HeCalc(9, 0, 1, p_out, 2, T_out, 1)
    state_out = {"h": h_out, "T": T_out, "p": p_out}
    return state_out

def RemoteCavity(T_in, p_in, m_dot, d_inner, l_pipe, Q_ex, N):
    """
    Function to calculate the temperature and pressure development for the cavity mockup.
    This functions combines the functions SimplePipe and FlowRestriction.
    The capillaries are modeled with the function SimplePipe.

    Parameters
    ----------
    Flow conditions
        T_in : K
            Inlet temperature in K, e.g. T_in = 10 #K.
        p_in : Pa
            Inlet pressure in Pa, e.g. p_in = 3e5 #Pa.
        m_dot : kg/s
            Mass flux in kg/s, e.g. m_dot = 0.5e-3 #kg/s
        d_inner : m
            Inner diameter of the capillary in m, e.g. d_inner = 0.003 #m
        l_pipe : m
            Length of the pipe in m, e.g. l_pipe = 3 #m
        Q_ex : W
            External heat load.

    Simulation specs
        N :
            Resolution of the calculation, i.e. number of increments, e.g. N = 100.

    Returns
    -------
    state_out : dict (h, T, p)
        Enthalpy, temperature, and pressure of the following state structured as a python dictionary.

    """

    ## Estimation of the influence of the arcs
    # Amount of 180° arcs: 5
    # Resistance coefficient for the 180° arc equal to 2*90° arc value according to VDI Heatatlas!
    f_arc = 2 * 1.3
    # Calculation according to VDI Heatatlas 2013
    # Assumption isoenthalpic flow
    state_Arc = FlowRestriction(T_in, p_in, m_dot, d_inner, 5*f_arc)
    p_Arc = state_Arc.get("p")
    T_Arc = state_Arc.get("T")

    ## Estimation of the external heat load on a compressible flow
    # Preparation of the variables to use the SimplePipe function
    # Heat transfer area of one pipe. Attention: d_inner is used!
    A_pipe = np.pi * d_inner * l_pipe #m²
    # Specific external heat load
    q_pipe = Q_ex/A_pipe #W/m²

    # Calling of the function SimplePipe
    state_out = SimplePipe(T_Arc, p_Arc, m_dot, d_inner, l_pipe, N, 0, q_pipe)
    #Transfer results
    p_out = state_out.get("p")
    T_out = state_out.get("T")
    h_out = state_out.get("h")
    state_out = {"h": h_out, "T": T_out, "p": p_out}

    return state_out


def Piping(T_in, p_in, m_dot, d_inner, l_pipe, f, epsilon_pipe, T_shield, N):
    """
    Function to calculate the temperature and pressure development inside the piping.
    This functions combines the functions SimplePipe and FlowRestriction.

    Parameters
    ----------
    Flow conditions
        T_in : K
            Inlet temperature in K, e.g. T_in = 10 #K.
        p_in : Pa
            Inlet pressure in Pa, e.g. p_in = 3e5 #Pa.
        m_dot : kg/s
            Mass flux in kg/s, e.g. m_dot = 0.5e-3 #kg/s
        d_inner : m
            Inner diameter of the capillary in m, e.g. d_inner = 0.003 #m
        l_pipe : m
            Length of the pipe in m, e.g. l_pipe = 3 #m
        f : -
            Friction loss coefficient of the combined flow restrictions.
        epsilon_pipe : -
            Emissivity coefficient of the pipe.
        T_shield : K
            Temperature of the LN2 shield in the cryostat.

    Simulation specs
        N :
            Resolution of the calculation, i.e. number of increments, e.g. N = 100.

    Returns
    -------
    state_out : dict (h, T, p)
        Enthalpy, temperature, and pressure of the following state structured as a python dictionary.

    """

    ## Estimation of the influence of the arcs
    # Calculation according to VDI Heatatlas 2013
    # Assumption isoenthalpic flow
    state_Arc = FlowRestriction(T_in, p_in, m_dot, d_inner, f)
    p_Arc = state_Arc.get("p")
    T_Arc = state_Arc.get("T")

    ## Estimation of the influence of thermal radiation on the compressible flow

    # Emission coefficent for an enclosed vessel
    # Assuming much bigger hot surface -> emissivity of hot surface doesnt matter anymore, just the cold one
    # Thus the simple equation can be used
    q_pipe = epsilon_pipe * sp.constants.Stefan_Boltzmann * (T_shield**4 - T_Arc**4) #W

    # Calling of the function SimplePipe
    state_out = SimplePipe(T_Arc, p_Arc, m_dot, d_inner, l_pipe, N, 0, q_pipe)
    #Transfer results
    p_out = state_out.get("p")
    T_out = state_out.get("T")
    h_out = state_out.get("h")
    state_out = {"h": h_out, "T": T_out, "p": p_out}

    return state_out


###############################################################################
####        #        #          #####     ######     ########   #          ####
####       # #      # #        #     #    #     #    #          #          ####
####      #   #    #   #      #       #   #      #   #          #          ####
####     #     #  #     #     #       #   #      #   #######    #          ####
####    #       ##       #     #     #    #     #    #          #          ####
####   #                  #     #####     ######     ########   ########   ####
###############################################################################
## Boundary conditions ##

## Flow conditions
# Mass flux
# m_dot = 0.47659e-3 #kg/s #Run1
m_dot = 0.59762e-3 #kg/s #Run2
# Nominal pressure
# p_nominal = 12.4135e5 #Pa #Run1
p_nominal = 18.7365e5 #Pa #Run2


## TS specs
# Shield temperature of the LN2 shield
T_shield = 130. #K

# Default inner diameter of the capillary/pipes
d_inner = 0.006 #m
# Emissivity of the pipes
epsilon_pipe = 0.05

## CryoCooler PT420 specs
# cooler_performance_data = 'data/PT420_cooler_performance.txt'
cooler_performance_data = 'data/PT420_cooler_performance.txt'
T_s1, T_s2 = T1T2_interp(cooler_performance_data)
# Diameter of the coldheads
d_coldhead_1st = 0.1295 #m
d_coldhead_2nd = 0.094 #m

## Cryocooler HEX specs
# Pipe diameter
d_i_1st = 0.016 #m
d_i_2nd = 0.006 #m
# Length of the HEX capillary around the coldheads
l_HEX_1st = 2.4 #m
l_HEX_2nd = 2.3 #m
# External heat load on the cryocooler stages from EH
# Run 1
# Q_1st_ex = 29.42627604#W
# Q_1st_ex = 21.38660321 #W
# Run 2
# Q_1st_ex = 29.78267861 #W
Q_1st_ex = 21.38660321 #W

Q_2nd_ex = 0 #9.9033446317665 #W

## CFHX specs
# Effectiveness CFHX
# Using a more conservative value, because we are in a different pressure regime than measured
eff_CFHX = 0.965

## MFM specs
Q_MFM = 0.0 #W
d_i_MFM = 0.006 #m
l_MFM = 0.7 #m
# 3 90° arcs
zeta_MFM = 3 * 1.3

## Cavity specs
# Q_Cav = 6.3 + 6.3612#W #Run1
Q_Cav = 6.3 + 5.1996#W #Run2
d_i_Cav = 0.008 #m
l_Cav = 1.45 #m
# 6 180° arcs
zeta_Cav = 6 * 2 * 1.3

## Client specs
l_Client = 0.8
d_i_Client = 0.006 #m
# 2 90° arcs
zeta_Client = 2*1.3

## Simulation specs
# Resolution
N = 200

## Initiaization
# Initialization of the values after the Cryofan
p_Cryofan_out = p_nominal #Pa
T_Cryofan_out = 54.0 #K
# Initialization of the cooling power of the two CryoCooler stages
Q_CC_1st = 30 #W
Q_CC_2nd = 4 #W
# For the convergence check of the deposited heat
Q_CC_1st_old = 20 #W
Q_CC_2nd_old = 4 #W
# Estimation of the state after Cav = CFHX LP in
p_CFHX_LP_in = p_Cryofan_out - 200 #Pa
T_CFHX_LP_in = 5.5 #K

## Simulation control
# Relaxation factor for stability
relax_factor = 0.20
# Convergence criterion
convergence_criterion_outer = 10e-5
convergence_criterion_inner = 10e-4


##############################################################################
## Calculation Single PT420 ##

# Convergence criterion
has_converged_outer = False

while not has_converged_outer:

    ##########################################################################
    ## Cryofan ##
    # p_Cryofan_out = 20e5 #Pa
    # T_Cryofan_out = 60.0 #K

    ##########################################################################
    ## Cryocooler 1st stage ##

    # Calculating the constant stage temperature depending on the cooling power of both stages
    T_CC_1st = T_s1(P1 = Q_CC_1st, P2 = Q_CC_2nd)
    # print("T_CC_1st: ", T_CC_1st)

    # Using the HEX_CryoCooler function to calculate the stage after the 1st stage depending on state 1
    state_CC_1st = HEX_CryoCooler(T_Cryofan_out, p_Cryofan_out, m_dot, d_i_1st, T_CC_1st, d_coldhead_1st, l_HEX_1st, Q_1st_ex, N)
    # Transfering the results
    p_CC_1st_out = state_CC_1st.get("p")
    T_CC_1st_out = state_CC_1st.get("T")
    Q_CC_1st = (1-relax_factor) * Q_CC_1st + relax_factor * state_CC_1st.get("Q_cooling")
    # Q_CC_1st = state_2.get("Q_cooling")

    ##########################################################################
    ## MFM / Remote cooling TS ##

    # Heat load for the MFM / TS
    # Specific heat load on the outside of the pipe
    q_MFM = (Q_MFM/(np.pi * d_i_MFM * l_MFM)) #W/m²
    state_MFM = SimplePipe(T_CC_1st_out, p_CC_1st_out, m_dot, d_i_MFM, l_MFM, 50, 0, q_MFM)
    p_MFM_out = state_MFM.get("p")
    T_MFM_out = state_MFM.get("T")

    # Piping between 1st stage and CFHX HP in
    state_CFHX_HP_in = Piping(T_MFM_out, p_MFM_out, m_dot, d_inner, 0.7, 3*1.3, epsilon_pipe, T_shield, 50)
    p_CFHX_HP_in = state_CFHX_HP_in.get("p")
    T_CFHX_HP_in = state_CFHX_HP_in.get("T")

    ##########################################################################
    ## CFHX ##
    # Need for iterativ calculation because the inlet of the CFHX depends on the outlet

    # Convergence criterion
    has_converged_inner = False

    while not has_converged_inner:

        state_CFHX = CFHX(m_dot, p_nominal, eff_CFHX, p_CFHX_HP_in, T_CFHX_HP_in, p_CFHX_LP_in, T_CFHX_LP_in)
        # Transfering the results
        p_CFHX_HP_out = state_CFHX.get("p_HP")
        T_CFHX_HP_out = state_CFHX.get("T_HP")

        ######################################################################
        ## Cryocooler 2nd stage ##

        # Piping between CFHX HP out and 2nd stage
        state_CC_2nd_in = Piping(T_CFHX_HP_out, p_CFHX_HP_out, m_dot, d_inner, 0.6, (2*1.3+4*1.3), epsilon_pipe, T_shield, 50)
        p_CC_2nd_in = state_CC_2nd_in.get("p")
        T_CC_2nd_in = state_CC_2nd_in.get("T")

        # Calculating the constant stage temperature depending on the cooling power of both stages
        T_CC_2nd = T_s2(P1 = Q_CC_1st, P2 = Q_CC_2nd)
        state_CC_2nd = HEX_CryoCooler(T_CC_2nd_in, p_CC_2nd_in, m_dot, d_i_2nd, T_CC_2nd, d_coldhead_2nd, l_HEX_2nd, Q_2nd_ex, N)
        # Transfering the results
        p_CC_2nd_out = state_CC_2nd.get("p")
        T_CC_2nd_out = state_CC_2nd.get("T")
        Q_CC_2nd = (1-relax_factor) * Q_CC_2nd + relax_factor * state_CC_2nd.get("Q_cooling")
        print("T_CC_2nd_out: ", T_CC_2nd_out, "K")


        ######################################################################
        ## Remote cooling Cavity ##

        # Piping of the client side
        state_Client_out = Piping(T_CC_2nd_out, p_CC_2nd_out, m_dot, d_inner, 0.8, 2*1.3, epsilon_pipe, T_shield, 50)
        p_Client_out = state_Client_out.get("p")
        T_Client_out = state_Client_out.get("T")

        # Remote cooling of the cavity mockup
        state_Cav = RemoteCavity(T_Client_out, p_Client_out, m_dot, d_i_Cav, l_Cav, Q_Cav, N)
        p_Cav_out = state_Cav.get("p")
        T_Cav_out = state_Cav.get("T")

        # Piping before CFHX LP in
        state_CFHX_LP_in = Piping(T_Cav_out, p_Cav_out, m_dot, d_inner, 0.55, 4*1.3, epsilon_pipe, T_shield, 50)

        # Check convergence criterion
        has_converged_inner = ( abs((state_CFHX_LP_in.get("p") - p_CFHX_LP_in)/state_CFHX_LP_in.get("p")) < convergence_criterion_inner
                            and abs((state_CFHX_LP_in.get("T") - T_CFHX_LP_in)/state_CFHX_LP_in.get("T")) < convergence_criterion_inner
                            and abs((Q_CC_2nd-Q_CC_2nd_old)/Q_CC_2nd) < convergence_criterion_inner )
        # Transfering the results
        p_CFHX_LP_in = state_CFHX_LP_in.get("p")
        T_CFHX_LP_in = state_CFHX_LP_in.get("T")
        Q_CC_2nd_old = Q_CC_2nd



    # End of inside loop
    ##########################################################################
    ## CFHX LP ##

    # Transfering the results from the CFHX earlier
    p_CFHX_LP_out = state_CFHX.get("p_LP")
    T_CFHX_LP_out = state_CFHX.get("T_LP")

    # Negligible piping in between the CFHX and the Cryofan
    p_Cryofan_in = p_CFHX_LP_out
    T_Cryofan_in = T_CFHX_LP_out

    ##########################################################################
    ## Cryofan ##
    state_Cryofan = CryoFan(m_dot, p_Cryofan_in, T_Cryofan_in, p_Cryofan_out, T_Cryofan_out)
    # Convergence check
    has_converged_outer = ( abs((state_Cryofan.get("p") - p_Cryofan_out)/state_Cryofan.get("p")) < convergence_criterion_outer
                        and abs((state_Cryofan.get("T") - T_Cryofan_out)/state_Cryofan.get("T")) < convergence_criterion_outer
                        and abs((Q_CC_1st-Q_CC_1st_old)/Q_CC_1st) < convergence_criterion_outer )
    Q_CC_1st_old = Q_CC_1st
    # Transfering the results
    p_Cryofan_out = state_Cryofan.get("p")
    T_Cryofan_out = state_Cryofan.get("T")
    print("State after Cryofan: ", state_Cryofan)

# End of the outside loop
##############################################################################

##############################################################################
## Postprocessing

# Calculation of the volume flow through Cryofan
Rho_Cryofan_in = hp.HeCalc(3,  0, 1, p_Cryofan_in, 2, T_Cryofan_in, 1) #kg/m³
Rho_Cryofan_out = hp.HeCalc(3,  0, 1, p_Cryofan_out, 2, T_Cryofan_out, 1) #kg/m³
Rho_Cryofan = 0.5 * (Rho_Cryofan_in + Rho_Cryofan_out)
vol_flow = m_dot / Rho_Cryofan * 3600 #m³/h
# print("Volume flow through Cryofan: ", vol_flow, " m³/h")
# Calculation of the pressure head in m
dp_system = p_Cryofan_out-p_Cryofan_in #Pa
press_head_system =  dp_system/(Rho_Cryofan*9.81) #m
print("Pressure head system: ", press_head_system, " m")

# Stage temperatures
print("T_CC_1st: ", T_CC_1st, " K")
print("T_CC_2nd: ", T_CC_2nd, " K")
# T_in and T_out for MFM and Cav
print("T_Cav_in: ", T_CC_2nd_out, " K")
print("T_Cav_out: ", T_Cav_out, " K")
print("T_MFM_in: ", T_CC_1st_out, " K")
print("T_MFM_out: ", T_MFM_out, " K")

##############################################################################