# -*- coding: utf-8 -*-
"""
Thermodynamic model to calculate the RADES-BabyIAXO cryostat.
All the components are strucutred into functions in the present file and then connected accordingly.
Different layouts are predefined and commented below: single cryocooler, two cryocoolers in series, two cryocoolers in parallel.



Created on Thu Aug 12 14:30:22 2021

@author: msiodlac
"""


import numpy as np
import scipy as sp
# import matplotlib.pyplot as plt
from cryocooler_performance import T1T2_interp

import sys
PATH = r'C:\Users\msiodlac\cernbox\Documents\0_cern_msiodlac\python_shared\hepak-py-wrap';
sys.path.insert(0, PATH)
import hepak_wrapper as hp
import csv


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

def T_func_MLI(T_hot, T_cold, pressure_vac, N_MLI):
    """
    Function to represent the temperature dependency of the heat loads.
    Here it is assumed that all heat loads follow the MLI model dependency
    First the estimated heat loads get divided by the reference value of the MLI model at the default temperatures and then they get multiplied by the different value according to the present temperatures.

    Parameters
    ----------
    T_hot : TYPE
        DESCRIPTION.
    T_cold : TYPE
        DESCRIPTION.
    pressure_vac : TYPE
        DESCRIPTION.
    N_MLI : TYPE
        DESCRIPTION.


    Returns
    -------
    q_MLI : TYPE
        DESCRIPTION.

    """

    ## Conduction and radiation contribution
    # Empirical model parameters
    alpha = 1.401*10**-4
    beta = 3.741*10**-9
    q_cond = (alpha/(N_MLI+1) * (T_hot+T_cold)/2 * (T_hot-T_cold)) #W/m²
    q_rad = (beta/(N_MLI+1) * (T_hot**4-T_cold**4)) #W/m²

    ## Residual gas contribution
    # Good assumption for helium according to Grohmann lecture and Ekin
    alpha_0 = 0.5 #-   ->conservative value
    kappa = 5/3 #-
    # Arithmetic mean for K according to Grohmann lecture
    T_mean = (T_hot + T_cold)/2.0 #K
    # "Constant" factor K
    K = (sp.constants.R / (8 * sp.constants.pi * 4*10**-3 * T_mean)) * (kappa+1)/(kappa-1) #W/(m²PaK)
    # Gas conduction
    q_gas = alpha_0 * K * pressure_vac * (T_hot - T_cold) #W/m²

    # Final equation with safety factor
    q_MLI = (q_cond + q_rad + q_gas) #W/m²

    return q_MLI

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


        # ## Heat and flow correlation according to Erwin Richter (Torsten's Book)
        # # Laminar regime
        # if Re < Re_prime:
        #     Nu_d = 3.66
        #     f = 64/Re
        # # Transitional regime
        # elif Re > Re_prime and Re < Re_primeprime:
        #     if Pr < 1:
        #         K_pr = 5 * ((2 + ((10/(Pr**2)) - 1)**0.5)**-1)
        #     elif Pr > 1:
        #         K_pr = (11/2) * (1 + (1 + 77/(4*Pr**2))**0.5)**-1
        #     C = 0.73
        #     K_k = ((d_inner/R_k)**0.25) * (1 + (2.79/((Re**(1/2))*((d_inner/R_k)**(1/4)))))
        #     Nu_d = C * (Re**(1/2)) * K_pr * K_k
        #     # See if this is adding any value, for dP in critical it is a very rough estimation.
        #     if Re > 3000:
        #         Zeta = (1.82*(np.log10(Re))-1.64)**-2
        #         f = Zeta
        #     elif Re < 3000:
        #         f = 64/Re
        # # Turbulent regime
        # elif Re > Re_primeprime:
        #     K_k = 1 + 1.8*(d_inner/R_k)
        #     Zeta = (1.82*(np.log10(Re))-1.64)**-2
        #     Nu_d = ( ((Zeta/8)*(Re - 1000)*Pr) / (1 + 12.7* ((Zeta/8)**0.5) * ((Pr**(2/3)) - 1)))  *K_k
        #     f = Zeta
        # else:
        #     raise Exception('No heat correlation for this Re regime.')


        ## Calculation

        # Heat transfer to the wall
        alpha = Nu*Lambda/d_inner
        dQ = alpha * A_heat * (T[i] - T_w)

        # ## Incompressible calculation
        # # Enthalpy going in
        # H_in = m_dot * Cp *T[i] #W
        # # Enthalpy going out
        # H_out = H_in - dQ
        # # Pressure drop in this increment
        # dp = f * ((Rho * u**2)/2) * (dx/d_inner)
        # # Update temperature and pressure for the next increment
        # T[i+1] = (H_out/m_dot) * 1/Cp
        # p[i+1] = p[i] - dp


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

    # # Enthalpy calculation
    # H_in = m_dot * hp.HeCalc(9, 0, 1, p[0], 2, T[0], 1) #W
    # H_out = m_dot * hp.HeCalc(9, 0, 1, p[N-1], 2, T[N-1], 1) #W
    # print("dH over HEX: ", H_in - H_out, " W")


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


    def Q_Noordenwind(vol_flow, Rho):
        # Pressure head as a fit function depending on the volume flow in m
        # press_head = -1.14257*vol_flow**3 - 5.61539*vol_flow**2 + 4.03585*vol_flow + 55.13512
        # Efficiency of the Noordenwind CryoFan in -.
        # Fit function of measured data for rpm = 21000.
        # Parameter: Volume flow needs to be in m³/h.
        efficiency_Noordenwind = 0.01 *(-9.8314 * vol_flow**4 + 4.8943e1 * vol_flow**3 - 9.1513e1 * vol_flow**2 + 8.8026e1 * vol_flow + 8.6729) #-
        # Dynamic loss of the Noordenwind CryoFan in W/rho.
        # Fit function of measured data for rpm = 21000.
        # Parameter: Volume flow needs to be in m³/h.
        dynLoss_Noordenwind = 6.5376e-3 * vol_flow**4 - 8.1741e-2 * vol_flow**3 + 1.8210e-1 * vol_flow**2 - 3.1761e-2 * vol_flow + 5.4652e-2 #W/rho

        # Friction loss
        Q_friction = dynLoss_Noordenwind * Rho #W
        # Dynamic heat load
        Q_dynamic = Q_friction/efficiency_Noordenwind - Q_friction #W
        # Static heat load
        # Using the given value for operation at 10 K to be conservative
        # Using a constant value from the datasheet for a working point at 30 K and 20 bara
        Q_static = 5.0 #W

        return Q_friction + Q_dynamic + Q_static


    def Q_AbsolutSystem(vol_flow, rpm, Rho, T_in):
        """
        Calculation of the total heat load from the Absolut System CryoFan.
        The heat load consists of the input heat load and the static heat load.

        Parameters
        ----------
        vol_flow :
            Mean volume flow at the CryoFan in m³/h.
        rpm :
            rpm of the CryoFan in 1/min.
        Rho :
            Mean density at the CryoFan in kg/m³.
        T_in : TYPE
            Inlet temperature at the CryoFan in K.

        Returns
        -------
        Q_CryoFan
            Total heat load from the Cryofan in W.
        """

        # Translation of the variables into the paramter space from Absolut System
        n = rpm / 60 #1/s
        Q_v = vol_flow / 3600 #m³/s
        flow_coeff = Q_v / n #m³/s /n

        # Input heat load
        # Input power as a function of the flow coefficient Q_v
        InPower = -1.86e+15 * flow_coeff**4 + 1.07e10 * flow_coeff**3 - 2.68e4 * flow_coeff**2 + 3.71e-2 * flow_coeff + 2.67e-10 #W/(rho n³)
        Q_in = InPower * n**3 * Rho #W/

        # Static heat load
        # Linear heat coefficient from the manufacturer
        # Using the inlet temperature to be conservative
        Q_static = 8e-3 * (300 - T_in) #W

        # print("vol_flow ", vol_flow)
        # print("Q_static ", Q_static)
        # print("Q_in ", Q_in)

        return Q_in + Q_static


    def Q_Cierzo(vol_flow, Rho):
        """
        Calculation of the total heat load from the Cierzo CryoFan.
        The heat load consists of the friction loss, the dynamic heat load and the static heat load.

        Parameters
        ----------
        vol_flow :
            Mean volume flow at the CryoFan in m³/h.
        Rho :
            Mean density at the CryoFan in kg/m³.

        Returns
        -------
        Q_CryoFan
            Total heat load from the Cryofan in W.
        """

        # Efficiency of the Cierzo CryoFan in -.
        # Fit function of measured data for rpm = 16000.
        # Parameter: Volume flow needs to be in m³/h.
        efficiency_Cierzo = 0.01 *(-3.1010e3 * vol_flow**4 + 4.2491e3 * vol_flow**3 - 2.5961e3 * vol_flow**2 + 7.2712e2 * vol_flow - 1.2495e1) #-

        # Dynamic loss of the Cierzo CryoFan in W/rho.
        # Fit function of measured data for rpm = 16000.
        # Parameter: Volume flow needs to be in m³/h.
        dynLoss_Cierzo = -1.4842e-1 * vol_flow**4 + 1.5072e-2 * vol_flow**3 - 4.9882e-2 * vol_flow**2 + 6.9258e-2 * vol_flow - 2.1721e-4 #W/rho


        # Friction loss
        Q_friction = dynLoss_Cierzo * Rho #W

        # Dynamic heat load
        Q_dynamic = Q_friction/efficiency_Cierzo - Q_friction #W

        # Static heat load
        # Using the given value for operation at 10 K to be conservative
        # Using a fixed value
        Q_static = 5.0 #W

        # print("vol_flow ", vol_flow)
        # print("Q_static ", Q_static)
        # print("Q_friction ", Q_friction)
        # print("Q_dynamic ", Q_dynamic)

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

    ## Noordenwind CryoFan
    # Call of the function for the Noordenwind CryoFan
    Q_CryoFan = Q_Noordenwind(vol_flow, Rho)

    ## AbsolutSystem CryoFan
    # Chose a rpm appropriate to the system flow curve
    # rpm = 9000
    # Call of the function for the Absolut System CryoFan
    # Q_CryoFan = Q_AbsolutSystem(vol_flow, rpm, Rho, T_in)

    ## Cierzo CryoFan
    # Call of the function for the Cierzo CryoFan
    # Q_CryoFan = Q_Cierzo(vol_flow, Rho)


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


def Remote_Cooling(T_in, p_in, m_dot, d_inner, l_pipe, N, BC_flag, BC_value):
    """
    Function to calculate the temperature and pressure distribution for a system consisting of flow-splitter, capillaries along a cooling surface and flow-merger.
    This functions combines the calculation for a single pipe with the flow resistance of the auxillary flow equipement.
    The capillaries are modeled with the function SimplePipe.
    You can chose between a fixed heat load or a fixed wall temperature as boundary condition which is passed on to the function Simple Pipe.

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
            Choose BC_flag = 0 and BC_value for adiabatic conditions!

    Returns
    -------
    state_out : dict (h, T, p)
        Enthalpy, temperature, and pressure of the following state structured as a python dictionary.

    """
    # See Powerpoint presentation for the legend of the numbers!

    # Split into how many pipes?
    split_number = 4
    m_dot_split = m_dot/split_number


    # Resistance coefficient for my 180° arc equal to 2*90° arc value according to VDI Heatatlas for my setup!
    f_arc = 2 * 1.3

    # Resistance coefficient for seperator and merger
    # VDI Heatatlas only has values for "Hosenrohr" seperator into two pipes and for larger diameters
    # d_i = 140 mm, f_separator = 0.035, f_merger = 0.08
    # Here conservative assumption!
    f_separator = 1
    f_merger = 1



    ## Calculation of state 1
    p_1 = p_in
    T_1 = T_in

    ## Calculation of state 2: Seperator
    # Calculation according to VDI Heatatlas 2013
    # Assumption isothermal flow
    state_2 = FlowRestriction(T_1, p_1, m_dot_split, d_inner, f_separator)
    p_2 = state_2.get("p")
    T_2 = state_2.get("T")

    ## Calculation of state 3: Cooling
    # Preparation of the variables to use the SimplePipe function
    # Heat transfer area of one pipe. Attention: d_inner is used!
    A_pipe = np.pi * d_inner * l_pipe #m²

    # Calculation specific to the chosen boundary condition
    if BC_flag == 0:
        # q_flux is fixed
        # Total heat load that all the capillaries together have to absorb
        Q_total = BC_value #W
        # Heat load that one single pipe needs to absorb.
        # There is now only one pipe into the one direction and no return pipe (See sketch).
        Q_pipe = Q_total/split_number #W
        q_pipe = Q_pipe / A_pipe
    else:
        # Wall temperature is fixed
        q_pipe = BC_value #K

    # Calling of the function SimplePipe
    state_3 = SimplePipe(T_2, p_2, m_dot_split, d_inner, l_pipe, N, BC_flag, q_pipe)
    p_3 = state_3.get("p")
    T_3 = state_3.get("T")

    ## Calculation of state 4: Arc
    # Calculation according to VDI Heatatlas 2013
    # Assumption isothermal flow
    state_4 = FlowRestriction(T_3, p_3, m_dot_split, d_inner, f_arc)
    p_4 = state_4.get("p")
    T_4 = state_4.get("T")

    ## Calculation of state 5: Merger
    # Calculation according to VDI Heatatlas 2013
    # Assumption isothermal flow
    state_5 = FlowRestriction(T_4, p_4, m_dot, d_inner, f_merger)
    p_5 = state_5.get("p")
    T_5 = state_5.get("T")

    ## Calculation of state 6: Adiabatic return
    # Assumption of adiabatic flow -> fixed q and q = 0
    # Return pipe of the pipe from state 3
    # Calling of the function SimplePipe
    state_6 = SimplePipe(T_5, p_5, m_dot, d_inner, l_pipe, N, 0, 0)

    # Get the results
    p_6 = state_6.get("p")
    T_6 = state_6.get("T")
    h_6 = state_6.get("h")

    # # Is the heat load due to friction relevant?
    # # And do we need to account for that everywhere we have a pressure drop?
    # Rho_in = hp.HeCalc(3,  0, 1, p_in, 2, T_in, 1) #kg/m³
    # Rho_out = hp.HeCalc(3,  0, 1, p_6, 2, T_6, 1) #kg/m³
    # vol_dot = m_dot/(0.5*(Rho_in+Rho_out))
    # print("Q_friction: ", (p_in-p_6)*vol_dot, " W")
    # # Last time I checked: 0.01 W -> negligible at this postion
    # # The overall friction loss is acounted for in the CryoFan


    ## Prepare the plot output of the results
    x_plot_3 = state_3.get("x_plot")
    T_plot_3 = state_3.get("T_plot")
    x_plot_6 = state_6.get("x_plot")
    T_plot_6 = state_6.get("T_plot")
    x_plot_6 = x_plot_6 + x_plot_3[len(x_plot_3)-1]
    x_plot = np.append(x_plot_3, x_plot_6)
    T_plot = np.append(T_plot_3, T_plot_6)

    # # Crosscheck the absorbed heat
    # Q_crosscheck = (state_3.get("q_pipe") + state_5.get("q_pipe")) * A_pipe * split_number
    # print("Crosscheck absorbed heat in W", Q_crosscheck)

    state_out = {"h": h_6, "T": T_6, "p": p_6, "x_plot": x_plot, "T_plot": T_plot}
    return state_out


###############################################################################
####        #        #          #####     ######     ########   #          ####
####       # #      # #        #     #    #     #    #          #          ####
####      #   #    #   #      #       #   #      #   #          #          ####
####     #     #  #     #     #       #   #      #   #######    #          ####
####    #       ##       #     #     #    #     #    #          #          ####
####   #                  #     #####     ######     ########   ########   ####
###############################################################################
# Name of the data file
# filename = "data_Series_10bar_0_2.csv"

## Boundary conditions ##

## Flow conditions
# Mass flux
m_dot = 0.5e-3 #kg/s
# Nominal pressure
p_nominal = 20e5 #Pa

# Inner diameter of the capillary/pipe
d_inner = 0.004 #m
d_i_1st = 0.008 #m


## CryoCooler PT420 specs
cooler_performance_data = 'data/PT420_cooler_performance.txt'
T_s1, T_s2 = T1T2_interp(cooler_performance_data)
# Diameter of the coldheads
d_coldhead_1st = 0.1295 #m
d_coldhead_2nd = 0.094 #m
# Length of the HEX capillary around the coldheads
l_HEX_1st = 2.0 #m
l_HEX_2nd = 2.0 #m

## CFHX specs
# Effectiveness CFHX epsilon
# Using a more conservative value, because we are in a different pressure regime than measured
epsilon = 0.965

## Heat load specs
# Vacuum pressure
pressure_vac = 1e-4 #Pa
# MLI layer
N_MLI_TS = 30
N_MLI_CM = 10
# Default temperatures
T_RT = 293 #K
T_50 = 50 #K
T_4_5 = 4.5 #K
# Mean temperatures
T_TS_mean = 50 #K
T_CM_mean = 4.5 #K
# Heat load on the thermal shield TS inside the bore @50 K and 4.5 K
Q_TS = 30.2 #W
# Heat load on the cold mass CM inside the bore @50 K and 4.5 K
Q_CM = 1.38 #W
# Heat load on the TS inside the service box @50 K and 4.5 K
Q_TS_SB = 7.1 #W
# Heat load on the thermalization plate inside the service box @50 K and 4.5 K
Q_ThermPlate = 0.31 #W
# Heat load pre factor to be multiplied with the MLI model
HL_factor_TS = Q_TS/T_func_MLI(T_RT, T_50, pressure_vac, N_MLI_TS)
HL_factor_CM = Q_CM/T_func_MLI(T_50, T_4_5, pressure_vac, N_MLI_CM)
HL_factor_TS_SB = Q_TS_SB/T_func_MLI(T_RT, T_50, pressure_vac, N_MLI_TS)
# HL_factor_ThermPlate = Q_ThermPlate/T_func_MLI(T_50, T_4_5, pressure_vac, N_MLI_CM)

## Remote cooling specs
# Length of one pipe inside the bore
# Here the layout is fixed with arcs, splitters, ...
l_pipe = 10.0 #m
# Data for the thermalization plate
# Here a straight pipe is assumed altough it will be in a Meander
l_plate = 0.5 #m

## Simulation specs
# Resolution
N = 500

## Initiaization
# Initialization of the values after the Cryofan
p_Cryofan_out = p_nominal #Pa
T_Cryofan_out = 64.0 #K
# Initialization of the cooling power of the two CryoCooler stages
Q_CC1_1st = 45 #W
Q_CC1_2nd = 4 #W
# For the convergence check of the deposited heat
Q_CC1_1st_old = 45 #W
Q_CC1_2nd_old = 4 #W
# Estimation of the state after CM = CFHX LP in
p_CM_out = p_Cryofan_out - 200 #Pa
T_CM_out = 5.5 #K

## Simulation control
# Relaxation factor for stability
relax_factor = 0.3
# Convergence criterion
# For Single PT420 Calculation: e-5 and e-4
convergence_criterion_outer = 10e-5
convergence_criterion_inner = 10e-4




"""
##############################################################################
## Calculation Default ##


## Loading of an old calculation: 0.6 g/s, Single PT420
p_Cryofan_out = 2000000.0
p_CC1_1st_out = 1999986.1618004139
p_CFHX_HP_out = 1999950.8632502756
p_Plate_out = 1999629.0336587504
p_CM_out = 1999378.3031910039
p_CFHX_LP_out = 1999342.7762619525
p_TS_out = 1997114.3062455435

T_Cryofan_out = 47.43078920826713
T_CC1_1st_out = 36.33899415216434
T_CFHX_HP_out = 9.144751197481604
T_CC1_2nd = 7.191674038664196
T_Plate_out = 7.8325794260355615
T_CM_out = 8.270243882405838
T_CFHX_LP_out = 35.25326509478185
T_TS_out = 45.773315564905026

Q_CC1_1st = 36.00894047782224
Q_CC1_2nd = 5.5112875423716075


# Convergence criterion
has_converged_outer = False

while not has_converged_outer:

    ##########################################################################
    ## Cryofan ##
    # p_Cryofan_out = 20e5 #Pa
    # T_Cryofan_out = 60.0 #K

    ##########################################################################
    ## CC1 1st stage ##

    # Calculating the constant stage temperature depending on the cooling power of both stages
    T_CC1_1st = T_s1(P1 = Q_CC1_1st, P2 = Q_CC1_2nd)
    # print("T_CC1_1st: ", T_CC1_1st)

    # Heat load for the TS service box with the current temperatures
    # Shifted temperature to be realistic
    # Here a dT of 10 K over the whole service box is assumed
    T_TS_SB_mean = T_CC1_1st + 10./2. #K
    Q_TS_SB = HL_factor_TS_SB * T_func_MLI(T_RT, T_TS_SB_mean, pressure_vac, N_MLI_TS)

    # Using the HEX_CryoCooler function to calculate the stage after the 1st stage depending on state 1
    state_CC1_1st = HEX_CryoCooler(T_Cryofan_out, p_Cryofan_out, m_dot, d_i_1st, T_CC1_1st, d_coldhead_1st, l_HEX_1st, Q_TS_SB, N)

    # Transfering the results
    p_CC1_1st_out = state_CC1_1st.get("p")
    T_CC1_1st_out = state_CC1_1st.get("T")
    Q_CC1_1st = (1-relax_factor) * Q_CC1_1st + relax_factor * state_CC1_1st.get("Q_cooling")
    # Q_CC1_1st = state_2.get("Q_cooling")

    ##########################################################################
    ## CFHX ##
    # Need for iterativ calculation because the inlet of the CFHX depends on the outlet

    # Convergence criterion
    has_converged_inner = False

    while not has_converged_inner:

        state_CFHX = CFHX(m_dot, p_nominal, epsilon, p_CC1_1st_out, T_CC1_1st_out, p_Plate_out, T_Plate_out)

        # Transfering the results
        p_CFHX_HP_out = state_CFHX.get("p_HP")
        T_CFHX_HP_out = state_CFHX.get("T_HP")

        ######################################################################
        ## CC1 2nd stage ##

        # Calculating the constant stage temperature depending on the cooling power of both stages
        T_CC1_2nd = T_s2(P1 = Q_CC1_1st, P2 = Q_CC1_2nd)
        # print("T_CC1_2nd: ", T_CC1_2nd)

        state_CC1_2nd = HEX_CryoCooler(T_CFHX_HP_out, p_CFHX_HP_out, m_dot, d_inner, T_CC1_2nd, d_coldhead_2nd, l_HEX_2nd, 0, N)

        # Transfering the results
        p_CC1_2nd_out = state_CC1_2nd.get("p")
        T_CC1_2nd_out = state_CC1_2nd.get("T")
        Q_CC1_2nd = (1-relax_factor) * Q_CC1_2nd + relax_factor * state_CC1_2nd.get("Q_cooling")

        ######################################################################
        ## Remote cooling CM ##

        # Heat load for the CM with the current temperatures
        # Arithmetic mean temperature between inlet and outlet for average temperature
        T_CM_mean = 0.5*(T_CC1_2nd_out + T_CM_out)
        Q_CM = HL_factor_CM * T_func_MLI(T_TS_mean, T_CM_mean, pressure_vac, N_MLI_CM)

        state_CM = Remote_Cooling(T_CC1_2nd_out, p_CC1_2nd_out, m_dot, d_inner, l_pipe, N, 0, Q_CM)
        # Transfering the results
        p_CM_out = state_CM.get("p")
        T_CM_out = state_CM.get("T")


        ######################################################################
        ## Cooling thermalization plate electrical equipment ##

        # Specific HL of the plate
        # Heat load is fixed and assumed to be constant
        q_plate = Q_ThermPlate/(np.pi * d_inner * l_plate) #W/m²

        state_Plate = SimplePipe(T_CC1_2nd_out, p_CC1_2nd_out, m_dot, d_inner, l_plate, N, 0, q_plate)


        # Check convergence criterion
        has_converged_inner = ( abs((state_Plate.get("p") - p_Plate_out)/state_Plate.get("p")) < convergence_criterion_inner
                            and abs((state_Plate.get("T") - T_Plate_out)/state_Plate.get("T")) < convergence_criterion_inner
                            and abs((Q_CC1_2nd-Q_CC1_2nd_old)/Q_CC1_2nd) < convergence_criterion_inner )

        # Transfering the results
        p_Plate_out = state_Plate.get("p")
        T_Plate_out = state_Plate.get("T")
        Q_CC1_2nd_old = Q_CC1_2nd



    # End of inside loop
    ##########################################################################
    ## CFHX LP ##

    # Transfering the results from the CFHX earlier
    p_CFHX_LP_out = state_CFHX.get("p_LP")
    T_CFHX_LP_out = state_CFHX.get("T_LP")

    ##########################################################################
    ## Remote cooling TS ##

    # Heat load for the TS with the current temperatures
    # Arithmetic mean temperature between inlet and outlet for average temperature
    T_TS_mean = 0.5*(T_CFHX_LP_out + T_TS_out)
    Q_TS = HL_factor_TS * T_func_MLI(T_RT, T_TS_mean, pressure_vac, N_MLI_TS)

    state_TS = Remote_Cooling(T_CFHX_LP_out, p_CFHX_LP_out, m_dot, d_inner, l_pipe, N, 0, Q_TS)
    # state_TS_out = Remote_Cooling_ReturnLine(T_CFHX_LP_out, p_CFHX_LP_out, m_dot, d_inner, l_pipe, N, 0, Q_TS)

    # Transfering the results
    p_TS_out = state_TS.get("p")
    T_TS_out = state_TS.get("T")

    ##########################################################################
    ## Piping in the SB ##

    # Pressure drop due to all the flow restrictions inside the system that are not included in remote cooling
	# Estimate how many arcs and turns there will be and then put them all together in one point
	# Circulator → CC1 1st: 90°arc, CC1 1st → CC2 1st: 90° arc, CC2 1st → CFHX: 2 90° arc, CFHX → CC1 2nd: 2 90° arc
	# CC1 2nd → CC2 2nd: 90° arc, CC2 2nd → CM: 2 90° arc, CM → Therm. Plate: 90° arc, Therm Plate Maender: 4 180° arc
	# Therm. Plate → CFHX: 2 90° arc, CFHX → TS: 2 90° arc, TS → Circulator: 90° arc
	# SUM: 15 90° arc + 4 180° arc = 23 90° arc = 23 * 1.3
    state_Piping_SB = FlowRestriction(T_TS_out, p_TS_out, m_dot, d_inner, 25*1.3)
    # Presure drop and HL due to thermal radiation of the TS in the service box
    # 2.5 m length and a pipe temperature of the CC1 2nd stage are assumed
    q_pipe_SB = (Q_therm_rad(2.45, (np.pi * d_inner * 2.5), T_TS_SB_mean, T_CC1_2nd_out, 0.05, 0.05)
                /(np.pi * d_inner * 2.5)) #W/m²
    state_Cryofan_in = SimplePipe(state_Piping_SB.get("T"), state_Piping_SB.get("p"), m_dot, d_inner, 2.5, 50, 0, q_pipe_SB)
    p_Cryofan_in = state_Cryofan_in.get("p")
    T_Cryofan_in = state_Cryofan_in.get("T")

    ##########################################################################
    ## Cryofan ##

    state_Cryofan = CryoFan(m_dot, p_Cryofan_in, T_Cryofan_in, p_Cryofan_out, T_Cryofan_out)

    # Convergence check
    has_converged_outer = ( abs((state_Cryofan.get("p") - p_Cryofan_out)/state_Cryofan.get("p")) < convergence_criterion_outer
                        and abs((state_Cryofan.get("T") - T_Cryofan_out)/state_Cryofan.get("T")) < convergence_criterion_outer
                        and abs((Q_CC1_1st-Q_CC1_1st_old)/Q_CC1_1st) < convergence_criterion_outer )
    Q_CC1_1st_old = Q_CC1_1st

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



# T_in and T_out for CM and TS
print("T_CM_in: ", T_CC1_2nd_out, " K")
print("T_CM_out: ", T_CM_out, " K")
print("T_TS_in: ", T_CFHX_LP_out, " K")
print("T_TS_out: ", T_TS_out, " K")

# Heat loads that are a function of temperature
print("Q_TS: ", Q_TS, " W")
print("Q_CM: ", Q_CM, " W")
print("Q_TS_SB: ", Q_TS_SB, " W")
print("Q_ThermPlate: ", Q_ThermPlate, " W")

##############################################################################
"""

## Volume of Helium Estimation
# Only taking into account the helium of the coldhead HX and the remote cooling

# V_remote_TS = np.pi/4 * d_inner**2 * (2 * 4 * l_pipe)
# V_remote_CM = np.pi/4 * d_inner**2 * (2 * 4 * l_pipe)
# V_1st = np.pi/4 * d_i_1st**2 * l_HEX_1st
# V_2nd = np.pi/4 * d_inner**2 * l_HEX_2nd

# V_total = V_remote_TS + V_remote_CM + V_1st + V_2nd






# pressure_head = -1.1426*vol_flow**3 - 5.6154*vol_flow**2 + 4.0359*vol_flow + 55.135
# print("pressure_head Cryofan at that v_dot: ", pressure_head)




## Plotting

# # Plotting of the CM temperature distribution

# x_plot_CM = state_6.get("x_plot")
# T_plot_CM = state_6.get("T_plot")

# x_plot_lin = [0, 20]
# T_plot_lin = [Temp_5, Temp_6]

# fig = plt.subplots(1, 1)
# plt.plot(x_plot_CM, T_plot_CM, 'r', label='T CM')
# # plt.plot(x_plot_lin, T_plot_lin, 'b--', label='Linear')
# plt.xlabel("length / m")
# plt.ylabel("temperature / K")

# plt.legend()
# plt.tight_layout()

# # Plotting of the TS temperature distribution

# x_plot_TS = state_9.get("x_plot")
# T_plot_TS = state_9.get("T_plot")

# T_plot_lin_TS = [Temp_8, Temp_9]

# fig = plt.subplots(1, 1)
# plt.plot(x_plot_TS, T_plot_TS, 'r', label='T TS')
# # plt.plot(x_plot_lin, T_plot_lin_TS, 'b--', label='Linear')
# plt.xlabel("length / m")
# plt.ylabel("temperature / K")

# plt.legend()
# plt.tight_layout()

# # -> Temperature distribution is pretty linear!

##############################################################################


"""
##############################################################################
## Calculation Single PT420 Variation ##

## Loading of an old calculation: 0.6 g/s, Single PT420
p_Cryofan_out = 2000000.0
p_CC1_1st_out = 1999986.1618004139
p_CFHX_HP_out = 1999950.8632502756
p_CM_out = 1999629.0336587504
p_Plate_out = 1999378.3031910039
p_CFHX_LP_out = 1999342.7762619525
p_TS_out = 1997114.3062455435
T_Cryofan_out = 47.43078920826713
T_CC1_1st_out = 36.33899415216434
T_CFHX_HP_out = 9.144751197481604
T_CC1_2nd = 7.191674038664196
T_CM_out = 7.8325794260355615
T_Plate_out = 8.270243882405838
T_CFHX_LP_out = 35.25326509478185
T_TS_out = 45.773315564905026
Q_CC1_1st = 36.00894047782224
Q_CC1_2nd = 5.5112875423716075


# Convergence criterion
has_converged_outer = False

while not has_converged_outer:

    ##########################################################################
    ## Cryofan ##
    # p_Cryofan_out = 20e5 #Pa
    # T_Cryofan_out = 60.0 #K

    ##########################################################################
    ## CC1 1st stage ##

    # Calculating the constant stage temperature depending on the cooling power of both stages
    T_CC1_1st = T_s1(P1 = Q_CC1_1st, P2 = Q_CC1_2nd)
    # print("T_CC1_1st: ", T_CC1_1st)
    # Heat load for the TS service box with the current temperatures
    # Shifted temperature to be realistic
    # Here a dT of 10 K over the whole service box is assumed
    T_TS_SB_mean = T_CC1_1st + 10./2. #K
    Q_TS_SB = HL_factor_TS_SB * T_func_MLI(T_RT, T_TS_SB_mean, pressure_vac, N_MLI_TS)
    # Using the HEX_CryoCooler function to calculate the stage after the 1st stage depending on state 1
    state_CC1_1st = HEX_CryoCooler(T_Cryofan_out, p_Cryofan_out, m_dot, d_i_1st, T_CC1_1st, d_coldhead_1st, l_HEX_1st, Q_TS_SB, N)
    # Transfering the results
    p_CC1_1st_out = state_CC1_1st.get("p")
    T_CC1_1st_out = state_CC1_1st.get("T")
    Q_CC1_1st = (1-relax_factor) * Q_CC1_1st + relax_factor * state_CC1_1st.get("Q_cooling")
    # Q_CC1_1st = state_2.get("Q_cooling")

    ##########################################################################
    ## CFHX ##
    # Need for iterativ calculation because the inlet of the CFHX depends on the outlet

    # Convergence criterion
    has_converged_inner = False

    while not has_converged_inner:

        state_CFHX = CFHX(m_dot, p_nominal, epsilon, p_CC1_1st_out, T_CC1_1st_out, p_Plate_out, T_Plate_out)
        # Transfering the results
        p_CFHX_HP_out = state_CFHX.get("p_HP")
        T_CFHX_HP_out = state_CFHX.get("T_HP")

        ######################################################################
        ## CC1 2nd stage ##

        # Calculating the constant stage temperature depending on the cooling power of both stages
        T_CC1_2nd = T_s2(P1 = Q_CC1_1st, P2 = Q_CC1_2nd)
        state_CC1_2nd = HEX_CryoCooler(T_CFHX_HP_out, p_CFHX_HP_out, m_dot, d_inner, T_CC1_2nd, d_coldhead_2nd, l_HEX_2nd, 0, N)
        # Transfering the results
        p_CC1_2nd_out = state_CC1_2nd.get("p")
        T_CC1_2nd_out = state_CC1_2nd.get("T")
        Q_CC1_2nd = (1-relax_factor) * Q_CC1_2nd + relax_factor * state_CC1_2nd.get("Q_cooling")

        ######################################################################
        ## Remote cooling CM ##

        # Heat load for the CM with the current temperatures
        # Arithmetic mean temperature between inlet and outlet for average temperature
        T_CM_mean = 0.5*(T_CC1_2nd_out + T_CM_out)
        Q_CM = HL_factor_CM * T_func_MLI(T_TS_mean, T_CM_mean, pressure_vac, N_MLI_CM)
        state_CM_out = Remote_Cooling(T_CC1_2nd_out, p_CC1_2nd_out, m_dot, d_inner, l_pipe, N, 0, Q_CM)
        # Transfering the results
        p_CM_out = state_CM_out.get("p")
        T_CM_out = state_CM_out.get("T")

        ######################################################################
        ## Cooling thermalization plate electrical equipment ##

        # Heat load of the thermalization plate
        # Assumed to be constant!
        # Q_CM = HL_factor_CM * T_func_MLI(T_RT, T_TS_mean, pressure_vac, N_MLI_CM)
        # Specific HL of the plate
        # Heat load is fixed and assumed to be constant
        q_plate = Q_ThermPlate/(np.pi * d_inner * l_plate) #W/m²
        state_Plate = SimplePipe(T_CM_out, p_CM_out, m_dot, d_inner, l_plate, N, 0, q_plate)

        # Check convergence criterion
        has_converged_inner = ( abs((state_Plate.get("p") - p_Plate_out)/state_Plate.get("p")) < convergence_criterion_inner
                            and abs((state_Plate.get("T") - T_Plate_out)/state_Plate.get("T")) < convergence_criterion_inner
                            and abs((Q_CC1_2nd-Q_CC1_2nd_old)/Q_CC1_2nd) < convergence_criterion_inner )
        # Transfering the results
        p_Plate_out = state_Plate.get("p")
        T_Plate_out = state_Plate.get("T")
        Q_CC1_2nd_old = Q_CC1_2nd


    # End of inside loop
    ##########################################################################
    ## CFHX LP ##

    # Transfering the results from the CFHX earlier
    p_CFHX_LP_out = state_CFHX.get("p_LP")
    T_CFHX_LP_out = state_CFHX.get("T_LP")

    ##########################################################################
    ## Remote cooling TS ##

    # Heat load for the TS with the current temperatures
    # Arithmetic mean temperature between inlet and outlet for average temperature
    T_TS_mean = 0.5*(T_CFHX_LP_out + T_TS_out)
    Q_TS = HL_factor_TS * T_func_MLI(T_RT, T_TS_mean, pressure_vac, N_MLI_TS)
    state_TS = Remote_Cooling(T_CFHX_LP_out, p_CFHX_LP_out, m_dot, d_inner, l_pipe, N, 0, Q_TS)
    # Transfering the results
    p_TS_out = state_TS.get("p")
    T_TS_out = state_TS.get("T")

    ##########################################################################
    ## Piping in the SB ##

    # Pressure drop due to all the flow restrictions inside the system that are not included in remote cooling
	# Estimate how many arcs and turns there will be and then put them all together in one point
	# Circulator → CC1 1st: 90°arc, CC1 1st → CC2 1st: 90° arc, CC2 1st → CFHX: 2 90° arc, CFHX → CC1 2nd: 2 90° arc
	# CC1 2nd → CC2 2nd: 90° arc, CC2 2nd → CM: 2 90° arc, CM → Therm. Plate: 90° arc, Therm Plate Maender: 4 180° arc
	# Therm. Plate → CFHX: 2 90° arc, CFHX → TS: 2 90° arc, TS → Circulator: 90° arc
	# SUM: 15 90° arc + 4 180° arc = 23 90° arc = 23 * 1.3
    state_Piping_SB = FlowRestriction(T_TS_out, p_TS_out, m_dot, d_inner, 25*1.3)
    # Presure drop and HL due to thermal radiation of the TS in the service box
    # 2.5 m length and a pipe temperature of the CC1 2nd stage are assumed
    q_pipe_SB = (Q_therm_rad(2.45, (np.pi * d_inner * 2.5), T_TS_SB_mean, T_CC1_2nd_out, 0.05, 0.05)
                /(np.pi * d_inner * 2.5)) #W/m²
    state_Cryofan_in = SimplePipe(state_Piping_SB.get("T"), state_Piping_SB.get("p"), m_dot, d_inner, 2.5, 50, 0, q_pipe_SB)
    p_Cryofan_in = state_Cryofan_in.get("p")
    T_Cryofan_in = state_Cryofan_in.get("T")

    ##########################################################################
    ## Cryofan ##
    state_Cryofan = CryoFan(m_dot, p_Cryofan_in, T_Cryofan_in, p_Cryofan_out, T_Cryofan_out)
    # Convergence check
    has_converged_outer = ( abs((state_Cryofan.get("p") - p_Cryofan_out)/state_Cryofan.get("p")) < convergence_criterion_outer
                        and abs((state_Cryofan.get("T") - T_Cryofan_out)/state_Cryofan.get("T")) < convergence_criterion_outer
                        and abs((Q_CC1_1st-Q_CC1_1st_old)/Q_CC1_1st) < convergence_criterion_outer )
    Q_CC1_1st_old = Q_CC1_1st
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
# T_in and T_out for CM and TS
print("T_CM_in: ", T_CC1_2nd_out, " K")
print("T_CM_out: ", T_CM_out, " K")
print("T_TS_in: ", T_CFHX_LP_out, " K")
print("T_TS_out: ", T_TS_out, " K")
# Heat loads that are a function of temperature
print("Q_TS: ", Q_TS, " W")
print("Q_CM: ", Q_CM, " W")
print("Q_TS_SB: ", Q_TS_SB, " W")
print("Q_ThermPlate: ", Q_ThermPlate, " W")
##############################################################################
"""





##############################################################################
## Calculation Two Cryocoolers Series ##

# Initialization of the cooling power of the two CryoCooler stages
Q_CC2_1st = float(Q_CC1_1st) #W
Q_CC2_2nd = float(Q_CC1_2nd) #W
# For the convergence check of the deposited heat
Q_CC2_1st_old = float(Q_CC1_1st_old) #W
Q_CC2_2nd_old = float(Q_CC1_2nd_old) #W

## Loading of an old calculation: 0.6 g/s, Series PT420
p_Cryofan_out = p_nominal/2000000.0 * 2000000.0
p_CC1_1st_out = p_nominal/2000000.0 * 1999965.9819977405
p_CFHX_HP_out = p_nominal/2000000.0 * 1999895.6769761513
p_Plate_out = p_nominal/2000000.0 * 1998753.0260555905
p_CM_out = p_nominal/2000000.0 * 1998786.6854423531
p_CFHX_LP_out = p_nominal/2000000.0 * 1998710.5655103405
p_TS_out = p_nominal/2000000.0 * 1994683.8950614152

T_Cryofan_out = 32.87885346721496
T_CC1_1st_out = 25.680754036159442
T_CFHX_HP_out = 6.269904305908413
T_CC1_2nd = 5.283034360592272
T_Plate_out = 5.018369709932516
T_CM_out = 4.890050618552861
T_CFHX_LP_out = 24.935341995913983
T_TS_out = 31.64677087371011

Q_CC1_1st = 16.005104660339025
Q_CC1_2nd = 2.77613079229747


# Convergence criterion
has_converged_outer = False

while not has_converged_outer:
    ##########################################################################
    ## Cryofan ##
    # p_Cryofan_out = 20e5 #Pa
    # T_Cryofan_out = 60.0 #K

    ##########################################################################
    ## CC1 1st stage ##

    # Calculating the constant stage temperature depending on the cooling power of both stages
    T_CC1_1st = T_s1(P1 = Q_CC1_1st, P2 = Q_CC1_2nd)
    # Using the HEX_CryoCooler function to calculate the stage after the 1st stage depending on state 1
    state_CC1_1st = HEX_CryoCooler(T_Cryofan_out, p_Cryofan_out, m_dot, d_i_1st, T_CC1_1st, d_coldhead_1st, l_HEX_1st, 0, N)
    # Transfering the results
    p_CC1_1st_out = state_CC1_1st.get("p")
    T_CC1_1st_out = state_CC1_1st.get("T")
    Q_CC1_1st = (1-relax_factor) * Q_CC1_1st + relax_factor * state_CC1_1st.get("Q_cooling")

    ##########################################################################
    ## CC2 1st stage ##

    # Calculating the constant stage temperature depending on the cooling power of both stages
    T_CC2_1st = T_s1(P1 = Q_CC2_1st, P2 = Q_CC2_2nd)

    ## Heat load of service box TS
    # Heat load for the TS service box with the current temperatures
    # Shifted temperature to be realistic
    # Here a dT of 10 K over the whole service box is assumed
    T_TS_SB_mean = T_CC2_1st + 10./2. #K
    Q_TS_SB = HL_factor_TS_SB * T_func_MLI(T_RT, T_TS_SB_mean, pressure_vac, N_MLI_TS)

    # Using the HEX_CryoCooler function to calculate the stage after the 1st stage depending on state 1
    state_CC2_1st = HEX_CryoCooler(T_CC1_1st_out, p_CC1_1st_out, m_dot, d_i_1st, T_CC2_1st, d_coldhead_1st, l_HEX_1st, Q_TS_SB, N)
    # Transfering the results
    p_CC2_1st_out = state_CC2_1st.get("p")
    T_CC2_1st_out = state_CC2_1st.get("T")
    Q_CC2_1st = (1-relax_factor) * Q_CC2_1st + relax_factor * state_CC2_1st.get("Q_cooling")

    ##########################################################################
    ## CFHX ##
    # Need for iterativ calculation because the inlet of the CFHX depends on the outlet

    # Convergence criterion
    has_converged_inner = False

    while not has_converged_inner:

        state_CFHX = CFHX(m_dot, p_nominal, epsilon, p_CC2_1st_out, T_CC2_1st_out, p_Plate_out, T_Plate_out)

        # Transfering the results
        p_CFHX_HP_out = state_CFHX.get("p_HP")
        T_CFHX_HP_out = state_CFHX.get("T_HP")

        ######################################################################
        ## CC1 2nd stage ##

        # Calculating the constant stage temperature depending on the cooling power of both stages
        T_CC1_2nd = T_s2(P1 = Q_CC1_1st, P2 = Q_CC1_2nd)
        state_CC1_2nd = HEX_CryoCooler(T_CFHX_HP_out, p_CFHX_HP_out, m_dot, d_inner, T_CC1_2nd, d_coldhead_2nd, l_HEX_2nd, 0, N)
        # Transfering the results
        p_CC1_2nd_out = state_CC1_2nd.get("p")
        T_CC1_2nd_out = state_CC1_2nd.get("T")
        Q_CC1_2nd = (1-relax_factor) * Q_CC1_2nd + relax_factor * state_CC1_2nd.get("Q_cooling")

        ######################################################################
        ## CC2 2nd stage ##

        # Calculating the constant stage temperature depending on the cooling power of both stages
        T_CC2_2nd = T_s2(P1 = Q_CC2_1st, P2 = Q_CC2_2nd)
        state_CC2_2nd = HEX_CryoCooler(T_CC1_2nd_out, p_CC1_2nd_out, m_dot, d_inner, T_CC2_2nd, d_coldhead_2nd, l_HEX_2nd, 0, N)
        # Transfering the results
        p_CC2_2nd_out = state_CC2_2nd.get("p")
        T_CC2_2nd_out = state_CC2_2nd.get("T")
        Q_CC2_2nd = (1-relax_factor) * Q_CC2_2nd + relax_factor * state_CC2_2nd.get("Q_cooling")

        ######################################################################
        ## Remote cooling CM ##

        # Heat load for the CM with the current temperatures
        # Arithmetic mean temperature between inlet and outlet for average temperature
        T_CM_mean = 0.5*(T_CC2_2nd_out + T_CM_out)
        Q_CM = HL_factor_CM * T_func_MLI(T_TS_mean, T_CM_mean, pressure_vac, N_MLI_CM)
        state_CM_out = Remote_Cooling(T_CC2_2nd_out, p_CC2_2nd_out, m_dot, d_inner, l_pipe, N, 0, Q_CM)
        # Transfering the results
        p_CM_out = state_CM_out.get("p")
        T_CM_out = state_CM_out.get("T")

        ######################################################################
        ## Cooling thermalization plate electrical equipment ##

        # Heat load of the thermalization plate
        # Assumed to be constant!
        # Q_CM = HL_factor_CM * T_func_MLI(T_RT, T_TS_mean, pressure_vac, N_MLI_CM)
        # Specific HL of the plate
        # Heat load is fixed and assumed to be constant
        q_plate = Q_ThermPlate/(np.pi * d_inner * l_plate) #W/m²
        state_Plate = SimplePipe(T_CM_out, p_CM_out, m_dot, d_inner, l_plate, N, 0, q_plate)

        # Check convergence criterion
        has_converged_inner = ( abs((state_Plate.get("p") - p_Plate_out)/state_Plate.get("p")) < convergence_criterion_inner
                            and abs((state_Plate.get("T") - T_Plate_out)/state_Plate.get("T")) < convergence_criterion_inner
                            and abs((Q_CC1_2nd-Q_CC1_2nd_old)/Q_CC1_2nd) < convergence_criterion_inner
                            and abs((Q_CC2_2nd-Q_CC2_2nd_old)/Q_CC2_2nd) < convergence_criterion_inner )
        # Transfering the results
        p_Plate_out = state_Plate.get("p")
        T_Plate_out = state_Plate.get("T")
        Q_CC1_2nd_old = Q_CC1_2nd
        Q_CC2_2nd_old = Q_CC2_2nd



    # End of inside loop
    ##########################################################################
    ## CFHX LP ##

    # Transfering the results from the CFHX earlier
    p_CFHX_LP_out = state_CFHX.get("p_LP")
    T_CFHX_LP_out = state_CFHX.get("T_LP")

    ##########################################################################
    ## Remote cooling TS ##

    # Heat load for the TS with the current temperatures
    # Arithmetic mean temperature between inlet and outlet for average temperature
    T_TS_mean = 0.5*(T_CFHX_LP_out + T_TS_out)
    Q_TS = HL_factor_TS * T_func_MLI(T_RT, T_TS_mean, pressure_vac, N_MLI_TS)
    state_TS = Remote_Cooling(T_CFHX_LP_out, p_CFHX_LP_out, m_dot, d_inner, l_pipe, N, 0, Q_TS)
    # Transfering the results
    p_TS_out = state_TS.get("p")
    T_TS_out = state_TS.get("T")

    ##########################################################################
    ## Piping in the SB ##

    # Pressure drop due to all the flow restrictions inside the system that are not included in remote cooling
	# Estimate how many arcs and turns there will be and then put them all together in one point
	# Circulator → CC1 1st: 90°arc, CC1 1st → CC2 1st: 90° arc, CC2 1st → CFHX: 2 90° arc, CFHX → CC1 2nd: 2 90° arc
	# CC1 2nd → CC2 2nd: 90° arc, CC2 2nd → CM: 2 90° arc, CM → Therm. Plate: 90° arc, Therm Plate Maender: 4 180° arc
	# Therm. Plate → CFHX: 2 90° arc, CFHX → TS: 2 90° arc, TS → Circulator: 90° arc
	# SUM: 15 90° arc + 4 180° arc = 23 90° arc = 23 * 1.3
    state_Piping_SB = FlowRestriction(T_TS_out, p_TS_out, m_dot, d_inner, 23*1.3)
    # Presure drop and HL due to thermal radiation of the TS in the service box
    # 2.5 m length and a pipe temperature of the CC1 2nd stage are assumed
    q_pipe_SB = (Q_therm_rad(2.45, (np.pi * d_inner * 2.5), T_TS_SB_mean, T_CC1_2nd_out, 0.05, 0.05)
                /(np.pi * d_inner * 2.5)) #W/m²
    state_Cryofan_in = SimplePipe(state_Piping_SB.get("T"), state_Piping_SB.get("p"), m_dot, d_inner, 2.5, 50, 0, q_pipe_SB)
    p_Cryofan_in = state_Cryofan_in.get("p")
    T_Cryofan_in = state_Cryofan_in.get("T")

    ##########################################################################
    ## Cryofan ##

    state_Cryofan = CryoFan(m_dot, p_Cryofan_in, T_Cryofan_in, p_Cryofan_out, T_Cryofan_out)

    # Convergence check
    has_converged_outer = (   abs((state_Cryofan.get("p") - p_Cryofan_out)/state_Cryofan.get("p")) < convergence_criterion_outer
                          and abs((state_Cryofan.get("T") - T_Cryofan_out)/state_Cryofan.get("T")) < convergence_criterion_outer
                          and abs((Q_CC1_1st-Q_CC1_1st_old)/Q_CC1_1st) < convergence_criterion_outer
                          and abs((Q_CC2_1st-Q_CC2_1st_old)/Q_CC2_1st) < convergence_criterion_outer )
    Q_CC1_1st_old = Q_CC1_1st
    Q_CC2_1st_old = Q_CC2_1st

    # Transfering the results
    p_Cryofan_out = state_Cryofan.get("p")
    T_Cryofan_out = state_Cryofan.get("T")

    print("State after Cryofan: ", state_Cryofan)

# End of the outside loop
##############################################################################
## Postprocessing ##

# Calculation of the volume flow through Cryofan
Rho_Cryofan_in = hp.HeCalc(3,  0, 1, p_Cryofan_in, 2, T_Cryofan_in, 1) #kg/m³
Rho_Cryofan_out = hp.HeCalc(3,  0, 1, p_Cryofan_out, 2, T_Cryofan_out, 1) #kg/m³
Rho_Cryofan = 0.5 * (Rho_Cryofan_in + Rho_Cryofan_out)
vol_flow = m_dot / Rho_Cryofan * 3600 #m³/h
# Calculation of the pressure head in m
dp_system = p_Cryofan_out-p_Cryofan_in #Pa
press_head_system =  dp_system/(Rho_Cryofan*9.81)
print("Pressure head system: ", press_head_system, " m")

# T_in and T_out for CM and TS
print("T_CM_in: ", T_CC2_2nd_out, " K")
print("T_CM_out: ", T_CM_out, " K")
print("T_TS_in: ", T_CFHX_LP_out, " K")
print("T_TS_out: ", T_TS_out, " K")

# Heat loads that are a function of temperature
print("Q_TS: ", Q_TS, " W")
print("Q_CM: ", Q_CM, " W")
print("Q_TS_SB: ", Q_TS_SB, " W")
print("Q_ThermPlate: ", Q_ThermPlate, " W")

##############################################################################
"""
#Save the data
# Data of the file
fields = ['Description', 'Value']
# Restructuring into rows
rows = np.zeros((24,2), dtype = float)
rows[0,1] = m_dot
rows[1,1] = p_nominal
rows[2,1] = p_Cryofan_in
rows[3,1] = press_head_system
rows[4,1] = vol_flow
rows[5,1] = T_CM_out
rows[6,1] = T_TS_out
rows[7,1] = T_CC1_1st_out
rows[8,1] = T_CC2_1st_out
rows[9,1] = T_CC1_2nd_out
rows[10,1] = T_CC2_2nd_out
rows[11,1] = T_CFHX_HP_out
rows[12,1] = T_CFHX_LP_out
rows[13,1] = T_Cryofan_in
rows[14,1] = T_Cryofan_out
rows[15,1] = T_Plate_out
rows[16,1] = Q_TS
rows[17,1] = Q_CM
rows[18,1] = Q_TS_SB
rows[19,1] = Q_ThermPlate
rows[20,1] = T_CC1_1st
rows[21,1] = T_CC2_1st
rows[22,1] = T_CC1_2nd
rows[23,1] = T_CC2_2nd

# Writing to csv file
with open('Results/' + filename, 'w', newline='') as csvfile:
    # Creating a csv writer object
    csvwriter = csv.writer(csvfile)
    # Writing the data
    csvwriter.writerow(fields)
    csvwriter.writerows(rows)
##############################################################################
"""

"""
##############################################################################
## Calculation Two Cryocoolers Parallel ##

# Initialization of the cooling power of the two CryoCooler stages
Q_CC2_1st = float(Q_CC1_1st) #W
Q_CC2_2nd = float(Q_CC1_2nd) #W
# For the convergence check of the deposited heat
Q_CC2_1st_old = float(Q_CC1_1st_old) #W
Q_CC2_2nd_old = float(Q_CC1_2nd_old) #W

## Loading of an old calculation: 0.6 g/s, Series PT420
p_Cryofan_out = 2000000.0
p_CC1_1st_out = 1999965.9819977405
p_CFHX_HP_out = 1999895.6769761513
p_Plate_out = 1998753.0260555905
p_CM_out = 1998786.6854423531
p_CFHX_LP_out = 1998710.5655103405
p_TS_out = 1994683.8950614152

T_Cryofan_out = 32.87885346721496
T_CC1_1st_out = 25.680754036159442
T_CFHX_HP_out = 6.269904305908413
T_CC1_2nd = 5.283034360592272
T_Plate_out = 5.018369709932516
T_CM_out = 4.890050618552861
T_CFHX_LP_out = 24.935341995913983
T_TS_out = 31.64677087371011

Q_CC1_1st = 16.005104660339025
Q_CC1_2nd = 2.77613079229747


# Convergence criterion
has_converged_outer = False

while not has_converged_outer:
    ##########################################################################
    ## Cryofan ##
    # p_Cryofan_out = 20e5 #Pa
    # T_Cryofan_out = 60.0 #K

    ##########################################################################
    ## CC1 1st stage ##

    # Calculating the constant stage temperature depending on the cooling power of both stages
    T_CC1_1st = T_s1(P1 = Q_CC1_1st, P2 = Q_CC1_2nd)

    ## Heat load of service box TS
    # Heat load for the TS service box with the current temperatures
    # Shifted temperature to be realistic
    # Here a dT of 10 K over the whole service box is assumed
    T_TS_SB_mean = T_CC1_1st + 10./2. #K
    Q_TS_SB = HL_factor_TS_SB * T_func_MLI(T_RT, T_TS_SB_mean, pressure_vac, N_MLI_TS)

    # Using the HEX_CryoCooler function to calculate the stage after the 1st stage depending on state 1
    # Parallel layout so both cryocoolers get treated exactly the same
    # Mass flow and external heat load is halfed
    state_CC1_1st = HEX_CryoCooler(T_Cryofan_out, p_Cryofan_out, 0.5*m_dot, d_i_1st, T_CC1_1st, d_coldhead_1st, l_HEX_1st, 0.5*Q_TS_SB, N)
    # Transfering the results
    p_CC1_1st_out = state_CC1_1st.get("p")
    T_CC1_1st_out = state_CC1_1st.get("T")
    Q_CC1_1st = (1-relax_factor) * Q_CC1_1st + relax_factor * state_CC1_1st.get("Q_cooling")

    ##########################################################################
    ## CC2 1st stage ##

    # Second cryocooler in parallel, thus exactly the same
    T_CC2_1st = T_CC1_1st
    p_CC2_1st_out = p_CC1_1st_out
    T_CC2_1st_out = T_CC1_1st_out
    Q_CC2_1st = Q_CC1_1st

    ##########################################################################
    ## CFHX ##
    # Need for iterativ calculation because the inlet of the CFHX depends on the outlet

    # Convergence criterion
    has_converged_inner = False

    while not has_converged_inner:

        state_CFHX = CFHX(m_dot, p_nominal, epsilon, p_CC1_1st_out, T_CC1_1st_out, p_Plate_out, T_Plate_out)

        # Transfering the results
        p_CFHX_HP_out = state_CFHX.get("p_HP")
        T_CFHX_HP_out = state_CFHX.get("T_HP")

        ######################################################################
        ## CC1 2nd stage ##

        # Parallel layout so both cryocoolers get treated exactly the same
        # Mass flow and external heat load is halfed
        # Calculating the constant stage temperature depending on the cooling power of both stages
        T_CC1_2nd = T_s2(P1 = Q_CC1_1st, P2 = Q_CC1_2nd)
        state_CC1_2nd = HEX_CryoCooler(T_CFHX_HP_out, p_CFHX_HP_out, 0.5*m_dot, d_inner, T_CC1_2nd, d_coldhead_2nd, l_HEX_2nd, 0.5*0, N)
        # Transfering the results
        p_CC1_2nd_out = state_CC1_2nd.get("p")
        T_CC1_2nd_out = state_CC1_2nd.get("T")
        Q_CC1_2nd = (1-relax_factor) * Q_CC1_2nd + relax_factor * state_CC1_2nd.get("Q_cooling")

        ######################################################################
        ## CC2 2nd stage ##

        # Second cryocooler in parallel, thus exactly the same
        T_CC2_2nd = T_CC1_2nd
        p_CC2_2nd_out = p_CC1_2nd_out
        T_CC2_2nd_out = T_CC1_2nd_out
        Q_CC2_2nd = Q_CC1_2nd

        ######################################################################
        ## Remote cooling CM ##

        # Heat load for the CM with the current temperatures
        # Arithmetic mean temperature between inlet and outlet for average temperature
        T_CM_mean = 0.5*(T_CC1_2nd_out + T_CM_out)
        Q_CM = HL_factor_CM * T_func_MLI(T_TS_mean, T_CM_mean, pressure_vac, N_MLI_CM)
        state_CM_out = Remote_Cooling(T_CC1_2nd_out, p_CC1_2nd_out, m_dot, d_inner, l_pipe, N, 0, Q_CM)
        # Transfering the results
        p_CM_out = state_CM_out.get("p")
        T_CM_out = state_CM_out.get("T")

        ######################################################################
        ## Cooling thermalization plate electrical equipment ##

        # Heat load of the thermalization plate
        # Assumed to be constant!
        # Q_CM = HL_factor_CM * T_func_MLI(T_RT, T_TS_mean, pressure_vac, N_MLI_CM)
        # Specific HL of the plate
        # Heat load is fixed and assumed to be constant
        q_plate = Q_ThermPlate/(np.pi * d_inner * l_plate) #W/m²
        state_Plate = SimplePipe(T_CM_out, p_CM_out, m_dot, d_inner, l_plate, N, 0, q_plate)

        # Check convergence criterion
        has_converged_inner = ( abs((state_Plate.get("p") - p_Plate_out)/state_Plate.get("p")) < convergence_criterion_inner
                            and abs((state_Plate.get("T") - T_Plate_out)/state_Plate.get("T")) < convergence_criterion_inner
                            and abs((Q_CC1_2nd-Q_CC1_2nd_old)/Q_CC1_2nd) < convergence_criterion_inner
                            and abs((Q_CC2_2nd-Q_CC2_2nd_old)/Q_CC2_2nd) < convergence_criterion_inner )
        # Transfering the results
        p_Plate_out = state_Plate.get("p")
        T_Plate_out = state_Plate.get("T")
        Q_CC1_2nd_old = Q_CC1_2nd
        Q_CC2_2nd_old = Q_CC2_2nd



    # End of inside loop
    ##########################################################################
    ## CFHX LP ##

    # Transfering the results from the CFHX earlier
    p_CFHX_LP_out = state_CFHX.get("p_LP")
    T_CFHX_LP_out = state_CFHX.get("T_LP")

    ##########################################################################
    ## Remote cooling TS ##

    # Heat load for the TS with the current temperatures
    # Arithmetic mean temperature between inlet and outlet for average temperature
    T_TS_mean = 0.5*(T_CFHX_LP_out + T_TS_out)
    Q_TS = HL_factor_TS * T_func_MLI(T_RT, T_TS_mean, pressure_vac, N_MLI_TS)
    state_TS = Remote_Cooling(T_CFHX_LP_out, p_CFHX_LP_out, m_dot, d_inner, l_pipe, N, 0, Q_TS)
    # Transfering the results
    p_TS_out = state_TS.get("p")
    T_TS_out = state_TS.get("T")

    ##########################################################################
    ## Piping in the SB ##

    # Pressure drop due to all the flow restrictions inside the system that are not included in remote cooling
	# Estimate how many arcs and turns there will be and then put them all together in one point
	# Circulator → CC1 1st: 90°arc, CC1 1st → CC2 1st: 90° arc, CC2 1st → CFHX: 2 90° arc, CFHX → CC1 2nd: 2 90° arc
	# CC1 2nd → CC2 2nd: 90° arc, CC2 2nd → CM: 2 90° arc, CM → Therm. Plate: 90° arc, Therm Plate Maender: 4 180° arc
	# Therm. Plate → CFHX: 2 90° arc, CFHX → TS: 2 90° arc, TS → Circulator: 90° arc
	# SUM: 15 90° arc + 4 180° arc = 23 90° arc = 23 * 1.3
    state_Piping_SB = FlowRestriction(T_TS_out, p_TS_out, m_dot, d_inner, 23*1.3)
    # Presure drop and HL due to thermal radiation of the TS in the service box
    # 2.5 m length and a pipe temperature of the CC1 2nd stage are assumed
    q_pipe_SB = (Q_therm_rad(2.45, (np.pi * d_inner * 2.5), T_TS_SB_mean, T_CC1_2nd_out, 0.05, 0.05)
                /(np.pi * d_inner * 2.5)) #W/m²
    state_Cryofan_in = SimplePipe(state_Piping_SB.get("T"), state_Piping_SB.get("p"), m_dot, d_inner, 2.5, 50, 0, q_pipe_SB)
    p_Cryofan_in = state_Cryofan_in.get("p")
    T_Cryofan_in = state_Cryofan_in.get("T")

    ##########################################################################
    ## Cryofan ##

    state_Cryofan = CryoFan(m_dot, p_Cryofan_in, T_Cryofan_in, p_Cryofan_out, T_Cryofan_out)

    # Convergence check
    has_converged_outer = (   abs((state_Cryofan.get("p") - p_Cryofan_out)/state_Cryofan.get("p")) < convergence_criterion_outer
                          and abs((state_Cryofan.get("T") - T_Cryofan_out)/state_Cryofan.get("T")) < convergence_criterion_outer
                          and abs((Q_CC1_1st-Q_CC1_1st_old)/Q_CC1_1st) < convergence_criterion_outer
                          and abs((Q_CC2_1st-Q_CC2_1st_old)/Q_CC2_1st) < convergence_criterion_outer )
    Q_CC1_1st_old = Q_CC1_1st
    Q_CC2_1st_old = Q_CC2_1st

    # Transfering the results
    p_Cryofan_out = state_Cryofan.get("p")
    T_Cryofan_out = state_Cryofan.get("T")

    print("State after Cryofan: ", state_Cryofan)

# End of the outside loop
##############################################################################
## Postprocessing ##

# Calculation of the volume flow through Cryofan
Rho_Cryofan_in = hp.HeCalc(3,  0, 1, p_Cryofan_in, 2, T_Cryofan_in, 1) #kg/m³
Rho_Cryofan_out = hp.HeCalc(3,  0, 1, p_Cryofan_out, 2, T_Cryofan_out, 1) #kg/m³
Rho_Cryofan = 0.5 * (Rho_Cryofan_in + Rho_Cryofan_out)
vol_flow = m_dot / Rho_Cryofan * 3600 #m³/h
# Calculation of the pressure head in m
dp_system = p_Cryofan_out-p_Cryofan_in #Pa
press_head_system =  dp_system/(Rho_Cryofan*9.81)
print("Pressure head system: ", press_head_system, " m")

# T_in and T_out for CM and TS
print("T_CM_in: ", T_CC1_2nd_out, " K")
print("T_CM_out: ", T_CM_out, " K")
print("T_TS_in: ", T_CFHX_LP_out, " K")
print("T_TS_out: ", T_TS_out, " K")

# Heat loads that are a function of temperature
print("Q_TS: ", Q_TS, " W")
print("Q_CM: ", Q_CM, " W")
print("Q_TS_SB: ", Q_TS_SB, " W")
print("Q_ThermPlate: ", Q_ThermPlate, " W")

##############################################################################
"""

"""
#Save the data
# Data of the file
fields = ['Description', 'Value']
# Restructuring into rows
rows = np.zeros((20,2), dtype = float)
# rows[0,0] = "mass flow"
rows[0,1] = m_dot
# rows[1,0] = "p nom"
rows[1,1] = p_nominal
# rows[1,0] = "p circ in"
rows[2,1] = p_Cryofan_in
# rows[2,0] = "press head"
rows[3,1] = press_head_system
# rows[3,0] = "vol flow"
rows[4,1] = vol_flow
# rows[4,0] = "T CM out"
rows[5,1] = T_CM_out
# rows[5,0] = "T TS out"
rows[6,1] = T_TS_out
# rows[6,0] = "T CC1 1st out"
rows[7,1] = T_CC1_1st_out
# rows[7,0] = "T CC1 2nd out"
rows[8,1] = T_CC1_2nd_out
# rows[8,0] = "T CFHX HP out"
rows[9,1] = T_CFHX_HP_out
# rows[9,0] = "T CFHX LP out"
rows[10,1] = T_CFHX_LP_out
rows[11,1] = T_Cryofan_in
rows[12,1] = T_Cryofan_out
# rows[11,0] = "T IVC out"
rows[13,1] = T_Plate_out
# rows[12,0] = "Q TS"
rows[14,1] = Q_TS
# rows[13,0] = "Q CM"
rows[15,1] = Q_CM
# rows[14,0] = "Q TS SB"
rows[16,1] = Q_TS_SB
# rows[15,0] = "Q_ThermPlate"
rows[17,1] = Q_ThermPlate
# rows[16,0] = "T CC1 1st "
rows[18,1] = T_CC1_1st
# rows[17,0] = "T CC1 2nd"
rows[19,1] = T_CC1_2nd


# Writing to csv file
with open('Results/' + filename, 'w', newline='') as csvfile:
    # Creating a csv writer object
    csvwriter = csv.writer(csvfile)
    # Writing the data
    csvwriter.writerow(fields)
    csvwriter.writerows(rows)
"""
