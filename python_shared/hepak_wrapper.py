# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 13:30:32 2020

@author: bnaydeno
"""


from ctypes import cdll, c_int, c_double, byref, c_long, WinDLL
import os, sys, platform

PYTHON_VERSION = sys.version_info

if PYTHON_VERSION[0] < 3:
    raise Exception("Version 3 of Python is needed")

def getLib():
    if platform.system() == "Windows":
        if PYTHON_VERSION[1] > 7:
            return WinDLL(os.path.join(
                os.path.dirname(__file__),
                'heprop.dll'), winmode = 0x8)
        else:
            return cdll.LoadLibrary(os.path.join(
                os.path.dirname(__file__),
                'heprop.so'))            
    elif platform.system() == "Darwin":
        return cdll.LoadLibrary(os.path.join(
            os.path.dirname(__file__),
            'heprop.dylib'))
    elif platform.system() == "Linux":
        return cdll.LoadLibrary(os.path.join(
            os.path.dirname(__file__),
            'heprop.so'))
    else:
         raise Exception('System unknown. Wrapper supports Windows, Darwin (Mac) and Linux.')

def h_pT(p, T):
    
    t = getLib()

    # OUTPUTS
    
    IDID = c_int(); # If IDID < 0, the calculation failed. Check Ranges. When IDID < -200, the reason for the failure is ambiguous.
    PRP =  ((c_double * 42) * 3)();
    
    # INPUTS
    JOUT = c_int(11000);  # Defined what parameters we want to be calculated. Use 11000
    J1 = c_int(1); # First variable choice
    VALU1 = c_double(p); # First variable value
    J2 = c_int(2); # Second variable choice
    VALU2 = c_double(T); # Second variable value
    JPRECS =  c_int(2); # 1 for modest, 2 for normal and 3 for hight precision.
    JUNITS = c_int(1); # 1 for SI units

    try:
        t.calc_(byref(IDID), byref(PRP), byref(JOUT), byref(J1), byref(VALU1), byref(J2), byref(VALU2), byref(JPRECS), byref(JUNITS))
    except:
        print("Unhandled Exception")
    
    
    if IDID.value < 0:
        raise Exception('The calculation has failed. Check IDID = '+ str(IDID.value) +' in HEPAK\'s documentation.')
    
    
    return PRP[0][9]

# Returns enthalpy for liquid phase in saturated conditions.
def h_p_sL(p):
    
    t = getLib()
    # OUTPUTS
    
    IDID = c_int(); # If IDID < 0, the calculation failed. Check Ranges. When IDID < -200, the reason for the failure is ambiguous.
    PRP =  ((c_double * 42) * 3)();
    
    # INPUTS
    JOUT = c_int(11000);  # Defined what parameters we want to be calculated. Use 11000
    J1 = c_int(1); # First variable choice
    VALU1 = c_double(p); # First variable value
    J2 = c_int(12); # Second variable choice
    VALU2 = c_double(0); # Second variable value
    JPRECS =  c_int(2); # 1 for modest, 2 for normal and 3 for hight precision.
    JUNITS = c_int(1); # 1 for SI units

    try:
        t.calc_(byref(IDID), byref(PRP), byref(JOUT), byref(J1), byref(VALU1), byref(J2), byref(VALU2), byref(JPRECS), byref(JUNITS))
    except:
        print("Unhandled Exception")
        
    if IDID.value < 0:
        raise Exception('The calculation has failed. Check IDID = '+ str(IDID.value) +' in HEPAK\'s documentation.')
    
    return PRP[1][9]


# Returns enthalpy for vapor phase in saturated conditions.
def h_p_sV(p):
    
    t = getLib()
    
    # OUTPUTS
    
    IDID = c_int(); # If IDID < 0, the calculation failed. Check Ranges. When IDID < -200, the reason for the failure is ambiguous.
    PRP =  ((c_double * 42) * 3)();
    
    # INPUTS
    JOUT = c_int(11000);  # Defined what parameters we want to be calculated. Use 11000
    J1 = c_int(1); # First variable choice
    VALU1 = c_double(p); # First variable value
    J2 = c_int(12); # Second variable choice
    VALU2 = c_double(0); # Second variable value
    JPRECS =  c_int(2); # 1 for modest, 2 for normal and 3 for hight precision.
    JUNITS = c_int(1); # 1 for SI units
 
    try:
        t.calc_(byref(IDID), byref(PRP), byref(JOUT), byref(J1), byref(VALU1), byref(J2), byref(VALU2), byref(JPRECS), byref(JUNITS))
    except:
       print("Unhandled Exception")
        
    if IDID.value < 0:
        raise Exception('The calculation has failed. Check IDID = '+ str(IDID.value) +' in HEPAK\'s documentation.')
    
    return PRP[2][9]


def T_ph(p,h):
    
    t = getLib()
    
    # OUTPUTS
    
    IDID = c_int(); # If IDID < 0, the calculation failed. Check Ranges. When IDID < -200, the reason for the failure is ambiguous.
    PRP =  ((c_double * 42) * 3)(); 
    
    # INPUTS
    JOUT = c_int(11000);  # Defined what parameters we want to be calculated. Use 11000
    J1 = c_int(1); # First variable choice
    VALU1 = c_double(p); # First variable value
    J2 = c_int(6); # Second variable choice
    VALU2 = c_double(h); # Second variable value
    JPRECS =  c_int(2); # 1 for modest, 2 for normal and 3 for hight precision.
    JUNITS = c_int(1); # 1 for SI units
    
    try:
        t.calc_(byref(IDID), byref(PRP), byref(JOUT), byref(J1), byref(VALU1), byref(J2), byref(VALU2), byref(JPRECS), byref(JUNITS))
    except:
        print("Unhandled Exception")

    if IDID.value < 0:
        raise Exception('The calculation has failed. Check IDID = '+ str(IDID.value) +' in HEPAK\'s documentation.')
        
    # print(IDID.value, '\n')
    # print(PRP[0][:], '\n')
    # print(PRP[1][:], '\n')
    # print(PRP[2][:], '\n')
    
    return PRP[0][2]


##### getPRP help

############################### INPUTS ########################################

# JOUT: AIM 
# C-------------------------------------- DEFINE PARAMETERS TO BE CALCULATED
# C
# C Output from the calculation is an index IDID and an array PROP,
# C which contains 41 different thermodynamic parameters.  These are
# C described in full detail below.  
# C
# C If you are you are interested in just a few of the 41 parameters and
# C wish to reduce calculation time, the 5-digit integer variable JOUT
# C controls blocks of the calculation, as follows:
# C
# C Pressure, temperature, density, volume, quality and PV/RT are always 
# C calculated, i.e., PROP (0-5) + (if IDID > 1) PROP(6-7).
# C
# C Each non-zero digit of JOUT specifies additional parameters:
# C
# *   If the 10000's digit > 0:   PROP(8-13)  = State Variables
# *   If the 1000's digit > 0 :   PROP(14-23) = Derivatives
# *   If the 100's digit > 0  :   PROP(25-26) = Thermal conductivity and
# C                               Viscosity and (if IDID > 1) PROP(29) =
# C                               surface tension.
# C                               Also, if Derivatives have been calculated,
# C                               PROP(27-28) = Diffusivity & Prandtl No.
# C                               will be calculated.
# *   If the 10's digit > 0   :   PROP(30-31) = Dielectric constant and
# C                               refractive index.   
# *       Note! PROP(40,0) = wavelength (micrometers) is a required input for
# *       the refractive index calculation! 
# C
# *   If the 1's digit > 0    :  PROP(8-23) + PROP(31-39) =  HeII properties.
# C       State and derivative properties will automatically be calculated
# C       when HeII properties are requested.
# C
# C  The following will result in calculation of all possible fluid properties.
# C
#       JOUT = 11111


# J1 and J2
#      JINn         identifies
# C     ---          ----------
# C      1         P  Pressure
# C      2         T  Temperature
# C      3         D  Density
# C      4         V  Specific volume
# C      5         S  Entropy
# C      6         H  Enthalpy
# C      8         G  Gibbs Energy
# C      9         X  Quality
# C     10        dT  delta-Temperature to the lambda line
# C     11         M  on the melting line
# C     12         2  on the saturation line, calculate both liquid and vapor
# C     13        sL  on the saturation line, calculate liquid phase only
# C     14        sV  on the saturation line, calculate vapor phase only
# C     15         L  on the lambda line
# 
#
#   Valid pairs of input variables are indicated in the following array:
#
# C
# *         P   T   D   V   S   H   U   G   X  dT   M   2  sL  sV   L 
# C   JINn  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
# C         -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# C  P  1|      x   x   x   x   x   x   x   x   x   x   x   x   x   x 
# C  T  2|  x       x   x   x           x   x       x   x   x   x   x 
# C  D  3|  x   x           x   x   x           x       x   x   x   x 
# C  V  4|  x   x           x   x   x           x       x   x   x   x 
# C  S  5|  x   x   x   x       x                       x   x   x      
# C  H  6|  x       x   x   x                           x   x  (x)    
# C  U  7|  x       x   x                               x   x  (x)    
# C  G  8|  x   x                                       x   x   x     
# C  X  9|  x   x                                                     
# C dT 10|  x       x   x                               x   x         
# C  M 11|  x   x                                                   x 
# C  2 12|  x   x   x   x   x   x   x   x       x                   x 
# C sL 13|  x   x   x   x   x   x   x   x       x                   x 
# C sV 14|  x   x   x   x   x  (x) (x)  x                             
# C  L 15|  x   x   x   x                           x   x   x 



# JPRECS - PRECISION
# C You may select from three precision criteria:
# C     Modest precision, fast calculations  ........  (1)
# C     Normal precision and speed           ........  (2)
# C     High precision, slow calculations    ........  (3)
    
############################### OUTPUTS #######################################

# PRP - Array structure and Units Selection
#                                          ------ Units selection -------        
# C   i   PROP(i,j)          units set> 1         2         3         4         
# C - -   -------------                                                         
# C * 0:  X = quality                  [-]       [-]       [-]       [-]        
# C         = vapor mass fraction
# C * 1:  P = Pressure                 Pa        kPa       kPa       psi        
# C * 2:  T = Temperature              K         K         K         R          
# C * 3:  D = Density                  kg/m3     kg/m3     kg/m3     lbm/ft3    
# C * 4:  V = Specific volume          m3/kg     cm3/g     cm3/g     ft3/lbm    
# C   5:  PV/RT                        [-]       [-]       [-]       [-]        
# C   6:  dP/dT (saturation line)      Pa/K      kPa/K     kPa/K     psi/R      
# C   7:  Latent Heat                  J/kg      J/g       J/mol     BTU/lb     
# C * 8:  S = Entropy                  J/kg-K    J/g-K     J/mol-K   BTU/lb-R   
# C * 9:  H = Enthalpy                 J/kg      J/g       J/mol     BTU/lb     
# C  10:  A = Helmholtz Energy         J/kg      J/g       J/mol     BTU/lb     
# C *11:  U = Internal Energy          J/kg      J/g       J/mol     BTU/lb     
# C *12:  G = Gibbs Energy             J/kg      J/g       J/mol     BTU/lb     
# C  13:  (reserved for Fugacity)      
# C  14:  Cp                           J/kg-K    J/g-K     J/mol-K   BTU/lb-R   
# C  15:  Cv                           J/kg-K    J/g-K     J/mol-K   BTU/lb-R   
# C  16:  gamma = Cp/Cv                [-]       [-]       [-]       [-]        
# C  17:  expansivity = (T/V)(dV/dT)   [-]       [-]       [-]       [-]        
# C  18:  Gruneisen parameter          [-]       [-]       [-]       [-]        
# C        = (V/Cv)(dP/dT) at constant V
# C  19:  isothermal compressibility   1/Pa      1/kPa     1/kPa     1/psi      
# C       = (1/D) dD/dP
# C  20:  velocity of sound            m/s       m/s       m/s       ft/s       
# C  21:  Joule-Thomson coeff          K/Pa      K/kPa     K/kPa     R/psi      
# C       = dT/dP at constant H 
# C  22:  dP/dD at constant T          Pa-m3/kg  kPa-m3/kg kPa-m3/kg psi-ft3/lb 
# C  23:  dP/dT at constant D          Pa/K      kPa/K     kPa/K     psi/R      
# C  24:  V*dH/dV at constant P        J/kg      J/g       J/mol     BTU/lb     
# C  25:  Viscosity                    Pa-s      uPa-s     uPa-s     lbm/ft-s   
# C  26:  Conductivity                 W/m-K     W/m-K     W/m-K     BTU/hrftR  
# C  27:  Thermal diffusivity          m2/s      cm2/s     m2/s      ft2/s      
# C  28:  Prandtl number               [-]       [-]       [-]       [-]        
# C  29:  Surface tension              N/m       dyne/cm   dyne/cm   lbF/ft     
# C  30:  Dielectric constant          [-]       [-]       [-]       [-]        
# C  31:  Refractive index             [-]       [-]       [-]       [-]        
# C *32 = isochoric dT to the lambda line
# C                                    [K]       [K]       [K]       [R]        
# C *33 = isobaric dT to the lambda line
# C                                    [K]       [K]       [K]       [R]        
# C  34 = superfluid density fraction  [-]       [-]       [-]       [-]        
# C  35 = second sound velocity        m2/s      m2/s      m2/s      ft2/s      
# C  36 = fourth sound velocity        m2/s      m2/s      m2/s      ft2/s      
# C  37 = Gorter-Mellink mutual friction constant       
# C                                    m-s/kg    cm-s/g    m-s/kg    m-s/kg     
# C  38 = SuperFluid Thermal Conductivity function    
# C       (see definition below)       W3/m5-kg  W3/cm5-g  W3/cm5-g  W3/m5-kg   
# C  39 = lambda line temperature (isochoric to the state point)
# C                                    K         K         K         K
# C  40 = wavelength (input)           micrometers (for all units sets)
# C  41 (reserved for internal use, no output data)
    

# IDID:  TEST FOR SUCCESS OR FAILURE
# C
# C First, test for success or failure of the calculation.
# C If IDID < 0, the calculation failed.  The most probable cause is that
# C     the input parameters are out of range.  The magnitude of IDID
# C     specifies the failure in greater detail, as follows:
# C IDID    -------------- meaning -------------------
# C   -11   Invalid input parameter(s), JIN1 and/or JIN2          
# C   -12   Invalid combination of input parameters JIN1 and JIN2      
# C  -101   Input Pressure <= zero.                                    
# C  -102   Input Pressure too high; out of range.                     
# C  -103   Input Temperature < 0.8 K; out of range.                   
# C  -104   Input Temperature > 1500 K; out of range.                  
# C  -105   Input Density <= zero.                                     
# C  -106   Input Density is outside of valid range.                   
# C  -107   Solid phase encountered.                                   
# C  -108   Entropy out of range.                                      
# C  -109   Enthalpy out of range.                                     
# C  -110   Internal Energy out of range.                              
# C  -111   Input Gibbs energy must be in liquid between 0.8 and 3.0 K.
# C  -112   Gibbs energy iteration requires 0.8 <= Temperature <=3.0 K.
# C  -113   Pressure out of range for Gibbs energy iteration.          
# C  -121   Saturation pressure outside of valid range.                
# C  -122   Saturation temperature out of range (0.8 to 5.1953 K)      
# C  -123   Saturation density out of range (0.000888 to 146.16 Kg/m3) 
# C  -124   Saturation entropy out of range (0.00471 to 23.936 J/Kg-K) 
# C  -125   Saturation enthalpy out of valid range for the liquid.     
# C  -127   Saturation internal energy outside of valid range for the liquid.   
# C  -128   Saturation Gibbs energy out of range.                      
# C  -129   Saturation quality must be between 0 and 1 (inclusive).    
# C  -131   Melting pressure out of range (25.328 to 1013.25 bars)     
# C  -132   Melting temperature > 14.0 K, pressure exceeds valid range.
# C  -141   Lambda pressure out of range (.050418 to 30.134 bars).     
# C  -142   Lambda temperature out of range (2.1768 to 1.7673 K).      
# C  -143   Lambda density out of range (146.15 to 179.83 Kg/m3).      
# C  -144   On input, absolute value of delta-T must be < 0.5 K        
# C  -161   Calculated temperature < 0.8 K                             
# C  -162   Calculated temperature > 1500. K                           
# C  -163   Calculated pressure > maximum valid pressure.              
# C  -164   Calculated density <= zero.                                
# C  -170   Input specific volume <= zero.                             
# C  -180   Indeterminant function in compressed liquid.      
# C  -201   Iteration failure with (P,T) input.  Out of range?
# C  -202   Iteration failure with (H,S) input.  Out of range?
# C  -203   Iteration failure with (P,G) input.  Out of range?      
# C  -204   Iteration failure with (D,P) input; compressed liquid?  
# C  -205   Failure with (D,S), (D,H) or (D,U) input.  Out of range?
# C  -206   Failure with (P,S), (P,H) or (P,U) input.  Out of range?
# C  -207   Unexpected iteration failure near the lambda line.  



# C-------- Fluid phase information associated with PROP (i,j) and IDID
# C
# C IDID      ------ location of output variables -------
# C    1      single phase fluid in PROP(i,0) only
# C    2      sat liquid and vapor in PROP(i,1) and PROP(i,2) only; 
# C    3      same as IDID=2 + mixture properties in PROP(i,0)
# C    4      sat liquid in PROP(i,1) only (JINn=13, see below);
# C    5      sat vapor  in PROP(i,2) only (JINn=14, see below);
# C    <0     indicates that an error has occurred, results are invalid.
# C
# C Thus, on output:
# C
# C PROP (i,0) refers to single phase fluid,   when IDID = 1
# C                      liquid-vapor mixtures when IDID = 3
# C            will be zero for other values of IDID
# C PROP (i,1) refers to saturated liquid when IDID = 2, 3, or 4
# C            will be zero for other values of IDID
# C PROP (i,2) refers to saturated vapor when IDID = 2, 3, or 5
# C            will be zero for other values of IDID

# def checkState(IDID, PROP):
#         if IDID == 1: # single-phase
#             return PROP[:][0]
#         elif IDID == 
            
#     return PROP
    

def getPRP(JOUT, J1, VALU1, J2, VALU2, JPRECS=2, JUNITS=1):
    
    t = getLib()

    # OUTPUTS
    
    IDID = c_long(); # If IDID < 0, the calculation failed. Check Ranges. When IDID < -200, the reason for the failure is ambiguous.
    PRP = ((c_double * 42) * 3)(); 
    
    # INPUTS
    JOUT = c_int(JOUT);  # Defined what parameters we want to be calculated. Use 11000
    J1 = c_int(J1); # First variable choice
    VALU1 = c_double(VALU1); # First variable value
    J2 = c_int(J2); # Second variable choice
    VALU2 = c_double(VALU2); # Second variable value
    JPRECS =  c_int(JPRECS); # 1 for modest, 2 for normal and 3 for hight precision.
    JUNITS = c_int(JUNITS); # 1 for SI units


    t.calc_(byref(IDID), byref(PRP), byref(JOUT), byref(J1), byref(VALU1), byref(J2), byref(VALU2), byref(JPRECS), byref(JUNITS))
    
    if IDID.value < 0:
        raise Exception('The calculation has failed. Check IDID = '+ str(IDID.value) +' in HEPAK\'s documentation.')
    
    return PRP, IDID.value


def translateInput(inputValue):
    # HEPAK manual has a mistake in that the Excel function uses numbers from
    # table III-2 as input values instead of the ones from table III-4. Hence,
    # this function translates numbers from table III-4 to table III-2 from 
    # the HEPAK manual providing the same user experience as the Excel function.
    
    # As for the user using HeCalc function, the inputs should be provided
    # according to the following table:

    #         P   T   D   V   S   H   U   G   X  dT   M   2  sL  sV   L 
    #         1   2   3   4   8   9   11  12  0  
                                                                     
    #  P  1|      x   x   x   x   x   x   x   x   x   x   x   x   x   x 
    #  T  2|  x       x   x   x           x   x       x   x   x   x   x 
    #  D  3|  x   x           x   x   x           x       x   x   x   x 
    #  V  4|  x   x           x   x   x           x       x   x   x   x 
    #  S  8|  x   x   x   x       x                       x   x   x      
    #  H  9|  x       x   x   x                           x   x  (x)    
    #  U 11|  x       x   x                               x   x  (x)    
    #  G 12|  x   x                                       x   x   x     
    #  X   |  x   x                                                     
    # dT   |  x       x   x                               x   x         
    #  M   |  x   x                                                   x 
    #  2   |  x   x   x   x   x   x   x   x       x                   x 
    # sL   |  x   x   x   x   x   x   x   x       x                   x 
    # sV   |  x   x   x   x   x  (x) (x)  x                             
    #  L   |  x   x   x   x                           x   x   x 


   #  Input variable 'dT' is interpreted as the isobaric temperature
   #  difference to the lambda line when paired with 1 (="P"), the isochoric
   #  temperature difference to the lambda line when paired with 3 (="D")
   #  or 4 (="V"), and the temperature difference from a point on the 
   #  saturation line to the lambda point when paired with ("2") or
   #  ("sL").  


    translateDict = {
        1 : 1,
        'P': 1,
        2: 2,
        'T': 2,
        3: 3,
        'D': 3,
        4: 4,
        'V': 4,
        8: 5,
        'S': 5,
        9: 6,
        'H': 6,
        11: 7,
        'U': 7,
        12: 8,
        'G': 8,
        0: 9,
        'X': 9,
        'dT': 10,
        'M': 11,
        '2': 12,
        'sL': 13,
        'sV': 14,
        'L': 15
    }
    
    if inputValue not in translateDict:
        raise Exception(inputValue + ' is not a valid input value. Try again.')
    
    return translateDict[inputValue]

# Similar functionality as Excel
def HeCalc(Index, Phase, Input1, Value1,  Input2, Value2, Unit):
    
    # Index is a number between 0 and 39 (inclusive) which denotes the fluid property to be returned
    # Phase is a number between 0 and 5 (inclusive) which allows the user to specify which fluid phase (single-phase, liquid, mixture, or vapor) property is to be returned
    Input1 = translateInput(Input1)
    Input2 = translateInput(Input2)

    t = getLib()

    # OUTPUTS
    
    IDID = c_long(); # If IDID < 0, the calculation failed. Check Ranges. When IDID < -200, the reason for the failure is ambiguous.
    PRP = ((c_double * 42) * 3)(); 
    
    # INPUTS
    JOUT = c_int(11111);  # Defined what parameters we want to be calculated. Use 11000
    J1 = c_int(Input1); # First variable choice
    VALU1 = c_double(Value1); # First variable value
    J2 = c_int(Input2); # Second variable choice
    VALU2 = c_double(Value2); # Second variable value
    JPRECS =  c_int(2); # 1 for modest, 2 for normal and 3 for hight precision.
    JUNITS = c_int(Unit); # 1 for SI units
    


    t.calc_(byref(IDID), byref(PRP), byref(JOUT), byref(J1), byref(VALU1), byref(J2), byref(VALU2), byref(JPRECS), byref(JUNITS))
    IDID = IDID.value
    

    if IDID < 0:
        raise Exception('The calculation has failed. Check IDID = '+ str(IDID) +' in HEPAK\'s documentation.')
    
    if Phase == 0:
        if IDID == 1:
            # single phase fluid in PROP(i,0) only
            return PRP[0][Index]
        elif IDID == 3:
            # mixture properties in PROP(i,0)
            return PRP[0][Index]
        elif IDID == 2:
            # saturated vapor properties 
            if Input1 == 14 or Input2 == 14:
                return PRP[2][Index]
            # saturated liquid properties 
            else:
                return PRP[1][Index]
        elif IDID == 4:
            return PRP[1][Index]
        elif IDID == 5:
            return PRP[2][Index]
            
    elif Phase == 1:
        if IDID == 1:
            # single phase fluid in PROP(i,0) only
            return PRP[0][Index]
        elif IDID == 3:
            # mixture properties in PROP(i,0)
            return PRP[0][Index]
        else:
            return "N/A"     
    elif Phase == 2:
        if IDID == 1:
            # single phase fluid in PROP(i,0) only
            return PRP[0][Index]
        else:
            return "N/A"
    elif Phase == 3:
        if IDID == 3:
            # mixture properties in PROP(i,0)
            return PRP[0][Index]
        else:
            return "N/A"      
    elif Phase == 4:
        # returns the saturated liquid properties 
        if IDID == 2 or IDID == 3 or IDID == 4:
            return PRP[1][Index]
        else:
            return "N/A"
    elif Phase == 5:
        # returns the saturated vapor properties 
        if IDID == 2 or IDID == 3 or IDID == 5:
            return PRP[2][Index]
        else:
            return "N/A"     
    else:
         raise Exception('Invalid Phase value. It should be between 0 and 5.')