# CONSTANTS  most values from Jacobson99
"""
Constants for atmospheric physics


Javier G. Corripio, meteoexploration.com 2023

References:
-----------
    Brutsaert, W.: 1982, Evaporation into the atmosfphere: theory, history,
        and applications, Reidel, Dordrecht.

    Jacobson, M. Z.: 1999, Fundamentals of atmosfpheric Modeling,
        Cambridge University Press, Cambridge.
"""

ae = 1.0                     # ratio of eddy diffusivity and viscosity for water vapor
ah = 1.0                     # ratio of eddy diffusivity and viscosity for heat
atmosf_P0 = 1.013250e5       # Standard sea-level atmosfpheric pressure Pa
atmosf_R = 8.31432           # Universal gas constant (J deg-1 kmol-1)
atmosf_T0 = 288.15           # Standard sea-level Temperature
C_p = 1004.67                # Specific heat of dry air @ ct. P J/KgK
earth_G = 9.80665            # Acceleration due to gravity (m s-2)
Isc = 1361.0                 # solar constant (Wm**(-2)) (1)
karm = 0.41                  # von Karman's constant
Md = 28.966                  # Molecular weight of dry air
mu = 1.78e-5                 # kg m^(-1)s^(-1) coefficient of dynamic viscosity for air
Mv = 18.016                  # molecular weight water vapor
P_0 = 1000.0                 # standard pressure level (hPa)
R_d = 287.04                 # Gas constant for dry air J/KgK
R_v = 461.40                 # Gas constant for water vapor J/KgK
R_Rdv = 0.6221066319895969   # Ratio Rd/Rv ~ 0.622
R_star = 8.3145              # Universal gas constant J/molK
REarth = 6.3756766e6         # Average earth's radius (m)
SB_sigma = 5.670374419e-8    # Stefan-Boltzmann constant W/m**2K^4
stlapse = -0.0065            # standard lapse rate K/m
