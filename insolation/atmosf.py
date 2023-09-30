from insolation import constant
import numpy as np



"""
atmosferic physics functions used by insolation

Javier G. Corripio, meteoexploration.com 2023

p2rho(Pz, Ta, RH)
    Computes density of air  at given pressure, temperature
    and relative humidity

rh2sh(RH, Ta, Pz, ice=0)
    Computes specific humidity from relative humidity at given temperature
    and pressure

wvapsat(TempK,ice=0)
    Computes saturated vapor pressure over water and over ice

z2p(z,P0=1.013250E5,T0=288.15)
    Computes pressure for a given altitude according to 
    US standard atmosfphere

"""

def p2rho(Pz, Ta, RH):
    """
        Computes density of air  at given pressure, temperature and relative humidity

    Parameters
    ----------
          Pz    float, actual pressure (hPa)
          Ta    float, air temperature (K)
          RH    float, relative humidity (%)

    Returns
    -------
        float, density of air (scalar or array)

    Nomenclature
    ------------   
          q_a       rho_v/rho  specific humidity air. Brutsaert 3.2
          e_0       partial pressure water vapour
          e_v_star  vapour pressure at saturation (hpa)
          rho_d     density dry air
          rho_v     density water vapour
          0.622 = 18.016/28.966 ratio of molecular weights of water and dry air

    References
    ----------
        Brutsaert, W.: 1982, Evaporation into the atmosfphere: theory, history,
            and applications, Reidel, Dordrecht.

        Jacobson, M. Z.: 1999, Fundamentals of atmosfpheric Modeling,
        Cambridge University Press, Cambridge.

    Examples
    --------
        from insolation import atmosf
        Pz = atmosf.z2p(2100)
        rho = atmosf.p2rho(Pz,275.,50.)
        print(rho)
        0.9930447027246551


    """

    e_v_star = wvapsat(Ta)                        # e_v* in hPa
    e_0 = e_v_star * RH/100.0                     # actual partial pressure of w. vapour,hPa
    P_d = Pz - e_0                                # partial pressure of dry air
    rho_d = 100. * P_d /(constant.R_d*Ta)         # Brutsaert 3.4 (*100  hPa to Pascal)
    rho_v = 0.622*100. * e_0/(constant.R_d*Ta)    # Brutsaert 3.5
    rho = rho_d + rho_v                           # density in kg/m^3
    return(rho)


def rh2sh(RH, Ta, Pz, ice=0):
    """
        Computes specific humidity from relative humidity at given temperature and pressure

    Parameters
    ----------
          RH    float, relative humidity (%)
          Ta    float, air temperature (K)
          Pz    float, actual pressure (hPa)

    Keywords
    --------
        ice   computes vapour pressure over ice

    Returns
    -------
        float, specific humidity, kg of water per  kg of DRY air

    Nomenclature
    ------------   
          q_a       rho_v/rho  specific humidity air. Brutsaert 3.2
          e_0       partial pressure water vapour
          e_v_star  vapour pressure at saturation (hpa)
          rho_d     density dry air
          rho_v     density water vapour
          0.622 = 18.016/28.966 ratio of molecular weights of water and dry air

    References
    ----------
        Brutsaert, W.: 1982, Evaporation into the atmosfphere: theory, history,
            and applications, Reidel, Dordrecht.

        Jacobson, M. Z.: 1999, Fundamentals of atmosfpheric Modeling,
            Cambridge University Press, Cambridge.

    Examples
    --------
        from insolation import atmosf
        RH = 65.
        Ta = 275.
        z = 2000.
        Pz = atmosf.z2p(z)
        q = atmosf.rh2sh(RH,Ta,Pz)
        print(q)
        0.0027881425515333042


    """

    e_v_star = wvapsat(Ta,ice)                  # Saturation vapour pressure in hPa
    e_0 = e_v_star * RH/100.0                   # actual partial pressure of water vapour, hPa
    P_d = Pz - e_0                              # partial pressure of dry air
    # q = (R_Rdv*e_0)/(p_d+R_Rdv*e_0)              # specific humidity. Jacobson99 (2.27)
    rho_d = 100.*P_d /(constant.R_d * Ta)       # dry air Brutsaert 3.4 (*100  hPa to Pascal)
    rho_v = 0.622*100.*e_0/(constant.R_d * Ta)  # density water vapor Brutsaert 3.5
    rho = rho_d + rho_v                         # density in kg/m^3
    q = rho_v/rho                               # specific humidity air, Brutsaert 3.2
    return(q)


def wvapsat(TempK,ice=0):
    """
    computes saturated vapor pressure over water and over ice

    Parameters
    ----------
        TempK float, temperature in Kelvin

    Keywords
    --------
        ice, set this keyword to compute aturated vapor pressure over ice

    Returns
    -------
        float, saturated vapor pressure in hPa 

    Notes
    -----
        Uses Lowe's polynomials

    References
    ----------
        Lowe, P. R.: 1977, An approximating polynomial for the computation of
            saturation vapor pressure, Journal of Applied Meteorology 16, 100-103.

    Examples
    --------
        from insolation import atmosf
        import numpy as np
        import matplotlib.pyplot as plt
        # Saturated vapor pressure at 0 deg C
        TaK = 273.15
        print(atmosf.wvapsat(TaK))
        6.1032389423162385
        # Saturated vapor pressure in the range of Earth's measured temperatures
        minmax_earth = np.arange(-89.2, 56.7)
        minmax_earthK = minmax_earth + 273.15
        wvaprng = atmosf.wvapsat(minmax_earthK)
        plt.plot(minmax_earth,wvaprng)
        plt.xlabel('Temperature [$^{\circ}$C]')
        plt.ylabel('Saturation Vapour Pressure [hPa]')
        plt.title("Saturation Vapor Pressure in Earth's range of measured temperatures")
        plt.show()

    """

    if (np.min(TempK) < 100.0):
        print('looks like temperature is in Celsius or Farenheit. I need Kelvin.')
        return(np.nan)
    # Lowe's polynomials for vapor pressure
    a0 = 6984.505294
    a1 = -188.9039310
    a2 = 2.133357675
    a3 = -1.288580973e-2
    a4 = 4.393587233e-5
    a5 = -8.023923082e-8
    a6 = 6.136820929e-11
    Temp = TempK
    if (ice == 1):
        Temp = TempK - 273.15
        a0 = 6.109177956
        a1 = 5.03469897e-1
        a2 = 1.886013408e-2
        a3 = 4.176223716e-4
        a4 = 5.824720280e-6
        a5 = 4.838803174e-8
        a6 = 1.838826904e-10
    wvap_sat = a0+Temp*(a1+Temp*(a2+Temp*(a3+Temp*(a4+Temp*(a5+Temp*a6)))))
    return(wvap_sat)



def z2p(z,P0=1.013250E5,T0=288.15):
    """
    Computes pressure for a given altitude according to US standard atmosfphere

    Parameters
    ----------
          z float, height a.s.l. (m)
          P0 float, reference sea level pressure (Pa)
          T0 float, reference sea level temperature (K)

    Returns
    -------  
        float (scaler or array), local pressure at height z m a.s.l. 

    References
    ----------
        US Standard Atmosfphere
        U.S. NOAA: 1976, U.S. standard atmosfphere, 1976, NOAA-S/T# 
        76-1562, U.S. National Oceanic and atmosfpheric Administration, 
        National Aeronautics and Space Administration, 
        United States Air Force, Washington. 227 pp.

    Notes
    -----
        pressure as a function of altitude US standard atmosfphere p. 12 & 3
        stlapse*1000 as in table 4, p. 3  temperature gradient is given in K/km
 
    Examples
    --------       
        from insolation import atmosf
        import numpy as np
        import matplotlib.pyplot as plt
        # Atmosferic pressure at 1450 m a.s.l.
        print(atmosf.z2p(1450))
        850.7871646008141
        # Plot the variation of pressure with altitude from sea level to the height of Everest
        z = np.arange(8848+1)
        P_z = atmosf.z2p(z)
        plt.plot(z, P_z)
        plt.xlabel('Altitude [m]')
        plt.ylabel('Pressure [hPa]')
        plt.title("Atmospheric Pressure from Sea Level to Everest")
        plt.grid(visible=True,linestyle="-.")
        plt.show()

    """

    H1 = (constant.REarth * z) /(constant.REarth + z)  # Geopotential height
    HB = 0.0
    zp = P0*(T0/(T0+constant.stlapse*(H1-HB)))**((constant.earth_G * constant.Md)/(constant.atmosf_R * constant.stlapse*1000))
    zp = zp/100.0  #Pa to hPa
    return(zp)
