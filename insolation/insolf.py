import numpy as np
from numpy import linalg as LA
from insolation import atmosf


"""
Functions to compute insolation over tilted surfaces or complex terrain.

Javier G. Corripio meteoexploration.com 2023

aspect(dem,dlxy,degrees = False):
    Calculates the aspect or orientation of the slope of every grid cell in a digital elevation model (DEM).

cgrad(dem, dlxy, cArea = False):
    Computes a unit vector normal to every grid cell in a digital elevation model.

daylength(latitude, longitude, jd, tmz):
    Computes duration of day light for a given latitude and Julian Day.

declination(jd):
    Computes declination from Julian Date.

doshade(dem, res, sun_v, num_sweeps=1):
    Computes cast shadows over a DEM
    python implementation from openamundsen package based on the f90 function in R insol
    by Florian Hanzer https://github.com/openamundsen/openamundsen

eqtime(jd):
    Computes the equation of time for a given Julian Date.

hourangle(jd,longitude,timezone):
    Hour angle, internal function for solar position.

hillshading(dem,dlxy,sunv):
    Computes the intensity of illumination over a surface (DEM) according to the position of the sun.

insolation(zenith,jd,height,visibility,RH,tempK,O3,alphag):
    Computes direct and diffuse solar irradiance perpendicular to the beam, for a given zenith angle, Julian Day, altitude and atmospheric conditions.

julian_day(y, m, d, h, min = 0, sec = 0):
    Computes Julian Day  (days np.since January 1, 4713 BCE at noon UTC).

normalvector(slope,aspect):
    Calculates a unit vector normal to a surface defined by slope inclination and slope orientation.

slope(dem,dlxy,degrees = False):
    Calculates the inclination of the slope of every grid cell in a digital elevation model (DEM). Zero is an horizontal surface.

sunpos(sunv):
    Returns a matrix of azimuth and zenith angles of the sun given the unit vectors from the observer to the direction of the sun.

sunr(jd):
    Calculates the Earth radius vector.

sunvector(jd,latitude,longitude,timezone):
    Calculates a unit vector in the direction of the sun from the observer position.

""" 


def aspect(dem,dlxy,degrees = False):
    """
    Calculates the aspect or orientation of the slope of every grid cell in a digital elevation model (DEM)

    Parameters
    ----------
        dem,  ndarray (2D, float), Digital Elevation Model with elevation data.
        dlxy, int, resolution of the DEM, pixel size along x and y coordinates.

    Keywords
    --------
        degrees, boolean, if True return degrees instead of radians.

    Returns
    -------
        aspect, 2D ndarray of floats of the same size as the input dem.
    
    Examples
    --------
        # Calculates the aspect of every grid cell of a DEM in the Pyrennes and plot it
        from insolation import insolf
        import rasterio
        import matplotlib.pyplot as plt
        demP = rasterio.open('https://meteoexploration.com/insol/data/dempyrenees.tif')
        dlxy = demP.res[0]
        dem = demP.read(1)
        demaspect = insolf.aspect(dem,dlxy,degrees = True)
        plt.imshow(demaspect, cmap='Greys_r')
        plt.colorbar()
        plt.title("Aspect of every grid cell of a DEM in the Pyrenees")
        plt.show()

    """

    demcgrad = cgrad(dem,dlxy)
    y = demcgrad[1,:,:]
    x = demcgrad[0,:,:]
    aspct = np.arctan2(y,x) + np.pi/2
    aspct[aspct < 0] = aspct[aspct < 0] + 2 * np.pi
    if (degrees):
        aspct = np.degrees(aspct)
    return(aspct)



def cgrad(dem, dlxy, cArea = False):
    """
    Computes a unit vector normal to every grid cell in a digital elevation model.

    Parameters
    ----------
        dem, ndarray (2D, float), Digital Elevation Model with elevation data
        dlxy, int, resolution of the DEM, pixel size along x and y coordinates.
        cArea, boolean, if True returns the surface area of every grid cell instead of the gradient.

    Keywords
    --------
        cArea, boolean, if True returns a 2D matrix with the surface area of every grid cell.

    Returns
    -------
        cellgrunit, ndarray float, a 3D matrix corresponding to the x, y, z coordinates of a unit vector 
        perpendicular to every grid cell. 
        If cArea is True, the result is a 2D matrix with the surface area of every grid cell.

    Examples
    --------
        # Visualize x, y z components of vectors normal to a DEM representing a regular pyramid
        from insolation import insolf
        import numpy as np
        import matplotlib.pyplot as plt
        ncols = 100
        nrows = 100
        # you can play with the height to see how the z-component becomes smaller as the height increases and the pyramid becomes steeper
        height = 500  
        nh = 2*height/ncols
        m = np.zeros([ncols,nrows])
        for i in np.arange(ncols):
            for j in np.arange(nrows):
                m[i,j]=height-max(abs(nh*i-height),abs(nh*j-height)) 

        grdm = insolf.cgrad(m,1)
        xcomponent = grdm[0,:,:]
        ycomponent = grdm[1,:,:]
        zcomponent = grdm[2,:,:]
        titletext = ['xcomponent','ycomponent','zcomponent']
        plt.imshow(m, cmap='Greys_r')
        plt.colorbar()
        plt.contour(m)
        plt.title('DEM of a pyramid')
        plt.show()
        print("press 'q' to continue")
        plt.show()
        for p in np.arange(3):
            # plt.imshow(grdm[p,:,:], cmap='Greys_r',vmin=-1, vmax=1)
            plt.imshow(grdm[p,:,:], cmap='Greys_r')
            plt.colorbar()
            plt.contour(m, cmap='Greys_r')
            plt.title(titletext[p])
            print("press 'q' to continue")
            plt.show()

        # Visualize x, y , z components of vector normal to surface for a DEM of the Pyrenees
        # To make the coordinates agree with solar position convention, vector coordinates are positive eastwards, southwards and upward.
        from insolation import insolf
        import rasterio
        import numpy as np
        import matplotlib.pyplot as plt
        # It would be faster if you download the DEM to your local computer
        demP = rasterio.open('https://meteoexploration.com/insol/data/dempyrenees.tif')
        dlxy = demP.res[0]
        dem = demP.read(1)
        demgrd = insolf.cgrad(dem,dlxy)
        xcomponent = demgrd[0,:,:]
        ycomponent = demgrd[1,:,:]
        zcomponent = demgrd[2,:,:]
        titletext = ['xcomponent','ycomponent','zcomponent']
        plt.imshow(dem, cmap='Greys_r')
        plt.colorbar()
        plt.contour(dem)
        plt.title('DEM of the Pyrenees')
        print("press 'q' to continue")
        plt.show()
        for p in np.arange(3):
            plt.imshow(demgrd[p,:,:], cmap='Greys_r')
            plt.colorbar()
            plt.title(titletext[p])
            print("press 'q' to contiqnue")
            plt.show()


    """

    dem = dem[::-1,...] # flip dem to conform python x, y order to to insol coordinate system convention
    rows, cols = dem.shape
    cellgr = np.zeros((3, rows, cols))
    mm = dem
    md = mm[:-1, 1:]
    mr = mm[1:, :-1]
    mrd = mm[1:, 1:]
    cellgr[1, :-1, :-1] = -0.5 * dlxy * (mm[:-1, :-1] + md - mr - mrd)  
    cellgr[0, :-1, :-1] = 0.5 * dlxy * (mm[:-1, :-1] - md + mr - mrd)
    cellgr[2, :-1, :-1] = dlxy * dlxy 
    cellgr[:,:,-1]=cellgr[:,:,-2]  
    cellgr[:,-1,:]=cellgr[:,-2,:]
    cellgr = np.flip(cellgr, axis=1) # flip results back
    cellArea = LA.norm(cellgr,axis=0)
    cellgrunit = np.divide(cellgr,cellArea)
    if cArea:
        return(cellArea)
    else:
        return(cellgrunit)




def daylength(latitude, longitude, jd, tmz):
    """ 
    Computes duration of day light for a given latitude and Julian Day.

    Parameters
    ----------
        latitude,  float, latitude in degrees and decimal fraction.
        longitude, float, longitude in degrees and decimal fraction.
        jd, float, Julian Day, see help(insolf.julian_day).
        tmz, float, Timezone, west of Greenwich is negative.

    Returns
    -------
        sunrise, sunset, duration of day light in hours, ndarray of floats of 
        shape (3, n) where n is length of input array.

    References
    ----------
        Corripio, J. G. (2003). Vectorial algebra algorithms for calculating 
        terrain parameters from DEMs and solar radiation modelling in mountainous
        terrain. International Journal of Geographical Information Science, 17(1), 1–23.  
        https://doi.org/10.1080/713811744

    Notes
    -----
        It considers sunrise and sunset as the time when the center of the sun passes above or below the horizon. 
        It does not take into account limb, summer time, atmospheric refraction or twilight.
        To add, refraction, dip due to observer altitude and more.

    Examples
    --------
        from insolation import insolf
        import numpy as np
        import matplotlib.pyplot as plt

        # Day length at Greenwich Observatory on the 21st of March
        insolf.daylength(51.4767,0.0,insolf.julian_day(2023,3,21,12),0)
        # array([ 6.10004075, 18.14060504, 12.04056429])

        # Daylength at Reykjavík on the 21st of June
        insolf.julian_day(2023,6,21,12)
        2460117.0
        insolf.daylength(64.12,-21.87, 2460117.0, 0)
        # array([ 3.26596303, 23.70863056, 20.44266753])

        # Daylength for the whole 2023 year at Reykjavík 
        jd2023noon = np.arange(insolf.julian_day(2023,1,1,12),insolf.julian_day(2023,12,31,12)+1)
        doy = jd2023noon - jd2023noon[0] + 1
        day_length = insolf.daylength(64.12, -21.87, jd2023noon, 0)
        plt.plot(doy,day_length[2,:])
        plt.xlabel('Day of the year')
        plt.ylabel('Day length [h]')
        plt.ylim(0,24)
        plt.title('Daylength over the year at Reykjavík')
        plt.show()

    """

    EqTime = eqtime(jd)
    delta = declination(jd)
    tanlatdel = -np.tan(np.radians(latitude)) * np.tan(np.radians(delta))
    tanlatdel = np.where(tanlatdel > 1, 1, tanlatdel)
    tanlatdel = np.where(tanlatdel < -1, -1, tanlatdel)
    omega = np.arccos(tanlatdel)
    daylen = (2*omega)/(2*np.pi/24.)
    stndmeridian = tmz*15.
    deltaLatTime = longitude-stndmeridian
    deltaLatTime = deltaLatTime * 24/360. 
    sunrise = 12*(1-omega/np.pi)-deltaLatTime-EqTime/60. 
    sunset = 12*(1+omega/np.pi)-deltaLatTime-EqTime/60.
    sunrise = np.where(omega == 0, np.nan, sunrise)
    sunset = np.where(omega == 0, np.nan, sunset)
    return(np.array([sunrise,sunset,daylen]))



def declination(jd):
    """
    Computes declination from Julian Date.

    Parameters
    ----------
        jd, float, Julian date
   
    Returns
    -------
        float, declination in degrees and decimal fraction

    References
    ----------
        https://gml.noaa.gov/grad/solcalc/calcdetails.html

        Meeus, J. 1999. Astronomical Algorithms. Willmann-Bell, Richmond, Virginia, USA.

        Reda, I. and Andreas, A. 2003. Solar Position Algorithm for Solar Radiation Applications. 55 pp.; 
        NREL Report No. TP-560-34302, Revised January 2008. https://www.nrel.gov/docs/fy08osti/34302.pdf

    Examples
    --------
        # Current time declination
        from insolation import insolf
        import datetime
        tt = datetime.datetime.now()
        jd = insolf.julian_day(tt.year,tt.month,tt.day,tt.hour+tt.minute/60.+tt.second/3600.)
        decl = insolf.declination(jd)

    """

    T = (jd - 2451545)/36525
    epsilon = (23 + 26/60 + 21.448/3600) - (46.815/3600) * T - (0.00059/3600) * T**2 + (0.001813/3600) * T**3
    L0 = 280.46645 + 36000.76983 * T + 0.0003032 * T**2
    M = 357.5291 + 35999.0503 * T - 0.0001559 * T**2 - 4.8e-07 * T**3
    e = 0.016708617 - 4.2037e-05 * T - 1.236e-07 * T**2
    C = (1.9146 - 0.004817 * T - 1.4e-05 * T**2) * np.sin(np.radians(M)) + \
    (0.019993 - 0.000101 * T) * np.sin(2 * np.radians(M)) + 0.00029 * np.sin(3 * np.radians(M))
    Theta = L0 + C
    v = M + C
    Omega = 125.04452 - 1934.136261 * T + 0.0020708 * T**2 + (T**3)/450000
    lambdad = Theta - 0.00569 - 0.00478 * np.sin(np.radians(Omega))
    delta = np.arcsin(np.sin(np.radians(epsilon)) * np.sin(np.radians(lambdad)))
    return(np.degrees(delta))



def doshade(dem, res, sun_v, num_sweeps=1):
    """ 
    Computes cast shadows over a DEM
    unorthodox way to call the function to avoid unnecessary imports (numba) if not used.
    python implementation from openamundsen package based on the f90 function in R insol
    by Florian Hanzer https://github.com/openamundsen/openamundsen
    

    Parameters
    ----------
        dem, ndrray, Digital elevation Model 2D array 
        res, int, resolution of the DEM (pixel size)
        sun_v, float array, unit vector in the direction of the sun
        num_sweeps

    Returns
    -------
        2D array of the same shade as dem with 1 for illuminated pixels and 0 for pixels in the shade

    References
    ----------
        Corripio, J. G. (2003). Vectorial algebra algorithms for calculating
        terrain parameters from DEMs and solar radiation modelling in mountainous
        terrain. International Journal of Geographical Information Science, 17(1),
        1–23. https://doi.org/10.1080/713811744

    Examples
    --------
        # Calculate and plot the cast shadows of a pyramid with the sun at the northwest and 15 degrees elevation
        from insolation import insolf
        import numpy as np
        import matplotlib.pyplot as plt
        # define the sun vector: northwest at 15 degrees elevation
        sv = insolf.normalvector(75,315)
        # define broadly a regular pyramid
        ncols = 100
        nrows = 100
        height = 50  
        nh = 2*height/ncols
        m = np.zeros([ncols,nrows])
        for i in np.arange(ncols):
            for j in np.arange(nrows):
                m[i,j]=height-max(abs(nh*i-height),abs(nh*j-height)) 

        # place it on a larger flat area
        mm = np.zeros([500,500])
        mm[200:300,200:300] = m
        # calculate and plot shadows 
        sh = insolf.doshade(mm, 1, sv)
        plt.imshow(sh, cmap='Greys_r')
        plt.colorbar()
        plt.contour(mm)
        plt.title('Shadow of a Pyramid with the sun at 315 azimuth, 15 elevation')
        print("press 'q' to continue")
        plt.show()

    """

    from insolation import shadows
    return(shadows.shadows(dem, res, sun_v, num_sweeps))




def eqtime(jd):
    """
    Computes the equation of time for a given Julian Date.

    Parameters
    ----------
        jd, float, Julian date.

    Returns
    -------
        float, equation of time in minutes.

    References
    ----------
        https://gml.noaa.gov/grad/solcalc/calcdetails.html

        Meeus, J. 1999. Astronomical Algorithms. Willmann-Bell, Richmond, Virginia, USA.

        Reda, I. and Andreas, A. 2003. Solar Position Algorithm for Solar Radiation Applications. 55 pp.; 
            NREL Report No. TP-560-34302, Revised January 2008. https://www.nrel.gov/docs/fy08osti/34302.pdf
    
    Examples
    --------
        from insolation import insolf
        import numpy as np
        import matplotlib.pyplot as plt

        # plot the equation of time for 2023 at daily intervals
        jd2023noon = np.arange(insolf.julian_day(2023,1,1,12),insolf.julian_day(2023,12,31,12)+1)
        plt.plot(insolf.eqtime(jd2023noon))
        plt.axhline(y=0.0, linestyle='-', color='k', linewidth=0.5)
        plt.xlabel('Day of the year')
        plt.ylabel('Equation of time [minutes]')
        plt.show()

        # Analema
        plt.plot(insolf.eqtime(jd2023noon), insolf.declination(jd2023noon))
        plt.show()

        # Analema viewed from Greenwich Observatory
        latGwch = 51.4791
        x = 180 + insolf.eqtime(jd2023noon) * 15/60
        y = 90 - latGwch + insolf.declination(jd2023noon)
        # plt.plot(x,y)
        plt.scatter(x,y,s=40, facecolors='none', edgecolors='k')
        plt.xlabel("Azimuth")
        plt.ylabel("Elevation")
        plt.ylim(0,90)
        plt.title('Equation of time for 2023 at daily intervals')
        plt.show()

        # Add the solstices and equinoxes (nearest day, see Meeus ch. 26 for more precision)
        decl = insolf.declination(jd2023noon)
        wintersolstice = np.argmin(decl)
        summersolstice = np.argmax(decl)
        # equinoxes when decl is closest to zero, find spring of=n first half of the year and then autumn (for N hemisphere)
        springequinox = np.argmin(np.abs(decl[0:180]))
        autumnequinox = 180 + np.argmin(np.abs(decl[180:365]))
        nodeseqx = ([springequinox,summersolstice,autumnequinox,wintersolstice])
        plt.scatter(x,y,s=40, facecolors='none', edgecolors='k')
        plt.scatter(x[nodeseqx],y[nodeseqx], s=80, facecolors='c', edgecolors='k')
        plt.xlabel("Azimuth")
        plt.ylabel("Elevation")
        plt.ylim(0,90)
        plt.title('Analema, equinoxes and solstices from Greenwich')
        plt.show()

    """

    jdc = (jd - 2451545.0)/36525.0
    sec = 21.448 - jdc*(46.8150 + jdc*(0.00059 - jdc*(0.001813)))
    e0 = 23.0 + (26.0 + (sec/60.0))/60.0  
    ecc = 0.016708634 - jdc * (0.000042037 + 0.0000001267 * jdc)
    oblcorr = e0 + 0.00256 * np.cos(np.radians(125.04 - 1934.136 * jdc))   
    y = (np.tan(np.radians(oblcorr)/2))**2
    l0 = 280.46646 + jdc * (36000.76983 + jdc*(0.0003032))
    l0 = (l0-360*(l0//360))%360
    rl0 = np.radians(l0)
    gmas = 357.52911 + jdc * (35999.05029 - 0.0001537 * jdc)
    gmas = np.radians(gmas)
    EqTime = y*np.sin(2*rl0)-2.0*ecc*np.sin(gmas)+4.0*ecc*y*np.sin(gmas)*np.cos(2*rl0)-0.5*y**2*np.sin(4*rl0)-1.25*ecc**2*np.sin(2*gmas)
    return(np.degrees(EqTime)*4)


def hourangle(jd,longitude,timezone):
    """
    Hour angle, internal function for solar position.

    Parameters
    ----------
          jd, float, Julian date
          longitude, float, longitude in decimal degrees
          timezone, float, timezone in hours, west of Greenwich is negative

    """

    hour = ((jd-np.floor(jd))*24+12) % 24
    eqtimeval = eqtime(jd)
    stndmeridian = timezone*15                  
    deltalontime = longitude-stndmeridian       
    deltalontime = deltalontime * 24.0/360.0    
    omegar = np.pi*( ( (hour + deltalontime + eqtimeval/60)/12.0 ) - 1.0) 
    return(omegar)



def hillshading(dem,dlxy,sunv):
    """
    Computes the intensity of illumination over a surface (DEM) according to the position of the sun.

    Parameters
    ----------
        dem, ndrray, Digital elevation Model 2D array 
        dlxy, int, resolution of the DEM (pixel size)
        sunv, float array, unit vector in the direction of the sun

    Returns
    -------
        hsh, ndarray of floats of the same dimensions as dem Illumination intensity ranging from 0 to 1

    References
    ----------
        Corripio, J. G. (2003). Vectorial algebra algorithms for calculating
        terrain parameters from DEMs and solar radiation modelling in mountainous
        terrain. International Journal of Geographical Information Science, 17(1),
        1–23. https://doi.org/10.1080/713811744        

    Examples
    --------
        from insolation import insolf
        import rasterio
        import numpy as np
        import matplotlib.pyplot as plt
        demP = rasterio.open('https://meteoexploration.com/insol/data/dempyrenees.tif')
        dlxy = demP.res[0]
        dem = demP.read(1)
        # Sun at azimuth 315, 35 degrees elevation
        sv = insolf.normalvector(55,315)
        grd = insolf.cgrad(dem,dlxy)
        hsh = grd[0,:,:]*sv[0] + grd[1,:,:]*sv[1] + grd[2,:,:]*sv[2]
        # remove negative incidence angles (self shading) 
        hsh = (hsh+np.abs(hsh))/2
        plt.imshow(hsh, cmap='Greys_r',vmin=0, vmax=1)
        plt.colorbar()
        plt.show()
        # Add cast shadows
        sh = insolf.doshade(dem, dlxy, sv)
        plt.imshow(hsh*sh, cmap='Greys_r',vmin=0, vmax=1)
        plt.colorbar()
        plt.show() 

    """

    demcgrad = cgrad(dem,dlxy)
    hsh = demcgrad[0,:,:]*sunv[0] + demcgrad[1,:,:]*sunv[1] + demcgrad[2,:,:]*sunv[2]
    hsh = (hsh + abs(hsh))/2.0      # set to zero self-shading pixels
    return(hsh) 



def insolation(zenith,jd,height,visibility,RH,tempK,O3,alphag):
    """
    Computes direct and diffuse solar irradiance perpendicular to the beam, for a given 
    zenith angle, Julian Day, altitude and atmospheric conditions.

    Parameters
    ----------
          zenith       Zenith angle in degrees
          jd           Julian Day
          height       Altitude above sea level
          visibility   Visibility [km]
          RH           Relative humidity [%]
          tempK        Air temperature [K]
          O3           Ozone thickness [m]
          alphag       Albedo of the surrounding terrain [0 to 1]

    Returns
    -------
        array-like, first element is direct radiation, second element is diffuse radiation.
                    If any input is an array First row is direct radiation and second row diffuse radiation

    References
    ----------
        Bird, R. E. and Hulstrom, R. L. (1981a) Review, evaluation and improvements of direct irradiance models, 
        Trans. ASME J. Solar Energy Eng. 103, 182-192.

        Bird, R. E. and Hulstrom, R. L. (1981b) A simplified clear sky model for direct 
        and diffuse insolation on horizontal surfaces, 
        Technical Report SERI/TR-642-761, Solar Research Institute, Golden, Colorado.
        
        Iqbal, M. (1983) An Introduction to Solar Radiation, Academic Press, Toronto.

    Examples
    --------
        from insolation import insolf
        import datetime
        import rasterio
        import numpy as np
        import matplotlib.pyplot as plt
        from tabulate import tabulate

        insolf.insolation(30,2456007,3200,28,60,278.15,0.02,0.2)
        # Daily insolation in Spring
        # Find nearest hour to sunrise and sunset
        TIMESTAMP = "2023-03-21 12:00:00"
        tt = datetime.datetime.strptime(TIMESTAMP,'%Y-%m-%d %H:%M:%S')
        jd = insolf.julian_day(tt.year,tt.month,tt.day,tt.hour+tt.minute/60.+tt.second/3600.)
        latitude = 42.675
        longitude = 0.033
        timezone = 1
        sr,ss,dl = insolf.daylength(42.675,0.033, jd, timezone)
        sr = np.ceil(sr)
        ss = np.floor(ss)
        sr2ss = np.arange(sr,ss)
        jdrng = insolf.julian_day(tt.year,tt.month,tt.day,sr2ss)
        sunv = insolf.sunvector(jdrng,latitude,longitude,timezone)
        azimuth,zenith = insolf.sunpos(sunv)
        height = 2000
        visibility = 60
        RH = 55
        tempK = 285.0
        O3 = 0.02
        alphag = 0.2
        Idir,Idiff = insolf.insolation(zenith,jdrng,height,visibility,RH,tempK,O3,alphag)
        Hnormal = np.array([0,0,1])
        IdirHoriz = Idir * np.dot(np.transpose(sunv),Hnormal)
        print("Hourly direct normal, direct horizontal and diffuse insolation at 43.7N, 0E, Spring")
        print(tabulate({"Hour": sr2ss,"I_direct_normal": Idir, "I_direct_horizontal": IdirHoriz, "I_diff": Idiff}, headers="keys"))
          Hour    I_direct_normal    I_direct_horizontal    I_diff
        ------  -----------------  ---------------------  --------
             8            487.849                83.0165   58.1299
             9            723.68                253.228    87.8123
            10            829.915               419.907   101.989
            11            884.963               555.478   110.169
            12            912.945               645.265   114.711
            13            922.682               680.859   116.373
            14            916.784               659.104   115.362
            15            893.716               581.95    111.559
            16            846.622               456.52    104.389
            17            756.851               295.762    92.0956
            18            567.277               122.506    68.1983
        # plot the results for the whole day
        dayh = np.arange(0,24,5/60)     # every 5 minutes
        jdrng = insolf.julian_day(tt.year,tt.month,tt.day,dayh)
        sunv = insolf.sunvector(jdrng,latitude,longitude,timezone)
        azimuth,zenith = insolf.sunpos(sunv)
        Idir,Idiff = insolf.insolation(zenith,jdrng,height,visibility,RH,tempK,O3,alphag)
        IdirHoriz = Idir * np.dot(np.transpose(sunv),Hnormal)
        plt.figure(figsize=(10,6))
        plt.plot(dayh,Idir,'r', label='Direct Normal')
        plt.plot(dayh,IdirHoriz,'y', label='Direct Horizontal')
        plt.plot(dayh,Idiff,'c', label='Diffuse')
        plt.legend(loc="upper left")
        plt.grid(visible=True,linestyle="-.")
        plt.xlabel('Time')
        plt.ylabel('Insolation Wm$^{-2}$')
        plt.title("Direct normal, direct horizontal and diffuse insolation at 43.7N, 0E, Spring")
        plt.ylim(0,1380)
        plt.show()

        # Alternatively, loop over every hour
        # Compute insolation on a DEM of the pyrenees:
        demP = rasterio.open('https://meteoexploration.com/insol/data/dempyrenees.tif')
        # demP = rasterio.open('dempyrenees.tif')
        dlxy = demP.res[0]
        dem = demP.read(1)
        I_tot = np.zeros(dem.shape)
        for i in np.arange(len(dayh)):
            sunv = insolf.sunvector(jdrng[i],latitude,longitude,timezone)
            azimuth,zenith = insolf.sunpos(sunv)
            hsh = insolf.hillshading(dem,dlxy,sunv)
            Idir,Idiff = insolf.insolation(zenith,jdrng[i],height,visibility,RH,tempK,O3,alphag)
            # Global insolation in MJ/m^2
            W2MJ = 3600e-6 # Watts/m^2 to MJ/m^2 hourly computation
            I_tot = I_tot + W2MJ*Idir*hsh + W2MJ*Idiff   # The diffuse fraction should consider the skyview factor for a propper calculation... to be implemented

        plt.imshow(I_tot, cmap='plasma')
        plt.colorbar()
        plt.title("Dayly insolation over the Pyrenees in Spring $MJ^{-2}$")
        plt.show()

    """

    Isc = 1361.0            # solar constant (Wm**(-2)) (1)
    zenith = np.where(zenith > 90, 90, zenith)
    theta = np.radians(zenith)
    ssctalb = 0.9  # single scattering albedo (aerosols)(Iqbal, 1983)
    Fc = 0.84      # ratio of forward to total energy scattered (Iqbal, 1983)
    Pz = atmosf.z2p(height)
    Mr = 1.0/(np.cos(theta)+0.15*((93.885-zenith)**(-1.253)))
    Ma = Mr*Pz/1013.25
    #** Use Lowe(1977) Lowe's polynomials for vapor pressure
    wvap_s =  atmosf.wvapsat(tempK)
    #Wprec = 0.493*(RH/100.0)*wvap_s/tempK   #precipitable water in cm Leckner (1978)
    Wprec = 46.5*(RH/100.0)*wvap_s/tempK  #Prata 1996
    rho2 = (1/sunr(jd))**2
    TauR = np.exp((-.09030*(Ma**0.84) )*(1.0+Ma-(Ma**1.01)) )
    TauO = 1.0-( ( 0.1611*(O3*Mr)*(1.0+139.48*(O3*Mr))**(-0.3035) ) - 0.002715*(O3*Mr)*( 1.0+0.044*(O3*Mr)+0.0003*(O3*Mr)**2 )**(-1))
    TauG = np.exp(-0.0127*(Ma**0.26))
    TauW = 1.0-2.4959*(Wprec*Mr)*( (1.0+79.034*(Wprec*Mr))**0.6828 + 6.385*(Wprec*Mr) )**(-1)
    TauA = ( 0.97-1.265*(visibility**(-0.66)) )**(Ma**0.9)   #Machler, 1983
    TauTotal = TauR*TauO*TauG*TauW*TauA   
    In = 0.9751*rho2*Isc*TauTotal
    tauaa = 1.0-(1.0-ssctalb)*(1.0-Ma+Ma**1.06)*(1.0-TauA)
    # tauaa = 1 - (1 - ssctalb) * (1 - Ma + Ma^1.06) * (1 - TauA)
    Idr = 0.79*rho2*Isc*np.cos(theta)*TauO*TauG*TauW*tauaa*0.5*(1.0-TauR)/(1.0-Ma+Ma**(1.02))
    tauas = (TauA)/tauaa
    Ida = 0.79*rho2*Isc*np.cos(theta)*TauO*TauG*TauW*tauaa*Fc*(1.0-tauas)/(1.0-Ma+Ma**1.02)
    alpha_atmos = 0.0685+(1.0-Fc)*(1.0-tauas)
    Idm = (In*np.cos(theta)+Idr+Ida)*alphag*alpha_atmos/(1.0-alphag*alpha_atmos)
    Id = Idr+Ida+Idm
    In = np.where(zenith >= 90, 0, In)   # Set In=0 if after sunset
    Id = np.where(zenith >= 90, 0, Id)
    return(np.array([In,Id]))





def julian_day(Y, m, d, H, min = 0, sec = 0):
   """
   Computes Julian Day (days since January 1, 4713 BCE at noon UTC)
   from 1990 edition of the U.S. Naval Observatory's Almanac for Computers Valid AD 1801 - 2099



   Parameters
   ----------
        Y (int)        Year, format yyyy
        m (short int)  Month number.
        d (short int)  Day-of-month.
        H (float)      Hour-of-day and decimal fraction.

   Returns
   -------
        float, Julian Day (days since January 1, 4713 BCE at noon UTC)

   References
   ----------

        https://adsabs.harvard.edu/full/1983IAPPP..13...16F
        https://aa.usno.navy.mil/software/novaspy_intro
        https://aa.usno.navy.mil/faq/JD_formula
        check results 
        https://aa.usno.navy.mil/data/JulianDate

   Examples
   --------
        from insolation import insolf
        import datetime
        TIMESTAMP = "2023-09-23 06:00:00"
        tt = datetime.datetime.strptime(TIMESTAMP,'%Y-%m-%d %H:%M:%S')
        print(insolf.julian_day(tt.year,tt.month,tt.day,tt.hour+tt.minute/60.+tt.second/3600.))
        2460210.75
        insolf.julian_day(2023,21,3,12)  # Check the month/day order, it will accept months > 12 without warning
        2460556.0
        insolf.julian_day(2024,9,3,12)
        2460557.0

   """

   H = H + min/60.0 + sec/3600.0
   tjd = 367*Y - (7*(Y+(m+9)//12))//4 + (275*m)//9 + d + 1721013.5 + H/24.0 
   return(tjd)



def normalvector(slope,aspect):
    """
    Calculates a unit vector normal to a surface defined by slope inclination and slope orientation.

    Parameters
    ----------
        slope (float)       slope inclination in degrees
        aspect (float)      slope orientation in degrees

    Returns
    -------
        float vector
    References
    ----------
        Corripio, J. G. (2003). Vectorial algebra algorithms for calculating
        terrain parameters from DEMs and solar radiation modelling in mountainous
        terrain. International Journal of Geographical Information Science, 17(1),
        1–23. https://doi.org/10.1080/713811744 

    Examples
    --------
        from insolation import insolf
        import numpy as np

        # horizontal surface
        insolf.normalvector(0,0)
        # array([ 0., -0.,  1.])

        # surface 45 degrees south
        insolf.normalvector(45,180)
        # array([8.65956056e-17, 7.07106781e-01, 7.07106781e-01])
        
    """

    sloper = np.radians(slope) 
    aspectr = np.radians(aspect)
    nvx = np.sin(aspectr)*np.sin(sloper)
    nvy = -np.cos(aspectr)*np.sin(sloper)
    nvz = np.cos(sloper)
    return(np.array([nvx,nvy,nvz]))
   


def slope(dem,dlxy,degrees = False):
    """
    Calculates the inclination of the slope of every grid cell in a digital elevation model (DEM). Zero is an horizontal surface.

    Parameters
    ----------
        dem,  ndarray (2D, float), Digital Elevation Model with elevation data.
        dlxy, int, resolution of the DEM, pixel size along x and y coordinates.

    Keywords
    --------
        degrees, boolean, if True return degrees instead of radians.

    Returns
    -------
        slp, 2D ndarray of floats of the same size as the input dem.
   
    Examples
    --------
        # Calculates the slope of every grid cell of a DEM in the Pyrennes and plot it
        from insolation import insolf
        import rasterio
        import matplotlib.pyplot as plt
        demP = rasterio.open('https://meteoexploration.com/insol/data/dempyrenees.tif')
        dlxy = demP.res[0]
        dem = demP.read(1)
        demslope = insolf.slope(dem,dlxy,degrees = True)
        plt.imshow(demslope, cmap='Greys')
        plt.colorbar()
        plt.show()

    """

    demcgrad = cgrad(dem,dlxy)
    z = demcgrad[2,:,:]
    slp = np.arccos(z)
    if (degrees):
        slp = np.degrees(slp)
    slp = slp
    return(slp)


def sunpos(sunv, degrees=True):
    """
    Returns a matrix of azimuth and zenith angles of the sun given the unit vectors from the observer to the direction of the sun.

    Parameters
    ----------
        sunv (float or array of floats)    coordinates x, y, z of the unit vector in the direction of the sun

    Returns
    -------
        array of azimuth and zenith angles

    References
    ----------
        Corripio, J. G. (2003). Vectorial algebra algorithms for calculating
        terrain parameters from DEMs and solar radiation modelling in mountainous
        terrain. International Journal of Geographical Information Science, 17(1),
        1–23. https://doi.org/10.1080/713811744 

    Examples
    --------
        from insolation import insolf
        import numpy as np
        import matplotlib.pyplot as plt

        # Julian Day hourly intervals at spring equinox
        jdrng = insolf.julian_day(2013,3,21,np.arange(7, 20))

        # sun path over the day at sprint equinox in the Pyrenees
        latitude = 42.675
        longitude = 0.033
        tmz = 1
        azimuth, zenith = insolf.sunpos(insolf.sunvector(jdrng,latitude,longitude,tmz))
        elevation = 90-zenith
        plt.scatter(azimuth,elevation,s=240, facecolors='y', edgecolors='y')
        plt.show()

    """

    azimuth = np.pi - np.arctan2(sunv[0],sunv[1] )
    zenith = np.arccos(sunv[2])
    if degrees:
        azimuth = np.degrees(azimuth)
        zenith = np.degrees(zenith)

    return(np.array([azimuth,zenith]))


def sunr(jd):
    """
    Calculates the Earth radius vector.

    Parameters
    ----------
        jd   Julian Day.

    Returns
    -------
        float Earth Radius Vector in Astronomical Units (AU).


    References
    ----------
        https://gml.noaa.gov/grad/solcalc/calcdetails.html

        Meeus, J. 1999. Astronomical Algorithms. Willmann-Bell, Richmond, Virginia, USA.

        Reda, I. and Andreas, A. 2003. Solar Position Algorithm for Solar Radiation Applications. 55 pp.; 
        NREL Report No. TP-560-34302, Revised January 2008. https://www.nrel.gov/docs/fy08osti/34302.pdf 

    Examples
    --------
        # Plot the variation of the earth radius vector over two years
        from insolation import insolf
        import numpy as np
        import matplotlib.pyplot as plt
        ndays = np.arange(730)
        jd_nexty = insolf.julian_day(2013,1,11,12) + ndays
        sun_r = insolf.sunr(jd_nexty)
        plt.plot(ndays,sun_r)
        plt.axhline(y=1.0, linestyle='-', color='k', linewidth=0.5)
        plt.show()

    """

    jdc=(jd - 2451545.0)/36525.0
    ecc=0.016708634-jdc*(0.000042037+0.0000001267*jdc)
    gmas=357.52911+jdc*(35999.05029-0.0001537*jdc)
    gmasr=np.radians(gmas)
    seqc=np.sin(gmasr)*(1.914602-jdc*(0.004817+0.000014*jdc))+np.sin(2*gmas)*(0.019993-0.000101*jdc)+ np.sin(3*gmasr)*0.000289
    sta=gmas+seqc
    sunrv=(1.000001018 * (1 - ecc**2)) / (1 + ecc * np.cos(np.radians(sta)))
    return(sunrv)



def sunvector(jd,latitude,longitude,timezone):
    """
    Calculates a unit vector in the direction of the sun from the observer position.

    Parameters
    ----------
        jd (float)          Julian date
        latitude (float)    Latitude of observer in degrees and decimal fraction
        longitude (float)   Longitude of observer in degrees and decimal fraction
        timezone (float)    Timezone in hours, west of Greenwich is negative

    Returns
    -------
        float vector or array of vectors.

    References
    ----------
        Corripio, J. G. (2003). Vectorial algebra algorithms for calculating
        terrain parameters from DEMs and solar radiation modelling in mountainous
        terrain. International Journal of Geographical Information Science, 17(1),
        1–23. https://doi.org/10.1080/713811744 

    Examples
    --------
        from insolation import insolf
        from mpl_toolkits import mplot3d
        import numpy as np
        import matplotlib.pyplot as plt
        # Current solar vector at Greenwich Observatory on the 21st of June, 2023 at midday
        latitude = 51.4778
        longitude = -0.0017
        tmz = 0
        jd = insolf.julian_day(2023,6,21,12)
        insolf.sunvector(jd,latitude,longitude,0)
        # Path of the sun over Greenwich in summer, z axis in degrees of elevation
        jdrng = insolf.julian_day(2023,6,21,np.arange(4,20.5,0.5))
        sv = insolf.sunvector(jdrng,latitude,longitude,0)
        fig = plt.figure()
        ax = plt.axes(projection ='3d')
        x = sv[0,:]
        y = sv[1,:]
        z = 90 - np.degrees(np.arccos(sv[2,:]))
        ax.plot3D(x, y, z, 'orange')
        ax.set_title('Sun position from Greenwich on the 21st of June')
        plt.show()

    """

    omegar = hourangle(jd,longitude,timezone)
    deltar = np.radians(declination(jd))
    lambdar = np.radians(latitude) 
    svx = -np.sin(omegar)*np.cos(deltar)
    svy = np.sin(lambdar)*np.cos(omegar)*np.cos(deltar)-np.cos(lambdar)*np.sin(deltar)
    svz = np.cos(lambdar)*np.cos(omegar)*np.cos(deltar)+np.sin(lambdar)*np.sin(deltar) 
    return(np.array([svx,svy,svz]))


