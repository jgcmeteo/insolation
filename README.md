
# insolation

A collection of python functions to compute insolation over inclined surfaces or complex terrain.

![insolation](https://www.meteoexploration.com/R/insol/figures/hillshadingpyrenees.jpg)


## Documentation

Full documentation at: https://www.meteoexploration.com/insol/python/index.html

## Overview

Included functions are:

`aspect(dem, dlxy, degrees = False)`  
&ensp;     Calculates the aspect or orientation of the slope of every grid cell in a digital elevation model (DEM).

`cgrad(dem, dlxy, cArea = False)`  
&ensp;   Computes a unit vector normal to every grid cell in a digital elevation model.

`daylength(latitude, longitude, jd, tmz)`  
&ensp;    Compute duration of day light for a given latitude and Julian Day.

`declination(jd)`  
    Computes declination from Julian Date.

`doshade(dem, res, sun_v, num_sweeps=1)`  
&ensp;    Computes cast shadows over a DEM python implementation from openamundsen package based on the f90 function in R insol by Florian Hanzer https://github.com/openamundsen/openamundsen

`eqtime(jd)`  
 &ensp;   Computes the equation of time for a given Julian Date.

`hourangle(jd,longitude, timezone)`  
 &ensp;   Hour angle, internal function for solar position.

`hillshading(dem, dlxy, sunv)`  
 &ensp;   Computes the intensity of illumination over a surface (DEM) according to the position of the sun.

`insolation(zenith, jd, height, visibility, RH, tempK, O3, alphag)`  
 &ensp;   Computes direct and diffuse solar irradiance perpendicular to the beam, for a given zenith angle, Julian Day, altitude and atmospheric conditions.

`julian_day(y, m, d, h, min = 0, sec = 0)`  
 &ensp;   Computes Julian Day  (days np.since January 1, 4713 BCE at noon UTC).

`normalvector(slope, aspect)`  
 &ensp;   Calculates a unit vector normal to a surface defined by slope inclination and slope orientation.

`slope(dem,dlxy, degrees = False)`  
&ensp;    Calculates the inclination of the slope of every grid cell in a digital elevation model (DEM). Zero is an horizontal surface.

`sunpos(sunv)`  
&ensp;    Returns a matrix of azimuth and zenith angles of the sun given the unit vectors from the observer to the direction of the sun.

`sunr(jd)`  
&ensp;   Calculates the Earth radius vector.

`sunvector(jd, latitude, longitude, timezone)`  
&ensp;    Calculates a unit vector in the direction of the sun from the observer position.


Included functions in **atmosf** are:


`p2rho(Pz, Ta, RH)`  
&ensp;    Computes density of air  at given pressure, temperature and relative humidity

`rh2sh(RH, Ta, Pz, ice=0)`  
&ensp;    Computes specific humidity from relative humidity at given temperature and pressure

`wvapsat(TempK,ice=0)`  
&ensp;    Computes saturated vapor pressure over water and over ice

`z2p(z,P0=1.013250E5,T0=288.15)`  
&ensp;    Computes pressure for a given altitude according to  US standard atmosfphere


## References

- Bird, R. E. and Hulstrom, R. L. (1981a) Review, evaluation and improvements of direct irradiance models, 
        Trans. ASME J. Solar Energy Eng. 103, 182-192.

- Bird, R. E. and Hulstrom, R. L. (1981b) A simplified clear sky model for direct 
        and diffuse insolation on horizontal surfaces, 
        Technical Report SERI/TR-642-761, Solar Research Institute, Golden, Colorado.

- Bourges, B.: 1985, Improvement in solar declination computation, Solar Energy 35(4), 367-369. 

- Brutsaert, W.: 1982, Evaporation into the atmosfphere: theory, history,
            and applications, Reidel, Dordrecht.

- Corripio, J. G. (2003). Vectorial algebra algorithms for calculating
        terrain parameters from DEMs and solar radiation modelling in mountainous
        terrain. International Journal of Geographical Information Science, 17(1),
        1–23. https://doi.org/10.1080/713811744 

- Danby, J. M. Eqn. 6.16.4 in Fundamentals of Celestial Mechanics, 2nd ed. Richmond, VA: Willmann-Bell, p. 207, 1988. 

- https://aa.usno.navy.mil/data/JulianDate 

- https://aa.usno.navy.mil/faq/JD_formula

- https://aa.usno.navy.mil/software/novaspy_intro

- https://adsabs.harvard.edu/full/1983IAPPP..13...16F

- https://gml.noaa.gov/grad/solcalc/calcdetails.html


- Iqbal, M. (1983) An Introduction to Solar Radiation, Academic Press, Toronto.

- Jacobson, M. Z.: 1999, Fundamentals of atmosfpheric Modeling,
        Cambridge University Press, Cambridge.

- Lowe, P. R.: 1977, An approximating polynomial for the computation of
        saturation vapor pressure, Journal of Applied Meteorology 16, 100-103.

- Meeus, J. 1999. Astronomical Algorithms. Willmann-Bell, Richmond, Virginia, USA.

- Reda, I. and Andreas, A. 2003. Solar Position Algorithm for Solar Radiation Applications. 55 pp.; 
        NREL Report No. TP-560-34302, Revised January 2008. https://www.nrel.gov/docs/fy08osti/34302.pdf

- US Standard Atmosfphere
        U.S. NOAA: 1976, U.S. standard atmosfphere, 1976, NOAA-S/T# 76-1562, U.S. National Oceanic and atmosfpheric Administration, 
        National Aeronautics and Space Administration, 
        United States Air Force, Washington. 227 pp.




## Requirements

- python 3 or higher
- numpy
- rasterio
- matplotlib




## Installation


pip install insolation

