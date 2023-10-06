# Celestial Navigation - Position Determination Python Script 
# By Chris Seymour, EIT
# Copyright Ulnooweg Education Centre, 2023, All rights reserved
###############################################################################
import time
import math
import numpy as np
from astropy import units as u
from astropy.units import UnitsError
from astropy.coordinates import Angle, Latitude, Longitude
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
import astropy.time as astime

# Options
decimalDegrees = True
testMode = True
numObs = 3 # might implement more than 3 observations in the future

# Constants
UNIX_epoch_JD = 2440587.5
deg2rad = math.pi / 180.0
rad2deg = 180.0 / math.pi
NavStars = (  # Sideral Hour Angle (SHA) [deg], declination (DEC) [deg]
    (319.0, 89.0),*
    (358.0, 29.0),
    (354.0, -42.0),
    (350.0, 56.0),
    (349.0, -18.0),
    (336.0, -57.0),
    (328.0, 23.0),
    (316.0, -40.0),
    (315.0, 4.0),
    (309.0, 50.0),
    (291.0, 16.0),
    (282.0, -8.0),
    (281.0, 46.0),
    (279.0, 6.0),
    (279.0, 29.0),
    (276.0, -1.0),
    (271.0, 7.0),
    (264.0, -53.0),
    (259.0, -17.0),
    (256.0, -29.0),
    (245.0, 5.0),
    (244.0, 28.0),
    (234.0, -59.0),
    (223.0, -43.0),
    (222.0, -70.0),
    (218.0, -9.0),
    (208.0, 12.0),
    (194.0, 62.0),
    (183.0, 15.0),
    (176.0, -17.0),
    (174.0, -63.0),
    (172.0, -57.0),
    (167.0, 56.0),
    (159.0, -11.0),
    (153.0, 49.0),
    (149.0, -60.0),
    (149.0, -36.0),
    (146.0, 19.0),
    (140.0, -61.0),
    (138.0, -16.0),
    (137.0, 74.0),
    (127.0, 27.0),
    (113.0, -26.0),
    (108.0, -69.0),
    (103.0, -16.0),
    (97.0, -37.0),
    (96.0, 13.0),
    (91.0, 51.0),
    (84.0, -34.0),
    (81.0, 39.0),
    (76.0, -26.0),
    (63.0, 9.0),
    (54.0, -57.0),
    (50.0, 45.0),
    (34.0, 10.0),
    (28.0, -47.0),
    (16.0, -30.0),
    (14.0, 15.0),
    )

# Definitions

rng = np.random.default_rng() # random number generator for testing

def haversine_np(lon1, lat1, lon2, lat2): # haversine function for error calculation
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)

    All args must be of equal length.    

    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2

    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367 * c
    return km

def ERAfromUTC(UTC):
    """
    Calculate the current Earth Rotation Angle (degrees) given the current UTC unix time
    """
    UT1 = UTC  # TODO implement UTC to UT1 calculation
    JD_UT1 = (UT1 / 86400.0) + 2440587.5  # UT1 to julian date
    T_u = (JD_UT1 - 2451545.0)  # from Sideral Time wiki page
    ERA = 360 * (0.7790572732640 + 1.00273781191135448 * T_u)  # [deg]
    return ERA

###############################################################################
# Main loop
###############################################################################

while True:
    
    if testMode:        
        # setup test observation location / frame
        testLat = Latitude(rng.uniform(-90,+90), unit=u.deg)
        testLon = Longitude(rng.uniform(-180,+180), unit=u.deg)
        testHeight = rng.uniform(-10, +10)
        testTime = time.time()
        testPos = EarthLocation.from_geodetic(lat=testLat, lon=testLon, height = testHeight, ellipsoid = 'WGS84')
        testFrame = AltAz( location=testPos, obstime = astime.Time(testTime, format='unix') )
        testTime = time.time()
        testERA = ERAfromUTC(testTime)
        
        # grab 3 non repeting stars that are visible from the test location
        testStars = [None for j in range(numObs)] # list of the stars used
        k = 0
        while k < numObs:
            candStar = rng.integers(low=0, high=57, endpoint=True) # pick a random candidate star
            if candStar in testStars: # if this is a repeat star
                continue # try again
            candStarLon = -Longitude( (NavStars[candStar][0]+testERA)%360, unit=u.deg, wrap_angle = 180*u.deg) # lon of the candidate star's sub-stellar point
            candStarLat = Latitude(NavStars[candStar][1], unit=u.deg) # lat of the candidate star's sub-stellar point
            deltaLon = abs(testLon - candStarLon)
            deltaLat = abs(testLat - candStarLat)
            if math.hypot(deltaLat.degree, deltaLon.degree) >= 90.0: #if the sub stellar point is more than 90 deg away
                continue # try again
            else:
                testStars[k] = candStar # use this star
                k = k+1
                continue
                
        # create empty observation objects
        #       [star number, elevation (Ho), time, GHA]
        obs = [[None for i in range(4)] for j in range(numObs)]
        # populate observations
        for i in range(numObs):
            
            # test star number
            obs[i][0] = testStars[i]
            
            
            # test observation time
            obs[i][2] = testTime
            
            # test GHA
            SHA = NavStars[ testStars[i] ][0] # lookup the star numbers & return SHA
            UT1 = obs[i][2]  # TODO implement UTC to UT1 calculation
            JD_UT1 = (UT1 / 86400.0) + 2440587.5  # UT1 to julian date
            T_u = (JD_UT1 - 2451545.0)  # from Sideral Time wiki page
            ERA = 360 * (0.7790572732640 + 1.00273781191135448 * T_u)  # [deg]
            GHA = (ERA + SHA) % 360.0  # TODO wrap from +180 to -180
            obs[i][3] = GHA
        
            # test elevation observation
            starObject = SkyCoord(ra=-SHA*u.deg, dec=NavStars[testStars[i]][1]*u.deg, frame='icrs')
            starObject = starObject.transform_to(testFrame)
            obs[i][1] = starObject.alt
            
        
    else: #get observations from user
        # obs = [star number, elevation (Ho), time, GHA]
        obs = [[None for i in range(4)] for j in range(numObs)]
        for i in range(numObs):  # for each observation
            # get star number from user & validate input
            while True:
                try:
                    starNum = int(input("Star Number: "))
                    assert starNum in range(len(NavStars))
                except(ValueError, AssertionError):
                        print("Invalid star number, try again")
                        continue
                obs[i][0] = starNum  
                break
            # get elevation from user and validate input
            while True:
                try:
                    Ho = Angle(input("Elevation (d = degrees, m = minutes): "))
                    assert Ho.is_within_bounds("-90d", "+90d")
                except(ValueError, TypeError, AssertionError, UnitsError):
                        print("Invalid elevation, try again. (remember your units!)")
                        continue
                obs[i][1] = Ho
                break
            # Get time
            obs[i][2] = time.time()
            # calculate GHA
            SHA = NavStars[starNum][0]  # lookup the star numbers & return SHA
            UT1 = obs[i][2]  # TODO implement UTC to UT1 calculation
            JD_UT1 = (UT1 / 86400.0) + 2440587.5  # UT1 to julian date
            T_u = (JD_UT1 - 2451545.0)  # from Sideral Time wiki page
            ERA = 360 * (0.7790572732640 + 1.00273781191135448 * T_u)  # [deg]
            GHA = (ERA + SHA) % 360.0  # TODO wrap from +180 to -180
            obs[i][3] = GHA
    
    
    # Setup and solve matrix equation
    A = np.zeros([3, 3])  # initalize A as square matrix
    B = np.zeros([3, 1])  # initalize B as column matrix
    for i in range(numObs):  # populate A and B
        star = obs[i][0]  # which star is used for this observation
        DEC = NavStars[star][1] * deg2rad  # lookup DEC, convert to rads
        Ho = obs[i][1] * deg2rad  # lookup elevation, convert to rads
        GHA = obs[i][3] * deg2rad  # lookup, and convert to rads
        A[i] = [np.cos(DEC) * np.cos(GHA), np.cos(DEC) * np.sin(GHA), np.sin(DEC)]
        B[i] = [np.sin(Ho)]
    try:
        # TODO: check for matrix invertability, e.g. if det(A)<0.1 raise a warning
        if np.linalg.det(A) < 0.1:
            input("Warning: solution is ill-defined. Press Enter to proceed.")
        A_inv = np.linalg.inv(A)  # TODO will transpose work here?
        X = np.dot(A_inv, B)  # solve matrix equation for cartesian position
    except (ValueError):
        input("No solution exists. Press Enter to try again.")
        continue
    # Convert to Lat / Lon
    x, y, z = X[0], X[1], X[2]
    lat = Latitude( np.arctan2(z, np.sqrt(x * x + y * y)), unit = u.rad)
    lon = - Longitude( np.arctan2(y, x), unit = u.rad, wrap_angle = 180*u.deg )
    
    # Print results    
    print("Determined position is:")
    if decimalDegrees:
        print("Latitude:  ", lat.to_string(unit = u.deg, decimal = True, precision = 4 ))
        print("Longitude: ", lon.to_string(unit = u.deg, decimal = True, precision = 4 ))
    else:
        print("Latitude:  ", lat.to_string(unit = u.deg, decimal = False, sep = 'dms', precision = 0 ))
        print("Longitude: ", lon.to_string(unit = u.deg, decimal = False, sep = 'dms', precision = 0 ))
        
    if testMode:
        print("test latitude = ", testLat.to_string(unit=u.deg, decimal=True, precision=4))
        print("test longitude = ", testLon.to_string(unit=u.deg, decimal=True, precision=4))
        posError = haversine_np(lon.degree[0], lat.degree[0], testLon.degree, testLat.degree)
        print("position error = ", posError)
    
    input("Press Enter to continue.\n")
    continue
