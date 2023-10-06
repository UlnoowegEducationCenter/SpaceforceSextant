# Celestial Navigation Position Deterimination (MLAT) - Web Version
# By Chris Seymour, EIT
# Copyright Ulnooweg Education Centre, 2023, All rights reserved
###############################################################################
import numpy as np
from numpy import genfromtxt
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import Angle, Latitude, Longitude, EarthLocation
from astropy.coordinates import AltAz, SkyCoord, golden_spiral_grid
from prettytable import PrettyTable


#### Options
outputDMS = False # default output is DD.dddd
testMode  = True  # runs an automated test instead of normal operation
numObs    = 4


#### Constants

# list of the navigational stars' RA and Dec.
NavStars = genfromtxt( 'NavigationalStars.csv', delimiter=',' ) # read the star list csv (from SIMBAD)
NavStars = np.delete( NavStars, 0, axis=0) # delete the label row
NavStars = np.delete( NavStars, (0,1), axis=1) # delete the name collumns
# ( (RA, DEC), ... ) row = almanac number w/ polaris = 0
NavStars = tuple(map(tuple,NavStars)) #convert from numpy array to tuple

# a standard atmosphere for refraction calculation purposes
atm_std = {
    "pressure"     : 101.325*u.kPa,
    "temperature"  : 15*u.deg_C,
    "RH"           : 60*u.pct,
    "wavelength"   : 550*u.nm
    }

# an atmosphere object to pass if you don't want to simulate refraction
atm_No ={
    "pressure"     : 0,
    "temperature"  : None,
    "RH"           : None,
    "wavelength"   : None
    }

# refraction correction factor table (0-5000' ASL), taken from 2023 nautical almanac
Ro =  tuple(Angle([ #sextant reading, refraction correction
       ['+90d00m', '00m'],
       ['+63d00m', '01m'],
       ['+33d00m', '02m'],
       ['+21d00m', '03m'],
       ['+16d00m', '04m'],
       ['+12d00m', '05m'],
       ['+10d00m', '06m'],
       ['+08d10m', '07m'],
       ['+06d50m', '08m'],
       ['+06d00m', '09m'],
       ['+05d20m', '10m'],
       ['+04d30m', '12m'],
       ['+03d30m', '14m'],
       ['+02d50m', '16m'],
       ['+02d20m', '18m'],
       ['+01d50m', '20m'],
       ['+01d12m', '25m'],
       ['+00d34m', '30m'],
       ['+00d06m', '35m'],
       ['-00d18m', '40m']
       ]))

f =  tuple([ #temperature, temperature correction factor
       [+47.0*u.deg_C, 0.9*u.dimensionless_unscaled],
       [+26.0*u.deg_C, 1.0*u.dimensionless_unscaled],
       [+05.0*u.deg_C, 1.1*u.dimensionless_unscaled],
       [-16.0*u.deg_C, 1.2*u.dimensionless_unscaled],
       [-37.0*u.deg_C, 1.3*u.dimensionless_unscaled],
       ])


#### Functions 

# get observations from usert
def getObservationsFromUser(numObs):
    # [star number (int), elevation (latitude), time (datetime)]
    userObs = [ [None, None, None] for i in range(numObs) ]
    for i in range(numObs):
        while True:
            print('-'*16)
            try:
                starNum = input("Star Number: ")
                starNum = int(starNum)
                assert starNum not in [userObs[i][0] for i in range(numObs)] # star not reused
                assert starNum in range(len(NavStars)) # is an actual nav star
            except (ValueError, AssertionError): # excepting errors from invalid star numbers
                print('Invalid star number, try again')
                continue
            userObs[i][0] = starNum
            break
        while True:
            try:
                obsEl = input("Elevation:   ")
                obsEl = Latitude( obsEl, unit = u.deg )
            except (ValueError, TypeError): # excepting errors from invalid elevations
                print('Invalid elevation, try again')
                continue
            userObs[i][1] = obsEl
            userObs[i][2] = Time.now()
            break
        i = i+1 # iterate to get next observation
    return userObs

# create fake observations given a test location and time
def getTestObservations(testLocation, testTime, numObs = 3, round2minutes = True, atm = atm_std, minAlt = 10*u.deg ):
    # [star number (int), elevation (latitude), time (datetime)]
    testObs = [ [None, None, None] for i in range(numObs) ]
    
    #create the test observation reference frame
    testFrame = AltAz(location          = testLocation, 
                      obstime           = testTime, 
                      pressure          = atm["pressure"], 
                      temperature       = atm["temperature"], 
                      relative_humidity = atm["RH"],
                      obswl             = atm["wavelength"]
                      )
    
    
    for i in range(numObs): #for each observation at this location
        while True: #pick a non-reused star  
            try:
                starNum = np.random.randint(low=0, high=len(NavStars))
                assert starNum not in [testObs[i][0] for i in range(numObs)] # star not reused
            except (AssertionError): # if the star is reused
                continue
            testObs[i][0] = starNum
            ra    = NavStars[starNum][0] # read right ascention of selected star
            dec   = NavStars[starNum][1] # read declination of selected star
            star  = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs') #create a skycoord of the selected star
            star  = star.transform_to(testFrame) #transform to the test reference frame
            obsEl = star.alt #extract the elevation (altitude) from the transformed skycoord
            if obsEl < minAlt: # if this star is below 10 degrees elevation
                continue # pick another star and try again
            else:
                break # use this star
        #print(round2minutes)
        #input(obsEl)
        if round2minutes == True: #account for precision of measurement
            obsEl = Latitude(obsEl.to_string(unit=u.deg, fields = 2),unit=u.deg)
            # obsEl = Latitude(obsEl.to_string(unit=u.deg, fields = 1),unit=u.deg)
        #input(obsEl)
        testObs[i][1] = obsEl
        testObs[i][2] = testTime
    return testObs

# estimate location based on multilateration
def multilateratePosition(obs):
    numObs = len(obs)
    assert numObs >= 3
    #obs = [star number (int), elevation (latitude), time (datetime)]
    A = np.zeros([numObs,3]) # initalize the A matrix
    for i in range(numObs): # populate the A matrix
        
        # unpack observation
        starNum = obs[i][0]
        Ho      = obs[i][1]
        obsTime = obs[i][2]
        
        # read star's declination
        dec = Latitude( NavStars[starNum][1], unit = u.deg )
        
        #calculate Greenwich Hour Angle (GHA) of star
        ERA = Time.earth_rotation_angle( obsTime, longitude = 0.0) # earth rotation angle at obstime, also called GHA_star in some books
        SHA = Longitude( 360 - NavStars[starNum][0], unit = u.deg ) # sideral hour angle from right ascention (from csv)
        GHA = Longitude( ERA + SHA, unit = u.deg ) # greenwich hour angle 
    
        A[i] = [np.cos(dec)*np.cos(GHA), np.cos(dec)*np.sin(GHA), np.sin(dec)] / np.sin(Ho)
    
    B = np.ones( (numObs,1) ) # create B vector
    
    # convert to least squares problem
    pseudoA = np.matmul( np.transpose(A), A)
    pseudoB = np.matmul( np.transpose(A), B)
    
    # solve
    X_star = np.linalg.solve(pseudoA,pseudoB)*u.km # solve Ax=B

    x, y, z = X_star[0], X_star[1], X_star[2] # unpack results
    estLat  =  Latitude( np.arctan2( z, np.hypot(x,y) )[0], unit = u.deg )
    estLon  = -Longitude( np.arctan2( y, x )[0], unit = u.deg, wrap_angle = 180*u.deg )
    
    obsLocation = EarthLocation.from_geodetic(lat = estLat, lon = estLon, height = 0.0, ellipsoid='WGS84')
    return obsLocation

# haversine function for error calculation
def haversine_np(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in degrees)

    All args must be of equal length.    

    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2

    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367 * c * u.km / u.rad
    return km

# refraction correction lookup
def refractionCorrection(obs, temp=atm_std["temperature"]):
    correctedObs = obs
    
    # find the temp factor that will be used for all observations
    tempFactor = 1.0; #initalizing temp correction factor as the default 1.0
    assert f[-1][0] <= temp <= f[0][0] #check temp is in the range we have data for
    for j in range(len(f)): # iterate though j
        if temp >= f[j][0]:
            tempFactor = f[j-1][1]
            break
        else:
            j = j+1
            continue
    
    for i in range(len(obs)): # for each observation, find and apply the correction
        correction = Angle('0.0d') # initalize correction to the default value of 0 
        assert Ro[-1][0] <= obs[i][1] <= Ro[0][0] # check altitude is in the range we have data for    
        for k in range(len(Ro)): # iterate through Ro
            if obs[i][1] >= Ro[k][0]: # if the alt is higher than the Ro entry
                correction = Ro[k-1][1]*tempFactor
                correctedObs[i][1] = obs[i][1] - correction # add the previous row's correction times the temp correction factor
                break #exit for loop after applying the correction
            else: # continue iterating
                j = j+1
                continue
    
    return correctedObs

# an alias for the numpy random number generator
rng = np.random.default_rng()

###############################################################################
# Main loop
###############################################################################


while testMode == False: # normal operating mode
    
    # Get observations
    obs = getObservationsFromUser(numObs)
    
    # Correct for refraction
    obs = refractionCorrection(obs)
    
    # multilaterate position
    pos = multilateratePosition(obs)
    
    # Output
    print('-'*16)
    if outputDMS:
        Lat_str = pos.lat.to_string( unit=u.deg, decimal=False, sep='dms', precision=0)
        Lon_str = pos.lon.to_string( unit=u.deg, decimal=False, sep='dms', precision=0)    
    else:
        Lat_str = pos.lat.to_string( unit=u.deg, decimal=True, precision=4)
        Lon_str = pos.lon.to_string( unit=u.deg, decimal=True, precision=4) 
        pass
    print('Your determined position is:\nLatitude:  ', Lat_str,'\nLongitude: ', Lon_str)
    input('Press Enter to continue')
    continue
    
while testMode == True: # test mode
    
    # test options
    numCases        = 100
    testTime        = Time.now()
    h_max           = 100 # meters
    h_min           = -10 # meters
    round2minutes   = True
    testAtm         = atm_std
    minAlt          = 5.0*u.deg
    simulateRefraction = True
    correctForRefraction = True
    numObs = 4
    
    # create numCases spaced test coordinates evenly spaced on the globe
    grid     = golden_spiral_grid(numCases)
    
    #turns off refraction in simulated observations if desired
    if simulateRefraction == False:
        testAtm = atm_No
        
              # [[est lat, true lat], [est lon, true lon], error ]
    results  = [[[ None, None ], [None, None], None] for i in range(numCases)] # initializing the results object
    testObs  = [[None, None, None] for i in range(numObs)] #initializing the test obs object
    for i in range(numCases):
        # pick a random height in our range
        h_test  = (h_max-h_min)*rng.random()+h_min
        #create an EarthLocation for the test case
        truePos = EarthLocation.from_geodetic(lat=grid[i].lat, 
                                              lon=grid[i].lon, 
                                              height=h_test, 
                                              ellipsoid = 'WGS84')
        
        # simulate observations from our test case location
        testObs = getTestObservations(truePos, 
                                      testTime, 
                                      numObs, 
                                      round2minutes,
                                      testAtm,
                                      minAlt
                                      )
        
        # Correct for refraction using air almanac tables
        if correctForRefraction:
            testObs = refractionCorrection(testObs)
        
        # estimate position using multilateration
        estPos  = multilateratePosition(testObs)
        
        # package results
        results[i][0] = [estPos.lat, truePos.lat]
        results[i][1] = [estPos.lon, truePos.lon]
        results[i][2] = haversine_np(truePos.lon, truePos.lat, estPos.lon, estPos.lat)
    
    # calculate mean and std. dev of error
    meanError = sum([results[i][2] for i in range(numCases)])/numCases
    stdDevError = np.std([results[i][2].value for i in range(numCases)])*u.km
    
    # setup results table
    print('\n Position Determination Simulated Test Results:\n') # setup results table print
    headers = ['Test Case','Estimated Latitude','True Latitude', 'Estimated Longitude', 'True Longitude', 'Total Position Error [km]']
    table = PrettyTable(headers)
    for i in range(numCases):
        if outputDMS:
            estLatStr  = results[i][0][0].to_string( unit=u.deg, decimal=False, sep='dms', precision=0)
            trueLatStr = results[i][0][1].to_string( unit=u.deg, decimal=False, sep='dms', precision=0)
            estLonStr  = results[i][1][0].to_string( unit=u.deg, decimal=False, sep='dms', precision=0)
            trueLonStr = results[i][1][1].to_string( unit=u.deg, decimal=False, sep='dms', precision=0)
            errorStr   = results[i][2].to_string( unit=u.km, precision=3) # round error to nearest meter
        else:
            estLatStr  = results[i][0][0].to_string( unit=u.deg, decimal=True, precision=4)
            trueLatStr = results[i][0][1].to_string( unit=u.deg, decimal=True, precision=4)
            estLonStr  = results[i][1][0].to_string( unit=u.deg, decimal=True, precision=4)
            trueLonStr = results[i][1][1].to_string( unit=u.deg, decimal=True, precision=4)
            errorStr   = results[i][2].to_string( unit=u.km, precision=3) # round error to nearest meter
        table.add_row( [ i, estLatStr, trueLatStr, estLonStr, trueLonStr, errorStr] )
    
    # output
    print(table)
    print("Mean Error = ", meanError.to_string( unit=u.km, precision = 3))
    print("Std. Deviation of Error = ", stdDevError.to_string( unit=u.km, precision = 3))
    break

    
    
    
    
        
        

