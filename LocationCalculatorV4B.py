# Celestial Navigation Location Deterimination Version 4B - Automated Test
# By Chris Seymour, EIT
# Copyright Ulnooweg Education Centre, 2023, All rights reserved
###############################################################################


#### import packages
import numpy as np
from sys import exit
from numpy import genfromtxt
from astropy.time import *
from astropy.coordinates import *
from astropy import units as u


#### Load Navigational Stars

# Read CSV of navigational stars
NavStars = genfromtxt( 'NavigationalStars.csv', delimiter=',' ) # read the star list csv
NavStars = np.delete( NavStars, 0, axis=0) # delete the label row
NavStars = np.delete( NavStars, (0,1,2), axis=1) # delete the name collumns
#row now corresponds to almanac number with the addition of 0 = polaris


#### Define Functions

# fn that creates a list of test cases
def makeTestCases(numCases):
    TestCaseNum = numCases # number of test cases
    Grid = golden_spiral_grid(TestCaseNum)
    HeightMax = 100 # max observation random height
    Heights = np.random.rand(TestCaseNum,1)*HeightMax
    TestCases = [None]*TestCaseNum
    for i in range(TestCaseNum):
        TestCases[i] = EarthLocation.from_geodetic(lat=Grid[i].lat, lon=Grid[i].lon, height = Heights[i], ellipsoid = 'WGS84')
    return TestCases

# fn that creates test observations from the test cases
def getTestObservations(TestCases,time):
    
    numObservationsPerLocation = 3
    numTestCases = len(TestCases)
    observations = [ [None for i in range(numObservationsPerLocation)] for j in range(numTestCases) ]
    
    for i in range(len(TestCases)): # for each test case
        testFrame = AltAz(location=TestCases[i],obstime=time) #create the observation reference frame
        testStars = [None for j in range(numObservationsPerLocation)] # initalize the test stars list
        for j in range(numObservationsPerLocation): 
            while True:
                candidateStar = np.random.randint(low=0, high=len(NavStars))  
                if candidateStar in testStars:
                    continue
                # elif visible==False: #calculate visibility
                #     continue
                else:
                    testStars[j]=candidateStar
                    break  
        
        for k in range(3): # for each sighting from test location
            testStar = testStars[k]
            ra   = NavStars[testStar][0] # read right ascention of selected star
            dec  = NavStars[testStar][1]
            Star = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
            Star = Star.transform_to(testFrame)
            alt  = Star.alt
            
            ERA = Time.earth_rotation_angle( time, longitude = 0) # earth rotation angle at the observation time
            SHA = Longitude( 360 - ra, unit = u.deg ) # sideral hour angle from right ascention (from csv)
            GHA = Longitude( ERA + SHA, unit = u.deg ) # greenwich hour angle 
            dec = Latitude( dec, unit = u.deg ) # declination from table
        
            observations[i][k] = [testStar, time, alt, GHA, dec]

    return observations

# Function to turn each observation into a row of the A matrix:
def obs2RowOfA(observation):
    
    Dec  = observation[4].rad # convert observation objects to floats in units of radians
    GHA  = observation[3].rad
    Ho   = observation[2].rad
    
    # create A matrix row 
    A_row = (np.cos(Dec)*np.cos(GHA), np.cos(Dec)*np.sin(GHA), np.sin(Dec) )
    A_row = A_row/np.sin(Ho)

    return A_row


#### Run Tests

numCases = 60
TestCases = makeTestCases(numCases)
TestObservations = getTestObservations(TestCases,Time.now()) ##### this one is fucked for some reason

print(TestObservations)

#### Calculate Position
DetPos = [ [None, None] ]*numCases
for i in range(numCases): #for every test case
    A = np.zeros( [3,3] )
    B = np.ones( (3,1) )
    
    for j in range(3):
        obs = TestObservations[i][j]
        A[j] = np.transpose( obs2RowOfA(obs) )
    
    X = np.linalg.solve(A,B)
    x, y, z = X[0], X[1], X[2]
    Lat = Latitude( np.arctan2( z, np.hypot(x,y) )[0]*u.rad, unit = u.deg )
    Lon = Longitude( np.arctan2( y, x )[0]*u.rad, unit = u.deg, wrap_angle = 180*u.deg )
    DetPos[i] = [ [Lat, Lon] ]

#### print outputs
for i in range(numCases):
    print('Test Case ',i+1,' Determined Position: ',DetPos[i],'Actual Position: ',TestCases[i])