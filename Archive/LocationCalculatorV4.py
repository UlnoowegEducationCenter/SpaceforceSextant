# Celestial Navigation Location Deterimination Version 4 (MLAT)
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


#### Load nav stars
   
# Read CSV of navigational stars
NavStars = genfromtxt( 'NavigationalStars.csv', delimiter=',' ) # read the star list csv
NavStars = np.delete( NavStars, 0, axis=0) # delete the label row
NavStars = np.delete( NavStars, (0,1,2), axis=1) # delete the name collumns
#row now corresponds to almanac number with the addition of 0 = polaris


#### Get observations from user

# get observation function
def getObservation(UsedStars):
    while True:
        print('-'*16)
        try:
            StarNum = input("Star Number: ")
            if StarNum == 'exit': # option to abort program execution
                exit("Program aborted")
            else:
                StarNum = int(StarNum)
                assert StarNum not in UsedStars
                assert StarNum in range(len(NavStars))
        except (ValueError, AssertionError): # excepting errors from invalid star numbers
            print('Invalid star number, try again')
            continue
        UsedStars.append(StarNum)
        break
    
    while True:
        try:
            obsEl = input("Elevation:   ")
            if obsEl == 'exit': # option to abort program execution
                exit("Program aborted")
            else:
                obsEl   = Latitude( obsEl, unit = u.deg )
        except (ValueError, TypeError): # excepting errors from invalid elevations
            print('Invalid elevation, try again')
            continue
        break
    
    # Calculate Greenwich Hour Angle (GHA) and declination
    ERA = Time.earth_rotation_angle( Time.now(), longitude = 0) # earth rotation angle at timenow
    SHA = Longitude( 360 - NavStars[StarNum,0], unit = u.deg ) # sideral hour angle from right ascention (from csv)
    GHA = Longitude( ERA + SHA, unit = u.deg ) # greenwich hour angle 
    dec = Latitude( NavStars[StarNum,1], unit = u.deg ) # declination from table
    
    observation = [StarNum, Time.now(), obsEl, GHA, dec ]
    
    return observation

# get 3 observations
StarsUsed = [] # initalize the list of used stars to prevent reuse
obs = list( ([0]*5, [0]*5, [0]*5) ) # initalize a list to store 3 observations
for i in range(3):
    obs[i] = getObservation(StarsUsed)


#### Calculate Position  

# Function to turn each observation into a row of the A matrix:
def obs2RowOfA(observation):
    Dec  = observation[4].rad # convert observation objects to floats in units of radians
    GHA  = observation[3].rad
    Ho   = observation[2].rad
    
    # create A matrix row 
    A_row = [np.cos(Dec)*np.cos(GHA), np.cos(Dec)*np.sin(GHA), np.sin(Dec)]
    A_row = A_row / np.sin(Ho)
    return A_row

# set up matrix equation
A = np.zeros([3,3]) # prealocate A
for i in range(3): # populate A
    A[i] = obs2RowOfA(obs[i])
B = np.ones( (3,1)) # create B

# solve matrix eqation
X = np.linalg.solve(A,B)

#convert cartesian to lat/lon
x, y, z = X[0], X[1], X[2]
Lat = Latitude( np.arctan2( z, np.hypot(x,y) )[0]*u.rad, unit = u.deg )
Lon = -Longitude( np.arctan2( y, x )[0]*u.rad, unit = u.deg, wrap_angle = 180*u.deg )

#### Output

# format lat/lon output strings
Lat_str = Lat.to_string( unit=u.deg, decimal=False, sep='dms', precision=0)
Lon_str = Lon.to_string( unit=u.deg, decimal=False, sep='dms', precision=0)

# main output
print('-'*16)
print('Your determined position is:\nLatitude:  ', Lat_str,'\nLongitude: ', Lon_str)

# Todo: option to add more values
