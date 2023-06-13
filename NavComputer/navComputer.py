# Celestial Navigation - Position Determination for Adafruit Feather RP2040
# By Chris Seymour, EIT
# Copyright Ulnooweg Education Centre, 2023, All rights reserved
###############################################################################


#### Import packages

#import board #TODO uncomment when running on device
import sys 
import time # TODO replace with RTC module
import math
import numpy as np # TODO replace with ulab.numpy


#### Constants
UNIX_epoch_JD = 2440587.5
deg2rad = math.pi / 180.
rad2deg = 180. / math.pi
NavStars = ( # tuple of stars Sideral Hour Angle (SHA) [deg], and declination (DEC) [deg], index is star number with polaris = 0 Sourced from wikipedia
    (319.,89.),
    (358.,29.),
    (354.,-42.),
    (350.,56.),
    (349.,-18.),
    (336.,-57.),
    (328.,23.),
    (316.,-40.),
    (315.,4.),
    (309.,50.),
    (291.,16.),
    (282.,-8.),
    (281.,46.),
    (279.,6.),
    (279.,29.),
    (276.,-1.),
    (271.,7.),
    (264.,-53.),
    (259.,-17.),
    (256.,-29.),
    (245.,5.),
    (244.,28.),
    (234.,-59.),
    (223.,-43.),
    (222.,-70.),
    (218.,-9.),
    (208.,12.),
    (194.,62.),
    (183.,15.),
    (176.,-17.),
    (174.,-63.),
    (172.,-57.),
    (167.,56.),
    (159.,-11.),
    (153.,49.),
    (149.,-60.),
    (149.,-36.),
    (146.,19.),
    (140.,-61.),
    (138.,-16.),
    (137.,74.),
    (127.,27.),
    (113.,-26.),
    (108.,-69.),
    (103.,-16.),
    (97.,-37.),
    (96.,13.),
    (91.,51.),
    (84.,-34.),
    (81.,39.),
    (76.,-26.),
    (63.,9.),
    (54.,-57.),
    (50.,45.),
    (34.,10.),
    (28.,-47.),
    (16.,-30.),
    (14.,15.)
    )


#### Get observations from user
numObs = 3
obs = [ [None for i in range(4)] for j in range(numObs) ] # star number, elevation (Ho), time, GHA
for i in range(numObs):
    
    # get star number from user & validate input
    while True:
        try:
            starNum = int(input("Star Number: "))
            #assert starNum not in obs[:][0] # TODO figure out why this isnt working - supposed to prevent reuse of star numbers
            assert starNum in range(len(NavStars))
        except (ValueError, AssertionError): # excepting errors from invalid star numbers
            print('Invalid star number, try again')
            continue
        obs[i][0] = starNum
        break
    
    # get elevation from user and validate input
    while True:
        try:
            Ho = float(input("Elevation: "))
            assert (-90.<= Ho) & (Ho <= 90.)
        except (ValueError, TypeError, AssertionError):
            print("invalid elevation, try again")
            continue
        obs[i][1] = Ho
        break
    
    # Get time
    t_unix = time.time() # TODO replace with RTC function
    obs[i][2] = t_unix
    
    # calculate GHA
    SHA = NavStars[starNum][0] # lookup the star number in the observation and return the SHA
    UT1 = obs[i][2] # universal time calculated from unix time (from observation via. RTC) #TODO correction factor from UTC to UT1
    JD_UT1 = (UT1/86400.) + 2440587.5 # UT1 to julian date
    T_u = JD_UT1 - 2451545.0 # this line and the next one are from wikipedia's page on Sideral Time
    ERA = 360*(0.7790572732640 + 1.00273781191135448*T_u) # earth rotation angle in degrees
    GHA = (ERA + SHA) % 360. # TODO make this work like longitude (wraps from +180 to -180)
    obs[i][3] = GHA


#### Setup and solve matrix equation

A = np.zeros([3,3]) # initalize A as square matrix
B = np.zeros([3,1]) # initalize B as column matrix
for i in range(numObs): #populate A and B
    starNum = obs[i][0] # which star is used for this observation
    DEC = NavStars[starNum][1]*deg2rad # lookup declination of star and convert to rads
    Ho = obs[i][1]*deg2rad # lookup elevation of observation and convert to rads
    GHA = obs[i][3]*deg2rad
    A[i] = [np.cos(DEC)*np.cos(GHA), np.cos(DEC)*np.sin(GHA), np.sin(DEC)]
    B[i] = [np.sin(Ho)]

try:
    #TODO: insert a check for matrix invertability, e.g. if det(A)<0.1 raise a warning
    A_inv = np.linalg.inv(A) # TODO investigate whether this needs to be inverse or will transpose work?
    print("A_inv is a ",type(A_inv),"\nA_inv=",A_inv) # TODO comment out in production code
    X = np.dot(A_inv,B) #solve matrix equation for cartesian position
except(ValueError):
    print('No solution exists, please try again.')
    sys.exit(1)

#### Convert to Lat / Lon and print
x, y, z = X[0], X[1], X[2]
lat = np.arctan2(z,np.sqrt(x*x+y*y))*rad2deg
lon = np.arctan2(y,x)*rad2deg

print("Latitude = ",lat)
print("Longitude = ",lon)