# Celestial Navigation - Position Determination for Adafruit Feather RP2040
# By Chris Seymour, EIT
# Copyright Ulnooweg Education Centre, 2023, All rights reserved
###############################################################################

###############################################################################
# Setup
###############################################################################

# Import packages
import time
import math
import numpy as np

# Constants
UNIX_epoch_JD = 2440587.5
deg2rad = math.pi / 180.0
rad2deg = 180.0 / math.pi
NavStars = (  # Sideral Hour Angle (SHA) [deg], declination (DEC) [deg]
    (319.0, 89.0),
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

###############################################################################
# Main loop
###############################################################################

while True:

    # Get observations from user
    numObs = 3
    obs = [[None for i in range(4)] for j in range(numObs)]
    # [star number, elevation (Ho), time, GHA]
    for i in range(numObs):  # for each observation
        # get star number from user & validate input
        while True:
            try:
                starNum = int(input("Star Number: "))
                assert starNum in range(len(NavStars))
            except (ValueError, AssertionError):
                print("Invalid star number, try again")
                continue
            obs[i][0] = starNum
            break

        # get elevation from user and validate input
        while True:
            try:
                Ho = float(input("Elevation: "))
                assert (-90.0 <= Ho) & (Ho <= 90.0)
            except (ValueError, TypeError, AssertionError):
                print("invalid elevation, try again")
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
        starNum = obs[i][0]  # which star is used for this observation
        DEC = NavStars[starNum][1] * deg2rad  # lookup DEC, convert to rads
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
        supervisor.reload()

    # Convert to Lat / Lon
    x, y, z = X[0], X[1], X[2]
    lat = np.arctan2(z, np.sqrt(x * x + y * y)) * rad2deg
    lon = np.arctan2(y, x) * rad2deg

    # Print results
    print("Determined position is:")
    print("Latitude:  ", str(lat[0]), " degrees")
    print("Longitude: ", str(lon[0]), " degrees")