import time
import math
import numpy as np

def calculate_position(star_observations):
    # Constants
    UNIX_epoch_JD = 2440587.5
    deg2rad = math.pi / 180.0
    rad2deg = 180.0 / math.pi
    NavStars = (  
        # Sideral Hour Angle (SHA) [deg], declination (DEC) [deg]
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

    numObs = len(star_observations)  # Number of observations
    obs = [[None for i in range(4)] for j in range(numObs)]  # [star number, elevation (Ho), time, GHA]

    for i in range(numObs):  # for each observation
        starNum, Ho, observation_time = star_observations[i]

        # Fill obs with data
        obs[i][0] = starNum
        obs[i][1] = Ho
        obs[i][2] = observation_time

        # Calculate GHA
        SHA = NavStars[starNum][0]  # lookup the star numbers & return SHA
        UT1 = obs[i][2]  # TODO implement UTC to UT1 calculation
        JD_UT1 = (UT1 / 86400.0) + UNIX_epoch_JD  # UT1 to Julian date
        T_u = (JD_UT1 - 2451545.0)  # from Sideral Time wiki page
        ERA = 360 * (0.7790572732640 + 1.00273781191135448 * T_u)  # [deg]
        GHA = (ERA + SHA) % 360.0  # TODO wrap from +180 to -180
        obs[i][3] = GHA

    A = np.zeros([numObs, numObs])  # initalize A as square matrix
    B = np.zeros([numObs, 1])  # initalize B as column matrix

    for i in range(numObs):  # populate A and B
        starNum = obs[i][0]  # which star is used for this observation
        DEC = NavStars[starNum][1] * deg2rad  # lookup DEC, convert to rads
        Ho = obs[i][1] * deg2rad  # lookup elevation, convert to rads
        GHA = obs[i][3] * deg2rad  # lookup, and convert to rads
        A[i] = [np.cos(DEC) * np.cos(GHA), np.cos(DEC) * np.sin(GHA), np.sin(DEC)]
        B[i] = [np.sin(Ho)]

    if np.linalg.det(A) < 0.1:  # check for matrix invertability
        raise ValueError("Solution is ill-defined")

    A_inv = np.linalg.inv(A)  
    X = np.dot(A_inv, B)  # solve matrix equation for cartesian position

    x, y, z = X[0], X[1], X[2]
    lat = np.arctan2(z, np.sqrt(x * x + y * y)) * rad2deg
    lon = np.arctan2(y, x) * rad2deg

    return lat[0], lon[0]  # return the estimated position

# Test the function
star_observations = [
    (3, 30, 1678378000),  # (star number, elevation, time)
    (5, 45, 1678378200),
    (7, 60, 1678378400)
]

latitude, longitude = calculate_position(star_observations)
print(f"Latitude: {latitude} degrees, Longitude: {longitude} degrees")
