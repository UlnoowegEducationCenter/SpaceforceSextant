from django.shortcuts import render
from .forms import ObservationForm

import time
import math
import numpy as np


def index_view(request):
    if request.method == 'POST':
        forms = [ObservationForm(request.POST, prefix=str(x)) for x in range(3)]  # Assuming 3 observations
        if all([form.is_valid() for form in forms]):
            star_observations = [(form.cleaned_data['star_number'], form.cleaned_data['elevation'], time.time()) for form in forms]
            try:
                latitude, longitude = calculate_position(star_observations, force=request.POST.get('force', False))
            except ValueError as e:
                print(e) 
                return render(request, 'error.html', {'error_message': str(e)})
            return render(request, 'result.html', {'latitude': latitude, 'longitude': longitude})
    else:
        forms = [ObservationForm(prefix=str(x)) for x in range(3)]  # Assuming 3 observations
    return render(request, 'index.html', {'forms': forms})

def calculate_position(star_observations, force=False):
    # Constants and data
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

    for i in range(len(star_observations)):
        SHA = NavStars[star_observations[i][0]][0]
        UT1 = star_observations[i][2]
        JD_UT1 = (UT1 / 86400.0) + 2440587.5
        T_u = (JD_UT1 - 2451545.0)
        ERA = 360 * (0.7790572732640 + 1.00273781191135448 * T_u)
        GHA = (ERA + SHA) % 360.0
        star_observations[i] = list(star_observations[i]) + [GHA]

    # Setup and solve matrix equation
    A = np.zeros([3, 3])
    B = np.zeros([3, 1])
    for i in range(len(star_observations)):
        DEC = NavStars[star_observations[i][0]][1] * deg2rad
        Ho = star_observations[i][1] * deg2rad
        GHA = star_observations[i][3] * deg2rad
        A[i] = [np.cos(DEC) * np.cos(GHA), np.cos(DEC) * np.sin(GHA), np.sin(DEC)]
        B[i] = [np.sin(Ho)]

    # if np.linalg.det(A) < 0.1:
    #    if not force:
    #        raise ValueError("Warning: solution is ill-defined.")

    # Process observations
    # if np.linalg.det(A) < 0.1:
    #    raise ValueError("Warning: solution is ill-defined.")
    A_inv = np.linalg.inv(A)
    X = np.dot(A_inv, B)

    x, y, z = X[0], X[1], X[2]
    lat = np.arctan2(z, np.sqrt(x * x + y * y)) * rad2deg
    lon = np.arctan2(y, x) * rad2deg

    return lat[0], lon[0]
