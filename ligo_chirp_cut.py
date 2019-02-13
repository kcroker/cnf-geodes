#!/usr/bin/python3

import scipy.integrate as integrate
import numpy as np
import scipy.interpolate as interpolate
import math
import sys

import units

# Adjust units to be in Mpc
# Planck 2018 median values
OmegaM = 0.315
OmegaL = 1 - OmegaM
H0 = 1
UnitAdjust = units.H0recip_MKS * units.c_MKS / units.m_per_Mpc

# Redshift range of the input LIGO map
a = 0
b = 1

zs = np.linspace(a, b, 1e3)
integrands = [math.sqrt( (1+z)**3 * OmegaM + OmegaL) for z in zs]

# Use Hogg (1999) notation
# Comoving distance to redshift z
discrete_DC = integrate.cumtrapz(integrands, zs)
DC = interpolate.interp1d(zs[1:], discrete_DC)

# Luminosity distance
discrete_DL = [(1+z)*DC(z)*UnitAdjust for z in zs[1:]]

# Now interpolate the inverse function to get redshift in terms of
# luminosity distance
z_of_DL = interpolate.interp1d(discrete_DL, zs[:-1])

detchirps_discrete = []

#
## 2D interpolator: include correction for mass ratio
#
qs = np.linspace(1e-3, 1, 25)
lookup = {}

# Read in all the chirps we know about
for line in open("data/luminosity_horiz_chirp.dat"):
    try:
        detchirp, luminosityd = map(float, line.split())

        # Make a dictionary so this becomes easy to get later
        lookup[detchirp] = luminosityd
        
        detchirps_discrete.append(detchirp)
    except Exception as e:
        print(e, file=sys.stderr)


def DLtoz_wrapper(dl):
    if dl < discrete_DL[0]:
        return 0.0
    else:
        return z_of_DL(dl)

flatz = []
chirpchirp = []
qq = []
for q in qs:
    for chirp in detchirps_discrete:
        try:
            # z_of_DL is determined entirely by Friedmann
            flatz.append(DLtoz_wrapper(lookup[chirp] * q * math.sqrt((1+q)/2.0)))
            chirpchirp.append(chirp)
            qq.append(q)
            #print("%.15e %.15e %.15e %.15e %.15e" % (q, chirp, flatz[-1], lookup[chirp], lookup[chirp] * q * math.sqrt((1+q)/2.0)), file=sys.stderr)
        except KeyError as e:
            print(e, file=sys.syderr)
            pass
        except ValueError as e:
            print(e, file=sys.stderr)
            pass
        
LIGOhorizon = interpolate.interp2d(chirpchirp, qq, flatz)

#
## 1D interpolator
#

# for line in open("data/luminosity_horiz_chirp.dat"):
#     try:
#         data = [float(x) for x in line.split()]
#         detectorchirps_discrete.append(data[0])
#         LIGOhorizon_discrete.append(z_of_DL(data[1]))
#     except Exception as e:
#         pass

# LIGOhorizon = interpolate.interp1d(detectorchirps_discrete, LIGOhorizon_discrete)

    
def horizon(m1, m2, z):
    try:
        zh = LIGOhorizon(source_chirp(m1, m2) * (1+z), float(m2)/m1)
        return zh
    except ValueError as e:
        print(e)
        
def source_chirp(m1, m2):
    return (m1*m2)**(3.0/5) / (m1 + m2)**(1.0/5)

def detector_chirp(z, m1, m2):
    detchirp = (1+z)*source_chirp(m1, m2)
    return detchirp

#
# DEBUG: take masses and output horizons
# 

# for line in sys.stdin:
#     try:
#         m1, m2 = map(float, line.split())
#         if m2 > m1:
#             m1, m2 = m2, m1

#         print(source_chirp(m1, m2), m2/m1)
#         print(LIGOhorizon(source_chirp(m1, m2), m2/m1))
#     except Exception as e:
#         print(e)
