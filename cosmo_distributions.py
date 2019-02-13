#!/usr/bin/python3
#
# cosmo_distributions.py
# Copyright(C) 2019 Kevin Croker
#
# Table backed stats distributions that work quickly.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#


import scipy.integrate as integrate
import scipy.interpolate as interpolate
import scipy.stats as st
import numpy as np
import math


# Distribution in formation redshift
OmegaM = 0.315
OmegaL = 1 - OmegaM
H0 = 1
R = 0.27

# Converts from astronomer units to \rho_cr/reciprocal hubble
UnitCorrection = 0.277

# From Madau & Dickinson 2014.
def psi(z):
    return UnitCorrection*0.015*(1+z)**2.7/(1 + ( (1+z)/2.9 )**5.6)
def H(z):
    return H0*math.sqrt(OmegaM*(1+z)**3 + OmegaL)
def dStardz(z):
    return (R - 1)*psi(z)/((1+z)*H(z))

# LIGO analyses sometimes consider things that are uniform in *time*, convert this to a statement
# in redshift
def dtdz(z):
    return 1/((1+z)*H(z))

#
# Construct the comoving volume correction
#
zs_rev = np.linspace(-0.2, 4.2, 1e3)
igrands = [1/H(z) for z in zs_rev]
discrete_DC = integrate.cumtrapz(igrands, zs_rev)
DC = interpolate.interp1d(zs_rev[1:], discrete_DC)

#
# Proportionality factors here don't matter, since we are going to build a weighted
# distribution in redshift from this, and then normalize it.  We're not trying to 
# predict absolute rates, only relative rates.
#
def comoving_shell_volume(z):
    return DC(z)**2  / H(z)

# for z in zs_rev[1:]:
#     print("%.15e %.15e %.15e" % (z, comoving_shell_volume(z)))

# exit(1)

# The only way to do this in a sane way is to interpolate it
# and then do the dirty
#
# https://stackoverflow.com/questions/51946345/stats-rv-continuous-slow-when-when-using-custom-pdf
#
# This is (dN_star/dz) dz = dN_star
# the distribution for formation times will be weighted relative to the number of objects expected
# to be produced at that time.  This means the comoving stellar density rate times a comoving volume shell,
# assuming fixed average stellar mass.
#
class madaudick_gen(st.rv_continuous):

    def _argcheck(self, a, b, res):
        try:
            None if self.res else None
        except AttributeError:
            # Set up everything
            self.res = res

            # Compute the correct extension on the bottom side and the right number steps,
            # so that the lower margin falls in the interpolated range for the integrated functions
            # which lose a domain point.
            N = round(1/res)
            zs = np.linspace(a+(a - b)*res, b-(a - b)*res, N+3)

            # Now when the top or bottom is chopped for the cdf or ppf
            # we still interpolate over [a,b] inclusive
            ys = [dStardz(z)*comoving_shell_volume(z) for z in zs]

            # Now use the above strategy, but normalize
            discrete_cdf1 = integrate.cumtrapz(ys, zs)
            self.d = discrete_cdf1[-1]

            # Normalize the CDF
            discrete_cdf1 = discrete_cdf1/self.d

            # The interpolated CDF and PPF are then normalized correctly
            self.cdf1 = interpolate.interp1d(zs[1:], discrete_cdf1)
            self.ppf1 = interpolate.interp1d(discrete_cdf1, zs[:-1])
            
        return (a > b) and (res < 1) and (res > 1e-10) 
    
    def _cdf(self, x, a, b, res):
        return self.cdf1(x)
    
    def _ppf(self, x, a, b, res):
        return self.ppf1(x)

    # Explicitly normalize the PDF
    def _pdf(self, x, a, b, res):
        return dStardz(x)/self.d

madaudick = madaudick_gen(a=4, b=0)

#
# Make the timeUniform distribution
#
class timeuniform_gen(st.rv_continuous):

    def _argcheck(self, a, b, res):
        try:
            None if self.res else None
        except AttributeError:
            # Set up everything
            self.res = res

            # Compute the correct extension on the bottom side and the right number steps,
            # so that the lower margin falls in the interpolated range for the integrated functions
            # which lose a domain point.
            N = round(1/res)
            zs = np.linspace(a+(a - b)*res, b-(a - b)*res, N+3)
#            print(zs)
#            exit(1)
            
            # So this is that the rate of production per comoving volume is constant
            ys = [dtdz(z)*comoving_shell_volume(z) for z in zs]

            # Now use the above strategy, but normalize
            discrete_cdf1 = integrate.cumtrapz(ys, zs)
            self.d = discrete_cdf1[-1]
            discrete_cdf1 = discrete_cdf1/self.d
            
            self.cdf1 = interpolate.interp1d(zs[1:], discrete_cdf1)
            self.ppf1 = interpolate.interp1d(discrete_cdf1, zs[:-1])
            
        return (a > b) and (res < 1) and (res > 1e-10) 
    
    def _cdf(self, x, a, b, res):
        return self.cdf1(x)
    
    def _ppf(self, x, a, b, res):
        return self.ppf1(x)

    def _pdf(self, x, a, b, res):
        return dtdz(x)/self.d

timeuniform = timeuniform_gen(a=4, b=0)

#
# Amazingly, this is not standardized??
# I feel like I'm unnecessarily repeating a bunch of stuff.
# And there are interpolation value errors every now and then, so there is an off-by-one
# in the interpolation code somewhere.
#

# (Silently fail even if the uncropped distribution is not normalizable, for speed.)
# Maybe we need to use generators.
class negplaw_gen(st.rv_continuous):

    # Use _argcheck to initialize and set things up?
    def _argcheck(self, a, b, alpha):
        try:
            None if self.alpha else None
        except AttributeError:
            self.a = a
            self.b = b
            self.alpha = alpha

            if b is None:
                self.d = (alpha - 1)/a**(1-alpha)
            else:
                self.d = (1-alpha)/(b**(1-alpha) - a**(1-alpha))

        # Do I need the (1-alpha) < 0 check?
        return (alpha > 0) and (a < b)

    # Define the _cdf
    def _cdf(self, x, a, b, alpha):
        return self.d/(1-alpha) * (x**(1-alpha) - a**(1-alpha))
    
    # This is the inverse of the _cdf
    # This bastard MUST be defined in order to get anything reasonable out of your distribution.
    # This is a somewhat awful function in general. Ugh.
    def _ppf(self, y, a, b, alpha):
        return (y/self.d *(1-alpha) + a**(1-alpha))**(1.0/(1-alpha))
    
    def _pdf(self, x, a, b, alpha):
        return self.d * x**(-self.alpha)

negplaw = negplaw_gen(a = 1.0, b = None)


# # Verify?
#
# import fastbin
#
# # Are separations uniform in log spaced bins?
# # (This won't quite be correct, because we are binning logrithmically)
# bins = np.logspace(np.log10(minR/10.0), np.log10(maxR*10))
# fastbin.dumpbinned(Ris, bins)
# # Seems to work.

# # Are the mass ratios uniform in normal scale?
# fastbin.dumpbinned(qs, [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1])
# # Seems to work.

# # Are the masses themselves powerlaw distributed?
# fastbin.dumpbinned(m2s, np.logspace(np.log10(minMass/10.0), np.log10(maxMass*10)))
# # Seems to work.

# # Are the redshifts distributed with respect to star formation rate?
# fastbin.dumpbinned(zis, np.linspace(0, 8, 64))
# # Seems to work.

# # Is the double uniform distribution for BH sane?
# fastbin.dumpbinned(m1s, np.linspace(0,20,100), sys.stdout)
# # Looks solid.

# # Is the uniform in time distribution (displayed in redshift) sane?
# q = timeuniform(6, 0, 1e-5)
# zs = q.rvs(size=cnt)
# fastbin.dumpbinned(zs, np.linspace(0,6,60), sys.stdout)
# # Smexy.
# exit(1)
# Verified.
