#!/usr/bin/python3
#
# odb_analysis_revisited.py
# Copyright(C) 2019 Kevin Croker
#
# Implements the GEODE mass function model given in Eqn. (18) of Paper III.
# An integral over the blueshifted Salpeter initial remnant distributions,
# weighted by Madau and Dickinson's SFR.  The normalization is not physically meaningful,
# but that doesn't matter.
#
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

import tabcalc as tc
import numpy as np
import math
import sys

# These are in \rho_cr and reciprocal hubble units
Omegam = 0.315
OmegaL = 1 - Omegam
H0 = 1
R = 0.27

# Converts from astronomer (Msol/Mpc^3; seconds) units to (\rho_cr; reciprocal hubble) units
# UnitCorrection = 0.277

# Lets use Msol/Mpc^3
UnitCorrection = 1

# From Madau & Dickinson 2014.
# Generalized to include the Planck linear model for DE evolution: w_de = c0 \pm (1-a)*c1  (I forgot which + or -, see the paper)
# (though this may introduce very small inconsistencies with the rest of the Dickinson model)
def psi(z):
    return UnitCorrection*0.015*(1+z)**2.7/(1 + ( (1+z)/2.9 )**5.6)
def H(z,c0,c1):
    return H0*math.sqrt(Omegam*(1+z)**3 + OmegaL*(1+z)**(3*(c0 + c1 + 1))*math.exp(-3*c1*z/(z + 1)))
def dStardz(z,c0,c1):
    if z > 20:
        return 0
    else:
        return (R - 1)*psi(z)/((1+z)*H(z,c0,c1))

# Return values for perfect dark energy, which should be a very good approximation
def dStarda(a):
    return -1/a**2 * dStardz(1.0/a - 1, -1, 0)
    
# And get how things look at the present day
# Get a domain of masses
def makeDistribution(params, max_mass, resolution):
    
    # Unpack the parameters
    alpha, ml, mh, k, a = params

    print("Powerlaw: alpha = %f, ml = %f, mh = %f, k = %f, a = %f" % params, file=sys.stderr)

    # Set mc = 0.1 for the Salpeter STELLAR distribution
    mc = 0.1
    
    # Get the prefactor
    prefactor = a**(k*(alpha - 1)) * (alpha - 2) / (mc*(ml**(1-alpha) - mh**(1-alpha)))
   
    masses = [m for m in np.linspace(ml, max_mass, resolution)]

    # Lists for the final values
    distribution = []
    final_masses = []

    # Regard z = 20 as onset of Pop III formation
    for mass in masses:
        # Make a linspace to integrate on
        amax = a*min(1, (mh/mass)**(1.0/k))
        amin = a*(ml/mass)**(1.0/k)
        a0_domain = [a0 for a0 in np.linspace(amin, amax, resolution)]

        # KC 9/30/18
        # Notice how the lower-limit works.  It goes far enough back in time
        # so that an object with the CUTOFF mass made at that time would grow
        # to the mass being interrogated by "mass"
        #
        # So for things at the high end of the spectrum, this just isn't possible.
        # That doesn't mean that such things don't exist, because objects far from the
        # cutoff can easily grow to that size.

        final_masses.append(mass)
        distribution.append( mass**(-alpha) * prefactor * (tc.TblBackedFunction(a0_domain, [dStarda(a0)*a0**(k*(1-alpha)) for a0 in a0_domain])).integrateTop() )

    finaltable = tc.TblBackedFunction(final_masses, distribution)
    finaltable.addComment("BH mass function, m^-|alpha| remnant distribution (cutoff at mc), a^k growth, at time a") 
    finaltable.addComment("|alpha| = %f, ml = %f, mh = %f, k = %f, a = %f" % params)

    return finaltable

# Test?
#   adomain = [a for a in np.linspace(0.02, 1, 1000)]
#   (tc.TblBackedFunction(adomain, [dStarda(a) for a in adomain])).writeout("dStarda_test.dat")
# Confirmed.

# Crude check for arguments
if len(sys.argv) < 6:
    print("Usage: ./odb_analysis_revisited.py <stellar powerlaw exponent> <remnant low mass cutoff> <remnanat high mass cutoff> <GEODE powerlaw growth> <obs. redshift>")
    exit(1);

# Usually Salpeter |alpha| 
alpha = float(eval(sys.argv[1]))

# Remnant low cutoff in solar masses
ml = float(eval(sys.argv[2]))

# Remnant high cutoff in solar masses
mh = float(eval(sys.argv[3]))

# Growth k
k = float(eval(sys.argv[4]))

# Get a scale factor from the redshift
a = 1.0/(1 + float(eval(sys.argv[5])))

# Get a distribution up to high masses
# 1e4 MSol with a resolution of 1e3
finaltable_highmass = makeDistribution((alpha, ml, mh, k, a), mh, 1e2)
finaltable_lowmass = makeDistribution((alpha, ml, mh, k, a), 1e2, 1e3)

# Assemble this into the full distribution
newdist = []

# Block 0
# Cut it at the low mass of the remnants, because they always grow heavier
newmasses = [m for m in np.linspace(ml, mh, 1e5)]
newdist = [ finaltable_lowmass(m) if m < 1e2 else finaltable_highmass(m) for m in newmasses]
combineddistro = tc.TblBackedFunction(newmasses, newdist)
combineddistro.addComment("GEODE distribution, for gap and cap Salpeter, with uniform resolution")
combineddistro.addComment("alpha = %f, ml = %f, mh = %f, k = %f, z = %f" % (alpha, ml, mh, k, float(eval(sys.argv[5]))))
combineddistro.writeout()

# Seems like it worked?
print("")

# Block 1
finalcdf = tc.indefinite_integral_top(combineddistro)
norm = finalcdf.y[-1]
finalcdf = finalcdf/finalcdf.y[-1]
finalcdf.addComment("CDF of the above GEODE mass function")
finalcdf.addComment("Normalized by (final) CDF value = %f, at mass %f" % (norm, finalcdf.x[-1]))
finalcdf.writeout()

# Print a newline
print("")

# Block 2
# Now we want to get the slope of these graphs, in logscale.
loglog_finaltable = tc.TblBackedFunction([np.log10(x) for x in finaltable_lowmass.x], [np.log10(y) for y in finaltable_lowmass.y])
loglog_finaltable.writeout()

# Print a newline
print("")

# Block 3
# Take a derivative
loglog_finaltable_diff = tc.derivative(loglog_finaltable)
loglog_finaltable_diff.addComment("Slope of GEODE mass function in log-log, with above parameters")
loglog_finaltable_diff.writeout()
