#!/usr/bin/python3
#
# geode-belczynski.py
# Copyright(C) 2019 Kevin Croker
#
# Completes the mergers of the DCOs formed in Belczynski and Mink (2015)'s data sets
# Only difference is that we:
#   1) track the blueshift effect in the GEODEs mass
#   2) sample redshifts according to some distribution and run their data sets anchored at these redshifts
#
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

import scipy.stats as st
import numpy as np
import math
import sys
import subprocess

# And the integration stub, so we don't need to shell out.
import geode_binary_module
import fastbin
import units
import cosmo_distributions

if len(sys.argv) < 4:
    print("Usage: ./%s <output file> <geode-geode | geode-ns | bh-bh> <madau | uniform>", file=sys.stderr)
    exit(1)
    
#################
#
# LOAD A MINK & BELCZYNSKI 2015 DATA SET
#
###################

# Read Mink & Belczynski format from stdin
# Ma (msol) | Mb (msol) | a (AU) | eccentricity
data = []
for datum in sys.stdin:

    Ma, Mb, Ri, ei = map(float, datum.split()[:4])

    # Place the larger one into Ma.  This will always be the GEODE
    if Mb > Ma:
        Ma,Mb = Mb,Ma

    # Switch from solar radii to AU.  This was a gnarly bug.
    Ri = float(Ri) * units.AU_per_Rsol
    
    # Our order is different, sorry
    data.append((Ri, ei, Ma, Mb))

################
#
# SPECIFICS FOR OUR SIMULATIONS
#
##################

# Enable buffering on status so we can watch?
status = open("current_status", "w", 1)

# We will check response to weighting in redshift and tim
cnt=100

zis = None
m1Function = None
m2Function = None
model = None
z_max = None

# Any larger and we are getting into low metallicity territory
# For consistency on the production runs, set everything to the same
# maximum redshift.
z_max = 4

if sys.argv[2] == "geode-ns":
    model = geode_binary_module.geode_ns_binary
    m1Function = lambda mi,scale: mi*scale
    m2Function = lambda mi,scale: mi
    
elif sys.argv[2] == "geode-geode":
    print("Interpreting data as geode-geode", file=status)
    model = geode_binary_module.geode_binary
    m1Function = lambda mi,scale: mi*scale
    m2Function = m1Function
    
elif sys.argv[2] == "bh-bh": 
    model = geode_binary_module.bh_binary
    m1Function = lambda mi,scale: mi
    m2Function = m1Function

else:
    print("Error: unknown DCO %s" % sys.argv[2], file=sys.stderr)
    exit(1)
    
print("Set z_max to %f" % z_max, file=status)

# Get a distribution of redshifts that is sane, given our model to investigate and aLIGO canonical
# mass cuts
#
# Check for systematics against LIGO canonical assumptions.
# In particular, we want to make sure the redshift maps have hotspots
# that are robust to the time distribution assumptions, so that it
# reflects the actual GEODE blueshift munging
redshiftDistribution = None 
if sys.argv[3] == 'madau':
    redshiftDistribution = cosmo_distributions.madaudick(z_max, 0, 1e-5)
elif sys.argv[3] == 'uniform':
    redshiftDistribution = cosmo_distributions.timeuniform(z_max, 0, 1e-5)
else:
    print("Error: unknown time distribution %s" % sys.argv[3], file=sys.stderr)
    exit(2)


###################
#
# RUN THE MERGERS
#
###################
#

# Because the lsoda integrator puts turds on my stdout
f = open(sys.argv[1], "w")

print("# Redshift (final) | Separation (final) | Eccentricity (final) | Primary (final) | Secondary (final) | Angular Momentum (final) | Redshift (formation) | Separation (formation) | Eccentricity (formation) | Primary (formation) | Secondary (formation)", file=f)

# Check robustness to sampling strategy
# Here we check each binary against a variety of redshifts
n = 0
hits = []
resolution = 1e-3

# Scale factor cutoff for no merger
geode_binary_module.no_merger_cutoff = 0.999

for dcnt, datum in enumerate(data):

    # Get a distinct redshift distribution
    try:
        zis = redshiftDistribution.rvs(size=cnt)
    except ValueError as e:
        # rvs asked for a redshift outside of the range.  wtf
        print("Pulled bad redshift range?  (setup looks fine.  wtf)", file=status)
        zis = redshiftDistribution.rvs(size=cnt)
    
    # Unpack
    Ri,ei,m1,m2 = datum

    print("Processing Belczynski & Mink binary #%d" % (dcnt+1), file=status)

    for zi in zis:

        # Get the scale factor to anchor at
        a_i = 1.0/(1 + zi)
        recipai3 = 1.0/a_i**3
        geode_binary_module.a_start = a_i

        # Make the grid here (this is slow again)
        alist = np.linspace(a_i, 1.0, 1.0/resolution)

        # Track number processed
        n += 1
            
        # Evolve
        a_f, L_f, R_f, e_f = geode_binary_module.evolve(alist, datum, model)
 
        # Suitably blueshift the masses
        # (Eliminate divisions within the loop)
        blueshift = a_f**3 * recipai3
        m1_f = m1Function(m1, blueshift)
        m2_f = m2Function(m2, blueshift)

        # Cut if it didn't merge.  (Don't attempt to do the LIGO cuts!)
        z_f = 1.0/a_f - 1
        if a_f > geode_binary_module.no_merger_cutoff:
            continue

        # The result, as outputted by the previous stub...
        result = (z_f, R_f, e_f, m1_f, m2_f, L_f, zi) + datum
    
        # Otherwise, we registered a hit
        print("Pair %.5d HIT! (rate = %.3f)" % (n, (len(hits)+1)/n), file=status)
        hits.append(result)

        # Output them as we find them
        print("%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e" % result, file=f)

    # Flush stdout
    f.flush()

# Append rate information at the end
print("# Rate: %f" % ((len(hits)+1)/n), file=f)

# ##################
# #
# # DEBUG ANALYZE THE HITS 
# #
# ##################

# # Advance to the next block
# print("", file=f)
# print("# Binned histogram of masses!", file=f)

# # Make a histogram of the final masses
# mfs = []
# for hit in hits:
#     mfs.append(hit[3])
#     mfs.append(hit[4])

# fastbin.dumpbinned(mfs, np.arange(0,115, 1), f)
