#!/usr/bin/python3
#
# fast2dbin.py
# Copyright(C) 2019 Kevin Croker
#
# Do awful and nasty things to coax gnuplot into making attractive heatmaps.
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


import sys
import numpy as np

chirps = []
merger_zs = []
for datum in sys.stdin:
    try:
        # I'm sure theres a pythonic way of slicing this...
        zf, a, b, m1, m2 = map(float, datum.split()[:5])

        chirps.append( (m1 * m2)**(3.0/5)/(m1 + m2)**(1.0/5) )
        merger_zs.append(zf)
    except Exception as e:
        print(e, file=sys.stderr)
    
# x shows up on y, y shows up on x... thx python dox
histo, xedges, yedges = np.histogram2d(merger_zs, chirps, normed=False,  bins=int(sys.argv[1]))

# Divide out histo by the total count
histo = histo / len(chirps)

# Since we are displaying features in redshift,
# even time-uniform behaviour will exhibit monotonic increase...

# Make some axes strings for gnuplot
# Eh.  For heatmaps we don't want bin centers anymore...
# For gnuplot matrix image plots, apparently we do...
binvars = [ (yedges[1:], yedges[0]), (xedges[1:], xedges[0]) ]
headers = []
minors = []
for n,herp in enumerate(binvars):
    binrs,binl = herp
    headers.append([])
    minors.append([])
    for r in binrs:
        headers[n].append(binl + (r - binl)/2.0)
        minors[n].append(binl)
        binl = r

# gnuplot is retarded and forgets to latex rowheaders and columnheaders when
# splotting a matrix type
# I'm so done with this program after this paper series omfg
for m,header in enumerate(headers):
    print("# (", end='')
    for n,value in enumerate(header):
        if m == 0:
            if (header[-1] - header[0]) < 10:
                print('\'\scalebox{0.5}{%.1f}\' %d, ' % (value, n), end='')
            else:
                print('\'\scalebox{0.5}{%.0f}\' %d, ' % (value, n), end='')
        else:
            print('\'\scalebox{0.5}{%.3f}\' %d, ' % (value, n), end='')
                        
    print(")")

# Now output minor tics so that we can draw a grid 
# on the box boundaries.
# Notice that this style uses row numbers...
for m,minor in enumerate(minors):
    print("# (", end='')
    for n,value in enumerate(minor):
        if n == len(minor)-1:
            print('"" %f 1' % (n + 0.5), end='')
        else:
            print('"" %f 1, ' % (n + 0.5), end='')
                        
    print(")")

# Print out max and min for y and x, so that we can 
# manually reconstruct the scaling.
print("# %.15e" % headers[1][0])
print("# %.15e" % headers[1][-1])
print("# %.15e" % headers[0][0])
print("# %.15e" % headers[0][-1])

#print("burn ", end='')
#for col in headers[0]:
#    print("\scalebox{0.5}{%.1f} " % col, end='')
#print('')

# Output the dirty
for n,row in enumerate(histo):

    # Compute the 
    #    rowstr = "\scalebox{0.5}{%.2f} " % headers[1][n]
    rowstr = ""
    for dat in row:
        rowstr = rowstr + " %.15e" % dat

    print(rowstr)


