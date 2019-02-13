#!/usr/bin/python3
#
# skript.py
# Copyright(C) 2019 Kevin Croker
#
# Adapted from https://stackoverflow.com/questions/31184103/how-to-improve-the-rendering-of-gradients-and-filled-elements-in-gnuplot
# to handle gnuplot's actual "set tables" output, with all its bs included.
# This takes those tables and produces an ascii stream that encodes parameters, which can be fed
# to "set object" to construct a filled polygon.  When you do this, the regions work fine.
# (Though you can still see some slight alpha munging across some graphs for some reason)
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

from __future__ import print_function
import sys
import math

f=open(sys.argv[1], 'r')
dats = []
for line in f:
    if line[0] == "#" or line[0] == "\n":
        continue

    # Otherwise, add it
    dats.append((line.strip()).split())

# Organize it (I'm sure theres a pythonic way to do this)
x = []
y = []
z = []
for dat in dats:
    if not len(dat) > 0:
        continue

    xm = float(dat[0])
    ym = float(dat[1])
    zm = float(dat[2])
    if not math.isnan(xm) and not math.isnan(ym) and not math.isnan(zm):
        x.append(xm)
        y.append(ym)
        z.append(zm)

# Coords?
if len(sys.argv) == 3:
    output = sys.argv[2]
    coords = ""
else:
    coords = sys.argv[2]
    output = sys.argv[3]

print('%s="' % output, end='')
for i in range(0,len(x)):
    if (i == 0):
        print('from %s {0},{1} '.format(x[i], y[i]) % coords, end='')
    else:
        print('to %s {0},{1} '.format(x[i], y[i]) % coords, end='')
for i in range(len(x)-1,-1,-1):
    print('to %s {0},{1} '.format(x[i], z[i]) % coords, end='')
print('"', end='')
