#!/usr/bin/python3
#
# fastbin.py
# Copyright(C) 2019 Kevin Croker
#
# Just some dirty to make histograms! 
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

import numpy as np
import sys

def dumpbinned(var, bins, f):
    if bins is None:
        hist, bin_bds = np.histogram(var)
    else:
        hist,bin_bds = np.histogram(var, bins)

    print(bin_bds, file=sys.stderr)
    binrs = bin_bds[1:]
    binl = bin_bds[0]

    # Normalize by the number of things to bin
    print("# N = %d" % (len(var)), file=f)
    for val,binr in zip(hist,binrs):
        print("%.15e %.15e" % (binl + (binr - binl)/2.0, val/len(var)), file=f)
        binl = binr
    print("", file=f)

def dumpbinned_raw(var, bins, f):
    if bins is None:
        hist, bin_bds = np.histogram(var)
    else:
        hist,bin_bds = np.histogram(var, bins)

    print(bin_bds, file=sys.stderr)
    binrs = bin_bds[1:]
    binl = bin_bds[0]

    # # Normalize by the number of things to bin
    # print("# N = %d" % (len(var)), file=f)
    # for val,binr in zip(hist,binrs):
    #     print("%.15e %.15e" % (binl + (binr - binl)/2.0, val/len(var)), file=f)
    #     binl = binr
    # print("", file=f)

    #    Unnormalized
    print("# N = %d" % (len(var)), file=f)
    binl = bin_bds[0]
    for val,binr in zip(hist,binrs):
        print("%.15e %.15e" % (binl + (binr - binl)/2.0, val), file=f)
        binl = binr
    print("", file=f)

