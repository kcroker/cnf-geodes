#!/usr/bin/python3

import fastbin
import numpy as np
import sys
import pdb


import ligo_chirp_cut

# sys.argv[1] is a string of columns to stack, ONE INDEXED
# sys.argv[2] is a string bin start, bin end, bin count
cols = list(map(int, sys.argv[1].split()))
start,end,step = tuple(map(float, sys.argv[2].split()))

norm = True
if len(sys.argv) > 3:
    nonorm = False

def bfloat(f):
    try:
        return float(f)
    except:
        return float('nan')
    
dats = []
for n,line in enumerate(sys.stdin):
    datum = line.split()
    try:
    
        if len(datum) < 1:
            continue
    
        if datum[0].startswith("#"):
            continue

        # Turn into floats
        datum = [bfloat(x) for x in datum]

        for col in cols:
            dats.append(datum[col-1])

    except Exception as e:
        print(e, file=sys.stderr)
    
bins = np.arange(start,end,step)
if norm:
    fastbin.dumpbinned(dats, bins, sys.stdout)
else:
    fastbin.dumpbinned_raw(dats, bins, sys.stdout)
    
# Now do mass sum
