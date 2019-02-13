#!/usr/bin/python3
#
# Cuts mergers on the ligo horizon, adjusted for unequal mass ratio
#
import ligo_chirp_cut
import sys
n = 0
for datum in sys.stdin:

    try:
        values = datum.split()
        
        if float(values[0]) > ligo_chirp_cut.horizon(float(values[3]), float(values[4]), float(values[0])):
            continue
        n = n + 1
        print(datum, end='')
    except Exception as e:
        print(e, file=sys.stderr)

# Print out the total
print("# %d" % n)
