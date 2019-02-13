#!/usr/bin/python3
#
# units.py
# Copyright(C) 2019 Kevin Croker
#
# Python "library" for table backed data.
# This is not meant to be memory efficient at all, but to be very easy to read.
# It was also a practice for me using python classes.
# At first, directly playing with cumtrapz felt awkward, but now its standard.
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
# 
import math
import numpy as np
import datetime
from scipy import integrate, misc, interpolate
import sys

#
# Current version only does uniform-spacing backed tables
#
# It shouldn't be too hard to add internal tracking of the numerical uncertainty
# in the backed function
#
# Current version assumes that the point locations coincide!
#
# Assume all functions have the exact same interpolation
#
def divcheck(numer, denom):
    if denom == 0:
        if numer == 0:
            return 0
        else:
            # derp
            return 87
    else:
        return numer/denom
   
class TblBackedFunction:

    def addComment(self, comment):
        # Append a thing to be used as a comment whenever writeout() is called
        self.comments.append(comment)

    def clearComments(self):
        # Ref an empty list
        comments = []
        
    def writeout(self, backing_file=""):
        # Attempt to open
        try:
            if not backing_file:
                # Assume they want a dump to stdio
                self.fhandle = sys.stdout
            else:
                self.fhandle = open(backing_file, 'w')
        except IOError as oops:
            print("Failed to open backing file " + backing_file)
            raise oops
        try:
            
            # First, write the comments as a header
            for comment in self.comments:
                self.fhandle.write("# %s\n" % comment)

            # Now write out the data
            for [xdat, ydat] in zip(self.x, self.y):
                self.fhandle.write("%.15e\t%.15e\n" % (xdat, ydat))

        except IOError as oops:
            print("Failed to write backing file %s. (Error text follows)" % backing_file, oops)

    def __mul__(self, other):
        tmpx = [x for x in self.x]
        if isinstance(other, self.__class__):
            N = range(len(self))
            tmpy = [self.y[n] * other.y[n] for n in N]
            return TblBackedFunction(tmpx, tmpy)
        else:
            tmpy = [q * other for q in self.y]
            return TblBackedFunction(tmpx, tmpy)

    
    # Define commutative multiplication
    __rmul__ = __mul__

    # Define convenient negation
    def __neg__(self):
        return -1.0*self

    def __sub__(self, other):
        N = range(len(self))
        tmpx = [self.x[n] for n in N]
        if isinstance(other, self.__class__):
            tmpy = [self.y[n] - other.y[n] for n in N]
            return TblBackedFunction(tmpx, tmpy)
        else:
            tmpy = [q - other for q in self.y]
            return TblBackedFunction(tmpx, tmpy)
            
    # Define commutative subtraction
    __rsub__ = __sub__

    def __truediv__(self, other):
        tmpx = [x for x in self.x]
        
        if isinstance(other, self.__class__):
            N = range(len(self))

            # Check for 0/0 at least
            tmpy = [divcheck(self.y[n],other.y[n]) for n in N]
            return TblBackedFunction(tmpx, tmpy)
        else:
            tmpy = [q / float(other) for q in self.y]
            return TblBackedFunction(tmpx, tmpy)

    # Define "commutative" division
    def __rtruediv__(self, atom):
        tmpx = [x for x in self.x]

        # Dunno if this will ever fire due to evaluation order
        if isinstance(atom, self.__class__):
            return __truediv__(self, atom)
        tmpy = [float(atom) / q for q in self.y]
        return TblBackedFunction(tmpx, tmpy)

    # Define addition
    def __add__(self, other):
        tmpx = [x for x in self.x]
        
        if isinstance(other, self.__class__):
            N = range(len(self))
            tmpy = [self.y[n] + other.y[n] for n in N]
            return TblBackedFunction(tmpx, tmpy)
        else:
            tmpy = [q + other for q in self.y]
            return TblBackedFunction(tmpx, tmpy)
            
    # Define commutative addition
    __radd__ = __add__

    # Define length 
    def __len__(self):
        return self.mylen

    def where(self, c):
        guess = math.floor((c-self.x[0])/self.resolution)
        if guess == len(self):
            guess -= 1

        # It will either be to the right or left of guess
        if self.x[guess] <= c:
            return guess
        elif self.x[guess] > c:
            return guess - 1
        else:
            return guess + 1

    def iwhere(self, n):
        return self.x[n]

    # Evaluation for arbitrary c via Neville interpolation of 5 bracketing points 
    def __call__(self, c):
        # Check bounds
        if c < self.xmin or c > self.xmax:
            raise Exception("Requested value not in range of backing table: ", c)

        where = self.where(c)

        if self.interpolator is None:
            self.interpolator = interpolate.interp1d(self.x, self.y)

        return self.interpolator(c)

    # Switch x and y
    def invert(self):
        tmp = [x for x in self.x]
        self.x = self.y
        self.y = tmp

        # Reset things
        self.xmin = min(self.x)
        self.xmax = max(self.x)
        self.ymin = min(self.y)
        self.ymax = max(self.y)
        self.mylen = len(self.x)
        
        # Subtract 0.5 so that we always floor right
        self.resolution = (self.xmax - self.xmin)/self.mylen

        # print("Data covers (", self.xmin, ", ", self.xmax, ") --> (", self.ymin, ", ", self.ymax, ") with resolution ", self.resolution)

    # Deep copy
    def copy(self):
        return TblBackedFunction(self.x, self.y)

    def integrateTop(self, option="derp"):
        return indefinite_integral_top(self, option).y[-1]

    def integrateBottom(self, option="derp"):
        return indefinite_integral_bottom(self, option).y[-1]

    def derivative(self):
        return derivative(self)
    
    # Constructor
    def __init__(self, backing_file, xcol=0, ycol=1, nonuniform=False):

        # Initialize the comments list
        self.comments = []
        self.addComment("Creation date (UTC): %s" % datetime.datetime.utcnow())
        
        # If its a np array, make it into a list real fast
        if isinstance(backing_file, np.ndarray):
            backing_file = [x for x in backing_file]

        # Make an interpolator
        self.interpolator = None
        
        # Is backing_file a list?
        if isinstance(backing_file, list):
            self.backing_file = "(lists)"

            # Alias for clarity
            x = backing_file
            y = xcol

            # Sanity check
            if not len(x) == len(y):
                raise Exception("Backing lists do not have the same length: ", len(x), len(y))

            # Assign the tables
            self.x = x
            self.y = y

            # XXX Add check to determine non-uniformity
            # and does neville work on non-uniform samples?
            
        elif isinstance(backing_file, str):
            self.backing_file = backing_file
            self.xcol = xcol
            self.ycol = ycol

            # Attempt to open
            try:
                self.fhandle = open(backing_file, 'r')
            except IOError as oops:
                print("Failed to open backing file " + backing_file)
                raise oops

            # Add a comment that we came from this file
            self.addComment("Data backed by: %s" % backing_file)
            
            # Backing file is open, parse it down
            data = []
            for line in self.fhandle:
                if line[0] == "#" or line[0] == "\n":
                    continue

                # Otherwise, add it
                data.append((line.strip()).split())
                
            # XXX Sorting was broken, fix it
            # XXX can't handle reverse order domain....
            
            # Make data
            self.x = []
            self.y = []

            # Add the data, skipping nans
            for datum in data:
                xm = float(datum[self.xcol])
                ym = float(datum[self.ycol])

                if not math.isnan(xm) and not math.isnan(ym):
                    self.x.append(xm)
                    self.y.append(ym)
        
            # Sort the data
            sorted_combo = sorted(zip(self.x, self.y), key=lambda q : q[0])

            for n,derp in enumerate(sorted_combo):
                self.x[n] = derp[0]
                self.y[n] = derp[1]
            
            # Close the backing file
            self.fhandle.close()
        else:
            raise Exception("Unsupported data spec type: ", type(backing_file))
        
        # Determine maximum, minimum, and resolution
        self.xmin = min(self.x)
        self.xmax = max(self.x)
        self.ymin = min(self.y)
        self.ymax = max(self.y)
        self.mylen = len(self.x)
        
        # Subtract 0.5 so that we always floor right
        self.resolution = (self.xmax - self.xmin)/self.mylen

        # print("Data covers (", self.xmin, ", ", self.xmax, ") --> (", self.ymin, ", ", self.ymax, ") with resolution ", self.resolution)

        # Verify that the data is uniform, and the step everywhere is equal to the 
        # resolution
        # ...

### Integration and differentiation routines 

#
# Right now, we need integral bottom
#
#   cumtrapz.
# 
# Compute all the cumulative integrals are tblf 
def indefinite_integral_bottom(tblf, upperbound="derp"):

    if upperbound == "derp":
        maxn = len(tblf.x) - 1
    else:
        maxn = tblf.where(upperbound)
        
    subx = (tblf.x[:maxn])[::-1]
    suby = [-1.0*z for z in (tblf.y[:maxn])[::-1]]

    cum = integrate.cumtrapz(suby, subx, initial=0)

    # Now return the list in correct order
    return TblBackedFunction(subx[::-1], cum[::-1])

def indefinite_integral_top(tblf, lowerbound="derp"):

    # Check for omitted argument
    if lowerbound == "derp":
        minn = 0
    elif lowerbound < tblf.x[0]:
        raise Exception("Lowerbound not included in the domain of the table")
    else:
        minn = tblf.where(lowerbound)
        
    subx = tblf.x[minn:]
    suby = tblf.y[minn:]

    cum = integrate.cumtrapz(suby, subx, initial=0)

    # Now return the list in correct order
    return TblBackedFunction(subx, cum)

# Remove NaN's for visualization purposes... use very carefully...
def stripNansInfs(tblf):
    newlist = []
    nanlist = []
    for datum in zip(tblf.x, tblf.y):
        if not math.isnan(datum[1]) and not math.isinf(datum[1]):
            newlist.append(datum)
        else:
            nanlist.append(datum)

    newx,newy = list(zip(*newlist))
    stripped = TblBackedFunction(list(newx), list(newy))
    stripped.addComment("NaN's were removed at: ")
    for nan in nanlist:
        stripped.addComment("%f %f" % (nan[0], nan[1]))

    return stripped

from scipy.signal import butter, lfilter, freqz

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def lowpass(tblf, cutoff, fs, order):
    xform = butter_lowpass_filter(tblf.y, cutoff, fs, order)
    return TblBackedFunction(tblf.x, xform)

    
# This should eventually be replaced with this hot smoothing thing that I use
# in downstream stuff...
def derivative(tblf):
    # Based on toth's awesome
    # x0=NaN
    # y0=NaN
    # plot 'test.dat' using (dx=$1-x0,x0=$1,$1-dx/2):(dy=$2-y0,y0=$2,dy/dx) w l notitle
    #
    # because np.gradient is fux0r3d
    
    x0 = tblf.x[0]
    y0 = tblf.y[0]
    newx = []
    newy = []
    
    for x,y in zip(tblf.x[1:], tblf.y[1:]):
        dx = x - x0
        x0 = x
        newx.append(x - dx/2.0)

        dy = y - y0
        y0 = y
        newy.append(dy/dx)
        
    return TblBackedFunction(newx, newy)

# UNIT TESTING STUB
def paulstubbs():

    f1 = TblBackedFunction("dStarda.dat", nonuniform=False)
    f1.writeout("idempotency_test.dat")

    for n in range(len(f1)):
        print("Backing list invertability N -> R -> N TEST: n = %d, actual doman = %f, inverse = %d" % (n, f1.iwhere(n), f1.where(f1.iwhere(n))))
        if not n == f1.where(f1.iwhere(n)):
            raise Exception("N->R->N invertability broken")

    print("PASSED N -> R -> N invertability test for backing list")

    for k in np.linspace(f1.xmin, f1.xmax, len(f1)*3):
        if f1.where(k) == len(f1) - 1:
            pass
        else:
            print("Backing list invertability R -> N -> R TEST: %.10e <= %.10e < %.10e" % (f1.iwhere(f1.where(k)), k, f1.iwhere(f1.where(k) + 1)))

        if f1.iwhere(f1.where(k)) <= k and k < f1.iwhere(f1.where(k) + 1):
            pass
        else:
            raise Exception("R->N->R invertability broken")
            
    print ("PASSED R -> N -> R invertability test for backing list")

    # Verify + invertability for equal resolution objects
    ftest = f1 - f1
    ftest.writeout("add_invert_check.dat")

    # Verify * invertability for equal resolution objects directly
    ftest = f1/f1
    ftest.writeout("mult_invert_check.dat")

    # Verify * invertability for atomic operations
    ftest = 10.0*f1/10.0 - f1
    ftest.writeout("atom_mult_invert_check.dat")

    # Verify + invertability for atomic operations
    ftest = 10.0+f1 - 10.0
    ftest.writeout("atom_add_invert_check.dat")

    # Verify division OF atoms
    ftest = 10/f1
    ftest.writeout("div_of_atom_check.dat")

    ftest = indefinite_integral_bottom(f1, 1)
    ftest.writeout("integral_of_dStarda.dat")
    
    ftest2 = -derivative(ftest)
    ftest2.writeout("Dintegral_of_dStarda.dat")

    # Verify the fundamental theorem of calculus
    ftest2 = tc.derivative(f1)
    dStarda_again = tc.indefinite_integral_top(ftest2, ftest2.xmin)
    dStarda_again.writeout("ftc_test.dat")

# # And write some more when you don't have a thesis to finish


#
# Note that the integrator will try to work with the underlying backing data
# 
# definite_integral_top(f, anchor, x) will then give another table backed function 
#   which is the integral of f, from anchor to x (the variable).
#
# definite_integral_bottom(f, x, anchor) does what you expect
#
# derivative(f) will give another backed function which is the derivative of f
#
# define all the operations to be pointwise for the object too.
#
        
        
