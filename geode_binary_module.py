#! /usr/bin/python3
#
# geode_binary_module.py
# Copyright(C) 2019 Kevin Croker
#
# Solves the orbital dynamics equations for DCOs, given the GEODE blueshift.
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

from scipy.integrate import odeint
import numpy as np
import sys
import math

# Assuming Planck 2018 flat cosmology
# and cosmological units.
OmegaM = 0.315
OmegaL = 1 - OmegaM

# CODATA G and c in mks, converted to MSol, AU, and H0^{-1}
# 67.4 km/(Mpc s) from Planck 2018
# H0Recip_mks = 1.0/(67.4 / 3.086e19)
H0recip_MKS = 4.5786e17 # seconds
G_MKS = 6.67408e-11
c_MKS = 2.998e8 # m/s
G_over_c2_MKS =  G_MKS/c_MKS**2  # m/kg
AU_per_meter = 1/1.496e11
kg_per_msol = 1.989e30 
secs_per_H0inv = 4.5786e17

# Things in astronomical units
G_over_c2_ASTRON = G_over_c2_MKS * AU_per_meter * kg_per_msol
c_ASTRON = c_MKS * AU_per_meter * secs_per_H0inv
G_ASTRON = G_over_c2_ASTRON * c_ASTRON**2

# So I don't compute them inline every time
G3_over_c5_ASTRON = (G_over_c2_ASTRON)**3 * c_ASTRON
sqrt_G7_over_c12_ASTRON = math.sqrt(G3_over_c5_ASTRON**2 * G_over_c2_ASTRON)

MyrToRecipH0 = (1.0/14516.9)

# We assumed the time unit here was H_0^{-1}, so this should be right
def H(a):
    return math.sqrt(OmegaM/a**3 + OmegaL)

# Rescale L so that the integrator does something
# (the angular momentum natural timescale is one orbital period, which is short af compared to a hubble time)
L_Rescale = 1e-10

# Make some global decls
m1 = None
m2 = None
a_start = None

no_merger_cutoff = 0.999

#
# Evolution equations due to GW emission are taken from Schutz, 1st ed. p. 250 Chapter 9, Exercise 46
# Originally from Peters 1964
#
# Since these equations are determined over the course of one orbit, the mass change due to blueshift can be treated
# adiabatically.  I.e. orbits are on the timescales of, at most 1000 of years (orbital timescales).
# The blueshift is on the order of millions of years (cosmological timescales).
#
# Since we do not track the actual orbit, we have to cutoff the merger when its going to happen in the next 10000 years or so.
# We do this by examining when periastron is < 1e-2.  
#
# (Note that R is used to indicate semi-major axis because a is used for RW scale factor)
#

class mergeException(Exception):
    pass

def geode_binary(state, a):
    L, R, e = state
    
    # Periastron occurs at R(1-e)
    # Stop the integration if R(1-e) < 10^-2 AU
    #
    # (This is robust to periastron or major axis as the cutoff)
    # Value chosen to keep the integrator from complaining and slowing down due to
    # unchecked stderr dumps
    if R*(1-e) < 1e-2 or a > no_merger_cutoff:
        # Leave by raising an exception.
        # This is my greatest hack ;)
        merge = mergeException()
        merge.state = (a, *state)
        raise merge
    else:
        # Follow L for sanity checks, even though we don't need to.
        dLda = -32.0/5 * (m1 * m2)**2*math.sqrt(m1 + m2)/(a * H(a) * R**(7.0/2)) * (a/a_start)**(15.0/2) * sqrt_G7_over_c12_ASTRON / (1-e**2)**2 * (1 + 7.0/8*e**2) * L_Rescale

        deda = (1/(H(a)*a) * (-304.0/15) * m1 * m2 * (m1 + m2) * (a/a_start)**9 * e/(R**4 * (1-e**2)**(5.0/2)))*(1 + 121.0/304 * e**2) * G3_over_c5_ASTRON        
        # dR gets two distinct effects: GW wave loss, GEODE blueshift
        # Note that this RAPIDLY transitions from blueshift dominated to GW dominated...

        # This amounts to the partial derivative holding e and L fixed, since we are dealing with the mass
        # The Schutz expression already considers that e and L are changing, with mass fixed :D
        # So we start this expression with derivatives of the masses.
        dRda = -9/a_start * ( (m1 + m2)/(m1*m2)**2 ) * (a_start/a)**10 * (L/L_Rescale)**2 / G_ASTRON / (1-e**2)

        # The 1/(1-e^2) in the denominator of the blueshift expression (from the usual expression for semi-major axis)
        # significantly shortens the lifetime
        dRda = dRda - 64.0/5 * m1*m2*(m1 + m2)/(a*H(a) * R**3) * (a/a_start)**9 * G3_over_c5_ASTRON / (1 - e**2)**(7.0/2) * (1 + 73.0/24*e**2 + 37.0/96 * e**4)
        
    # Return the values
    return [dLda, dRda, deda]

def geode_ns_binary(state, a):
    L, R, e = state
        
    if R*(1-e) < 1e-2 or a > no_merger_cutoff:
        # Leave by raising an exception.
        # super hax
        merge = mergeException()
        merge.state = (a, *state)
        raise merge
    else:
        # Notice that the heavier primary will always be the BH, based on the CE studies done by Dominik et al. (2012)
        dLda = -32.0/5 * (m2*m1*(a/a_start)**3)**2 * math.sqrt(m2 + m1*(a/a_start)**3)/(a * H(a) * R**(7.0/2)) * sqrt_G7_over_c12_ASTRON / (1-e**2)**2 * (1 + 7.0/8*e**2) * L_Rescale

        # "Eggshentrishitty" --Connery
        deda = (1/(H(a)*a)) * (-304.0/15) * e/( (1-e**2)**(5.0/2) * R**4) * (m2*m1*(a/a_start)**3) * (m1*(a/a_start)**3 + m2) * (1 + 121.0/304 * e**2) * G3_over_c5_ASTRON
        
        # dR gets two distinct effects: GW wave loss, GEODE blueshift
        # Note that this RAPIDLY transitions from blueshift dominated to GW dominated...
        dRda = -1.0/(a_start*(m1*m2)**2) * (3*m1*(a_start/a)**4 + 6*m2*(a_start/a)**7) * (L/L_Rescale)**2 / G_ASTRON / (1-e**2)
        dRda = dRda - 64.0/5 * (a/a_start)**3 * (m1*m2) * ( m1*(a/a_start)**3 + m2 ) / (a*H(a) * R**3) / (1 - e**2)**(7.0/2) * (1 + 73.0/24*e**2 + 37.0/96 * e**4) * G3_over_c5_ASTRON
        
    return [dLda, dRda, deda]

def bh_binary(state, a):
    L, R, e = state
    
    if R*(1-e) < 1e-2 or a > no_merger_cutoff:
        # Leave by raising an exception.
        merge = mergeException()
        merge.state = (a, *state)
        raise merge
    else:
        # These are the standard equations, expressed in scale factor
        deda = (1/(H(a)*a) * (-304.0/15) * m1 * m2 * (m1 + m2) * e/(R**4 * (1-e**2)**(5.0/2)))*(1 + 121.0/304 * e**2) * G3_over_c5_ASTRON 
        dLda = -32.0/5 * (m1 * m2)**2*math.sqrt(m1 + m2)/(a * H(a) * R**(7.0/2)) * sqrt_G7_over_c12_ASTRON / (1-e**2)**2 * (1 + 7.0/8*e**2) * L_Rescale
        dRda = - 64.0/5 * m1*m2*(m1 + m2)/(a*H(a) * R**3) / (1 - e**2)**(7.0/2) * (1 + 73.0/24*e**2 + 37.0/96 * e**4) * G3_over_c5_ASTRON
    
    # Return the values
    return [dLda, dRda, deda]

def evolve(alist, paramset, model):

    global m1, m2

    # Set the function parameters
    R_i, e_i, m1, m2  = paramset

    # Initial conditions
    init_state = [m1 * m2 * math.sqrt(G_ASTRON*R_i*(1-e_i**2)/(m1 + m2))*L_Rescale, R_i, e_i]
    
    # Integrate
    state = None
    try:
        odeint(model, init_state, alist)
    except mergeException as merged:
        state = merged.state

    # Return
    return state
