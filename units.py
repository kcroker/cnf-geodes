#!/usr/bin/python3
#
# units.py
# Copyright(C) 2019 Kevin Croker
#
# Utility class, easily updatable.  Not that we had this bright idea early enough...
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


import math

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
Rsol_per_AU = 1/0.00464913034
AU_per_Rsol = 0.00464913034
m_per_Mpc = 3.086e22

# Check this one...
MyrToRecipH0 = (1.0/14516.9)

# Things in astronomical units
G_over_c2_ASTRON = G_over_c2_MKS * AU_per_meter * kg_per_msol
c_ASTRON = c_MKS * AU_per_meter * secs_per_H0inv
G_ASTRON = G_over_c2_ASTRON * c_ASTRON**2

# So I don't compute them inline every time
G3_over_c5_ASTRON = (G_over_c2_ASTRON)**3 * c_ASTRON
sqrt_G7_over_c12_ASTRON = math.sqrt(G3_over_c5_ASTRON**2 * G_over_c2_ASTRON)
