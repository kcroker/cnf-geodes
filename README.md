# cnf-geodes
Simulation, analysis, and visualization pipeline for "Implications of symmetry and pressure in Friedmann cosmology. III. The black hole mass function"

## Overview

This is a set of small analysis tools, which grew organically during the development of Paper III.
I have attempted to be as Pythonic as possible, especially when writing in bash ;)
But seriously, please be forgiving with respect to inefficiencies in the implementation.
I have also attempted to carve out all vestigial code, but may have missed things here or there.
In particular, `units.py` contains lots of stuff that is duplicated elsewhere (but everything should be consistent).

The tools split into three categories:
1.) Computation: e.g. a model GEODE theoretical mass function, a detector sensitivity cut
2.) Solving the Keplerian orbital dynamics of binary systems, possibly with GEODEs
3.) Auxiliary scripts which process outputs from 1.) and 2.) into the visualizations present in Paper III

Most tools are written in the UNIX style: the accept line based input on STDIN and produce line based output on STDOUT, complaining to STDERR.

## Data Products

All data products are space-delimited ASCII tables to double precision.
Where not obvious, the format and contents of the tables are given as a comment in the table headers.
Units for the main data outputs are
+ **Mass**: Solar mass
+ **Distance**: AU

## Usage



