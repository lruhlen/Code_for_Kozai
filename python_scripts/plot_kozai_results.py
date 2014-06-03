#============================================================
# Import modules
#============================================================
import re, os, fnmatch, sys, easygui as eg
from numpy import *
from matplotlib import * 
from pylab import *
import asciitable, string, datetime, shutil, atpy

# Declare constant values
secPerYear = 3600.0 * 24.0 * 365.0
rJup = 7.1492e9 
rSun = 6.955e10
LSun = 3.9e33

# Select a directory/run to analyze
savedirBase = eg.diropenbox(msg='Pick a run directory to examine')

# Load the atmos summary data
atmos = atpy.Table(savedirBase+'/atmos_data.txt',type='ascii')

# Plot the atmos data
plt.figure(1)
plot(atmos.Time/secPerYear,atmos.R/rJup,'o',markersize=10)
plt.xlabel('Time (years)')
plt.ylabel('Radius of tau = 2/3 surface (Rjup)')
show()

plt.figure(2)
plot(atmos.Teff,atmos.L/LSun,'o',markersize=10)
plt.xlabel('Teff (K)')
plt.xlim([1e3,0.0])
plt.ylabel('Luminsity (Lsun)')
plt.yscale('log')
show()
