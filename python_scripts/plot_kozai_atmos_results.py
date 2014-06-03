#============================================================
# Import modules
#============================================================
import re, os, fnmatch, sys, easygui as eg
from numpy import *
from matplotlib import * 
from pylab import *
import numpy, matplotlib, pylab
import asciitable, string, datetime, shutil, atpy, random

def set_colors(NumColors=8):
    #    import matplotlib, random
    #rcParams['lines.linewidth'] = 2
    cc = [0] * NumColors
    neg_colors = ['blue','green','pink','orange','black','red','BlueViolet']
    
    if NumColors > 0:
        for item in range(NumColors):
            cc[item] = ( (item) * 256 / NumColors)
            rcParams['axes.color_cycle'] = list(cm.hsv(cc))
    else:
        print 'You FOOL!  You must pass a POSITIVE number to set_colors!'


#============================================================
def plot_atmos(files=''):
    # Declare constant values
    secPerYear = 3600.0 * 24.0 * 365.0
    rJup = 7.1492e9 
    rSun = 6.955e10
    LSun = 3.9e33
    mJup = 1.89813e30
    msize = 6

    if files == '':
        # Select a directory/run to analyze
        savedirBase = eg.diropenbox(msg='Pick a run directory to examine',default='/Volumes/Data/Work/Research/BodenheimerCode/Code_for_Kozai/outputs/2013/')
        files = [savedirBase+'/']
    
    for item in files:
        #############################################################
        # Make sure the run actually has some output models to plot...
        # (start here where I get back)
        #############################################################
        if len(os.listdir(item)) <= 3:
            print item+' has no (or few) converged models, so it is not getting plotted.'
            continue

        
        # Find the kozai timescale for this run
        infile = open(item+"/full_run_output.txt",'r')
        for line in infile:
            if (line.find('EFAC') != -1):
                efac = line.split()[-1]
                sigmasquared = line.split()[-3]
                jnought = line.split()[-5]
            if (line.find('KZPE') != -1):
                tKozai = line.split()[-1]                 
            if (line.find('Zmass') != -1):
                Zmass = line.split()[-1]
                Zmass = round( (float(Zmass)/mJup),3)
                break
        infile.close()

        #pltlabel = savedirBase.split('/')[-1]
        #pltlabel = eg.enterbox(msg='Enter the Kozai timescale for this run')
        #pltlabel = tKozai+' yrs'
        #pltlabel = efac
        #pltlabel = str(sqrt(float(sigmasquared)))
        #pltlabel = jnought
        titlestr=""
        #titlestr = eg.enterbox(msg='What do you want to title this set of plots?')
        #pltlabel = eg.enterbox(msg="Enter identifying info about this data")
        #        pltlabel = string.join(item.split('/')[-3:],'/')
        pltlabel = str(Zmass)+' Mjup'

	# Load the atmos summary data
        atmos = atpy.Table(item+'/atmos_data.txt',type='ascii')

	
        # Plot the atmos data
        #plt.close("all")
        plt.figure(1)
        plot(atmos.Time * float(tKozai),atmos.R/rJup,'-o',markersize=msize,label=pltlabel,alpha=0.75)
        plt.grid(True,which="both")
        plt.xlabel('Years')
        plt.xscale('log')
        plt.ylabel('Radius of tau = 2/3 surface (Rjup)')
        plt.title(titlestr)
        plt.legend(loc="best",prop={'size':10},ncol=1)
        plt.ylim([0,4])
        plt.xlim([0,2e8])
        
#        plt.figure(2)
#        plot(atmos.Time* float(tKozai),atmos.L,'-o',markersize=msize,label=pltlabel,alpha=0.75)
#        plt.grid(True,which="both")
#        plt.xlabel('Years')
#        plt.xscale("linear")
#        plt.ylabel('erg/s')
#        plt.yscale('log')
#        plt.title(titlestr)
#        plt.legend(loc="best",prop={'size':10},ncol=3)
#       
#        plt.figure(3)
#        plot(atmos.Time* float(tKozai),atmos.Teff,'-o',markersize=msize,label=pltlabel,alpha=0.75)
#        plt.grid(True,which="both")
#        plt.xlabel('Years')
#        plt.yscale('log')
#        plt.ylabel('Teff (K)')
#        plt.title(titlestr)
#        plt.legend(loc="best",prop={'size':10},ncol=3)
#       
#        plt.figure(4)
#        plot(atmos.Teff,atmos.L/LSun,'o',markersize=8,label=pltlabel,alpha=0.5)
#        plt.grid(True,which="both")
#        plt.xlabel('Teff (K)')
#        plt.xlim([4e3,0.0])
#        plt.ylabel('Luminsity (Lsun)')
#        plt.yscale('log')
#        plt.title(titlestr)
#        plt.legend(loc="best",prop={'size':10},ncol=3)
#
    show()
    return 0
