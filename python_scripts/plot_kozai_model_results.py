#============================================================
# Import modules
#============================================================
import re, os, fnmatch, sys, easygui as eg
from numpy import *
from matplotlib import * 
from pylab import *
import asciitable, string, datetime, shutil, atpy, random

def set_colors(NumColors=8):
    #    import matplotlib, random
    cc = [0] * NumColors
    for item in range(NumColors):
        cc[item] = ( (item) * 256 / NumColors)
    rcParams['axes.color_cycle'] = list(cm.hsv(cc))
    rcParams['lines.linewidth'] = 2
#============================================================

def plot_models(dir=''):
    # Declare constant values
    secPerYear = 3600.0 * 24.0 * 365.0
    rJup = 7.1492e9
    mJup = 1.89813e30
    rSun = 6.955e10
    LSun = 3.9e33

    if dir == '':
        # Select a directory/run to analyze
        savedirBase = eg.diropenbox(msg='Pick a run directory to examine',default='/Volumes/Data/Work/Research/BodenheimerCode/Code_for_Kozai/outputs/2013/')
        dir = [savedirBase]

		
    # Generate a list of all the models in that dir
    filelist = [ ]
    for file in os.listdir(dir):
        if fnmatch.fnmatch(file,'model_*.txt'):
            filelist.append(file)
            
    # Sort the file list by model number
    temp = []
    for item in filelist:
        foo = item.split('_')[-1]
        foo = foo.split('.')[0]
        temp.append(int(foo))
    temp.sort()
    i=0
    nummods = size(temp)
    for item in temp:
        i = i+1
        col = cm.spring(i*250/nummods)
        # Load the model summary data
        model = atpy.Table(dir+'/model_'+str(item)+'.txt',type='ascii')
	
        # Plot the data
        foo = 'model_'+str(item)
	
        plt.figure(5)
        plt.grid(True)
        plt.xlabel('Mass (Mjup)')
        plt.ylabel('Pressure')
        plot(model.M/mJup,model.P,label=foo,alpha=0.9,color=col)
        #        legend(loc="best",prop={'size':10},ncol=4)
	
        plt.figure(6)
        plt.grid(True)
        plt.xlabel('Mass (Mjup)')
        plt.ylabel('Radius (Rjup)')
        plot(model.M/mJup,model.R/rJup,label=foo,alpha=0.9,color=col)
        #        legend(loc="best",prop={'size':10},ncol=4)
        
        plt.figure(7)
        plt.grid(True)
        plt.xlabel('Mass (Mjup)')
        plt.ylabel('Luminosity (erg/s)')
        plot(model.M/mJup,model.L,label=foo,alpha=0.9,color=col)
        #        legend(loc="best",prop={'size':10},ncol=4)
        
	
        plt.figure(8)
        plt.grid(True)
        plt.xlabel('Mass (Mjup)')
        plt.ylabel('Temperature (K)')
        plot(model.M/mJup,model.T,label=foo,alpha=0.9,color=col)
        #        legend(loc="best",prop={'size':10},ncol=4)
	    
	
        #	show()

    return 0
