#============================================================
# Import modules
#============================================================
#import re, os, fnmatch, sys, easygui as eg
#from numpy import *
#from matplotlib import * 
#from pylab import *
#import asciitable, string, datetime, shutil, atpy, random

#============================================================
def plot_batch(data,color_num=0): # Note: 'data' has to be a dataframe type (as in, what search_index.py loads into the selected_data variable)

    # Declare constant values
    secPerYear = 3600.0 * 24.0 * 365.0
    rJup = 7.1492e9 
    rSun = 6.955e10
    LSun = 3.9e33
    mJup = 1.89813e30
    AU = 1.49597871e13

    msize = 6
    
    color_list = ['black','darkviolet','aquamarine','green','slategray','yellow','darkorange','magenta','red','brown','indigo','lightcoral','olivedrab','slateblue','mediumvioletred','gold','firebrick','dodgerblue','deeppink','darkorchid','darkcyan']
    batch_color = color_list[mod(color_num,len(color_list))]


    files = data.Path.values
    files = files.tolist()

    plt.figure(1)

        
    for item in files:
        if len(os.listdir(item)) <= 2:
            print item+' has no (or few) converged models, so it is not getting plotted.'
            continue
        if not os.path.exists(item+"/atmos_data.txt"):
            print item+"/atmos_data.txt does not exist.  Skipping to next file."
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
                break
        infile.close()

        titlestr=""
	
	# Load the atmos summary data
        atmos = atpy.Table(item+'/atmos_data.txt',type='ascii')

	
        # Plot the atmos data
        plot(atmos.Time * float(tKozai),atmos.R/rJup,'-',color=batch_color)
        #        plt.legend(loc="best",prop={'size':10},ncol=3)


    # Annotate/legend the batch as a whole
    plt.grid(True,which="both")
    plt.xlabel('Years')
    plt.xscale('log')
    plt.ylabel('Radius of tau = 2/3 surface (Rjup)')
    plt.title(titlestr)

    plt_label = '~'+str(round(data.Zmass[0]/mJup,3)) + ' Mjup'
    plot(None,'-',color=batch_color,label=plt_label)
    plt.legend(loc="best",prop={'size':10})
        
    show()
    return 0
