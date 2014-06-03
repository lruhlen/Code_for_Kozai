#============================================================
# Import modules
#============================================================
import re, os, fnmatch, sys, easygui as eg
import numpy, matplotlib, pylab
import asciitable, string, datetime, shutil, atpy
from update_index import *
#============================================================

def parse_results(SourceFile=''):
    #============================================================
    # Determine where to store the results
    #============================================================
    date = datetime.datetime.today()
    month = date.strftime("%B")
    day = date.strftime("%d")
    year = date.strftime("%Y")
    secPerYear = 3600.0 * 24.0 * 365.0

    homedir = os.getcwd()
    savedirBase = "/Volumes/Data/Work/Research/BodenheimerCode/Code_for_Kozai/outputs/"+year+"/"+month+"/"+month+"_"+day+"_"+year+"/"
    tmp = [0]
    versionNum=[0]

    if not os.path.exists(savedirBase): # If the save-directory does not already exist, create it
        os.makedirs(savedirBase)
    else:  # See how many versions, if any, have already been written here
        filelist = [ ]
        for file in os.listdir(savedirBase):
            if fnmatch.fnmatch(file,'*v[0-9]*'):
                filelist.append(file)
                # See how many versions already exist
            for name in filelist:
                prefix = string.split(name,'v')
                tmp.append(int(prefix[1]))
    
    # Figure out which version number we're writing to file, now
    versionNum = str(max(tmp)+1)
    
    # Make the directory for our current version
    savedirBase = savedirBase+"v"+versionNum+"/"
    os.makedirs(savedirBase)
    
    # Copy the full run results file to that directory
    if SourceFile == '':
        SourceFile =  eg.fileopenbox(msg='Select the file you want to parse',default='/Volumes/Data/Work/Research/BodenheimerCode/Code_for_Kozai/outputs/')

    shutil.move(SourceFile,savedirBase+'full_run_output.txt')

    # Move into the savedirBase folder
    os.chdir(savedirBase)

    InModel = 0
    InAtmos = 0
    currentTime = ''
    # model number, time, top of atmos radius, top of atmos luminosity, top of atmos temperature
    atmosOutFile = open('atmos_data.txt','w')
    atmosOutFile.write('ModelNum    Time    R     L     Teff\n')
    atmosOutFile.close()
    
    # Go through the full run outputs file
    infile = open("full_run_output.txt",'r')
    for line in infile:
        # Skip over blank lines
        if (line == "\n"):
            continue  
        # Find the kozai timescale for this run
        if (line.find('KZPE') != -1):
            tKozai = float(line.split()[1]) * secPerYear
        # Find the start of each model
        if (line.find('Start of model') != -1):
            InModel = 1
            # Figure out which model number this one is
            modelNum = line.split()[-1]
            modelOutFile = open('model_'+modelNum+'.txt','w')
            continue
        # Find the end of each model
        if (line.find('End of model') != -1):
            InModel = 0
            modelOutFile.close()
            continue
        # Find the start of each atmos
        if (line.find('Start of atmos') != -1):
            InAtmos = 1
            atmosNum = line.split()[-1]
            atmosOutFile = open('atmos_data.txt','a')
            continue
        # Find the end of each atmos
        if (line.find('End of atmos') != -1):
            InAtmos=0
            atmosOutFile.close()
            continue
    
        # Write the model data to its own file
        if (InModel == 1):
            if (line.find('TIME:') != -1):
                currentTime= line.split('TIME:')[-1].split()[0]
                continue
            else:
                # Strip out the 'are we convecting?' asterisks
                line = line.replace('*','')
                modelOutFile.write(line)
    
        # Write the atmos data to... its own file?  Or collate results in a single file?
        if (InAtmos == 1) and (line.find('Top') != -1):
            # Manipulate/parse the line string to extract/refashion its content for
            # printing to the atmos output file
            foo = line.split(":")[-1].split()
            Teff = foo[1]
            R = foo[3]
            L = foo[5]
            # We want the following info, in this order:
            # model number, time, top of atmos radius, top of atmos luminosity, top of atmos temperature
            currentTime.replace('D','E')
            print 'currentTime is ',currentTime
            currentTime = str(float(currentTime)/tKozai)
            print 'currentTime is ',currentTime            
            print ' '
            line = atmosNum+'    '+currentTime+'    '+R+'     '+L+'      '+Teff+'\n'
    
            # Write the line to file
            atmosOutFile.write(line)
    
    
    # Change back to the original directory    
    os.chdir(homedir)
    
    # Update the simulation index with this run's info
    #run update_index.py
    print 'Updating the index file'
    foo = update_index_main(savedirBase)

    return 0
