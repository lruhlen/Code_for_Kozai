def parse_kozai_results_on_cluster(input_filename, raw_output):
    #============================================================
    import datetime, os, fnmatch, sys, numpy, matplotlib, pylab, string,shutil
    #============================================================

    #============================================================
    # Convert the raw_output string to an array
    #============================================================
    raw_output = raw_output.split('\n')

    #============================================================
    # Determine where to store the results
    #============================================================
    date = datetime.datetime.today()
    month = date.strftime("%B")
    day = date.strftime("%d")
    year = date.strftime("%Y")
    secPerYear = 3600.0 * 24.0 * 365.0

    homedir = "/home/sgeadmin"
    savedirBase = "/home/sgeadmin/outputs/"+input_filename+"/"
    
    # Make the directory for our current version
    if (not os.path.exists(savedirBase)):
        os.makedirs(savedirBase)
    else:
        print 'woo' 
    
    # Copy the .mod file into the directory
    modfileName = '/home/sgeadmin/outputs/'+input_filename+'.mod'
    shutil.move(modfileName,savedirBase)

    # Move into the savedirBase folder
    os.chdir(savedirBase)

    InModel = 0
    InAtmos = 0
    currentTime = ''
    # model number, time, top of atmos radius, top of atmos luminosity, top of atmos temperature
    atmosOutFile = open('atmos_data.txt','w')
    atmosOutFile.write('ModelNum    Time    R     L     Teff\n')
    atmosOutFile.close()

    full_output_record = open('full_run_output.txt','w')
    
    # Go through the full run outputs file
    for line in raw_output:
        full_output_record.write(line)
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
    
    
    full_output_record.close()

    # Change back to the original directory    
    os.chdir(homedir)

    return len(raw_output)
