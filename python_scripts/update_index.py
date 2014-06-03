#============================================================
# Import modules
#============================================================
import re, os, fnmatch, sys, easygui as eg
from numpy import *
from matplotlib import * 
from pylab import *
import asciitable, string, datetime, shutil, atpy

#============================================================
# Define useful functions
#============================================================
def modification_date(filename):
    #    t = os.path.getmtime(filename)
    t = filename.split('/')
    index = 2 + t.index(filter(lambda x:x.startswith('201'),t)[0])
    t = t[index]
    return datetime.datetime.strptime(t,'%B_%d_%Y')
#============================================================
def extract_value(line,key,splitval):
    loc = line.find(key) + len(key)
    line = line[loc:]
    loc = line.find(splitval) + len(splitval)
    line = line[loc:]
    val = line.split()[0]
    return val
#============================================================
def get_update_list(base_path=''):
    if base_path == '':
        base_path = eg.diropenbox(msg='Select directory to crawl',default ='/Volumes/Data/Work/Research/BodenheimerCode/Code_for_Kozai/outputs/2013/')         
        #    result = []
    for (path,dirs,files) in os.walk(base_path):
        if 'full_run_output.txt' in files:
            update_index_main(path)
            #      result.append(path)    
    return 0

#============================================================
# Main function
#============================================================
def update_index_main(dirname=''):
    # Declare variables relevant to the index file
    path_to_index_file = '/Volumes/Data/Work/Research/BodenheimerCode/Code_for_Kozai/outputs/SIMULATION_INDEX.txt'
    header = ['DateCreated','Path','STARTING MODEL','BINARY OUTPUT FILE','MODA','NMOD','NRIT','ITMN','ITMX','JADD','JSUB','NATM','Atmx','Atmn','dTAX','L/H','dLmx','dLmn','dXmx','dXmn','dPmx','dPmn','Crad','Cwrk','dZmx','dZmn','dZdt','epsP','epsR','epsL','epsT','SminP','SminR','SminL','SminT','SmaxP','SmaxR','SmaxL','SmaxT','dTIM','FACT','dTMN','dTMX','CHMN','CHMX','XX','YY', 'JNOU', 'SIGM', 'KZPE', 'ZSTA', 'ECCN', 'SEMI', 'VALN', 'QVAL','Zmass','TMAX','TWRT']

    header_val = {}
    for item in header:
        header_val[item] = 'NaN'

    # Check if the index file exists.  If not, create it.
    if (os.path.isfile(path_to_index_file) == False):
        # Make the file
        index_file = open(path_to_index_file,'w')
        index_file.write(string.join(header,'\t'))
        index_file.write('\n')
    # Otherwise, open the existing index file in write/append mode
    else:
        index_file = open(path_to_index_file,'a')

    # If an input filename wasn't provided, prompt the user to select one:
    if (dirname == ''):
        # Ask the user which result to parse and add to the index
        name_of_file_to_add = eg.diropenbox(msg='Which run do you want to add to the index?',default ='/Volumes/Data/Work/Research/BodenheimerCode/Code_for_Kozai/outputs/') 
    else:
        name_of_file_to_add = dirname


    # Get the DateCreated and Path values for this file
    date_run_was_created = modification_date(name_of_file_to_add)
    ##########
    header_val['DateCreated'] = date_run_was_created.strftime('%Y-%m-%d')
    ##########
    run_path = name_of_file_to_add

    # Open the run file for parsing
    file_to_add = open(name_of_file_to_add+'/full_run_output.txt','r')
    header_val['Path'] = name_of_file_to_add
    stop_string = 'COMRD: END OF INPUT DATA'

    # Read the input params (file header, sort of), and parse them
    for line in file_to_add:
        # Ignore comment lines in the input file
        if (line.find('*') > -1):
            continue
        # Stop reading in the file once you reach the designated stopping point
        elif (line.find(stop_string) > -1):
            break
        else:
            #------------------------------------------------
            ## Special Parsing Cases ##
            #------------------------------------------------
            # If the line contains input or output file info, parse like this:
            if (line.find(header[2]) > -1 ):
                header_val[header[2]] = re.split(':\s+[0-9]+',line.strip())[-1].split()[-1]        
            elif (line.find(header[3]) > -1):
                header_val[header[3]] = re.split(':\s+[0-9]+',line.strip())[-1].split()[-1]
                
            # If the line contains SMIN or SMAX file info, parse like this:
            elif (line.find('SMIN') > -1):
                foo = line.split('=')[-1].split()
                for i in range(4):
                    header_val[header[31+i]] = foo[i]
            elif (line.find('SMAX') > -1):
                foo = line.split('=')[-1].split()
                for i in range(4):
                    header_val[header[35+i]] = foo[i]

            #------------------------------------------------
            ## Normal Parsing Case ##
            #------------------------------------------------
            else:
                for item in header:
                    if (line.find(item) > -1):
                        foo = extract_value(line,item,'=')
                        header_val[item] = foo

                
        # Close the input params file
    file_to_add.close()  

    # Get the dictionary entries in the right (header) order, and convert the values to a single line/string.
    foo = []
    for item in header:
        foo.append(header_val[item])
    foo = string.join(foo,'\t')

    # Write the run info to the index file
    index_file.write(foo+'\n')

    # Have something here that eliminates duplicates from the full_data index...?
    # May want to have a seperate program for doing that...

    # Close the index file
    index_file.close()
    return 0
