###############################################
# NOTE: YOU NEED TO COPY AND PASTE THESE
# FUNCTIONS DIRECTLY INTO THE PYTHON TERMINAL
# WINDOW.
#
# Typing 'run batch_commands.py'
# doesn't work because of some sort of issue
# with importing and using the os module.
###############################################

#----------------------------------------------
from parse_kozai_results import parse_results
from multiprocessing import Pool
#from os import *
import os

#----------------------------------------------
def batch_parse(num_runs=0):
    if num_runs == 0:
        print 'You fool!  You need to tell me how many runs there are to parse!'
        return 1
    else:
        batch_output_dir = '/Volumes/Data/Work/Research/BodenheimerCode/Code_for_Kozai/outputs/batch_outputs/'
        for i in range(1,num_runs+1,1):
            file_to_parse = batch_output_dir+'run'+str(i)+'.txt'
            parse_results(file_to_parse)
            print 'Parsed ',file_to_parse
        return 0
    
#----------------------------------------------
def simulation_syscall(input_filename=''):
    #    os.system("echo "+foo)
    homedir = os.getcwd()
    outdir = '/Volumes/Data/Work/Research/BodenheimerCode/Code_for_Kozai/outputs/batch_outputs/'
    indir = '/Volumes/Data/Work/Research/BodenheimerCode/Code_for_Kozai/inputs/batch_inputs/'
    rundir = '/Volumes/Data/Work/Research/BodenheimerCode/Code_for_Kozai'

    filename_base = input_filename.split('/')[-1]
    output_filename = outdir+filename_base
    input_filename = indir+filename_base

    os.chdir(rundir)
    #    os.system("./full_tidal_heating < "+input_filename+" > "+output_filename)
    os.system("./full_tidal_heating_time_limited < "+input_filename+" > "+output_filename)
    #    os.system("growlnotify -m  \" "+output_filename + " completed\" -s \"Background script notification\" &")
    print 'Done running '+input_filename+" simulation."
    foo = output_filename.split('.txt')[-1]+'.mod'
    os.system("rm "+input_filename+" "+foo)
    os.chdir(homedir)
    print filename_base
    
    return filename_base

#----------------------------------------------
def batch_run(num_proc=6):
    pool = Pool(num_proc)
    indir = '/Volumes/Data/Work/Research/BodenheimerCode/Code_for_Kozai/inputs/batch_inputs/'

    files = os.listdir(indir)
    pool.map(simulation_syscall,files)
    pool.close()
    pool.join()

    num = len(files)
    batch_parse(num)

    
    #   print 'Done with this simulation batch run!'
    return 0
