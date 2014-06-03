#============================================================
# Import modules
#============================================================
import re, os, fnmatch, sys, easygui as eg
import pandas as pd
import asciitable, string, datetime, shutil, atpy, time
from numpy import *
from matplotlib import * 
from pylab import *
from datetime import date

def search_index():
    #============================================================
    # Define constants
    #============================================================
    mJup = 1.89813e30
    AU = 1.49597871e13

    #============================================================
    # Read in the index file (data) as (to) a python table
    #============================================================
    path_to_index_file = '/Volumes/Data/Work/Research/BodenheimerCode/Code_for_Kozai/outputs/SIMULATION_INDEX.txt'
    #index_data = atpy.Table(path_to_index_file,type='ascii',masked=True)
    index_data = atpy.Table(path_to_index_file,type='ascii')
    dates = pd.to_datetime(index_data.DateCreated)

    #============================================================
    # Get a list of all the parameter/column names available to sort on
    #============================================================
    param_list = index_data.columns.keys
    kozai_indices = (0,47,48,49,50,51,52,54,55,56,57)
    kozai_params = map(param_list.__getitem__,kozai_indices)
    kozai_params.append('MinNumResults')

    other_var_indices = []
    for i in range(len(param_list)):
        if i not in kozai_indices:
            other_var_indices.append(i)
    other_params = map(param_list.__getitem__,other_var_indices)

    #============================================================
    ## Load the index data into a data frame (for database querying/searching later on)
    #============================================================
    full_data = pd.DataFrame(index_data[::],index=dates)
    full_data = full_data.sort()

    #============================================================
    # Create the 'selected data' copy of the full data dataframe.
    #============================================================
    selected_data = full_data

    #============================================================
    # Open a GUI window that displays all the params available to sort on
    #============================================================
    sort_params1 = []
    sort_params2=[]
    sort_params1 = eg.multchoicebox(msg='Pick the KOZAI params to sort on',choices = kozai_params)
    #sort_params2 = eg.multchoicebox(msg='Pick the other params to sort on',choices = other_params)
    sort_params = sort_params1 + sort_params2

    #============================================================
    # Ask the user for min/max or 'contains string' criteria on all the selected fields
    #============================================================
    for item in range(len(sort_params)):
        
        if sort_params[item] == 'Path':
            #        thing = raw_input('File path contains: ')
            thing = eg.enterbox('File path contains: ')
            selected_data = selected_data[selected_data.Path.str.contains(thing)]
        
        elif sort_params[item] == 'DateCreated':
            minDate = eg.enterbox(msg='Enter the earliest date (yyyy-mm-dd):')
            maxDate = eg.enterbox(msg='Enter latest date (yyyy-mm-dd):',default=minDate)
            selected_data = selected_data[minDate:maxDate]
            print 'minDate: ',minDate,' maxDate: ',maxDate

        elif sort_params[item] == 'Zmass':
            foo = sort(selected_data.drop_duplicates('Zmass').Zmass.values)
            min = eg.indexbox(msg = 'Minimum mass value, in units of Mjup',choices = foo/mJup)
            max = min + eg.indexbox(msg = 'Maximum mass value, in units of Mjup',choices = foo[min::]/mJup)
            min = foo[min]
            max = foo[max]
            selected_data = selected_data[(selected_data[sort_params[item]] <= max) & (selected_data[sort_params[item]] >= min)] 
            print min," <= Zmass <= ",max

        elif sort_params[item] == 'SEMI':
            foo = sort(selected_data.drop_duplicates('SEMI').SEMI.values)
            min = eg.indexbox(msg='Minimum semi-major axis value, in AU',choices=foo/AU)
            max = eg.indexbox(msg='Maximum semi-major axis value, in AU',choices=foo/AU)
            min = foo[min]
            max = foo[max]
            selected_data = selected_data[(selected_data[sort_params[item]] <= max) & (selected_data[sort_params[item]] >= min)] 
            print min," <= SEMI <= ",max

        elif sort_params[item] == 'MinNumResults':
            min = eg.enterbox(msg='Enter the minimum number of results each simulation should have',default=1)
            values = map(check_num_sims_in_dir,selected_data.Path.values, [1*int(min)]*int(len(selected_data)) )
            selected_data = selected_data[values]
            print 'Minimum number of converged models required: ',min
        else:
            foo = sort(selected_data.drop_duplicates(sort_params[item])[sort_params[item]].values)
            min = eg.indexbox(msg='Minimum '+sort_params[item]+' value',choices=foo)
            max = eg.indexbox(msg='Maximum '+sort_params[item]+' value',choices=foo)
            min = foo[min]
            max = foo[max]
            selected_data = selected_data[(selected_data[sort_params[item]] <= max) & (selected_data[sort_params[item]] >= min)] 
            print min," <= ",sort_params[item]," <= ",max



            result = []
            #            print '\n+++++++++++++++++++++++++++'

    if len(selected_data) < 1:
        print 'No simulations match your requirements'
        print '+++++++++++++++++++++++++++\n'
    else:
        print str(len(selected_data))+ ' simulations meet these requirements'
        #    print 'Enter \'print result\' to see them!'
        #    print '+++++++++++++++++++++++++++'
        #    result = selected_data.Path.values
        #    result = result.tolist()

    return selected_data


def check_num_sims_in_dir(filepath,min):
    if len(os.listdir(filepath)) <= (min+2):
        #        print filepath+' has no (or few) converged models, so it is not getting plotted.'
        return False
    if not os.path.exists(filepath+"/atmos_data.txt"):
        #        print item+"/atmos_data.txt does not exist.  Skipping to next file."
        return False
    return True
