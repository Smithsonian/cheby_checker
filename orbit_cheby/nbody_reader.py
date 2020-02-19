# -*- coding: utf-8 -*-
# orbit_cheby/orbit_cheby/nbody_reader

'''
    --------------------------------------------------------------
    sifter's / mpchecker's orbit_cheby module.
    
    Jan 2020
    Matt Payne & Margaret Pan & Mike Alexandersen
    
    This module provides functionalities to read the output 
    from nbody simulations
    
    Obviously this presupposes that the result of nbody files are
    in some standard format
    
    One standard format currently implemented below is JSON
    
    --------------------------------------------------------------
    '''


# Import third-party packages
# --------------------------------------------------------------
import sys, os
import numpy as np
import json
import itertools


# Import neighboring packages
# --------------------------------------------------------------


# Define top-line parameters
# --------------------------------------------------------------

# Field names
object_name     = 'unpacked_designation'
time_fieldname  = 'MJD_TDB' # To be clarified if this is correct
coord_names     = ['x','y','z','vx','vy','vz']
# This generates a 'triangular' list of 21-combinations of coord-names
covar_names     = [ "_".join([coord_names[i], coord_names[j]]) for i in range(len(coord_names)) for j in range(i,len(coord_names))  ]
# Number of fields
nCovars     = len( covar_names )
nComponents = nCovars + len(coord_names)
nFields     = nComponents + 1

# Functions to read nbody results
# --------------------------------------------------------------

def parse_nbody_json( json_filepath ):
    '''
        Read a json file and return a dictionary
    '''
    # Simple read of json -> dictionary
    if os.path.isfile( json_filepath ):
        with open(json_filepath) as f:
            nbody_dict =  json.load(f)

    # Check the dictionary is valid
    nbody_dict = check_nbody_dict( nbody_dict )

    # Convert to numpy arrays
    lists = [coord_names, covar_names, [time_fieldname] ]
    for item in itertools.chain(*lists):
        nbody_dict[item] = np.asarray(nbody_dict[item])

    return nbody_dict

def parse_nbody_txt( text_filepath ):
    '''
        Read a text file
        
        returns:
        --------
        name: str
         - name of object being processed: hopefully is an UNPACKED DESIGNATION
        array:
         - time, state, triangular-cov
    '''
    if os.path.isfile( text_filepath ):
        
        # Simple read of file -> array
        array =  np.loadtxt(text_filepath)

        # Get name
        with open(text_filepath) as fh:
            # remove '\n' from the end, and '# ' from the start
            name = fh.readline().strip()[2:]

    
    return name, array



# Functions to check nbody validity
# --------------------------------------------------------------

def check_nbody_dict( nbody_dict ):
    '''
        Check whether the input json has the expected
        structure / variables

    '''
    # is this actually a dictionary ?
    assert isinstance( nbody_dict , dict )

    # are all the required keys present?
    lists = [coord_names, covar_names, [object_name] , [time_fieldname]]
    for item in itertools.chain(*lists):
        assert item in nbody_dict
    
    # are all the lists of the required length ?
    lists = [coord_names, covar_names]
    for item in itertools.chain(*lists):
        assert len(nbody_dict[item]) == len(nbody_dict[time_fieldname])

    return nbody_dict

# Functions to create json
# --------------------------------------------------------------

def create_nbody_json( json_filepath ):
    '''
        Convenience function to create a json file of coordinates
        Is only for the purposes of data exploration 
        
        *** NOTE THAT IN PRACTICE MJP DOES NOT LIKE THIS APPROACH ***
        *** IT SEEMS TO BE TOO TIME-INTENSIVE ***********************
        
    '''
    # Dictionary to hold things ...
    data_dict = {object_name : '2022 AA' }
    
    # Fake an array of 20k times
    data_dict[time_fieldname] = np.arange(40000, 60000, 1)
    
    # Fake some slowly varying coordinate data
    r = np.random.uniform(low=2.5, high=3.5)
    p_days = r**(3./2.) * 365.25
    data_dict['x'] = r*np.cos( 2.*np.pi * data_dict[time_fieldname]/p_days )
    data_dict['y'] = r*np.sin( 2.*np.pi * data_dict[time_fieldname]/p_days )
    
    # Doing this to check whether zeros cause problems for cheby
    data_dict['z'] = np.zeros(len(data_dict[time_fieldname]))
    
    # Doing this to see if the # coefficients is stable
    for n, v in enumerate(['vx','vy','vz']):
        data_dict[v] = 0.1*n*data_dict['x']
    for n, c in enumerate(covar_names):
        data_dict[c] = 0.01*n*data_dict['x']

    # Json hates arrays, so convert to list
    for k,v in data_dict.items():
        if isinstance(v, np.ndarray):
            data_dict[k] = v.tolist()

    # Save to file as json
    with open(json_filepath, 'w') as json_file:
        json.dump(data_dict, json_file, sort_keys=True)

def create_nbody_txt( text_filepath ):
    '''
        Convenience function to create a text-file of coordinates
        ****This is only for the purposes of data exploration****
        
        #  I am constructing an array that is ~20,000 rows long
        # Each row stores 1-time, 6-posn&vel, and 21-CoVarCoeffs
        #
        #    ------ 28 ------>>>
        #   |         [[ t, x, y, z, vx, vy, vz, x_x, x_y, ... ]
        #   |         ...
        # 20,000      ...
        #   |         ...
        #   |
        #   V         ...
        #   V         ]
        #

    '''
    
    # times
    times = np.arange(40000, 60000, 1)
    
    # Set up empty array and put times into zeroth element
    a = np.zeros( ( len(times) , nFields)  )
    a[:,0] = times
    
    # Fake some slowly varying coordinate data
    r = np.random.uniform(low=2.5, high=3.5)
    p_days = r**(3./2.) * 365.25
    a[:,1] = r*np.cos( 2.*np.pi * times/p_days )
    a[:,2] = r*np.sin( 2.*np.pi * times/p_days )
    
    # Leave z == 0 to check whether zeros cause problems for cheby
    #a[3]
    
    # Populate vel & CoVar components
    for n, v in enumerate(['vx','vy','vz']):
        a[:,4+n]  = 0.1*n*a[:,1]
    for n, c in enumerate(covar_names):
        a[:,7+n]  = 1.e-9*n*a[:,1]

    # Save to file
    np.savetxt(text_filepath , a , header='2022 AA' )

