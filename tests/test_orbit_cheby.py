# -*- coding: utf-8 -*-
# sifter/tests/test_base

'''
    --------------------------------------------------------------
    tests of orbit_cheby's base class
    
    Jan 2020
    Matt Payne
    
    --------------------------------------------------------------
'''


# Import third-party packages
# --------------------------------------------------------------
import sys
import os
import numpy as np
import pytest

# Import neighboring packages
# --------------------------------------------------------------
# Can't get the damn thing to work from symbolic link in ...
# '/Users/matthewjohnpayne/opt/anaconda3/lib/python3.7/site-packages/orbit_cheby/__init__.py'
# So just jacking this in for now ...
sys.path.append( os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'orbit_cheby') )
#os.chdir('/Users/matthewjohnpayne/Envs/orbit_cheby/orbit_cheby')

import orbit_cheby
import nbody_reader




# Convenience data / functions to aid testing
# --------------------------------------------------------------
# Set up a filepath (file will be created during testing)
test_filepath        = os.path.join(os.path.dirname(os.getcwd() ), 'dev_data', '2022AA_demo.txt')
name, times, states  = nbody_reader.parse_nbody_txt( test_filepath )

def create_single_MSC_from_array() :
    return orbit_cheby.MSC_Loader(FROM_ARRAY = True ,
                                 unpacked_provisional_designations = name,
                                 times_TDB = times,
                                 statearray = states).MSCs[0]

def create_list_of_MSCs_from_arrays(n=5) :
    arrays, names = [],[]
    for i in range(1,n):
        name, times, states  = nbody_reader.parse_nbody_txt( text_filepath )
        names.append(name+"_"+str(i))
        arrays.append(states)
    states_3D = np.stack(arrays, axis=2)

    return orbit_cheby.MSC_Loader(FROM_ARRAY = True ,
                                     unpacked_provisional_designations = names,
                                     times_TDB = times,
                                     statearray = states_3D).MSCs


# Actual tests ...
# --------------------------------------------------------------

'''
@pytest.mark.parametrize(('test_filepath'), [test_filepath])
def test_text_file_creation(test_filepath):
    
    # Remove test file if it exists
    if os.path.isfile(test_filepath):
        os.remove(test_filepath)
    
    # Use convenience func in nbody_reader to create a text file
    nbody_reader.create_nbody_txt(test_filepath)
'''






def test_create_empty_MSC():
    
    # Initialize the multi_sector_cheby
    result = orbit_cheby.MSC()

    # Check the result is as expected
    assert isinstance( MSC , orbit_cheby.MSC )
    assert 'sector_coeffs' in MSC.__dict__
    assert isinstance( MSC.sector_coeffs , list )


