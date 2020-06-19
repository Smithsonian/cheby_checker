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
from orbit_cheby import nbody_reader



# Convenience data / functions to aid testing
# --------------------------------------------------------------
# Set up a filepath (file will be created during testing)
test_filepath = os.path.join(os.path.dirname(os.getcwd() ), 'dev_data', '2022AA_demo.txt')


# Actual tests ...
# --------------------------------------------------------------

@pytest.mark.parametrize(('test_filepath'), [test_filepath])
def test_text_file_creation(test_filepath):
    
    # Remove test file if it exists
    if os.path.isfile(test_filepath):
        os.remove(test_filepath)

    # Use convenience func in nbody_reader to create a text file
    nbody_reader.create_nbody_txt(test_filepath)

    # Check that the test file has been created
    assert os.path.isfile(text_filepath)



@pytest.mark.parametrize(('test_filepath'), [test_filepath])
def test_text_file_creation(test_filepath):
    
    # Use convenience func in nbody_reader to create a text file
    nbody_reader.create_nbody_txt(test_filepath)
    
    # Parse the text file
    result = nbody_reader.parse_nbody_txt(test_filepath)

    # Check that the parsed result is as expected
    assert len(result) == 2
    name, a = result
    assert isinstance(name, str), isinstance(a, np.ndarray)
    assert a.shape == (20000,28)





