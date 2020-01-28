# -*- coding: utf-8 -*-
# sifter/tests/test_base

'''
    --------------------------------------------------------------
    tests of sifter's base class
    
    Jan 2020
    Matt Payne & Mike Alexandersen
    
    --------------------------------------------------------------
'''


# Import third-party packages
# --------------------------------------------------------------
import sys, os
import numpy as np

# Import neighboring packages
# --------------------------------------------------------------
sys.path.append( os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'sifter') )
import precalc



def test_instantiation():
    assert isinstance( precalc.Base() , precalc.Base)

def test_variables():
    B = precalc.Base()
    
    # Test necessary variables are available and have the expected values ...
    assert 'sifter' in B._fetch_data_directory()
    assert 'HP_nside' in B.__dict__
    assert 'HP_order' in B.__dict__
    assert 'npix' in B.__dict__
    assert isinstance(B.npix, np.int64)




# Call the tests while developing ...
test_instantiation()
test_variables()
