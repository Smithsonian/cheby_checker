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
import pytest
import numpy as np

# Import neighboring packages
# --------------------------------------------------------------
from sifter import precalc




def test_instantiation():
    assert isinstance( precalc.Base() , precalc.Base)

def test_variables():
    B = precalc.Base()
    
    # Test necessary variables are available and have the expected values ...
    assert 'sifter' in B._fetch_data_directory()
    assert 'HP_nside' in B.__dict__
    assert 'HP_order' in B.__dict__
    assert 'HP_npix' in B.__dict__
    assert isinstance(B.HP_npix, np.int64)



# End of file.
