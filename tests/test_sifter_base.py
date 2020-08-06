# -*- coding: utf-8 -*-
# cheby_checker/tests/test_sifter_base

'''
    --------------------------------------------------------------
    tests of sifter's base class

    Jan 2020
    Matt Payne & Mike Alexandersen

    --------------------------------------------------------------
'''


# Import third-party packages
# --------------------------------------------------------------
import sys
import os
import numpy as np

# Import neighboring packages
# --------------------------------------------------------------
sys.path.append(os.path.dirname(os.path.dirname(
                os.path.realpath(__file__))))
from cheby_checker import sifter_precalc as precalc


def test_instantiation():
    '''Test that the empty class instantiates correctly.'''
    assert isinstance(precalc.Base(), precalc.Base)


def test_variables():
    '''
    Test that necessary variables are available
    and have the expected values/types...
    '''
    B = precalc.Base()
    
    assert 'sifter' in B._fetch_data_directory()
    assert 'HP_nside' in B.__dict__
    assert 'HP_order' in B.__dict__
    assert 'HP_npix' in B.__dict__
    assert isinstance(B.HP_npix, np.int64)


# End of file.
