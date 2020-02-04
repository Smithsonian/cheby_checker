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
import sys, os
import numpy as np

# Import neighboring packages
# --------------------------------------------------------------
sys.path.append( os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'orbit_cheby') )
import precalc



def test_instantiation():
    pass # assert isinstance( precalc.Base() , precalc.Base)




# Call the tests while developing ...
#test_instantiation()
print('All tests passed')
