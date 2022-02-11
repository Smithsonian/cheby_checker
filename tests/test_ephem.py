"""
	Tests of the ephem module
	Currently (20200804) incomplete
"""

# Import third-party packages
# --------------------------------------------------------------
import numpy as np
import sys, os
from astropy_healpix import healpy

# Import neighboring packages
# --------------------------------------------------------------
from cheby_checker import ephem

# Files / Directories
# --------------------------------------------------------------
HEAD_DIR = os.path.dirname(os.path.realpath(os.getcwd()))
DATA_DIR = os.path.join(HEAD_DIR, 'dev_data')
filenames = [os.path.join(DATA_DIR, file)
              for file in ['30101.eq0_horizons', '30102.eq0_horizons']]


# Tests ...
# --------------------------------------------------------------
def convenience_data_generation():
    ''' Get data into database : See Demonstrate_EndToEnd_Orbit_Precalc.ipynb'''

    # Use the Sim-object approach to run a simulation
    Sim = mpc_nbody.NbodySim(filenames[0], 'eq')
    Sim(tstep=20, trange=1000)

    # Initialize an MSC object
    M = orbit_cheby.MSC_Loader(NbodySim = Sim).MSCs[0]

    # Declare a "PreCalc" object
    P = precalc.PreCalc()

    # Do the upsert
    P.upsert( MSCs , observatoryXYZ)
    


def test_instantiation():
    # Instantiate :
    # NB(1) As written, need to instantiate with variables ...
    # ... __init__(self, designations, times, observatoryXYZ = None , obsCode = None)
    # NB(2) As written, queried designation needs to be in database
    #E = ephem.Ephem()
    assert isinstance(E,ephem.Ephem) , f'E is not of the expected type:{type(E)}'


