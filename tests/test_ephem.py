"""
	Tests of the ephem module
	Currently (20200804) incomplete
"""

# Import third-party packages
# --------------------------------------------------------------
import os

# Import neighboring packages
# --------------------------------------------------------------
import pytest

from cheby_checker import ephem, nbody

# Files / Directories
# --------------------------------------------------------------
HEAD_DIR = os.path.dirname(os.path.realpath(os.getcwd()))
DATA_DIR = os.path.join(HEAD_DIR, 'dev_data')
filenames = [os.path.join(DATA_DIR, file)
              for file in ['30101.eq0_horizons', '30102.eq0_horizons']]


# Tests ...
# --------------------------------------------------------------
def convenience_data_generation():
    """ Get data into database : See Demonstrate_EndToEnd_Orbit_Precalc.ipynb """

    # Use the Sim-object approach to run a simulation
    # TODO: Same fix as in test_convenience_functions.py
    # Sim = mpc_nbody.NbodySim(filenames[0], 'eq')
    Sim = nbody.NbodySim()
    Sim._parse_orbfit_json(filenames[0])
    Sim(tstep=20, trange=1000)

    # Initialize an MSC object
    M = orbit_cheby.MSC_Loader(NbodySim = Sim).MSCs[0]

    # Declare a "PreCalc" object
    P = precalc.PreCalc()
    # Do the upsert
    P.upsert(MSCs , observatoryXYZ)
    

@pytest.mark.skip(reason="designations, times not defined")
def test_instantiation():
    E = ephem.Ephem(designations, times)
    assert isinstance(E, ephem.Ephem), f'E is not of the expected type:{type(E)}'
