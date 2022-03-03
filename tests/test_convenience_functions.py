# -*- coding: utf-8 -*-
# sifter/tests/test_orbit_cheby

"""
--------------------------------------------------------------
tests of orbit_cheby's base class

Jan 2020
Matt Payne

--------------------------------------------------------------
"""

# Import third-party packages
# --------------------------------------------------------------
import os
import pytest

# Import neighboring packages
# --------------------------------------------------------------
from cheby_checker import nbody


# Constants & Test Data
# -----------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'dev_data')
FLAT_FILES = [os.path.join(DATA_DIR, '2022AA_demo.txt') , os.path.join(DATA_DIR, 'simulation_states.dat')]
orbfit_filenames = [os.path.join(DATA_DIR, file) for file in ['30101.eq0_horizons', '30102.eq0_horizons']]


# Convenience data / functions to aid testing
# --------------------------------------------------------------
@pytest.mark.parametrize(('orbfit_file') , orbfit_filenames)
def test_Sim(orbfit_file):
    """
    A convenience function to return a simulation object
    Proper testing of mpc_nbody is done elsewhere (test_run_nbody.py)
    """
    # First, let's initiate the class with an input file:
    Sim = nbody.NbodySim() # 'eq' was an arg to this, may be deprecated per current nbody.py
    Sim._parse_orbfit_json(orbfit_file)

    # Now run the integrator, by calling the object.
    Sim(tstep=20, trange=1000)
    
    # Quick tests ...
    for attrib in ["output_times", "output_vectors"]:
        assert attrib in Sim.__dict__

    return Sim
