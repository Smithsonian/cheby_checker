# -*- coding: utf-8 -*-
# mpc_nbody/tests/test_run_nbody.py

"""
----------------------------------------------------------------------------
tests for mpc_nbody's run_nbody module.

Mar 2020
Mike Alexandersen & Matthew Payne & Matthew Holman

----------------------------------------------------------------------------
"""

# import third-party packages
# -----------------------------------------------------------------------------
import sys
import os
import numpy as np
import pytest
from astroquery.jplhorizons import Horizons

# Import neighbouring packages
# -----------------------------------------------------------------------------
sys.path.append(os.environ['REBX_DIR'])
from examples.ephem_forces import ephem_forces
from tests.test_nbody_NbodySim import compare_xyzv

# Constants & Test Data
# -----------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__))), 'dev_data')

@pytest.mark.skip(reason="archaic")
def test_nbody_vs_Horizons(tstart, tstep, trange, geocentric,
                           targets, id_type, threshold_xyz, threshold_v):
    """
    Test that putting input from Horizons in gives Horizons consistent output.
    """
    centre = '500' if geocentric else '500@0'
    
    # Make the single array with 6 elements for each particle.
    horizons_in = []
    for i, targi in enumerate(targets):
        horizons_in = np.concatenate([horizons_in,
                                        nice_Horizons(targi, centre, tstart, id_type[i])])
                                        
    print('Pre-Integration...')
    for _ in horizons_in: print('%18.16f'%_)
    print( tstart, tstep, trange, geocentric)

    # Run nbody integrator
    (   input_vectors,
        input_covariances,
        input_n_particles,
        output_times,
        output_vectors,
        output_covariance,
        output_n_times,
        output_n_particles
     ) = mpc_nbody.run_nbody(horizons_in, tstart, tstep, trange, init_covariances = None, geocentric=geocentric )
    print('Finished integration')
    assert False
    # Check ~5 time steps (or less if there are many)
    for j in sorted(set(np.linspace(0, output_n_times - 1, 5).astype(int)) ):
        
        # Get Horizons positions for that time and compare
        for i, targi in enumerate(targets):
            
            # Get the states we want to compare
            print(j, ' : output_times[j]',output_times[j])
            horizons_xyzv = nice_Horizons(targi, centre, output_times[j],
                                          id_type[i])
            mpc_xyzv      = output_vectors[j, i, :]

            print('output_vectors',output_vectors.shape )
            # Check whether position/v within threshold.
            error, good_tf = compare_xyzv(horizons_xyzv, mpc_xyzv,
                                          threshold_xyz, threshold_v)
            if np.all(good_tf):
                print('Awesome!')
            else:
                print(f'Time, timestep: {output_times[j]:}, {j:}')
                print(f'Horizons : {horizons_xyzv:}')
                print(f'N-body   : {mpc_xyzv:}')
                print(f'Position off by [au]: {error[:3]:}')
                print(f'Velocity off by [au/day]: {error[3:6]:}')
            assert np.all(good_tf)
    assert output_n_particles == len(targets)
    assert np.all(input_vectors == horizons_in)
    assert input_n_particles == output_n_particles
    ### This should get refactored to use is_nbody_output_good_enough !!!



# Non-test helper functions
# -----------------------------------------------------------------------------
def nice_Horizons(target, centre, epochs, id_type):
    """
    Only require the inputs I actually want to vary.
    Return in the format I actually want, not an astropy table.
    """
    horizons_table  = Horizons(target, centre, epochs=epochs, id_type=id_type)
    horizons_vector = horizons_table.vectors(refplane='earth')
    horizons_xyzv   = horizons_vector['x', 'y', 'z', 'vx', 'vy', 'vz']
    print('nice_Horizons : ', np.array(list(horizons_xyzv.as_array()[0])))
    return np.array(list(horizons_xyzv.as_array()[0]))

#test_nbody_vs_Horizons(2458850.0, 1.0, 10., 0, ['2020 CD3'], ['smallbody'], 1e-10, 1e-11)

try:  # Import ephem_forces from wherever REBX_DIR is set to live
    sys.path.append(os.environ['REBX_DIR'])
    from examples.ephem_forces.ephem_forces import integration_function
except (KeyError, ModuleNotFoundError):
    from reboundx.examples.ephem_forces.ephem_forces import integration_function

tstart, tstep, trange, geocentric, n_particles, init_states = \
    2458850.0, 0.1, 10.0, 0, 1, \
    np.array([-0.18297081,  0.8869136,   0.39472742, -0.01719277, -0.0028644,  -0.00117345])
# TODO: Not sure if these last three (required) params are correct. See also test_malloc_reboundx.py
n_var = 6
invar_part = np.zeros(6, dtype=int)
invar = np.identity(6)
print('tstart,    tstep,    trange,    geocentric,    n_particles,    init_states')
print(tstart,    tstep,    trange,    geocentric,    n_particles,    init_states)
times, states, var, var_ng, status = integration_function(
      tstart, tstart + trange, tstep,
      geocentric, n_particles, init_states,
      n_var, invar_part, invar
    )
