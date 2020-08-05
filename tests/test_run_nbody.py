# -*- coding: utf-8 -*-
# mpc_nbody/tests/test_run_nbody.py

'''
----------------------------------------------------------------------------
tests for mpc_nbody's run_nbody module.

Mar 2020
Mike Alexandersen & Matthew Payne & Matthew Holman

----------------------------------------------------------------------------
'''

# import third-party packages
# -----------------------------------------------------------------------------
import sys
import os
import numpy as np
import pytest
from astroquery.jplhorizons import Horizons

# Import neighbouring packages
# -----------------------------------------------------------------------------
try:  # Import ephem_forces from whereever REBX_DIR is set to live
    sys.path.append(os.environ['REBX_DIR'])
    from examples.ephem_forces import ephem_forces
except (KeyError, ModuleNotFoundError):
    from reboundx.examples.ephem_forces import ephem_forces

sys.path.append(os.path.dirname(os.path.dirname(
                os.path.realpath(__file__))))
from tests.test_parse_input import is_parsed_good_enough, compare_xyzv
from cheby_checker import mpc_nbody
from cheby_checker.parse_input import ParseElements


# Constants & Test Data
# -----------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__))), 'dev_data')


# Tests
# -----------------------------------------------------------------------------

def test_initialize_integration_function():
    '''
    If we put ANYTHING into the ephem_forces.integration_function,
    will it work or crash and burn?
    Most likely if there is a problem, it'll cause pytest to crash entirely,
    so might as well start with this.

    *** MJP ***
    This test is insufficient. It does not capture the problem related to the
    hard-coded ephem file, which returns a message
    "could not load DE430 file, fool!"

    *** MA ***
    Indeed, when that file is missing, the C code prints that insult,
    then causes a complete crash, exiting python without a python error.
    As such, it is impossible to guard against it, even with a
    try:
      stuff
    except:
      print('Oh no, you probably have files missing!')
    as it exits python entirely (somehow).
    When run interactively, the C prints the insult, but when run inside the
    pytest, it doesn't get printed (because pytest suppresses outputs) and
    python just crashes and leaves you with no test output.
    So in principle it's clear that it didn't work, but it doesn't actually
    give useful information.
    For now, I have found out that the two files needed are actually smaller
    than the GitHub 100 MB file limit, so I'll just commit them into the repo.
    '''
    tstart, tstep, trange = 2456184.7, 20.0, 600
    geocentric = 0
    n_particles = 1
    instates = np.array([-3.1, 2.7, 3.6, -0.006, -0.004, -0.002])
    (times, states, n_out, n_particles_out
     ) = ephem_forces.integration_function(tstart, tstep, trange, geocentric,
                                           n_particles, instates)
    assert n_particles_out == n_particles
    assert isinstance(n_particles_out, int)
    assert isinstance(n_out, int)
    assert isinstance(states, np.ndarray)
    assert isinstance(times, np.ndarray)


# A @pytest.mark.parametrize basically defines a set of parameters that
# the test will loop through.
@pytest.mark.parametrize(
    ('vectors', 'tstart', 'tstep', 'trange'),
    [
     (np.array([-2.093834952466475E+00, 1.000913720009255E+00,
                4.197984954533551E-01, -4.226738336365523E-03,
                -9.129140909705199E-03, -3.627121453928710E-03]),
      2456117.641933589, 20, 600),
     (ParseElements(os.path.join(DATA_DIR, '30101.eq0_horizons'), 'eq'),
      2456117.641933589, 20, 600),
     ([ParseElements(os.path.join(DATA_DIR, '30101.eq0_horizons'), 'eq'),
       ParseElements(os.path.join(DATA_DIR, '30102.eq0_horizons'), 'eq')],
      2456117.641933589, 20, 600)
      ])
def test_run_nbody(vectors, tstart, tstep, trange):
    '''
    Test whether the run_nbody function works correctly.
    This test for now only tests that the function doesn't crash and burn,
    not the actual output.
    '''
    (input_vectors, input_n_particles,
     output_times, output_vectors, output_n_times, output_n_particles
     ) = mpc_nbody.run_nbody(vectors, tstart, tstep, trange, geocentric=False)
    assert isinstance(input_vectors, np.ndarray)
    assert isinstance(input_n_particles, int)
    assert isinstance(output_times, np.ndarray)
    assert isinstance(output_vectors, np.ndarray)
    assert isinstance(output_n_times, int)
    assert isinstance(output_n_particles, int)
    assert len(output_times) == output_n_times
    assert output_n_particles in np.shape(output_vectors)
    assert output_n_times in np.shape(output_vectors)
    assert (6 in np.shape(output_vectors)) | (27 in np.shape(output_vectors))


# Splitting the parameters into two @pytest.mark.parametrize statements
# essentially makes it a nested loop (so all combinations are tested).
@pytest.mark.parametrize(
    ('tstart', 'tstep', 'trange', 'geocentric', 'targets', 'id_type'),
    [
     (2456117.641933589, 20.0, 600, 0, ['30101'], ['smallbody']),
     (2456184.7528431923, 20.0, 600, 0, ['30102'], ['smallbody']),
     (2456142.5, 20.0, 60, 0, ['30101', '30102'],
      ['smallbody', 'smallbody']),
      ])
@pytest.mark.parametrize(
    ('threshold_xyz', 'threshold_v'),
    [
     (1e-10, 1e-11),  # 1e-10 au ~ 15m, 1e-11 au/day ~ 1.5 m/day
     (5e-11, 2e-13),  # 5e-11 au ~ 7.5m, 2e-13 au/day ~ 30 mm/day
      ])
def test_nbody_vs_Horizons(tstart, tstep, trange, geocentric,
                           targets, id_type, threshold_xyz, threshold_v):
    '''
    Test that putting input from Horizons in gives Horizons consistent output.
    '''
    centre = '500' if geocentric else '500@0'
    # Make the single array with 6 elements for each particle.
    horizons_in = []
    for i, targi in enumerate(targets):
        horizons_in = np.concatenate([horizons_in, nice_Horizons(targi, centre,
                                      tstart, id_type[i])])
    # Run nbody integrator
    (input_vectors, input_n_particles,
     output_times, output_vectors, output_n_times, output_n_particles
     ) = mpc_nbody.run_nbody(horizons_in, tstart, tstep, trange, geocentric)
    # Check 20 time steps (or less if there are many)
    for j in set(np.linspace(0, output_n_times - 1, 20).astype(int)):
        # Get Horizons positions for that time and compare
        for i, targi in enumerate(targets):
            horizons_xyzv = nice_Horizons(targi, centre, output_times[j],
                                          id_type[i])
            mpc_xyzv = output_vectors[j, i, :]
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
    print(input_vectors, horizons_in, np.all(input_vectors == horizons_in))
    assert np.all(input_vectors == horizons_in)
    print(input_n_particles, output_n_particles,
          input_n_particles == output_n_particles)
    assert input_n_particles == output_n_particles
    ### This should get refactored to use is_nbody_output_good_enough !!!


def test_NbodySim_empty():
    '''
    Test the mpc_nbody.NbodySim class. Test empty initialization.
    '''
    assert isinstance(mpc_nbody.NbodySim(), mpc_nbody.NbodySim)


@pytest.mark.parametrize(
    ('data_file', 'file_type', 'holman_ic_test_file', 'nbody_test_file'),
    [
     pytest.param('30101.ele220', 'ele220', 'holman_ic_30101', 'nbody_30101',
                  marks=pytest.mark.xfail(reason='Not implemented yet.')),
     pytest.param('30102.ele220', 'ele220', 'holman_ic_30102', 'nbody_30102',
                  marks=pytest.mark.xfail(reason='Not implemented yet.')),
     ('30101.eq0_horizons', 'eq', 'holman_ic_30101_horizons',
      'nbody_30101_horizons'),
     ('30102.eq0_horizons', 'eq', 'holman_ic_30102_horizons',
      'nbody_30102_horizons'),
      ])
def test_NbodySim(data_file, file_type, holman_ic_test_file, nbody_test_file):
    '''
    Test the mpc_nbody.NbodySim class. Test empty initialization.
    '''
    Sim = mpc_nbody.NbodySim(os.path.join(DATA_DIR, data_file), file_type,
                             save_parsed=True)
    Sim(tstep=20, trange=600, save_output=True)
    is_parsed_good_enough(os.path.join(DATA_DIR, holman_ic_test_file))
    is_nbody_output_good_enough(Sim.output_times, Sim.output_vectors,
                                target=data_file[:5])


# Non-test helper functions
# -----------------------------------------------------------------------------

def is_nbody_output_good_enough(times, data, target='30102'):
    '''
    Helper function for determining whether the saved output from an nbody
    integration is good enough.
    '''
    # Check 20 timesteps (or less if there are many)
    some_times = np.linspace(0, len(times) - 1, 20).astype(int)
    for j in set(some_times):
        # Get Horizons positions for that time and compare
        horizons_xyzv = nice_Horizons(target, '500@0', times[j],
                                      'smallbody')
        mpc_xyzv = data[j, 0, :]
        # Check whether position/v within threshold.
        error, good_tf = compare_xyzv(horizons_xyzv, mpc_xyzv,
                                      5e-11, 2e-13)  # 7.5m, 30 mm/day
        if np.all(good_tf):
            print('Awesome!')
        else:
            print(f'Time, timestep: {times[j]:}, {j:}')
            print(f'Horizons : {horizons_xyzv:}')
            print(f'N-body   : {mpc_xyzv:}')
            print(f'Position off by [au]: {error[:3]:}')
            print(f'Velocity off by [au/day]: {error[3:6]:}')
        assert np.all(good_tf)


def nice_Horizons(target, centre, epochs, id_type):
    '''
    Only require the inputs I actually want to vary.
    Return in the format I actually want, not an astropy table.
    '''
    horizons_table = Horizons(target, centre, epochs=epochs, id_type=id_type)
    horizons_vector = horizons_table.vectors(refplane='earth')
    horizons_xyzv = horizons_vector['x', 'y', 'z', 'vx', 'vy', 'vz']
    return np.array(list(horizons_xyzv.as_array()[0]))


# End
