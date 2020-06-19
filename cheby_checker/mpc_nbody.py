# -*- coding: utf-8 -*-
# mpc_nbody/mpc_nbody/mpc_nbody.py

'''
----------------------------------------------------------------------------
mpc_nbody's wrapper module that calls the parser to parse an orbit,
and the reboundx/examples/ephem_forces n-body integrator.

Apr 2020
Mike Alexandersen & Matthew Payne & Matthew Holman

This module provides functionalities to
(a)
(b)
(c)

The output is then either saved to a file or passed to the orbit_cheby routine.
----------------------------------------------------------------------------
'''


# Import third-party packages
# -----------------------------------------------------------------------------
import sys
import os
import numpy as np


# Import neighbouring packages
# -----------------------------------------------------------------------------
try:  # Import ephem_forces from whereever REBX_DIR is set to live
    sys.path.append(os.environ['REBX_DIR'])
    from examples.ephem_forces.ephem_forces import integration_function
except (KeyError, ModuleNotFoundError):
    from reboundx.examples.ephem_forces.ephem_forces import integration_function

sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import parse_input

# Default for caching stuff using lru_cache
# -----------------------------------------------------------------------------


# Constants and stuff
# -----------------------------------------------------------------------------
DATA_PATH = os.path.realpath(os.path.dirname(__file__))
DATA_DIR = os.path.join(os.path.dirname(DATA_PATH), 'dev_data')
au_km = 149597870.700  # This is now a definition


# Data classes/methods
# -----------------------------------------------------------------------------
class NbodySim():
    '''
    Class for containing all of the N-body related stuff.
    '''

    def __init__(self, input_file=None, filetype=None, save_parsed=False):
        #If input filename provided, process it:
        if isinstance(input_file, str) & isinstance(filetype, str):
            self.pparticle = parse_input.ParseElements(input_file, filetype,
                                                       save_parsed=save_parsed)
        else:
            print("Keywords 'input_file' and/or 'filetype' missing; "
                  "initiating empty object.")
            self.pparticle = None
        self.geocentric = False  # Can be changed to something like
        #self.pparticle.geocentric if ParseElements gains knowledge.
        self.input_vectors = None
        self.input_n_particles = None
        self.output_times = None
        self.output_vectors = None
        self.output_n_times = None
        self.output_n_particles = None
        self.time_parameters = None

    def __call__(self, tstart=None, vectors=None, tstep=20, trange=600,
                 save_output=None, verbose=False):
        if vectors is None:
            vectors = self.pparticle
            if vectors is None:
                raise TypeError("If you didn't parse a particle from an input "
                                "file, you must supply 'vectors'.")
        if tstart is None:
            try:
                tstart = self.pparticle.time.tdb.jd
            except AttributeError:
                print("If you didn't parse a particle from an input file, "
                      "you must supply a 'tstart' value.")
                raise TypeError("If you didn't parse a particle from input "
                                "file, you must supply a 'tstart' value.")
        (self.input_vectors, self.input_n_particles, self.output_times,
         self.output_vectors, self.output_n_times, self.output_n_particles
         ) = run_nbody(vectors, tstart, tstep, trange, self.geocentric, verbose)
        print(f'###!!!{type(self.output_times):}!!!###' if verbose else '')
        self.time_parameters = [tstart, tstep, trange]
        if save_output is not None:
            if isinstance(save_output, str):
                self.save_output(output_file=save_output)
            else:
                self.save_output()

    def save_output(self, output_file='simulation_states.dat'):
        """
        Save all the outputs to file.

        Inputs:
        -------
        output_file : string, filename to write elements to.

        The file is overwritten if it already exists.
        """
        outfile = open(output_file, 'w')
        outfile.write(f'#Input vectors: [')
        for coo in self.input_vectors:
            outfile.write(f'{coo} ')
        outfile.write(f'\b]')
        outfile.write(f'\n#Input N_particles: {self.input_n_particles:}')
        outfile.write('\n#Start time, timestep, time range: '
                      f'{self.time_parameters:}')
        outfile.write(f'\n#Output N_times: {self.output_n_times:}')
        outfile.write(f'\n#Output N_particles: {self.output_n_particles:}')
        outfile.write('\n#')
        outfile.write('\n#Time               ')
        for j in np.arange(self.output_n_particles):
            outfile.write('x                  y                  '
                          'z                   dx                  '
                          '  dy                    dz                  ')
        for i, timei in enumerate(self.output_times):
            outfile.write(f'\n{timei:} ')
            for j in np.arange(self.output_n_particles):
                for k in np.arange(6):
                    outfile.write(f'{self.output_vectors[i, j, k]:} ')
        outfile.write('\n#End')

# Functions
# -----------------------------------------------------------------------------


def run_nbody(input_vectors, tstart, tstep, trange, geocentric=False,
              verbose=False):
    '''
    Run the nbody integrator with the parsed input.

    Input:
    ------
    input_vectors = Either ParseElements object,
                    list of ParseElements objects,
                    or numpy array of elements.
    tstart = float, Julian Date at start of integration.
    tstep = float or integer, major time step of integrator.
    trange = float or integer, rough total time of integration.
    geocentric = boolean, use geo- (True) or heliocentric (False)

    Output:
    -------
    reparsed_input = numpy array, input elements, reparsed into array
    n_particles = integer, the input number of particles
    times = numpy array, all the output times (including sub-steps)
    output_vectors = numpy array, output elements of dimensions
                                  (n_times, n_particles_out, 6)
    n_times = integer, number of time outputs
    n_particles_out = integer, number of output particles (different why?)
    '''
    # First get input (3 types allowed) into a useful format:
    reparsed_input, n_particles = _fix_input(input_vectors, verbose)
    # Now run the nbody integrator:
    (times, output_vectors, n_times, n_particles_out
     ) = integration_function(tstart, tstep, trange, geocentric,
                              n_particles, reparsed_input)
    return(reparsed_input, n_particles, times, output_vectors,
           n_times, n_particles_out)


def _fix_input(pinput, verbose=False):
    '''
    Convert the input to a useful format.

    Input:
    ------
    pinput = Either ParseElements object,
             list of ParseElements objects,
             or numpy array of elements.

    Output:
    -------
    reparsed = numpy array of elements.
    len(reparsed)//6 = integer, number of particles.
    '''
    if isinstance(pinput, parse_input.ParseElements):
        print('###!!!ONE!!!###' if verbose else '')
        els = pinput.barycentric_equatorial_cartesian_elements
        reparsed = np.array([els[i] for i in ['x_BaryEqu', 'y_BaryEqu',
                                              'z_BaryEqu', 'dx_BaryEqu',
                                              'dy_BaryEqu', 'dz_BaryEqu']])
    elif isinstance(pinput, list):
        print('###!!!TWO!!!###' if verbose else '')
        if isinstance(pinput[0], parse_input.ParseElements):
            print('###!!!TWO.5!!!###' if verbose else '')
            reparsed = []
            for particle in pinput:
                els = particle.barycentric_equatorial_cartesian_elements
                _ = [reparsed.append(els[i])
                     for i in ['x_BaryEqu', 'y_BaryEqu', 'z_BaryEqu',
                               'dx_BaryEqu', 'dy_BaryEqu', 'dz_BaryEqu']]
            reparsed = np.array(reparsed)
        else:
            reparsed = np.array(pinput)
    elif isinstance(pinput, np.ndarray):
        if (len(np.shape(pinput)) == 1) & (len(pinput) % 6 == 0):
            print('###!!!THREE!!!###' if verbose else '')
            reparsed = pinput
    else:
        raise(TypeError('"pinput" not understood.\n'
                        'Must be ParseElements object, '
                        'list of ParseElements or numpy array.'))
    return reparsed, len(reparsed) // 6

# End
