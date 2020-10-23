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

sys.path.append(
                os.path.dirname(os.path.dirname(
                                                os.path.realpath(__file__))))
from cheby_checker import parse_input


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

    def __init__(self, input=None, filetype='eq', save_parsed=False):

        # Expected class variables
        self.pparticle          = None # This will be ParseElements object ...
        self.geocentric         = False
        self.input_vectors      = None
        self.input_n_particles  = None
        self.output_times       = None
        self.output_vectors     = None
        self.output_n_times     = None
        self.output_n_particles = None
        self.time_parameters    = None

        #If input filename provided, process it using ParseElements :
        try:
            assert input
            self.pparticle = parse_input.ParseElements( input       = input,
                                                        filetype    = filetype,
                                                        save_parsed = save_parsed)
        except Error as e:
            print(f"Instantiating empty object.")
            print(f"Reason:{e}")
            


    def __call__(self,
                tstart      ,
                tstep       ,
                trange      ,
                epoch       =None,
                vectors     =None,
                covariances =None,
                save_output =None,
                verbose     =False):
        '''
        Allows the NbodySim object to be called as a function
         - Performs integration by calling *integration_function* within *run_nbody*
        '''
        
        # Get the initial conditions for the particle integration
        if vectors is None or epoch is None :
            try:
                epoch       = self.pparticle.time.tdb.jd
                vectors     = self.pparticle.bary_eq_vec
                covariances = self.pparticle.bary_eq_cov
            except AttributeError as e :
                print(f "Invalid ParseElements ? : {e}")

        if vectors is None or epoch is None :
            raise TypeError("If you didn't parse a particle from an input "
                                "file, you must supply 'vectors' and 'epoch'.")

        # Call the routine that runs the nbody integration
        (self.epoch,
        self.input_vectors,
        self.input_covariances,
        self.input_n_particles,
        self.output_times,
        self.output_vectors,
        self.output_covariances,
        self.output_n_times,
        self.output_n_particles) = self.run_nbody(  epoch
                                                    vectors,
                                                    covariances
                                                    tstart,
                                                    tstep,
                                                    trange,
                                                    geocentric  = self.geocentric,
                                                    verbose     = verbose)
        print(f'###!!!{type(self.output_times):}!!!###' if verbose else '')
        
        # Save output if requested
        if save_output is not None:
            if isinstance(save_output, str):
                self.save_output(output_file=save_output)
            else:
                self.save_output()


    # --------------------------------------------
    # Functions to run integration & parse results
    # --------------------------------------------

    def run_nbody(  self,
                    epoch
                    input_states,
                    input_covariances,
                    tstart,
                    tstep,
                    trange,
                    init_covariances = None,
                    geocentric=False,
                    verbose=False):
        '''
        Run the nbody integrator with the parsed input.

        Input:
        ------
        epoch :          float
         - Common epoch at which all input vectors are defined
        input_states :   np.ndarray
         - 2D, shape = (N_particles, 6)
         - elements (xyzuvw) for N_particles to be integrated
        input_covariances :   np.ndarray or None
         - 3D, shape = (N_particles, 6,6)
        tstart = float,
         - Julian Date at start of integration.
        tstep = float or integer,
         - major time step of integrator.
        trange = float or integer,
         - rough total time of integration.
        geocentric = boolean,
         - use geo- (True) or heliocentric (False)

        Output:
        -------
        reparsed_input      = numpy array
         - input elements, reparsed into array
        n_particles         = integer,
         - the input number of particles
        times               = numpy array,
         - all the output times (including sub-steps)
        states              = numpy array,
         - output elements of dimensions (n_times, n_particles_out, 6)
        output_covariance   = numpy array,
         - output elements of dimensions (n_times, n_particles_out*6, 6)
        n_times             = integer,
         - number of time outputs
        n_particles_out     = integer,
         - number of output particles (different why?)
        '''
        
        # Get the number of particles
        assert input_states.ndim == 2 and input_states.shape[1] == 6
        n_particles = input_states.shape[0]
        
        # Now run the nbody integrator:
        times, output_vectors, n_times, n_particles_out = integration_function( tstart,
                                                                                tstep,
                                                                                trange,
                                                                                geocentric,
                                                                                n_particles,
                                                                                input_states)
                                                                                
        # Split the output vectors
        final_states, partial_derivatives = _split_integration_output(output_vectors , n_times, n_particles_out )
                                  
        # Calculate the covariance matrix (at each timestep) from the \partial X / \partial X_0 data
        if init_covariances is not None:
            final_covariance_arrays = self._get_covariance_from_tangent_vectors(init_covariances, partial_derivatives )
        else:
            final_covariance_arrays = None
            
        return( epoch,
                init_states,
                init_covariances,
                n_particles,
                times,
                final_states,
                final_covariance_arrays,
                n_times,
                n_particles_out)

    def _split_integration_output(self, output_from_integration_function , n_times, n_particles_out ):
        '''
        Get the states and partials out of the array returned from the integration_function
        
        MJP : 20200902 : Reminder to self.
        I may want to make the state array be 3D      : (n_times, n_particles, 6 )
        I may want to make the covariance array be 4D : (n_times, n_particles, 6,6 )
        '''
        original_shape = output_from_integration_function.shape
        assert original_shape == (n_times, 7*n_particles_out, 6)
        
        # Make a handy mask : True => CoVariance Components, False => State components
        mask = np.ones_like(output_from_integration_function, dtype=bool)
        mask[:,0:n_particles_out,:] = False

        # State (xyzuvw) components:
        states = output_from_integration_function[~mask].reshape(original_shape[0],-1,6)

        # Partial Deriv Components
        # MJP 20200902: *** UNCLEAR WHETHER THESE WILL NEED TO BE TRANSPOSED ***
        partials = output_from_integration_function[mask].reshape(original_shape[0],n_particles_out,6,6)

        return states, partials
        
    def _get_covariance_from_tangent_vectors(self, init_covariances, partial_derivatives ):
        '''
        
        Follow Milani et al 1999
        Gamma_t = [partial X / partial X_0] Gamma_0 [partial X / partial X_0]^T
        '''
        assert init_covariances.ndim == 3 and partial_derivatives.ndim == 4, \
            f'init_covariances.ndim={init_covariances.ndim} , partial_derivatives.ndim={partial_derivatives.ndim}'
        
        # Take the inverse of the covariance matrix to get the normal matrix
        # NB: We make a stack of identical matricees to use in the matrix multiplication below
        Gamma0      = np.linalg.inv(init_covariances)
        GammaStack0 = np.tile( Gamma0, (partial_derivatives.shape[0],1,1,1) )

        # We need each of the individual pd arrays to be individually transposed
        # NB, tuple fixes dimensions 0 & 1 , while indicates that dimensions 2 & 3 will be swapped/transposed
        pds_transposed  = partial_derivatives.transpose( (0,1,3,2) )

        # Do matrix multiplication: using the pd's to get the CoV as a func of time
        # NB matmul/@ automagically knows how to work on a stack of matricees
        # - https://stackoverflow.com/questions/34142485/difference-between-numpy-dot-and-python-3-5-matrix-multiplication
        GammaStack_t    = pds_transposed @ GammaStack0 @ partial_derivatives
        
        # Magically, np.linalg.inv also knows how to deal with a stack of arrays/matrices
        # - https://stackoverflow.com/questions/11972102/is-there-a-way-to-efficiently-invert-an-array-of-matrices-with-numpy
        CoV_t           = np.linalg.inv( GammaStack_t )
        
        return CoV_t
            

    # -----------------------------------
    # Misc Funcs ...
    # -----------------------------------

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
                      f'{tstart, tstep, trange}')
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


    """
    def _fix_input(pinput, init_covariances, verbose=False):
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
        if isinstance(pinput, parse_input.ParseElements) :
            print('###!!!ONE!!!###' if verbose else '')
            init_states      = pinput.bary_eq_vec
            # We wrap in an extra array because we want the covariance input to be 3D: shape = (n_particles,6,6)
            init_covariances = np.array( [pinput.bary_eq_cov] )

        elif isinstance(pinput, list) :
            print('###!!!TWO!!!###' if verbose else '')
            if np.all( [ isinstance(_, parse_input.ParseElements) for _ in pinput] ):
                print('###!!!TWO.5!!!###' if verbose else '')
                # We flatten here because the *integration_function* demands 1D input
                init_states         = np.array( [p.bary_eq_vec for p in pinput] ).flatten()
                init_covariances    = np.array( [p.bary_eq_cov for p in pinput] )
            else:
                init_states         = np.array(pinput)
                init_covariances    = np.array(init_covariances) if init_covariances is not None else None
                
        elif isinstance(pinput, np.ndarray):
            if (pinput.ndim == 1) & (len(pinput) % 6 == 0):
                print('###!!!THREE!!!###' if verbose else '')
                init_states         = pinput
                init_covariances    = np.array(init_covariances) if init_covariances is not None else None

        else:
            raise(TypeError('"pinput" not understood.\n'
                            'Must be ParseElements object, '
                            'list of ParseElements, or numpy array.'))
                            
        return init_states, init_covariances, len(init_states) // 6
    """



# End
