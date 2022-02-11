# -*- coding: utf-8 -*-

'''
    --------------------------------------------------------------
    tests of orbit_cheby's base class
    ../tests/test_orbit_cheby_horizons.py
    
    Concentrating on tests that use the exact same starting positions
    as horizons so that I can double-check the accuracy of ..
    ... the rebound-to-cheby conversions
    ... the XYZ-to-RADEC conversions

    Feb 2022
    Matt Payne
    
    --------------------------------------------------------------
'''


# Import third-party packages
# --------------------------------------------------------------
import sys
import os
import numpy as np
import pytest
import json
import glob
import math
from datetime import datetime
from astropy.time import Time

# Import MPC packages
# -----------------------------------------------------------------------------
from mpc_orb.parse import MPCORB

# Import neighboring packages
# --------------------------------------------------------------
sys.path.append(os.environ['REBX_DIR'])
from examples.ephem_forces import ephem_forces

from checby_checker import nbody
from checby_checker import orbit_cheby
from checby_checker import cheby_checker
from checby_checker import obs_pos
from checby_checker import coco
from checby_checker import convenience_Horizons as Horizons

this_dir = os.path.abspath(os.path.dirname( __file__ ))
repo_dir = os.path.abspath(os.path.dirname( this_dir ))
data_dir = os.path.join(repo_dir, 'dev_data')
json_dir = os.path.join(data_dir, 'mpc_orb_jsons')
std_json_dir = os.path.join(json_dir, 'standard_mp')


# Constants & Test Data
# -----------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'dev_data')

FLAT_FILES = [  os.path.join(DATA_DIR, '2022AA_demo.txt') ,
                os.path.join(DATA_DIR, 'simulation_states.dat')]
orbfit_filenames = [os.path.join(DATA_DIR, file) for file in ['30101.eq0_horizons', '30102.eq0_horizons']]
mpc_orb_json_files = glob.glob(std_json_dir + '/*.json' )


# Convenience data / functions to aid testing
# --------------------------------------------------------------


def convenience_call_to_nbody_run_mpcorb( json_filepath_or_list, tstart=2459200 , tstop=2459295 ):
    ''' 
        A convenience function to return a simulation object
        Proper testing of mpc_nbody is done elsewhere (test_nbody_run.py)
    '''

    # Instantiate
    N = nbody.NbodySim()
    #N.verbose = True
    
    # Now run the integrator
    # NB: tstart = 2459200 = First date in Sector # 600 (check using B.map_JD_to_sector_number_and_sector_start_JD(2459200, B.standard_MJDmin )
    #     tstop = 2459295 = :Last day in Sector # 602   (check using B.map_JD_to_sector_number_and_sector_start_JD(2459295, B.standard_MJDmin )
    # This is useful to know for later testing ...
    mpcorb_list = [ json_filepath_or_list ] if isinstance(json_filepath_or_list,str) else json_filepath_or_list
    N.run_mpcorb( tstart = tstart , tstop = tstop , mpcorb_list = mpcorb_list )

    # Quick tests ...
    for attrib in ["output_times", "output_states"]:
        assert attrib in N.__dict__
        assert isinstance( N.__dict__[attrib] , np.ndarray )

    return N
    
def convenience_call_to_get_MSCs_from_file_via_nbody_run( json_filepath_or_list, tstart=2459200 , tstop=2459295):
    '''
        A convenience function to return a list of MSCs
        Starts from mpc_orb_json file(s)
        Runs NBody integration via nbody.NbodySim.run_mpcorb
        Uses MSC_Loader to do all of the work to declare and populate a list of MSC objects
        Tests of MSC_Loader are performed below (e.g. *test_create_loader* & *test_loader_from_nbodysim*)
        
    '''

    # Instantiate Nbody Object, Run NBody Simulation, and return populated Object
    N = convenience_call_to_nbody_run_mpcorb(json_filepath_or_list , tstart=tstart , tstop=tstop)
    
    # Use the MSC_Loader to do all of the work to declare and populate a list of MSC objects
    MSCs = orbit_cheby.MSC_Loader(NbodySim = N).MSCs
    
    return MSCs
    


def similar_xyzuvw(xyzv0, xyzv1, threshold_xyz=1e-13, threshold_v=1e-14): # 15 mm, 1.5 mm/day
    '''
    Calculate the difference between two sets of cartesian coordinates.
    '''
    if isinstance(xyzv0, list):
        xyzv0 = np.array(xyzv0)
    if isinstance(xyzv1, list):
        xyzv1 = np.array(xyzv1)
    error = xyzv0 - xyzv1
    if len(error) == 3:
        good_tf = np.abs(error) < np.array([threshold_xyz] * 3 )
    if len(error) == 6:
        good_tf = np.abs(error) < np.array([threshold_xyz] * 3 + [threshold_v] * 3)

    return np.all(good_tf), error
    
    
def similar_angles( ang1_deg, ang2_deg, threshold_arcsec = 1.0):
    '''
    '''
    err_deg    = ang1_deg-ang2_deg
    err_arcsec = err_deg * 3600.
    good= np.abs(err_arcsec) < np.array([threshold_arcsec])
    return np.all(good), err_arcsec

def read_rwo_for_mpcorb_json( mpc_orb_json_file ):
    '''
    Given an mpcorb_json file in the dev_data/... directory, find the associated
    orbfit RWO file, and then parse the file, extracting arrays of ...:
     - times
     - RA
     - DEC
     - RA_resid
     - DEC_resid
     - obsCode
    '''
    
    # Get RWO name from JSON name
    json_dir, json_file = os.path.split(mpc_orb_json_file)
    rwo_filepath        = os.path.join(json_dir ,json_file.replace(".json",".rwo") )
    
    # Read data from RWO file
    with open(rwo_filepath) as json_file:
        data = json.load(json_file)

    # Get arrays of interest
    ymd , ra, dec, ra_resid_arcsec, dec_resid_arcsec, obscode = [],[],[],[],[],[]
    for _ in data["optical_list"]:
        ymd.append( ( int(_["year"]),    int(_["month"]),   float(_["day"]) ) )
        ra.append(  ( float(_["ra_hrs"]),  float(_["ra_min"]) , float(_["ra_sec"]) ) )
        dec.append( ( float(_["dec_deg"]), float(_["dec_min"]), float(_["dec_sec"]) ) )
        ra_resid_arcsec.append(  float(_["ra_resid"])  )
        dec_resid_arcsec.append( float(_["dec_resid"]) )
        obscode.append( _["obscode"])

    # Convert stupid fucking formats
    def convert_time(yy, mm, dd):
        #Convert fraction of the day in hh:mm:ss
        dd_split = math.modf(dd)
        dd_int = int(dd_split[1])
        hours = dd_split[0]
        hours = hours*24
        hours_split = math.modf(hours)
        hours_int = int(hours_split[1])
        minutes = hours_split[0]
        minutes = minutes*60
        minutes_split = math.modf(minutes)
        minutes_int = int(minutes_split[1])
        seconds = minutes_split[0]
        seconds = seconds*60
        seconds_split = math.modf(seconds)
        seconds_int = int(seconds_split[1])
        seconds = int(str(seconds-int(seconds))[2:8])
        #Get MJD
        times = Time(datetime(yy,mm,dd_int,hours_int,minutes_int,seconds_int,seconds),scale='utc')
        jd = times.tdb.jd
        return dd_int, hours_int, minutes_int, seconds_int, jd
    def convert_ra(hrs,mins,secs):
        return 15.0 * (float(hrs) + float(mins)/60.0 + float(secs)/3600.0)
    def convert_dec(degs, mins, secs):
        return float(degs) + float(mins)/60.0 + float(secs)/3600.0

    times, ra_deg, dec_deg = [],[],[]
    for YMD, RA, DEC in zip(ymd , ra, dec):
        times.append( convert_time( *YMD)[-1] )
        ra_deg.append( convert_ra( *RA) )
        dec_deg.append( convert_dec( *DEC ) )
        
    # Make into arrays
    # Return
    return np.array(times) , np.array(ra_deg) , np.array(dec_deg) , np.array(ra_resid_arcsec) , np.array(dec_resid_arcsec) , np.array(obscode)






# Actual tests
# -----------------------------------------------------------------


def test_accuracy_cartesians_A(  ):
    '''
        Test the ACCURACY of the cheby coefficients
        
        Use a similar approach to that taken in the
        test_nbody_run.test_production_integration_function_wrapper_D function
         - Start from exactly the same position as Horizons.
        
        Here we focus on 519104 == 2010 LE128
    '''

    # Define the variables that will be used in the query
    target  = '519104' # Asteroid #123456
    centre  = '500@0'  # <<-- Barycentric
    epochs  = ['2458850.0','2458880.0']
    id_type = 'smallbody'
    refplane='earth' # <<--Equatorial

    # Call the *nice_Horizons* function to get the cartesian states at the first time
    # This is returning EQUATORIAL BERYCENTRIC coordinates
    horizons_zero = Horizons.nice_Horizons(target, centre, epochs[0], id_type, refplane=refplane)
    
    # Call the production_integration_function_wrapper
    # - We already tested this in (e.g.) test_nbody_run.test_production_integration_function_wrapper_D
    # - But there's no harm in repeating the test on a different object here
    # The integration is performed in EQUATORIAL BERYCENTRIC coordinates
    tstart = epoch = epochs[0]
    tstop  = epochs[-1]
    instates = np.array([ horizons_zero ])
    outtimes, states, partial_derivatives_wrt_state, partial_derivatives_wrt_NG, return_value = ephem_forces.production_integration_function_wrapper(   float(tstart),
                                                                float(tstop),
                                                                float(epoch),
                                                                instates,
                                                                non_grav_dict_list = None,
                                                                geocentric = 0 ,
                                                                epsilon = 1e-10,
                                                                tstep_min = 0.02,
                                                                tstep_max = 32.)
                                                                
    # Now call horizons again at some of the output times at which the *production_integration_function_wrapper()* produced output
    for n, t in enumerate(outtimes):
    
        # Check horizons at the simulation-time output
        h = Horizons.nice_Horizons(target, centre, t, id_type, refplane=refplane)

        # Check similarity: threshold_xyz=1e-12, threshold_v=1e-13 : # 15 cm, 1.5 cm/day
        similar_bool , error = similar_xyzuvw(h, states[n][0], threshold_xyz=1e-11, threshold_v=1e-12)
        assert similar_bool
         

    # ----------------------------------------------------------------
    # Now transition to testing the accuracy of the cheby coefficients
    # - The first dozen-or-so lines are required to do set-up
    # - In practice these will all be handled by convenience functions
    # ----------------------------------------------------------------
    
    # Instantiate an MSC object
    M=orbit_cheby.MSC()

    # ## ### #### #####
    # Run the *from_coord_arrays* function that we want to test ...
    # NB1: This attempts to load from arrays ...
    # NB2: Need to do array slicing that will be taken care of within MSC_Loader in practical operation
    # ## ### #### #####
    primary_unpacked_provisional_designation = target
    M.from_coord_arrays(primary_unpacked_provisional_designation, outtimes , states[:,0,:] )

    # check that the expected attributes have been set
    for attr in ["TDB_init", "TDB_final", "sector_init", "sector_final", "sector_coeffs"]:
        assert attr in M.__dict__


    # Check the settings used for the cheby-fitting
    assert M.minorder == 5 , 'The default is expected to be *5*, albeit there is little reasoning for why ...'

    # Check that the sector_coeff dict is of the correct shape
    for sector_number , cheb_coeffs in M.sector_coeffs.items():

        # The # of coefficients for each sector should be > minorder
        #      ( if order = 7 , N_coefficients = 8)
        # NB> v.shape = ( N_coefficients , N_coords  ) , and N_coords = 6
        assert cheb_coeffs.shape[0] > M.minorder
        assert cheb_coeffs.shape[1] , N.output_states[:,0,:].shape[1]
            
            
    # -------------------------------------------------------------------
    # *** TESTING THE ACCURACY OF THE *evaluate_components* FUNCTION ***
    # Evaluate cheby at given time: returns array of shape (N_times, N_components)
    chebyEval = M.evaluate_components(outtimes[:-1])  # Not evaluating the last one as the sector is not supported
    assert chebyEval.shape == ( len(outtimes[:-1]), 6 ) # No cov supplied, so N_components = 6

    # At each time, that the cheby-evaluation is "close-enough" to the input state
    for n, t in enumerate(outtimes[:-1]):             # Not evaluating the last one as the sector is not supported
        similarBool, err_arr = similar_xyzuvw(states[n][0] , chebyEval[n] , threshold_xyz=1e-10, threshold_v=1e-11)
        assert similarBool , f'inputState={inputState} , chebyEval={chebyEval} , err_arr={err_arr}'
    # --------------------------------------------------------

    
    # -------------------------------------------------------------------
    # *** TESTING THE ACCURACY OF THE *generate_XYZ* FUNCTION ***
    # Evaluate cheby at given time: returns array of shape (N_times, 6)
    chebyEval = M.generate_XYZ(outtimes[:-1])  # Not evaluating the last one as the sector is not supported
    assert chebyEval.shape == ( len(outtimes[:-1]), 3 )

    # At each time, that the cheby-evaluation is "close-enough" to the input state
    for n, t in enumerate(outtimes[:-1]):             # Not evaluating the last one as the sector is not supported
        similarBool, err_arr = similar_xyzuvw(states[n][0][:3] , chebyEval[n] , threshold_xyz=1e-10, threshold_v=1e-11)
        assert similarBool , f'inputState={inputState} , chebyEval={chebyEval} , err_arr={err_arr}'
    # --------------------------------------------------------
    

    # -------------------------------------------------------------------
    # *** TESTING THE ACCURACY OF THE *generate_XYZUVW* FUNCTION ***
    # Evaluate cheby at given time: returns array of shape (N_times, 6)
    chebyEval = M.generate_XYZUVW(outtimes[:-1])  # Not evaluating the last one as the sector is not supported
    assert chebyEval.shape == ( len(outtimes[:-1]), 6 )

    # At each time, that the cheby-evaluation is "close-enough" to the input state
    for n, t in enumerate(outtimes[:-1]):  # Not evaluating the last one as the sector is not supported
        similarBool, err_arr = similar_xyzuvw(states[n][0] , chebyEval[n] , threshold_xyz=1e-10, threshold_v=1e-11)
        assert similarBool , f'inputState={inputState} , chebyEval={chebyEval} , err_arr={err_arr}'
    # --------------------------------------------------------
    


def test_accuracy_cartesians_B(  ):
    '''
        Test the ACCURACY of the cheby coefficients
        
        EXACTLY THE SAME AS *test_accuracy_cartesians_A* EXCEPT THAT
        HERE WE FOCUS ON 719 == Albert
    '''

    # Define the variables that will be used in the query
    target  = '719' # Asteroid #123456
    centre  = '500@0'  # <<-- Barycentric
    epochs  = ['2458850.0','2458880.0']
    id_type = 'smallbody'
    refplane='earth' # <<--Equatorial

    # Call the *nice_Horizons* function to get the cartesian states at the first time
    # This is returning EQUATORIAL BERYCENTRIC coordinates
    horizons_zero = Horizons.nice_Horizons(target, centre, epochs[0], id_type, refplane=refplane)
    
    # Call the production_integration_function_wrapper
    # - We already tested this in (e.g.) test_nbody_run.test_production_integration_function_wrapper_D
    # - But there's no harm in repeating the test on a different object here
    # The integration is performed in EQUATORIAL BERYCENTRIC coordinates
    tstart = epoch = epochs[0]
    tstop  = epochs[-1]
    instates = np.array([ horizons_zero ])
    outtimes, states, partial_derivatives_wrt_state, partial_derivatives_wrt_NG, return_value = ephem_forces.production_integration_function_wrapper(   float(tstart),
                                                                float(tstop),
                                                                float(epoch),
                                                                instates,
                                                                non_grav_dict_list = None,
                                                                geocentric = 0 ,
                                                                epsilon = 1e-10,
                                                                tstep_min = 0.02,
                                                                tstep_max = 32.)
                                                                
    # Now call horizons again at some of the output times at which the *production_integration_function_wrapper()* produced output
    for n, t in enumerate(outtimes):
    
        # Check horizons at the simulation-time output
        h = Horizons.nice_Horizons(target, centre, t, id_type, refplane=refplane)

        # Check similarity: threshold_xyz=1e-12, threshold_v=1e-13 : # 15 cm, 1.5 cm/day
        similar_bool , error = similar_xyzuvw(h, states[n][0], threshold_xyz=1e-11, threshold_v=1e-12)
        assert similar_bool
         

    # ----------------------------------------------------------------
    # Now transition to testing the accuracy of the cheby coefficients
    # - The first dozen-or-so lines are required to do set-up
    # - In practice these will all be handled by convenience functions
    # ----------------------------------------------------------------
    
    # Instantiate an MSC object
    M=orbit_cheby.MSC()

    # ## ### #### #####
    # Run the *from_coord_arrays* function that we want to test ...
    # NB1: This attempts to load from arrays ...
    # NB2: Need to do array slicing that will be taken care of within MSC_Loader in practical operation
    # ## ### #### #####
    primary_unpacked_provisional_designation = target
    M.from_coord_arrays(primary_unpacked_provisional_designation, outtimes , states[:,0,:] )

    # check that the expected attributes have been set
    for attr in ["TDB_init", "TDB_final", "sector_init", "sector_final", "sector_coeffs"]:
        assert attr in M.__dict__


    # Check the settings used for the cheby-fitting
    assert M.minorder == 5 , 'The default is expected to be *5*, albeit there is little reasoning for why ...'

    # Check that the sector_coeff dict is of the correct shape
    for sector_number , cheb_coeffs in M.sector_coeffs.items():

        # The # of coefficients for each sector should be > minorder
        #      ( if order = 7 , N_coefficients = 8)
        # NB> v.shape = ( N_coefficients , N_coords  ) , and N_coords = 6
        assert cheb_coeffs.shape[0] > M.minorder
        assert cheb_coeffs.shape[1] , N.output_states[:,0,:].shape[1]
            
            
    # -------------------------------------------------------------------
    # *** TESTING THE ACCURACY OF THE *evaluate_components* FUNCTION ***
    # Evaluate cheby at given time: returns array of shape (N_times, N_components)
    chebyEval = M.evaluate_components(outtimes[:-1])  # Not evaluating the last one as the sector is not supported
    assert chebyEval.shape == ( len(outtimes[:-1]), 6 ) # No cov supplied, so N_components = 6

    # At each time, that the cheby-evaluation is "close-enough" to the input state
    for n, t in enumerate(outtimes[:-1]):             # Not evaluating the last one as the sector is not supported
        similarBool, err_arr = similar_xyzuvw(states[n][0] , chebyEval[n] , threshold_xyz=1e-10, threshold_v=1e-11)
        assert similarBool , f'inputState={inputState} , chebyEval={chebyEval} , err_arr={err_arr}'
    # --------------------------------------------------------

    
    # -------------------------------------------------------------------
    # *** TESTING THE ACCURACY OF THE *generate_XYZ* FUNCTION ***
    # Evaluate cheby at given time: returns array of shape (N_times, 6)
    chebyEval = M.generate_XYZ(outtimes[:-1])  # Not evaluating the last one as the sector is not supported
    assert chebyEval.shape == ( len(outtimes[:-1]), 3 )

    # At each time, that the cheby-evaluation is "close-enough" to the input state
    for n, t in enumerate(outtimes[:-1]):             # Not evaluating the last one as the sector is not supported
        similarBool, err_arr = similar_xyzuvw(states[n][0][:3] , chebyEval[n] , threshold_xyz=1e-10, threshold_v=1e-11)
        assert similarBool , f'inputState={inputState} , chebyEval={chebyEval} , err_arr={err_arr}'
    # --------------------------------------------------------
    

    # -------------------------------------------------------------------
    # *** TESTING THE ACCURACY OF THE *generate_XYZUVW* FUNCTION ***
    # Evaluate cheby at given time: returns array of shape (N_times, 6)
    chebyEval = M.generate_XYZUVW(outtimes[:-1])  # Not evaluating the last one as the sector is not supported
    assert chebyEval.shape == ( len(outtimes[:-1]), 6 )

    # At each time, that the cheby-evaluation is "close-enough" to the input state
    for n, t in enumerate(outtimes[:-1]):  # Not evaluating the last one as the sector is not supported
        similarBool, err_arr = similar_xyzuvw(states[n][0] , chebyEval[n] , threshold_xyz=1e-10, threshold_v=1e-11)
        assert similarBool , f'inputState={inputState} , chebyEval={chebyEval} , err_arr={err_arr}'
    # --------------------------------------------------------
    


def test_accuracy_RADEC_A(  ):
    '''
        Test the ACCURACY of the RADEC coefficients
        
        Use a similar approach to that taken in the
        test_nbody_run.test_production_integration_function_wrapper_D function
         - Start from exactly the same position as Horizons.
        
        Here we focus on 519104 == 2010 LE128
    '''

    # Define the variables that will be used in the query
    target  = '519104' # Asteroid #123456
    centre  = '500@0'  # <<-- Barycentric
    epochs  = ['2458850.0','2458880.0']
    id_type = 'smallbody'
    refplane='earth' # <<--Equatorial

    # Call the *nice_Horizons* function to get the cartesian states at the first time
    # This is returning EQUATORIAL BERYCENTRIC coordinates
    horizons_zero = Horizons.nice_Horizons(target, centre, epochs[0], id_type, refplane=refplane)
    
    # Call the production_integration_function_wrapper
    # - We already tested this in (e.g.) test_nbody_run.test_production_integration_function_wrapper_D
    # - But there's no harm in repeating the test on a different object here
    # The integration is performed in EQUATORIAL BERYCENTRIC coordinates
    tstart = epoch = epochs[0]
    tstop  = epochs[-1]
    instates = np.array([ horizons_zero ])
    outtimes, states, partial_derivatives_wrt_state, partial_derivatives_wrt_NG, return_value = ephem_forces.production_integration_function_wrapper(   float(tstart),
                                                                float(tstop),
                                                                float(epoch),
                                                                instates,
                                                                non_grav_dict_list = None,
                                                                geocentric = 0 ,
                                                                epsilon = 1e-10,
                                                                tstep_min = 0.02,
                                                                tstep_max = 32.)
                                                                
    # Now call horizons again at some of the output times at which the *production_integration_function_wrapper()* produced output
    for n, t in enumerate(outtimes):
    
        # Check horizons at the simulation-time output
        h = Horizons.nice_Horizons(target, centre, t, id_type, refplane=refplane)

        # Check similarity: threshold_xyz=1e-12, threshold_v=1e-13 : # 15 cm, 1.5 cm/day
        similar_bool , error = similar_xyzuvw(h, states[n][0], threshold_xyz=1e-11, threshold_v=1e-12)
        assert similar_bool

         

    # ----------------------------------------------------------------
    # Now transition to testing the accuracy of the RADEC coefficients
    # - The first dozen-or-so lines are required to do set-up
    # - In practice these will all be handled by convenience functions
    # ----------------------------------------------------------------
    
    # Instantiate an MSC object
    M=orbit_cheby.MSC()

    # ## ### #### #####
    # Run the *from_coord_arrays* function that we want to test ...
    # NB1: This attempts to load from arrays ...
    # NB2: Need to do array slicing that will be taken care of within MSC_Loader in practical operation
    # ## ### #### #####
    primary_unpacked_provisional_designation = target
    M.from_coord_arrays(primary_unpacked_provisional_designation, outtimes , states[:,0,:] )

    # check that the expected attributes have been set
    for attr in ["TDB_init", "TDB_final", "sector_init", "sector_final", "sector_coeffs"]:
        assert attr in M.__dict__


    # Check the settings used for the cheby-fitting
    assert M.minorder == 5 , 'The default is expected to be *5*, albeit there is little reasoning for why ...'

    # Check that the sector_coeff dict is of the correct shape
    for sector_number , cheb_coeffs in M.sector_coeffs.items():

        # The # of coefficients for each sector should be > minorder
        #      ( if order = 7 , N_coefficients = 8)
        # NB> v.shape = ( N_coefficients , N_coords  ) , and N_coords = 6
        assert cheb_coeffs.shape[0] > M.minorder
        assert cheb_coeffs.shape[1] , N.output_states[:,0,:].shape[1]
            
            

    # We will need some observatory positions to provide as inputs to the RA,Dec function ...
    # I know from test_obs_pos.test_get_barycentric_A() that the observatory positions are consistent with JPL at the ~15m (1e-10) level
    # But remember that the obs_pos code takes UTC not TDB
    obsCode = 'F51'
    t = Time( outtimes[:-1] , format='jd', scale='tdb' )
    t = t.utc  # <<-- This converts to utc (from tdb, above)
    helio_eq_observatoryXYZ = np.array( [obs_pos.ObsPos().get_heliocentric_equatorial_xyz(jd , obsCode=obsCode) for jd in t.jd ] )
    bary_eq_observatoryXYZ  = coco.equatorial_helio2bary( helio_eq_observatoryXYZ ,outtimes[:-1])
    
    # -------------------------------------------------------------------
    # *** TESTING THE ACCURACY OF THE *generate_RaDec* FUNCTION ***
    chebyEval_XYZ   = M.generate_XYZ(outtimes[:-1])  # Not evaluating the last one as the sector is not supported
    chebyEval_Unit  = M.generate_UnitVector(outtimes[:-1], bary_eq_observatoryXYZ)
    chebyEval_RADEC = M.generate_RaDec(outtimes[:-1], bary_eq_observatoryXYZ )  # Not evaluating the last one as the sector is not supported
    assert chebyEval_XYZ.shape == ( len(outtimes[:-1]), 3 )
    assert chebyEval_Unit.shape == ( len(outtimes[:-1]), 3 )
    assert chebyEval_RADEC.shape == ( len(outtimes[:-1]), 2 )
    # -------------------------------------------------------------------



    # At each time, test that the cheby-evaluation is "close-enough" to the input state
    for n, tdb in enumerate(outtimes[:-1]):             # Not evaluating the last one as the sector is not supported
        #print(n,tdb)
        
        # ---- 1 : --------------------------------
        # Repeating the check from *test_obs_pos.test_get_barycentric_A* that the observatory positions agree with horizons ...
        # OBSERVATORY POSITIONS AGREE AT THE 1e-10 => 15 m LEVEL
        # Query Horizons for barycentric posn of F51
        # Won't work directly, so have to put these hacks in place ...
        # (a) get posn of heliocentre w.r.t. barycenter, and then add to posn of F51 w.r.t. heliocenter
        # (b) setting the target as the Sun, and the center as the observatory)
        # - as such, we need to multiply the values by -1
        target = '10' ; cent = 'F51' ; id_type = 'majorbody' ; refplane = 'earth'
        refplane='earth'
        horizons_helio_eq = Horizons.nice_Horizons(target, cent, str(tdb), id_type, refplane=refplane )
        horizons_helio_eq = -1*horizons_helio_eq
        b2h = Horizons.nice_Horizons(target, '@0', str(tdb), id_type, refplane=refplane )
        horizons_bary_eq_obs = b2h + horizons_helio_eq
        similarBool, err_arr = similar_xyzuvw(horizons_bary_eq_obs[:3] , bary_eq_observatoryXYZ[n] , threshold_xyz=1e-10, threshold_v=1e-10)
        '''
        print('---'*22, 1)
        print('horizons_bary_eq_obs=', horizons_bary_eq_obs)
        print('bary_eq_observatoryXYZ[n]=', bary_eq_observatoryXYZ[n])
        print('err',err_arr)
        '''
        assert similarBool , f'horizons_bary_eq_obs[:3]={horizons_bary_eq_obs[:3]} , bary_eq_observatoryXYZ[n]={bary_eq_observatoryXYZ[n]} , err_arr={err_arr}'


        # ---- 2 : --------------------------------
        # Repeating the check from *test_accuracy_cartesians_B* that the cheby-eval of the XYZ agrees with horizons ...
        # BARYCENTRIC EQUATORIAL POSITIONS OF TTHE OBJECT AGREE AT THE 1e-12 => 15cm LEVEL
        target = '519104' ; cent = '500@0' ; id_type = 'smallbody' ; refplane = 'earth'
        h_bary_eq = Horizons.nice_Horizons(target, cent, tdb, id_type, refplane=refplane)
        similarBool, err_arr = similar_xyzuvw(h_bary_eq[:3] , chebyEval_XYZ[n] , threshold_xyz=1e-12, threshold_v=1e-12)
        '''
        print('---'*22, 2)
        print('h_bary_eq=', h_bary_eq)
        print('chebyEval_XYZ[n]=', chebyEval_XYZ[n])
        print('err',err_arr)
        '''
        assert similarBool , f'h_bary_eq[:3]={h_bary_eq[:3]} , chebyEval_XYZ[n]={chebyEval_XYZ[n]} , err_arr={err_arr}'


        # ---- 3 : --------------------------------
        # Doing the RA-DEC comparison with horizons
        target = '519104' ; cent = 'F51' ; id_type = 'smallbody' ; refplane = 'earth'
        
        # ARRGGHH: Fucking assumed time format changes between vectors & ephemerides calls !!!
        tmp_t = Time( tdb , format='jd', scale='tdb' )
        tmp_t = tmp_t.utc  # <<-- This converts to utc (from tdb, above)
        utc = tmp_t.jd
        h_RADEC = Horizons.nice_Horizons_radec(target, cent, utc, id_type, refplane=refplane)
        similarBool, err_arr = similar_angles(  h_RADEC[0], chebyEval_RADEC[n], threshold_arcsec = 1.0)
        '''
        print('---'*22, 4)
        print('h_RADEC=', h_RADEC)
        print('chebyEval_RADEC[n]=', chebyEval_RADEC[n])
        '''
        print(similarBool, err_arr)
        assert similarBool , f'h_RADEC[n]={h_RADEC[n]} , chebyEval_RADEC[0]={chebyEval_RADEC[0]} , err_arr={err_arr}'
  

     



"""
def test_accuracy_RADEC_A(  ):
    '''
        Test the ACCURACY of the RA,Dec calculated from cheby coefficients
        
        Use a similar approach to that taken in the
        test_nbody_run.test_production_integration_function_wrapper_D function
         - Start from exactly the same position as Horizons.
        
        Here we focus on 519104 == 2010 LE128
    '''

    # Define the variables that will be used in the query
    target  = '519104' # Asteroid #123456
    centre  = '500@0'  # <<-- Barycentric
    epochs  = ['2458850.0','2458880.0']
    id_type = 'smallbody'
    refplane='earth' # <<--Equatorial

    # Call the *nice_Horizons* function to get the cartesian states at the first time
    # This is returning EQUATORIAL BERYCENTRIC coordinates
    horizons_zero = Horizons.nice_Horizons(target, centre, epochs[0], id_type, refplane=refplane)
    print('horizons_zero=',horizons_zero)
    
    # Call the production_integration_function_wrapper
    # - We already tested this in (e.g.) test_nbody_run.test_production_integration_function_wrapper_D
    # - But there's no harm in repeating the test on a different object here
    # The integration is performed in EQUATORIAL BERYCENTRIC coordinates
    tstart = epoch = epochs[0]
    tstop  = epochs[-1]
    instates = np.array([ horizons_zero ])
    outtimes, states, partial_derivatives_wrt_state, partial_derivatives_wrt_NG, return_value = ephem_forces.production_integration_function_wrapper(   float(tstart),
                                                                float(tstop),
                                                                float(epoch),
                                                                instates,
                                                                non_grav_dict_list = None,
                                                                geocentric = 0 ,
                                                                epsilon = 1e-10,
                                                                tstep_min = 0.02,
                                                                tstep_max = 1.)
                                                                
    # Now call horizons again at some of the output times at which the *production_integration_function_wrapper()* produced output
    for n, t in enumerate(outtimes):
    
        # Check horizons at the simulation-time output
        h = Horizons.nice_Horizons(target, centre, t, id_type, refplane=refplane)

        # Check similarity: threshold_xyz=1e-12, threshold_v=1e-13 : # 15 cm, 1.5 cm/day
        similar_bool , error = similar_xyzuvw(h, states[n][0], threshold_xyz=1e-12, threshold_v=1e-13)
        assert similar_bool
         

    # ----------------------------------------------------------------
    # Now transition to testing the accuracy of the RA,Dec components
    # calculated using chebys
    # - The first dozen-or-so lines are required to do set-up
    # - In practice these will all be handled by convenience functions
    # ----------------------------------------------------------------
    
    # Instantiate an MSC object
    M=orbit_cheby.MSC()

    # ## ### #### #####
    # Run the *from_coord_arrays* function that we want to test ...
    # NB1: This attempts to load from arrays ...
    # NB2: Need to do array slicing that will be taken care of within MSC_Loader in practical operation
    # ## ### #### #####
    M.from_coord_arrays(primary_unpacked_provisional_designation, outtimes , states )

    # check that the expected attributes have been set
    for attr in ["TDB_init", "TDB_final", "sector_init", "sector_final", "sector_coeffs"]:
        assert attr in M.__dict__

    # check that the "supported_sector_numbers" are as expected
    # NB: convenience_call_to_nbody_run_mpcorb uses
    #     tstart = 2459200 = First date in Sector # 600 (check using B.map_JD_to_sector_number_and_sector_start_JD(2459200, B.standard_MJDmin )
    #     tstop = 2459295 = :Last day in Sector # 602   (check using B.map_JD_to_sector_number_and_sector_start_JD(2459295, B.standard_MJDmin )
    # Hence expected_supported_sector_numbers = [600,601,602]
    B = cheby_checker.Base()
    assert M.sector_init  == B.map_JD_to_sector_number_and_sector_start_JD(N.tstart, B.standard_MJDmin )[0]
    assert M.sector_final == B.map_JD_to_sector_number_and_sector_start_JD(N.tstop, B.standard_MJDmin )[0]
    assert M.TDB_init     == B.map_JD_to_sector_number_and_sector_start_JD(N.tstart, B.standard_MJDmin )[1]
    assert M.TDB_final    == B.map_JD_to_sector_number_and_sector_start_JD(N.tstop, B.standard_MJDmin )[1] + B.sector_length_days - B.epsilon

    # Check the settings used for the cheby-fitting
    #print('B.standard_MJDmin =', B.standard_MJDmin )
    assert M.minorder == 5 , 'The default is expected to be *5*, albeit there is little reasoning for why ...'

    # Check that the sector_coeff dict is of the correct shape
    for sector_number , cheb_coeffs in M.sector_coeffs.items():

        # The # of coefficients for each sector should be > minorder
        #      ( if order = 7 , N_coefficients = 8)
        # NB> v.shape = ( N_coefficients , N_coords  ) , and N_coords = 6
        assert cheb_coeffs.shape[0] > M.minorder
        assert cheb_coeffs.shape[1] , N.output_states[:,0,:].shape[1]
            
            
    # Check the accuracy of the RA,DEC coordinates ...
    for n, t in enumerate(outtimes): # Loop over the times from the NBody integration output ...
        
        # Query Horizons for their RA,Dec
"""
        

"""
def test_evaluate_components_A():
    '''
        Use MSC function to evaluate cheby components at set of supplied times
        
        Here I am checking that the *evaluate_components* function
        returns coordinate values similar to those from rebounx:
         - So in essence I am testing the round-trip of calling "generate_cheb_for_sector" + "evaluate_components"
        
        ***Need to test the hell out of this by doing more tests below, as it underlies everything***
    '''

    # The file that we will use as the source of the data
    mpc_orb_json_filepath = mpc_orb_json_files[0]
    
    # Do a local read of the json (using MPCORB) as it's useful to be able to read the name of the object
    MO = MPCORB(mpc_orb_json_filepath)
    primary_unpacked_provisional_designation = MO.designation_data["unpacked_primary_provisional_designation"]
    
    # Instantiate Nbody Object, Run NBody Simulation, and return populated Object
    N = convenience_call_to_nbody_run_mpcorb(mpc_orb_json_filepath)
    
    # Instantiate MSC object & populate using data from Nbody Object ...
    M=orbit_cheby.MSC()
    M.from_coord_arrays(primary_unpacked_provisional_designation, N.output_times , N.output_states[:,0,:] )

    # Indicees for data supported by this MSC:
    init, final = M.get_valid_range_of_dates()
    ind = np.where( (N.output_times >= init) & (N.output_times <= final ) )

    # Make the call to the low-level function, *evaluate_components*
    # NB: returns array of shape (N_times, N_components)
    evaluatedComponents      = M.evaluate_components( N.output_times[ind] )
    evaluatedComponentsSlice = M.evaluate_components( N.output_times[ind] , component_slice_spec=slice(0,3))

    # Make a call to the higher level function, *generate_XYZ*
    # NB returned result is of shape (len(times_tdb) , 3)
    XYZs = M.generate_XYZ( N.output_times[ind] )

    # Check that the shape of the returned components is as expected ...
    assert N.output_states[:,0,:].shape   == evaluatedComponents.shape
    assert evaluatedComponentsSlice.shape == XYZs.shape

    # Check that the values in evaluatedComponentsSlice == those in XYZs
    assert np.all( evaluatedComponentsSlice == XYZs )

    # Check that the values in evaluatedComponentsSlice are all close to the input state
    for n, comp in enumerate( N.output_states[:,0,:] ):
    
        # Input xyz is the first 3 components of N.output_states[:,0,:]
        input_xyz = comp[:3]
        output_xyz = evaluatedComponentsSlice[n]
        similarBool, err_arr = similar_xyzuvw(input_xyz, output_xyz , threshold_xyz=1e-11, threshold_v=1e-11)
        assert similarBool , f'input_xyz={input_xyz} , output_xyz={output_xyz} , err_arr={err_arr}'



def test_generate_RaDec_A():
    '''
    
    Test the generate of RaDec coordinates from an MSC (Cheby) Object
    
    Note (1) that under the hood this calls *generate_UnitVector*, so in
    order to compare to the RAs & DECs used in the orbit fitting, so have
    to call generate_RaDec which relies on generate_UnitVector
    
    Note (2) that at this point I am switching to using MSCs generated using the
    convenience function *convenience_call_to_get_MSCs_from_file_via_nbody_run( json_filepath_or_list)*
     - This uses MSC_Loader under the hood
     - MSC_Loader is tested below ... see, e.g., *test_create_loader*  &  *test_loader_from_nbodysim*
    
    '''

    # The file that we will use as the source of the data
    mpc_orb_json_filepath = mpc_orb_json_files[0]
    print(mpc_orb_json_filepath)
    
    # Getting MSCs via convenience call to MSC_Loader, NBodySim, etc
    B = cheby_checker.Base()
    MSCs = convenience_call_to_get_MSCs_from_file_via_nbody_run( mpc_orb_json_filepath , tstart=B.standard_MJDmin , tstop=B.standard_MJDmax )
    MSC = MSCs[0]
    
    # Read in orbfit RWO file of positions and residuals for the orbit in question
    times_tdb_jd , ra_deg , dec_deg , ra_resid_arcsec , dec_resid_arcsec , obscode =  read_rwo_for_mpcorb_json( mpc_orb_json_filepath )
    
    # We will need some observatory positions to provide as inputs ...
    n_samples = 300
    helio_eq_observatoryXYZ = np.array( [obs_pos.ObsPos().get_heliocentric_equatorial_xyz(jd , obsCode=obsCode) for jd,obsCode in zip(times_tdb_jd[-n_samples:],obscode[-n_samples:])] )

    # conversion to barycentric equatorial coordinates
    bary_eq_observatoryXYZ = coco.equatorial_helio2bary( helio_eq_observatoryXYZ ,times_tdb_jd[-n_samples:])

    # Simulate RA,Dec observations at the times of the actual observations from the rwo file
    # *** This is the function we are testing ***
    RaDEC = MSC.generate_RaDec( times_tdb_jd[-n_samples:] , observatoryXYZ=bary_eq_observatoryXYZ )
    
    def get_orbfit_calculated( ra_observed_deg, dec_observed_deg , ra_resid_arcsec , dec_resid_arcsec  ):
        return  ra_observed_deg - (ra_resid_arcsec/3600.)/np.cos(np.radians(dec_observed_deg)) , dec_observed_deg - dec_resid_arcsec/3600.
        
    # Check whether the simulated RA,Dec positions are close to the observed RA,Dec
    for n, t in enumerate(times_tdb_jd[-n_samples:]):
        # Back-calculate where orbfit must have calculaetd the object to be ...
        orbfit_calc_ra, orbfit_calc_dec  = get_orbfit_calculated(   ra_deg[-n_samples:][n],
                                                                    dec_deg[-n_samples:][n],
                                                                    ra_resid_arcsec[-n_samples:][n],
                                                                    dec_resid_arcsec[-n_samples:][n])
        # Check whether the cheby-calc is similar to the orbfit calc ...
        good, err_arcsec = similar_angles( np.array( [orbfit_calc_ra, orbfit_calc_dec] ), RaDEC[n], threshold_arcsec = 2.0)
        assert good, f'Calc RA,Dec disagreed by >2": err_arcsec:{err_arcsec} , orbfit:{orbfit_calc_ra, orbfit_calc_dec  } , cheby:{RaDEC[n]} '
        print(n, t , good , err_arcsec)
    
def test_generate_RaDec_B():
    '''
    
    Test the generate of RaDec coordinates from an MSC (Cheby) Object
    
    Structured as per test_generate_RaDec_A, but now testing MULTIPLE input json files
    
    '''

    # Loop over different input json files ...
    for mpc_orb_json_filepath in mpc_orb_json_files:
        print(mpc_orb_json_filepath)
        
        # Getting MSCs via convenience call to MSC_Loader, NBodySim, etc
        B = cheby_checker.Base()
        MSCs = convenience_call_to_get_MSCs_from_file_via_nbody_run( mpc_orb_json_filepath , tstart=B.standard_MJDmin , tstop=B.standard_MJDmax )
        MSC = MSCs[0]
        
        # Read in orbfit RWO file of positions and residuals for the orbit in question
        times_tdb_jd , ra_deg , dec_deg , ra_resid_arcsec , dec_resid_arcsec , obscode =  read_rwo_for_mpcorb_json( mpc_orb_json_filepath )
        
        # We will need some observatory positions to provide as inputs ...
        n_samples = 10
        helio_eq_observatoryXYZ = np.array( [obs_pos.ObsPos().get_heliocentric_equatorial_xyz(jd , obsCode=obsCode) for jd,obsCode in zip(times_tdb_jd[-n_samples:],obscode[-n_samples:])] )

        # conversion to barycentric equatorial coordinates
        bary_eq_observatoryXYZ = coco.equatorial_helio2bary( helio_eq_observatoryXYZ ,times_tdb_jd[-n_samples:])

        # Simulate RA,Dec observations at the times of the actual observations from the rwo file
        # *** This is the function we are testing ***
        RaDEC = MSC.generate_RaDec( times_tdb_jd[-n_samples:] , observatoryXYZ=bary_eq_observatoryXYZ )
        
        def get_orbfit_calculated( ra_observed_deg, dec_observed_deg , ra_resid_arcsec , dec_resid_arcsec  ):
            return  ra_observed_deg - (ra_resid_arcsec/3600.)/np.cos(np.radians(dec_observed_deg)) , dec_observed_deg - dec_resid_arcsec/3600.
            
        # Check whether the simulated RA,Dec positions are close to the observed RA,Dec
        for n, t in enumerate(times_tdb_jd[-n_samples:]):
            # Back-calculate where orbfit must have calculaetd the object to be ...
            orbfit_calc_ra, orbfit_calc_dec  = get_orbfit_calculated(   ra_deg[-n_samples:][n],
                                                                        dec_deg[-n_samples:][n],
                                                                        ra_resid_arcsec[-n_samples:][n],
                                                                        dec_resid_arcsec[-n_samples:][n])
                                                                        
            # Check whether the cheby-calc is similar to the orbfit calc ...
            # Note that at least one of the test points differs by ~2" in RA
            good, err_arcsec = similar_angles( np.array( [orbfit_calc_ra, orbfit_calc_dec] ), RaDEC[n], threshold_arcsec = 2.0)
            assert good, f'Calc RA,Dec disagreed by >2": err_arcsec:{err_arcsec} , orbfit:{orbfit_calc_ra, orbfit_calc_dec  } , cheby:{RaDEC[n]} '
            
        
"""


#test_accuracy_cartesians_B()
test_accuracy_RADEC_A()
