# -*- coding: utf-8 -*-

'''
    --------------------------------------------------------------
    tests of orbit_cheby's base class
    sifter/tests/test_orbit_cheby

    Jan 2020
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
from cheby_checker import nbody
from cheby_checker import orbit_cheby
from cheby_checker import cheby_checker
from cheby_checker import obs_pos
from cheby_checker import coco

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

# Actual tests : MSC Object [ tests of MSC_Loader object are below]
# -----------------------------------------------------------------



def test_from_coord_arrays_B(  ):
    '''
        Use MSC functionality to load from numpy arrays
        
        Significant reliance on *generate_cheb_for_sector()* under-the-hood
         - We test the ACCURACY of the returned coefficients, so we are testing
           both *from_coord_arrays*  &  *from_coord_arrays*
         
        Here we test whether various quanities are populated as expected
        when we pass in data from MANY SINGLE-PARTICLE integrations
        
        Follows tthe same pattern as test_from_coord_arrays_A, above
        
    '''

    # Loop over the files that we will use as the source of the data
    for mpc_orb_json_filepath in mpc_orb_json_files:
        print('\n','--'*22)
        print(mpc_orb_json_filepath)
        # Instantiate Nbody Object, Run NBody Simulation, and return populated Object
        N = convenience_call_to_nbody_run_mpcorb(mpc_orb_json_filepath)
        #print('N.output_times=',N.output_times)
        
        # Do a local read of the json (using MPCORB) as it's useful to be able to read the name of the object
        M = MPCORB(mpc_orb_json_filepath)
        unpacked_primary_provisional_designation = M.designation_data["unpacked_primary_provisional_designation"]

        # Instantiate MSC object
        M=orbit_cheby.MSC()

        # ## ### #### #####
        # Run the *from_coord_arrays* function that we want to test ...
        # NB1: This attempts to load from arrays ...
        # NB2: Need to do array slicing that will be taken care of within MSC_Loader in practical operation
        # ## ### #### #####
        M.from_coord_arrays(unpacked_primary_provisional_designation, N.output_times , N.output_states[:,0,:] )

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
            
            
        # Check the accuracy of the sector_coeff dict coefficients ...
        for n, t in enumerate(N.output_times): # Loop over the times from the NBody integration output ...

            # Input state ...
            inputState = N.output_states[:,0,:][n]
            
            # Evaluate cheby at given time ...
            # NB: returns array of shape (N_times, N_components)
            chebyEval = M.evaluate_components([t])[0]
            
            # Check that the cheby-evaluation is "close-enough" to the input state
            # NB1: threshold_xyz=1e-11 => ~150cm == 1.5m
            # NB2: **** 2022-02-07 *** : For reasons related to lack of control over the number of output points, I have changed
            #      min number of points-per-secttor to 5 (instead of 7) and worsened the allowed accuracy on this tast
            #      from 1e-11 to 1e-10 (1,500cm). I would like to go back to the more stringent consttraint. But in order
            #      to go back, I/we need better control over the output points returned from the production_integration_function_wrapper
            similarBool, err_arr = similar_xyzuvw(inputState,chebyEval , threshold_xyz=1e-10, threshold_v=1e-11)
            assert similarBool , f'inputState={inputState} , chebyEval={chebyEval} , err_arr={err_arr}'



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
    unpacked_primary_provisional_designation = MO.designation_data["unpacked_primary_provisional_designation"]
    
    # Instantiate Nbody Object, Run NBody Simulation, and return populated Object
    N = convenience_call_to_nbody_run_mpcorb(mpc_orb_json_filepath)
    
    # Instantiate MSC object & populate using data from Nbody Object ...
    M=orbit_cheby.MSC()
    M.from_coord_arrays(unpacked_primary_provisional_designation, N.output_times , N.output_states[:,0,:] )

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
            
        
