# -*- coding: utf-8 -*-
# sifter/sifter/precalc.py

'''
    --------------------------------------------------------------
    sifter's precalculation module.

    Jan 2020
    Matt Payne & Mike Alexandersen

    This module provides functionalities to
    (a) create & save new precalculations
    (b) load extant precalculations

    *WRITE MORE STUFF*

    --------------------------------------------------------------
    '''


# Import third-party packages
# --------------------------------------------------------------
import sys
import os
import getpass
import warnings
from astropy_healpix import HEALPix
from astropy.coordinates import SkyCoord
from astropy import units as u


# Different machines set up differently ...
# ... adding paths to force stuff to work while developing
if getpass.getuser() in ['matthewjohnpayne']:
    sys.path.append('/Users/matthewjohnpayne/Envs/mpcvenv/obs80/')
    from obs80 import parse80
elif getpass.getuser() in ['malexand']:
    sys.path.append('/data/mpc/share/apps/python_libs/mpcpp/obs80/')
    from mpcpp.obs80.obs80 import parse80
else:
    from obs80 import parse80


# Import neighboring packages
# --------------------------------------------------------------
sys.path.append(os.path.join( os.path.dirname(os.path.dirname(os.path.realpath(__file__) )), 'sifter'))
import sql


# Default for caching stuff using lru_cache
# -------------------------------------------------------------


# Data classes/methods
# -------------------------------------------------------------


class Base():
    '''
        Parent class to hold some file/directory definitions & methods
    '''

    def __init__(self):

        # Filing ...
        #self.HP     = '_healpix.txt'

        # Healpix ...
        self.HP_nside = 16
        self.HP_order = 'nested'
        self.HPix = HEALPix(nside=self.HP_nside, order=self.HP_order)
        self.HP_npix = self.HPix.npix

        # sqlite database specs ...
        self.db_filename = 'sifter.db'

    def _fetch_data_directory(self, ):
        '''
            Returns the default path to the directory where data will be
            downloaded.

            By default, this method will return ~/.sifter_data/data
            and create this directory if it does not exist.

            If the directory cannot be accessed or created, then it returns
            the local directory (".")

            Returns
            -------
            data_dir : str
            Path to location of `data_dir` where data (FITs files) will be
            downloaded
        '''

        data_dir = os.path.join(os.path.expanduser('~'), '.sifter_data')

        # If it doesn't exist, make a new data directory
        if not os.path.isdir(data_dir):

            try:
                os.mkdir(data_dir)

            # downloads locally if OS error occurs
            except OSError:
                warnings.warn('Warning: unable to create {}. '
                              'Download directory set to be the current '
                              'working directory instead.'.format(data_dir))
            data_dir = '.'

        return data_dir



class Tracklets(Base):
    '''
        Class to facilitate "precalculations" on a list of tracklets
        '''
    
    def __init__(self, observations=None):
        
        # Give access to "Base" methods
        super().__init__()
        
        # connect to db
        self.conn = sql.create_connection(sql.fetch_db_filepath())
        
        # if observations supplied, process them ...
        if observations is not None:
            self.save_tracklets(*self.parse_observation_lists(observations))

    def parse_observation_lists(self, list_of_observation_pairs):
    '''
        read observational input (probably be in obs80-string formats)
        
        Inputs:
        -------
        list_of_observation_pairs : list-of-tuples/lists?
        - each in obs80 format ???
        
        Returns
        -------
        list of tracklet dictionaries
        - specified as per "parse_observation_pair" function
        '''
            
            return [self.parse_observation_pair(observation_pair)
                    for observation_pair in list_of_observation_pairs]
        
    def parse_observation_pair(self, observation_pair):
        '''
            read observational input (probably be in obs80-string formats)
            
            Inputs:
            -------
            observation_pair : list
            - list of strings containing obs80 lines
            
            Returns
            -------
            JD: integer
            - date of first observation
            HP: integer
            - healpix of first observations
            tracklet_name: string
            - unique name for tracklet; 26 characters
            tracklet_dict: dictionary
            - container for all data
            - contains:
            - JD: integer; date of first observation
            - HP: integer;  Healpix of first observation
            - JD2: integer; date of last observation
            - HP2: integer; Healpix of last observations
            - RoM: astropy Quantity; Rate of motion (angle per time)
            - AoM: astropy Quantity; Angle of motion, positive,
            measured from East towards North.
            - tracklet_name: string; Unique ID for the tracklet
            - observations: list of strings; the input obs80 lines
            - should be everything needed for subsequent detailed calculations
            '''
        # Check number of observations given.
        nobs = len(observation_pair)
        if nobs == 0:
            raise RuntimeError("Received zero observations. Can't real.")
        if nobs == 1:
            print("### WARNING ###\nReceived 1 observation.")
            print("Will proceed, but RoM and AoM will be 0.\n### WARNING ###")
        elif nobs > 2:
            print("Received more than 2 observations.")
            print("Only using first and last one.")
                    else:
            pass
        # Parse obs80 lines
        parsed = [obs for obs in parse80(observation_pair)]

        # Convert nteger julian date
        JDfloat = parsed[0].jdutc
        JDfloat2 = parsed[-1].jdutc
        JD = round(JDfloat)
        JD2 = round(JDfloat2)
        # *** WILL WE EVER NEED INTEGER JD FOR LAST OBS??? ***
        
        # Calculate Rate of Motion and Angle of Motion
        Coord = SkyCoord(parsed[0].ra * u.hourangle, parsed[0].dec * u.deg)
        Coord2 = SkyCoord(parsed[-1].ra * u.hourangle, parsed[-1].dec * u.deg)
        Delta_JD = (JDfloat2 - JDfloat) * u.day
        RoM = Coord.separation(Coord2) / Delta_JD
        AoM = Coord.position_angle(Coord2)
        
        # Find the healpix that the coordinates are in.
        HP = self.HPix.lonlat_to_healpix(Coord.ra, Coord.dec)
        HP2 = self.HPix.lonlat_to_healpix(Coord.ra, Coord.dec)
        # *** WILL WE EVER NEED HP FOR LAST OBS??? ***
        
        # Create a unique ID for the tracklet. This is done by combining
        # the 5 digit number or 6-7 character temporary designation,
        # the julian date of first obs (up to 14 characters),
        # and the 3 digit observatory code. 26 characters total
        namestr = parsed[0].num if parsed[0].num != '' else parsed[0].desig
        namestr = namestr.ljust(7, '_')
        jdstr = str(JDfloat).ljust(14, '_')
        codstr = parsed[0].cod.ljust(3, '_')
        tracklet_name = '{}_{}_{}'.format(namestr, jdstr, codstr)
        
        # We should also pre-calc and save the observatory location
        # for each observation
        print('We should also pre-calc and save the observatory location '
              'for each observation')
            
        # put everything of use in the tracklet_dictionary
        tracklet_dictionary = {'JD': JD,
                                'HP': HP,
                                'JD2': JD2,
                                'HP2': HP2,
                                'tracklet_name': tracklet_name,
                                'RoM': RoM,
                                'AoM': AoM,
                                'observations': observation_pair
                                }

    return tracklet_dictionary


    def save_tracklets(self, tracklet_dictionary_list):
        '''
            This should use the results from parse_observations
            and store them appropriately in a nice file/database structure.
            
            Inputs:
            -------
            JD_list : list-of-integers
            - day
            
            HP_list : list-of-integers
            - healpix
            
            tracklet_name_list: list-of-strings ?
            - unique identifier for tracklet
            
            tracklet_dictionary_list: list-of-dictionaries
            - all data that we want to save for each tracklet
            
            Returns
            -------
        
        '''
        JD = [tracklet_dic['JD'] for tracklet_dic in tracklet_dictionary_list]
        HP = [tracklet_dic['HP'] for tracklet_dic in tracklet_dictionary_list]
        tracklet_name = [tracklet_dic['tracklet_name'] for tracklet_dic in tracklet_dictionary_list]
        
        # upload data (making equal-length lists)
        return sql.upsert_tracklets(self.conn, JD, HP, tracklet_name,
                                    tracklet_dictionary_list)

    def delete_tracklets(self, tracklet_name_list):
        '''
        We need some method to remove tracklets
        '''
        return sql.delete_tracklets(self.conn, tracklet_name_list)


# End of file


"""
class Tracklet(Base):
    '''
        Used to perform "precalculations" on an individual tracklet
    '''

    def __init__(self, observations=None):

        # Give access to "Base" methods
        super().__init__()

        # connect to db
        self.conn = sql.create_connection(sql.fetch_db_filepath())

        # if observations supplied, process them ...
        if observations is not None:
            self.save_tracklet(self.parse_observations(observations))

    def parse_observations(self, observation_pair):
        '''
            read observational input (probably be in obs80-string formats)

            Inputs:
            -------
            observation_pair : list
             - list of strings containing obs80 lines

            Returns
            -------
            tracklet_dict: dictionary
             - container for all data
             - contains:
                 - JD: integer; date of first observation
                 - HP: integer;  Healpix of first observation
                 - JD2: integer; date of last observation
                 - HP2: integer; Healpix of last observations
                 - RoM: astropy Quantity; Rate of motion (angle per time)
                 - AoM: astropy Quantity; Angle of motion, positive,
                                          measured from East towards North.
                 - tracklet_name: string; Unique ID for the tracklet
                 - observations: list of strings; the input obs80 lines
             - should be everything needed for subsequent detailed calculations
            '''
        # Check number of observations given.
        nobs = len(observation_pair)
        if nobs == 0:
            raise RuntimeError("Received zero observations. Can't real.")
        if nobs == 1:
            print("### WARNING ###\nReceived 1 observation.")
            print("Will proceed, but RoM and AoM will be 0.\n### WARNING ###")
        elif nobs > 2:
            print("Received more than 2 observations.")
            print("Only using first and last one.")
        else:
            pass
        # Parse obs80 lines
        parsed = [obs for obs in parse80(observation_pair)]

        # Convert nteger julian date
        JDfloat = parsed[0].jdutc
        JDfloat2 = parsed[-1].jdutc
        JD = round(JDfloat)
        JD2 = round(JDfloat2)
        # *** WILL WE EVER NEED INTEGER JD FOR LAST OBS??? ***

        # Calculate Rate of Motion and Angle of Motion
        Coord = SkyCoord(parsed[0].ra * u.hourangle, parsed[0].dec * u.deg)
        Coord2 = SkyCoord(parsed[-1].ra * u.hourangle, parsed[-1].dec * u.deg)
        Delta_JD = (JDfloat2 - JDfloat) * u.day
        RoM = Coord.separation(Coord2) / Delta_JD
        AoM = Coord.position_angle(Coord2)

        # Find the healpix that the coordinates are in.
        HP = self.HPix.lonlat_to_healpix(Coord.ra, Coord.dec)
        HP2 = self.HPix.lonlat_to_healpix(Coord.ra, Coord.dec)
        # *** WILL WE EVER NEED HP FOR LAST OBS??? ***

        # Create a unique ID for the tracklet. This is done by combining
        # the 5 digit number or 6-7 character temporary designation,
        # the julian date of first obs (up to 14 characters),
        # and the 3 digit observatory code. 26 characters total
        namestr = parsed[0].num if parsed[0].num != '' else parsed[0].desig
        namestr = namestr.ljust(7, '_')
        jdstr = str(JDfloat).ljust(14, '_')
        codstr = parsed[0].cod.ljust(3, '_')
        tracklet_name = '{}_{}_{}'.format(namestr, jdstr, codstr)

        # we should also pre-calc and save the observatory location
        # for each observation
        print('We should also pre-calc and save the observatory location'
              ' for each observation.')

        # put everything of use in the tracklet_dictionary
        tracklet_dictionary = {'JD': JD,
                               'HP': HP,
                               'JD2': JD2,
                               'HP2': HP2,
                               'tracklet_name': tracklet_name,
                               'RoM': RoM,
                               'AoM': AoM,
                               'observations': observation_pair
                               }

        return tracklet_dictionary
        #return JD, HP, tracklet_name, tracklet_dictionary

    def save_tracklet(self, tracklet_dict):
        '''
            This should use the results from parse_observations
            and store them appropriately in a nice file/database structure.

            Inputs:
            -------
            JD : integer
            - day

            HP : integer
            - healpix

            tracklet_name: string ?
             - unique identifier for tracklet

            tracklet_dict:
             - all data that we want to save

            Returns
            -------

        '''

        # upload data
        return sql.upsert_tracklet(self.conn,
                                   tracklet_dict['JD'],
                                   tracklet_dict['HP'],
                                   tracklet_dict['tracklet_name'],
                                   tracklet_dict)

    def delete_tracklet(self, tracklet_name):
        '''
            We need some method to remove a tracklet
        '''
        return sql.delete_tracklet(self.conn, tracklet_name)

"""
