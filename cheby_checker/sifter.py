# -*- coding: utf-8 -*-
# cheby_checker/cheby_checker/sifter.py

"""
--------------------------------------------------------------
sifter module.

Jan 2020
Matt Payne & Mike Alexandersen

This module provides overall access to sifter:
sifter searches the ITF & database for tracklets that match an orbit.


*WRITE MORE STUFF*

--------------------------------------------------------------
"""

# Import third-party packages
# --------------------------------------------------------------
import sys, os
import argparse


# Import neighboring packages
# --------------------------------------------------------------
from . import sifter_query as query  # Yay, python3 relative import :-)

# Set default search parameters
# --------------------------------------------------------------
default_dict = {
    'cons_arcsecs' : 2.0,
    'max_residuals_arcsecs' : 500.0,
    'orbit_file' : 'cheby.json',
}

def sifter():
    # -------- PARSE ANY COMMAND-LINE ARGUMENTS ----------------
    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument('-s', '--screen', dest='print_to_screen', action="store_true",
                            help="e.g.: sifter.py -s. Prints output to screen. \
                                    Default = False.")

    arg_parser.add_argument('-c', '--consistency', dest='cons_arcsecs', default = default_dict['cons_arcsecs'],
                        help="e.g.: sifter.py -c 1.0 . Sets 'consistency' in arc-seconds. \
                                    Default value = default_dict['cons_arcsecs'] = %r ." % default_dict['cons_arcsecs'])

    arg_parser.add_argument('-m', '--maxres', dest='max_residuals_arcsecs', default = default_dict['max_residuals_arcsecs'],
                        help="e.g.: sifter.py -m 1000.0 . Sets 'maximum residuals' in arc-seconds. \
                                    Default value = default_dict['max_residuals_arcsecs'] = %r ." % default_dict['max_residuals_arcsecs'])

    arg_parser.add_argument('-o', '--orbfile', dest='orbit_file', default = default_dict['orbit_file'],
                        help="e.g.: sifter.py -o K20B12Q.json : Sets name of the input file that contains the orbital information, \
                                    N.B. Currently assumes will be JSON format and contain CHEBYSHEV COEFFS. \
                                    Default value = default_dict['orbit_file'] = %r ." % default_dict['orbit_file'])

    # Add an optional additional argument:
    arg_parser.add_argument('arg', nargs='?', default=None)

    # Parse the arguments:
    args = arg_parser.parse_args()

    # **DEBUGGING** Print out variable names & values
    for k,v in args.__dict__.items():
        print(k,v)

    # -------- EXECUTE QUERY (TO SEARCH ITF) ----------------
    #
    # (i) Check the input file of cheby-coeffs exists
    #  - Could/Should extend this to check its actually a valid file ...
    # asset os.path.isfile( args.orbit_file ), \
    #   'Could not proceed with query as %r does not appear to exist ... \n ... perhaps supply a full filepath?' % args.orbit_file
    #
    # (ii) Read the input file of cheby-coeffs into a dictionary
    #with open( args.orbit_file ) as json_file:
    #   cheby_dict = json.load(json_file)

    # (iii) Run the main sifter-query to see if the input orbit matches anything in the ITF
    #query.query( cheby_dict )


if __name__ == '__main__':
    sifter()
