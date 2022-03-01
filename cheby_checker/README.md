# 2020-08-04
# MJP
# Summary of WIP and To-Do

# To Date / WIP 

- In rough order of usage / data-flow, the following modules should exist ...

### obs_pos.py

- Functionality to allow calculation of observatory positions. 
- It's not at all efficient, but I want to have the function outside of any of the other modules
- Initially written to be dependent on naval obs code (), but later could be shifted to wis.py

[//]: # (TODO: Check efficiency)

 ### nbody

 - Wrapper around REBOUND 
 - See notebooks/Demonstate_Functionality_mpc_nbody.ipynb
 
 ### orbit_cheby
 
 - The main work-horse of the package, based around MSC (Multi Sector Cheby) class. 
 - Contains functionalities to fit MSC to nbody data, evaluate MSCs for ephemeris (see below), calculate RA,Dec from Cartesians, etc etc
 - See notebooks/Demonstrate_Orbital_Chebyshev_Functionality.ipynb [ *I suspect that this notebook is incomplete / has only been shown to explicitly work for the first ~half of the cells*]. 
 - explicit tests do *not* seem to exist, so should be extracted from the demo notebook.
 
 ### precalc & sql 
 
 - The precalc routines perform the necessary precalculations (e.g. nightly healpix locations) and then save the data into the sqlite db
 - The sql module handles all of the reads & writes to/from the sqlite db. 
 - notebooks/Demonstrate_SQLandPreCalc.ipynb 
 
 ### ephem
 
 - Very light wrapper around MSC class. Generates an ephemeris (sky-posn, etc) for given object(s).
 - Function(s) sketched-out but not tested / implemented in any way
 
### mpchecker2
 
- Contains stub functions to do *all* checking-related work (pCheck, MPChecker, CheckID, ...)
- Functions sketched-out but not tested / implemented in any way
- A notebook with stub sections in it exists, notebooks/http://localhost:8867/notebooks/Demonstrate_EndToEnd_Orbit_Checking.ipynb, but with little content

[//]: # (TODO: check this out)

### data_classes

- Classes to hold ad-hoc data (detections, pointing, etc) used by various check functions
- Some tests written in tests/towards_tests_of_data_classes.py

[//]: # (TODO: tests may be under a different module now.)

 ### sifter, sifter_precalc, sifter_sql, sifter_query:

- files from the formerly independent `sifter` repository, now fully merged into here. 
- sifter_precalc parses observations, identifies tracklets and splits observations into seperate tracklets, figures out where the observer is, deals with coordinate frame conversion, calculates the important info about a tracklet (integer night, HealPix, etc) and saves it to a dictionary. 
- sifter_sql creates a database for storing the dictionaries from sifter_precalc
- sifter_query can be used to query the database from sifter_sql such that you can ask "on a given night, what tracklets are in HealPix X?"
- sifter kind of ties it all together (WIP)

# To Do 

- Edit sifter to use the stand-alone obs_pos for figuring out where observatories are, and remove the duplicated code from sifter_precalc.py [MA]

- Develop tests for the many many functionalities within orbit_cheby [MA/MJP]

- Add functionality from "playing_with_convex_hulls_and_ellipsoid_representations" to allow generation of ellipsoid boundary & calculation of overlap of convex hulls. This should probably be added to MSC (boundary) and Detections/Residuals (overlap). [MJP]

- Develop test / demo of the end-to-end flow of data from orbit-fit, through nbody, into sifter_precalc, and then the database. [MA/MJP]

- Develop test / demo of the end-to-end extraction of data from the database, through ephem, into pcheck [MJP]

- Develop tests of functionalities in mpchecker. Start with pCheck (There are many more lower-level functionalities that need to have proper tests written for them, but I suspect that getting the pcheck & ephem codes working and tested will prompt me to get a whole lot of lower-level code sorted too) . [MJP]

[//]: # (Todo: Check status of the above.)

- Modularize cheby_checker and propagate these modules to unit tests, too.
- Run tests as they are (on a Docker container) before making them unit testable.
- Determine SQLite vs. Postgres