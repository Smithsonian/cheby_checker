# cheby_checker tests go in this directory

## 2022 : MJP: The following have been checked to pass tests 

### test_convenience_Horizons.py
 - convenience_Horizons.py only depends on EXTERNAL packages
 - pytest test_convenience_Horizons.py
 - passes within containerized environment 

### test_obs_pos.py
 - obs_pos.py depends on INTERNAL MPC_library.py: eventually want to change dependancy to wis.py
 - pytest test_obs_pos.py
 - passes within containerized environment

### test_cheby_checker.py
 - cheby_checker.py only depends on EXTERNAL packages
 - pytest test_cheby_checker.py
 - passes within containerized environment

### test_sql.py :
- sql.py has an INTERNAL dependance on the cheby_checker.py module (above)
- pytest test_sql.py
- passes within containerized environment

### test_coco.py :
- coco.py has an INTERNAL dependance on the cheby_checker.py module (above)
- coco.py also depends on INTERNAL MPC_library.py: should shift the rotn & ecliptic angle variables to coco & Base respectively.
- pytest test_coco.py
*** PROBLEM WITH TESTS IN CONTAINER!!!***
   ---  The version on laptop gives different answers depending on whether the center is '399' or 'Geocenter'
   >>> Horizons('399', '500@0', epochs=2458937.000000000, id_type='majorbody').vectors(refplane='earth')
            ['x',
            -0.9931022752569104 
   >>> Horizons('Geocenter', '500@0', epochs=2458937.000000000, id_type='majorbody').vectors(refplane='earth')
            ['x',
            -0.9931015634277218 

   --- The version in the container gives the same answer for both
   >>> Horizons('399', '500@0', epochs=2458937.000000000, id_type='majorbody').vectors(refplane='earth')
            'x', ...
            -0.9931022752569104 ...

   --- The version on the website ...
        $$SOE
        2458937.000000000 = A.D. 2020-Mar-28 12:00:00.0000 TDB [del_T=     69.185671 s]
        X =-9.931022752569104E-01 Y =-1.208197618935286E-01 Z =-5.232300575907833E-02
        ...
        $$EOE

