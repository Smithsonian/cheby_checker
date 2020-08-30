# Imports / set-up
import numpy as np ;
import os,sys ;
sys.path.append(os.environ['REBX_DIR']) ;
from examples.ephem_forces.ephem_forces import integration_function ;

def test_malloc():
  ''' early version had malloc problem : check for that '''

  tstart=2456117.641933589 ;
  tstep=20 ;
  trange=1000 ;
  geocentric=False ;
  n_particles=1 ;
  reparsed_input=np.array([-2.0938349524664743,1.0009137200092553,0.41979849545335507,-0.004226738336365523, -0.009129140909705197, -0.0036271214539287102])

 
  # Call repeatedly to check for random crash with malloc ...
  for i in range(1000):
    integration_function(tstart, tstep, trange, geocentric,n_particles, reparsed_input)
