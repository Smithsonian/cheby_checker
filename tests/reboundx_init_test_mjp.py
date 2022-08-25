'''
MJP : 2021-07-16

Copying some simple run commands from MJH's "integrator_example.ipynb" notebook

I just want to have a very simple test/demo that I can use to see whether a basic installation
has worked or not ...

'''

# Do basic required imports
import numpy as np
import time 
import sys, os 
sys.path.append(os.environ['REBX_DIR'])
from examples.ephem_forces import ephem_forces


instates = np.array([[3.338875349745594E+00, -9.176518281675284E-01, -5.038590682977396E-01, 2.805663319000732E-03, 7.550408687780768E-03, 2.980028206579994E-03]])
n_particles = 1

tstart, trange = 2458849.5, 2000
epoch = tstart
tend = tstart + trange
print('tstart,tend,epoch,instates',tstart,tend,epoch,instates)

# Do the run(s) ...
start_time = time.time()
times, states, var, var_ng, status = ephem_forces.production_integration_function_wrapper(tstart, tend, epoch, instates)
print('after 1st run ... ellapsed_time = ', time.time() - start_time )
times, states, var, var_ng, status = ephem_forces.production_integration_function_wrapper(tstart, tend, epoch, instates)
print('after 2nd run ...ellapsed_time = ', time.time() - start_time)
times, states, var, var_ng, status = ephem_forces.production_integration_function_wrapper(tstart, tend, epoch, instates)
print('after 3rd and final run ...ellapsed_time = ', time.time() - start_time)

