import sys
import os

# Globally add the reboundx directory to the path, because it can't be pip installed.
sys.path.append(os.environ['REBX_DIR'])
