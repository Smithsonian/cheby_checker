# Might be useful for "test_sockets" to be able to remote-launch a server
import sys,os

# Import neighboring packages
# --------------------------------------------------------------
sys.path.append(os.path.dirname(os.path.dirname(
    os.path.realpath(__file__))))
from cheby_checker import sockets

# Run socket server 
# --------------------------------------------------------------
S = sockets.Server()
S._demo_listen(startup_func=True)
