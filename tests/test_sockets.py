# Import third-party packages
# --------------------------------------------------------------
import numpy as np
import threading
import socket

# Import neighboring packages
# --------------------------------------------------------------
from cheby_checker import sockets

# Some tests of general sockets methods
# (tests of cheby-specific sockets classes are below)
# --------------------------------------------------------------
def run_fake_server(HOST = '127.0.0.1', PORT = 65432):
    # Run a server to listen for a connection and then close it
    # NB Don't run as-is : put in a thread, or something, to background-it
    server_sock = socket.socket()
    server_sock.bind((HOST,PORT))
    server_sock.listen(0)
    server_sock.accept()
    server_sock.close()
    assert True


def run_fake_client(HOST = '127.0.0.1', PORT = 65432):
    # This is our fake test client that is just going to attempt to connect and disconnect
    fake_client = socket.socket()
    fake_client.settimeout(2)
    fake_client.connect((HOST,PORT))
    fake_client.close()


def test_threading_server_and_client(HOST = '127.0.0.1', PORT = 65431):
    """ Adapted from https://www.devdungeon.com/content/unit-testing-tcp-server-client-python """

    # Start fake server in background thread
    server_thread = threading.Thread(target=run_fake_server, args=(HOST,PORT))
    server_thread.start()

    # Test the clients basic connection and disconnection
    # *** If the above server is not running, then this will not connect ***
    run_fake_client(HOST=HOST, PORT=PORT)

    # Ensure server thread ends
    server_thread.join()


# Tests of cheby-specific sockets classes
# --------------------------------------------------------------
#def test_server_instantiation():
#    S = sockets.Server()
#    assert isinstance(S,sockets.Server)
#    C = sockets.Client()
#    assert isinstance(C,sockets.Client)
#    return True
def test_demo_client_server_connect():
    
    # launch client
    C = sockets.Client()
    
    # Send default "Hello World" message
    # & collect returned signal from server
    received = C._demo_client(VERBOSE=True)

    # Check content of demo message
    # is as expected (hard-coded demo content)
    assert isinstance(received, dict)
    
    ERR_STR = "\n".join( [ "*** received = %r " %  received,
    "This is likely to be caused if/when a server is NOT RUNNING.",
    "A server needs to be launched BEFORE running these tests.",
    "It is likely that I could do so as part of the test process.",
    "But I haven't worked out how to do that yet."
    " (https://realpython.com/testing-third-party-apis-with-mock-servers/) "
    "For now, the simplest way to get a demo server running is to execute the following python command",
    "python sockets_server_starter_DO_NOT_DELETE.py &"
    "Then the pytests can be run.",
    "Try to remember to kill the server afterwards ... "
    ] )
    assert "msg" in received, ERR_STR
    
    assert received["msg"] == 'Hello World!'
    assert "received" in received
    assert received["received"] == True


def test_demo_big_message_exchange():
    
    # launch client
    C = sockets.Client()
    
    # Send large message
    # & collect returned signal from server
    n = int(3e3)
    received = C._demo_client(data_dict = { "np" : np.random.random_sample( (n,n) ) }, VERBOSE=True)
    
    # check that all is as expected
    assert isinstance(received, dict)
    assert "received" in received
    assert received["received"] == True
    assert "np" in received
    assert isinstance( received["np"], np.ndarray )
    assert received["np"].shape == (n,n)
