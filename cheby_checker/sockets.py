# -*- coding: utf-8 -*-
# cheby_checker/cheby_checker/ephem

"""
--------------------------------------------------------------
cheby_checker's / sockets module.

Aug 2020
Matt Payne

This module is intended to handle the server-client connections
for cheby-checker

It is currently developmental while MJP gets his head around the
sockets methods

It is intended to set up a server (e.g. on marsden) that will
listen for requests for data from clients (either on marsden,
and/or on some other machines). It is expected that this will
be requests for ephemeris / mpchecker / pcheck / ...

--------------------------------------------------------------
"""


# Import third-party packages
# --------------------------------------------------------------
import threading
import socket
import pickle
import numpy as np
import struct


# Import neighboring packages
# --------------------------------------------------------------
#try:
#    import orbit_cheby
#except ImportError:
#    from . import orbit_cheby


# Socket-Server-Related Object Definitions
# --------------------------------------------------------------
class SharedSocket(object):
    """
    Primarily used to provide shared methods
    to send & receive messages of arbitrary length/size
    https://stackoverflow.com/questions/17667903/python-socket-receive-large-amount-of-data
    """

    default_server_host = '' # '127.0.0.1'
    default_server_port = 54321

    def __init__(self,):
        pass

    def send_msg(self, sock, msg):
        # Prefix each message with a 4-byte length (network byte order)
        msg = struct.pack('>I', len(msg)) + msg
        sock.sendall(msg)

    def recv_msg(self, sock):
        # Read message length and unpack it into an integer
        raw_msglen = self.recvall(sock, 4)
        if not raw_msglen:
            return None
        msglen = struct.unpack('>I', raw_msglen)[0]
        # Read the message data
        return self.recvall(sock, msglen)

    def recvall(self, sock, n):
        # Helper function to recv n bytes or return None if EOF is hit
        data = bytearray()
        while len(data) < n:
            packet = sock.recv(n - len(data))
            if not packet:
                return None
            data.extend(packet)
        return data


class Server(SharedSocket):
    """
    Set up a server that will listen for clients

    usage:
    ------
    Server().listen()


    """

    def __init__(self, host=None, port=None):
        
        self.host = host if host is not None else self.default_server_host
        self.port = port if port is not None else self.default_server_port
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        #  associate the socket with a specific network interface and port number
        self.sock.bind((self.host, self.port))


    def _demo_listen(self, startup_func = None ):
        """
        Demo function to illustrate how to set-up server
        and to allow tests of various functionalities
        """
    
        # Allow for doing a lengthy job on start-up
        # E.g. loading a big ephemeris file
        if startup_func:
            self._a_function_that_takes_a_few_seconds_to_evaluate(i=5)
            
        # listen() enables a server to accept() connections
        self.sock.listen(5)
        print('\nI are listen')
        while True :
            
            # accept() blocks and waits for an incoming connection.
            # One thing that’s imperative to understand is that we now have a
            # new socket object from accept(). This is important since it’s the
            # socket that you’ll use to communicate with the client. It’s distinct
            # from the listening socket that the server is using to accept new
            # connections:
            client, address = self.sock.accept()
            client.settimeout(60)
            
            # Either of the below work ...
            #self.listenToClient(client,address)
            threading.Thread(target = self._listenToClient,args = (client,address)).start()

    def _listenToClient(self, client, address):
        while True:
            try:
                data        = self.recv_msg(client)
                if data:
                
                    # Could insert some function here to do
                    # useful work: E.g. call cheby_checker
                    
                    # Set the response to echo back the data
                    # but with the addition of a recieved boolean
                    data_string = self._add_received_to_data(data)
                    self.send_msg(client, data_string)
                    
                else:
                    raise error('Client disconnected')
            except:
                client.close()
                return False
                
    def _add_received_to_data(self, data):
        data_object             = pickle.loads(data)
        data_object["received"] = True
        return                  pickle.dumps(data_object)
        
    def _a_function_that_takes_a_few_seconds_to_evaluate(self, i=1):
        print("\nDoing a lengthy calculation")
        for _ in range(i):
            n           = int(1e4)
            expected    = (n,n)
            actual      = np.random.random_sample((n,n)).shape
        return True

class Client(SharedSocket):
    """
    Convenience class & method(s) for connecting to server
    """

    def __init__(self, host=None, port=None):
        self.server_host = host if host is not None else self.default_server_host
        self.server_port = port if port is not None else self.default_server_port
            
    def _demo_client(self, data_dict = {'msg' : 'Hello World!' } , VERBOSE = False ):

        #try:
            # Pickle the provided data dictionary
            data_string = pickle.dumps(data_dict)
            if VERBOSE:
                print('message sent (before pickling)     :',data_dict)
                
            # Create a socket objects
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:

                s.settimeout(111)
                # Connect to the server
                s.connect((self.server_host, self.server_port))
                # Send data to the server
                self.send_msg(s, data_string)
                # Read the reply from the server
                data_string = self.recv_msg(s)
                data_dict   = pickle.loads(data_string)
                if VERBOSE:
                    print('message received (after unpickling):',data_dict)
        #except Exception as e:
            #data_dict = {'ERROR': e}
            
            return data_dict


if __name__ == "__main__":
    Server().listen(startup_func = True )
