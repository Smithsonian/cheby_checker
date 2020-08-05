"""
	Tests of the data_classes module
	Currently (20200804) incomplete
	No tests of Detections & Residuals
"""

# Import third-party packages
# --------------------------------------------------------------
import numpy as np
import sys, os
from astropy_healpix import healpy

# Import neighboring packages
# --------------------------------------------------------------
sys.path.append(os.path.dirname(os.path.dirname(
    os.path.realpath(__file__))))
from cheby_checker import data_classes


# -------------------------- Vectorial --------------------------
# NB This is the parent class for Pointings & Detections
# ---------------------------------------------------------------
def test_vectorial():
    # Instantiate (empty) Vectorial object and check type
    V = data_classes.Vectorial()
    assert isinstance(V,data_classes.Vectorial) , f'V is not of the expected type:{type(V)}'
    assert len(V.__slots__) == len( ['iama', 'arr', *data_classes.var_names['Convenience'], *data_classes.var_names['Vectorial']] )

    # Check Vectorial methods
    # (i) Create empty 'arr' of specified size
    V = data_classes.Vectorial()
    assert not hasattr(V, 'arr')
    n = 3
    V._custom_array_init( n )
    assert hasattr(V, 'arr')
    assert V.arr.shape == (n,len(data_classes.var_names['Vectorial']) )

    # (ii) Create 'arr' containing input data (1d input)
    V = data_classes.Vectorial()
    assert not hasattr(V, 'arr')
    inputdata = np.random.random_sample( len(data_classes.var_names['Vectorial']) )
    V._custom_array_init( inputdata )
    assert hasattr(V, 'arr')
    assert V.arr.shape == (1,len(data_classes.var_names['Vectorial']) )
    assert np.all(V.arr[0] == inputdata), f'V.arr[0] = {V.arr[0]}, inputdata = {inputdata}'

    # (iii) Create 'arr' containing input data (2d input)
    V = data_classes.Vectorial()
    assert not hasattr(V, 'arr')
    n = 3
    inputdata = np.random.random_sample( (n, len(data_classes.var_names['Vectorial'])) )
    V._custom_array_init( inputdata )
    assert hasattr(V, 'arr')
    assert V.arr.shape == (n,len(data_classes.var_names['Vectorial']) )
    assert np.all(V.arr == inputdata)

    # (iv) Demonstrate overwritten __getattribute__ ...
    # Here we show that the named attributes are equal to slices of the input array
    for var,num in data_classes.var_maps['Vectorial'].items():
        assert np.all(getattr(V,var) == inputdata[:,num])
    # Here we show that the convenience variables are as expected ...
    for var in data_classes.var_maps['Convenience']:
        if var == 'pos':
            assert np.all(getattr(V,var) == np.array( [V.pos1,V.pos2,V.pos3] ) )
        if var == 'unit_vec':
            assert np.all(getattr(V,var) == np.array( [ np.cos(V.dec)*np.cos(V.ra),
                                                        np.cos(V.dec)*np.sin(V.ra),
                                                        np.sin(V.dec)] ) )
            
        if var == 'HP':
            assert np.all(getattr(V,var) == healpy.vec2pix(V.HP_nside, V.unit_vec[0], V.unit_vec[1], V.unit_vec[2],
                                                            nest=True if V.HP_order=='nested' else False ))
        if var == 'HPlist':
            pass
        
        
    
    
    
# -------------------------- Pointings --------------------------
def test_pointings():
    # (i) Instantiate Pointing object containing input data (1d input) and check
    # NB (1) Has to have data to instantiate
    # NB (2) This also demonstrates that Pointings inherits from Vectorial
    inputdata = np.random.random_sample( len(data_classes.var_names['Pointings']) )
    P = data_classes.Pointings(inputdata)
    assert isinstance(P,data_classes.Vectorial)  , f'P is not of the expected type:{type(P)}'
    assert isinstance(P,data_classes.Pointings)  , f'P is not of the expected type:{type(P)}'
    assert len(P.__slots__) == len( [*data_classes.var_names['Pointings']] )
    assert hasattr(P, 'arr')
    assert P.arr.shape == (1,len(data_classes.var_maps['Pointings']) )
    assert np.all(P.arr[0] == inputdata)


    # (ii) Create Pointings containing input data (2d input)
    n = 3
    inputdata = np.random.random_sample( (n, len(data_classes.var_maps['Pointings'])) )
    P = data_classes.Pointings(inputdata)
    assert hasattr(P, 'arr')
    assert P.arr.shape == (n,len(data_classes.var_maps['Pointings']) )
    assert np.all(P.arr == inputdata)

    # (iv) Demonstrate  __getattribute__ ...
    # Here we show that the named attributes are equal to slices of the input array
    for var,num in data_classes.var_maps['Pointings'].items():
        assert np.all(getattr(P,var) == inputdata[:,num])

    # Here we show that the convenience variables are available (via Vectorial) and have the expected values expected ...
    for var in data_classes.var_maps['Convenience']:
        if var == 'pos':
            assert np.all(getattr(P,var) == np.array( [P.pos1,P.pos2,P.pos3] ) )
        if var == 'unit_vec':
            assert np.all(getattr(P,var) == np.array( [ np.cos(P.dec)*np.cos(P.ra),
                                                        np.cos(P.dec)*np.sin(P.ra),
                                                        np.sin(P.dec)] ) )
            
        if var == 'HP':
            assert np.all(getattr(P,var) == healpy.vec2pix(P.HP_nside, P.unit_vec[0], P.unit_vec[1], P.unit_vec[2],
                                                            nest=True if P.HP_order=='nested' else False ))
        if var == 'HPlist':
            pass
        


# -------------------------- Detections --------------------------
def test_detections():
    pass



# -------------------------- Residuals ---------------------------
def test_residuals():
    pass 
