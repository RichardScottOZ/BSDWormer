import numpy as np
#import scipy as sp
from osgeo import gdalnumeric
from osgeo import gdal
#from osgeo import ogr
#from osgeo import osr
#from osgeo import gdal_array
from osgeo import gdalconst
from matplotlib import pyplot as plt
import os.path
import sys
import os
#from matplotlib import image as mpimg

class UpwardCont(object):
    """Overview of whole worming method.
    
    Usage example (and Doctest)
    
    >>> foo = UpwardCont()
    >>> assert isinstance(foo,UpwardCont)
    >>> assert (foo.spatial_grid == None)
    >>> grid = np.zeros((512,512),dtype=complex)
    >>> foo.setSpatialGrid(grid)
    >>> assert np.allclose(foo.spatial_grid,grid)
    >>> assert foo.spatial_grid.dtype == np.dtype(complex) 
    >>> foo.setkxkyGrid(grid)
    >>> assert np.allclose(foo.kxky_grid,grid)
    >>> assert foo.kxky_grid.dtype == np.dtype(complex)
    >>> foo.spatial_grid[0,0] = 1
    >>> print(foo.kxky_grid[255,255])  

    """

    def __init__(self):
        self.spatial_grid = None
        
    def setSpatialGrid(self,grid):
        """Setter for spatial_grid
        """
        self.spatial_grid = grid
    
    def setkxkyGrid(self,grid):
        """Setter for kxky_grid (wavenumber domain)
        """
        self.kxky_grid = grid
    
    def simpleFFT(self,spatial_grid):
        """ Perform a simple FFT without pre-conditioning
            Input: complex; Output: complex
        """        
        self.kxky_grid = np.fft.fft2(spatial_grid)

    
if __name__ == '__main__':
    import doctest
    doctest.testmod()    