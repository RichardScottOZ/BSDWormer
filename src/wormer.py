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
#from matplotlib import image as mpimg

class Wormer(object):
    """Overview of whole worming method.
    
    Usage example (and Doctest)
    
    >>> foo = Wormer()
    >>> assert isinstance(foo,Wormer)
    >>> assert (foo.base_grid == None)
    >>> grid = np.zeros((512,512),np.float)
    >>> foo.setBaseGrid(grid)
    >>> assert np.allclose(foo.base_grid,grid)
    >>> filename = os.path.abspath('../test_data/SuartBasin/PaddedSuartBasin.ers')
    >>> foo.importGdalRaster(filename)
    >>> assert foo.gdal_input_filename == filename
    >>> assert foo.base_grid.shape == (2048,1536)
    >>> foo.viewRaster(foo.base_grid)
    
    """
    
    def __init__(self):
        self.base_grid = None
        
    def setBaseGrid(self,grid):
        """Setter for base_grid
        """
        self.base_grid = grid
        
    def importGdalRaster(self,gdal_filename):
        """Import a GDAL raster into a numpy array stashed away as an attribute.
        """
        self.gdal_input_filename = gdal_filename
        self.setBaseGrid(np.array(gdalnumeric.LoadFile(self.gdal_input_filename)))
        
    def viewRaster(self,numpy_grid):
        plt.imshow(numpy_grid)
        plt.show()

    

if __name__ == '__main__':
    import doctest
    doctest.testmod()