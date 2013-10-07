import numpy as np
import scipy as sp
from osgeo import gdalnumeric
from osgeo import gdal
#from osgeo import ogr
#from osgeo import osr
#from osgeo import gdal_array
from osgeo import gdalconst

class Wormer(object):
    """Class to be overview of whole worming method
    
    >>> foo = Wormer()
    >>> assert isinstance(foo,Wormer)
    >>> assert (foo.base_grid == None)
    
    """
    
    def __init__(self):
        self.base_grid = None
        
    def setBaseGrid(self,grid):
        """Setter for base_grid
        >>> foo = Wormer()
        >>> grid = np.zeros((512,512),np.float)
        >>> foo.setBaseGrid(grid)
        >>> assert np.allclose(foo.base_grid,grid)
        """
        self.base_grid= grid
        
    def importGdalRaster(self,gdal_filename):
        """Import a GDAL raster into a numpy array
        >>> foo = Wormer()
        >>> filename = '../test_data/SuartBasin/suart_basin.ers'
        >>> foo.importGdalRaster(filename)
        >>> assert foo.gdal_input_filename == filename
        >>> assert foo.base_grid.shape == (1440,960)
        """
        self.gdal_input_filename = gdal_filename
        self.setBaseGrid(np.array(gdalnumeric.LoadFile(self.gdal_input_filename))

    
if __name__ == "__main__":
    import doctest
    doctest.testmod()