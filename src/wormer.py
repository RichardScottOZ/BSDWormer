import numpy as np
import scipy as sp

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
        """Setter for base_grid
        >>> foo = Wormer()
        >>> filename = '../test_data/SuartBasin/suart_basin.ers'
        >>> foo.importGdalRaster(filename)
        >>> assert foo.gdal_input_filename == filename
        """
        self.gdal_input_filename = gdal_filename

    
if __name__ == "__main__":
    import doctest
    doctest.testmod()