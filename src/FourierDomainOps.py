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

class UpwardCont(object):
    """Overview of whole worming method.
    
    Usage example (and Doctest)
    
    >>> foo = UpwardCont()
    >>> assert isinstance(foo,UpwardCont)
    >>> assert (foo.kx_ky_grid == None)

    """

    def __init__(self):
        self.kx_ky_grid = None



if __name__ == '__main__':
    import doctest
    doctest.testmod()    