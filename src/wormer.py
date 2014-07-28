import numpy as np
from osgeo import gdalnumeric
from osgeo import gdal
from osgeo import gdalconst
from matplotlib import pyplot as plt
import os.path
import FourierDomainGrid as GRID
import FourierDomainOps as FDO
from Utility import isclose, viewRaster
from scipy.ndimage.measurements import label


class Wormer(object):
    """Overview of whole worming method.
    
    Usage example (and Doctest)
    
    >>> foo = Wormer()
    >>> assert isinstance(foo,Wormer)
    >>> assert (foo.base_grid == None)
    >>> grid = np.zeros((512,512),np.float)
    >>> foo.setBaseGrid(grid)
    >>> assert np.allclose(foo.base_grid,grid)
    >>> filename = os.path.abspath('../test_data/SuratBasin/PaddedSuratBasin.ers')
    >>> foo.importGdalRaster(filename)
    >>> assert foo.gdal_input_filename == filename
    >>> assert foo.base_grid.shape == (2048,1536)
    >>> #foo.viewRaster(foo.base_grid)
    >>> viewRaster(foo.wormLevel(dz=-10.*foo.dy))
    
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
        self.ds = gdal.Open(self.gdal_input_filename,gdalconst.GA_ReadOnly) 
        self.geomat = self.ds.GetGeoTransform()
        # FIXME: the following assumes North is up in the image
        self.dx = self.geomat[1]
        self.dy = self.geomat[5]
        
    def exportNewGdalRaster(self, narray, gdal_filename, fmt='ERS'):
        """Export a numpy array into a GDAL image, using the initial image as a prototype.
        """
        gdalnumeric.SaveArray( narray, gdal_filename, format = fmt, prototype = self.ds )
        
    def wormLevel(self,dz):
        # import pdb; pdb.set_trace()
        fdg = GRID.FourierDomainGrid(dx=self.dx, dy=self.dy)
        fdg.setSpatialGrid(self.base_grid)
        fdg.setHatGrid(fdg.simpleFFT(fdg.spatial_grid.astype(np.complex)))
        fdg.buildWavenumbers(fdg.spatial_grid)
        up_fdg = GRID.FourierDomainGrid(dx=self.dx, dy=self.dy)
        fdo = FDO.FourierDomainOps(fdg)
        fdo.buildModK()
        fdo.buildUpwardContinuationOp(dz)
        up_fdg.setHatGrid(fdg.hat_grid * fdo.F_up)
        up_fdg.setSpatialGrid(up_fdg.simpleIFFT(up_fdg.hat_grid))
        bar = fdo.CannyEdgeDetect(up_fdg)
        return bar
    
    def collectWorms(self,level_image, structure=np.ones((3,3))):
        """Takes level_image and labels all of the connected things.
        'connected' is defined by whatever is touching within the
        structuring element -- which defaults to all pixels one
        step away in primary or diagonal directions.
        Returns an array holding the labels, and the number of labels.
        """
        (lbls,n_labels) = label(level_image > -100.,structure=structure)
        return (lbls,n_labels)

        
        
if __name__ == '__main__':
    import doctest
    doctest.testmod()