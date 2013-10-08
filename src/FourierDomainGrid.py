import numpy as np

class FourierDomainGrid(object):
    """Deal with Fourier Domain entities and Operators.
    
    Usage example (and Doctest)
    
    >>> foo = FourierDomainGrid()
    >>> assert isinstance(foo,FourierDomainGrid)
    
    >>> assert (foo.spatial_grid == None)
    
    >>> grid = np.zeros((512,512),dtype=complex)
    >>> foo.setSpatialGrid(grid)
    >>> assert np.allclose(foo.spatial_grid,grid)
    >>> assert foo.spatial_grid.dtype == np.dtype(complex) 
    
    >>> foo.setHatGrid(grid)
    >>> assert np.allclose(foo.hat_grid,grid)
    >>> assert foo.hat_grid.dtype == np.dtype(complex)
    
    # Testing that the FFT of a delta function at the origin is 1 in the transform domain
    # The origin of coordinates from the viewpoint of the FFT is [0,0]
    >>> foo.spatial_grid[0,0] = 1.+0.j
    >>> foo.setHatGrid(foo.simpleFFT(foo.spatial_grid))
    >>> assert np.allclose(foo.hat_grid,(1.+0j))    
    >>> assert np.allclose(foo.simpleIFFT(foo.hat_grid),foo.spatial_grid)
    
    # Testing that the IFT of a delta function at the origin is 1 in the REAL domain
    # The origin of coordinates from the viewpoint of the IFT is [0,0]
    >>> hat_grid = np.zeros((512,512),dtype=complex)
    >>> hat_grid[0,0] = (1.+0.j)
    >>> foo.setSpatialGrid(foo.simpleIFFT(hat_grid))
    >>> assert np.allclose(foo.spatial_grid,(1.+0j)/(512.*512.))    
    >>> assert np.allclose(foo.simpleFFT(foo.spatial_grid),hat_grid)
    
    """

    def __init__(self):
        self.spatial_grid = None
        
    def wavenumbers_x(self,shape ):
        pass
        
    def setSpatialGrid(self,grid):
        """Setter for spatial_grid
        """
        self.spatial_grid = grid
    
    def setHatGrid(self,grid):
        """Setter for hat_grid (wavenumber domain)
        """
        self.hat_grid = grid
    
    def simpleFFT(self,spatial_grid):
        """ Perform a simple FFT without pre-conditioning
            Input: complex; Output: complex
        """        
        return np.fft.fft2(spatial_grid)

    def simpleIFFT(self,hat_grid):
        """ Perform a simple inverse FFT without pre-conditioning
            Input: complex; Output: complex
        """        
        return np.fft.ifft2(hat_grid)
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()    