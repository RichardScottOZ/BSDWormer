import numpy as np

class FourierDomainOps(object):
    """Deal with Fourier Domain entities and Operators.
    
    Usage example (and Doctest)
    
    >>> import FourierDomainGrid as GRID
    >>> bar = GRID.FourierDomainGrid()
    >>> bar.setSpatialGrid(np.ones((512,256),dtype=np.complex))
    >>> bar.setHatGrid(bar.simpleFFT(bar.spatial_grid))
    >>> bar.buildWavenumbers(bar.spatial_grid,dx=1.0,dy=1.0)
    >>> foo = FourierDomainOps(bar)
    >>> assert isinstance(foo,FourierDomainOps)
    >>> foo.buildModK()
    >>> assert foo.modk.shape == foo.fdg.grid_shape
    >>> assert foo.modk[0,0] == 0.0
    >>> assert np.allclose(foo.modk[0,-1],np.abs(foo.fdg.kx[-1]))
    >>> assert np.allclose(foo.modk[-1,0],np.abs(foo.fdg.ky[-1]))
    >>> assert np.allclose(foo.modk[-1,-1],np.sqrt(foo.fdg.ky[-1]**2 + foo.fdg.kx[-1]**2))
    >>> assert np.allclose(foo.modk[5,-6],np.sqrt(foo.fdg.ky[5]**2 + foo.fdg.kx[-6]**2))

    """
    
    def __init__(self, fourier_domain_grid):
        self.fdg = fourier_domain_grid
    
    def buildModK(self):
        kx = self.fdg.kx
        ky = self.fdg.ky
        # ky is 1D as is kx. We are building something of shape
        # len(ky) rows with len(kx) cols by the call to outer.
        # The value being summed are the squares of the values in kx, ky
        # and the two squared values are being summed and inserted into the 
        # result appropriately. Then we do an elementwise sqrt to 
        # form the L2 modulus. 
        self.modk = np.sqrt(np.add.outer(ky*ky,kx*kx))

    
    def UpwardContinuationOp(self,delta_z):
        """Fourier Transform of the Upward Continuation Operator
        """
        self.F_up = np.exp(-delta_z*self.modk)
    
     
if __name__ == '__main__':
    import doctest
    doctest.testmod()         