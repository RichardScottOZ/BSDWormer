import numpy as np

class FourierDomainOps(object):
    """Deal with Fourier Domain entities and Operators.
    
    Usage examples (and Doctests)
    
    >>> import FourierDomainGrid as GRID
    >>> bar = GRID.FourierDomainGrid(dx=1.0,dy=1.0)
    >>> bar.setSpatialGrid(np.ones((512,256),dtype=np.complex))
    >>> bar.setHatGrid(bar.simpleFFT(bar.spatial_grid))
    >>> bar.buildWavenumbers(bar.spatial_grid)
    >>> foo = FourierDomainOps(bar)
    >>> assert isinstance(foo,FourierDomainOps)
    >>> foo.buildModK()
    >>> assert foo.modk.shape == foo.fdg.grid_shape
    >>> assert foo.modk[0,0] == 0.0
    >>> assert np.allclose(foo.modk[0,-1],np.abs(foo.fdg.kx[-1]))
    >>> assert np.allclose(foo.modk[-1,0],np.abs(foo.fdg.ky[-1]))
    >>> assert np.allclose(foo.modk[-1,-1],np.sqrt(foo.fdg.ky[-1]**2 + foo.fdg.kx[-1]**2))
    >>> assert np.allclose(foo.modk[5,-6],np.sqrt(foo.fdg.ky[5]**2 + foo.fdg.kx[-6]**2))
    
    # Test the x gradient for a delta function somewhere on the grid.
    >>> bar = GRID.FourierDomainGrid(dx=1.0,dy=1.0)
    >>> bar.setSpatialGrid(np.zeros((512,256),dtype=np.complex))
    >>> bar.spatial_grid[255,127] = 1.0
    >>> bar.setHatGrid(bar.simpleFFT(bar.spatial_grid))
    >>> bar.buildWavenumbers(bar.spatial_grid)
    >>> foo = FourierDomainOps(bar)
    >>> foo.buildDxOp()
    >>> dxBar = GRID.FourierDomainGrid(dx=1.0,dy=1.0)
    >>> dxBar.setHatGrid(bar.hat_grid*foo.F_dxOp)
    >>> dxBar.setSpatialGrid(dxBar.simpleIFFT(dxBar.hat_grid))
    >>> assert np.allclose(dxBar.spatial_grid[255,125],0.0)
    >>> assert np.allclose(dxBar.spatial_grid[255,126],1.0)
    >>> assert np.allclose(dxBar.spatial_grid[255,127],0.0)
    >>> assert np.allclose(dxBar.spatial_grid[255,128],-1.0)
    >>> assert np.allclose(dxBar.spatial_grid[255,129],0.0)
    >>> assert np.allclose(dxBar.spatial_grid[254,:],0.0)
    >>> assert np.allclose(dxBar.spatial_grid[256,:],0.0)

    # Test the y gradient for a delta function somewhere on the grid.
    >>> bar = GRID.FourierDomainGrid(dx=1.0,dy=1.0)
    >>> bar.setSpatialGrid(np.zeros((512,256),dtype=np.complex))
    >>> bar.spatial_grid[255,127] = 1.0
    >>> bar.setHatGrid(bar.simpleFFT(bar.spatial_grid))
    >>> bar.buildWavenumbers(bar.spatial_grid)
    >>> foo = FourierDomainOps(bar)
    >>> foo.buildDyOp()
    >>> dyBar = GRID.FourierDomainGrid(dx=1.0,dy=1.0)
    >>> dyBar.setHatGrid(bar.hat_grid*foo.F_dyOp)
    >>> dyBar.setSpatialGrid(dyBar.simpleIFFT(dyBar.hat_grid))
    >>> assert np.allclose(dyBar.spatial_grid[253,127],0.0)
    >>> assert np.allclose(dyBar.spatial_grid[254,127],1.0)
    >>> assert np.allclose(dyBar.spatial_grid[255,127],0.0)
    >>> assert np.allclose(dyBar.spatial_grid[256,127],-1.0)
    >>> assert np.allclose(dyBar.spatial_grid[257,127],0.0)
    >>> assert np.allclose(dyBar.spatial_grid[:,126],0.0)
    >>> assert np.allclose(dyBar.spatial_grid[:,128],0.0)
    
    # Test for upward continuation
    >>> bar = GRID.FourierDomainGrid(dx=1.0,dy=1.0)
    >>> bar.setSpatialGrid(np.zeros((512,256),dtype=np.complex))
    >>> bar.spatial_grid[255,127] = 1.0
    >>> bar.setHatGrid(bar.simpleFFT(bar.spatial_grid))
    >>> bar.buildWavenumbers(bar.spatial_grid)
    >>> foo = FourierDomainOps(bar)
    >>> foo.buildModK()
    >>> foo.buildUpwardContinuationOp(20.0)
    >>> upcontBar = GRID.FourierDomainGrid(dx=1.0,dy=1.0)
    >>> upcontBar.setHatGrid(bar.hat_grid * foo.F_up)
    >>> upcontBar.setSpatialGrid(upcontBar.simpleIFFT(upcontBar.hat_grid))
    >>> assert np.allclose(np.sum(upcontBar.spatial_grid),1.0)
    >>> foo.buildUpwardContinuationOp(10.0)
    >>> upcontBar = GRID.FourierDomainGrid(dx=1.0,dy=1.0)
    >>> upcontBar.setHatGrid(bar.hat_grid * foo.F_up)
    >>> upcontBar.setSpatialGrid(upcontBar.simpleIFFT(upcontBar.hat_grid))
    >>> assert np.allclose(np.sum(upcontBar.spatial_grid),1.0)
    >>> foo.buildUpwardContinuationOp(5.0)
    >>> upcontBar = GRID.FourierDomainGrid(dx=1.0,dy=1.0)
    >>> upcontBar.setHatGrid(bar.hat_grid * foo.F_up)
    >>> upcontBar.setSpatialGrid(upcontBar.simpleIFFT(upcontBar.hat_grid))
    >>> assert np.allclose(np.sum(upcontBar.spatial_grid),1.0)
    >>> foo.buildUpwardContinuationOp(1.0)
    >>> upcontBar = GRID.FourierDomainGrid(dx=1.0,dy=1.0)
    >>> upcontBar.setHatGrid(bar.hat_grid * foo.F_up)
    >>> upcontBar.setSpatialGrid(upcontBar.simpleIFFT(upcontBar.hat_grid))
    >>> assert np.allclose(np.sum(upcontBar.spatial_grid),1.0)
    

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

    
    def buildUpwardContinuationOp(self,delta_z):
        """Fourier Transform of the Upward Continuation Operator
        
        FIXME: Find an analytic expression for some upward continuation, 
        and write a test that checks it.
        
        The tests in the doctest above only exercise this loosely.
        We are passing a bunch of necessary tests. But these tests are not 
        sufficient tests.
        """
        self.F_up = np.exp(-delta_z*self.modk)
        
    def buildDxOp(self):
        """The analytic expression for a derivative is available in closed form in the 
        Fourier Domain. Unfortunately, it exponentially amplifies the shortest 
        wavelengths in the signal, leading to extreme amplification of the 
        content that usually isn't known very well at all. The resulting 
        derivatives are 'noisy' in the extreme.
        As is well known, another approach is to estimate gradients as central 
        differences in the Fourier domain on the
        discrete grid (or, perhaps using higher order finite difference stencils 
        for increased accuracy).
        In this routine, we are implementing a central difference in x 
        explicitly using the Fourier shift theorem. Denoting the transform pair
        f <=> F, the shift theorem states:
          F{f(x-x_0)}(kx) = exp(-2 pi j kx x_0) F(kx)
        
        We write the central difference on a grid spacing of \Delta x as:
        df/dx \approx [ f( x + \Delta x) - f( x - \Delta x) ]/(2 \Delta x)
        
        Employing the shift theorem, in the Fourier domain this becomes:
        F[df/dx](kx) \approx 
          F(kx) [ exp(2 pi j kx \Delta x) - exp(-2 pi j kx \Delta x) ] / (2 \Delta x)
        Factoring out the F(kx), the operator becomes the rest of the expression.

        Now, that expression has the functional form of a sine. Writing 
        y = 2 pi kx \Delta x, From Euler's formula we find:
        
        sin(y) = (exp(j y) - exp(-j y))/2j
        so the derivative operator reduces to:
        2j/(\Delta x) sin(2 pi kx \Delta x)
        This expression has the major benefit of not amplifying the short 
        wavelengths.
        """
        kx = self.fdg.kx
        dx = self.fdg.dx
        self.F_dxOp = ((2.*(0.+1.j))/dx)*np.sin(2.*np.pi*kx[np.newaxis,:]*dx)
    
    def buildDyOp(self):
        """
        See the comment for the DxOp immediately above. Everything is the same 
        except that we substitute y for x
        """
        ky = self.fdg.ky
        dy = self.fdg.dy
        self.F_dyOp = ((2.*(0.+1.j))/dy)*np.sin(2.*np.pi*ky[:,np.newaxis]*dy)
    
     
if __name__ == '__main__':
    import doctest
    doctest.testmod()         