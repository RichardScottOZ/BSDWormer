import numpy as np
import FourierDomainGrid as GRID
from Utility import isclose

class FourierDomainOps(object):
    """Deal with Fourier Domain entities and Operators.
    
    Usage examples (and Doctests)
    
    >>> bar = GRID.FourierDomainGrid(dx=1.0,dy=1.0)
    >>> bar.setSpatialGrid(np.ones((512,256),dtype=np.complex))
    >>> bar.setHatGrid(bar.simpleFFT(bar.spatial_grid))
    >>> bar.buildWavenumbers(bar.spatial_grid)
    >>> foo = FourierDomainOps(bar)
    >>> assert isinstance(foo,FourierDomainOps)
    >>> assert np.allclose(bar.spatial_grid,foo.fdg.spatial_grid)
    >>> assert foo.F_dxOp is None
    >>> assert foo.F_dyOp is None
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
    
    # Test buildGradVector
    >>> bar = GRID.FourierDomainGrid(dx=1.0,dy=1.0)
    >>> bar.setSpatialGrid(np.zeros((512,256),dtype=np.complex))
    >>> bar.spatial_grid[255,127] = 1.0
    >>> bar.setHatGrid(bar.simpleFFT(bar.spatial_grid))
    >>> foo = FourierDomainOps(bar) #We have hat_grid already
    >>> (x_vect,y_vect) = foo.buildGradVector(bar)
    >>> assert isinstance(x_vect, GRID.FourierDomainGrid)
    >>> assert isinstance(y_vect, GRID.FourierDomainGrid)
    >>> assert x_vect.hat_grid.shape == (512,256)
    >>> assert y_vect.hat_grid.shape == (512,256)
    >>> assert x_vect.spatial_grid.shape == (512,256)
    >>> assert y_vect.spatial_grid.shape == (512,256)
    >>> del foo
    >>> del bar
    
    >>> bar = GRID.FourierDomainGrid(dx=1.0,dy=1.0)
    >>> bar.setSpatialGrid(np.zeros((512,256),dtype=np.complex))
    >>> bar.spatial_grid[255,127] = 1.0
    >>> foo = FourierDomainOps(bar) 
    >>> (x_vect,y_vect) = foo.buildGradVector(bar)
    >>> assert x_vect.hat_grid.shape == (512,256)
    >>> assert y_vect.hat_grid.shape == (512,256)
    >>> assert x_vect.spatial_grid.shape == (512,256)
    >>> assert y_vect.spatial_grid.shape == (512,256)    
    >>> assert np.allclose(x_vect.spatial_grid[255,125],0.0)
    >>> assert np.allclose(x_vect.spatial_grid[255,126],1.0)
    >>> assert np.allclose(x_vect.spatial_grid[255,127],0.0)
    >>> assert np.allclose(x_vect.spatial_grid[255,128],-1.0)
    >>> assert np.allclose(x_vect.spatial_grid[255,129],0.0)
    >>> assert np.allclose(x_vect.spatial_grid[254,:],0.0)
    >>> assert np.allclose(x_vect.spatial_grid[256,:],0.0)
    >>> assert np.allclose(y_vect.spatial_grid[253,127],0.0)
    >>> assert np.allclose(y_vect.spatial_grid[254,127],1.0)
    >>> assert np.allclose(y_vect.spatial_grid[255,127],0.0)
    >>> assert np.allclose(y_vect.spatial_grid[256,127],-1.0)
    >>> assert np.allclose(y_vect.spatial_grid[257,127],0.0)
    >>> assert np.allclose(y_vect.spatial_grid[:,126],0.0)
    >>> assert np.allclose(y_vect.spatial_grid[:,128],0.0)
    
    # Test buildMod2DVect
    >>> mod_grad = foo.buildMod2DVect(x_vect,y_vect)
    >>> assert np.allclose(1.0, mod_grad[255,126])
    >>> assert np.allclose(1.0, mod_grad[255,128])
    >>> assert np.allclose(1.0, mod_grad[254,127])
    >>> assert np.allclose(1.0, mod_grad[256,127])
    >>> assert np.allclose(0.0, mod_grad[255,127])
    
    # Test buildNormed2DVect
    >>> (unit_x,unit_y) = foo.buildUnit2DVect(x_vect,y_vect,mod_grad)
    >>> assert np.allclose(unit_x.spatial_grid[255,125],0.0)
    >>> assert np.allclose(unit_x.spatial_grid[255,126],1.0)
    >>> assert np.allclose(unit_x.spatial_grid[255,127],0.0)
    >>> assert np.allclose(unit_x.spatial_grid[255,128],-1.0)
    >>> assert np.allclose(unit_x.spatial_grid[255,129],0.0)
    >>> assert np.allclose(unit_x.spatial_grid[254,:],0.0)
    >>> assert np.allclose(unit_x.spatial_grid[256,:],0.0)
    >>> assert np.allclose(unit_y.spatial_grid[253,127],0.0)
    >>> assert np.allclose(unit_y.spatial_grid[254,127],1.0)
    >>> assert np.allclose(unit_y.spatial_grid[255,127],0.0)
    >>> assert np.allclose(unit_y.spatial_grid[256,127],-1.0)
    >>> assert np.allclose(unit_y.spatial_grid[257,127],0.0)
    >>> assert np.allclose(unit_y.spatial_grid[:,126],0.0)
    >>> assert np.allclose(unit_y.spatial_grid[:,128],0.0)

    >>> (unit_x,unit_y) = foo.buildUnit2DVect(x_vect,y_vect)
    >>> assert np.allclose(unit_x.spatial_grid[255,125],0.0)
    >>> assert np.allclose(unit_x.spatial_grid[255,126],1.0)
    >>> assert np.allclose(unit_x.spatial_grid[255,127],0.0)
    >>> assert np.allclose(unit_x.spatial_grid[255,128],-1.0)
    >>> assert np.allclose(unit_x.spatial_grid[255,129],0.0)
    >>> assert np.allclose(unit_x.spatial_grid[254,:],0.0)
    >>> assert np.allclose(unit_x.spatial_grid[256,:],0.0)
    >>> assert np.allclose(unit_y.spatial_grid[253,127],0.0)
    >>> assert np.allclose(unit_y.spatial_grid[254,127],1.0)
    >>> assert np.allclose(unit_y.spatial_grid[255,127],0.0)
    >>> assert np.allclose(unit_y.spatial_grid[256,127],-1.0)
    >>> assert np.allclose(unit_y.spatial_grid[257,127],0.0)
    >>> assert np.allclose(unit_y.spatial_grid[:,126],0.0)
    >>> assert np.allclose(unit_y.spatial_grid[:,128],0.0)

    """
    
    def __init__(self, fourier_domain_grid):
        self.fdg = fourier_domain_grid
        self.F_dxOp = None
        self.F_dyOp = None
    
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
        self.F_up = np.exp(-2.*np.pi*delta_z*self.modk)
        
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
        
    def buildGradVector(self,fdg):
        """Builds the 2D gradient vector of a FDG.
        """
        try:
            assert isinstance(np.ndarray,fdg.hat_grid)
        except:
            fdg.setHatGrid(fdg.simpleFFT(fdg.spatial_grid))
            
        try:
            assert isinstance(np.ndarray,fdg.kx)
            assert isinstance(np.ndarray,fdg,ky)
        except:
            fdg.buildWavenumbers(fdg.hat_grid)   

        if self.F_dxOp == None:
            self.buildDxOp()
        if self.F_dyOp == None:
            self.buildDyOp()
        
        dx_vect = GRID.FourierDomainGrid(dx=fdg.dx,dy=fdg.dy)
        dx_vect.setHatGrid(self.F_dxOp*fdg.hat_grid)
        dx_vect.setSpatialGrid(dx_vect.simpleIFFT(dx_vect.hat_grid))
        
        dy_vect = GRID.FourierDomainGrid(dx=fdg.dx,dy=fdg.dy)
        dy_vect.setHatGrid(self.F_dyOp*fdg.hat_grid)
        dy_vect.setSpatialGrid(dy_vect.simpleIFFT(dy_vect.hat_grid))
        return (dx_vect,dy_vect)
    
    def buildMod2DVect(self,fdg_x,fdg_y):
        """This code has a numerical flaw in that it will return
        zero modulus in the diagonal FD cells (pixels) for a delta
        function input. This is because the gradient operators are 
        built for 1-D gradients only, and have zero entries in those
        diagonal pixels. 
        We will barge on, but be aware that some class of subtle numerical
        gradient artifacts are introduced via our strategy.
        Maybe we can fix this in future versions of the code???
        """
        return np.sqrt(fdg_x.spatial_grid**2 + fdg_y.spatial_grid**2)

    def buildUnit2DVect(self,fdg_x,fdg_y,mod=None):
        """Deal with normalizing vectors with possibly zero lengths.
        """
        if mod == None:
            mod = self.buildMod2DVect(fdg_x,fdg_y)
        mask = isclose(0.,mod).astype(int)
        normed_x = GRID.FourierDomainGrid(dx=fdg_x.dx,dy=fdg_x.dy)
        normed_x.setSpatialGrid(np.choose(mask,[fdg_x.spatial_grid/mod,0.0]))
        normed_y = GRID.FourierDomainGrid(dx=fdg_y.dx,dy=fdg_y.dy)
        normed_y.setSpatialGrid(np.choose(mask,[fdg_y.spatial_grid/mod,0.0]))
        return (normed_x,normed_y)
        
    def CannyEdgeDetect(self,fdg):
        """A 2D Canny edge detector using upward continuation
        as a blurring operation, instead of a Gaussian (which would
        be appropriate for a Diffusion equation problem).
        Implementing: scalar-product(grad(norm(grad(f))),grad(f)/norm(grad(f))) = 0 
        """
        (grad_x,grad_y) = self.buildGradVector(fdg) # vector
        norm_grad = self.buildMod2DVect(grad_x,grad_y)   # scalar
        unit_x,unit_y = self.buildUnit2DVect(grad_x,grad_y)
        (grad_of_norm_x,grad_of_norm_y) = self.buildGradVector(norm_grad)
        inner_product = grad_of_norm_x*unit_x + grad_of_norm_y*unit_y

if __name__ == '__main__':
    import doctest
    doctest.testmod()         