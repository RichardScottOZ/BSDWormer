import numpy as np
from matplotlib import pyplot as plt
import sys
sys.path += ['/Users/frank/Documents/Src/Git Stuff/BSDWormer/src/Py3Vtk/pyvtk/build/lib/pyvtk']
import pyvtk as PV


def viewRaster(numpy_grid):
    plt.imshow(numpy_grid)
    plt.show()

def writeVtkImage(filename,image,origin,spacing):
    image_points = [val for val in np.ravel(image)]
    image_coords = [image.shape[1],image.shape[0],1]
    vtk = PV.VtkData(PV.StructuredPoints(image_coords,
                                         origin=origin,
                                         spacing=spacing),
                     'Image',
                     PV.PointData(PV.Scalars(image_points,name='field vals')),
                     'Image Data'
                    )
    vtk.tofile(filename+"_image.vtk")

def isclose(a, b, rtol=1.e-5, atol=1.e-8, check_invalid=True):
    """Similar to numpy.allclose, but returns a boolean array.
    See numpy.allclose for an explanation of *rtol* and *atol*.
    
    # A few quick tests...
    >>> assert np.any(isclose(0.300001, np.array([0.1, 0.2, 0.3, 0.4])))

    >>> x = np.array([0.1, np.nan, np.inf, -np.inf])
    >>> y = np.array([0.1000001, np.nan, np.inf, -np.inf])
    >>> assert np.all(isclose(x, y))

    >>> x = np.array([0.1, 0.2, np.inf])
    >>> y = np.array([0.101, np.nan, 0.2])
    >>> assert not np.all(isclose(x, y))
    """

    def within_tol(x, y, atol, rtol):
        return np.less_equal(np.abs(x-y), atol + rtol * np.abs(y))
    x = np.array(a, copy=False)
    y = np.array(b, copy=False)
    if not check_invalid:
        return within_tol(x, y, atol, rtol)
    xfin = np.isfinite(x)
    yfin = np.isfinite(y)
    if np.all(xfin) and np.all(yfin):
        return within_tol(x, y, atol, rtol)
    else:
        # Avoid subtraction with infinite/nan values...
        cond = np.zeros(np.broadcast(x, y).shape, dtype=np.bool)
        mask = xfin & yfin
        cond[mask] = within_tol(x[mask], y[mask], atol, rtol)
        # Inf and -Inf equality...
        cond[~mask] = (x[~mask] == y[~mask])
        # NaN equality...
        cond[np.isnan(x) & np.isnan(y)] = True
        return cond

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    
