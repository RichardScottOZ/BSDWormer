import numpy as np
from matplotlib import pyplot as plt
import sys
sys.path += ['/Users/frank/Documents/Src/Git Stuff/BSDWormer/src/Py3Vtk/pyvtk/build/lib/pyvtk']
import pyvtk as PV
from osgeo import gdal,ogr,osr



def viewRaster(numpy_grid):
    plt.imshow(numpy_grid)
    plt.show()
    
def writeGDALRasterFromNumpyArray(dst_filename,array,geotransform,proj):
    driver = gdal.GetDriverByName('GTiff')
    y_pixels,x_pixels = array.shape

    dataset = driver.Create(
        dst_filename,
        x_pixels,
        y_pixels,
        1,
        gdal.GDT_Float32, )

    dataset.SetGeoTransform(geotransform)  
    dataset.SetProjection(proj)
    dataset.GetRasterBand(1).WriteArray(array)
    dataset.FlushCache()  # Write to disk.

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

def writeVtkWorms(filename,points,lines,vals,single_level=True):
    if not single_level:
        # Warning. This code is buggy and produces a messed up vtk file
        # It displays, but there are all kinds of problems in it
        # FIXME: cure this stanza or junk it...
        points_list = [p for v in points.values() for p in v]
        lines_list = [p for v in lines.values() for p in v]
        vals_list = [p for v in vals.values() for p in v]
    else:
        points_list = points
        lines_list = lines
        vals_list = vals
    
    vtk = PV.VtkData(PV.PolyData(points=points_list,
                                 lines=lines_list),
                     'Worm Segments',
                     PV.PointData(PV.Scalars(np.fabs(vals_list),name='mutliscale edge magnitudes')),
                     'Edge Magnitudes'
                    )
    vtk.tofile(filename+"_worms.vtk")
    
def writeVtkWormLevels(filename,points,lines,vals):
    assert len(points) == len(lines)
    assert len(lines) == len(vals)
    for level in range(1,len(points)+1):
        pts = points[level]
        lns = lines[level]
        vs = vals[level]
        name = filename + '_level_{0}'.format(level)
        writeVtkWorms(name,pts,lns,vs)

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
    
"""The following code comes from 
<https://gis.stackexchange.com/questions/57834/how-to-get-raster-corner-coordinates-using-python-gdal-bindings>
and supposedly comes from the metageta project <https://code.google.com/p/metageta/>
which is MIT licensed but with the "MetaGETA name substituted for MIT..
I reproduce the text of the MIT license here to (hopefully) remain in 
compliance with the terms of that license.

MetaGETA license

Copyright (c) 2013 Australian Government, Department of the Environment

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

def GetExtent(gt,cols,rows):
    ''' Return list of corner coordinates from a geotransform

        @type gt:   C{tuple/list}
        @param gt: geotransform
        @type cols:   C{int}
        @param cols: number of columns in the dataset
        @type rows:   C{int}
        @param rows: number of rows in the dataset
        @rtype:    C{[float,...,float]}
        @return:   coordinates of each corner
    '''
    ext=[]
    xarr=[0,cols]
    yarr=[0,rows]

    for px in xarr:
        for py in yarr:
            x=gt[0]+(px*gt[1])+(py*gt[2])
            y=gt[3]+(px*gt[4])+(py*gt[5])
            ext.append([x,y])
            print(x,y)
        yarr.reverse()
    return ext

def ReprojectCoords(coords,src_srs,tgt_srs):
    ''' Reproject a list of x,y coordinates.

        @type geom:     C{tuple/list}
        @param geom:    List of [[x,y],...[x,y]] coordinates
        @type src_srs:  C{osr.SpatialReference}
        @param src_srs: OSR SpatialReference object
        @type tgt_srs:  C{osr.SpatialReference}
        @param tgt_srs: OSR SpatialReference object
        @rtype:         C{tuple/list}
        @return:        List of transformed [[x,y],...[x,y]] coordinates
    '''
    trans_coords=[]
    transform = osr.CoordinateTransformation( src_srs, tgt_srs)
    for x,y in coords:
        x,y,z = transform.TransformPoint(x,y)
        trans_coords.append([x,y])
    return trans_coords

""" This is converted into a doctest...
    >>> raster=r'somerasterfile.tif'
    >>> ds=gdal.Open(raster)

    >>> gt=ds.GetGeoTransform()
    >>> cols = ds.RasterXSize
    >>> rows = ds.RasterYSize
    >>> ext=GetExtent(gt,cols,rows)

    >>> src_srs=osr.SpatialReference()
    >>> src_srs.ImportFromWkt(ds.GetProjection())
    >>> #tgt_srs=osr.SpatialReference()
    >>> #tgt_srs.ImportFromEPSG(4326)
    >>> tgt_srs = src_srs.CloneGeogCS()

    >>> geo_ext=ReprojectCoords(ext,src_srs,tgt_srs)
"""

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    
