import numpy as np
#import numpy.ma as ma
#from osgeo import gdalnumeric
#from osgeo import gdal
#from osgeo import gdalconst
from matplotlib import pyplot as plt
import os.path
#import FourierDomainGrid as GRID
#import FourierDomainOps as FDO
#from Utility import isclose, viewRaster, writeVtkImage, writeVtkWorms, GetExtent, ReprojectCoords
import wormer as w
#import networkx as nx
#from scipy import spatial
from WriteWormsToPostGIS import PostGISWriter as PGW
#from geometry import GeoTransformToGCPs,MapToPixel
import sys
job = w.Wormer()
filename = os.path.abspath('../test_data/AppBasinGPFA/ADK_BGA_UTM_2500.ers')
job.importGdalRaster(filename)
job.importExternallyPaddedRaster('../test_data/AppBasinGPFA/ADK_BGA_UTM_2500_padded.ers')
#raw_input()
#job.padded_grid = job.base_grid
#nodata = -65.
#foo.base_grid = np.where(foo.base_grid < nodata,nodata,foo.base_grid)
#job.buildPaddedRaster((512,1024))
#job._georeferencePaddedGrid()
# HACK! We already have a padded PSG image from Oasis Montaj.
#job.padded_geotransform = job.geomat
pgw = PGW(srid=32618,db='postgresql://frank@localhost:5432/frank')
for dz in range(1,15):
    #if dz == 0:
    #    dz = 0.001
    delta_z = 500.
    #dzm = -dz*job.dy
    # Bloody images. job.dy is negative. Grrrr. So we need to have a positive dzm.
    dzm = dz * delta_z
    job.wormLevelAsPoints(dz=dzm)
    job.buildWormSegs(dz=dz,clipped=True,nodata_in_worm_image=-100.,from_image=False)
    job.buildLevelForVTK(dz=dz,delta_z_in_units=dzm)
    #pgw.addWormLayer(job,dz)
    #pgw.addWormPoints(job,dz)
    pgw.addWormsAtHeightToDB(job,dz,srid=32618)
# This thing does some SRID postprocessing on the bloody PostGIS database. 
# Why? Because I couldn't figure out any other way to make it work with geoalchemy2...
pgw.cleanUpDatabase()
