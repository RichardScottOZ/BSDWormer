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
from scipy import spatial
import networkx as nx




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
        self.worm_image = fdo.CannyEdgeDetect(up_fdg)
        return self.worm_image
    
    def buildGraph(self,seg,vals):
        self.G = nx.Graph()
        tree = spatial.KDTree(seg)
        for i,s in enumerate(seg):
            neighbors = tree.query(s,k=3,distance_upper_bound=1.42)
            self.G.add_node(i,pos=tuple(s),val=vals[i])
            n1 = neighbors[1][1]
            self.G.add_edge(i,n1)
            n2 = neighbors[1][2]
            self.G.add_edge(i,n2)
        
    def buildWormSegs(self):
        worm_points = np.argwhere(self.worm_image > -100)
        worm_vals = self.worm_image[worm_points[:,0],worm_points[:,1]]
        self.buildGraph(worm_points,worm_vals)
        mst = nx.minimum_spanning_tree(self.G)
        num_nodes = worm_points.shape[0]
        visited = np.zeros(num_nodes+1,dtype=np.bool_)
        self.segs = []
        for nd in range(num_nodes+1):                       # Loop over all nodes
            if not visited[nd]:                             # Check to see if we've been here before
                dfst = nx.dfs_tree(mst,nd)                  # Build a Depth First Search Tree rooted at nd 
                visited[nd] = True                          # Mark nd as visited
                bingo = []                                  # Start a new list to hold results; will be included as element in master list
                last_node = nd                              # Initialize last_node
                for (source,dest,d) in nx.dfs_labeled_edges(dfst,nd):
                    # Run through the labeled edges of a Depth First Search 
                    if d['dir'] != 'forward':
                        # Backtracking or out of tree edge; not interested
                        continue
                    if source == dest:
                        # Root of tree
                        last_node = dest                    # Initialize last node
                        visited[dest] = True                # and mark dest visited
                        continue
                    if source >= num_nodes or dest >= num_nodes:
                        # FIXME: HACK!!! Got to figure out what's going on here
                        # Although, this does seem to help with fixing
                        # the "slash through everything" edges...
                        # A mystery.
                        continue
                    if last_node != source:
                        # We've switched branches!
                        self.segs += [bingo]                # Add the accumulated edges to segs
                        bingo = [(source,dest)]             # and initialize a new bingo list with this edge
                    else:
                        # We're tracking a single branch (last_node == source)
                        bingo += [(source,dest)]            # Add the new edge to the bingo list
                    last_node = dest
                    visited[dest] = True
                self.segs += [bingo]                        # Add the last bingo list
        
        
if __name__ == '__main__':
    import doctest
    doctest.testmod()