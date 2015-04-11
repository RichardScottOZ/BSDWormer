import numpy as np
from osgeo import gdalnumeric
from osgeo import gdal
from osgeo import gdalconst
#from matplotlib import pyplot as plt
#import os.path
import FourierDomainGrid as GRID
import FourierDomainOps as FDO
#from Utility import isclose, viewRaster
#from scipy.ndimage.measurements import label
from scipy import spatial
import networkx as nx
from Utility import writeVtkWorms, writeGDALRasterFromNumpyArray
from FftUtils import mk_apod_mask, embed_data
from osgeo.gdal import GCP, GCPsToGeoTransform
from geometry import MapToPixel, PixelToMap, CellSize, ExtentToGCPs, GeoTransformToGCPs


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
    
    """
    
    def __init__(self):
        self.base_grid = None
        self.all_points = {}
        self.all_lines = {}
        self.all_vals = {}
        self.G = {}
        
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
        try:
            self.no_data_value = self.ds.GetRasterBand(1).GetNoDataValue()
        except:
            self.no_data_value = None
        # FIXME: the following assumes North is up in the image
        self.dx = self.geomat[1]
        self.dy = self.geomat[5]
        
    def exportNewGdalRaster(self, narray, gdal_filename, fmt='ERS'):
        """Export a numpy array into a GDAL image, using the initial image as a prototype.
        """
        gdalnumeric.SaveArray( narray, gdal_filename, format = fmt, prototype = self.ds )
        
    def wormLevelAsImage(self,dz):
        # import pdb; pdb.set_trace()
        fdg = GRID.FourierDomainGrid(dx=self.dx, dy=self.dy)
        fdg.setSpatialGrid(self.padded_grid)
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
    
    def buildGraph(self,seg,vals,geotransform,dz):
        self.G[dz] = nx.Graph()
        tree = spatial.KDTree(seg)
        for i,s in enumerate(seg):
            neighbors = tree.query(s,k=3,distance_upper_bound=1.42)
            geographic = PixelToMap(s[1],s[0],geotransform)
            self.G[dz].add_node(i,pos=tuple(s),val=vals[i],geog_pos=geographic)
            n1 = neighbors[1][1]
            self.G[dz].add_edge(i,n1)
            n2 = neighbors[1][2]
            self.G[dz].add_edge(i,n2)
        # remove duplicate edges 
        # <http://networkx.lanl.gov/reference/generated/networkx.generators.degree_seq.configuration_model.html>
        self.G[dz] = nx.Graph(self.G[dz])
        
    def buildWormSegs(self,nodata_in_worm_image = -100, clipped=True, log_vals = True, dz=None):            
        if clipped:
            # The coords returned by np.argwhere are setting the start of 
            # the clipping region to coordinates (0,0)
            # Deal with that by pulling vals from the clipped region...
            #clip_img = self.worm_image[self.padded_slice_y,self.padded_slice_x]
            clip_img = np.where(self.externally_sized_mask > self.no_data_value,
                                self.worm_image,self.no_data_value)[self.padded_slice_y,self.padded_slice_x]
                
            worm_points = np.argwhere(clip_img > nodata_in_worm_image )
            worm_vals = clip_img[worm_points[:,0],worm_points[:,1]]
            gt = self.geomat
        else:
            worm_points = np.argwhere(self.worm_image > nodata_in_worm_image)
            worm_vals = self.worm_image[worm_points[:,0],worm_points[:,1]]
            gt = self.padded_geotransform
        worm_vals *= dz
        if log_vals:
            worm_vals = np.log10(worm_vals)
        self.buildGraph(worm_points,worm_vals,gt,dz)
        mst = nx.minimum_spanning_tree(self.G[dz])
        num_nodes = worm_points.shape[0]
        visited = np.zeros(num_nodes+1,dtype=np.bool_)
        self.segs = []
        for nd in range(num_nodes):                       # Loop over all nodes
            if not visited[nd]:                             # Check to see if we've been here before
                dfst = nx.dfs_tree(mst,nd)                  # Build a Depth First Search Tree rooted at nd 
                visited[nd] = True                          # Mark nd as visited
                bingo = []                                  # Start a new list to hold results; will be included as element in master list
                last_node = nd                              # Initialize last_node
                for (source,dest,d) in nx.dfs_labeled_edges(dfst,nd):
                    # Run through the labeled edges of a Depth First Search
                    if visited[source] and visited[dest]:
                        # We've already traversed this edge for some reason; bail out of further processing
                        continue
                    if d['dir'] != 'forward':
                        # Backtracking or out-of-spanning-tree edge; not interested
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
                
        
    
    def buildLevelForVTK(self,dz,invert_z = True):
        """Computes the VTK representation of worms for a single level."""
        # Organize the worm data so that the VTK writing stuff will swallow it...
        G_dz = dz
        max_valid_node = self.G[G_dz].number_of_nodes()-1
        #max_valid_node = self.G[G_dz].number_of_nodes()
        all_points = [self.G[G_dz].node[p]['geog_pos'] for p in range(max_valid_node)] # The last one is empty; in Geographic coords

        if invert_z:
            dz = -dz
            
        all_points = [(e[0],e[1],dz*self.dx) for e in all_points] # Make sure we get East and North correct
        if self.all_points.get(G_dz) == None:
            base_point_num = 0
            self.all_points[G_dz] = all_points
        else:
            base_point_num = len(self.all_points[G_dz]) 
            self.all_points[G_dz] += all_points
        all_vals = [self.G[G_dz].node[p]['val'] for p in range(max_valid_node)] # The last one is empty
        all_lines = []
        for s in self.segs:
            if s == []:
                continue
            if s[0][0] >= max_valid_node or s[0][1] >= max_valid_node:
                continue
            # We need to int-ify things here so we don't trip
            # over an isinstance(foo,int) with foo as numpy.int64 in the
            # VTK writer checks since that fails... Go figure.
            ln = [int(s[0][0]+base_point_num),int(s[0][1]+base_point_num)]
            for edge in s[1:]:
                if edge[1] >= max_valid_node:
                    break
                # Ditto...
                ln += [int(edge[1]+base_point_num)]
            all_lines += [ln]
            
        if self.all_vals.get(G_dz) == None:
            self.all_vals[G_dz] = all_vals
        else:
            self.all_vals[G_dz] += all_vals
            
        if self.all_lines.get(G_dz) == None:
            self.all_lines[G_dz] = all_lines
        else:
            self.all_lines[G_dz] += all_lines


    def buildPaddedRaster(self,padded_shape,rolloff_size = 100,pad_type='hann'):
        """Returns an 'apodized' padded version of base_grid.
        
        #FIXME: Turn this into (doc?)tests...
        
        assert(np.allclose(baz[slice_y,slice_x],foo.base_grid-np.mean(foo.base_grid)))
        assert(np.allclose((baz*bar)[slice_y,slice_x],foo.base_grid-np.mean(foo.base_grid)))
        
        # Deal with the top left corner, and surrounding pixels
        top_left = (slice_y.start,slice_x.start)
        print(top_left)
        assert(bar[top_left]==1.0)
        nw_top_left = (slice_y.start-1,slice_x.start-1)
        print(nw_top_left)
        assert(bar[nw_top_left] < 1.0)
        n_top_left = (slice_y.start-1,slice_x.start)
        print(n_top_left)
        assert(bar[n_top_left] < 1.0)
        w_top_left = (slice_y.start,slice_x.start-1)
        print(w_top_left)
        assert(bar[w_top_left] < 1.0)
        
        # Deal with the top right corner, and surrounding pixels
        top_right = (slice_y.start,slice_x.stop)
        print(top_right)
        assert(bar[top_right]==1.0)
        e_top_right = (slice_y.start,slice_x.stop+1)
        print(e_top_right)
        assert(bar[e_top_right]<1.0)
        ne_top_right = (slice_y.start-1,slice_x.stop+1)
        print(ne_top_right)
        assert(bar[ne_top_right]<1.0)
        n_top_right = (slice_y.start-1,slice_x.stop)
        print(n_top_right)
        assert(bar[n_top_right]<1.0)
        
        # Deal with the bottom left corner, and surrounding pixels
        bot_left = (slice_y.stop,slice_x.start)
        print(bot_left)
        assert(bar[bot_left]==1.0)
        sw_bot_left = (slice_y.stop+1,slice_x.start-1)
        print(sw_bot_left)
        assert(bar[sw_bot_left] < 1.0)
        s_bot_left = (slice_y.stop+1,slice_x.start)
        print(s_bot_left)
        assert(bar[s_bot_left] < 1.0)
        w_bot_left = (slice_y.stop,slice_x.start-1)
        print(w_bot_left)
        assert(bar[w_bot_left] < 1.0)
        
        # Deal with the bottom right corner, and surrounding pixels
        bot_right = (slice_y.stop,slice_x.stop)
        print(bot_right)
        assert(bar[bot_right]==1.0)
        e_bot_right = (slice_y.stop,slice_x.stop+1)
        print(e_bot_right)
        assert(bar[e_bot_right]<1.0)
        se_bot_right = (slice_y.stop+1,slice_x.stop+1)
        print(se_bot_right)
        assert(bar[se_bot_right]<1.0)
        s_bot_right = (slice_y.stop+1,slice_x.stop)
        print(s_bot_right)
        assert(bar[s_bot_right]<1.0)
        """
        ysz,xsz = self.base_grid.shape
        rolloff_sz = rolloff_size
        bar = mk_apod_mask(masksz=padded_shape,
                   apodsz=(ysz+rolloff_sz,xsz+rolloff_sz),
                   wsize=rolloff_sz,
                   apod_f = pad_type)
        baz,p_slice_y,p_slice_x = embed_data(self.base_grid,
                                             big_shape=padded_shape,
                                             pad_size=(rolloff_sz//2))
        self.padded_grid = baz*bar
        self.padded_slice_y = p_slice_y
        self.padded_slice_x = p_slice_x
        
    def _georeferencePaddedGrid(self):
        """Builds a spatial reference system appropriate for a padded grid."""
        ypix_max,xpix_max = self.base_grid.shape
        gcps = GeoTransformToGCPs(self.geomat,xpix_max,ypix_max)
        padded_gcps = []
        for gcp in gcps:
            if gcp.GCPPixel == 0.0:
                x_pad = self.padded_slice_x.start
            else:
                x_pad = self.padded_slice_x.stop
        
            if gcp.GCPLine == 0.0:
                y_pad = self.padded_slice_y.start
            else:
                y_pad = self.padded_slice_y.stop
        
            padded_gcp = GCP(gcp.GCPX,gcp.GCPY,gcp.GCPZ,x_pad,y_pad)
            padded_gcps += [padded_gcp]
        # This is a GDAL geotransform for the padded grid
        # in the Spatial Reference System of the original dataset
        self.padded_geotransform = GCPsToGeoTransform(padded_gcps)
        
    def _writePaddedGridGeotiff(self,filename):
        """Writes the padded grid out as a geotiff raster."""
        writeGDALRasterFromNumpyArray(filename,
                                      self.padded_grid,
                                      self.padded_geotransform,
                                      self.ds.GetProjection())
        
    def importExternallyPaddedRaster(self,gdal_filename):
        """Import a GDAL raster of the externally padded image 
           (e.g. by preprocessing with Oasis Montaj) into a numpy array stashed away as an attribute.
        """
        self.gdal_external_padded_filename = gdal_filename
        self.padded_grid = np.array(gdalnumeric.LoadFile(self.gdal_external_padded_filename))
        ep_ds = gdal.Open(self.gdal_external_padded_filename,gdalconst.GA_ReadOnly) 
        ep_geomat = ep_ds.GetGeoTransform()
        # FIXME: the following assumes North is up in the image
        #ep_dx = ep_geomat[1]
        #ep_dy = ep_geomat[5]
        #ep_linemax,ep_pixmax = self.padded_grid.shape
        internal_linemax,internal_pixmax = self.base_grid.shape
        internal_corner_gcps = GeoTransformToGCPs(self.geomat,internal_pixmax,internal_linemax)
        # FIXME: This is a terrapixel image (1 million by 1 million)
        # that _should_ hold us for a while as an outer bound of what we can import into 
        # any machine with "reasonable" memory.
        # However, eventually something more general is desirable...
        internal_pixmin = 1000000
        internal_pixmax = -1000000
        internal_linemin = 1000000
        internal_linemax = -1000000
        for gcp in internal_corner_gcps:
            px,py = MapToPixel(gcp.GCPX,gcp.GCPY,ep_geomat)
            internal_pixmin = min(internal_pixmin,px)
            internal_pixmax = max(internal_pixmax,px)
            internal_linemin = min(internal_linemin,py)
            internal_linemax = max(internal_linemax,py)
        #print(internal_pixmin,internal_pixmax,internal_linemin,internal_linemax)
        self.padded_slice_x = slice(internal_pixmin,internal_pixmax,1)
        self.padded_slice_y = slice(internal_linemin,internal_linemax,1)
        #print(self.padded_slice_x,self.padded_slice_y)
        self.padded_geotransform = ep_geomat
        if self.no_data_value != None:
            self.externally_sized_mask = np.ones(self.padded_grid.shape,np.float32)
            self.externally_sized_mask[:,:] = self.no_data_value
            self.externally_sized_mask[self.padded_slice_y,self.padded_slice_x] = self.base_grid
        

        
        
if __name__ == '__main__':
    import doctest
    doctest.testmod()