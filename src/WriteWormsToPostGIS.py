"""
Created on Jan 23, 2015

 @author: frank
"""
from sqlalchemy import create_engine, func, inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, Float, ForeignKey
from geoalchemy2 import Geometry
from geoalchemy2.elements import WKTElement, WKBElement
from sqlalchemy.orm import sessionmaker, relationship, backref
from sqlalchemy.sql import text
from io import StringIO
from math import atan2, pi

# This is the base of all PostGIS table names for this project
# With a little luck, all of this "by hand" construction of tablenames
# will get fixed in the worming code shortly, but for now, let's keep on doing this.
basename = 'ADK_CEUS_PSG_1250'
layer_name = basename
points_name = basename + '_points'
levels_name = basename + '_levels'
levels_points_name = basename + '_levels_points'


def ToMLSZM_WKT(foo,height,SRID=4326):
    dummy = StringIO()
    dummy.write('SRID=%d; MULTILINESTRINGZM('%SRID)
    first_line = True
    for l in foo.all_lines[height]:
        if first_line:
            first_line = False
        else:
            dummy.write(',')
        dummy.write('(')
        last_point_idx = len(l)-1
        for j,lp in enumerate(l):
            # No bloody comma on the end
            dummy.write('%g %g %g %g'%(foo.all_points[height][lp][0], 
                                       foo.all_points[height][lp][1], 
                                       foo.all_points[height][lp][2], 
                                       foo.all_vals[height][lp]))
            if j != last_point_idx:
                dummy.write(',')
        dummy.write(')')
    dummy.write(')')
    ret = dummy.getvalue()
    dummy.close()
    return ret

def makeLine(pt1,pt2):
    geom = func.ST_AsEWKB(func.ST_MakeLine(pt1,pt2))
    return geom

Base = declarative_base()

class WormLayer(Base):
    __tablename__ = basename
    id = Column(Integer, primary_key=True)
    height = Column(Float)
    geom = Column(Geometry('MULTILINESTRINGZM'))
    
class WormPoint(Base):
    __tablename__ = points_name
    # Primary Key. Boring.
    worm_point_id = Column(Integer, primary_key=True, index=True)
    # An id from the worming code. I don't remember if it was unique, so I didn't use it as a PK.
    vtk_id = Column(Integer,index=True)
    # Coordinates of the point in some 'native' CRS
    x = Column(Float)
    y = Column(Float)
    z = Column(Float)
    # The scalar value of the magnitude of the horizontal gradient
    grad = Column(Float)
    # The height of upward continuation from which the grad and coordinates were drawn.
    height = Column(Float)
    # A PostGIS point geometry, in the native CRS
    pt = Column(Geometry('POINTZM'),index=True)
    # A duplicate of pt in WGS84 coordinates; converted by PostGIS at write-time
    wgs84_pt = Column(Geometry('POINTZM'),index=True)
    # Database magic that links entries in this table with entries in another table
    level = relationship('WormLevel', secondary=levels_points_name)
    
class WormLevel(Base):
    __tablename__ = levels_name
    # A PK, there are only ~10 entries in this table, so it's tiny, so no index.
    worm_level_id = Column(Integer, primary_key=True)
    # The actual level (prob in meters, but potentially varies...)
    level = Column(Float)
    # Database magic that links entries in this table with entries in another table
    point = relationship('WormPoint', secondary=levels_points_name)
    
class WormLevelPoints(Base):
    __tablename__ = levels_points_name
    # This table has a "composite primary key" composed of the first 2 ForeignKey entries and the internal primary key
    # This is the level_id in the external table
    worm_level_id = Column(Integer, ForeignKey(levels_name + '.worm_level_id'), primary_key=True)
    # This is the point id of the END point of a line segment.
    point_id = Column(Integer, ForeignKey(points_name + '.worm_point_id'), primary_key=True)
    # In addition to participating in a composite primary key, this field is 
    # a unique-within-a-level index for worm segments. 
    worm_seg_id = Column(Integer,primary_key=True,index=True)
    # Database magic that links entries in this table with entries in another table
    worm_level = relationship(WormLevel, backref=backref("worm_point_assoc"))
    # Database magic that links entries in this table with entries in another table
    worm_point = relationship(WormPoint, backref=backref("worm_level_assoc"))
    # This is an index number internal to each worm segment, numbering the edges
    # FIXME (maybe) This terminology needs to be cleaned up.
    seg_sequence_num = Column(Integer)
    # This holds the PostGIS geometry structure for a single edge, in some native CRS.
    line_segmt = Column(Geometry('LINESTRINGZM'),index=True)
    # This scalar gradient value is derived from the average of the point grads on either end of the edge
    # Currently, the upstream code is doing that for the LOG(value), so this is in fact now
    # sqrt(grad(pt1) * grad(pt2))
    line_grad = Column(Float)
    # The azimuth of the edge in degrees East of North.
    azimuth = Column(Float)
    # This is the point ID in the points table of the starting point of an edge
    # FIXME (maybe) this could and probably should be an actual relation into the points table, for ease of retrieval.
    start_point_id = Column(Integer)
    # This is a duplicate of line_segmt but explicitly stored in wgs84.
    wgs84_line_segmt = Column(Geometry('LINESTRINGZM'),index=True)

    



class PostGISWriter(object):
    """
    Encapsulates the method for writing worms to XYZM points in PostGIS.
    Also writes the x, y, z, and m values to a different table.
    """

    def __init__(self, db="postgresql://frank:f00bar@localhost:5433/frank", srid=4326):
        """
        Constructor
        """
        self.engine = create_engine('%s'%db, echo=False)
        Session = sessionmaker(bind=self.engine)
        self.session = Session()
        self.connect = self.engine.connect()
        self.srid = srid
        if not self.engine.dialect.has_table(self.engine.connect(), layer_name):
            WormLayer.__table__.create(self.engine)
        if not self.engine.dialect.has_table(self.engine.connect(), points_name):
            WormPoint.__table__.create(self.engine)
        if not self.engine.dialect.has_table(self.engine.connect(), levels_name):
            WormLevel.__table__.create(self.engine)
        if not self.engine.dialect.has_table(self.engine.connect(), levels_points_name):
            WormLevelPoints.__table__.create(self.engine)
        
        
    def addWormLayer(self,layers,height):
        worm_layer = WormLayer(height=height,
                               geom=ToMLSZM_WKT(layers,height=height,SRID=self.srid))
        self.session.add(worm_layer)
        self.session.commit()

    def addWormPoints(self,layers,height):
        for l in layers.all_lines[height]:
            for lp in l:
                x = layers.all_points[height][lp][0]
                y = layers.all_points[height][lp][1]
                z = layers.all_points[height][lp][2]
                grad = layers.all_vals[height][lp]
                pt_wkt = 'POINTZM(%g %g %g %g)'%(x,y,z,grad)
                worm_point = WormPoint(vtk_id=lp,
                                       x=x,
                                       y=y,
                                       z=z,
                                       grad=grad,
                                       height=height,
                                       pt=pt_wkt)
                self.session.add(worm_point)
        self.session.commit()

    def addWormsAtHeightToDB(self,layers,height, srid=4326):
        
        # Build the level record in the database
        worm_level = WormLevel(level=height)
        self.session.add(worm_level)
        self.session.commit()
        
        # Load all of the points records into the database
        points = {}
        for pt_id,pt in enumerate(layers.all_points[height]):
            x = pt[0]
            y = pt[1]
            z = pt[2]
            grad = layers.all_vals[height][pt_id]
            pt_wkt = 'POINTZM(%g %g %g %g)'%(x,y,z,grad)
            pt_ewkt = WKTElement(pt_wkt,srid=srid)
            wgs84_pt = func.ST_Transform(pt_ewkt,4326)
            worm_point = WormPoint(vtk_id=pt_id,
                                   x=x,
                                   y=y,
                                   z=z,
                                   grad=grad,
                                   height=height,
                                   pt=pt_ewkt,
                                   wgs84_pt = wgs84_pt
                                   )
            self.session.add(worm_point)
            points[pt_id] = worm_point
            
        self.session.commit()

        # Grab the table object so that we can do simple inserts instead of through the 
        # molasses-slow declarative_base route for inserts...
        wlp_table = inspect(WormLevelPoints).mapped_table
        
        # And finally, load all of the association table records
        connect = self.session.connection()
        r1 = connect.execute(wlp_table.select())
        for seg_id,l in enumerate(layers.all_lines[height]):
            start_pt = points[l[0]]
            sp = layers.all_points[height][start_pt.vtk_id]
            for seq_num,lp in enumerate(l):
                # seq_num == 0 is a special case
                # sp is always equal to ep, due to a vagary of
                # the network/graph constructing algorithm.
                # Let's just not even bother putting them into the database
                if seq_num == 0:
                    continue
                end_pt = points[lp]
                ep = layers.all_points[height][end_pt.vtk_id]
                line_grad = (start_pt.grad + end_pt.grad)/2.0
                # FIXME (Maybe). This assumes that all coords are in a UTM or
                # otherwise Cartesian coordinate system. This should be true
                # for the worms, but BEWARE if you happen to be
                # performing your worm job with geographic coordinates.
                # We are calculating the azimuth (East from North)
                # from the start point to the end point.
                dEastings = ep[0] - sp[0]
                dNorthings = ep[1] - sp[1]
                azimuth = atan2(dEastings,dNorthings)*360./(2.*pi)
                sgmt_wkt = 'LINESTRINGZM(%g %g %g %g, %g %g %g %g)'%(sp[0],
                                                                     sp[1],
                                                                     sp[2],
                                                                     start_pt.grad,
                                                                     ep[0],
                                                                     ep[1],
                                                                     ep[2],
                                                                     end_pt.grad)
                sgmt_ewkt = WKTElement(sgmt_wkt,srid=srid)
                #wgs84_sgmt = func.ST_Transform(sgmt_ewkt,4326)
                connect.execute(wlp_table.insert(),
                                worm_level_id = worm_level.worm_level_id,
                                point_id = end_pt.worm_point_id,
                                seg_sequence_num = seq_num,
                                worm_seg_id = seg_id,
                                line_segmt = sgmt_ewkt,
                                line_grad = line_grad,
                                azimuth = azimuth,
                                start_point_id = start_pt.worm_point_id
                                )
                start_pt = end_pt
                sp = ep
        
        self.session.commit()
        
                
    def cleanUpDatabase(self):
        connect = self.session.connection()
        connect.execute(text('UPDATE "'+levels_points_name+'" SET wgs84_line_segmt = ST_Transform(line_segmt,4326);'))
        connect.execute(text("SELECT UpdateGeometrySRID('"+levels_points_name+"','wgs84_line_segmt',4326);"))
        self.session.commit()
        self.session.close()
        

    def _rollbackBadDatabaseTransaction(self):
        self.session.rollback()
        