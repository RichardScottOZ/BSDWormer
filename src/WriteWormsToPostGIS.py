'''
Created on Jan 23, 2015

 @author: frank
'''
from sqlalchemy import create_engine, func, inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, Float, ForeignKey
from geoalchemy2 import Geometry
from geoalchemy2.elements import WKTElement, WKBElement
from sqlalchemy.orm import sessionmaker, relationship, backref
from sqlalchemy.sql import text
from io import StringIO
from math import atan2, pi


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
    __tablename__ = 'AppBasinMergedBGA2500'
    id = Column(Integer, primary_key=True)
    height = Column(Float)
    geom = Column(Geometry('MULTILINESTRINGZM'))
    
class WormPoint(Base):
    __tablename__ = 'AppBasinMergedBGA2500_points'
    worm_point_id = Column(Integer, primary_key=True, index=True)
    vtk_id = Column(Integer,index=True)
    x = Column(Float)
    y = Column(Float)
    z = Column(Float)
    grad = Column(Float)
    height = Column(Float)
    pt = Column(Geometry('POINTZM'),index=True)
    wgs84_pt = Column(Geometry('POINTZM'),index=True)
    level = relationship('WormLevel', secondary='AppBasinMergedBGA2500_levels_points')
    
class WormLevel(Base):
    __tablename__ = 'AppBasinMergedBGA2500_levels'
    worm_level_id = Column(Integer, primary_key=True)
    level = Column(Float)
    point = relationship('WormPoint', secondary='AppBasinMergedBGA2500_levels_points')
    
class WormLevelPoints(Base):
    __tablename__ = 'AppBasinMergedBGA2500_levels_points'
    worm_level_id = Column(Integer, ForeignKey('AppBasinMergedBGA2500_levels.worm_level_id'), primary_key=True)
    point_id = Column(Integer, ForeignKey('AppBasinMergedBGA2500_points.worm_point_id'), primary_key=True)
    worm_seg_id = Column(Integer,primary_key=True,index=True)
    worm_level = relationship(WormLevel, backref=backref("worm_point_assoc"))
    worm_point = relationship(WormPoint, backref=backref("worm_level_assoc"))
    seg_sequence_num = Column(Integer)
    line_segmt = Column(Geometry('LINESTRINGZM'),index=True)
    line_grad = Column(Float)
    azimuth = Column(Float)
    start_point_id = Column(Integer)
    wgs84_line_segmt = Column(Geometry('LINESTRINGZM'),index=True)

    



class PostGISWriter(object):
    '''
    Encapsulates the method for writing worms to XYZM points in PostGIS.
    Also writes the x, y, z, and m values to a different table.
    '''

    def __init__(self, db='postgresql://frank@localhost/frank', srid=4326):
        '''
        Constructor
        '''
        self.engine = create_engine('%s'%db, echo=False)
        Session = sessionmaker(bind=self.engine)
        self.session = Session()
        self.connect = self.engine.connect()
        self.srid = srid
        if not self.engine.dialect.has_table(self.engine.connect(), "AppBasinMergedBGA2500"):
            WormLayer.__table__.create(self.engine)
        if not self.engine.dialect.has_table(self.engine.connect(), "AppBasinMergedBGA2500_points"):
            WormPoint.__table__.create(self.engine)
        if not self.engine.dialect.has_table(self.engine.connect(), "AppBasinMergedBGA2500_levels"):
            WormLevel.__table__.create(self.engine)
        if not self.engine.dialect.has_table(self.engine.connect(), "AppBasinMergedBGA2500_levels_points"):
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
        connect.execute(text('UPDATE "AppBasinMergedBGA2500_levels_points" SET wgs84_line_segmt = ST_Transform(line_segmt,4326);'))
        connect.execute(text("SELECT UpdateGeometrySRID('AppBasinMergedBGA2500_levels_points','wgs84_line_segmt',4326);"))
        self.session.commit()
        self.session.close()
        

    def _rollbackBadDatabaseTransaction(self):
        self.session.rollback()
        