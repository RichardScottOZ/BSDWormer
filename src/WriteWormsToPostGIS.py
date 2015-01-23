'''
Created on Jan 23, 2015

@author: frank
'''
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Float
from geoalchemy2 import Geometry
from sqlalchemy.orm import sessionmaker
from io import StringIO

def ToWKT(foo,SRID=4326):
    dummy = StringIO()
    dummy.write('SRID=%d; MULTILINESTRINGZM('%SRID)
    first_line = True
    for l in foo.all_lines:
        if first_line:
            first_line = False
        else:
            dummy.write(',')
        dummy.write('(')
        last_point_idx = len(l)-1
        for j,lp in enumerate(l):
            # No bloody comma on the end
            dummy.write('%g %g %g %g'%(foo.all_points[lp][0], foo.all_points[lp][1], foo.all_points[lp][2], foo.all_vals[lp]))
            if j != last_point_idx:
                dummy.write(',')
        dummy.write(')')
    dummy.write(')')
    ret = dummy.getvalue()
    dummy.close()
    return ret

Base = declarative_base()

class WormLayer(Base):
    __tablename__ = 'worm_layer'
    id = Column(Integer, primary_key=True)
    height = Column(Float)
    geom = Column(Geometry('MULTILINESTRINGZM'))


class PostGISWriter(object):
    '''
    Encapsulates the method for writing worms to XYZM points in PostGIS. 
    '''

    def __init__(self, db='frank', tablename='worm_layer', srid=4326):
        '''
        Constructor
        '''
        self.engine = create_engine('postgresql://frank@localhost/%s'%db, echo=False)
        Session = sessionmaker(bind=self.engine)
        self.session = Session()
        self.srid = srid
        WormLayer.__tablename__ = tablename
        WormLayer.__table__.create(self.engine)
        
    def addWormLayer(self,layer,height):
        worm_layer = WormLayer(height,geom=ToWKT(layer,SRID=self.srid))
        self.session.add(worm_layer)
        self.session.commit()
        
    def _rollbackBadDatabaseTransaction(self):
        self.session.rollback()
        