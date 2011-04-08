import os, sys

THIS_DIR = os.path.abspath( os.path.split(__file__)[0] )

try:
  import physical
  a = physical.unit.m
  del a
except:
  PHYS_DIR = os.path.abspath( THIS_DIR + '/../../../../physical/python' )
  sys.path.append( PHYS_DIR )

CHIMP_DATA = os.path.abspath( THIS_DIR + '/../../data/particledb.xml' )

