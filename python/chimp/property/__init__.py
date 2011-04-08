#!/usr/bin/env python

from lxml.builder import E, ET
from lxml import etree

import os, sys

THIS_DIR = os.path.abspath( os.path.split(__file__)[0] )
from ..config import CHIMP_DATA

import physical
from physical.unit import *
from physical.constant import * # pi == constant.pi
from physical import unit       # so we can get unit.pi
from physical import element


xmlDb = etree.parse( CHIMP_DATA )
xmlDb.xinclude()

def getMass(P):
    """
    Returns the mass of a Particle if known.
    """
    _m = xmlDb.xpath("//Particle[@name='{name}']/mass/node()".format(name=P))[0]
    # translate from my c++ language to python
    _m = eval( _m.replace('::', '.') )
    return _m.coeff

def getCharge(P):
    """
    Returns the charge of a Particle if known.
    """
    _m = xmlDb.xpath( "//Particle[@name='{name}']/charge/node()" \
                      .format(name=P) )[0]
    # translate from my c++ language to python
    _m = eval( _m.replace('::', '.') )
    return _m.coeff

def getPolarizability(P):
    """
    Returns the Polarizability of a Particle if known.
    """
    _m = xmlDb.xpath( "//Particle[@name='{name}']/polarizability/node()" \
                      .format(name=P) )[0]
    # translate from my c++ language to python
    _m = eval( _m.replace('::', '.') )
    return _m.coeff

