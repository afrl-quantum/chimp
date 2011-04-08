
from lxml.builder import E as _E

class DATA:
    def __init__(self, reference, xscale=1, yscale=1, data=list()):
        self.reference = reference
        self.xscale = xscale
        self.yscale = yscale
        self.data = data


    def getXML(self):
        x_cs = _E.cross_section( {'model'  : 'data',
                                 'xscale' : str(self.xscale),
                                 'yscale' : str(self.yscale) } )

        vals = [ _E.val({'x':str(d[0]), 'y':str(d[1])}) for d in self.data ]
        x_cs.extend( vals )

        if self.reference:
          x_cs.extend( [self.reference.getXML()] )

        return x_cs


class VHS:
    def __init__(self, reference, value, T_ref, visc_T_law):
        self.reference = reference
        self.value = value
        self.T_ref = T_ref
        self.visc_T_law = visc_T_law

    def getXML(self):
        x_cs = _E.cross_section(
            {'model' : 'vhs'},
            _E.value( str(self.value) ),
            _E.T_ref( str(self.T_ref) ),
            _E.visc_T_law( str(self.visc_T_law) )
        )

        if self.reference:
          x_cs.extend( [self.reference.getXML()] )

        return x_cs



class LotzParameter:
    def __init__(self, P, q, a, b, c):
        self.P = P
        self.q = q
        self.a = a
        self.b = b
        self.c = c

    def getXML(self):
        return _E.LotzParameter(
            _E.P( self.P ),
            _E.q( self.q ),
            _E.a( self.a ),
            _E.b( self.b ),
            _E.c( self.c ),
        )

class Lotz:
    def __init__(self, reference, xscale=1, yscale=1, params=list()):
        self.reference = reference
        self.xscale = xscale
        self.yscale = yscale
        self.params = params


    def getXML(self):
        x_cs = _E.cross_section( {'model'  : 'lotz'} )

        x_params = [ p.getXML() for p in self.params ]
        x_cs.extend( x_params )

        if self.reference:
          x_cs.extend( [self.reference.getXML()] )

        return x_cs

