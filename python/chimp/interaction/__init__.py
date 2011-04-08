
from lxml.builder import E as _E
import cross_section
from chimp.property import getMass

def sortPartialEq(partial):
    terms = []
    for i in partial.items():
        m = getMass(i[0])
        terms.append( ( m, i[0], i[1] ) )
    terms.sort()
    return terms

def formatPartial(partial):
    S = sortPartialEq(partial)
    strs = []
    for term in S:
        if term[2] > 1: # if n > 1
            strs.append( str(term[2]) + ' ' + term[1] )
        else:
            strs.append( term[1] )
    return ' + '.join(strs)

def xmlPartial(partial):
    S = sortPartialEq(partial)
    terms = []
    for term in S:
        if len(terms):
            terms[-1].tail = ' + '

        if term[2] > 1: # if n > 1
            terms.append( _E.T( _E.n(str(term[2])), ' ', _E.P(term[1]) ) )
        else:
            terms.append( _E.T( _E.P(term[1]) ) )
    return terms


class Equation:
    def __init__(self, lhs=dict(), rhs=dict()):
        self.lhs = lhs
        self.rhs = rhs

    def __str__(self):
        return formatPartial(self.lhs) + '  -->  ' + formatPartial(self.rhs)

    def getXML(self):
        In = _E.In()
        it = xmlPartial(self.lhs)
        In.extend( it )

        Out = _E.Out()
        Out.extend( xmlPartial(self.rhs) )
        return _E.Eq( In, '  -->  ', Out )


class Interaction:
    def __init__(self, **kwargs):
        self.Eq = Equation( **kwargs )
        self.data = []

    def getXML(self):
        i = _E.Interaction( self.Eq.getXML() )
        dataNodes = [ d.getXML() for d in self.data ]
        i.extend( dataNodes )
        return i


class VSSInteraction(Interaction):
    def __init__(self, reference, vss_param_inv, **kwargs):
        Interaction.__init__(self, **kwargs)
        self.reference = reference
        self.vss_param_inv = vss_param_inv

    def getXML(self):
        x = Interaction.getXML(self)
        x.set('model', 'vss_elastic')
        x.extend( [_E.vss_param_inv( str(self.vss_param_inv) )] )
        if self.reference:
          x.extend( [self.reference.getXML()] )
        return x

