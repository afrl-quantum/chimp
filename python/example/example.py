#!/usr/bin/env python

##### FIX PYTHON PATH FOR THIS EXAMPLE #####
import sys
sys.path.append( '../' )
##### FIX PYTHON PATH FOR THIS EXAMPLE #####

from numpy import *
import chimp
from chimp.interaction import cross_section, Interaction
from lxml.builder import ET, E

interactions = E.Interactions()


######
R = chimp.Reference("""
      This data was cooked right out of my head!
        --Spencer
      """, bibtex="""
@Article{TheKey,
    author = {Bob T. Builder and Barbie Queue },
    title = {Stuff that is important},
    journal = {Chem. Rev. A},
    year = {1995},
    volume = {51},
    number = {},
    month = {},
    pages = {834-99},
    keywords = {PACS. 30.30.-s, 99.80.Pj, 08.20.Jp }
}     """)
rawdata = loadtxt('dummy-data.txt')
D = cross_section.DATA(reference=R, xscale='eV', yscale='nm^2', data=rawdata)
I = Interaction( lhs={'e^-':1,'H2':1}, rhs={'e^-':2,'H2^+':1} )

I.data.append( D )

interactions.append( I.getXML() )
#####


######
R = chimp.Reference("""
      This data was cooked right out of my head!
        --Spencer
      """, bibtex="""
@Article{TheDuke,
    author = {John Wayne},
    title = {True grit is not made from grits},
    journal = {Holywood. Rev. B},
    year = {1960},
    volume = {3},
    number = {},
    month = {},
    pages = {1-10},
}     """)
rawdata = loadtxt('dummy-data2.txt')
D = cross_section.DATA( reference=R,
                        xscale='1/K_B',
                        yscale='1e-12*barns',
                        data=rawdata )
I = Interaction( lhs={'Hg':1,'e^-':1}, rhs={'Hg^+':1,'e^-':2} )

I.data.append( D )

interactions.append( I.getXML() )
#####

print ET.tostring( E.Olson( interactions ), pretty_print=True )
