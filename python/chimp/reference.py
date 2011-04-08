from lxml.builder import E as _E

class Reference:
  def __init__(self, text, bibtex=None):
    self.text = text
    self.bibtex = bibtex

  def getXML(self):
    x = _E.reference( self.text )

    if self.bibtex:
      x.extend( [ _E.bibtex(self.bibtex) ] )

    return x

