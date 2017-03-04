import numpy as np

class XSF():
    def __init__(self,
                 structure = None):
        self.structure = structure

    def structureString(self):
        ''' Assume P1 cell - this string can express a structure '''
        from ase.data import atomic_numbers

        outString  = "CRYSTAL\n"

        outString += "PRIMVEC\n"
        outString += "\n".join([" ".join(map(str, self.structure.unitCell.vectors[i])) for i in xrange(3)])

        outString += "\nPRIMCOORD\n"
        outString += "%s %s\n"%(len([x for x in self.structure.speciesList if x.core.lower()[0] == 'c']), 1)

        #coordinates are in Angstoms
        if self.structure.speciesList[0].cartCoord is None:
            self.structure.setCartCoord()

#        outString += "\n".join([" ".join(map(str, [atomic_numbers[x.element]] + list(x.cartCoord))) for x in self.structure.speciesList if x.core.lower() == 'c'])
        for x in self.structure.speciesList:
            outString += "%s %s %s %s\n"%(atomic_numbers[x.element], x.cartCoord[0], x.cartCoord[1], x.cartCoord[2])

        return outString

    def dataGrid3DString(self, 
                         dataGrid   = None,
                         gridPoints = np.zeros(3),
                         centralPt  = np.zeros(3)):

        ''' N.B. the grid should be in Fortran style (column major) 
            Assume Gamma centred (this is what the 0,0,0 means??) 
            grid points is number along each axis e.g. (5,4,20) 
            dataGrid is array(np.prod(gridPoints)) '''

        if dataGrid is not None:
            self.dataGrid3D = dataGrid

        outString  = "BEGIN_BLOCK_DATAGRID_3D\n"
        outString += "Pointless comment\n"
        outString += "BEGIN_DATAGRID_3D_this_is_3Dgrid\n"

        # grid info (npts each axis, centering, vectors)
        outString += " ".join([str(x) for x in gridPoints]) + "\n"
        outString += " ".join([str(x) for x in centralPt]) + "\n"
        outString += "\n".join([" ".join(map(str, self.structure.unitCell.vectors[i])) for i in xrange(3)])

        # just write the whole data out (in column major)
        outString += "\n".join([str(x) for x in self.dataGrid3D])

        outString += "END_DATAGRID_3D\n"
        outString += "END_BLOCK_DATAGRID_3D"

        return outString

    def xsfString(self, **kwargs):
#                  dataGrid   = None,
#                  gridPoints = np.zeros(3),
#                  centralPt  = np.zeros(3)):

        return self.structureString() + "\n" + self.dataGrid3DString(**kwargs) #dataGrid   = dataGrid,
                                                                     #gridPoints = gridPoints,
#                         centralPt  = np.zeros(3)):
