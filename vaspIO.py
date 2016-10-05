import numpy as np
import fortranformat as ff

fortranWriter= ff.FortranRecordWriter('(3f15.10)')

class vaspInput():
    def __init__(self,
                 parentStructure    = None,
#                 aseStructure       = None,
                 fileName           = None,
                 title              = 'vasp input',
                 unitCell           = None,
                 natms              = None,
                 units              = 'eV',
                 temperature        = None,
                 integrator         = None,
                 moleculeOrder      = None,
                 ):
        self.fileName          = fileName
        self.title             = title
        self.natms             = natms
        self.parentStructure   = parentStructure
#        self.aseStructure      = aseStructure
        self.moleculeOrder     = moleculeOrder
        self.units             = units
        self.temperature       = temperature
        self.integrator        = integrator
        if self.parentStructure:
            self.unitCell      = self.parentStructure.unitCell
            self.speciesList   = self.parentStructure.speciesList

    def alphabeticalMoleculeOrder(self):
        self.moleculeOrder = sorted(set([x.element for x in self.speciesList]))

    def requiredOrderedAtoms(self):
        for ele in self.moleculeOrder:
            for at in self.speciesList:
                if at.element == ele and at.core.lower() == 'core':
                    yield at

    def atomicPositionsString(self, fractional):
        ''' if not fractional, use cartesian  '''

        if not self.moleculeOrder:
            self.alphabeticalMoleculeOrder()

        if fractional:
            positionsList = [x.fracCoord for x in self.requiredOrderedAtoms()]

        return "\n".join([fortranWriter.write(x) for x in positionsList])

    def numberSpecies(self):
        ''' Just numbers of atoms in self.moleculeOrder '''
        for m in self.moleculeOrder:
            yield len([x for x in self.speciesList if x.element == m and x.core.lower() == 'core'])

#    @classmethod
    def generatePOSCARString(self,
                             scalingFactor = 1.,
                             direct        = 'direct',
                             speciesOrder  = None):

        #order is important to VASP- must match POTCAR
        if not self.moleculeOrder:
            self.alphabeticalMoleculeOrder()

        outString  = self.title + "\n"
        outString += str(scalingFactor) + "\n"

        outString += "\n".join([fortranWriter.write(x) for x in self.unitCell.vectors]) + "\n"

        outString += " ".join(map(str, [_ for _ in self.numberSpecies()])) + "\n"

        if direct.lower()[0] == 'd':
            outString += "Direct\n"
            outString += self.atomicPositionsString(True)

        return outString

class VaspXML():
    def __init__(self):
        pass

    @classmethod
    def createStructure(cls,
#                        **kwargs):
                        xmlString = None,
                        structureIndex     = None):
        ''' This runs through all the calculations run to find a particular index, and makes a structure corresponding to the input for this '''
        import xmltodict
        from coreClasses import Structure, Species, UnitCell
        doc = xmltodict.parse(xmlString)
        nAtoms = int(doc['modeling']['atominfo']['atoms'])
#        nTypes = int(doc['modeling']['atominfo']['types'])
        positions = np.array([str(x).split() for x in doc['modeling']['calculation'][structureIndex]['structure']['varray']['v']], dtype='float64')
        elements = map(str, [doc['modeling']['atominfo']['array'][0]['set']['rc'][i]['c'][0] for i in xrange(nAtoms)])

        return Structure(unitCell    = UnitCell(vectors =  np.array([str(x).split() for x in doc['modeling']['calculation'][structureIndex]['structure']['crystal']['varray'][0]['v']], dtype='float64')),
                         speciesList = [Species(element = elements[i], fracCoord = positions[i]) for i in xrange(nAtoms)])

    @classmethod
    def initialStructure(cls,
                         xmlString):
         return cls.createStructure(xmlString      = xmlString,
                                    structureIndex = 0)

    @classmethod
    def lastAttemptedStructure(cls,
                               xmlString):
         return cls.createStructure(xmlString      = xmlString,
                                    structureIndex = -1)

    @classmethod
    def finalStructure(cls,
                         xmlString):
        ''' xml file has two structures outside of the calculation list'''
        import xmltodict
        from coreClasses import Structure, Species, UnitCell
        doc = xmltodict.parse(xmlString)
        nAtoms = int(doc['modeling']['atominfo']['atoms'])
        positions = np.array([str(x).split() for x in doc['modeling']['structure'][1]['varray']['v']], dtype='float64')
        elements = map(str, [doc['modeling']['atominfo']['array'][0]['set']['rc'][i]['c'][0] for i in xrange(nAtoms)])

        return Structure(unitCell    = UnitCell(vectors =  np.array([str(x).split() for x in doc['modeling']['structure'][1]['crystal']['varray'][0]['v']], dtype='float64')),
                         speciesList = [Species(element = elements[i], fracCoord = positions[i]) for i in xrange(nAtoms)])
#         return cls.createStructure(xmlString      = xmlString,
#                                    structureIndex = -1)
