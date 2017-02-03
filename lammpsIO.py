import numpy as np


class StageInstructions():
    ''' Put these in stageList of LAMMPSInput '''


class LAMMPSInput():
    ''' parentStructure is the input
        stageList contains instructions for each step
        Common variables:
            T #temperature
            N #number of steps (production) '''

    def __init__(self,
#                 parentStructure    = None,
#                 fileSuffix         = 'Null',
#                 title              = 'Null',
#                 stageList          = [],
#                 uniqueSpeciesOrder = None,
#                 units              = 'metal',
#                 lammpsVariables    = None,
#                 potentials         = None,
                 **kwargs
                 ):

        # defaults
        # this code is being reviewed 160117 - why should input have suffix- easier to write multiple types file from one class?
        self.fileSuffix         = 'Null'
        self.title              = 'Null'
        self.parentStructure    = None
        self.stageList          = []
        self.uniqueSpeciesOrder = None
        self.units              = 'metal'
        self.lammpsVariables    = {}
        self.moleculeTypes      = None

        for key in kwargs:
            setattr(self, key, kwargs[key])

    def defineUniqueSpeciesOrder(self):
        self.uniqueSpeciesOrder = self.parentStructure.uniqueSpecies(['element', 'core'])
        return self.uniqueSpeciesOrder

    def dataWrite(self, fileName=None):

        if not fileName:
            fileName = 'data.' + fileSuffix

        with open(fileName, 'w') as outf:
            outf.write(self.dataString)

    def numberBonds(self):
        ''' For now, just return the number of core-shell bonds in the parentStructure '''
        return len([x for x in self.parentStructure.speciesList if x.core[:4].lower() == 'shel'])

    def numberAngles(self):
        ''' For now 0 '''
        return 0

    def numberDihedrals(self):
        ''' For now 0 '''
        return 0

    def unitCellDefinition(self):
        # is this a 'right angled' cell, with conventional setting (diagonal)?
        if np.array_equal(np.dot(self.parentStructure.unitCell.vectors, np.identity(3)),
                          self.parentStructure.unitCell.vectors):
            return "0.0  %s  xlo xhi\n0.0  %s  ylo yhi\n0.0  %s  zlo zhi\n"%(self.parentStructure.unitCell.vectors[0,0], 
                                                                             self.parentStructure.unitCell.vectors[1,1], 
                                                                             self.parentStructure.unitCell.vectors[2,2])
        else:
            raise Exception('Implement triclinic cells')

    def labelElementWRTMolecules(self):
        ''' 0 if not a 'molecule' and 1,2,3 etc otherwise '''
        moleculeCounter = 1
        moleculeList    = []
        for i in xrange(len(self.uniqueSpeciesOrder)):
            if self.moleculeTypes is not None and self.uniqueSpeciesOrder[i].element in self.moleculeTypes:
                moleculeList.append((self.uniqueSpeciesOrder[i].element, moleculeCounter))
            else:
                moleculeList.append((self.uniqueSpeciesOrder[i].element, 0))
        return dict(moleculeList)

    def bondInformation(self, cut= 1.e-1):
        ''' String, list of bonds - only works wrt shell model so far (no actual molecules) 
            Use cartesian coordinates '''

        # Ignore if no shells
        if self.moleculeTypes is None:
            return ""

        outString = "\n\nBonds\n\n"

        moleculeCounter = 1
        for i, at1 in enumerate(self.parentStructure.speciesList):

            if at1.element not in self.moleculeTypes or at1.core.lower()[0] != 'c':
                continue

            for j, at2 in enumerate(self.parentStructure.speciesList):

                if at2.core.lower()[0] != 's' or at2.element != at1.element:
                    continue

                if np.linalg.norm(at1.cartCoord - at2.cartCoord) < cut:
                    outString += "%s  %s  %s  %s\n"%(moleculeCounter,
                                                     self.elementToMoleculeDict[at1.element],
                                                     i+1,
                                                     j+1)
                    moleculeCounter += 1
                    continue

        # make sure that all shells have been found
        assert(len([x for x in self.parentStructure.speciesList if x.core.lower()[0] == 's']) == moleculeCounter - 1)
        return outString

    def dataString(self, blurb='data file for LAMMPS'):
        ''' String for the data file '''

        from setTools import sameElementByAttributes, subsetByAttributes

        outString  = blurb + "\n\n"
        outString += str(len(self.parentStructure.speciesList)) + "   atoms\n"
        outString += str(self.numberBonds()) + "   bonds\n"
        outString += str(self.numberAngles()) + "   angles\n"
        outString += str(self.numberDihedrals()) + "   dihedrals\n\n"

        if not self.uniqueSpeciesOrder:
            self.defineUniqueSpeciesOrder()
        self.atomDict = dict([(x[1].element + x[1].core, x[0] + 1) for x in enumerate(self.uniqueSpeciesOrder)])

        # detect shells, and declare as molecule
        if any([x.core.lower()[0] == 's' for x in self.uniqueSpeciesOrder]):
            self.moleculeTypes = [x.element for x in self.uniqueSpeciesOrder if x.core.lower()[0] == 's']
        self.elementToMoleculeDict = self.labelElementWRTMolecules()

        outString += str(len(self.uniqueSpeciesOrder)) + "   atom types\n"
        # for now just stick to simple atoms and shells- ignore molecules!
        outString += str(len([x for x in self.uniqueSpeciesOrder if x.core[0].lower() == 's'])) + "   bond types\n" 
        outString += "0   angle types\n"
        outString += "0   dihedral types\n\n"

        outString += self.unitCellDefinition()
        outString += "\n"
#        outString += "UNIT CELL\n\n"

        assert(self.uniqueSpeciesOrder[0].mass > 0.)
        outString += "Masses\n\n"
        for i, at in enumerate(self.uniqueSpeciesOrder):
            outString += str(i + 1) + " " + str(at.mass) + "\n"
        outString += "\n"

        outString += "Atoms\n\n"
        assert(self.parentStructure.speciesList[0].charge is not None)
        # atom number, molecule type (0 if not), atom index, charge, cart coords
        for i, at in enumerate(self.parentStructure.speciesList):
            outString += "%s   %s   %s   %s  %s  %s  %s\n" %(i+1,
                                                             self.elementToMoleculeDict[at.element],
                                                             self.atomDict[at.element + at.core],
                                                             at.charge,
                                                             at.cartCoord[0], at.cartCoord[1], at.cartCoord[2])

        outString += self.bondInformation()

        return outString

    def vdwInformationFIELD(self, attributeList = ['element', 'core']):
        ''' modified from DLPOLY input (hence odd name)- make general if possible '''
        from coreClasses import VBuckingham
        from setTools import sameElementByAttributes

        # only include the potentials which are required for this structure
        uniqueSpeciesList = [x for x in self.parentStructure.uniqueSpecies(['element', 'core'])]
        if not any([x.core[0].lower() == 's' for x in uniqueSpeciesList]):
            attributeList = ['element']
        print attributeList
        twoBodyPots   = [x for x in self.potentials if isinstance(x, VBuckingham)
                                                    and any([sameElementByAttributes(x.species1, y, attributeList) for y in uniqueSpeciesList])
                                                    and any([sameElementByAttributes(x.species2, y, attributeList) for y in uniqueSpeciesList])]
#        print [(x.element, x.core) for x in uniqueSpeciesList], [(x.species1.element, x.species1.core) for x in self.potentials], twoBodyPots, 'e'

#        outString = "vdW %s\n"%len(twoBodyPots)
        outString = "pair_coeff  *  *  0.0  1.0  0.0\n"
        for p in twoBodyPots:
            indices = sorted([self.atomDict[p.species1.element + p.species1.core],
                              self.atomDict[p.species2.element + p.species2.core]])
            outString += "%s  %s  %s  %s  %s\n"%(indices[0],
                                                 indices[1],
                                                 p.A,
                                                 p.rho,
                                                 p.C6)
        return outString


    def inputString(self):
        outString = "# custom LAMMPS script\n\n"
        
        outString += "units %s\n\n" % self.units

        outString += "\n".join(["variable %s equal %s" %(x[0], x[1]) for x in self.lammpsVariables.items()])
        outString += "\n"

        outString += "#potentials\n"
        outString += self.vdwInformationFIELD()

        return outString


class LAMMPSOutput():

    def __init__(self):
        pass

    @staticmethod
    def structureFromDataFile(inFile,
                              speciesDict   = None, 
                              maintainOrder = False):

        ''' If speciesDict is None, work this all out from file --- do this later if needed 
            Orthorhombic lattice
            Cartesian coordinates '''

        from coreClasses import Structure, Species, UnitCell
        import copy

        with open(inFile, 'r') as inFi:
            inData = inFi.read()

        #assume square box at this stage
        print 'Assuming orthorhombic lattice'
        bounds = []
        
        unitCellData = copy.deepcopy(inData)
        for l in unitCellData.split('\n'):
            if any([x in l for x in ['xlo', 'ylo', 'zlo']]) :
                bounds.append(map(float, l.split()[:2]))

        assert(len(bounds) == 3)

        atomList = []
        # For now assume Atoms is last bit, and just take this
        atomData = copy.deepcopy(inData)

        for l in atomData.split('Atoms')[-1].split('Velocities')[0].split("\n"):

            ls = l.split(' ')
            try:
                species = copy.deepcopy(speciesDict[int(ls[2])])
                species.cartCoord = np.array(ls[4:7], dtype='float64')
                atomList.append((int(ls[0]), species))
            except:
                pass

        #atomList -> just list of species, with order maintained or not
        if maintainOrder:
            atomList = [x[1] for x in atomList]
        else:
            atomList = [x[1] for x in sorted(atomList, key=lambda y:y[0])]
            

        return Structure(unitCell = UnitCell(vectors = np.diag([max(x) - min(x) for x in bounds])),
                         speciesList = atomList)

    @staticmethod
    def logDataToDict(dataFile, cleanLogFirst=True):
        ''' If log file is cleaned beforehand, no need to cut out the useful bit '''
        with open(dataFile, 'r') as inFile:
            lines = inFile.readlines()
            if cleanLogFirst:
                usefulData  = False
                usefulLines = []
                header      = None
                for l in lines:
                    if l.startswith('Step Temp'):
                        usefulData = True
                        header     = l
                    elif l.startswith('Loop'): #better - anything not a number
                        usefulData = False
                    elif usefulData:
                        usefulLines.append(l)
                usefulLines = [header] + usefulLines
            else:
                usefulLines = lines

            keys = usefulLines[0].split()
            data = [ [] for _ in xrange(len(keys))]
            for l in usefulLines[1:]:
                for ix, x in enumerate(l.split()):
                    data[ix].append(x)

        return dict(zip(keys, data))

