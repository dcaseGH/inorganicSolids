#Core classes to act as interfaces between programs and perform simple manipulations
import numpy as np

class XYZFile:
    ''' Some of these may be better being staticmethods and some as instances 
        Beware of things where there may be lines with a character and then unit cells '''
    def __init__(self,
                 fileName          = None,
                 stringData        = None,
                 maxStructures     = None, 
                 linesPerStructure = None):
#                 nStructures       = None):

        self.fileName          = fileName
        self.stringData        = stringData
        self.maxStructures     = maxStructures
        self.linesPerStructure = linesPerStructure

    @staticmethod
    def returnString(inString, swapDict = {}):
        ''' Quick way to rewrite files with dummy labels (use swapDict) '''
        outLines = []
        for line in inString.split("\n"):
            if len(line) == 0:
                continue
            splitLine = line.split()

            if splitLine[0] in swapDict.keys():
                splitLine[0] = swapDict[splitLine[0]]

            outLines.append(" ".join(splitLine))
        return "\n".join(outLines)


    def setLinesPerStructure(self):
        ''' The offset will depend on whether there are unit cell vectors etc '''
        counter = 0
        header = True
        for line in open(self.fileName, 'r'):
            if header and self.standardSpeciesLine(line):
                header = False
            elif not header and not self.standardSpeciesLine(line):
                break
            counter += 1
        self.linesPerStructure = counter
        return self.linesPerStructure

    @staticmethod
    def standardSpeciesLine(testLine):
        ''' A standard species line has 4 things [element] [x] [y] [z]
            The element must begin with a character- the x,y,z are numbers '''
        splitLine = testLine.split()
        if len(splitLine) != 4:
            return False
        # Deal with cases where atom has a numerical index at another time
#        print 'MASSIVE HACK IN CORECLASSES.XYZ';return True
        if splitLine[0][0].isalpha() and all([x.replace('.', '', 1).replace('e', '', 1).replace('-', '', 1).replace('E', '', 1).isdigit() for x in splitLine[1:]]):
            return True
        return False

    @staticmethod
    def nStructures(fileName,
                    linesPerStructure):
        ''' calculate linesPerStructure manually or with quick subroutine 
            watch for blank lines '''
        num_lines = float( sum(1 for line in open(fileName, 'r')) )
        nStructures = num_lines / float(linesPerStructure)
        try:
            assert(nStructures.is_integer())
        except:
            raise Exception('Error XYZFile.nStructures')
        return int(nStructures)

    @staticmethod
    def returnXYZStrings(fileName,
                         linesPerStructure,
                         maxNStructures    = None,
                         selectList        = None,
                         returnIndex       = False):

        ''' Not fully tested or commented '''
        counter           = 0
        yieldCounter      = 0
        structureCounter  = 0
        lines   = []

        # define max structures if selectList to make this quicker
        if maxNStructures is None and selectList is not None:
            maxNStructures = len(selectList)

        for line in open(fileName, 'r'):

            counter += 1
            lines.append(line)

            if len(lines) == linesPerStructure:

                #do not allow blank lines at end of string
                lines[-1].replace('\n', '')
                if not selectList or structureCounter in selectList:
                    if returnIndex:
                        yield (yieldCounter, "".join(lines))
                    else:
                        yield "".join(lines)
                    yieldCounter += 1 

                lines = []
                structureCounter += 1

            # test this!!!! GeneratorExit ???
            if maxNStructures and yieldCounter >= maxNStructures:
                raise StopIteration


    @staticmethod
    def xyzStringToCellVectors(xyzString):
        ''' Cell info on second line-- assume orthorhombic for now- easy to add if need  '''
        cellInfo = xyzString.split("\n")[1].split()

        #change this if need to go beyond orthorhombic - assume second line is only this cell info
#        assert(len(cellInfo) == 3)
        if len(cellInfo) == 3:
            return np.diag(map(float, cellInfo))
        elif len(cellInfo) == 6:
            print 'warning, better to calculate angles etc and turn to cell'
            return np.array(map(float, [cellInfo[0], cellInfo[3], cellInfo[4], 0., cellInfo[1], cellInfo[5], 0., 0., cellInfo[2]])).reshape((3,3))

    @staticmethod
    def xyzStringToSpeciesList(xyzString, speciesDict=None, offset = 0):
        ''' Clean all this up at some point 
            Atom Dict is in case of labelling atoms 1,2,3 not Fe, S, O etc - also can add other info if with in the atomDict 
            atm the species dict should have str for keys '''
        from copy import deepcopy
        speciesList = []
        for l in xyzString.split("\n")[offset:]:
            if XYZFile.standardSpeciesLine(l):
                if speciesDict:
#                    tempAtom = deepcopy(speciesDict[int(l.split()[0])])
                    tempAtom = deepcopy(speciesDict[l.split()[0]])
                    tempAtom.cartCoord = np.array(l.split()[1:], dtype='float64')
                else:
                    tempAtom = Species(element   = l.split()[0],
                                       cartCoord = np.array(l.split()[1:], dtype='float64') )
                speciesList.append(tempAtom)
        return speciesList

    ####### THIS IS NOT READY !!!!!!!!!!!!!!!#########
    @staticmethod
    def xyzList(self, inFile, returnObjects = 'xyzString'):#, containsUnitCell = True):
        ''' Read file and break into list of xyzStrings, PMGStructures, or Structures 
            Either give a string, or dummy instance of either structure type-- do ASE if ever needed'''

        print "XYZ may or may not have unit cell - should work it out automatically"
        structureList = []

        topBufferLines = 2
        counter        = 0
        lines          = []
        for line in open(historyFileName):
            counter += 1
            if counter <= topBufferLines:
                nAtoms = line.split()[-1]
                continue
            exit()
            timestepBufferLines = int(nAtoms) * 2 + 4
            lines.append(line)

            if len(lines) == timestepBufferLines:

                # this should just be the quickest thing to break up the file
                if isInstance(returnObjects, str):
                    structureList.append("\n".join(lines))
                    if maxStructures and len(structureList)>= maxStructures:
                        return structureList

                    lines = []
                    continue

                latticeVectors = np.array([x.split() for x in lines[1:4]], dtype='float64')
                species, coords = [], []

                offset = 4
                for ix, x in enumerate(lines[offset:]):
                    if ix%2==0:

                        if x == '' or (ignoreShells and 'shl' in x.split()[0]):
                            continue
                        elif selectOnlySpecie and x.split()[0].replace('_', '') != selectOnlySpecie:
                            continue
                        else:
                            species.append(x.split()[0].replace('_', ''))
                            coords.append( np.array(lines[ix+1+offset].split(), dtype='float64') )

                if selectOnlySpecie and addOrigin:
                    species.append('X')
                    coords.append(np.zeros(3))

                if isInstance(returnObject, PMGS ):
                    structureList.append( PMGS(PMGL(latticeVectors),
                                               species,
                                               coords,
                                               coords_are_cartesian = True)
                                          )
                elif isInstance(returnObject, Structure):
                    structureList.append(Structure(unitCell    = UnitCell(vectors = latticeVectors),
                                                   speciesList = [Species(element   = species[i],
                                                                          cartCoord = coords[i]) for i in xrange(species)])
                                         )
                    
                if maxStructures and len(structureList)>= maxStructures:
                    return structureList

                lines = []

        return structureList



class Species:

    def __init__(self, 
                 label          = None,
                 potentialLabel = None,
                 element        = None,
                 mass           = None,
                 charge         = None,
                 fracCoord      = None,
                 cartCoord      = None,
                 cartVelocity   = None,
                 cartForce      = None,
                 core           = 'core'
                 ):

        self.label          = label
        self.potentialLabel = potentialLabel
        self.element        = element
        self.mass           = mass
        self.charge         = charge
        self.fracCoord      = fracCoord
        self.cartCoord      = cartCoord
        self.core           = core
        self.cartVelocity   = cartVelocity
        self.cartForce      = cartForce

    def toASEAtom(self):
        from ase import Atom
        return Atom(str(self.element), self.cartCoord)

    @classmethod
    def initFromASEAtom(cls, aseAtom):
        return cls(cartCoord = aseAtom.position, 
                    element  =  aseAtom.symbol)

    @classmethod
    def initFromPMGSite(cls, pmgSite):
        ''' Strictly PMG has _sites and so on ... 
            _species sits on _site '''

        return cls(cartCoord = pmgSite._coords,
                   fracCoord = pmgSite._fcoords,
                   element   = pmgSite._species._data.keys()[0].symbol) 

    @staticmethod
    def elementFromPMGSpecies(pmgSpecies):
        ''' Should return element string from an instance of a PMG _species '''
        return pmgSpecies._data.keys()[0].symbol

    def dlpolyLabelStandard(self):
        # for the moment try to use just this standard labelling scheme
        outString = self.element
        if len(self.element) < 2:
            outString += "_"
        if self.core[:4] == 'shel':
            outString += "-shl"
        return outString

    def setThermalVelocity(self, temp):
        ''' Velocity is normal with mu=0, sigma = (kT/m)**0.5
            because of shells, atm the mass is standard for core and v is 0 for shell '''

        if self.core[:4] == 'shel':
            self.cartVelocity = np.zeros(3)
            return self.cartVelocity

        import random
        from scipy.stats import norm
        from ase.data import atomic_masses, chemical_symbols

        # per mol
        sigma = (6.0221367e23 * 1.380658e-23 * temp / atomic_masses[chemical_symbols.index(self.element)]) ** 0.5
        self.cartVelocity = np.array([norm.ppf(random.random(), scale=sigma) for _ in xrange(3)])
        return self.cartVelocity

    def atomicValenceElectrons(self):
         from hardcode import atomicValenceElectrons
         return atomicValenceElectrons[self.element]

    def defaultCharge(self):
        ''' Use this v. sparingly 
            Do not assume charges are standard- always specify except for quick plots etc '''
        from hardcode import defaultCharges
        return defaultCharges[self.element]

class Defect():#Species):
    # Do super classes later- python 2 vs 3 issues (and I did it wrong the first time)
    def __init__(self,
                 defectType = None,
                 species    = None):
#        super(Species, self).__init__(species)
        self.species    = species
        self.defectType = defectType

    def charge(self):
        # minus the charge....
        return None


class UnitCell:
    # add vectors and stuff as needed
    def __init__(self,
                 angles  = None,
                 lengths = None,
                 vectors = None):

        self.angles     = angles
        self.lengths    = lengths
        self.vectors    = vectors
        self.invVectors = None

        if self.vectors is None and self.angles is not None and self.lengths is not None:
            self.vectors = self.calculateVectors(self.lengths, self.angles)

        if self.vectors is not None and self.angles is None and self.lengths is None:
            self.lengths = np.array(map(np.linalg.norm, [self.vectors[i] for i in xrange(3)]))
            self.angles  = (180. / np.pi) * np.array(map(np.arccos,[np.dot(self.vectors[1], self.vectors[2]) / (self.lengths[1] * self.lengths[2]),
                                                                    np.dot(self.vectors[0], self.vectors[2]) / (self.lengths[0] * self.lengths[2]),
                                                                    np.dot(self.vectors[0], self.vectors[1]) / (self.lengths[0] * self.lengths[1])]))

    @staticmethod
    def generalLengthsAnglesMatrix(matrix, anglesDegrees = True):
        ''' If have matrix in unusual form, can get angles and lengths
            Useful e.g. with LAMMPS cells 
            matrix is numpy array with indices (fractional, cartesian) 
            '''

        if anglesDegrees:
            convFactor = 180./np.pi
        else:
            convFactor = 1.

        lengths = np.array([np.linalg.norm(matrix[x, :]) for x in xrange(3)])
        angles  = convFactor * np.array(map(np.arccos, [np.dot(matrix[1], matrix[2]) / (lengths[1] * lengths[2]),
                                                        np.dot(matrix[0], matrix[2]) / (lengths[0] * lengths[2]),
                                                        np.dot(matrix[0], matrix[1]) / (lengths[0] * lengths[1])]))

        return lengths, angles

    def calculateVectors(self, lengths, angles):
        """
        stolen from pymatgen- will possibly use their code in future
        """
        a = lengths[0]
        b = lengths[1]
        c = lengths[2]
        alpha_r = angles[0] * np.pi / 180.
        beta_r  = angles[1] * np.pi / 180.
        gamma_r = angles[2] * np.pi / 180.
        val = (np.cos(alpha_r) * np.cos(beta_r) - np.cos(gamma_r))\
            / (np.sin(alpha_r) * np.sin(beta_r))
        #Sometimes rounding errors result in values slightly > 1.
        val = min([max([-1., val]), 1.])
        gamma_star = np.arccos(val)
        vector_a = [a * np.sin(beta_r), 0.0, a * np.cos(beta_r)]
        vector_b = [-b * np.sin(alpha_r) * np.cos(gamma_star),
                    b * np.sin(alpha_r) * np.sin(gamma_star),
                    b * np.cos(alpha_r)]
        vector_c = [0.0, 0.0, float(c)]
        self.vectors = np.array([vector_a, vector_b, vector_c])
        return self.vectors

    def setInvVectors(self):
        self.invVectors = np.linalg.inv(self.vectors)

    def createGrid(self, targetSeparation = 1., includeExtraPoint = True, returnMGrid = False, returnInfo = False, columnMajor = False):
        ''' Return a grid with roughly targetSeparation along axes (use meshgrid for certain things) 
            Normal use is (Npts,3) numpy array
            Alternatively, return 3 (nx,ny,nz) arrays- see mgrid vectorization '''

        if includeExtraPoint:
            extraPt = 1
        else:
            extraPt = 0

        nPts = np.floor(self.lengths / targetSeparation)
        nPtsI = np.array(np.floor(self.lengths / targetSeparation), dtype=int)

        if returnMGrid:
            print 'not implemented yet';exit()
            if includeExtraPoint:
                return np.mgrid[0:1:nPtsI[0]*1j,
                                0:1:nPtsI[1]*1j,
                                0:1:nPtsI[2]*1j]
            #mgrid must be orthogonal??
#                return [np.dot(x, self.vectors) for x in np.mgrid[0:1:nPtsI[0]*1j,
#                                                                  0:1:nPtsI[1]*1j,
#                                                                  0:1:nPtsI[2]*1j]]

        #delete this
        if returnInfo:
            print "Grid = (%s, %s, %s)"%(nPtsI[0] + extraPt,
                                         nPtsI[1] + extraPt,
                                         nPtsI[2] + extraPt)

        if columnMajor:
            outpoints = np.dot([(i,j,k)/nPts for k in xrange(nPtsI[2] + extraPt)
                                             for j in xrange(nPtsI[1] + extraPt)
                                             for i in xrange(nPtsI[0] + extraPt)], self.vectors)
        else:
            outpoints = np.dot([(i,j,k)/nPts for i in xrange(nPtsI[0] + extraPt)
                                             for j in xrange(nPtsI[1] + extraPt)
                                             for k in xrange(nPtsI[2] + extraPt)], self.vectors)
           
        if returnInfo:
            return (nPtsI + np.array([extraPt, extraPt, extraPt]), outpoints)

        return outpoints

    @classmethod
    def fromXYZFile(cls, xyzFile, transposeVectors = False):
        if transposeVectors:
            return cls(vectors = XYZFile.xyzStringToCellVectors(open(xyzFile.fileName, 'r').read()).T)
        else:
            return cls(vectors = XYZFile.xyzStringToCellVectors(open(xyzFile.fileName, 'r').read()))

    def vertices(self, fractional = True):
        ''' vertices of the unit cell - property?? '''

        fracVertices = np.array([(i,j,k) for i in xrange(2) for j in xrange(2) for k in xrange(2)])

        if not fractional:
            return np.dot(fracVertices, self.vectors)
        else:
            return fracVertices

    @classmethod
    def randomCell(cls,
                   angles              = np.array([90., 90., 90]), 
                   lengths             = None,
                   angleRandomNumbers  = None,
                   lengthRandomNumbers = None,
                   targetVolume        = None):
        ''' Generate a random cell 
            see organic structure generation (ask Dave) for details '''

        if angleRandomNumbers is not None:
            asin_val = np.arcsin(2.0 * angleRandomNumbers - 1.0 ) / np.pi
            minAngle, delta = 60., 120. - 60.
            _angles = (0.5 + asin_val)*delta + minAngle
        else:
            _angles = angles

        if lengths is not None:
            _lengths = lengths
        else:
            d2r = np.pi/180.
            from scipy.stats import norm
            vStar = (1.0+ 2.0 * np.cos(_angles[0]*d2r) * np.cos(_angles[1]*d2r) * np.cos(_angles[2]*d2r) - np.cos(_angles[0]*d2r)**2 - np.cos(_angles[1]*d2r)**2 - np.cos(_angles[2]*d2r)**2 )**0.5
            _sigma = 2.           
            _lengths = np.array([norm.ppf(x) * _sigma + (targetVolume/vStar)**(1./3.) for x in np.clip(lengthRandomNumbers, 0.01, 0.99)])
#            if rand_vec.unit_cell_lengths[0] > 0.99:
#                temp_lengths [ 0 ] = norm.ppf(0.99) * sd_param * self.length_bounds[0][0] + mean0
#            elif rand_vec.unit_cell_lengths[0]<0.01:
#                temp_lengths [ 0 ] = norm.ppf(0.01) * sd_param * self.length_bounds[0][0] + mean0
#            else:
#                temp_lengths [ 0 ] = norm.ppf(rand_vec.unit_cell_lengths[0]) * sd_param * self.length_bounds[0][0] + mean0            _sigma   = 0.1
#            _lengths = 
#            _lengths = np.random.normal(targetVolume ** (1./3.),
#                                        _sigma,
#                                        3)

        return cls(angles  = np.array(_angles),
                   lengths = np.array(_lengths))

class Potential:
    ''' At the moment, just a holder '''
    # no need to __init__ ???

    @staticmethod
    def adaptRemoveShells(potsIn):
        ''' Input a list of pots, change VBuckingham species.core -> core  '''
        from copy import deepcopy
        outPots = []
        for x in potsIn:
            if isinstance(x, VBuckingham):
                x.species1.core = 'core'
                x.species2.core = 'core'
                outPots.append(x)
        return outPots

class VSpring:
    ''' Assume that spring between species1 core and shel  '''
    def __init__(self,
                 species1   = None,
                 K          = None):
        self.species1 = species1
        self.K        = K

    def stringForm(self):
        if self.species1.potentialLabel:
            l1 = self.species1.potentialLabel
        else:
            l1 = self.species1.element

        return "spring\n%s %s"%(l1, self.K)

class VBuckingham:
    ''' Follow 1.4.3 GULP manual (rho is length (angstrom))
        fitA etc are flags (ignore for now) '''
    def __init__(self,
                 species1   = None,
                 species2   = None,
                 A          = None,
                 rho        = None,
                 C6         = None,
                 name       = 'buck',
                 toBeFitted = False,
                 cutMin     = 0.,
                 cutMax     = 12.,
                 fitA       = 0,
                 fitRho     = 0,
                 fitC6      = 0):

        self.species1 = species1
        self.species2 = species2
        self.A        = A
        self.fitA     = fitA
        self.rho      = rho
        self.fitRho   = fitRho
        self.C6       = C6
        self.fitC     = fitC6
        self.cutMin   = cutMin
        self.cutMax   = cutMax

    def energy(self, r, chargeProduct = None):
        ''' Units are same as A or C6, r is rho^{-1} (what does this mean???)
            Set chargeProduct to 'auto' if want to use defaults (check before using, e.g. Fe 2/3+ and so on) '''
        from hardcode import hartrees2eV, bohr2angstrom

        if type(chargeProduct) == float:
            return self.A * np.exp(-r / self.rho) - self.C6 * r **-6             
        elif type(chargeProduct) == str and chargeProduct.lower() == 'auto':
            return self.species1.defaultCharge() * self.species2.defaultCharge() * hartrees2eV * (r / bohr2angstrom)**-1 +\
                   self.A * np.exp(-r / self.rho) - self.C6 * r **-6 
        else:
            return self.A * np.exp(-r / self.rho) - self.C6 * r **-6 


    def stringForm(self):
        ''' specify potential labels if you want different types of same atom 
            N.B. THIS IS GULP STRING FORM '''
        if self.species1.potentialLabel:
            l1 = self.species1.potentialLabel
        else:
            l1 = self.species1.element

        if self.species2.potentialLabel:
            l2 = self.species2.potentialLabel
        else:
            l2 = self.species2.element
            
        return "buck\n" + " ".join([l1] + [self.species1.core] +
                                   [l2] + [self.species2.core] +
                                   [str(self.A)] + [str(self.rho)] + [str(self.C6)] +
                                   [str(self.cutMin)] + [str(self.cutMax)])

    def dlpolyString(self):
        return " ".join([self.species1.dlpolyLabelStandard(),
                         self.species2.dlpolyLabelStandard(),
                         "buck",
                         str(self.A),
                         str(self.rho),
                         str(self.C6)])

    @staticmethod
    def removeSpeciesFromListBucks(listPots, species, attributeList = ['element']):
        ''' listPots can be non buckingham, but just remove elements which are 
            buckingham and contain species '''

        from setTools import sameElementByAttributes

        return [x for x in listPots if x.__class__.__name__ != 'VBuckingham' 
                                    or (not sameElementByAttributes(x.species1, species, attributeList)
                                    and not sameElementByAttributes(x.species2, species, attributeList))]


    @staticmethod
    def latexTableFromListBucks(listPots, fileName = None):
        ''' From a list of potentials, make a latex table from the Buckingham pots 
            if fileName, write to this file, otherwise return a string '''
        stringOut  = r'\begin{tabular}{lrrr}' + "\n"
        stringOut += r'\toprule' + "\n"
        stringOut += r'Species &       A &       $\rho$ &       C$_{6}$ \\' + "\n"
        stringOut += r'\midrule' + "\n"

        for p in listPots:
            if p.__class__.__name__ == 'VBuckingham':
                a=        " & ".join(map(str, [p.species1.element + "-" + p.species2.element,
                                               p.A,
                                               p.rho,
                                               p.C6]))
                stringOut += a + r' \\' + "\n"

        stringOut += r'\bottomrule' + "\n"
        stringOut += r'\end{tabular}' + "\n"

        if fileName is not None:
            with open(fileName, 'w') as outf:
                outf.write(stringOut)

        return stringOut

    @staticmethod
    def plotListBuckinghams(listPots, nPts = 100, xBounds = [1., 8.], figsize=(10,10),
                            speciesList = None, matchAttributes = ['element']):
        ''' Make a matplotlib output for each VBuckingham '''
        import matplotlib.pyplot as plt
        from setTools import subsetByAttributes

        listBucks = [x for x in listPots if x.__class__.__name__ == 'VBuckingham' and
                                            subsetByAttributes([x.species1], speciesList, matchAttributes) and
                                            subsetByAttributes([x.species2], speciesList, matchAttributes)]

        xPts      = [min(xBounds) + (max(xBounds) - min(xBounds)) * x/nPts for x in xrange(nPts)]
        plt.subplots(figsize=figsize)
        for i, lb in enumerate(listBucks):
            plt.subplot(len(listBucks), 1, i)
            plt.plot(xPts, [lb.energy(x, chargeProduct='auto') for x in xPts])
            plt.title(lb.species1.element + "_" + lb.species2.element)

        plt.show()
#        return listBucks

    def fitToLJ(self,
                fitBounds    = np.array([1., 2.5]),
                fitPoints    = 100,
                initialGuess = np.array([10., 0.])):
        ''' return an LJ potential with fitted parameters '''
        from scipy.optimize import fmin

#        outLJPot = 
        from scipy.optimize import fmin

        testPoints = np.array([np.min(fitBounds) + (np.max(fitBounds) - np.min(fitBounds)) * x / float(fitPoints) for x in xrange(fitPoints)])

        def squareDifference(constants):
            ''' constants are np.array([c12, c6]) '''
            if constants[0] < 0. or constants[1] > 0.:
                penalty=1000.
            else:
                penalty = 0.
            return sum([(self.energy(r) - VLennardJones.r6r12Energy(r, constants[0], constants[1]))**2 for r in testPoints]) + penalty
#            return sum([(self.energy(r) - VLennardJones.r6r12Energy(r, constants[0], constants[1]) / self.energy(r))**2 for r in testPoints])

        print testPoints
        outputMinimization = fmin(squareDifference, initialGuess)
        print self.energy(2.), VLennardJones.r6r12Energy(2., outputMinimization[0], outputMinimization[1])
        print outputMinimization, "c12 c6"
#        print self.energy(1.), VLennardJones.r6r12(1., 
#        outLJPot.c6 = 
        return VLennardJones(species1 = self.species1,
                             species2 = self.species2,
                             c6       = outputMinimization[1],
                             c12      = outputMinimization[0])

class VLennardJones:
    # be clear about conventions- note LAMMPS may be different (and note signs etc)
    def __init__(self,
                 species1 = None,
                 species2 = None,
                 c6       = None,
                 c12      = None):

        self.species1 = species1
        self.species2 = species2
        self.c6       = c6
        self.c12      = c12

    def energy(self, r):
        return self.c12 * r**(-12) - self.c6 * r**(-6)

    @staticmethod
    def r6r12Energy(r, c12, c6):
        return c12 * r**(-12) - c6 * r**(-6)

    def returnEpsilonSigma(self):
        print self.c6, self.c12
        epsilon = self.c6 **2 / (4. * self.c12)
        sigma   = (self.c12 / self.c6)**(1./6.)
        return np.array([epsilon, sigma])

class VThreeBody:
    ''' Follow 1.4.3 GULP manual (rho is inverse length)
        fitA etc are flags (ignore for now) '''
    def __init__(self,
                 species1   = None,
                 species2   = None,
                 species3   = None,
                 K          = None,
                 theta0     = None,
                 cut12      = None,
                 cut13      = None,
                 cut23      = None):

        self.species1 = species1
        self.species2 = species2
        self.species3 = species3
        self.K        = K
        self.theta0   = theta0
        self.cut12    = cut12
        self.cut13    = cut13
        self.cut23    = cut23

    def stringForm(self):
        ''' specify potential labels if you want different types of same atom '''
        if self.species1.potentialLabel:
            l1 = self.species1.potentialLabel
        else:
            l1 = self.species1.element

        if self.species2.potentialLabel:
            l2 = self.species2.potentialLabel
        else:
            l2 = self.species2.element

        if self.species3.potentialLabel:
            l3 = self.species3.potentialLabel
        else:
            l3 = self.species3.element
            
        return "three\n" + " ".join([l1] + [self.species1.core] +
                                    [l2] + [self.species2.core] +
                                    [l3] + [self.species3.core] +
                                    [str(self.K)] + [str(self.theta0)] +
                                    [str(self.cut12)] + [str(self.cut13)] + [str(self.cut23)])


# buck coulomb.. inherit buckingham and coulomb

class SymmetryGroup:
    ''' Most likely a space group (point groups are for molecules)
        elementList is going to be in pmg form (so use element.affine_matrix)
        self.elementList[0].__class__ = <class 'pymatgen.core.operations.SymmOp'> '''
    def __init__(self,
                 labelHM     = None,
                 number      = None,
                 elementList = []):
        self.labelHM     = labelHM
        self.number      = number
        self.elementList = elementList

    @classmethod
    def fromCif(cls, fileName):
        ''' wrapper- could be useful to use pmg stuff '''

        from pymatgen.core import Structure as PMGS
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        analyzer = SpacegroupAnalyzer(PMGS.from_file(fileName))
        return cls(labelHM     = analyzer.get_space_group_symbol(),
                   number      = analyzer.get_space_group_number(),
                   elementList = analyzer.get_space_group())

    def generateUniquePoints(self, testPts, cut = 1.e-5, boundUnitCell = False):
        outList = []

        for t in testPts:
            for e in self.elementList:
#            for x in e.operate_multi(testPts):

                x = e.operate(t)
                if boundUnitCell:
                    testPoint = x - np.floor(x)
                else:
                    testPoint = x

                if not any([np.linalg.norm(testPoint - y) < cut for y in outList]):
                    outList.append(testPoint)

        return np.array(outList)

class Structure:
    def __init__(self,
                 unitCell      = None,
                 speciesList   = [],
                 symmetryGroup = None):
        self.unitCell      = unitCell
        self.speciesList   = speciesList
        self.symmetryGroup = symmetryGroup

    def cartCoords(self):
        ''' List of cart coords - remove shells later if needed'''

        if self.speciesList[0].cartCoord is None:
            self.setCartCoord()

        return np.array([x.cartCoord for x in self.speciesList])

    def fracCoords(self):
        ''' List of frac coords - remove shells later if needed'''

        if self.speciesList[0].fracCoord is None:
            self.setFracCoord()

        return np.array([x.fracCoord for x in self.speciesList])


    @classmethod
    def subSection(cls, parentStructure, limits):
        ''' Returns new structure which is cut-out
            just take the bit between limits (in frac coords)
            should use a P1 cell 
            limits are np.array([[lower bound vector], [upper bound vector]])'''

        lengthReductionRatio = limits[1] - limits[0]
        invRatio             = np.array([1./x for x in lengthReductionRatio])

        newSpeciesList = [x for x in parentStructure.speciesList if all([x.fracCoord[i] < limits[1][i] for i in xrange(3)]) and
                                                                            all([x.fracCoord[i] > limits[0][i] for i in xrange(3)])]

        for i in xrange(len(newSpeciesList)):
            tempFC = newSpeciesList[i].fracCoord * invRatio
            newSpeciesList[i].fracCoord = tempFC - np.floor(tempFC)

        return cls(unitCell    = UnitCell(vectors = np.array([parentStructure.unitCell.vectors[i] * (limits[1][i] - limits[0][i])  for i in xrange(3)])),
                   speciesList = newSpeciesList)

    @classmethod
    def cutOutParallelepiped(cls, parentStructure, parallelepiped, origin = np.zeros(3), periodicTestCutoff = None):
        ''' Rather than just take sub section (above), can allow a general parallelepiped to be cut out
            Should be same as above for orthorhombic cell and no rotation (not checked yet) 
            parallelepiped is np.array().shape = (3,3)
            periodicTestCutoff is in \\A - it will look to see whether another atom is found just outside box on other side '''

        from copy import deepcopy

        invVectors     = np.linalg.inv(parallelepiped)
        newSpeciesList = []

        parentStructure.setFracCoord()
        for i in xrange(len(parentStructure.speciesList)):
#            fracCoord = np.dot(parentStructure.speciesList[i].cartCoord, invVectors)
            fracCoord = np.dot(parentStructure.speciesList[i].cartCoord - origin, invVectors)
            if all([x >= 0. and x < 1. for x in fracCoord]):
                newSpeciesList.append(Species(element   = parentStructure.speciesList[i].element,
                                              core      = parentStructure.speciesList[i].core,
                                              charge    = parentStructure.speciesList[i].charge,
                                              fracCoord = fracCoord,
                                              cartCoord = np.dot(fracCoord, parallelepiped)))
#                parentStructure.speciesList[i].fracCoord = fracCoord
#                parentStructure.speciesList[i].cartCoord = np.dot(fracCoord, parallelepiped)
#                newSpeciesList.append(deepcopy(parentStructure.speciesList[i]))

        #note that at this stage returning an unrotated cell 050517
        return cls(unitCell    = UnitCell(vectors = parallelepiped),
                   speciesList = newSpeciesList)

    def withinPeriodicEnvelopment(self, bigCell, dist, distVec = None):
        ''' quasi periodicity is tested- are there similar atoms to ones near edges one cell translation away?
            make big cell to include things that were not in self 
            distVec is vector of length dist along (1,1,1) - (0,0,0) '''

        # if not given a vector to pick fractional cutoffs, use this
        if distVec is None:
            distVec = np.ones(3) * 3**-0.5 * dist
        fractionalCutoffs = np.vstack([np.dot(distVec, np.linalg.inv(self.unitCell.vectors)),
                                       1. - np.dot(distVec, np.linalg.inv(self.unitCell.vectors))])

        # total if statements = 6 faces, 12 edges, 8 vertices = 26 (labelled f,c,v below)
        for at in self.speciesList:
            #fx0
            if at.fracCoord[0] < fractionalCutoffs[0, 0]:
                if not any([np.linalg.norm(at.cartCoord - x.cartCoord + self.unitCell.vectors[0]) < dist for x in bigCell.speciesList if x.element == at.element]):
                    return False
                #exy0
                if at.fracCoord[1] < fractionalCutoffs[0, 1]:
                    if not any([np.linalg.norm(at.cartCoord - x.cartCoord + self.unitCell.vectors[0] + self.unitCell.vectors[1]) < dist for x in bigCell.speciesList if x.element == at.element]):
                        return False
                    #vxyz0
                    if at.fracCoord[2] < fractionalCutoffs[0, 2]:
                        if not any([np.linalg.norm(at.cartCoord - x.cartCoord + self.unitCell.vectors[0] + self.unitCell.vectors[1] + self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                            return False
                    #vxyz1
                    if at.fracCoord[2] > fractionalCutoffs[1, 2]:
                        if not any([np.linalg.norm(at.cartCoord - x.cartCoord + self.unitCell.vectors[0] + self.unitCell.vectors[1] - self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                            return False

                #exy1
                if at.fracCoord[1] > fractionalCutoffs[1, 1]:
                    if not any([np.linalg.norm(at.cartCoord - x.cartCoord + self.unitCell.vectors[0] - self.unitCell.vectors[1]) < dist for x in bigCell.speciesList if x.element == at.element]):
                        return False
                    #vxzy2
                    if at.fracCoord[2] < fractionalCutoffs[0, 2]:
                        if not any([np.linalg.norm(at.cartCoord - x.cartCoord + self.unitCell.vectors[0] - self.unitCell.vectors[1] + self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                            return False
                    #vxyz3
                    if at.fracCoord[2] > fractionalCutoffs[1, 2]:
                        if not any([np.linalg.norm(at.cartCoord - x.cartCoord + self.unitCell.vectors[0] - self.unitCell.vectors[1] - self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                            return False

                #exz0
                if at.fracCoord[2] < fractionalCutoffs[0, 2]:
                    if not any([np.linalg.norm(at.cartCoord - x.cartCoord + self.unitCell.vectors[0] + self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                        return False
                #exz1
                if at.fracCoord[2] > fractionalCutoffs[1, 2]:
                    if not any([np.linalg.norm(at.cartCoord - x.cartCoord + self.unitCell.vectors[0] - self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                        return False

            #fx1
            if at.fracCoord[0] > fractionalCutoffs[1, 0]:
                if not any([np.linalg.norm(at.cartCoord - x.cartCoord - self.unitCell.vectors[0]) < dist for x in bigCell.speciesList if x.element == at.element]):
                    return False
                #exy2
                if at.fracCoord[1] < fractionalCutoffs[0, 1]:
                    if not any([np.linalg.norm(at.cartCoord - x.cartCoord - self.unitCell.vectors[0] + self.unitCell.vectors[1]) < dist for x in bigCell.speciesList if x.element == at.element]):
                        return False
                    #vxyz4
                    if at.fracCoord[2] < fractionalCutoffs[0, 2]:
                        if not any([np.linalg.norm(at.cartCoord - x.cartCoord - self.unitCell.vectors[0] + self.unitCell.vectors[1] + self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                            return False
                    #vxyz5
                    if at.fracCoord[2] > fractionalCutoffs[1, 2]:
                        if not any([np.linalg.norm(at.cartCoord - x.cartCoord - self.unitCell.vectors[0] + self.unitCell.vectors[1] - self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                            return False

                #exy3
                if at.fracCoord[1] > fractionalCutoffs[1, 1]:
                    if not any([np.linalg.norm(at.cartCoord - x.cartCoord - self.unitCell.vectors[0] - self.unitCell.vectors[1]) < dist for x in bigCell.speciesList if x.element == at.element]):
                        return False
                    #vxzy6
                    if at.fracCoord[2] < fractionalCutoffs[0, 2]:
                        if not any([np.linalg.norm(at.cartCoord - x.cartCoord - self.unitCell.vectors[0] - self.unitCell.vectors[1] + self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                            return False
                    #vxyz7
                    if at.fracCoord[2] > fractionalCutoffs[1, 2]:
                        if not any([np.linalg.norm(at.cartCoord - x.cartCoord - self.unitCell.vectors[0] - self.unitCell.vectors[1] - self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                            return False

                #exz2
                if at.fracCoord[2] < fractionalCutoffs[0, 2]:
                    if not any([np.linalg.norm(at.cartCoord - x.cartCoord - self.unitCell.vectors[0] + self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                        return False
                #exz3
                if at.fracCoord[2] > fractionalCutoffs[1, 2]:
                    if not any([np.linalg.norm(at.cartCoord - x.cartCoord - self.unitCell.vectors[0] - self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                        return False

            #fy0
            if at.fracCoord[1] < fractionalCutoffs[0, 1]:
                if not any([np.linalg.norm(at.cartCoord - x.cartCoord + self.unitCell.vectors[1]) < dist for x in bigCell.speciesList if x.element == at.element]):
                    return False
                #eyz0
                if at.fracCoord[2] < fractionalCutoffs[0, 2]:
                    if not any([np.linalg.norm(at.cartCoord - x.cartCoord + self.unitCell.vectors[1] + self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                        return False
                #eyz1
                if at.fracCoord[2] > fractionalCutoffs[1, 2]:
                    if not any([np.linalg.norm(at.cartCoord - x.cartCoord + self.unitCell.vectors[1] - self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                        return False
            #fy1
            if at.fracCoord[1] > fractionalCutoffs[1, 1]:
                if not any([np.linalg.norm(at.cartCoord - x.cartCoord - self.unitCell.vectors[1]) < dist for x in bigCell.speciesList if x.element == at.element]):
                    return False
                #eyz2
                if at.fracCoord[2] < fractionalCutoffs[0, 2]:
                    if not any([np.linalg.norm(at.cartCoord - x.cartCoord - self.unitCell.vectors[1] + self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                        return False
                #eyz3
                if at.fracCoord[2] > fractionalCutoffs[1, 2]:
                    if not any([np.linalg.norm(at.cartCoord - x.cartCoord - self.unitCell.vectors[1] - self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                        return False

            #fz0
            if at.fracCoord[2] < fractionalCutoffs[0, 2]:
                if not any([np.linalg.norm(at.cartCoord - x.cartCoord + self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                    return False
            #fz1
            if at.fracCoord[2] > fractionalCutoffs[1, 2]:
                if not any([np.linalg.norm(at.cartCoord - x.cartCoord - self.unitCell.vectors[2]) < dist for x in bigCell.speciesList if x.element == at.element]):
                    return False

        return True


    def changeUnitCell(self, inLimits, retainUnitCell = False):
        ''' Changes this unit cell, such that atoms are now between new limits 
            N.B. only an origin at 0,0,0 is used for definition of unit cell atm 
            and this subroutine assumes that input self is in usual cell convention - [0,1)^3 (I believe) 

            example use: myStructureInstance.changeUnitCell(np.array([[-1,2], [-1,2], [-1,2]]),
                                                            retainUnitCell = True) '''

        #transpose if given rows <-> columns
        if inLimits.shape == (2,3):
            inLimits = inLimits.T

        # ensure min first
        limits = np.array([[min(l), max(l)] for l in inLimits])

        # upper limit is ceil because start in 0,1 and xrange stops before this number- to things cancel 
        possibleTranslations = np.array([[a, b, c] for a in xrange(int(np.floor(limits[0, 0])), int(np.ceil(limits[0, 1])))
                                                   for b in xrange(int(np.floor(limits[1, 0])), int(np.ceil(limits[1, 1])))
                                                   for c in xrange(int(np.floor(limits[2, 0])), int(np.ceil(limits[2, 1])))])

        if self.speciesList[0].fracCoord is None:
            self.setFracCoord()

        newSpeciesList = []
        for ns, s in enumerate(self.speciesList):
            for p in possibleTranslations:
                newPt = s.fracCoord + p
                if newPt[0] >= limits[0, 0] and newPt[0] < limits[0, 1] and\
                   newPt[1] >= limits[1, 0] and newPt[1] < limits[1, 1] and\
                   newPt[2] >= limits[2, 0] and newPt[2] < limits[2, 1]:
                    newSpeciesList.append(Species(element   = s.element,
                                                  core      = s.core,
                                                  fracCoord = newPt))

        self.speciesList = newSpeciesList
        self.setCartCoord()
        # catch errors hopefully with this
        if not retainUnitCell:
            self.unitCell    = None

        return self

    def lllReduce(self):
        ''' this code is taken from PyMatGen and uses its functionality
            see pmg.core.surface for an example (hence name slab) 
            run this as structure = structure.lllReduce() '''

        slab = self.toPMGStructure()
        lll_slab = slab.copy(sanitize=True)
        mapping = lll_slab.lattice.find_mapping(slab.lattice)

        return self.fromPMGStructure(lll_slab)

    def uniqueSpecies(self, attrList):
        ''' Return list, worry about changes if need be'''
        from setTools import subsetByAttributes
        from copy import deepcopy
        tempList = []
        for x in self.speciesList:
            if not subsetByAttributes([x], tempList, attrList):
                tempList.append(deepcopy(x))
        return tempList

    def setFracCoord(self):
        ''' CHECK CORRECT WAY AROUND!!!!!!!!!!!!!!
            assumes have lattice (vectors) set '''
#        if not self.unitCell.invVectors:
        if self.unitCell.invVectors is None:
            self.unitCell.setInvVectors()

        for i in xrange(len(self.speciesList)):
            self.speciesList[i].fracCoord = np.dot(self.speciesList[i].cartCoord, self.unitCell.invVectors)

    def setCartCoord(self):
        ''' CHECK CORRECT WAY AROUND!!!!!!!!!!!!!!
            assumes have lattice (vectors) set '''
        for i in xrange(len(self.speciesList)):
            self.speciesList[i].cartCoord = np.dot(self.speciesList[i].fracCoord, self.unitCell.vectors)

    def resetFracCoord(self):
        ''' set all Frac Coord to be [0,1)^3 '''
        for i in xrange(len(self.speciesList)):
            self.speciesList[i].fracCoord -= np.floor(self.speciesList[i].fracCoord)

        # if there are cart coords, reset these too
#        if self.speciesList[0].cartCoord is not None: Untested Change 150317
        if self.speciesList[0].cartCoord is None:
            self.setCartCoord()            
        print "Please set good tests for all these- consider different programs and their conventions"

    def applyFracDisplacement(self, fracDisplacement):
        ''' Move all species '''
        # calc Frac Coord if don't already have it
        if self.speciesList[0].fracCoord is None:
            self.setFracCoord()

        # move fracCoord
        for i in xrange(len(self.speciesList)):
            self.speciesList[i].fracCoord += fracDisplacement

        # if using cartCoord, recalc it
        if self.speciesList[0].cartCoord is not None:
            self.setCartCoord()


    def displaceAllSpecies(self, cartDisplacement = None, fracDisplacement = None):
        ''' Do cartDisplacement if needed 
            frac takes precedence if give both (DONT DO THIS!!) '''
        if fracDisplacement is not None:
            self.applyFracDisplacement(fracDisplacement)
            return True


    def addShells(self, newSpecies, displacement = np.zeros(3), cartesian = False):
        ''' can add displacement to a shell to avoid energy catastrophe
            note, either works on cartesian or fractional- possible bugs if both needed '''
        from setTools import sameElementByAttributes
        from copy import deepcopy as dc
        newSpecies.core = 'shel'
        tempList        = []
        for x in self.speciesList:
            if sameElementByAttributes(x, newSpecies, ['element']):
                if cartesian:
                    newSpecies.cartCoord = x.cartCoord + displacement
                else:
                    newSpecies.fracCoord = x.fracCoord + displacement
                tempList.append(dc(newSpecies))
        self.speciesList += tempList
        return len(self.speciesList)

    def removeShells(self):
        self.speciesList = [x for x in self.speciesList if x.core.lower()[0] != 's']
        return len(self.speciesList)

    #@property??
    def density(self, numberDensity = False):
        from ase.data import atomic_masses, chemical_symbols
        totalMass = sum([atomic_masses[chemical_symbols.index(t.element)] for t in self.speciesList])
        volume    = np.linalg.det(self.unitCell.vectors)
        if numberDensity:
            return float(len(self.speciesList)) / volume
        return totalMass / volume

    def addStandardCharges(self):
        ''' Assume all shells already have charge set '''
        uniqueSpecies  = self.uniqueSpecies(['element', 'core'])
        uniqueElements = self.uniqueSpecies(['element'])
        shellCharges = dict([(x.element, x.charge) for x in uniqueSpecies if x.core.lower()[0] == 's'])

        for at in uniqueElements:
            if any([(x.element == at.element and x.core.lower()[0] == 's') for x in uniqueSpecies]):
                self.setCharges([Species(element = at.element,
                                        core    = at.core,
                                        charge  = at.defaultCharge() - shellCharges[at.element])],
                                onlySetSome = True)
            else:
                self.setCharges([Species(element = at.element,
                                        core    = at.core,
                                        charge  = at.defaultCharge())],
                                onlySetSome = True)

        return self.setCharges([])

    def setCharges(self, templateSpecies, onlySetSome = False):
        ''' Pass in some templateSpecies- anything in speciesList that matches these gets the template's charge 
            use onlySetSome if this crashes because some are left unset '''
        from setTools import sameElementByAttributes
        for t in templateSpecies:
            for xs, s in enumerate(self.speciesList):
                if sameElementByAttributes(t, s, ['element', 'core']):
                    self.speciesList[xs].charge = t.charge
        if onlySetSome:
            return None
        return sum(x.charge for x in self.speciesList)

    def setASEMasses(self, templateSpecies):
        ''' Pass in some templateSpecies- anything in speciesList that matches these gets mass from ase.data  '''
        from setTools import sameElementByAttributes
        from ase.data import atomic_masses, chemical_symbols
        for t in templateSpecies:
            for xs, s in enumerate(self.speciesList):
                if sameElementByAttributes(t, s, ['element', 'core']):
                    self.speciesList[xs].mass = atomic_masses[chemical_symbols.index(t.element)]
        return #sum(x.charge for x in self.speciesList)

    def changeMassesForShells(self, templateSpecies, massShell = 0.):
        ''' Take some mass from core and add to shell where matches templateSpecies 
            Mass must have already been set for core '''
        from setTools import sameElementByAttributes
        if not  isinstance(templateSpecies, list):
            templateSpecies = [templateSpecies]
        for t in list(templateSpecies):
            for xs, s in enumerate(self.speciesList):
#                if sameElementByAttributes(t, s, ['element', 'core']):
                if sameElementByAttributes(t, s, ['element']):
                    if s.core == 'core':
                        self.speciesList[xs].mass -= massShell
                    elif s.core[:4] == 'shel':
                        self.speciesList[xs].mass = massShell
        return sum(x.charge for x in self.speciesList)

#    def expandAllSymmetry(self):
#        ''' Convert this into a P1 cell 
#            Use PMG?? '''
#
#        if self.symmetryGroup is not None:
#            self.symmetryGroup = SymmetryGroup(number = 1)
# 
#        raise Exception('Better to just read a cif with PMG or ASE for now')

#    @classmethod
    def superCell(self, superCell):
        ''' returns a superCell - input (1,1,1) for eg (tuple of ints (fracs??))'''
        return Structure.fromASEStructure( self.toASEStructure() * superCell )

    def toASEStructure(self):
        from ase import Atoms

        if self.speciesList[0].cartCoord is None:
            self.setCartCoord()

        return Atoms([x.toASEAtom() for x in self.speciesList],
                     cell = self.unitCell.vectors)
        
    @classmethod
    def fromASEStructure(cls, aseStructure):
        ''' Reading in an ASE structure can help make supercells
            the ASE cif reader doesn't read much information- just the basics '''
        return cls(unitCell    = UnitCell(vectors = aseStructure._cell),
                   speciesList = [Species.initFromASEAtom(x) for x in aseStructure])

    @classmethod
    def fromGULPOutput(cls, gulpOutputString):
        ''' Makes P1 cell from gulp output
            assumes that you've run an opti calculation (the comp writes info at the bottom and can be used, but see other subroutines, not this)
            please check this if running in a novel manner '''
        # test that can be turned into PMG structure or something and same as the custom gulp thing below
        from gulpIO import OutputFileGULP
        gulpOutputDict = OutputFileGULP().readOutputOpti(gulpOutputString)

        return cls(unitCell    = gulpOutputDict['unitCell'],
                   speciesList = gulpOutputDict['speciesList'])
    
    def toPMGStructure(self):
        ''' pymatgen structure- assume P1 cell '''
        from pymatgen.core import Structure as PMGS
        from pymatgen.core import Lattice   as PMGL

        # only concerned with core atoms (no shells)
        coreSpeciesIndices = [i for i in range(len(self.speciesList)) if self.speciesList[i].core == 'core']
        outLattice   = PMGL.from_parameters(self.unitCell.lengths[0],
                                            self.unitCell.lengths[1],
                                            self.unitCell.lengths[2],
                                            self.unitCell.angles[0],
                                            self.unitCell.angles[1],
                                            self.unitCell.angles[2])
        return PMGS.from_spacegroup(1,
                                    outLattice,
                                    [self.speciesList[i].element   for i in coreSpeciesIndices],
                                    [self.speciesList[i].fracCoord for i in coreSpeciesIndices])       

    @classmethod
    def fromPMGStructure(cls, pmgStructure):
        return cls(unitCell    = UnitCell(vectors = pmgStructure._lattice._matrix),
                   speciesList = [Species().initFromPMGSite(x) for x in pmgStructure._sites])

    @classmethod
    def fromXYZ(cls, xyzFileName, unitCell=None, transposeVectors=False):
        ''' Transpose vectors may be important if xyz from LAMMPS with lower triangular lattice '''
        xyzF = XYZFile(fileName = xyzFileName)

        if not unitCell:
            unitCell = UnitCell.fromXYZFile(xyzF, transposeVectors=transposeVectors)

        speciesList = XYZFile.xyzStringToSpeciesList(open(xyzF.fileName, 'r').read())

        return cls(unitCell    = unitCell,
                   speciesList = speciesList)

    @classmethod
    def fromCIF(cls, cifName, expandFullCell = False):
        ''' cifName is just the name of the file (inc. location if different directory)
            expandFullCell will return a P1 cell and use PMG to expand all symmetry operators 
            N.B. two methods dont give same order, and no cartCoords from default method '''

        if expandFullCell:
            try:
                from pymatgen.core import Structure as PMGS
                return cls.fromPMGStructure(PMGS.from_file(cifName))
            except ValueError:
                print 'Failure for pymatgen - try using ASE'
                from ase.io import read as ASERead
                return cls.fromASEStructure(aseStructure = ASERead(cifName))

        from cifIO import getUnitCellCIF, getSpeciesListCIF, getSymmetryGroupCIF
        try:
            symmetryGroup = getSymmetryGroupCIF(cifName)
        except:
            print "Couldn't determine symmetry group"
            symmetryGroup = None
        return cls(unitCell      = getUnitCellCIF     (cifName),
                   speciesList   = getSpeciesListCIF  (cifName),
                   symmetryGroup = symmetryGroup)


    def inputCIF(self, cifName):
        from cifIO import getUnitCellCIF, getSpeciesListCIF, getSymmetryGroupCIF
        self.unitCell      = getUnitCellCIF     (cifName)
        self.speciesList   = getSpeciesListCIF  (cifName)
        self.symmetryGroup = getSymmetryGroupCIF(cifName)

    def toCif(self, cifName):
        ''' Only write P1 at the moment- quickCif = True '''
        from cifIO import writeCIF
        writeCIF(cifName, self, quickCif = True)

    def toXYZ(self, xyzName, includeShells=False, includeCellVectors=True, reorderAlphabetically=True):
        with open(xyzName, 'w') as outf:
            outf.write(self.xyzString(includeShells=includeShells,
                                      includeCellVectors=includeCellVectors,
                                      reorderAlphabetically=reorderAlphabetically))

    def removeListIndices(self, deletionsList):
        ''' Deletes the species with these indices (it will do reversing for you) '''

        for i in deletionsList[::-1]:
            del(self.speciesList[i])

        return len(self.speciesList)

    def removeListSpecies(self, speciesToMatch, attributes = ['element']):
        ''' Removes a list of matching species from self.speciesList
            Returns this list separately '''

        from setTools import sameElementByAttributes

        separatedList = []
        deletionsList = []

        for indexS, s in enumerate(self.speciesList):
            if sameElementByAttributes(s, speciesToMatch, attributes):
                separatedList.append(self.speciesList[indexS])
                deletionsList.append(indexS)

        for i in deletionsList[::-1]:
            del(self.speciesList[i])

        return separatedList

    def removeSpeciesSobol(self, species, nRemove, offset = 0, swapSpecies = None, cartesian = True):
        ''' Remove nRemove species, that match character of species, according to closest fractional 
            coordinates to those from sobol list 
            Will only work well if atoms are sensibly bounded within lattice vectors (as looks in [0,1)^3 . LatticeVectors)
            Deletions are in order 
            if give a swapSpecies , deleted species -> this new one '''

        from setTools import sameElementByAttributes
        from sobol_lib import i4_sobol
   
        #Assume cartesian coordinates are wanted at this stage- easy to change if needed
        attributesList = [k for k in species.__dict__.keys() if species.__dict__[k] is not None]

#        matchingElements = [x for x in self.speciesList if sameElementByAttributes(x, species, attributesList)]
        deletedElements  = []

        for n in xrange(offset, offset + nRemove):
            sobVec, newSeed = i4_sobol(3, n)
            if cartesian:
                removeCentre = np.dot(sobVec, self.unitCell.vectors)
                distance, savedPosition = np.inf, None
#                for im, m in enumerate(matchingElements):
                for im, m in enumerate(self.speciesList):
                    if not sameElementByAttributes(m, species, attributesList):
                        continue
                    newDistance = np.linalg.norm(m.cartCoord - removeCentre)
                    if newDistance < distance and im not in deletedElements:
                        distance, savedPosition = newDistance, im
            deletedElements.append(savedPosition)

        # either replace or delete the selected species (only change the element - add charge/mass another time??)
        for i in deletedElements:
            if swapSpecies is not None:
                self.speciesList[i].element = swapSpecies.element
            else:
                del(self.speciesList[i])

        return len(self.speciesList)

    def customDisplacementsGULP(self, gulpOutput, originalCIF):
        ''' This is a subroutine to extract information from a gulpOutput
            In the original CIF are some atoms. This returns cell lengths (1st 3 elements)
            and then nearest neighbour distances for the atoms in the CIF
            called custom as this has little control, testing, weighting etc '''
        # Import pymatgen stuff to find nearest neighbours
        from pymatgen.core import Structure as PMGS
        from pymatgen.core import Lattice   as PMGL

        # only concerned with core atoms (no shells)
        coreSpeciesIndices = [i for i in range(len(self.speciesList)) if self.speciesList[i].core == 'core']
        outLattice   = PMGL.from_parameters(*list(gulpOutput['newLengths']) + list(gulpOutput['newAngles']))
        outStructure = PMGS.from_spacegroup(self.symmetryGroup.number,
                                            outLattice,
                                            [self.speciesList[i].element for i in coreSpeciesIndices],
                                            [gulpOutput['newCoords'][i] for i in coreSpeciesIndices])

        sitesNewStructure = []

        for x in [gulpOutput['newCoords'][i] for i in coreSpeciesIndices]:
            for s in outStructure._sites:
                if np.allclose(s._fcoords, x):
                    sitesNewStructure.append(s)

        newDisps = np.array(list(gulpOutput['newLengths']) +
                            [sorted(outStructure.get_neighbors(s, 4.), key = lambda x:x[1])[0][1] for s in sitesNewStructure])

        sitesOriginalStructure = []
        originalStructure = PMGS.from_file(originalCIF)

        for x in [self.speciesList[i].fracCoord for i in coreSpeciesIndices]:
            for s in originalStructure._sites:
                if np.allclose(s._fcoords, x, atol=1e-1):
                    sitesOriginalStructure.append(s)

        originalDisps = np.array(list(originalStructure.lattice._lengths) +
                          [sorted(originalStructure.get_neighbors(s, 4.), key = lambda x:x[1])[0][1] for s in sitesOriginalStructure])

        return originalDisps, newDisps

#    def closestSpeciesInListMatchIndex(self, targetSpecies, targetCoordinate, cartesian = False, matchAttributes = ['element']):
#        ''' Return the index of the closest Species (of correct type) '''
#
#        from setTools import sameElementByAttributes
#
#        if cartesian:
#            raise Exception("not implemented cartesian coord in Structure.closestSpeciesMatchIndex")
#
#        if self.speciesList[0].fracCoord is None:
#            self.resetFracCoords()
#
#        for iat, at in enumerate(self.speciesList):
#             if not sameElementByAttributes(targetSpecies, at, matchAttributes) :
#                 continue
#             disp =  np.linalg.norm(targetCoordinate - at.fracCoord)
#             if disp < minDisp and iat not in outList:
#                 minDisp, index = disp, iat

    def speciesMatchIndices(self,
                            targetSpecies   = None,
                            matchAttributes = ['element']):
        ''' Return list of indices of matches to targetSpecies '''

        from setTools import sameElementByAttributes

        return [x[0] for x in enumerate(self.speciesList) if sameElementByAttributes(x[1], targetSpecies, matchAttributes)]

    def checkOverlappingAtoms(self,
                              cutoff,
                              adjustAtoms = False):
        ''' confirm that cartCoords do not overlap (quicker check than radialDistributionFunction
            no PBC '''
        cartCoords = self.cartCoords()
        cutSquare  = cutoff**2
        foundAny   = False
        swapIndices = []
        eta        = 1.e-6 # numerical consideration to ensure push atoms apart enough
        for i1 in xrange(cartCoords.shape[0]):
            for i2 in xrange(cartCoords.shape[0]):
                if i1 == i2:
                    continue
                distVec = cartCoords[i1]- cartCoords[i2]
                if np.dot(distVec, distVec) < cutSquare:
                    if adjustAtoms:
                        foundAny = True
                        #move second atom (so can test later) (remember to move cartCoords[i2])
                        self.speciesList[i2].cartCoord = self.speciesList[i2].cartCoord - distVec * (cutoff - np.linalg.norm(distVec) + eta) / np.linalg.norm(distVec)
                        cartCoords[i2] = cartCoords[i2] - distVec * (cutoff - np.linalg.norm(distVec) + eta) / np.linalg.norm(distVec)
                    else:
                        return True
        return foundAny

    def removeOverlappingAtoms(self,
                               cutoff,
                               maxIt = 10):
        ''' Loop through checkOverlappingAtoms, pushing atoms to cutoff length apart 
            Return true if structure can be made to have no overlapping atoms '''
        for i in xrange(maxIt):
            if not self.checkOverlappingAtoms(cutoff, adjustAtoms = True):
                return True
        return False

    def indexNearestSpeciesToPoint(self, pt):
        ''' Return index of closest species to a pt (cartesian) (use faster routines if ever needed) '''
        bestDist = np.inf
        for i in xrange(len(self.speciesList)):
            newDist = np.linalg.norm(pt - self.speciesList[i].cartCoord)
            if newDist < bestDist:
                bestDist, closestI = newDist, i
        return closestI, bestDist

    def nearestNeighbourDistance(self,
                                 speciesListIndex = None, 
                                 targetSpecies    = None,
                                 cutoffRadius     = None):
        ''' Nearest neighbour (of particular type) is first instance in RDF 
            centred at speciesListIndex '''

        return self.radialDistributionFunction(speciesListIndex = speciesListIndex,
                                               targetSpecies    = targetSpecies,
                                               cutoffRadius     = cutoffRadius)[0]

    def radialDistributionFunction(self, 
                                   speciesListIndex     = None, 
                                   targetSpecies        = None,
                                   cutoffRadius         = None,
                                   returnVectors        = False,
                                   returnSpeciesList    = False):
        ''' Give a species list index, and get displacement list of nearest neighbours from this (\AA)
            targetSpecies is optional, and only the element is used at this stage 
            if returnVectors these are cartesian vectors 
            if returnSpeciesList this is a list of Species, with positions WRT central atom '''

        # ensure that fractional coordinates are available
#        if self.speciesList[0].fracCoord is None:
        if any(x.fracCoord is None for x in self.speciesList):
            self.setFracCoord()

        # cutoffRadius is for looking for neighbours- 4. is only really appropriate for close things
        if not cutoffRadius:
            cutoffRadius = 4.

        # start by temporarily making a PMG structure- use this functionality and watch fracCoord bounds
        centralFracCoord = np.mod(self.speciesList[speciesListIndex].fracCoord, 1)
        tempPMGStructure = self.makePMGStructure() 
        try:
            centralSite = [s for s in tempPMGStructure._sites if np.allclose(s._fcoords, centralFracCoord, atol=1.e-6)][0]
        except:
            print "Issues for coordinate at ", centralFracCoord #for some reason np.mod doesn't work at fracCoord = 1.
            tempOut = []
            for i in xrange(3):
                if abs(centralFracCoord[i] - 1.) < 0.0000001:
                    tempOut.append(0.)
                else:
                    tempOut.append(centralFracCoord[i])
            centralFracCoord = np.array(tempOut)
            centralSite = [s for s in tempPMGStructure._sites if np.allclose(s._fcoords, centralFracCoord, atol=1.e-6)][0]

        tempPMGNeighbours = sorted(tempPMGStructure.get_neighbors(centralSite, cutoffRadius), key = lambda x:x[1])
        if targetSpecies:
            listNeighbours = [x for x in tempPMGNeighbours if x[0]._species._data.keys()[0].symbol == targetSpecies.element]
        else:
            listNeighbours = tempPMGNeighbours

        if returnVectors:
            return [x[0]._coords - centralSite._coords for x in listNeighbours]
        elif returnSpeciesList:
            return [Species(element = Species.elementFromPMGSpecies(x[0]._species),
                            cartCoord = x[0]._coords - centralSite._coords) for x in listNeighbours]
#            print "NOt implaMEETNNTNTED yet";exit()
        else:
#            return [x[1] for x in tempPMGNeighbours]
            return [x[1] for x in listNeighbours]

    def numberValenceElectrons(self):
        return sum([x.atomicValenceElectrons() for x in self.speciesList if x.core == 'core'])

    def composition(self, considerShellsSeparately = False, denominator = None):
        ''' Return dict with elements and number of times they appear (ignore the shells stuff unless needed 
            denominator divides all numbers if present '''

        from setTools import subsetByAttributes

        if denominator:
            return dict(zip([x.element for x in self.uniqueSpecies(['element'])],
                            [float(len([x for x in self.speciesList if x.element == y.element and x.core[:4].lower() == 'core'])) / float(denominator) for y in self.uniqueSpecies(['element'])]))

        return dict(zip([x.element for x in self.uniqueSpecies(['element'])],
                        [len([x for x in self.speciesList if x.element == y.element and x.core[:4].lower() == 'core']) for y in self.uniqueSpecies(['element'])]))

    def compositionString(self, considerShellsSeparately = False, denominator = None):
        ''' Return a string with the composition '''

        composition = self.composition(considerShellsSeparately = considerShellsSeparately, denominator = denominator)

        return "; ".join(map(str, [" ".join(map( str, [composition[k], k])) for k in composition.keys()]))

    def makePMGStructure(self):
        from pymatgen.core import Structure as PMGS
        from pymatgen.core import Lattice   as PMGL

        outLattice   = PMGL.from_parameters(*list(self.unitCell.lengths) + list(self.unitCell.angles))

        if self.symmetryGroup is None or not self.symmetryGroup.number:
            print "symmetryGroup.number is unknown -> P1 cell"
            self.symmetryGroup = SymmetryGroup(number = 1)

        outStructure = PMGS.from_spacegroup(self.symmetryGroup.number,
                                            outLattice,
                                            [x.element   for x in self.speciesList if x.core[:4] == 'core'],
                                            [x.fracCoord for x in self.speciesList if x.core[:4] == 'core'])
        return outStructure

    def xyzString(self, includeShells=False, includeCellVectors=False, reorderAlphabetically=True):
        ''' make .xyz format, for PLUMED for e.g. '''
        import copy

        #print a dummy structure with atoms reordered
        newStructure = copy.deepcopy(self)
        if reorderAlphabetically:
            print "Reordering species alphabetically in xyzString"
            newStructure.speciesList = sorted([x for x in self.speciesList], key = lambda x: x.element)

        if newStructure.speciesList[0].cartCoord is None:
            newStructure.setCartCoord()

        outString = str(len([x for x in newStructure.speciesList if x.core[:4].lower() == 'core' or includeShells])) + "\n"

        # possibly include cell lengths- periodicity may not be possible in Plumed???
        if includeCellVectors:
            outString += " ".join(map(str, list(newStructure.unitCell.vectors.ravel()))) + "\n"
        else:
            outString += " Not included cell vectors\n"

        outString += "\n".join([" ".join([x.element] + map(str, x.cartCoord)) for x in newStructure.speciesList if x.core[:4].lower() == 'core' or includeShells])

        # the \n is essential for PLUMED
        return outString + "\n"

    def optimiseGULP(self, 
                     potentials = None, 
                     addShells = False, 
                     addStandardCharges = True, 
                     changeNonStandardCharges = False,
#                     returnGULPDict = True,  #always return a dict
                     timeout = 10000., 
                     keepShell = False):
        ''' If addShell set it to a list ie [Species(element = 'O',
                                                     charge  = -2.86)]
            timeout in seconds (can set to None) 
            returns (Structure, dictionary) '''
        import copy
        from hardcode           import GULP_EXE
        from subprocessHandling import RunCommandNew
        from gulpIO             import InputFileGULP, OutputFileGULP
        

        returnGulpDict = {'initialEnergy'   : None,
                          'finalEnergy'     : None,
                          'initialStructure': copy.deepcopy(self),
                          'finalStructure'  : None,
                          'phonons'         : None,
                          'phononsPositive' : None,
                          'completeSuccess' : False}

        if addShells is not False:
            if type(addShells) is not type([None]):
                addShells = list(addShells)
            for s in addShells:
                self.addShells(Species(element= s.element,
                                       charge = s.charge))

        if addStandardCharges:
            self.addStandardCharges()
        if changeNonStandardCharges:
            self.setCharges(changeNonStandardCharges)

        assert(abs(self.setCharges([])) < 1.e-5)

        inFile = InputFileGULP(fileName = '',
                               title    = 'generic title',
                               parentStructure = self,
                               keywords   = ['conp', 'opti', 'phonon'],
                               potentials = potentials)

        runner = RunCommandNew(GULP_EXE, inFile.stringForm())
        try:
            runner.run(timeout=timeout)
        except:
            return (None, returnGulpDict)

        try:
            returnGulpDict['phonons']         = OutputFileGULP().phononFrequencies(runner.output)
            returnGulpDict['phononsPositive'] = OutputFileGULP().phononFrequenciesPositive(runner.output)
            returnGulpDict['latticeEnergies'] = [x for x in OutputFileGULP.findAllLatticeEnergies(runner.output)]
            returnGulpDict['initialEnergy']   = returnGulpDict['latticeEnergies'][0]
            returnGulpDict['finalEnergy']     = returnGulpDict['latticeEnergies'][-1]
        except:
            return (None, returnGulpDict)

        self = Structure.fromGULPOutput(runner.output)

        if not keepShell:
            self.removeShells()

        returnGulpDict['completeSuccess'] = True
        return (self, returnGulpDict)
