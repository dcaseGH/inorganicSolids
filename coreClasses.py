#Perhaps make this good if needed

import numpy as np

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

    @classmethod
    def initFromASEAtom(cls, aseAtom):
        return cls(cartCoord = aseAtom.position, 
                    element  =  aseAtom.symbol)

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

    def energy(self, r):
        ''' Units are same as A or C6, r is rho^{-1} '''
        return self.A * np.exp(-r / self.rho) - C6 * r **-6 

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
    def __init__(self,
                 labelHM     = None,
                 number      = None,
                 elementList = []):
        self.labelHM     = labelHM
        self.number      = number
        self.elementList = []

class Structure:
    def __init__(self,
                 unitCell      = None,
                 speciesList   = [],
                 symmetryGroup = None):
        self.unitCell      = unitCell
        self.speciesList   = speciesList
        self.symmetryGroup = symmetryGroup

    def uniqueSpecies(self, attrList):
        ''' Return list, worry about changes if need be'''
        from setTools import subsetByAttributes
        tempList = []
        for x in self.speciesList:
            if not subsetByAttributes([x], tempList, attrList):
                tempList.append(x)
        return tempList

    def setFracCoord(self):
        ''' CHECK CORRECT WAY AROUND!!!!!!!!!!!!!!
            assumes have lattice (vectors) set '''
        if not self.unitCell.invVectors:
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
        if self.speciesList[0].cartCoord is not None:
            self.setCartCoord()            
        print "Please set good tests for all these- consider different programs and their conventions"

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

    def setCharges(self, templateSpecies):
        ''' Pass in some templateSpecies- anything in speciesList that matches these gets the template's charge  '''
        from setTools import sameElementByAttributes
        for t in templateSpecies:
            for xs, s in enumerate(self.speciesList):
                if sameElementByAttributes(t, s, ['element', 'core']):
                    self.speciesList[xs].charge = t.charge
        return sum(x.charge for x in self.speciesList)

    def setASEMasses(self, templateSpecies):
        ''' Pass in some templateSpecies- anything in speciesList that matches these gets mass from ase.data  '''
        from setTools import sameElementByAttributes
        from ase.data import atomic_masses, chemical_symbols
        for t in templateSpecies:
            for xs, s in enumerate(self.speciesList):
                if sameElementByAttributes(t, s, ['element', 'core']):
                    self.speciesList[xs].mass = atomic_masses[chemical_symbols.index(t.element)]
        return sum(x.charge for x in self.speciesList)

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

#    @classmethod
    def inputCIF(self, cifName):
        from cifIO import getUnitCellCIF, getSpeciesListCIF, getSymmetryGroupCIF
        self.unitCell      = getUnitCellCIF     (cifName)
        self.speciesList   = getSpeciesListCIF  (cifName)
        self.symmetryGroup = getSymmetryGroupCIF(cifName)

    def removeSpeciesSobol(self, species, nRemove, offset = 0):
        ''' Remove nRemove species, that match character of species, according to closest fractional 
            coordinates to those from sobol list 
            Will only work well if atoms are sensibly bounded within lattice vectors (as looks in [0,1)^3 . LatticeVectors)
            Deletions are in order '''

        from setTools import sameElementByAttributes
        from sobol_lib import i4_sobol
   
        #Assume cartesian coordinates are wanted at this stage- easy to change if needed
        cartesian = True
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

        for i in deletedElements:
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
