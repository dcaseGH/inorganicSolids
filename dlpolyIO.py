import numpy as np
import fortranformat as ff

'''
imcon meaning
0 no periodic boundaries
1 cubic boundary conditions
2 orthorhombic boundary conditions
3 parallelepiped boundary conditions
4 truncated octahedral boundary conditions
5 rhombic dodecahedral boundary conditions
6 x-y parallelogram boundary conditions with no periodicity in the z direction
7 hexagonal prism boundary conditions
'''


class dlpolyInput():
    ''' see dlpoly classic documentation
        levcfg = 
        imcom  = 
        natms can be left None, as number atoms calculable
        engcfg is energy-- what to do ???? '''
    def __init__(self,
                 parentStructure    = None,
                 aseStructure       = None,
                 fileName           = None,
                 title              = None,
                 unitCell           = None,
                 levcfg             = None,
                 imcon              = None,
                 natms              = None,
                 engcfg             = None,
                 units              = 'eV',
                 temperature        = None,
                 integrator         = None,
#                 keywords           = [],
                 potentials         = [],
                 moleculeOrder      = None,
#                 defectSpecies      = [],
#                 defectDetails      = {},
#                 outputInstructions = {}
                 ):
        self.fileName          = fileName
        self.title             = title
        self.natms             = natms
        self.engcfg            = engcfg
        self.imcon             = imcon
        self.levcfg            = levcfg
        self.parentStructure   = parentStructure
        self.aseStructure      = aseStructure
        self.potentials        = potentials
        self.moleculeOrder     = moleculeOrder
        self.units             = units
        self.temperature       = temperature
        self.integrator        = integrator
        #these are just pointers (i think) to make it quicker to get information
        if self.parentStructure:
            self.unitCell      = self.parentStructure.unitCell
            self.speciesList   = self.parentStructure.speciesList
            self.symmetryGroup = self.parentStructure.symmetryGroup
        else:
            self.unitCell      = None
            self.speciesList   = None
            self.symmetryGroup = None
        # can specify unitCell independently (why?????)
#        if unitCell:
#            self.unitCell  = unitCell

    def numberFormatter(self, q):
        ''' Always return string- perhaps reformat strings later  '''
        if type(q) == float:
            if abs(q) < 1.e-5:
                return '{:20}'.format(0.0)
            return '{:20}'.format(q)#.replace('e', 'E')
        elif type(q) == int:
            return '{:^10}'.format(q)
        else:
            return '{:8}'.format(q)

#    @classmethod
    def initFromASE(self, 
                    cifFile   = None,
                    superCell = None):
        ''' Atm initialize an instance of dlployInput, then run this to fill in details
            supercell is tuple of ints, i.e. (1,1,1)
            ASE has lots of stuff- e.g. velocity distn etc, maybe use later '''
        from coreClasses import Species, Structure, UnitCell
        # consider issues wiht supercell... primative cells etc... keep simple
        from ase.io import read as ASERead
        #self.parentStructure should be a core Structure instance, and ASE structure should be some
        self.aseStructure = ASERead(cifFile)
#        self.parentStructure = ASERead(cifFile)
        if superCell:
#            self.parentStructure = self.parentStructure * superCell
            self.aseStructure = self.aseStructure * superCell
#        self.parentStructure = Structure(unitCell    = self.aseStructure._cell,
#                                         speciesList = [Species.initFromASEAtom(x) for x in self.aseStructure])
        self.parentStructure = Structure(unitCell    = UnitCell(vectors = self.aseStructure._cell),
                                         speciesList = [Species.initFromASEAtom(x) for x in self.aseStructure])

        # these should just be pointers.. should still get functionality of structure????
        # CHANGE UNIT CELL TO BE UNIT CELL VECTORS
        self.unitCell    = self.parentStructure.unitCell
        self.speciesList = self.parentStructure.speciesList
        return

    @classmethod
    def initFromCONFIG(cls, configFile, levcfg = None):
        ''' Make an input from REVCON, for e.g. because you want to change things '''
        from coreClasses import Structure, Species, UnitCell
        structure = Structure()

        def atomsInConfig(linesWithAtoms, levcfgIn, levcfg, natms):
            ''' possible that levcfgIn and levcfg differ, but shouldn't matter
                only take information needed '''

            outList = []

            for iat in xrange(natms):
                newAt = Species(element = linesWithAtoms[iat * (levcfgIn + 2)].split()[0].split('_')[0])
                if 'sh' in linesWithAtoms[iat * (levcfgIn + 2)].split()[0]:
                    newAt.core = 'shel'
                newAt.cartCoord = np.array(map(float, linesWithAtoms[iat * (levcfgIn + 2) + 1].split()))
                if levcfgIn > 0 and levcfg > 0:
                    newAt.cartVelocity = np.array(map(float, linesWithAtoms[iat * (levcfgIn + 2) + 2].split()))
                if levcfgIn > 1 and levcfg > 1:
                    newAt.cartVelocity = np.array(map(float, linesWithAtoms[iat * (levcfgIn + 2) + 3].split()))
                outList.append(newAt)

            return outList

        with open(configFile, 'r') as cFile:
            lines = cFile.readlines()
            title = lines[0]
            levcfgIn = int(lines[1].split()[0])
            imcon    = int(lines[1].split()[1])
            natms    = int(lines[1].split()[2])
            engcfg   = float(lines[1].split()[3])

            if levcfg is None:
                levcfg = levcfgIn

            structure.unitCell = UnitCell(vectors = np.array([map(float, lines[2].split()),
                                                              map(float, lines[3].split()),
                                                              map(float, lines[4].split())]))

            structure.speciesList = atomsInConfig(lines[5:], levcfgIn, levcfg, natms)

        return cls(title           = title,
                   parentStructure = structure,
                   levcfg          = levcfg,
                   imcon           = imcon,
                   engcfg          = engcfg)


    def atomicInformation(self, atom, atomIndex):
        ''' Prepare the information prior to calling this  '''
        aString = " ".join([self.numberFormatter(atom.dlpolyLabelStandard()), self.numberFormatter(atomIndex), "\n"])
        line = ff.FortranRecordWriter('(3f20.10)')
        aString +=  line.write(map(float, atom.cartCoord.tolist())) + "\n"
        if self.levcfg > 0:
            aString += line.write(map(float, atom.cartVelocity.tolist())) + "\n"
        if self.levcfg > 1:
            aString += " ".join(list(atom.cartForce) + ["\n"])
        return aString

    def checkCorrectOrdering(self):
        ''' DL_POLY is very fussy about the order
            at this stage assume molecules = atoms and
            all input has been turned into correct format (self.speciesList contains Atoms etc)
            Run this before making any output  '''
        #only ever set this if confident that self.speciesList is ordered
        #try this at some point: all(l[i] <= l[i+1] for i in xrange(len(l)-1)), or use index in self.moleculeOrder in lambda fn
        if self.moleculeOrder:
            return True

        #cartesian cutoff- possible risk if ever put two atoms too close
        displacementCutoff = 1.e-2

        #sort alphabetically by element
        self.moleculeOrder = sorted(set([x.element for x in self.speciesList]))
        shellsList =  sorted(set([x.element for x in self.speciesList if x.core[:4] == 'shel']))
        print 'checking molecule order -> ' , self.moleculeOrder
        #put shells next to correct atom
        self.speciesList = sorted([x for x in self.speciesList], key = lambda x: x.element)

        newSpeciesList  = []
        # this loop rearranges so that shells are next to correct core
        # beware moving inside loop
        for iat1, at1 in enumerate(self.speciesList):
            if at1.core == 'core':
                newSpeciesList.append(at1)
                if at1.element in shellsList:
                    for iat2, at2 in enumerate(self.speciesList):
                        if at2.core[:4] == 'shel' and at2.element == at1.element and \
                           np.linalg.norm(at1.cartCoord - at2.cartCoord) < displacementCutoff:
                            newSpeciesList.append(at2) 
#                else:
#                    continue
        if len(newSpeciesList) != len(self.speciesList):
            print "NOT FOUND ALL SPECIES IN REORDERING", len(newSpeciesList) != len(self.speciesList);exit()
        self.speciesList = newSpeciesList

        def rearrangeBondedList(previousSpeciesList, centralElement, adjacentElement, bondCutoff):
            ''' This rearranges list so that central element and all the atoms bonded to it are at the top
                e.g. this is used to take PO4 units and rearrange them 
                Assume that the above business of moving order wrt core/shel is complete '''
            bondedList, nonBondedList = [], []
            for iat1, at1 in enumerate(previousSpeciesList):
                if at1.element == centralElement:
                    bondedList.append(at1)#previousSpeciesList[iat1])
                    for iat2, at2 in enumerate(previousSpeciesList):
                        if at2.element == adjacentElement and \
                           np.linalg.norm(at1.cartCoord - at2.cartCoord) < bondCutoff:
                            bondedList.append(at2)#previousSpeciesList[iat2])
                elif at1.element != adjacentElement:
                    nonBondedList.append(at1)#previousSpeciesList[iat1])

            if len(previousSpeciesList) != len(bondedList) + len(nonBondedList):
                raise ValueError("Error rearrangeBondedList %s ne %s %s")%(len(previousSpeciesList), len(bondedList), len(nonBondedList))

            return bondedList + nonBondedList

        def movedBondedAtomsWRTPeriodicBoundaries(previousSpeciesList, unitCell, centralElement, adjacentElement, bondCutoff):
            outList = []
            displacements = []
            for i1 in xrange(-1, 2):
                for i2 in xrange(-1, 2):
                    for i3 in xrange(-1, 2):
                        displacements.append(np.dot(np.array((i1, i2, i3)), unitCell))

            for iat1, at1 in enumerate(previousSpeciesList):
                if at1.element != adjacentElement:
                    outList.append(previousSpeciesList[iat1])
                for iat2, at2 in enumerate([x for x in previousSpeciesList if x.element == centralElement]):
                    for d in displacements:
                        if np.linalg.norm(at1.cartCoord + d + at2.cartCoord) < bondCutoff:
                            previousSpeciesList[iat1].cartCoord += d
                            outList.append(previousSpeciesList[iat1])
                            break
            return outList
        return True

    def stringFormCONFIG(self):
        ''' Maintain rubbish fortran names throughout for consistency
            This file is just the cell, atoms, periodicity etc  '''

        self.checkCorrectOrdering()

        inString = str(self.title) + "\n"

        # basic commands
        if not self.natms:     
            self.natms = len([x for x in self.speciesList if x.core == 'core'])
        inString += " ".join(map(self.numberFormatter, [self.levcfg, self.imcon, self.natms, self.engcfg, "\n"])).replace("None", "")

        # print unit cell
        formatter = ff.FortranRecordWriter('(3f20.10)')
        for i in xrange(3):
#            inString += " ".join(map(self.numberFormatter, self.unitCell[i].tolist())) + "\n"
#            inString += " ".join(map(self.numberFormatter, self.unitCell.vectors[i].tolist())) + "\n"
            inString += formatter.write(self.unitCell.vectors[i].tolist()) + "\n"

        for iat, atom in enumerate(self.speciesList):
            if self.levcfg > 0 and atom.cartVelocity is None:
                self.speciesList[iat].setThermalVelocity(self.temperature)
            inString += self.atomicInformation(self.speciesList[iat], iat + 1)

        return inString

    def molecularInformationFIELD(self):
        ''' Writes basic stuff about molecules
            ASSUME THAT MOLECULES ARE JUST ATOMS '''

        from coreClasses import VSpring
        springPots   = [x for x in self.potentials if isinstance(x, VSpring)]
        def correctSpringK(pots, at1, at2):
            ''' should probably generalize to all potentials and move to core classes?  '''
            for p in pots:
                if p.species1.element == at1.element and p.species1.element == at2.element:
                    return p.K
            #raise error if not in list?
            return None

        outString = "molecular types %s\n" % len(self.moleculeOrder)
        uniqueSpeciesList = [x for x in self.parentStructure.uniqueSpecies(['element', 'core'])]# if x.core == 'core']
        for a in self.moleculeOrder:
            aAsAtom = [x for x in uniqueSpeciesList if x.element == a and x.core == 'core'][0]
            outString += a + "\n"
            outString += "nummols %s\n" %(len([x for x in self.speciesList if x.element == a and x.core == 'core']))

            if any([x.element == a and x.core[:4] == 'shel' for x in uniqueSpeciesList]):
                outString += "atoms 2\n"
                outString += " ".join([aAsAtom.dlpolyLabelStandard(), str(aAsAtom.mass), str(aAsAtom.charge), "\n"])
                shellAtom = [x for x in uniqueSpeciesList if x.element == a and x.core[:4] == 'shel'][0]
                outString += " ".join([shellAtom.dlpolyLabelStandard(), str(shellAtom.mass), str(shellAtom.charge),  "\n"])
                # two units in shell, model = 1 (adiabatic (cheap) model)
                outString += "shell 1 1\n"
                #parts 1 and two of the molecule are interacting with shell strength K
                outString += "1 2 %s\n"%(correctSpringK(springPots, aAsAtom, shellAtom))
            else:
                outString += "atoms 1\n"
                #sitnam (site) name, weight (atomic site mass), chge [charge], nrept (repeat counter), ifrz (frozen), igrp (neutral group charge)????
                outString += " ".join([aAsAtom.dlpolyLabelStandard(), str(aAsAtom.mass), str(aAsAtom.charge), "\n"])
            outString += "finish\n"
        return outString

    def vdwInformationFIELD(self):
        ''' If you give superfluous potentials, perhaps DLPOLY wont like it???   '''
        from coreClasses import VBuckingham
        twoBodyPots   = [x for x in self.potentials if isinstance(x, VBuckingham)]

        outString = "vdW %s\n"%len(twoBodyPots)
        for p in twoBodyPots:
            outString += p.dlpolyString() + "\n"
        return outString

    #cant decorate as don't know self... class method??  
#    @checkCorrectOrdering(self)
    def stringFormFIELD(self):
        ''' Stuff to do with force-fields/energy etc  '''
        # call this here if can't decorate
        self.checkCorrectOrdering()
        inString  = self.title + "\n"
        inString += " ".join(['Units', str(self.units), "\n"])
        #inString += something to do with neutral/charged stuff... leave this for now

        # do this because already checked moleculeOrder        
        inString += self.molecularInformationFIELD()

        inString += self.vdwInformationFIELD()

        inString += "close\n"
        return inString       

    def stringFormCONTROL(self):
        '''   '''
        outString  = self.title + "\n\n"
        outString += "integration " + self.integrator + "\n\n"
        outString += "temperature " + str(self.temperature) + "\n"
        outString += "pressure " + str(self.pressure) + "\n"
        outString += "ensemble " + self.ensemble + "\n"
        outString += "\n"
        outString += "\n"
        outString += "\n"
        outString += "\n"
        outString += "\n"
        outString += "\n"

        return outString + 'finish'

class dlpolyJob():
    ''' Has information to run one DL_POLY classic job
        Mainly needs CONTROL, FIELD and CONFIG somehow- either reads the string or learns file name
        time in seconds (int or float??)
        Either give a file (address of) or string 
        Can technically set dlpolyInput to be a blank thing and just copy files around, but must set a time'''
    def __init__(self,
                 dlpolyInput   = None,
                 configString  = None,
                 controlString = None,
                 fieldString   = None,
                 configFile    = None,
                 controlFile   = None,
                 fieldFile     = None,
                 time          = None,
                ):

        if dlpolyInput and time:
            self.time        = float(time)
            self.dlpolyInput = dlpolyInput
        else:
            return False

        if configString:
            self.configString  = configString
        elif configFile:
            self.configFile    = configFile
        else:
            self.configString  = self.dlpolyInput.stringFormCONFIG()

        if fieldString:
            self.fieldString   = fieldString
        elif fieldFile:
            self.fieldFile     = fieldFile
        else:
            self.fieldString   = self.dlpolyInput.stringFormFIELD()

        if controlString:
            self.controlString = controltring
        elif controlFile:
            self.controlFile   = controlFile
        else:
            self.controlString = self.dlpolyInput.stringFormCONTROL()

        return True



class dlpolyOutput():
    def __init__(self,
                 inputInstance = None,
                 outputDF      = None,
                 outputString  = None
                ):
        self.inputInstance = inputInstance
        self.outputDF      = outputDF
        self.outputString  = outputString

    def reformat(self):
        ''' Need something to reformat numbers/units stuff '''
        pass

    def dataFrameOutput(self, *args, **kwargs):#outputString = None):
        import pandas as pd
        from hardcode import timeStringToFloat
        ''' Generate a pandas data frame from OUTPUT
            args can be in any order- pd will make a dictionary like thing 
            outputString must be the actual string, so open('OUTPUT', 'r').read() 
            All args are the outputs in the header--- they are floats
        '''

        errorDetected = False
        if not self.outputString:
            self.outputString = kwargs['outputString']

        resultsDict = {}        

        # time variables are used to ID each data frame- also record average as we go
        timeVariables = ['step', 'time', 'cpu-time']
        timeFormat    = [int, timeStringToFloat, timeStringToFloat]

        for k in timeVariables:
            resultsDict[k] = []

        for k in args:
            resultsDict[k] = []
            resultsDict[k + 'Average'] = []

        header = None
        dataLines = self.outputString.split("\n")
        for il in xrange(len(dataLines)):

            if 'run terminating' in dataLines[il]:
                break

            # formatting of OUTPUT is to make a header, and then ordered table below- get indices for table
            if header is None and '--------------------------------' in dataLines[il]:
                header = [dataLines[il + il2].split() for il2 in xrange(2, 5)]
                # replace cpu time for cpu_time- everything needs _ otherwise -> 2 strings
                header[2][1] = 'cpu_time'            
                # the first element was 'cpu' so delete it
                del(header[2][0])
                header = np.array(header, dtype=str)
                headerIndices = {}
                for x in args:
                    headerIndices[x] = zip(*np.where(header == x))[0]
#                headerIndices = [zip(*np.where(header == x)) for x in args]

            if isinstance(header, np.ndarray) and '--------------------------------' in dataLines[il]:

                # sometimes the header is repeated, so skip this (assume still the same)
                if '--------------------------------' in dataLines[il + 1]:
                    continue
                if 'step' in dataLines[il + 2]:
                    continue
                if 'run terminating' in dataLines[il + 3]:
                    break

                try:
                    # follow same thing as above when making header
                    dataInstant = np.array([dataLines[il + il2].split() for il2 in xrange(1, 4)], dtype=str)
                    dataAverage = [dataLines[il + il2].split() for il2 in xrange(5, 8)]
                    dataAverage[2].insert(0, 'None')
                    dataAverage = np.array(dataAverage, dtype=str)

                    # assume that timeVariables ordered as defined
                    for ik, k in enumerate(timeVariables):
                        resultsDict[k].append( timeFormat[ik]( dataInstant[ik, 0] ) )
                    for i in args:
                        resultsDict[i].append( float( dataInstant[ headerIndices[i] ] ) )
                        resultsDict[i + 'Average'].append( float( dataAverage[ headerIndices[i] ] ) )
                except:
                    errorDetected = True

        if errorDetected:
            print "Latent error detected, perhaps not all data good/available"

        self.outputDF = pd.DataFrame(resultsDict)
        return self.outputDF
