from coreClasses import Structure
import numpy as np

class OutputFileGULP():
    def __init__(self,
                 fileName        = None):
        self.fileName = fileName

    #def makeOutputStructure()
    @classmethod
    def latticeEnergy(self, inputString):
        return

    @classmethod
    def readOutputCompOpti(self, inputString):
        ''' If you run comp and opti it will make a table of changes - read this '''
        stringTable = inputString.split('Parameter   Initial value   Final value   Difference    Units      Percent')[1].split('Volume')[1]
        newCoords  = []
        newLengths = []
        newAngles  = []
        rows = stringTable.split('\n')       

        for r in rows:
            if 'Fractional' in r:
                newCoords.append( float(r.split()[3]))
            if 'Angstroms'  in r:
                newLengths.append(float(r.split()[2]))
            if 'Degrees'    in r:
                newAngles.append( float(r.split()[2]))

        return {'newCoords' : np.array(newCoords).reshape((len(newCoords) / 3, 3)),
                'newAngles' : np.array(newAngles),
                'newLengths': np.array(newLengths)}

    @classmethod
    def phononFrequencies(self, inputString):
        phonString = inputString.split('Frequencies (cm-1) [NB: Negative implies an imaginary mode]:\n\n')[1].split('\n\n\n')[0]
        return np.array(map( float, phonString.replace('\n', '').split() ))

    @classmethod
    def defectEnergy(self, inputString):
#        return float(inputString.split('Total defect energy ')[1].split('eV')[0].replace('=', ''))
        defectEnergies = []
        for l in inputString.split("\n"):
            if 'Total defect energy' in l:
                defectEnergies.append(float(l.split('Total defect energy ')[1].split('eV')[0].replace('=', '')))
        return defectEnergies                


class InputFileGULP(Structure):
    # see gulp webpage or http://www.nsccs.ac.uk/si_gulp.php for information on file type                                                        
    ''' '''
    def __init__(self,
                 parentStructure    = None,
                 fileName           = None,
                 title              = None,
                 keywords           = [],
                 potentials         = [],
                 defectSpecies      = [],
                 defectDetails      = {},
                 outputInstructions = {}
                 ):
        self.fileName      = fileName
        self.title         = title
        if parentStructure:
            self.unitCell      = parentStructure.unitCell
            self.speciesList   = parentStructure.speciesList
            self.symmetryGroup = parentStructure.symmetryGroup
        else:
            self.unitCell      = None
            self.speciesList   = None
            self.symmetryGroup = None
        self.keywords      = keywords
        self.potentials    = potentials
        self.defectSpecies = defectSpecies
        self.defectDetails = defectDetails
        self.outputInstructions = outputInstructions

    def potentialInputString(self, matchAttributes = ['element', 'core']):

        from itertools   import combinations
        from setTools    import sameSetByAttributes
        from coreClasses import VBuckingham, VSpring, VThreeBody

        outString = ""
        twoBodyPots   = [x for x in self.potentials if isinstance(x, VBuckingham)]
        springs       = [x for x in self.potentials if isinstance(x, VSpring)]
        threeBodyPots = [x for x in self.potentials if isinstance(x, VThreeBody)]
        
        # consider unique pairs and also unique items self-interaction
#        for pair in combinations(self.uniqueSpecies(matchAttributes), 2):
#            for pot in twoBodyPots:
##            for pot in self.potentials:
#                if sameSetByAttributes(pair, (pot.species1, pot.species2), matchAttributes):
#                    outString += pot.stringForm() + "\n"
#        for pair in [(x, x) for x in self.uniqueSpecies(matchAttributes)]:
##            print pair[0].element, pair[1].element
#            for pot in twoBodyPots:
##            for pot in self.potentials:
#                if sameSetByAttributes(pair, (pot.species1, pot.species2), matchAttributes):
#                    outString += pot.stringForm() + "\n"
        #these are hacked- only add needed ones

        for p in twoBodyPots:
            outString += p.stringForm() + "\n"
        for p in springs:
            outString += p.stringForm() + "\n"
        for p in threeBodyPots:
            outString += p.stringForm() + "\n"

        return outString

    def defectsInputString(self, coordinates = 'frac'):
        # one mad day may want to give cart coords
        outString = ""

        for d in self.defectSpecies:
            outString += d.defectType + "\n"
            outString += " ".join(map(str, [d.species.element] + list(d.species.fracCoord))) + "\n"

        outString += "centre\n"
        outString += " ".join(map(str, self.defectDetails['centre'])) + "\n"
        outString += "size\n"
        outString += str(self.defectDetails['innerRadius']) + " " + str(self.defectDetails['outerRadius']) + "\n"

        return outString

    def outputInstructionsString(self):
        return "\n".join(["output " + k + " " + self.outputInstructions[k] for k in self.outputInstructions])

    def stringForm(self, coordinates = 'frac'):
        ''' return a string which looks like an input file
            please note outputInstructions DOES NOT RELATE TO .GOUT
            for this, direct the output of however you run the calcn (this if for cif/xyz movies '''

        inString  = " ".join([x for x in self.keywords]) + "\n"
        if self.title:
            inString += "title \n %s \nend \n" %(self.title)
        inString += "cell\n"
        if self.unitCell:
            inString += " ".join([str(x) for x in self.unitCell.lengths]) + " "
            inString += " ".join([str(x) for x in self.unitCell.angles]) + "\n"

        if self.speciesList and coordinates.lower() == 'frac':
            inString += "fractional\n"
            inString += "\n".join([" ".join([a.element] + [a.core] + [str(x) for x in a.fracCoord]) for a in self.speciesList]) + "\n"

        if self.symmetryGroup and self.symmetryGroup.labelHM:
            inString += "space\n"
            # the labels from cifs also include H for hexagonal on the end- gulp doesn't like this (I'm not 100pc sure)
            # changed by dhc 110816 to remove the S at the end of some things
            if self.symmetryGroup.labelHM[-1] in ['H']:
                inString += self.symmetryGroup.labelHM[:-1] + "\n"
            elif self.symmetryGroup.labelHM[-1] in ['S']:
                inString += self.symmetryGroup.labelHM[:-1] + "\n"
            else:
                inString += self.symmetryGroup.labelHM + "\n"

        if self.speciesList:
            inString += "species\n"
            # technically should be getting potential label not element
            inString += "\n".join([" ".join([a.element] + [a.core] + [str(a.charge)]) for a in self.uniqueSpecies(['element', 'core', 'charge'])]) + "\n"

        if any([x.lower()[:4] == 'defe' for x in self.keywords]):
            inString += self.defectsInputString()
        inString += self.potentialInputString()

        if len(self.outputInstructions.keys()) > 0:
            inString += self.outputInstructionsString()

        return inString

    def writeFile(self, fileName):
        with open(fileName, 'w') as outf:
            outf.write(self.stringForm())