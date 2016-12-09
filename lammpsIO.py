import numpy as np


class StageInstructions():
    ''' Put these in stageList of LAMMPSInput '''


class LAMMPSInput():
    ''' parentStructure is the input
        stageList contains instructions for each step '''

    def __init__(self,
                 parentStructure    = None,
                 fileSuffix         = 'Null',
                 title              = 'Null',
                 stageList          = [],
                 uniqueSpeciesOrder = None
                 ):

        self.fileSuffix         = fileSuffix
        self.title              = title
        self.parentStructure    = parentStructure
        self.stageList          = stageList
        self.uniqueSpeciesOrder = uniqueSpeciesOrder

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

    def dataString(self, blurb='data file for LAMMPS'):
        ''' String for the data file '''
        outString  = blurb + "\n\n"
        outString += str(len(self.parentStructure.speciesList)) + "   atoms\n"
        outString += str(self.numberBonds()) + "   bonds\n"
        outString += str(self.numberAngles()) + "   angles\n"
        outString += str(self.numberDihedrals()) + "   dihedrals\n\n"

        if not self.uniqueSpeciesOrder:
            self.defineUniqueSpeciesOrder()

        outString += str(len(self.uniqueSpeciesOrder)) + "   atom types\n"
        # for now just stick to simple atoms and shells- ignore molecules!
        outString += str(len([x for x in self.uniqueSpeciesOrder if x.core[0].lower() == 's'])) + "   bond types\n" 
        outString += "0   angle types\n"
        outString += "0   dihedral types\n\n"

        outString += "UNIT CELL\n\n"

        assert(self.uniqueSpeciesOrder[0].mass > 0.)
        outString += "Masses\n\n"
        for i, at in enumerate(self.uniqueSpeciesOrder):
            outString += str(i) + " " + at.mass + "\n"


        return outString
