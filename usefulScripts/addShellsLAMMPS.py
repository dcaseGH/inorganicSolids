''' Script to read final output of a LAMMPS run with no shells, and add the shells '''
# ONE DAY DO THIS    run as python [thisScript.py] [input data file name] [output data file name] '''

from lammpsIO import LAMMPSOutput, LAMMPSInput
from coreClasses import Species

import sys

# for now all inputs are hard-coded into this script
inputFile   = sys.argv[1]#''
outputFile  = sys.argv[2]#''
speciesDict = {1: Species(element='B'),
               2: Species(element='Li'),
               3: Species(element='O'),
               4: Species(element='V')}
# make lists work if ever needed
shellList = ['O']
qOShell   = -2.86

structure = LAMMPSOutput.structureFromDataFile(inputFile, speciesDict = speciesDict) 
print "Input yields structure with %s species"%(len(structure.speciesList))

for element in shellList:
    structure.addShells(Species(element='O', charge = qOShell), cartesian=True)

assert(abs(structure.setCharges([Species(element = 'V',  charge = 5.),
                                 Species(element = 'Li', charge = 1.),
                                 Species(element = 'B',  charge = 3.),
                                 Species(element = element,  charge = -2. - qOShell)]
                                )) < 0.001)

#add standard masses (what I always use)
structure.setASEMasses(structure.uniqueSpecies(['element']))
structure.changeMassesForShells(Species(element='O'), massShell=1.)

lammpsInput = LAMMPSInput(parentStructure = structure)
with open(outputFile, 'w') as outf:
    outf.write(lammpsInput.dataString())
