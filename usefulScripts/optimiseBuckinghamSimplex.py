from coreClasses import Structure, Species
from vaspIO import VaspXML
from cifIO import writeCIFFromSpeciesList
from savedPotentials import currentPotentialSet141016
from cifIO import quickCifString, writeCIF
from gulpIO import InputFileGULP, OutputFileGULP
from hardcode import GULP_EXE
from subprocessHandling import RunCommand

import copy
import numpy as np
from scipy.optimize import fmin

xstalLTPExp = Structure().fromCIF('/Users/cmdc2-extra/Work/icsdFiles/EntryWithCollCode95979.cif', expandFullCell = True)

# expand to P1 cell- makes easier to compare and allows sym breaking to true minimum

qOCore = -2.86

xstalLTPExp.addShells(Species(element = 'O',  charge = qOCore))
xstalLTPExp.setCharges([Species(element = 'Ti', core='core', charge = 4.),
                        Species(element = 'Al', core='core', charge = 3.),
                        Species(element = 'Li', core='core', charge = 1.),
                        Species(element = 'P',  core='core', charge = 5.),
                        Species(element = 'O',  core='core', charge = -2 - qOCore)])

xstalLTPExp.resetFracCoord()

meanDisp = []
el = 'Ti'

potentialTiIndex = [x for x in enumerate(currentPotentialSet141016) if x[1].species1.element == 'Ti'][0][0]

print "index %s" % potentialTiIndex 

for si in xstalLTPExp.speciesMatchIndices(targetSpecies = Species(element = el)):
    print currentPotentialSet141016[potentialTiIndex].A,\
          currentPotentialSet141016[potentialTiIndex].rho, el, si, xstalLTPExp.nearestNeighbourDistance(speciesListIndex = si)
    meanDisp.append(xstalLTPExp.nearestNeighbourDistance(speciesListIndex = si))

oldTiODisp = np.mean(meanDisp)

print "oldTiODisp = %s" %oldTiODisp
#exit()

def attemptNewParametersTiO(inputList):
    TiOA, TiOrho = inputList[0], inputList[1]

    potentials = currentPotentialSet141016
    potentials[potentialTiIndex].A = TiOA
    potentials[potentialTiIndex].rho = TiOrho
    inFile = InputFileGULP(fileName = '',
                           title    = 'use to fit buckingham',
                           parentStructure = xstalLTPExp,
                           #                       keywords   = ['conp', 'opti', 'comp', 'phonon'],                                                                                                                                                                               
                           keywords   = ['conp', 'opti', 'phonon'],
                           potentials = potentials)

    #inFile.writeFile('LTP_opt_fit.gin')

    try:
        runner = RunCommand(GULP_EXE, inFile.stringForm())
        runner.run(timeout=10)

        with open('temp261016.gout', 'w') as outf:
            outf.write(runner.output)
        xstalLTPGULP = Structure().fromGULPOutput(runner.output)
        el = 'Ti'
        meanDisp = []
        for si in xstalLTPGULP.speciesMatchIndices(targetSpecies = Species(element = el)):
            print potentials[potentialTiIndex].A, potentials[potentialTiIndex].rho, el, si, xstalLTPGULP.nearestNeighbourDistance(speciesListIndex = si)
            meanDisp.append(xstalLTPGULP.nearestNeighbourDistance(speciesListIndex = si))
        assert(len(meanDisp) == 12)
        return abs(np.mean(meanDisp) - oldTiODisp)

    except:
        return np.inf



print fmin(attemptNewParametersTiO, ( currentPotentialSet141016[potentialTiIndex].A,  currentPotentialSet141016[potentialTiIndex].rho))

#writeCIF('xstalLTPGULP.cif', xstalLTPGULP)

#print 'gulp LTP'
#print xstalLTPGULP.unitCell.lengths
#for el in ['Al', 'P', 'Ti', 'Li']:
#
