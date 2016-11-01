from coreClasses import Structure, Species
from vaspIO import VaspXML
from cifIO import quickCifString, writeCIF
from gulpIO import InputFileGULP, OutputFileGULP
from hardcode import GULP_EXE
from subprocessHandling import RunCommand

from cifIO import writeCIFFromSpeciesList

import numpy as np
import copy

from savedPotentials import currentPotentialSet141016, currentPotentialSet271016

#pots = currentPotentialSet141016
pots = currentPotentialSet271016
qOCore = -2.86

for i in ['0']:#map(str, xrange(6)):
    
    LATPDict = {}
#    xstalFF   = Structure().fromCIF('LATP_final_' + i + '.cif')
    xstalPreFF = Structure().fromCIF('LATP_initial_'+ i + '.cif')
    xstalVasp = VaspXML().finalStructure(open('vasprun_LATP_' + i + '.xml', 'r').read())
    xstalExp  = Structure().fromCIF('EntryWithCollCode427619_hacked.cif')

    LATPDict['A'], LATPDict['B'], LATPDict['C'] = [None, None, None], [None, None, None], [None, None, None]
    LATPDict['Ti'], LATPDict['Li'], LATPDict['Li_{i}'], LATPDict['Al'], LATPDict['P'] = [None, None, None], [None, None, None], [None, None, None], [None, None, None], [None, None, None]

    xstalPreFF.addShells(Species(element = 'O',  charge = qOCore))
    xstalPreFF.setCharges([Species(element = 'Ti', core='core', charge = 4.),
                           Species(element = 'Al', core='core', charge = 3.),
                           Species(element = 'Li', core='core', charge = 1.),
                           Species(element = 'P',  core='core', charge = 5.),
                           Species(element = 'O',  core='core', charge = -2 - qOCore)])
    xstalPreFF.resetFracCoord()

    inFile = InputFileGULP(fileName = '',
                       title    = 'use to fit buckingham',
                       parentStructure = xstalPreFF,
#                       keywords   = ['conp', 'opti', 'phonon'],
                       keywords   = ['conp', 'opti'],
                       potentials = pots)

    from subprocessHandling import RunCommandSafe
    runner = RunCommandSafe(GULP_EXE, inFile.stringForm())
    runner.run()

    xstalFF = Structure().fromGULPOutput(runner.output)
    writeCIF('LATP_0_FFmin.cif', xstalFF)

    print 'xstalFF'
    print xstalFF.unitCell.lengths
    LATPDict['A'][0], LATPDict['B'][0], LATPDict['C'][0] = xstalFF.unitCell.lengths[0], xstalFF.unitCell.lengths[1], xstalFF.unitCell.lengths[2]
    for el in ['Al', 'P', 'Ti', 'Li']:
        meanList = []
        for si in xstalFF.speciesMatchIndices(targetSpecies = Species(element = el)):
            print el, si, xstalFF.nearestNeighbourDistance(speciesListIndex = si)
            if si != 108:
                meanList.append(xstalFF.nearestNeighbourDistance(speciesListIndex = si))
            else:
                LATPDict['Li_{i}'][0] = xstalFF.nearestNeighbourDistance(speciesListIndex = si)
        print np.mean(meanList)
        LATPDict[el][0] = np.mean(meanList)

    print 'xstalVasp'
    print xstalVasp.unitCell.lengths
    LATPDict['A'][1], LATPDict['B'][1], LATPDict['C'][1] = xstalVasp.unitCell.lengths[0], xstalVasp.unitCell.lengths[1], xstalVasp.unitCell.lengths[2]
    for el in ['Al', 'P', 'Ti', 'Li']:
        meanList = []
        for si in xstalVasp.speciesMatchIndices(targetSpecies = Species(element = el)):
            print el, si, xstalVasp.nearestNeighbourDistance(speciesListIndex = si)
            if si != 7:
                meanList.append(xstalVasp.nearestNeighbourDistance(speciesListIndex = si))
            else:
                LATPDict['Li_{i}'][1] = xstalVasp.nearestNeighbourDistance(speciesListIndex = si)
        print np.mean(meanList)
        LATPDict[el][1] = np.mean(meanList)

    print 'xstalExp'
    print xstalExp.unitCell.lengths
    LATPDict['A'][2], LATPDict['B'][2], LATPDict['C'][2] = xstalExp.unitCell.lengths[0], xstalExp.unitCell.lengths[1], xstalExp.unitCell.lengths[2]
    for el in ['Al', 'P', 'Ti', 'Li']:
        for si in xstalExp.speciesMatchIndices(targetSpecies = Species(element = el)):
            if si in [4, 5]:
                print el, si, xstalExp.radialDistributionFunction(speciesListIndex = si)[1]
                LATPDict[el][2] = xstalExp.radialDistributionFunction(speciesListIndex = si)[1]
            elif si == 6:
                LATPDict['Li_{i}'][2] = xstalExp.radialDistributionFunction(speciesListIndex = si)[1]
            else:
                print el, si, xstalExp.nearestNeighbourDistance(speciesListIndex = si)
                LATPDict[el][2] = xstalExp.nearestNeighbourDistance(speciesListIndex = si)

print LATPDict#;exit()

print 'Al2TiO5'
ATODict = {}

ATODict['A'], ATODict['B'], ATODict['C'] = [None, None, None], [None, None, None], [None, None, None]
ATODict['Ti'], ATODict['Al'] = [None, None, None], [None, None, None]

xstalATOVasp = VaspXML().finalStructure(open('AlTiO_0_K422.xml', 'r').read())
print 'vasp ATO'
print xstalATOVasp.unitCell.lengths
ATODict['A'][1], ATODict['B'][1], ATODict['C'][1] = xstalATOVasp.unitCell.lengths[0], xstalATOVasp.unitCell.lengths[1], xstalATOVasp.unitCell.lengths[2]
for el in ['Al', 'Ti']:
    meanList = []
    for si in xstalATOVasp.speciesMatchIndices(targetSpecies = Species(element = el)):
        print el, si, xstalATOVasp.nearestNeighbourDistance(speciesListIndex = si)
        meanList.append(xstalATOVasp.nearestNeighbourDistance(speciesListIndex = si))
    ATODict[el][1] = xstalATOVasp.nearestNeighbourDistance(speciesListIndex = si)

# Minimise ATO with ForceField

#xstalATOExp = VaspXML().initialStructure(open('AlTiO_0_K422.xml', 'r').read())#= Structure().fromCIF('')
xstalATOExp  = Structure().fromCIF('/Users/cmdc2-extra/Work/icsdFiles/MyBaseFileNameCollCode161695.cif')

writeCIF('xstalATOExp.cif', xstalATOExp)

print 'exp ATO'
print xstalATOExp.unitCell.lengths
ATODict['A'][2], ATODict['B'][2], ATODict['C'][2] = xstalATOExp.unitCell.lengths[0], xstalATOExp.unitCell.lengths[1], xstalATOExp.unitCell.lengths[2]
#xstalATOExp['A'][1], xstalATOVasp['B'][1], xstalATOVasp['C'][1] = xstalATOVasp.unitCell.lengths[0], xstalATOVasp.unitCell.lengths[1], xstalATOVasp.unitCell.lengths[2]
for el in ['Al', 'Ti']:
    meanList = []
    for si in xstalATOExp.speciesMatchIndices(targetSpecies = Species(element = el)):
        print el, si, xstalATOExp.nearestNeighbourDistance(speciesListIndex = si)
#        meanList.append(xstalATOExp.nearestNeighbourDistance(speciesListIndex = si))
        meanList.append(xstalATOExp.radialDistributionFunction(speciesListIndex = si)[1])
    ATODict[el][2] = np.mean(meanList)

#exit()

# should really deepcopy xstalATOVasp

#


xstalATOVasp.addShells(Species(element = 'O',  charge = qOCore))
xstalATOVasp.setCharges([Species(element = 'Ti', core='core', charge = 4.),
                         Species(element = 'Al', core='core', charge = 3.),
                         Species(element = 'O',  core='core', charge = -2 - qOCore)])
xstalATOVasp.resetFracCoord()

inFile = InputFileGULP(fileName = '',
                       title    = 'use to fit buckingham',
                       parentStructure = xstalATOVasp,
                       keywords   = ['conp', 'opti', 'comp', 'phonon'],
                       potentials = pots)#currentPotentialSet141016)

runner = RunCommand(GULP_EXE, inFile.stringForm())
runner.run(timeout=10)

xstalATOGULP = Structure().fromGULPOutput(runner.output)
with open('temp261016ATO.gout', 'w') as outf:
    outf.write(runner.output)

writeCIF('xstalGULP.cif', xstalATOGULP)

print 'gulp ATO'
print xstalATOGULP.unitCell.lengths
ATODict['A'][0], ATODict['B'][0], ATODict['C'][0] = xstalATOGULP.unitCell.lengths[0], xstalATOGULP.unitCell.lengths[1], xstalATOGULP.unitCell.lengths[2]
for el in ['Al', 'Ti']:
    meanList = []
    for si in xstalATOGULP.speciesMatchIndices(targetSpecies = Species(element = el)):
        print el, si, xstalATOGULP.nearestNeighbourDistance(speciesListIndex = si)
        meanList.append(xstalATOGULP.nearestNeighbourDistance(speciesListIndex = si))
    ATODict[el][0] = np.mean(meanList)

print "ATODict =", ATODict

LTPDict = {}
LTPDict['A'], LTPDict['B'], LTPDict['C'] = [None, None, None], [None, None, None], [None, None, None]
LTPDict['Ti'], LTPDict['Li'], LTPDict['P'] = [None, None, None], [None, None, None], [None, None, None]

# get it in a P1 cell- save trouble
xstalLTPExp = Structure().fromCIF('/Users/cmdc2-extra/Work/icsdFiles/EntryWithCollCode95979.cif', expandFullCell = True)
print 'xstalLTPExp'
print xstalLTPExp.unitCell.lengths
LTPDict['A'][2], LTPDict['B'][2], LTPDict['C'][2] = xstalLTPExp.unitCell.lengths[0], xstalLTPExp.unitCell.lengths[1], xstalLTPExp.unitCell.lengths[2]
for el in ['P', 'Ti', 'Li']:
    meanList = []
    for si in xstalLTPExp.speciesMatchIndices(targetSpecies = Species(element = el)):
        print el, si, xstalLTPExp.nearestNeighbourDistance(speciesListIndex = si)
        meanList.append(xstalLTPExp.nearestNeighbourDistance(speciesListIndex = si))
    print np.mean(meanList)
    LTPDict[el][2] = np.mean(meanList)
 
xstalLTPPreGulp = copy.deepcopy(xstalLTPExp)

# expand to P1 cell- makes easier to compare and allows sym breaking to true minimum
#xstalLTPPreGulp.expandAllSymmetry()

xstalLTPPreGulp.addShells(Species(element = 'O',  charge = qOCore))
xstalLTPPreGulp.setCharges([Species(element = 'Ti', core='core', charge = 4.),
                            Species(element = 'Al', core='core', charge = 3.),
                            Species(element = 'Li', core='core', charge = 1.),
                            Species(element = 'P',  core='core', charge = 5.),
                            Species(element = 'O',  core='core', charge = -2 - qOCore)])
xstalLTPPreGulp.resetFracCoord()

inFile = InputFileGULP(fileName = '',
                       title    = 'use to fit buckingham',
                       parentStructure = xstalLTPPreGulp,
#                       keywords   = ['conp', 'opti', 'comp', 'phonon'],
                       keywords   = ['conp', 'opti', 'phonon'],
                       potentials = pots)# currentPotentialSet141016)

inFile.writeFile('LTP_opt_fit.gin')

runner = RunCommand(GULP_EXE, inFile.stringForm())
runner.run(timeout=10)

with open('temp261016.gout', 'w') as outf:
    outf.write(runner.output)
xstalLTPGULP = Structure().fromGULPOutput(runner.output)

writeCIF('xstalLTPGULP.cif', xstalLTPGULP)

print 'gulp LTP'
print xstalLTPGULP.unitCell.lengths
LTPDict['A'][0], LTPDict['B'][0], LTPDict['C'][0] = xstalLTPGULP.unitCell.lengths[0], xstalLTPGULP.unitCell.lengths[1], xstalLTPGULP.unitCell.lengths[2]
for el in ['P', 'Ti', 'Li']:
    meanList = []
    for si in xstalLTPGULP.speciesMatchIndices(targetSpecies = Species(element = el)):
        print el, si, xstalLTPGULP.nearestNeighbourDistance(speciesListIndex = si)
        meanList.append(xstalLTPGULP.nearestNeighbourDistance(speciesListIndex = si))
    LTPDict[el][0] = np.mean(meanList)

xstalLTPVasp = VaspXML().finalStructure(open('vasprun_LTP_6.xml', 'r').read())
print 'vasp LTP'
print xstalLTPVasp.unitCell.lengths
LTPDict['A'][1], LTPDict['B'][1], LTPDict['C'][1] = xstalLTPVasp.unitCell.lengths[0], xstalLTPVasp.unitCell.lengths[1], xstalLTPVasp.unitCell.lengths[2]
for el in ['Li', 'Ti', 'P']:
    meanList = []
    for si in xstalLTPVasp.speciesMatchIndices(targetSpecies = Species(element = el)):
        print el, si, xstalLTPVasp.nearestNeighbourDistance(speciesListIndex = si)
        meanList.append(xstalLTPVasp.nearestNeighbourDistance(speciesListIndex = si))
    LTPDict[el][1] = xstalLTPVasp.nearestNeighbourDistance(speciesListIndex = si)



print LTPDict
