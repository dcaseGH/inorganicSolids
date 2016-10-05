import numpy as np
from sobol_lib import i4_sobol

def sobolInterstitialSites(npts = 1, seedOffset = 0, inputFile = None, newFracCoord = None, cheshireMatrix = np.identity(3)):
    ''' Untested- if get odd cheshireMatrix, do tests'''
    # != 3 for some screw axes e.g. Pmm2 etc
    nDim = 3
    return np.array([np.dot(cheshire, i4_sobol(nDim, seed)[0]) for seed in xrange(seedOffset, seedOffset + npts)])


def calcEnergyDefect(inputFile           = None,
                     gulpCommand         = None,
                     timeout             = 10):

    ''' Run general calcn but return information about defects
        inputFile must be a new instance of a gulp input file '''

    from gulpIO import OutputFileGULP
    from subprocessHandling import RunCommand
    from coreClasses import Defect
    import time

    ''' input file should know all about defect '''
    runner = RunCommand(gulpCommand, inputFile.stringForm())
    try:
        runner.run(timeout = timeout)
    except:
        print "GULP failed for %s insertion"%newFracCoord

    try:
        defectEnergies = OutputFileGULP.defectEnergy(runner.output)
        return {'gulpInput': inputFile,
                'initialDefectEnergy': defectEnergies[0],
                'finalDefectEnergy': defectEnergies[1],
                'timeStamp': time.time(),
                'valid': True
                }
    except:
        return {'gulpInput': inputFile, 
                'initialDefectEnergy': None,
                'finalDefectEnergy': None,
                'timeStamp': time.time(),
                'valid': False
                }

def calcEnergyInterstitialSite(inputFile           = None,
                               interstitialSpecies = None,
                               newFracCoord        = None,
                               cheshireMatrix      = np.identity(3),
                               gulpCommand         = None,
                               timeout             = 10):

    ''' Set the input file up to run the calculation that you want
        i.e. it has structure, potentials etc and keywords for defect calcn 
        interstitial species is instance of Species
        newFracCoord is np.array '''

    from gulpIO import OutputFileGULP
    from subprocessHandling import RunCommand
    from coreClasses import Defect
    import time

    ''' Put one extra atom at newCoord '''
    interstitialSpecies.fracCoord = newFracCoord
    inputFile.defectSpecies = [Defect(species = interstitialSpecies,
                                      defectType = 'interstitial')]
    inputFile.defectDetails['centre'] = newFracCoord

    runner = RunCommand(gulpCommand, inputFile.stringForm())
    try:
        runner.run(timeout = timeout)
    except:
        print "GULP failed for %s insertion"%newFracCoord

    try:
        defectEnergies = OutputFileGULP.defectEnergy(runner.output)
        return {'gulpInput': inputFile,
                'initialDefectEnergy': defectEnergies[0],
                'finalDefectEnergy': defectEnergies[1],
                'timeStamp': time.time(),
                'valid': True
                }
    except:
        return {'gulpInput': inputFile, 
                'initialDefectEnergy': None,
                'finalDefectEnergy': None,
                'timeStamp': time.time(),
                'valid': False
                }

