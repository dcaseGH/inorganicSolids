import unittest
#from gulpBasics import *
#from inputOutput import *
import numpy as np

''' Tests only work on dhcs computer atm- make portable if others ever want them  '''

class TestLiTiPO4(unittest.TestCase):
    def test_wholeMethodGettingDisplacements(self):
        ''' This is data used when fitting the LiTiPO4 potential 
            potentials are initial guess (from Islam's papers)
            structures checked in Mercury '''
        from coreClasses import Structure, VBuckingham, VSpring, VThreeBody, Species
        from gulpIO import InputFileGULP, OutputFileGULP
        from hardcode import GULP_EXE, CIF_FILE_DIRECTORY
        from subprocessHandling import RunCommand

        originalCIF = CIF_FILE_DIRECTORY + 'EntryWithCollCode95979.cif'
        xstal = Structure()
        xstal.inputCIF(originalCIF)
        qTiCore = 4.0
        qOCore  = -2.86
        xstal.addShells(Species(element = 'Ti', charge = qTiCore))
        xstal.addShells(Species(element = 'O',  charge = qOCore))
        xstal.setCharges([Species(element = 'Ti', core='core', charge = 4 - qTiCore),
                          Species(element = 'P',  core='core', charge = 5.),
                          Species(element = 'Li', core='core', charge = 1.),
                          Species(element = 'O',  core='core', charge = -2 - qOCore)])

        currentParam = {'TiOA': 5111.7,
                        'TiOrho': 0.26,
                        'LiOA': 632.1018,
                        'LiOrho': 0.2906,
                        'OOA': 22764.30,
                        'OOrho': 0.1490,
                        'OOC6': 27.89}

        inFile = InputFileGULP(fileName = 'LiTiP04',
                               title    = 'use to fit buckingham',
                               parentStructure = xstal,
                               keywords   = ['conp', 'opti', 'comp', 'phonon'],
                               potentials = [
                VBuckingham(species1 = Species(element = 'O', core='shel'),
                            species2 = Species(element = 'O', core='shel'),
                            A        = currentParam['OOA'],
                            rho      = currentParam['OOrho'],
                            C6       = currentParam['OOC6']),
                VBuckingham(species1 = Species(element = 'P'),
                            species2 = Species(element = 'O', core='shel'),
                            A        = 897.2648,
                            rho      = 0.3577,
                            C6       = 0.),
                VBuckingham(species1 = Species(element = 'Li'),
                            species2 = Species(element = 'O', core='shel'),
                            A        = currentParam['LiOA'],
                            rho      = currentParam['LiOrho'],
                            C6       = 0.),
                VBuckingham(species1 = Species(element = 'Ti', core='shel'),
                            species2 = Species(element = 'O',  core='shel'),
                            A        = currentParam['TiOA'],
                            rho      = currentParam['TiOrho'],
                            C6       = 0.0),
                VThreeBody( species1 = Species(element = 'P'),
                            species2 = Species(element = 'O', core='shel'),
                            species3 = Species(element = 'O', core='shel'),
                            K        = 1.322626,
                            theta0   = 109.47,
                            cut12    = 2.,
                            cut13    = 2.,
                            cut23    = 4.),
                VSpring(    species1 = Species(element = 'O'),
                            K        = 65.),
                VSpring(    species1 = Species(element = 'Ti'),
                            K        = 3140000.0)
                ])

        runner = RunCommand(GULP_EXE, inFile.stringForm())
        runner.run(timeout=10)
        gulpOutput = OutputFileGULP.readOutputCompOpti(runner.output)

        ''' Check that it's reading fractional coordinates and changing system  '''
        ''' Checked by hand for example that P atom moved 0.2903, 0 , 0.25 -> 0.29565, 0., 0.25 in Mercury '''
        np.testing.assert_array_almost_equal(np.array([ 0.29565, 0 , 0.25 ]),
                                             gulpOutput['newCoords'][2],
                                             decimal = 5)
#        print runner.output

        originalDisps, newDisps = xstal.customDisplacementsGULP(gulpOutput, originalCIF)

        ''' Check by hand - a=8.5110, b=8.5110, c=20.843, LiO=2.266,
                            TiO=1.880, PO=1.527 '''
        np.testing.assert_array_almost_equal(np.array([  8.511     ,   8.511     ,  20.843     ,   2.26551777,
                                                         1.88011998,   1.52709768,   1.5272784 ,   1.52709768]),
                                             originalDisps)
        ''' Check by hand - a=b=8.10874 c = 21.363, LiO = 2.278,
                            TiO=1.834, PO = 1.538 and 1.534 (distorted under minimisation) '''
        np.testing.assert_array_almost_equal(np.array([  8.108746  ,   8.108746  ,  21.363248  ,   2.2784803 ,
                                                         1.83389411,   1.53400978,   1.53400978,   1.53820873]),
                                             newDisps)


if __name__ == '__main__':
    unittest.main()

