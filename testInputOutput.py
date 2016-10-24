import unittest
#from gulpBasics import *
#from inputOutput import *
import numpy as np

''' Tests only work on dhcs computer atm- make portable if others ever want them  '''

class TestCIFFileReading(unittest.TestCase):
    def test_tidyStringNumber(self):
        from cifIO import tidyLabelWithDigits
        self.assertTrue(tidyLabelWithDigits('Ti3') == 'Ti')
        self.assertTrue(tidyLabelWithDigits('Ti') == 'Ti')
        self.assertTrue(tidyLabelWithDigits('1111111Ti3') == '')

    def test_readCIF(self):
        self.assertTrue(True)

    def test_getUnitCellCIF(self):
        from cifIO import getUnitCellCIF
        cell = getUnitCellCIF('/Users/cmdc2-extra/Work/icsdFiles/EntryWithCollCode193041.cif')
        np.testing.assert_array_equal(cell.angles,  np.array([90., 90., 120.]))
        np.testing.assert_array_equal(cell.lengths, np.array([8.39187, 8.39187, 22.73297]))

    def test_getSpeciesListCIF(self):
        from cifIO import getSpeciesListCIF
        speciesList = getSpeciesListCIF('/Users/cmdc2-extra/Work/icsdFiles/EntryWithCollCode193041.cif')
        self.assertTrue(len(speciesList) == 11)
        self.assertTrue(speciesList[0].label   == 'Li1')
        self.assertTrue(speciesList[2].element == 'Ti')
        np.testing.assert_array_equal(speciesList[-1].fracCoord, np.array([0.9126, 0.1442, 0.2995]))

    def test_getSymmetryGroup(self):
        from cifIO import getSymmetryGroupCIF
        symmetryGroup = getSymmetryGroupCIF('/Users/cmdc2-extra/Work/icsdFiles/EntryWithCollCode193041.cif')
        self.assertTrue(symmetryGroup.labelHM == 'R -3 H')
        self.assertTrue(symmetryGroup.number  == 148)

class TestGULPFileReading(unittest.TestCase):
    def test_readOutputCompOpti(self):
        from gulpIO import OutputFileGULP
        outputFile = OutputFileGULP()
        outputDict = outputFile.readOutputCompOpti(open('testFiles/LTP.gout', 'r').read())
        np.testing.assert_array_equal(outputDict['newAngles'] , np.array([90., 90., 120.]))
        np.testing.assert_array_equal(outputDict['newLengths'], np.array([7.266685, 7.266685, 19.999402]))
        np.testing.assert_array_equal(outputDict['newCoords'] , np.array([[0., 0., 0.],
                                                                          [0., 0., 0.146214],
                                                                          [0.279077, 0., 0.250000],
                                                                          [0.156684, 0.963000, 0.199072],
                                                                          [0.203860, 0.173778, 0.100096]]))


    def test_readOutputCompOpti(self):
        from gulpIO import OutputFileGULP
        np.testing.assert_almost_equal(OutputFileGULP().latticeEnergy(open('testFiles/exampleLATP.gout', 'r').read()),
                                       -4887.67113227)

        allEnergies = OutputFileGULP().listAllLatticeEnergies(open('testFiles/exampleLATP.gout', 'r').read())
        np.testing.assert_almost_equal(allEnergies[0], 
                                       -4887.67113227)
        np.testing.assert_almost_equal(allEnergies[1], 
                                       -4987.99610200)

    def test_makeStructureGULP(self):
        from coreClasses import Structure
        testStructure = Structure().fromGULPOutput(open('testFiles/exampleLATP.gout', 'r').read())
        np.testing.assert_array_almost_equal(testStructure.unitCell.vectors,
                                             np.array([[ 8.233124,   -0.000000,   -0.000000],
                                                       [-4.116562,    7.130095,   -0.000000],
                                                       [-0.000000,   -0.000000,   20.985017]]))

        np.testing.assert_array_almost_equal(testStructure.speciesList[-1].fracCoord,
                                             np.array([0.524320,    0.831077,    0.756371]))
        self.assertTrue(testStructure.speciesList[-1].core[:4] == 'shel')
        self.assertTrue(testStructure.speciesList[-1].element == 'O')


class TestGULPFileWriting(unittest.TestCase):
    def test_outputCIF(self):
        ''' simply test that a cif file can be made
            ideally, this would be extended to run GULP in a 
            pipe, and return the file as a string '''
        from gulpIO import InputFileGULP
        infile = InputFileGULP()
        infile.outputInstructions = {'cif': 'test.cif'}
        self.assertTrue(any([l == 'output cif test.cif' for l in infile.stringForm().split("\n")]))


class TestDLPOLYOutput(unittest.TestCase):
    def testDataFrameOutput(self):
        ''' Assume have pandas '''
        from dlpolyIO import dlpolyOutput
        import pandas as pd
        egOutput  = dlpolyOutput(outputString = open('testFiles/OUTPUT310816', 'r').read())
        dataFrame = egOutput.dataFrameOutput('temp_tot', 'temp_shl')
#        print dataFrame.loc[:10]
        self.assertEqual(dataFrame.loc[0]['step'], 1)
        self.assertEqual(dataFrame.loc[0]['cpu-time'], 0.36)
        self.assertEqual(dataFrame.loc[0]['temp_shlAverage'], 2661.0)

class TestDLPOLYInput(unittest.TestCase):
    def testInitFromConfig(self):
        '''  '''
        from dlpolyIO import dlpolyInput
        newInput = dlpolyInput().initFromCONFIG('testFiles/exampleCONFIG', levcfg = 0)
        self.assertEqual(newInput.levcfg, 0)
        self.assertEqual(len(newInput.parentStructure.speciesList), 9000)
        self.assertTrue( newInput.speciesList[0].cartVelocity is None )
        np.testing.assert_array_almost_equal(newInput.speciesList[0].cartCoord,
                                             np.array([-0.2470917286E-01,   -0.2350931424,         10.40699720]))

if __name__ == '__main__':
    unittest.main()

