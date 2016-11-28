import unittest
#from gulpBasics import *
#from inputOutput import *
import numpy as np

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


if __name__ == '__main__':
    unittest.main()

