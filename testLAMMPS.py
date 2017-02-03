import unittest
import numpy as np

class TestLAMMPSOutput(unittest.TestCase):
    def testDataFrameOutput(self):

        from lammpsIO import LAMMPSOutput
        from coreClasses import Species
        print 'aitnian'

        struc = LAMMPSOutput.structureFromDataFile('testFiles/example_data.lmp',
                                                   speciesDict = {1: Species(element = 'B',
                                                                             mass    = 10.811,
                                                                             charge  = 3.),
                                                                  2: Species(element = 'Li',
                                                                             mass    = 6.941,
                                                                             charge  = 1.),
                                                                  3: Species(element = 'O',
                                                                             mass    = 15.9994,
                                                                             charge  = -2.),
                                                                  4: Species(element = 'V',
                                                                             mass    = 50.9415,
                                                                             charge  = 5.)})
        self.assertEqual(len(struc.speciesList), 3000)
        #these are from the log file- unit cell lengths
        np.testing.assert_array_almost_equal(struc.unitCell.lengths,
                                             np.array([32.909985,    31.364635,     33.54656]))

if __name__ == '__main__':
    unittest.main()
